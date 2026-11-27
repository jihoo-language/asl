use std::cell::{Cell, OnceCell, RefCell};
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::odd_ext;
use crate::prime_linear::{add_scaled, Basis, Coordinate};
use crate::rv_as_two::Cover;
use crate::rv_cartier_arithmetic::{retain, Binomial, Powers};

type Column = Rc<[Coordinate]>;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Cartier,
    Iterate(usize),
    Translate(i128),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Cartier => "rvASTwoCartier",
            Self::Iterate(_) => "rvASTwoCartierIter",
            Self::Translate(_) => "rvASTwoTranslate",
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Basis,
    Rank(usize),
    StableRank,
    ANumber,
    OrdinaryRanks,
    HasseWitt,
}

struct Row {
    low: i128,
    high: i128,
    start: usize,
}

struct Images {
    ranks: Vec<usize>,
    stable: Basis,
}

#[derive(Clone, Copy)]
enum Action {
    Cartier,
    CartierDifference,
    DeckDifference,
}

pub(crate) struct Cartier {
    prime: Prime,
    size: usize,
    rows: Vec<Row>,
    binomial: Binomial,
    polynomial: Vec<(i128, u32)>,
    powers: Powers,
    columns: RefCell<BTreeMap<usize, Column>>,
    images: OnceCell<Images>,
    cells: Cell<usize>,
}

impl Cartier {
    pub(crate) fn build(model: &mut Model, cover: &Cover) -> Result<Self> {
        let prime = cover.prime();
        let p = i128::from(prime.get());
        let (d, e) = cover.pole_orders();
        let (d, e) = (d as i128, e as i128);
        let mut rows = Vec::new();
        let mut size = 0usize;
        for j in 0..=p - 2 {
            model.spend(1)?;
            // All rows include dx/x times y^j. The lower and upper
            // bounds are the regularity inequalities at zero and infinity.
            let low = -((p - 1) * (e + 1) - e * j).div_euclid(p);
            let high = ((p - 1) * (d + 1) - 2 * p - d * j).div_euclid(p);
            if low > -1 || high < -1 {
                return Err(basis_error());
            }
            rows.push(Row {
                low,
                high,
                start: size,
            });
            let width = usize::try_from(high - low + 1).map_err(|_| basis_error())?;
            size = size.checked_add(width).ok_or_else(basis_error)?;
        }
        if size != cover.genus() {
            return Err(basis_error());
        }
        let binomial = Binomial::build(model, prime, (p - 1) as usize)?;
        let polynomial = cover.laurent_terms(model)?;
        let cells = 10 + 5 * rows.len() + 3 * polynomial.len() + binomial.cells();
        check_storage(cells)?;
        model.allocate(cells)?;
        Ok(Self {
            prime,
            size,
            rows,
            binomial,
            polynomial,
            powers: Powers::default(),
            columns: RefCell::default(),
            images: OnceCell::new(),
            cells: Cell::new(cells),
        })
    }

    pub(crate) fn exponent(&self, coordinate: usize) -> (i128, usize) {
        let j = self.rows.partition_point(|row| row.start <= coordinate) - 1;
        let row = &self.rows[j];
        (row.low + (coordinate - row.start) as i128, j)
    }

    fn coordinate(&self, i: i128, j: usize) -> Option<usize> {
        let row = self.rows.get(j)?;
        (row.low..=row.high)
            .contains(&i)
            .then(|| row.start + (i - row.low) as usize)
    }

    pub(crate) fn generator(&self, model: &mut Model, i: i128, j: usize) -> Result<Value> {
        let index = self.coordinate(i, j).ok_or_else(|| {
            Error::plain(format!(
                "x^{i} y^{j} dx is outside the two-pole regular differential basis"
            ))
        })?;
        odd_ext::persist(model, vec![(index, 1)])
    }

    fn basis(&self, model: &mut Model) -> Result<Vec<(i128, usize)>> {
        check_storage(self.size.saturating_mul(3))?;
        model.spend(self.size)?;
        Ok(self
            .rows
            .iter()
            .enumerate()
            .flat_map(|(j, row)| (row.low..=row.high).map(move |i| (i, j)))
            .collect())
    }

    fn column(&self, model: &mut Model, index: usize) -> Result<Column> {
        if let Some(column) = self.columns.borrow().get(&index) {
            model.spend(1)?;
            return Ok(column.clone());
        }
        let (i, j) = self.exponent(index);
        let image = self.binomial.cartier_monomial(
            model,
            i,
            j,
            |model, k| {
                self.powers.get(
                    model,
                    self.prime,
                    self.polynomial.iter().copied(),
                    k,
                    &self.cells,
                )
            },
            &self.cells,
        )?;
        let image = image
            .into_iter()
            .map(|(i, j, c)| {
                let dest = self.coordinate(i, j).ok_or_else(basis_error)?;
                Ok((dest, c))
            })
            .collect::<Result<Vec<_>>>()?;
        retain(model, &self.cells, 5 + 2 * image.len())?;
        let image: Column = image.into();
        self.columns.borrow_mut().insert(index, image.clone());
        Ok(image)
    }

    fn step(&self, model: &mut Model, support: &[Coordinate]) -> Result<Vec<Coordinate>> {
        let mut result = Vec::new();
        for &(index, coefficient) in support {
            let image = self.column(model, index)?;
            result = add_scaled(model, self.prime, &result, &image, coefficient)?;
        }
        Ok(result)
    }

    fn translate(
        &self,
        model: &mut Model,
        support: &[Coordinate],
        scalar: u32,
    ) -> Result<Vec<Coordinate>> {
        let mut result = Vec::new();
        for &(index, coefficient) in support {
            let (i, j) = self.exponent(index);
            let image = self.binomial.translate(model, j, scalar, |k| {
                self.coordinate(i, k).ok_or_else(basis_error)
            })?;
            result = add_scaled(model, self.prime, &result, &image, coefficient)?;
        }
        Ok(result)
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        operation: Operation,
        source: Column,
    ) -> Result<Value> {
        if source.is_empty() {
            return Ok(odd_ext::zero());
        }
        let iterations = match operation {
            Operation::Cartier => 1,
            Operation::Iterate(0) => return Ok(Value::OddClass(source)),
            Operation::Iterate(count) => count,
            Operation::Translate(c) => {
                let scalar = c.rem_euclid(i128::from(self.prime.get())) as u32;
                return if scalar == 0 {
                    Ok(Value::OddClass(source))
                } else {
                    let result = self.translate(model, &source, scalar)?;
                    odd_ext::persist(model, result)
                };
            }
        };
        // By degree g the image is stable. The deck difference has one
        // (p−1)-block there, with fixed vector dx/x. Cartier commutes with
        // the deck action and fixes dx/x, hence C−1 is nilpotent and C^p=1
        // on this image. This keeps arbitrarily high actual iterates short.
        let count = if iterations > self.size {
            self.size + (iterations - self.size) % self.prime.get() as usize
        } else {
            iterations
        };
        let mut result = self.step(model, &source)?;
        for _ in 1..count {
            if result.is_empty() {
                break;
            }
            result = self.step(model, &result)?;
        }
        odd_ext::persist(model, result)
    }

    fn image_basis(
        &self,
        model: &mut Model,
        source: Option<&Basis>,
        action: Action,
        scratch: usize,
    ) -> Result<Basis> {
        let units = if source.is_some() { 0 } else { self.size };
        check_storage(
            self.cells
                .get()
                .saturating_add(scratch)
                .saturating_add(units),
        )?;
        model.spend(units)?;
        let inputs = source
            .into_iter()
            .flat_map(|basis| basis.vectors(self.prime))
            .chain((0..units).map(|i| vec![(i, 1)]));
        let mut result = Basis::default();
        for input in inputs {
            let moved = match action {
                Action::Cartier | Action::CartierDifference => self.step(model, &input)?,
                Action::DeckDifference => self.translate(model, &input, 1)?,
            };
            let image = match action {
                Action::Cartier => moved,
                _ => add_scaled(model, self.prime, &moved, &input, self.prime.negate(1))?,
            };
            check_storage(self.cells.get() + scratch + result.cells() + 2 * image.len())?;
            result.insert(model, self.prime, image)?;
            check_storage(self.cells.get() + scratch + result.cells())?;
        }
        Ok(result)
    }

    fn rank(&self, model: &mut Model, power: usize) -> Result<usize> {
        if power == 0 {
            return Ok(self.size);
        }
        if let Some(data) = self.images.get() {
            return Ok(data
                .ranks
                .get(power)
                .copied()
                .unwrap_or_else(|| data.stable.rank()));
        }
        let mut ranks = vec![self.size];
        let mut current: Option<Basis> = None;
        for _ in 0..power.min(self.size) {
            let scratch = current.as_ref().map_or(0, |basis| basis.cells()) + ranks.len();
            let next = self.image_basis(model, current.as_ref(), Action::Cartier, scratch)?;
            let rank = next.rank();
            if Some(&rank) == ranks.last() {
                if rank != self.prime.get() as usize - 1 {
                    return Err(Error::plain(
                        "two-pole Cartier stable image has unexpected dimension",
                    ));
                }
                retain(model, &self.cells, 5 + ranks.len() + next.cells())?;
                self.images.get_or_init(|| Images {
                    ranks,
                    stable: next,
                });
                return Ok(rank);
            }
            ranks.push(rank);
            current = Some(next);
        }
        if power >= self.size {
            return Err(Error::plain("two-pole Cartier image failed to stabilize"));
        }
        Ok(*ranks.last().unwrap())
    }

    fn ordinary(&self, model: &mut Model) -> Result<&Basis> {
        self.rank(model, self.size)?;
        Ok(&self.images.get().unwrap().stable)
    }

    pub(crate) fn ordinary_generator(&self, model: &mut Model, index: usize) -> Result<Value> {
        let stable = self.ordinary(model)?;
        if index >= stable.rank() {
            return Err(Error::plain(format!(
                "ordinary generator index {index} is outside dimension {}",
                stable.rank()
            )));
        }
        model.spend(1 + index)?;
        let vector = stable.vectors(self.prime).nth(index).unwrap();
        odd_ext::persist(model, vector)
    }

    fn difference_ranks(
        &self,
        model: &mut Model,
        stable: &Basis,
        action: Action,
    ) -> Result<Vec<usize>> {
        let mut ranks = vec![stable.rank()];
        let mut current: Option<Basis> = None;
        while *ranks.last().unwrap() != 0 {
            let scratch = current.as_ref().map_or(0, |basis| basis.cells()) + ranks.len();
            let next = self.image_basis(
                model,
                Some(current.as_ref().unwrap_or(stable)),
                action,
                scratch,
            )?;
            if next.rank() >= *ranks.last().unwrap() {
                return Err(Error::plain(
                    "two-pole ordinary difference failed to be nilpotent",
                ));
            }
            ranks.push(next.rank());
            current = Some(next);
        }
        Ok(ranks)
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, query: Query) -> Result<String> {
        match query {
            Query::Basis => Ok(format!("rvASTwoBasis[{name}] = {:?}", self.basis(model)?)),
            Query::Rank(k) => Ok(format!(
                "rvASTwoCartierRank[{name}, {k}] = {}",
                self.rank(model, k)?
            )),
            Query::StableRank => Ok(format!(
                "rvASTwoStableRank[{name}] = {}",
                self.ordinary(model)?.rank()
            )),
            Query::ANumber => Ok(format!(
                "rvASTwoaNumber[{name}] = {}",
                self.size - self.rank(model, 1)?
            )),
            Query::OrdinaryRanks => {
                let stable = self.ordinary(model)?;
                let deck = self.difference_ranks(model, stable, Action::DeckDifference)?;
                let cartier = self.difference_ranks(model, stable, Action::CartierDifference)?;
                Ok(format!(
                    "rvASTwoOrdinaryRanks[{name}] = T-1 {deck:?}, C-1 {cartier:?}"
                ))
            }
            Query::HasseWitt => {
                let stable = self.ordinary(model)?;
                self.difference_ranks(model, stable, Action::CartierDifference)?;
                // The actual stable restriction is unipotent, while C on
                // the complement quotient is nilpotent. This gives det(1−TC).
                let dimension = stable.rank();
                model.spend(1 + dimension)?;
                let coefficients: Vec<_> = (0..=dimension)
                    .map(|i| {
                        let coefficient = self.binomial.choose(dimension, i);
                        if i.is_multiple_of(2) {
                            coefficient
                        } else {
                            self.prime.negate(coefficient)
                        }
                    })
                    .collect();
                Ok(format!("rvASTwoHasseWitt[{name}] = {coefficients:?}"))
            }
        }
    }
}

fn basis_error() -> Error {
    Error::plain("two-pole Cartier calculation disagrees with the regular differential basis")
}
