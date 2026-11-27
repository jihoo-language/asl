use std::cell::{Cell, RefCell};
use std::collections::{btree_map::Entry, BTreeMap};
use std::rc::Rc;

use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::odd_ext;
use crate::prime_linear::{self, add_scaled, Basis, Coordinate};
use crate::rv_as_odd::Cover;

pub(crate) type Polynomial = Rc<[(u128, u32)]>;
pub(crate) type Column = Rc<[Coordinate]>;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Cartier,
    Iterate(usize),
    Translate(i128),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Cartier => "rvASOddCartier",
            Self::Iterate(_) => "rvASOddCartierIter",
            Self::Translate(_) => "rvASOddTranslate",
        }
    }
}

pub(crate) struct Cartier {
    prime: Prime,
    // Basis positions in each row of fixed y-exponent, including the end.
    offsets: Vec<usize>,
    factorials: Vec<u32>,
    inverse_factorials: Vec<u32>,
    powers: RefCell<Vec<Polynomial>>,
    columns: RefCell<BTreeMap<usize, Column>>,
    cells: Cell<usize>,
}

impl Cartier {
    pub(crate) fn build(model: &mut Model, cover: &Cover) -> Result<Self> {
        let prime = cover.prime();
        let p = prime.get() as u128;
        let d = cover.degree() as u128;
        let mut offsets = vec![0usize];
        if let Some(bound) = (2 * cover.genus() as u128).checked_sub(2) {
            for r in 0..=p - 2 {
                model.spend(1)?;
                let Some(x_bound) = bound.checked_sub(d * r) else {
                    break;
                };
                let count = usize::try_from(x_bound / p + 1).map_err(|_| basis_error())?;
                offsets.push(
                    offsets
                        .last()
                        .unwrap()
                        .checked_add(count)
                        .ok_or_else(basis_error)?,
                );
            }
        }
        if offsets.last() != Some(&cover.genus()) {
            return Err(basis_error());
        }
        // The Čech overlap basis can have y-exponent p-1.
        let max_r = (prime.get() - 1) as usize;
        let mut factorials = vec![1; max_r + 1];
        let mut inverse_factorials = vec![1; max_r + 1];
        model.spend(2 * max_r)?;
        for r in 1..=max_r {
            factorials[r] = prime.mul(factorials[r - 1], r as u32);
        }
        inverse_factorials[max_r] = prime_linear::invert(model, prime, factorials[max_r])?;
        for r in (1..=max_r).rev() {
            inverse_factorials[r - 1] = prime.mul(inverse_factorials[r], r as u32);
        }
        let cells = 9 + offsets.len() + factorials.len() + inverse_factorials.len();
        model.allocate(cells)?;
        Ok(Self {
            prime,
            offsets,
            factorials,
            inverse_factorials,
            powers: RefCell::new(vec![Rc::from([(0, 1)])]),
            columns: RefCell::new(BTreeMap::new()),
            cells: Cell::new(cells),
        })
    }

    fn size(&self) -> usize {
        *self.offsets.last().unwrap()
    }

    pub(crate) fn exponent(&self, coordinate: usize) -> (usize, usize) {
        let r = self.offsets.partition_point(|&offset| offset <= coordinate) - 1;
        (coordinate - self.offsets[r], r)
    }

    pub(crate) fn coordinate(&self, b: usize, r: usize) -> Option<usize> {
        let start = *self.offsets.get(r)?;
        let end = *self.offsets.get(r.checked_add(1)?)?;
        (b < end - start).then(|| start + b)
    }

    pub(crate) fn choose(&self, r: usize, j: usize) -> u32 {
        self.prime.mul(
            self.factorials[r],
            self.prime
                .mul(self.inverse_factorials[j], self.inverse_factorials[r - j]),
        )
    }

    fn retain(&self, model: &mut Model, cells: usize) -> Result<()> {
        let total = self.cells.get() + cells;
        check_storage(total)?;
        model.allocate(cells)?;
        self.cells.set(total);
        Ok(())
    }

    pub(crate) fn generator(&self, model: &mut Model, b: usize, r: usize) -> Result<Value> {
        let coordinate = self.coordinate(b, r).ok_or_else(|| {
            Error::plain(format!(
                "x^{b} y^{r} dx is outside the regular differential basis"
            ))
        })?;
        odd_ext::persist(model, vec![(coordinate, 1)])
    }

    pub(crate) fn basis(&self, model: &mut Model) -> Result<Vec<(usize, usize)>> {
        check_storage(self.size().saturating_mul(2))?;
        model.spend(self.size())?;
        Ok(self
            .offsets
            .windows(2)
            .enumerate()
            .flat_map(|(r, row)| (0..row[1] - row[0]).map(move |b| (b, r)))
            .collect())
    }

    pub(crate) fn power(
        &self,
        model: &mut Model,
        cover: &Cover,
        exponent: usize,
    ) -> Result<Polynomial> {
        let mut powers = self.powers.borrow_mut();
        while powers.len() <= exponent {
            let previous = powers.last().unwrap();
            let terms = cover.terms();
            model.spend(1usize.saturating_add(previous.len().saturating_mul(terms.len())))?;
            let mut next = BTreeMap::new();
            for &(i, a) in previous.iter() {
                for &(j, b) in terms {
                    let value = self.prime.mul(a, b);
                    model.spend(1 + next.len().max(1).ilog2() as usize)?;
                    match next.entry(i + j as u128) {
                        Entry::Occupied(mut entry) => {
                            let sum = self.prime.add(*entry.get(), value);
                            if sum == 0 {
                                entry.remove();
                            } else {
                                *entry.get_mut() = sum;
                            }
                        }
                        Entry::Vacant(entry) => {
                            entry.insert(value);
                            check_storage(self.cells.get() + 8 * next.len())?;
                        }
                    }
                }
            }
            self.retain(model, 3 + 4 * next.len())?;
            powers.push(next.into_iter().collect());
        }
        Ok(powers[exponent].clone())
    }

    pub(crate) fn column(
        &self,
        model: &mut Model,
        cover: &Cover,
        coordinate: usize,
    ) -> Result<Column> {
        if let Some(column) = self.columns.borrow().get(&coordinate) {
            model.spend(1)?;
            return Ok(column.clone());
        }
        let (b, r) = self.exponent(coordinate);
        let result = self
            .monomial_image(model, cover, b as i128, r)?
            .into_iter()
            .map(|(b, r, coefficient)| {
                let b = usize::try_from(b).map_err(|_| basis_error())?;
                let destination = self.coordinate(b, r).ok_or_else(basis_error)?;
                Ok((destination, coefficient))
            })
            .collect::<Result<Vec<_>>>()?;
        self.retain(model, 5 + 2 * result.len())?;
        let result: Column = result.into();
        self.columns.borrow_mut().insert(coordinate, result.clone());
        Ok(result)
    }

    pub(crate) fn monomial_image(
        &self,
        model: &mut Model,
        cover: &Cover,
        b: i128,
        r: usize,
    ) -> Result<Vec<(i128, usize, u32)>> {
        let p = i128::from(self.prime.get());
        let mut result = Vec::new();
        // Expand y^r = (y^p - f)^r. This formula on Laurent monomials
        // gives both the regular Cartier matrix and the Čech chart maps.
        for j in 0..=r {
            let power = self.power(model, cover, r - j)?;
            let mut scalar = self.choose(r, j);
            if !(r - j).is_multiple_of(2) {
                scalar = self.prime.negate(scalar);
            }
            model.spend(1 + power.len())?;
            for &(e, coefficient) in power.iter() {
                // d fits in usize and r<p is bounded by the cyclotomic order.
                let exponent = e as i128 + b;
                if exponent.rem_euclid(p) == p - 1 {
                    result.push((
                        exponent.div_euclid(p),
                        j,
                        self.prime.mul(coefficient, scalar),
                    ));
                    check_storage(self.cells.get() + 4 * result.len())?;
                }
            }
        }
        Ok(result)
    }

    fn step(
        &self,
        model: &mut Model,
        cover: &Cover,
        source: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let mut result = Vec::new();
        for &(index, scalar) in source {
            let column = self.column(model, cover, index)?;
            result = add_scaled(model, self.prime, &result, &column, scalar)?;
        }
        Ok(result)
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        cover: &Cover,
        operation: Operation,
        source: Column,
    ) -> Result<Value> {
        if source.is_empty() {
            return Ok(odd_ext::zero());
        }
        let genus = self.size();
        match operation {
            Operation::Iterate(0) => Ok(Value::OddClass(source)),
            Operation::Cartier | Operation::Iterate(_) => {
                let iterations = match operation {
                    Operation::Iterate(count) => count,
                    _ => 1,
                };
                let mut result = self.step(model, cover, &source)?;
                for _ in 1..iterations.min(genus) {
                    if result.is_empty() {
                        break;
                    }
                    result = self.step(model, cover, &result)?;
                }
                if iterations >= genus && !result.is_empty() {
                    return Err(nilpotence_error());
                }
                odd_ext::persist(model, result)
            }
            Operation::Translate(scalar) => {
                let c = scalar.rem_euclid(i128::from(self.prime.get())) as u32;
                if c == 0 {
                    return Ok(Value::OddClass(source));
                }
                let mut result = Vec::new();
                for &(index, coefficient) in source.iter() {
                    let (b, r) = self.exponent(index);
                    let mut power = 1;
                    let mut column = Vec::with_capacity(r + 1);
                    for j in (0..=r).rev() {
                        let index = self.coordinate(b, j).ok_or_else(basis_error)?;
                        let scalar = self.prime.mul(self.choose(r, j), power);
                        column.push((index, scalar));
                        power = self.prime.mul(power, c);
                    }
                    model.spend(column.len())?;
                    column.reverse();
                    result = add_scaled(model, self.prime, &result, &column, coefficient)?;
                }
                odd_ext::persist(model, result)
            }
        }
    }

    pub(crate) fn rank(&self, model: &mut Model, cover: &Cover, power: usize) -> Result<usize> {
        let genus = self.size();
        if power == 0 || genus == 0 {
            return Ok(genus);
        }
        let power = power.min(genus);
        check_storage(self.cells.get().saturating_add(genus.saturating_mul(5)))?;
        model.spend(genus)?;
        let mut basis: Vec<_> = (0..genus).map(|i| vec![(i, 1)]).collect();
        let mut basis_cells = genus * 5;
        for _ in 0..power {
            let mut span = Basis::default();
            let mut next = Vec::new();
            let mut next_cells = 0;
            for vector in &basis {
                let image = self.step(model, cover, vector)?;
                check_storage(
                    self.cells.get() + basis_cells + next_cells + span.cells() + 2 * image.len(),
                )?;
                let rank = span.rank();
                span.insert(model, self.prime, image.clone())?;
                if span.rank() != rank {
                    next_cells += 3 + 2 * image.len();
                    check_storage(self.cells.get() + basis_cells + next_cells + span.cells())?;
                    next.push(image);
                }
            }
            if next.is_empty() {
                return Ok(0);
            }
            if next.len() == basis.len() {
                return Err(nilpotence_error());
            }
            basis = next;
            basis_cells = next_cells;
        }
        Ok(basis.len())
    }
}

fn basis_error() -> Error {
    Error::plain(
        "odd Artin–Schreier Cartier calculation disagrees with the regular differential basis",
    )
}

fn nilpotence_error() -> Error {
    Error::plain("odd Artin–Schreier Cartier operator has an unexpected nonzero stable image")
}
