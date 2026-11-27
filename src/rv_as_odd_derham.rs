use std::cell::{Cell, RefCell};
use std::collections::{btree_map::Entry, BTreeMap, VecDeque};
use std::rc::Rc;

use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::odd_ext;
use crate::prime_linear::{add_scaled, Basis, Coordinate};
use crate::rv_as_odd::Cover;
use crate::rv_as_odd_cartier::{Cartier, Column};

// Laurent x-exponents need room for the original usize degree times p.
// The cyclotomic order bounds p, and every y exponent is less than p.
type Monomial = (i128, usize);
type Laurent = BTreeMap<Monomial, u32>;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Include,
    Frobenius,
    Verschiebung,
    Translate(i128),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Include => "rvASOddDeRhamIncl",
            Self::Frobenius => "rvASOddF",
            Self::Verschiebung => "rvASOddV",
            Self::Translate(_) => "rvASOddDeRhamTranslate",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Map {
    Frobenius,
    Verschiebung,
}

#[derive(Default)]
struct Triple {
    u: Laurent,
    v: Laurent,
    overlap: Laurent,
}

pub(crate) struct DeRham {
    prime: Prime,
    genus: usize,
    dimension: usize,
    cartier: Rc<Cartier>,
    columns: RefCell<BTreeMap<(Map, usize), Column>>,
    cells: Cell<usize>,
}

impl DeRham {
    pub(crate) fn build(model: &mut Model, cover: &Cover) -> Result<Self> {
        let genus = cover.genus();
        let dimension = genus.checked_mul(2).ok_or_else(dimension_error)?;
        Ok(Self {
            prime: cover.prime(),
            genus,
            dimension,
            cartier: cover.cartier(model)?,
            columns: RefCell::new(BTreeMap::new()),
            cells: Cell::new(0),
        })
    }

    pub(crate) fn generator(&self, model: &mut Model, index: usize) -> Result<Value> {
        if index >= self.dimension {
            return Err(Error::plain(format!(
                "rvASOddDeRhamGen index {index} is out of range for dimension {}",
                self.dimension
            )));
        }
        odd_ext::persist(model, vec![(index, 1)])
    }

    fn p(&self) -> usize {
        self.prime.get() as usize
    }

    fn gap(&self, index: usize) -> Monomial {
        let (b, r) = self.cartier.exponent(index);
        (-(b as i128) - 1, self.p() - 1 - r)
    }

    fn lift(&self, model: &mut Model, cover: &Cover, h: Monomial) -> Result<(Laurent, Laurent)> {
        let mut u = Laurent::new();
        let mut v = Laurent::new();
        for (term, coefficient) in self.derivative(model, cover, &Laurent::from([(h, 1)]))? {
            if term.0 >= 0 {
                u.insert(term, coefficient);
            } else {
                v.insert(term, self.prime.negate(coefficient));
            }
        }
        Ok((u, v))
    }

    fn derivative(&self, model: &mut Model, cover: &Cover, input: &Laurent) -> Result<Laurent> {
        let p = self.p() as i128;
        let mut output = Laurent::new();
        for (&(b, r), &coefficient) in input {
            let x_coefficient = self.prime.mul(coefficient, b.rem_euclid(p) as u32);
            add_term(model, self.prime, &mut output, (b - 1, r), x_coefficient)?;
            if r != 0 {
                // The affine equation gives dy = -f'(x) dx.
                let y_coefficient = self.prime.negate(self.prime.mul(coefficient, r as u32));
                model.spend(cover.terms().len())?;
                for &(power, c) in cover.terms() {
                    if power != 0 {
                        let df = self.prime.mul(self.prime.reduce(power), c);
                        add_term(
                            model,
                            self.prime,
                            &mut output,
                            (b + power as i128 - 1, r - 1),
                            self.prime.mul(y_coefficient, df),
                        )?;
                    }
                }
            }
        }
        Ok(output)
    }

    fn triple(
        &self,
        model: &mut Model,
        cover: &Cover,
        coordinates: &[Coordinate],
    ) -> Result<Triple> {
        let mut result = Triple::default();
        for &(index, coefficient) in coordinates {
            if index < self.genus {
                let (b, r) = self.cartier.exponent(index);
                let term = (b as i128, r);
                add_term(model, self.prime, &mut result.u, term, coefficient)?;
                add_term(model, self.prime, &mut result.v, term, coefficient)?;
            } else {
                let h = self.gap(index - self.genus);
                let (u, v) = self.lift(model, cover, h)?;
                add_form(model, self.prime, &mut result.u, &u, coefficient)?;
                add_form(model, self.prime, &mut result.v, &v, coefficient)?;
                add_term(model, self.prime, &mut result.overlap, h, coefficient)?;
            }
        }
        Ok(result)
    }

    fn decode(
        &self,
        model: &mut Model,
        cover: &Cover,
        mut triple: Triple,
    ) -> Result<Vec<Coordinate>> {
        let mut exact_u = Laurent::new();
        let mut exact_v = Laurent::new();
        let mut result = Vec::new();
        for (term @ (b, r), coefficient) in triple.overlap {
            model.spend(1)?;
            if b >= 0 {
                exact_u.insert(term, coefficient);
            } else if self.p() as i128 * b + cover.degree() as i128 * r as i128 <= 0 {
                exact_v.insert(term, coefficient);
            } else {
                // The remaining monomials are precisely the gap basis. The
                // sign is fixed by U - V = dh and h = x^(-b-1)y^(p-1-r).
                let x = usize::try_from(-b - 1).map_err(|_| basis_error())?;
                let index = self
                    .cartier
                    .coordinate(x, self.p() - 1 - r)
                    .ok_or_else(basis_error)?;
                result.push((self.genus + index, coefficient));
                let (u, v) = self.lift(model, cover, term)?;
                let negative = self.prime.negate(coefficient);
                add_form(model, self.prime, &mut triple.u, &u, negative)?;
                add_form(model, self.prime, &mut triple.v, &v, negative)?;
            }
        }
        let du = self.derivative(model, cover, &exact_u)?;
        let dv = self.derivative(model, cover, &exact_v)?;
        add_form(model, self.prime, &mut triple.u, &du, self.prime.negate(1))?;
        add_form(model, self.prime, &mut triple.v, &dv, 1)?;
        model.spend(triple.u.len() + triple.v.len())?;
        if triple.u != triple.v {
            return Err(Error::plain(
                "odd Artin–Schreier Čech charts do not glue after reduction",
            ));
        }
        for ((b, r), coefficient) in triple.u {
            let x = usize::try_from(b).map_err(|_| basis_error())?;
            let index = self.cartier.coordinate(x, r).ok_or_else(basis_error)?;
            result.push((index, coefficient));
        }
        check_storage(self.cells.get() + 2 * result.len())?;
        model.spend(result.len() * (1 + result.len().max(1).ilog2() as usize))?;
        result.sort_unstable();
        Ok(result)
    }

    fn pth_power_gap(&self, model: &mut Model, cover: &Cover, (b, r): Monomial) -> Result<Laurent> {
        let mut result = Laurent::new();
        let x = b * self.p() as i128;
        // h^p uses y^p = y+f; it has zero differential in characteristic p.
        for j in 0..=r {
            let scale = self.cartier.choose(r, j);
            let coefficients = self.cartier.power(model, cover, r - j)?;
            for &(power, coefficient) in coefficients.iter() {
                add_term(
                    model,
                    self.prime,
                    &mut result,
                    (x + power as i128, j),
                    self.prime.mul(coefficient, scale),
                )?;
            }
        }
        Ok(result)
    }

    fn cartier_form(&self, model: &mut Model, cover: &Cover, form: &Laurent) -> Result<Laurent> {
        let mut result = Laurent::new();
        for (&(b, r), &coefficient) in form {
            for (x, j, c) in self.cartier.monomial_image(model, cover, b, r)? {
                add_term(
                    model,
                    self.prime,
                    &mut result,
                    (x, j),
                    self.prime.mul(coefficient, c),
                )?;
            }
        }
        Ok(result)
    }

    fn column(
        &self,
        model: &mut Model,
        cover: &Cover,
        operation: Map,
        index: usize,
    ) -> Result<Column> {
        if index < self.genus {
            return match operation {
                Map::Frobenius => Ok([].into()),
                Map::Verschiebung => self.cartier.column(model, cover, index),
            };
        }
        if let Some(column) = self.columns.borrow().get(&(operation, index)) {
            model.spend(1)?;
            return Ok(column.clone());
        }
        let gap = self.gap(index - self.genus);
        let triple = match operation {
            Map::Frobenius => Triple {
                overlap: self.pth_power_gap(model, cover, gap)?,
                ..Triple::default()
            },
            Map::Verschiebung => {
                let (u, v) = self.lift(model, cover, gap)?;
                Triple {
                    u: self.cartier_form(model, cover, &u)?,
                    v: self.cartier_form(model, cover, &v)?,
                    ..Triple::default()
                }
            }
        };
        let output = self.decode(model, cover, triple)?;
        let cells = 6 + 2 * output.len();
        let total = self.cells.get() + cells;
        check_storage(total)?;
        model.allocate(cells)?;
        self.cells.set(total);
        let output: Column = output.into();
        self.columns
            .borrow_mut()
            .insert((operation, index), output.clone());
        Ok(output)
    }

    fn image(
        &self,
        model: &mut Model,
        cover: &Cover,
        operation: Map,
        input: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let mut output = Vec::new();
        for &(index, scalar) in input {
            let column = self.column(model, cover, operation, index)?;
            output = add_scaled(model, self.prime, &output, &column, scalar)?;
        }
        Ok(output)
    }

    fn translate(&self, model: &mut Model, input: &Laurent, scalar: u32) -> Result<Laurent> {
        let mut result = Laurent::new();
        for (&(b, r), &coefficient) in input {
            let mut power = coefficient;
            for j in (0..=r).rev() {
                add_term(
                    model,
                    self.prime,
                    &mut result,
                    (b, j),
                    self.prime.mul(power, self.cartier.choose(r, j)),
                )?;
                power = self.prime.mul(power, scalar);
            }
        }
        Ok(result)
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        cover: &Cover,
        operation: Operation,
        input: Column,
    ) -> Result<Value> {
        let coordinates = match operation {
            Operation::Include => return Ok(Value::OddClass(input)),
            Operation::Frobenius => self.image(model, cover, Map::Frobenius, &input)?,
            Operation::Verschiebung => self.image(model, cover, Map::Verschiebung, &input)?,
            Operation::Translate(scalar) => {
                let scalar = scalar.rem_euclid(self.p() as i128) as u32;
                if scalar == 0 {
                    return Ok(Value::OddClass(input));
                }
                let Triple { u, v, overlap } = self.triple(model, cover, &input)?;
                let triple = Triple {
                    u: self.translate(model, &u, scalar)?,
                    v: self.translate(model, &v, scalar)?,
                    overlap: self.translate(model, &overlap, scalar)?,
                };
                self.decode(model, cover, triple)?
            }
        };
        odd_ext::persist(model, coordinates)
    }

    pub(crate) fn pair(
        &self,
        model: &mut Model,
        cover: &Cover,
        left: &[Coordinate],
        right: &[Coordinate],
    ) -> Result<Value> {
        let first = self.triple(model, cover, left)?;
        let second = self.triple(model, cover, right)?;
        Ok(Value::PrimeScalar(
            self.pair_triples(model, &first, &second)?,
        ))
    }

    pub(crate) fn gram(&self, model: &mut Model, cover: &Cover) -> Result<Vec<Vec<u32>>> {
        let dimension = self.dimension;
        let mut cells = dimension.saturating_mul(dimension.saturating_add(9));
        check_storage(cells)?;
        model.spend(dimension * dimension)?;
        let mut triples = Vec::with_capacity(dimension);
        for i in 0..dimension {
            let triple = self.triple(model, cover, &[(i, 1)])?;
            cells += 8 * (triple.u.len() + triple.v.len() + triple.overlap.len());
            check_storage(cells)?;
            triples.push(triple);
        }
        let mut matrix = vec![vec![0; dimension]; dimension];
        for i in 0..dimension {
            for j in i + 1..dimension {
                let value = self.pair_triples(model, &triples[i], &triples[j])?;
                matrix[i][j] = value;
                matrix[j][i] = self.prime.negate(value);
            }
        }
        Ok(matrix)
    }

    fn pair_triples(&self, model: &mut Model, first: &Triple, second: &Triple) -> Result<u32> {
        let a = self.residue(model, &first.u, &second.overlap)?;
        let b = self.residue(model, &first.overlap, &second.v)?;
        Ok(self.prime.add(a, self.prime.negate(b)))
    }

    fn residue(&self, model: &mut Model, left: &Laurent, right: &Laurent) -> Result<u32> {
        let mut result = 0;
        let top = self.p() - 1;
        let (source, target) = if left.len() <= right.len() {
            (left, right)
        } else {
            (right, left)
        };
        for (&(b, r), &coefficient) in source {
            model.spend(2 * (1 + target.len().max(1).ilog2() as usize))?;
            let x = -1 - b;
            let mut complement = target.get(&(x, top - r)).copied().unwrap_or(0);
            if r == top {
                // y^(2p-2) contributes y^(p-1) after reducing y^p.
                complement = self
                    .prime
                    .add(complement, target.get(&(x, top)).copied().unwrap_or(0));
            }
            result = self
                .prime
                .add(result, self.prime.mul(coefficient, complement));
        }
        Ok(result)
    }
}

fn add_term(
    model: &mut Model,
    prime: Prime,
    target: &mut Laurent,
    monomial: Monomial,
    coefficient: u32,
) -> Result<()> {
    model.spend(1 + target.len().max(1).ilog2() as usize)?;
    if coefficient == 0 {
        return Ok(());
    }
    match target.entry(monomial) {
        Entry::Vacant(entry) => {
            entry.insert(coefficient);
            check_storage(8 * target.len())?;
        }
        Entry::Occupied(mut entry) => {
            let sum = prime.add(*entry.get(), coefficient);
            if sum == 0 {
                entry.remove();
            } else {
                *entry.get_mut() = sum;
            }
        }
    }
    Ok(())
}

fn add_form(
    model: &mut Model,
    prime: Prime,
    target: &mut Laurent,
    source: &Laurent,
    scalar: u32,
) -> Result<()> {
    for (&term, &coefficient) in source {
        add_term(model, prime, target, term, prime.mul(scalar, coefficient))?;
    }
    Ok(())
}

fn dimension_error() -> Error {
    Error::plain("odd Artin–Schreier de Rham dimension is out of range")
}

fn basis_error() -> Error {
    Error::plain("odd Artin–Schreier Čech form is outside the expected differential basis")
}

struct FlagStep {
    space: Basis,
    image_rank: usize,
}

struct Flag {
    steps: BTreeMap<usize, FlagStep>,
    pending: VecDeque<usize>,
    cells: usize,
}

impl Flag {
    fn include(
        &mut self,
        model: &mut Model,
        prime: Prime,
        space: Basis,
        retained: usize,
    ) -> Result<()> {
        let rank = space.rank();
        if let Some(existing) = self.steps.get(&rank) {
            return contained(model, prime, &space, &existing.space);
        }
        self.cells += 6 + space.cells();
        check_storage(retained + self.cells)?;
        self.steps.insert(
            rank,
            FlagStep {
                space,
                image_rank: 0,
            },
        );
        self.pending.push_back(rank);
        Ok(())
    }
}

impl DeRham {
    fn flag_image(
        &self,
        model: &mut Model,
        cover: &Cover,
        input: &Basis,
        retained: usize,
    ) -> Result<Basis> {
        let mut result = Basis::default();
        for vector in input.vectors(self.prime) {
            let image = self.image(model, cover, Map::Verschiebung, &vector)?;
            result.insert(model, self.prime, image)?;
            check_storage(self.cells.get() + retained + result.cells())?;
        }
        Ok(result)
    }

    fn flag_preimage(
        &self,
        model: &mut Model,
        cover: &Cover,
        target: &Basis,
        retained: usize,
    ) -> Result<Basis> {
        let mut quotient = Basis::default();
        for vector in target.vectors(self.prime) {
            quotient.insert(model, self.prime, vector)?;
            check_storage(self.cells.get() + retained + quotient.cells())?;
        }
        let mut result = Basis::default();
        // Solve the F-column relations modulo target. Target generators
        // carry zero source witness, so dependencies give F^{-1}(target).
        for i in 0..self.dimension {
            let image = self.column(model, cover, Map::Frobenius, i)?;
            if let Some(vector) =
                quotient.insert_with_witness(model, self.prime, image.to_vec(), vec![(i, 1)])?
            {
                result.insert(model, self.prime, vector)?;
            }
            check_storage(self.cells.get() + retained + quotient.cells() + result.cells())?;
        }
        Ok(result)
    }

    pub(crate) fn final_type(&self, model: &mut Model, cover: &Cover) -> Result<Vec<usize>> {
        check_storage(self.dimension.saturating_mul(4))?;
        let mut flag = Flag {
            steps: BTreeMap::new(),
            pending: VecDeque::new(),
            cells: 0,
        };
        flag.include(model, self.prime, Basis::default(), self.cells.get())?;
        let mut ambient = Basis::default();
        for i in 0..self.dimension {
            ambient.insert(model, self.prime, vec![(i, 1)])?;
            check_storage(self.cells.get() + flag.cells + ambient.cells())?;
        }
        flag.include(model, self.prime, ambient, self.cells.get())?;
        while let Some(rank) = flag.pending.pop_front() {
            let step = flag.steps.get_mut(&rank).unwrap();
            let image = self.flag_image(model, cover, &step.space, flag.cells)?;
            let preimage =
                self.flag_preimage(model, cover, &step.space, flag.cells + image.cells())?;
            step.image_rank = image.rank();
            flag.include(
                model,
                self.prime,
                image,
                self.cells.get() + preimage.cells(),
            )?;
            flag.include(model, self.prime, preimage, self.cells.get())?;
        }
        let mut previous = &flag.steps[&0];
        let mut previous_rank = 0;
        let mut sequence = Vec::with_capacity(self.genus);
        for (&rank, step) in flag.steps.iter().skip(1) {
            contained(model, self.prime, &previous.space, &step.space)?;
            let width = rank - previous_rank;
            let rise = step.image_rank.checked_sub(previous.image_rank);
            if !matches!(rise, Some(0)) && rise != Some(width) {
                return Err(Error::plain(
                    "odd Artin–Schreier canonical EO flag has invalid Verschiebung ranks",
                ));
            }
            for dimension in previous_rank + 1..=rank.min(self.genus) {
                let increment = if rise == Some(width) {
                    dimension - previous_rank
                } else {
                    0
                };
                sequence.push(previous.image_rank + increment);
            }
            previous = step;
            previous_rank = rank;
        }
        Ok(sequence)
    }
}

fn contained(model: &mut Model, prime: Prime, lower: &Basis, upper: &Basis) -> Result<()> {
    for vector in lower.vectors(prime) {
        if !upper.reduce(model, prime, vector)?.is_empty() {
            return Err(Error::plain(
                "odd Artin–Schreier canonical EO subspaces are not nested",
            ));
        }
    }
    Ok(())
}
