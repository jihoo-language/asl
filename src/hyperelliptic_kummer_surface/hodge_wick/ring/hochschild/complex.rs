use crate::error::Error;

use super::{
    basis_error, linear, storage, work_limit, Algebra, Cochain, Grade, Key, Model, Rational, Result,
};

#[derive(Default)]
struct Axis(Vec<Key>);

impl Axis {
    fn new(model: &mut Model, algebra: &Algebra, grade: &Grade) -> Result<Self> {
        if !grade.possible() {
            return Ok(Self::default());
        }
        let arity = usize::try_from(grade.arity).map_err(|_| work_limit())?;
        storage(3usize.saturating_mul(arity.saturating_add(1)))?;
        let mut word = vec![0; arity];
        let mut prefix = vec![0i128; arity + 1];
        let mut depth = 0;
        let mut cells = 0;
        let mut keys = Vec::new();
        loop {
            if depth == arity {
                let weight = (prefix[depth] + grade.shift) as usize;
                let outputs = algebra.weight_range(weight);
                model.spend(outputs.len().saturating_mul(arity + 1))?;
                for output in outputs {
                    cells += arity + 8;
                    storage(cells)?;
                    keys.push(Key {
                        word: word.clone(),
                        output,
                    });
                }
                if depth == 0 {
                    break;
                }
                depth -= 1;
                continue;
            }
            let remaining = (arity - depth - 1) as i128;
            let top = algebra.top() as i128;
            let mut found = false;
            for i in word[depth] + 1..algebra.len() {
                model.spend(1)?;
                let weight = prefix[depth] + algebra.weights[i] as i128 + grade.shift;
                if weight + remaining * top < 0 {
                    continue;
                }
                if weight + remaining > top {
                    break;
                }
                prefix[depth + 1] = prefix[depth] + algebra.weights[i] as i128;
                word[depth] = i;
                depth += 1;
                if depth < arity {
                    word[depth] = 0;
                }
                found = true;
                break;
            }
            if !found {
                if depth == 0 {
                    break;
                }
                depth -= 1;
            }
        }
        Ok(Self(keys))
    }

    fn locate(&self, key: &Key) -> Result<usize> {
        self.0
            .binary_search(key)
            .map_err(|_| Error::plain("Hochschild coefficient lies outside the homogeneous slice"))
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn coordinates(&self, value: &Cochain) -> Result<Vec<Rational>> {
        let mut row = vec![Rational::zero(); self.len()];
        for (key, &coefficient) in &value.entries {
            row[self.locate(key)?] = coefficient;
        }
        Ok(row)
    }

    fn cochain(&self, model: &mut Model, row: &[Rational]) -> Result<Cochain> {
        let mut result = Cochain::default();
        for (key, &coefficient) in self.0.iter().zip(row) {
            if !coefficient.is_zero() {
                result.insert(model, key.clone(), coefficient)?;
            }
        }
        Ok(result)
    }
}

pub(super) fn differential_terms(
    model: &mut Model,
    algebra: &Algebra,
    key: &Key,
    mut emit: impl FnMut(&mut Model, Key, Rational) -> Result<()>,
) -> Result<()> {
    let arity = key.word.len();
    if arity == 0 {
        // R is commutative.
        return Ok(());
    }
    model.spend(algebra.len() + arity)?;
    for x in 1..algebra.len() {
        for &(output, coefficient) in algebra.product(x, key.output) {
            let mut left = Vec::with_capacity(arity + 1);
            left.push(x);
            left.extend_from_slice(&key.word);
            emit(model, Key { word: left, output }, coefficient)?;
            let mut right = key.word.clone();
            right.push(x);
            let coefficient = if arity.is_multiple_of(2) {
                coefficient.negate()
            } else {
                coefficient
            };
            emit(
                model,
                Key {
                    word: right,
                    output,
                },
                coefficient,
            )?;
        }
    }
    for (slot, &input) in key.word.iter().enumerate() {
        for &(a, b, coefficient) in algebra.splits(input) {
            let mut word = Vec::with_capacity(arity + 1);
            word.extend_from_slice(&key.word[..slot]);
            word.extend([a, b]);
            word.extend_from_slice(&key.word[slot + 1..]);
            let coefficient = if slot.is_multiple_of(2) {
                coefficient.negate()
            } else {
                coefficient
            };
            emit(
                model,
                Key {
                    word,
                    output: key.output,
                },
                coefficient,
            )?;
        }
    }
    Ok(())
}

pub(super) struct Complex {
    axis: Axis,
    equations: linear::Echelon,
    boundaries: linear::Echelon,
}

impl Complex {
    pub(super) fn new(model: &mut Model, algebra: &Algebra, grade: &Grade) -> Result<Self> {
        let axis = Axis::new(model, algebra, grade)?;
        if axis.len() == 0 {
            let equations = linear::rref(model, Vec::new(), 0)?;
            let boundaries = linear::rref(model, Vec::new(), 0)?;
            return Ok(Self {
                axis,
                equations,
                boundaries,
            });
        }
        let next = Axis::new(model, algebra, &grade.with_arity_offset(1)?)?;
        let previous = Axis::new(model, algebra, &grade.with_arity_offset(-1)?)?;
        let mut differential = linear::matrix(next.len(), axis.len())?;
        for (column, key) in axis.0.iter().enumerate() {
            differential_terms(model, algebra, key, |model, term, coefficient| {
                let coordinate = &mut differential[next.locate(&term)?][column];
                *coordinate = coordinate.add(model, coefficient)?;
                Ok(())
            })?;
        }
        let equations = linear::rref(model, differential, axis.len())?;
        let mut boundaries = linear::matrix(previous.len(), axis.len())?;
        for (row, key) in boundaries.iter_mut().zip(previous.0) {
            differential_terms(model, algebra, &key, |model, term, coefficient| {
                let coordinate = &mut row[axis.locate(&term)?];
                *coordinate = coordinate.add(model, coefficient)?;
                Ok(())
            })?;
        }
        let boundaries = linear::rref(model, boundaries, axis.len())?;
        for row in &boundaries.rows {
            if !annihilates(model, &equations, row)? {
                return Err(Error::plain(
                    "Hochschild differential failed the square-zero check",
                ));
            }
        }
        Ok(Self {
            axis,
            equations,
            boundaries,
        })
    }

    /// Choose the preimage with the free coordinates of the linear system zero.
    pub(super) fn primitive(
        &self,
        model: &mut Model,
        algebra: &Algebra,
        grade: &Grade,
        boundary: &Cochain,
    ) -> Result<Cochain> {
        let previous = Axis::new(model, algebra, &grade.with_arity_offset(-1)?)?;
        let mut differential = linear::matrix(self.axis.len(), previous.len())?;
        for (column, key) in previous.0.iter().enumerate() {
            differential_terms(model, algebra, key, |model, term, coefficient| {
                let coordinate = &mut differential[self.axis.locate(&term)?][column];
                *coordinate = coordinate.add(model, coefficient)?;
                Ok(())
            })?;
        }
        let rhs = self.axis.coordinates(boundary)?;
        let solution = linear::solve(model, differential, previous.len(), &rhs)?
            .ok_or_else(|| Error::plain("Hochschild boundary has no exact primitive"))?;
        previous.cochain(model, &solution)
    }

    pub(super) fn dimensions(&self) -> [usize; 4] {
        let chains = self.axis.len();
        let cycles = chains - self.equations.rows.len();
        let boundaries = self.boundaries.rows.len();
        [chains, cycles, boundaries, cycles - boundaries]
    }

    fn reduce(&self, model: &mut Model, coordinates: &mut [Rational]) -> Result<()> {
        for (&pivot, row) in self.boundaries.pivots.iter().zip(&self.boundaries.rows) {
            linear::add_multiple(model, coordinates, row, coordinates[pivot].negate())?;
        }
        Ok(())
    }

    pub(super) fn project(&self, model: &mut Model, value: &Cochain) -> Result<Cochain> {
        let mut coordinates = self.axis.coordinates(value)?;
        if !annihilates(model, &self.equations, &coordinates)? {
            return Err(Error::plain(
                "normalized Hochschild cochain is not a cocycle",
            ));
        }
        self.reduce(model, &mut coordinates)?;
        self.axis.cochain(model, &coordinates)
    }

    pub(super) fn basis(&self, model: &mut Model, index: usize) -> Result<Cochain> {
        if index >= self.dimensions()[3] {
            return Err(basis_error());
        }
        let mut cycles = self.equations.kernel(self.axis.len())?;
        for row in &mut cycles {
            self.reduce(model, row)?;
        }
        let representative = linear::rref(model, cycles, self.axis.len())?;
        if representative.rows.len() != self.dimensions()[3] {
            return Err(Error::plain(
                "Hochschild cycle quotient has unexpected rank",
            ));
        }
        self.axis.cochain(model, &representative.rows[index])
    }
}

fn annihilates(
    model: &mut Model,
    equations: &linear::Echelon,
    coordinates: &[Rational],
) -> Result<bool> {
    for equation in &equations.rows {
        let mut test = Rational::zero();
        for (&form, &value) in equation.iter().zip(coordinates) {
            let product = form.multiply(model, value)?;
            test = test.add(model, product)?;
        }
        if !test.is_zero() {
            return Ok(false);
        }
    }
    Ok(true)
}
