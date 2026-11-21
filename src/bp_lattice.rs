use crate::bp_smith::Matrix;
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear;

// Every change replaces B_i by B_i - c_i B_k (i != k) and B_k by p B_k.
// This retains a basis of the exact kernel of the next digit congruence.
struct Change {
    pivot: usize,
    coefficients: Box<[(usize, u32)]>,
}

pub(crate) struct Lattice {
    scales: Vec<i128>,
    changes: Vec<Change>,
    prime: Prime,
    cells: usize,
}

impl Lattice {
    pub(crate) fn cells(&self) -> usize {
        self.cells
    }

    pub(crate) fn build(
        model: &mut Model,
        prime: Prime,
        weights: &[i128],
        filtration: i128,
        common: i128,
        constraint: Option<(&Matrix, &[i128], i128)>,
    ) -> Result<Self> {
        model.spend(weights.len())?;
        let scales = weights
            .iter()
            .map(|&weight| power(prime, (filtration - weight).max(0) - common))
            .collect::<Result<Vec<_>>>()?;
        let cells = 2 * scales.len();
        let mut result = Self {
            scales,
            changes: Vec::new(),
            prime,
            cells,
        };
        if let Some((differential, target_weights, threshold)) = constraint {
            for (index, &weight) in target_weights.iter().enumerate() {
                let order = threshold - weight - common;
                if order > 0 {
                    result.impose(model, differential.row(index), order)?;
                }
            }
        }
        Ok(result)
    }

    fn impose(&mut self, model: &mut Model, row: &[i128], order: i128) -> Result<()> {
        model.spend(row.len())?;
        if row.iter().all(|&coefficient| coefficient == 0) {
            return Ok(());
        }
        let mut row = row
            .iter()
            .zip(&self.scales)
            .map(|(&entry, &scale)| multiply(entry, scale))
            .collect::<Result<Vec<_>>>()?;
        let prime = i128::from(self.prime.get());
        for change in &self.changes {
            Self::apply_to_row(model, prime, &mut row, change)?;
        }
        let mut power = 1;
        for digit in 0..order {
            model.spend(row.len())?;
            debug_assert!(row.iter().all(|value| value % power == 0));
            let residue = |value: i128| (value / power).rem_euclid(prime) as u32;
            if let Some((pivot, leading)) = row
                .iter()
                .map(|&value| residue(value))
                .enumerate()
                .find(|&(_, value)| value != 0)
            {
                let inverse = prime_linear::invert(model, self.prime, leading)?;
                model.spend(row.len())?;
                let coefficients = row
                    .iter()
                    .enumerate()
                    .filter_map(|(index, &value)| {
                        let coefficient = self.prime.mul(residue(value), inverse);
                        (index != pivot && coefficient != 0).then_some((index, coefficient))
                    })
                    .collect::<Box<_>>();
                let change = Change {
                    pivot,
                    coefficients,
                };
                self.cells += 5 + 2 * change.coefficients.len();
                if self.cells > MAX_PIVOT_CELLS {
                    return Err(Error::plain("BP lattice basis exceeds the storage limit"));
                }
                Self::apply_to_row(model, prime, &mut row, &change)?;
                self.changes.push(change);
            }
            if digit + 1 < order {
                power = multiply(power, prime)?;
            }
        }
        Ok(())
    }

    fn apply_to_row(
        model: &mut Model,
        prime: i128,
        row: &mut [i128],
        change: &Change,
    ) -> Result<()> {
        model.spend(1 + change.coefficients.len())?;
        let pivot = row[change.pivot];
        for &(index, coefficient) in &change.coefficients {
            row[index] = subtract(row[index], multiply(pivot, i128::from(coefficient))?)?;
        }
        row[change.pivot] = multiply(pivot, prime)?;
        Ok(())
    }

    pub(crate) fn encode(
        &self,
        model: &mut Model,
        mut values: Vec<i128>,
    ) -> Result<Option<Vec<i128>>> {
        model.spend(values.len())?;
        for (value, &scale) in values.iter_mut().zip(&self.scales) {
            if *value % scale != 0 {
                return Ok(None);
            }
            *value /= scale;
        }
        let prime = i128::from(self.prime.get());
        for change in &self.changes {
            model.spend(1 + change.coefficients.len())?;
            let mut pivot = values[change.pivot];
            for &(index, coefficient) in &change.coefficients {
                pivot = add(pivot, multiply(values[index], i128::from(coefficient))?)?;
            }
            if pivot % prime != 0 {
                return Ok(None);
            }
            values[change.pivot] = pivot / prime;
        }
        Ok(Some(values))
    }

    pub(crate) fn generator(&self, model: &mut Model, index: usize) -> Result<Vec<i128>> {
        let mut values = vec![0; self.scales.len()];
        model.spend(values.len())?;
        values[index] = 1;
        self.decode(model, values)
    }

    pub(crate) fn decode(&self, model: &mut Model, mut values: Vec<i128>) -> Result<Vec<i128>> {
        let prime = i128::from(self.prime.get());
        for change in self.changes.iter().rev() {
            model.spend(1 + change.coefficients.len())?;
            let mut pivot = multiply(values[change.pivot], prime)?;
            for &(index, coefficient) in &change.coefficients {
                pivot = subtract(pivot, multiply(values[index], i128::from(coefficient))?)?;
            }
            values[change.pivot] = pivot;
        }
        model.spend(values.len())?;
        for (value, &scale) in values.iter_mut().zip(&self.scales) {
            *value = multiply(*value, scale)?;
        }
        Ok(values)
    }
}

pub(crate) fn power(prime: Prime, exponent: i128) -> Result<i128> {
    i128::from(prime.get())
        .checked_pow(u32::try_from(exponent).map_err(|_| overflow())?)
        .ok_or_else(overflow)
}

fn add(left: i128, right: i128) -> Result<i128> {
    left.checked_add(right).ok_or_else(overflow)
}

fn subtract(left: i128, right: i128) -> Result<i128> {
    left.checked_sub(right).ok_or_else(overflow)
}

fn multiply(left: i128, right: i128) -> Result<i128> {
    left.checked_mul(right).ok_or_else(overflow)
}

fn overflow() -> Error {
    Error::plain("BP filtered lattice arithmetic exceeds the signed i128 range")
}
