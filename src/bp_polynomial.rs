use std::cmp::Ordering;
use std::collections::{btree_map::Entry, BTreeMap};

use crate::error::{Error, Result};
use crate::model::Model;

const MAX_POLYNOMIAL_CELLS: usize = 200_000;
type Monomial = Box<[(usize, usize)]>;

#[derive(Clone, Default)]
pub(crate) struct Polynomial {
    terms: BTreeMap<Monomial, i128>,
    cells: usize,
}

impl Polynomial {
    pub(crate) fn integer(value: i128) -> Self {
        if value == 0 {
            return Self::default();
        }
        Self {
            terms: [(Box::new([]) as Monomial, value)].into(),
            cells: 3,
        }
    }

    pub(crate) fn variable(index: usize) -> Self {
        Self {
            terms: [(Box::new([(index, 1usize)]) as Monomial, 1)].into(),
            cells: 5,
        }
    }

    pub(crate) fn cells(&self) -> usize {
        self.cells
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(self.cells + rhs.cells)?;
        Ok(self.terms.cmp(&rhs.terms))
    }

    fn insert(
        &mut self,
        model: &mut Model,
        monomial: Monomial,
        coefficient: i128,
        subtract: bool,
    ) -> Result<()> {
        if coefficient == 0 {
            return Ok(());
        }
        let cells = 3 + 2 * monomial.len();
        model.spend(cells * (1 + self.terms.len().max(1).ilog2() as usize))?;
        match self.terms.entry(monomial) {
            Entry::Occupied(mut entry) => {
                let next = if subtract {
                    entry.get().checked_sub(coefficient)
                } else {
                    entry.get().checked_add(coefficient)
                }
                .ok_or_else(overflow)?;
                if next == 0 {
                    entry.remove();
                    self.cells -= cells;
                } else {
                    *entry.get_mut() = next;
                }
            }
            Entry::Vacant(entry) => {
                let coefficient = if subtract {
                    coefficient.checked_neg().ok_or_else(overflow)?
                } else {
                    coefficient
                };
                if cells > MAX_POLYNOMIAL_CELLS - self.cells {
                    return Err(Error::plain("BP polynomial exceeds the term storage limit"));
                }
                self.cells += cells;
                entry.insert(coefficient);
            }
        }
        Ok(())
    }

    pub(crate) fn add_assign(
        &mut self,
        model: &mut Model,
        rhs: &Self,
        subtract: bool,
    ) -> Result<()> {
        model.spend(rhs.cells)?;
        for (monomial, &coefficient) in &rhs.terms {
            self.insert(model, monomial.clone(), coefficient, subtract)?;
        }
        Ok(())
    }

    pub(crate) fn scale(&mut self, model: &mut Model, scalar: i128) -> Result<()> {
        model.spend(self.cells)?;
        if scalar == 0 {
            *self = Self::default();
        } else {
            for coefficient in self.terms.values_mut() {
                *coefficient = coefficient.checked_mul(scalar).ok_or_else(overflow)?;
            }
        }
        Ok(())
    }

    pub(crate) fn divisible(&self, model: &mut Model, positive: i128) -> Result<bool> {
        model.spend(self.terms.len())?;
        Ok(self
            .terms
            .values()
            .all(|coefficient| coefficient % positive == 0))
    }

    pub(crate) fn divide_exact(&mut self, model: &mut Model, positive: i128) -> Result<()> {
        if !self.divisible(model, positive)? {
            return Err(Error::plain("BP polynomial is not exactly divisible"));
        }
        model.spend(self.terms.len())?;
        for coefficient in self.terms.values_mut() {
            *coefficient /= positive;
        }
        Ok(())
    }

    pub(crate) fn multiply(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let mut result = Self::default();
        for (left, &a) in &self.terms {
            for (right, &b) in &rhs.terms {
                model.spend(1 + left.len() + right.len())?;
                let monomial = multiply_monomials(left, right)?;
                result.insert(
                    model,
                    monomial,
                    a.checked_mul(b).ok_or_else(overflow)?,
                    false,
                )?;
            }
        }
        Ok(result)
    }

    pub(crate) fn power(&self, model: &mut Model, mut exponent: usize) -> Result<Self> {
        model.spend(self.cells)?;
        if exponent == 0 {
            return Ok(Self::integer(1));
        }
        let mut result = Self::integer(1);
        let mut base = self.clone();
        loop {
            model.spend(1)?;
            if exponent & 1 != 0 {
                result = result.multiply(model, &base)?;
            }
            exponent >>= 1;
            if exponent == 0 {
                return Ok(result);
            }
            base = base.multiply(model, &base)?;
        }
    }

    pub(crate) fn max_used_index(
        &self,
        model: &mut Model,
        width: usize,
        slot: usize,
    ) -> Result<usize> {
        model.spend(self.cells)?;
        let start = slot * width;
        Ok(self
            .terms
            .keys()
            .flat_map(|monomial| monomial.iter())
            .filter_map(|&(variable, _)| variable.checked_sub(start).filter(|&index| index < width))
            .max()
            .map_or(0, |index| index + 1))
    }

    pub(crate) fn substitute(&self, model: &mut Model, images: &[Self]) -> Result<Self> {
        let mut result = Self::default();
        for (monomial, &coefficient) in &self.terms {
            let mut term = Self::integer(coefficient);
            for &(variable, exponent) in monomial.iter() {
                let power = images[variable].power(model, exponent)?;
                term = term.multiply(model, &power)?;
            }
            result.add_assign(model, &term, false)?;
        }
        Ok(result)
    }

    pub(crate) fn write(
        &self,
        mut out: impl FnMut(&str) -> Result<()>,
        width: usize,
    ) -> Result<()> {
        if self.is_zero() {
            return out("0");
        }
        for (index, (monomial, &coefficient)) in self.terms.iter().enumerate() {
            if index == 0 {
                if coefficient < 0 {
                    out("-")?;
                }
            } else {
                out(if coefficient < 0 { " - " } else { " + " })?;
            }
            let magnitude = coefficient.unsigned_abs();
            if magnitude != 1 || monomial.is_empty() {
                out(&magnitude.to_string())?;
                if !monomial.is_empty() {
                    out("*")?;
                }
            }
            for (i, &(variable, exponent)) in monomial.iter().enumerate() {
                if i != 0 {
                    out("*")?;
                }
                let height = variable % width + 1;
                let slot = variable / width;
                if slot == 0 {
                    out(&format!("v{height}"))?;
                } else {
                    out(&format!("t{slot}_{height}"))?;
                }
                if exponent != 1 {
                    out(&format!("^{exponent}"))?;
                }
            }
        }
        Ok(())
    }
}

fn multiply_monomials(
    mut left: &[(usize, usize)],
    mut right: &[(usize, usize)],
) -> Result<Monomial> {
    let mut result = Vec::with_capacity(left.len() + right.len());
    while let (Some((&(i, a), next_left)), Some((&(j, b), next_right))) =
        (left.split_first(), right.split_first())
    {
        match i.cmp(&j) {
            Ordering::Less => {
                result.push((i, a));
                left = next_left;
            }
            Ordering::Equal => {
                result.push((
                    i,
                    a.checked_add(b)
                        .ok_or_else(|| Error::plain("BP monomial exponent overflow"))?,
                ));
                left = next_left;
                right = next_right;
            }
            Ordering::Greater => {
                result.push((j, b));
                right = next_right;
            }
        }
    }
    result.extend_from_slice(left);
    result.extend_from_slice(right);
    Ok(result.into_boxed_slice())
}

fn overflow() -> Error {
    Error::plain("BP coefficient exceeds the signed i128 range")
}
