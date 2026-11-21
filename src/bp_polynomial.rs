use std::cmp::Ordering;
use std::collections::{btree_map::Entry, BTreeMap};

use crate::error::{Error, Result};
use crate::model::Model;

const MAX_POLYNOMIAL_CELLS: usize = 200_000;
pub(crate) type Monomial = Box<[(usize, usize)]>;

#[derive(Clone, Copy)]
pub(crate) struct Quotient {
    pub(crate) prime: u32,
    pub(crate) killed_variables: usize,
    pub(crate) cutoff: Option<(usize, usize)>,
}

impl Quotient {
    fn retains(self, monomial: &Monomial) -> bool {
        monomial
            .first()
            .is_none_or(|&(index, _)| index >= self.killed_variables)
            && self.cutoff.is_none_or(|(variable, power)| {
                monomial
                    .binary_search_by_key(&variable, |&(index, _)| index)
                    .ok()
                    .is_none_or(|index| monomial[index].1 < power)
            })
    }
}

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

    pub(crate) fn from_terms<'a>(
        model: &mut Model,
        terms: impl IntoIterator<Item = (&'a Monomial, i128)>,
    ) -> Result<Self> {
        let mut result = Self::default();
        for (monomial, coefficient) in terms {
            result.insert(model, monomial.clone(), coefficient, false)?;
        }
        Ok(result)
    }

    pub(crate) fn terms(&self) -> impl Iterator<Item = (&Monomial, i128)> {
        self.terms
            .iter()
            .map(|(monomial, &coefficient)| (monomial, coefficient))
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

    pub(crate) fn divide_exact(&mut self, model: &mut Model, positive: i128) -> Result<()> {
        debug_assert!(positive > 0);
        if !self.divide_if_exact(model, positive as u128)? {
            return Err(Error::plain("BP polynomial is not exactly divisible"));
        }
        Ok(())
    }

    pub(crate) fn divide_if_exact(&mut self, model: &mut Model, divisor: u128) -> Result<bool> {
        debug_assert!(divisor != 0);
        if divisor == 1 {
            return Ok(true);
        }
        model.spend(self.terms.len())?;
        if self
            .terms
            .values()
            .any(|value| value.unsigned_abs() % divisor != 0)
        {
            return Ok(false);
        }
        model.spend(self.terms.len())?;
        for value in self.terms.values_mut() {
            // Dividing by at least two also handles the magnitude of i128::MIN.
            let magnitude = value.unsigned_abs() / divisor;
            *value = value.signum() * magnitude as i128;
        }
        Ok(true)
    }

    pub(crate) fn reduce(&mut self, model: &mut Model, quotient: Quotient) -> Result<()> {
        model.spend(self.cells)?;
        let modulus = i128::from(quotient.prime);
        self.terms.retain(|monomial, coefficient| {
            *coefficient = coefficient.rem_euclid(modulus);
            let keep = *coefficient != 0 && quotient.retains(monomial);
            if !keep {
                self.cells -= 3 + 2 * monomial.len();
            }
            keep
        });
        Ok(())
    }

    pub(crate) fn divide_variable_power(
        self,
        model: &mut Model,
        variable: usize,
        power: usize,
    ) -> Result<Option<Self>> {
        let mut quotient = Self::default();
        for (monomial, coefficient) in self.terms {
            model.spend(1 + monomial.len())?;
            let Ok(index) = monomial.binary_search_by_key(&variable, |&(index, _)| index) else {
                return Ok(None);
            };
            let Some(remaining) = monomial[index].1.checked_sub(power) else {
                return Ok(None);
            };
            let mut monomial = monomial.into_vec();
            if remaining == 0 {
                monomial.remove(index);
            } else {
                monomial[index].1 = remaining;
            }
            quotient.insert(model, monomial.into_boxed_slice(), coefficient, false)?;
        }
        Ok(Some(quotient))
    }

    pub(crate) fn multiply(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        self.multiply_in(model, rhs, None)
    }

    pub(crate) fn multiply_in(
        &self,
        model: &mut Model,
        rhs: &Self,
        quotient: Option<Quotient>,
    ) -> Result<Self> {
        let cutoff = quotient.and_then(|q| q.cutoff);
        let mut result = Self::default();
        for (left, &a) in &self.terms {
            for (right, &b) in &rhs.terms {
                model.spend(1 + left.len() + right.len())?;
                let Some(monomial) = multiply_monomials(left, right, cutoff)? else {
                    continue;
                };
                result.insert(
                    model,
                    monomial,
                    a.checked_mul(b).ok_or_else(overflow)?,
                    false,
                )?;
            }
        }
        if let Some(quotient) = quotient {
            result.reduce(model, quotient)?;
        }
        Ok(result)
    }

    pub(crate) fn power(&self, model: &mut Model, exponent: usize) -> Result<Self> {
        self.power_in(model, exponent, None)
    }

    pub(crate) fn power_in(
        &self,
        model: &mut Model,
        mut exponent: usize,
        quotient: Option<Quotient>,
    ) -> Result<Self> {
        model.spend(self.cells)?;
        if exponent == 0 {
            return Ok(Self::integer(1));
        }
        let mut result = Self::integer(1);
        let mut base = self.clone();
        loop {
            model.spend(1)?;
            if exponent & 1 != 0 {
                result = result.multiply_in(model, &base, quotient)?;
            }
            exponent >>= 1;
            if exponent == 0 {
                return Ok(result);
            }
            base = base.multiply_in(model, &base, quotient)?;
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

    pub(crate) fn substitute(
        &self,
        model: &mut Model,
        images: &[Self],
        quotient: Option<Quotient>,
    ) -> Result<Self> {
        let mut result = Self::default();
        for (monomial, &coefficient) in &self.terms {
            let mut term = Self::integer(coefficient);
            for &(variable, exponent) in monomial.iter() {
                let power = images[variable].power_in(model, exponent, quotient)?;
                term = term.multiply_in(model, &power, quotient)?;
            }
            result.add_assign(model, &term, false)?;
        }
        Ok(result)
    }

    // The surviving variable indices must be mapped in strictly increasing order.
    pub(crate) fn reindex(
        &self,
        model: &mut Model,
        image: impl Fn(usize) -> Option<usize>,
    ) -> Result<Self> {
        model.spend(self.cells)?;
        let mut output = Self::default();
        for (monomial, &coefficient) in &self.terms {
            let mapped: Option<Monomial> = monomial
                .iter()
                .map(|&(index, power)| image(index).map(|index| (index, power)))
                .collect();
            if let Some(mapped) = mapped {
                output.insert(model, mapped, coefficient, false)?;
            }
        }
        Ok(output)
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
    cutoff: Option<(usize, usize)>,
) -> Result<Option<Monomial>> {
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
                if let Some((variable, power)) = cutoff {
                    if i == variable && (a >= power || b >= power - a) {
                        return Ok(None);
                    }
                }
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
    Ok(Some(result.into_boxed_slice()))
}

fn overflow() -> Error {
    Error::plain("BP coefficient exceeds the signed i128 range")
}
