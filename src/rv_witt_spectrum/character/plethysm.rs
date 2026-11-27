use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;

use super::super::norm::power;
use super::super::{add, multiply, zero, Integer};
use super::table::{
    accumulate, divide, factorial, from_power_basis, power_basis, CharacterOracle, Partition,
    MAX_WEIGHT,
};
use super::{output_weight, overflow, Character};

impl Character {
    pub(crate) fn plethysm(&self, model: &mut Model, inner: &Self) -> Result<Self> {
        let weight = output_weight(self.weight, inner.weight)?;
        if self.weight == 0 {
            return self.scale(model, 1);
        }
        if self.terms.is_empty() || inner.terms.is_empty() {
            return Ok(Self::zero(weight));
        }
        let first_inner = *inner.terms.values().next().unwrap();
        if self.weight == 1 {
            return inner.scale(model, *self.terms.values().next().unwrap());
        }
        if inner.weight == 1 && first_inner == 1 {
            return self.scale(model, 1);
        }
        if inner.weight == 0 {
            return self.at_constant(model, first_inner);
        }
        let n = bounded(weight)?;
        let (a, b) = (self.weight as u8, inner.weight as u8);
        let mut oracle = CharacterOracle::default();
        let inner_values = inner.class_values(model, &mut oracle, b)?;
        let inner_terms = power_basis(model, b, inner_values)?;
        let outer_values = self.class_values(model, &mut oracle, a)?;
        let outer_terms = power_basis(model, a, outer_values)?;
        let inner_factorial = factorial(b);
        // This is the wreath product order, so it divides (ab)!.
        let denominator = factorial(a) * inner_factorial.pow(u32::from(a));
        let mut power_basis = BTreeMap::new();
        for (cycles, coefficient) in outer_terms {
            let scale = inner_factorial.pow(u32::from(a) - cycles.len() as u32);
            let coefficient = multiply(model, coefficient, Integer::from_i128(scale as i128))?;
            let mut products = BTreeMap::from([(Vec::new(), coefficient)]);
            for r in cycles {
                let mut next = BTreeMap::new();
                for (left, &alpha) in &products {
                    for (right, &beta) in &inner_terms {
                        model.spend(1 + left.len() + right.len())?;
                        let partition = merge_cycles(left, right, r);
                        let product = multiply(model, alpha, beta)?;
                        accumulate(&mut next, partition, product)?;
                    }
                }
                products = next;
            }
            model.spend(products.len())?;
            for (partition, value) in products {
                accumulate(&mut power_basis, partition, value)?;
            }
        }
        from_power_basis(model, n, denominator, &power_basis, &mut oracle)
    }

    fn at_constant(&self, model: &mut Model, value: i128) -> Result<Self> {
        if value == 1 || value == -1 {
            let coefficient = self
                .terms
                .iter()
                .find_map(|(shape, &coefficient)| {
                    let rows = shape.rows();
                    if (value == 1 && rows.len() == 1) || (value == -1 && rows[0] == 1) {
                        Some(coefficient)
                    } else {
                        None
                    }
                })
                .unwrap_or(0);
            model.spend(self.terms.len())?;
            let coefficient = if value == -1 && !self.weight.is_multiple_of(2) {
                coefficient.checked_neg().ok_or_else(overflow)?
            } else {
                coefficient
            };
            return Self::constant(coefficient);
        }
        let n = bounded(self.weight)?;
        let mut oracle = CharacterOracle::default();
        let values = self.class_values(model, &mut oracle, n)?;
        let terms = power_basis(model, n, values)?;
        let value = Integer::from_i128(value);
        let mut sum = zero();
        for (cycles, coefficient) in terms {
            let substitute = power(model, value, cycles.len() as u128)?;
            sum = add(sum, multiply(model, coefficient, substitute)?)?;
        }
        Self::constant(divide(model, sum, factorial(n))?)
    }
}

fn merge_cycles(left: &[u8], right: &[u8], factor: u8) -> Partition {
    let mut merged = Vec::with_capacity(left.len() + right.len());
    let mut left = left.iter().copied().peekable();
    let mut right = right.iter().map(|&part| part * factor).peekable();
    while let (Some(&a), Some(&b)) = (left.peek(), right.peek()) {
        merged.push(if a >= b {
            left.next().unwrap()
        } else {
            right.next().unwrap()
        });
    }
    merged.extend(left);
    merged.extend(right);
    merged
}

fn bounded(n: u128) -> Result<u8> {
    if n > u128::from(MAX_WEIGHT) {
        return Err(Error::plain(
            "Schur plethysm exceeds the weight limit of 16 boxes",
        ));
    }
    Ok(n as u8)
}
