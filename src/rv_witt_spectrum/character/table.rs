use std::collections::{BTreeMap, HashMap};

use crate::error::{Error, Result};
use crate::model::Model;
use crate::wide_integer::Natural;

use super::super::{add, exact_divide, multiply, zero, Integer, SchurShape};
use super::{overflow, Character};

pub(super) const MAX_WEIGHT: u8 = 16;
const MEMO_LIMIT: usize = 200_000;
pub(super) type Partition = Vec<u8>;
pub(super) type Coefficients = BTreeMap<Partition, Integer>;

impl Character {
    pub(super) fn class_values(
        &self,
        model: &mut Model,
        oracle: &mut CharacterOracle,
        n: u8,
    ) -> Result<Coefficients> {
        let input: Vec<_> = self
            .terms
            .iter()
            .map(|(shape, &value)| {
                (
                    beta_numbers(shape.rows().iter().map(|&i| i as u8)),
                    Integer::from_i128(value),
                )
            })
            .collect();
        let mut terms = BTreeMap::new();
        for partition in partitions(model, n)? {
            let code = pack(&partition);
            let mut value = zero();
            for &(beta, coefficient) in &input {
                let irreducible = Integer::from_i128(oracle.value(model, beta, code)?);
                value = add(value, multiply(model, coefficient, irreducible)?)?;
            }
            if !value.is_zero() {
                terms.insert(partition, value);
            }
        }
        Ok(terms)
    }
}

// Store the Frobenius characteristic in the power basis with common
// denominator n! rather than one denominator z_ρ for each class.
pub(super) fn power_basis(
    model: &mut Model,
    n: u8,
    mut values: Coefficients,
) -> Result<Coefficients> {
    let factorial = factorial(n);
    for (cycle_type, value) in &mut values {
        let scale = Integer::from_i128((factorial / centralizer(cycle_type)) as i128);
        *value = multiply(model, *value, scale)?;
    }
    Ok(values)
}

pub(super) fn from_power_basis(
    model: &mut Model,
    n: u8,
    denominator: u128,
    power_basis: &Coefficients,
    oracle: &mut CharacterOracle,
) -> Result<Character> {
    let packed: Vec<_> = power_basis
        .iter()
        .map(|(cycle_type, &value)| (pack(cycle_type), value))
        .collect();
    let mut result = Character::zero(u128::from(n));
    for rows in partitions(model, n)? {
        let beta = beta_numbers(rows.iter().copied());
        let mut sum = zero();
        for &(code, power) in &packed {
            let irreducible = Integer::from_i128(oracle.value(model, beta, code)?);
            sum = add(sum, multiply(model, power, irreducible)?)?;
        }
        let coefficient = divide(model, sum, denominator)?;
        if coefficient != 0 {
            let shape = SchurShape::new(rows.into_iter().map(u128::from).collect())?;
            result.terms.insert(shape, coefficient);
        }
    }
    Ok(result)
}

pub(super) fn divide(model: &mut Model, value: Integer, denominator: u128) -> Result<i128> {
    exact_divide(model, value, Natural::from_small(denominator))?
        .signed()
        .ok_or_else(overflow)
}

pub(super) fn accumulate(
    terms: &mut Coefficients,
    partition: Partition,
    value: Integer,
) -> Result<()> {
    if !value.is_zero() {
        let sum = add(terms.get(&partition).copied().unwrap_or_else(zero), value)?;
        if sum.is_zero() {
            terms.remove(&partition);
        } else {
            terms.insert(partition, sum);
        }
    }
    Ok(())
}

#[derive(Default)]
pub(super) struct CharacterOracle(HashMap<(u64, u128), i128>);

impl CharacterOracle {
    pub(super) fn value(&mut self, model: &mut Model, beads: u64, classes: u128) -> Result<i128> {
        model.spend(1)?;
        if classes == 0 {
            return Ok(i128::from(beads == (1u64 << MAX_WEIGHT) - 1));
        }
        if let Some(&cached) = self.0.get(&(beads, classes)) {
            return Ok(cached);
        }
        let length = (classes & 31) as u32;
        let mut result = 0i128;
        let mut sources = beads >> length;
        while sources != 0 {
            model.spend(1)?;
            let target = sources.trailing_zeros();
            sources &= sources - 1;
            let target_bit = 1u64 << target;
            if beads & target_bit != 0 {
                continue;
            }
            let source_bit = target_bit << length;
            let between = (source_bit - target_bit) ^ target_bit;
            let value = self.value(model, beads ^ source_bit ^ target_bit, classes >> 5)?;
            // Moving one beta-number down by length removes a rim hook.
            // Its height counts the beads strictly between the endpoints.
            let contribution = if (between & beads).count_ones().is_multiple_of(2) {
                value
            } else {
                -value
            };
            result = result.checked_add(contribution).ok_or_else(overflow)?;
        }
        if self.0.len() == MEMO_LIMIT {
            return Err(Error::plain(
                "Schur character table exceeds the storage limit",
            ));
        }
        self.0.insert((beads, classes), result);
        Ok(result)
    }
}

pub(super) fn beta_numbers(rows: impl Iterator<Item = u8>) -> u64 {
    rows.chain(std::iter::repeat(0))
        .take(MAX_WEIGHT as usize)
        .enumerate()
        .fold(0, |set, (i, row)| {
            set | (1 << (usize::from(row + MAX_WEIGHT - 1) - i))
        })
}

pub(super) fn pack(cycle_type: &[u8]) -> u128 {
    cycle_type
        .iter()
        .rev()
        .fold(0, |packed, &part| (packed << 5) | u128::from(part))
}

pub(super) fn partitions(model: &mut Model, n: u8) -> Result<Vec<Partition>> {
    fn visit(
        model: &mut Model,
        remaining: u8,
        max: u8,
        prefix: &mut Partition,
        output: &mut Vec<Partition>,
    ) -> Result<()> {
        if remaining == 0 {
            model.spend(prefix.len() + 1)?;
            output.push(prefix.clone());
            return Ok(());
        }
        for row in (1..=remaining.min(max)).rev() {
            prefix.push(row);
            visit(model, remaining - row, row, prefix, output)?;
            prefix.pop();
        }
        Ok(())
    }
    let mut values = Vec::new();
    visit(model, n, n, &mut Vec::new(), &mut values)?;
    Ok(values)
}

pub(super) fn centralizer(partition: &[u8]) -> u128 {
    partition
        .chunk_by(|a, b| a == b)
        .map(|run| {
            let count = run.len() as u8;
            u128::from(run[0]).pow(u32::from(count)) * factorial(count)
        })
        .product()
}

pub(super) fn factorial(n: u8) -> u128 {
    (1..=n).map(u128::from).product()
}
