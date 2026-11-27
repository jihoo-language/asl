use std::collections::{BTreeMap, HashMap};

use crate::error::{Error, Result};
use crate::model::Model;
use crate::wide_integer::Natural;

use super::super::norm::power;
use super::super::{add, exact_divide, multiply, zero, Integer, SchurShape};
use super::{output_weight, overflow, Character};

const MAX_WEIGHT: u8 = 16;
const MEMO_LIMIT: usize = 200_000;
type Partition = Vec<u8>;
type Coefficients = BTreeMap<Partition, Integer>;

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
        let inner_terms = inner.class_terms(model, &mut oracle, b)?;
        let outer_terms = self.class_terms(model, &mut oracle, a)?;
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
        let terms = self.class_terms(model, &mut oracle, n)?;
        let value = Integer::from_i128(value);
        let mut sum = zero();
        for (cycles, coefficient) in terms {
            let substitute = power(model, value, cycles.len() as u128)?;
            sum = add(sum, multiply(model, coefficient, substitute)?)?;
        }
        Self::constant(divide(model, sum, factorial(n))?)
    }

    // In the power basis the class value on ρ is divided by z_ρ.
    // Store every coefficient with a common denominator n!.
    fn class_terms(
        &self,
        model: &mut Model,
        oracle: &mut CharacterOracle,
        n: u8,
    ) -> Result<Coefficients> {
        let factorial = factorial(n);
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
                let scale = Integer::from_i128((factorial / centralizer(&partition)) as i128);
                terms.insert(partition, multiply(model, value, scale)?);
            }
        }
        Ok(terms)
    }
}

fn from_power_basis(
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

fn divide(model: &mut Model, value: Integer, denominator: u128) -> Result<i128> {
    exact_divide(model, value, Natural::from_small(denominator))?
        .signed()
        .ok_or_else(overflow)
}

fn accumulate(terms: &mut Coefficients, partition: Partition, value: Integer) -> Result<()> {
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

#[derive(Default)]
struct CharacterOracle(HashMap<(u64, u128), i128>);

impl CharacterOracle {
    fn value(&mut self, model: &mut Model, beads: u64, classes: u128) -> Result<i128> {
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
                "Schur plethysm character table exceeds the storage limit",
            ));
        }
        self.0.insert((beads, classes), result);
        Ok(result)
    }
}

fn beta_numbers(rows: impl Iterator<Item = u8>) -> u64 {
    rows.chain(std::iter::repeat(0))
        .take(MAX_WEIGHT as usize)
        .enumerate()
        .fold(0, |set, (i, row)| {
            set | (1 << (usize::from(row + MAX_WEIGHT - 1) - i))
        })
}

fn pack(cycle_type: &[u8]) -> u128 {
    cycle_type
        .iter()
        .rev()
        .fold(0, |packed, &part| (packed << 5) | u128::from(part))
}

fn partitions(model: &mut Model, n: u8) -> Result<Vec<Partition>> {
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

fn centralizer(partition: &[u8]) -> u128 {
    partition
        .chunk_by(|a, b| a == b)
        .map(|run| {
            let count = run.len() as u8;
            u128::from(run[0]).pow(u32::from(count)) * factorial(count)
        })
        .product()
}

fn factorial(n: u8) -> u128 {
    (1..=n).map(u128::from).product()
}

fn bounded(n: u128) -> Result<u8> {
    if n > u128::from(MAX_WEIGHT) {
        return Err(Error::plain(
            "Schur plethysm exceeds the weight limit of 16 boxes",
        ));
    }
    Ok(n as u8)
}
