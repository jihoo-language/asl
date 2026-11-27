use crate::error::{Error, Result};
use crate::model::Model;
use crate::wide_integer::Natural;

use super::{Factor, Gram, Integer, State, COLORS, LIMBS};

pub(super) fn pair(model: &mut Model, gram: &Gram, left: &State, right: &State) -> Result<i128> {
    let mut sum = Integer::from_i128(0);
    for (left, a) in left.terms() {
        for (right, b) in right.terms() {
            model.spend(1 + left.len() + right.len())?;
            let weight = monomial_pair(model, gram, left, right)?;
            if !weight.is_zero() {
                let term = multiply(model, weight, Integer::product(a, b))?;
                sum = sum.add(term).map_err(|_| intermediate_range())?;
            }
        }
    }
    sum.signed()
        .ok_or_else(|| Error::plain("Nakajima pairing exceeds signed 128-bit range"))
}

fn monomial_pair(
    model: &mut Model,
    gram: &Gram,
    left: &[Factor],
    right: &[Factor],
) -> Result<Integer> {
    let left = left.chunk_by(|a, b| a.mode.length == b.mode.length);
    let right = right.chunk_by(|a, b| a.mode.length == b.mode.length);
    if !left.clone().map(profile).eq(right.clone().map(profile)) {
        return Ok(Integer::from_i128(0));
    }
    let mut contracted = Vec::new();
    for (rows, columns) in left.zip(right) {
        let graph = Contractions {
            gram,
            rows,
            columns,
        };
        let mut counts = [0; COLORS];
        for (i, column) in columns.iter().enumerate() {
            counts[i] = column.exponent;
        }
        let coefficient = graph.coefficient(model, 0, &mut counts)?;
        if coefficient.is_zero() {
            return Ok(coefficient);
        }
        contracted.push((coefficient, columns));
    }
    let mut result = Integer::from_i128(1);
    for (coefficient, columns) in contracted {
        let (mode, count) = profile(columns);
        let signed_mode = if mode % 2 == 0 { -mode } else { mode };
        result = multiply(model, result, coefficient)?;
        let mode_factor = power(model, signed_mode, count)?;
        result = multiply(model, result, mode_factor)?;
        // A monomial coefficient counts allocations to colors. The factorial
        // of each right multiplicity assigns the distinguishable copies of that color.
        for column in columns {
            let assignments = factorial(model, column.exponent)?;
            result = multiply(model, result, assignments)?;
        }
    }
    Ok(result)
}

fn profile(factors: &[Factor]) -> (i128, u128) {
    (
        factors[0].mode.length,
        factors.iter().map(|factor| factor.exponent).sum(),
    )
}

struct Contractions<'a> {
    gram: &'a Gram,
    rows: &'a [Factor],
    columns: &'a [Factor],
}

#[derive(Clone, Copy)]
struct Slot {
    ordered_index: usize,
    gram_entry: i128,
}

impl Contractions<'_> {
    // Take the coefficient of the remaining right monomial in the product of
    // powers of the linear forms given by the remaining rows of the Gram matrix.
    fn coefficient(
        &self,
        model: &mut Model,
        row: usize,
        counts: &mut [u128; COLORS],
    ) -> Result<Integer> {
        let Some(factor) = self.rows.get(row) else {
            return Ok(Integer::from_i128(1));
        };
        model.spend(1 + self.columns.len())?;
        let mut slots = Vec::with_capacity(self.columns.len());
        let mut available = 0;
        for (i, column) in self.columns.iter().enumerate() {
            let entry = self.gram.entry(factor.mode.color, column.mode.color);
            if entry != 0 && counts[i] != 0 {
                available += counts[i];
                slots.push(Slot {
                    ordered_index: i,
                    gram_entry: entry,
                });
            }
        }
        if available < factor.exponent {
            return Ok(Integer::from_i128(0));
        }
        self.distribute(model, row, counts, &slots, available, factor.exponent)
    }

    fn distribute(
        &self,
        model: &mut Model,
        row: usize,
        counts: &mut [u128; COLORS],
        slots: &[Slot],
        available: u128,
        needed: u128,
    ) -> Result<Integer> {
        if needed == 0 {
            return self.coefficient(model, row + 1, counts);
        }
        let (slot, rest) = slots.split_first().expect("positive capacity for row");
        let capacity = counts[slot.ordered_index];
        let remaining = available - capacity;
        let mut sum = Integer::from_i128(0);
        for amount in needed.saturating_sub(remaining)..=needed.min(capacity) {
            model.spend(1)?;
            counts[slot.ordered_index] = capacity - amount;
            let next = self.distribute(model, row, counts, rest, remaining, needed - amount);
            counts[slot.ordered_index] = capacity;
            let tail = next?;
            if !tail.is_zero() {
                let ways = binomial(model, needed, amount)?;
                let pairing = power(model, slot.gram_entry, amount)?;
                let term = multiply(model, ways, pairing)?;
                let term = multiply(model, term, tail)?;
                sum = sum.add(term).map_err(|_| intermediate_range())?;
            }
        }
        Ok(sum)
    }
}

fn multiply(model: &mut Model, left: Integer, right: Integer) -> Result<Integer> {
    model.spend(1 + left.magnitude.bits().div_ceil(32) * right.magnitude.bits().div_ceil(32))?;
    left.multiply_wide(right).map_err(|_| intermediate_range())
}

fn power(model: &mut Model, value: i128, mut exponent: u128) -> Result<Integer> {
    let mut result = Integer::from_i128(1);
    let mut factor = Integer::from_i128(value);
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = multiply(model, result, factor)?;
        }
        exponent >>= 1;
        if exponent != 0 {
            factor = multiply(model, factor, factor)?;
        }
    }
    Ok(result)
}

fn factorial(model: &mut Model, exponent: u128) -> Result<Integer> {
    let mut result = Integer::from_i128(1);
    // This loop stops at the exact arithmetic bound even for huge multiplicities.
    for i in 2..=exponent {
        model.spend(LIMBS)?;
        result = result
            .multiply(i as i128)
            .map_err(|_| intermediate_range())?;
    }
    Ok(result)
}

fn binomial(model: &mut Model, n: u128, k: u128) -> Result<Integer> {
    // Reserve four further limbs for multiplication by a single u128 before
    // division, so the width test applies to the binomial itself.
    let mut value = Natural::<LIMBS>::from_small(1);
    for j in 1..=k.min(n - k) {
        model.spend(LIMBS)?;
        let extended = value.resize::<{ LIMBS + 4 }>().expect("widening a natural");
        let weighted = extended
            .multiply_u128(n - j + 1)
            .expect("product fits the wider natural");
        let (quotient, remainder) = weighted.divide(model, j)?;
        debug_assert_eq!(remainder, 0);
        value = quotient.resize().ok_or_else(intermediate_range)?;
    }
    Ok(Integer {
        magnitude: value,
        negative: false,
    })
}

fn intermediate_range() -> Error {
    Error::plain("Nakajima Wick pairing intermediate exceeds 512 bits")
}
