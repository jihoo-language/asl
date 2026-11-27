use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_ore::gcd;

use super::{
    add, multiply, multiply_factor, overflow, storage, Cutoff, Integer, Product, Seed, TERM_LIMIT,
};

type Laurent = BTreeMap<i128, Integer>;

pub(super) fn coefficient(model: &mut Model, n: u128, degree: u128, r: i128) -> Result<Integer> {
    if n == 0 {
        return Err(Error::plain("generalized Kummer length must be positive"));
    }
    if n == 1 {
        return Ok(Integer::from_i128(i128::from(degree == 0 && r == 0)));
    }
    let index = n - 1;
    let discriminant_bound = index
        .checked_mul(degree)
        .and_then(|nd| nd.checked_mul(4))
        .and_then(|nd| {
            index
                .checked_mul(index)
                .and_then(|square| square.checked_add(nd))
        });
    if (degree == 0 && r.unsigned_abs() > index)
        || discriminant_bound.is_some_and(|bound| r.unsigned_abs() > bound.isqrt())
    {
        return Ok(Integer::from_i128(0));
    }
    let n = usize::try_from(n).map_err(|_| storage())?;
    let degree = usize::try_from(degree).map_err(|_| storage())?;
    // Bound the q and Laurent tables together before expanding any factors.
    let primitive = n.checked_mul(degree).ok_or_else(storage)?;
    let radius = n.saturating_add(primitive);
    let rows = degree
        .saturating_add(1)
        .saturating_mul(radius.saturating_mul(2).saturating_add(1));
    if rows.saturating_add(primitive.saturating_add(1)) > TERM_LIMIT {
        return Err(storage());
    }
    let cutoff = Cutoff::new(model, 0, primitive as u128)?;
    let seed = jacobi_seed(model, cutoff)?;
    let mut quotient = Seed::new();
    for d in 0..=degree {
        let mut residual = hecke_row(model, n, d, &seed)?;
        for j in 1..=d {
            for (r, a) in seed_row(&seed, j) {
                for (&(_, s), &e) in quotient.range((d - j, i128::MIN)..=(d - j, i128::MAX)) {
                    let term = multiply(model, a, e)?.negate();
                    add(model, &mut residual, r + s, term)?;
                }
            }
        }
        // A0 = (y-1)^2 / y. Both simple divisions check their remainders.
        let once = divide_by_y_minus_one(model, &residual)?;
        let twice = divide_by_y_minus_one(model, &once)?;
        for (r, value) in twice {
            add(model, &mut quotient, (d, r + 1), value)?;
        }
    }
    Ok(quotient
        .get(&(degree, r))
        .copied()
        .unwrap_or_else(|| Integer::from_i128(0)))
}

fn jacobi_seed(model: &mut Model, cutoff: Cutoff) -> Result<Product> {
    let mut seed = Product::from([
        ((0, 0, -1), Integer::from_i128(1)),
        ((0, 0, 0), Integer::from_i128(-2)),
        ((0, 0, 1), Integer::from_i128(1)),
    ]);
    for j in 1..=cutoff.q {
        for (r, exponent) in [(-1, 2), (0, -4), (1, 2)] {
            multiply_factor(
                model,
                &mut seed,
                cutoff,
                (0, j, r),
                Integer::from_i128(exponent),
            )?;
        }
    }
    Ok(seed)
}

fn seed_row(seed: &Product, degree: usize) -> impl Iterator<Item = (i128, Integer)> + '_ {
    seed.range((0, degree, i128::MIN)..=(0, degree, i128::MAX))
        .map(|(&(_, _, r), &value)| (r, value))
}

fn hecke_row(model: &mut Model, n: usize, d: usize, base: &Product) -> Result<Laurent> {
    let common = gcd(n as u128, d as u128) as usize;
    let mut row = Laurent::new();
    for a in 1..=common.isqrt() {
        model.spend(1)?;
        if !common.is_multiple_of(a) {
            continue;
        }
        let other = common / a;
        add_hecke_divisor(model, &mut row, base, n, d, a)?;
        if other != a {
            add_hecke_divisor(model, &mut row, base, n, d, other)?;
        }
    }
    Ok(row)
}

fn add_hecke_divisor(
    model: &mut Model,
    row: &mut Laurent,
    base: &Product,
    n: usize,
    d: usize,
    a: usize,
) -> Result<()> {
    let quotient = n / a;
    let seed_degree = quotient * (d / a);
    let scale = (quotient as i128)
        .checked_pow(3)
        .and_then(|cube| cube.checked_mul(n as i128))
        .ok_or_else(overflow)?;
    for (r, value) in seed_row(base, seed_degree) {
        let term = multiply(model, value, Integer::from_i128(scale))?;
        add(model, row, r * a as i128, term)?;
    }
    Ok(())
}

fn divide_by_y_minus_one(model: &mut Model, input: &Laurent) -> Result<Laurent> {
    let Some((&first, _)) = input.first_key_value() else {
        return Ok(Laurent::new());
    };
    let (&last, _) = input.last_key_value().unwrap();
    let mut output = Laurent::new();
    let mut carry = Integer::from_i128(0);
    model.spend(usize::try_from(last - first + 1).map_err(|_| storage())?)?;
    for r in first..=last {
        if let Some(&coefficient) = input.get(&r) {
            carry = carry.add(coefficient.negate()).map_err(|_| overflow())?;
        }
        if r != last {
            add(model, &mut output, r, carry)?;
        }
    }
    if !carry.is_zero() {
        return Err(Error::plain(
            "generalized Kummer elliptic genus Laurent division is inexact",
        ));
    }
    Ok(output)
}
