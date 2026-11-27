use crate::error::{Error, Result};
use crate::model::Model;

use super::{
    add, coefficient_overflow, degree_limit, determinant_from_ghosts, integrity, multiply, one,
    power_sums, require_degree, zero, Integer, Spectrum, MAX_DEGREE,
};

pub(crate) fn cycles(n: u128, k: u128) -> Result<(u128, u128)> {
    if n == 0 || k == 0 {
        return Err(Error::plain(
            "tensor-induction degree and ghost index must be positive",
        ));
    }
    let (mut g, mut rest) = (n, k);
    while rest != 0 {
        (g, rest) = (rest, g % rest);
    }
    // Rotation by k positions has g cycles. A factor traverses each cycle
    // k/g times before returning to its original tensor position.
    Ok((g, k / g))
}

pub(crate) fn power(model: &mut Model, mut value: Integer, mut exponent: u128) -> Result<Integer> {
    if exponent == 0 {
        return Ok(one());
    }
    if value.is_zero() {
        return Ok(zero());
    }
    if value.magnitude.small() == Some(1) {
        return Ok(Integer::from_i128(
            if value.negative && !exponent.is_multiple_of(2) {
                -1
            } else {
                1
            },
        ));
    }
    let mut product = one();
    loop {
        if !exponent.is_multiple_of(2) {
            product = multiply(model, product, value)?;
        }
        exponent /= 2;
        if exponent == 0 {
            return Ok(product);
        }
        value = multiply(model, value, value)?;
    }
}

impl Spectrum {
    pub(crate) fn norm(&self, model: &mut Model, n: u128) -> Result<Self> {
        if n == 0 {
            return Err(Error::plain("Witt spectrum norm degree must be positive"));
        }
        if n == 1 {
            return Ok(self.clone());
        }
        let a = self.numerator.len() - 1;
        let b = self.denominator.len() - 1;
        if a + b <= 1 {
            if b == 1 && n.is_multiple_of(2) {
                let coefficient = self.denominator[1]
                    .checked_neg()
                    .ok_or_else(coefficient_overflow)?;
                return Ok(Self {
                    numerator: Box::new([1, coefficient]),
                    denominator: Box::new([1]),
                });
            }
            return Ok(self.clone());
        }

        let (even, odd) = ranks(a, b, n)?;
        let rank = even.max(odd);
        // Precompute the ordinary power sums of the two input factors. At
        // most k/g terms of each sequence are needed for the k-th ghost.
        let cyclic = (1..=rank)
            .map(|k| cycles(n, k as u128))
            .collect::<Result<Vec<_>>>()?;
        let limit = cyclic.iter().map(|&(_, i)| i as usize).max().unwrap_or(0);
        let pos = power_sums(model, &self.numerator, limit)?;
        let neg = power_sums(model, &self.denominator, limit)?;
        model.spend(rank + even + odd)?;
        let mut even_ghosts = Vec::with_capacity(even);
        let mut odd_ghosts = Vec::with_capacity(odd);
        for (k, (g, i)) in cyclic.into_iter().enumerate() {
            let left = pos[i as usize];
            let right = neg[i as usize];
            let difference = add(left, right.negate())?;
            let supertrace = power(model, difference, g)?;
            let (even_trace, odd_trace) = if odd == 0 || (n / g).is_multiple_of(2) {
                (supertrace, zero())
            } else if even == 0 {
                (zero(), supertrace.negate())
            } else {
                let trace = power(model, add(left, right)?, g)?;
                (
                    half(add(trace, supertrace)?)?,
                    half(add(trace, supertrace.negate())?)?,
                )
            };
            if k < even {
                even_ghosts.push(even_trace);
            }
            if k < odd {
                odd_ghosts.push(odd_trace);
            }
        }
        let p = determinant_from_ghosts(model, &even_ghosts)?;
        let q = determinant_from_ghosts(model, &odd_ghosts)?;
        if p[even].is_zero() || q[odd].is_zero() {
            return Err(integrity());
        }
        Self::canonical(model, p, q)
    }

    pub(super) fn norm_ghost(&self, model: &mut Model, n: u128, k: u128) -> Result<Integer> {
        let (g, i) = cycles(n, k)?;
        let ghost = self.ghost(model, i)?;
        power(model, ghost, g)
    }
}

fn ranks(a: usize, b: usize, n: u128) -> Result<(usize, usize)> {
    let exponent = u32::try_from(n).map_err(|_| degree_limit())?;
    let total = (a + b)
        .checked_pow(exponent)
        .filter(|&total| total <= 2 * MAX_DEGREE)
        .ok_or_else(degree_limit)?;
    let imbalance = a
        .abs_diff(b)
        .checked_pow(exponent)
        .ok_or_else(degree_limit)?;
    if !(total - imbalance).is_multiple_of(2) {
        return Err(Error::plain(
            "tensor-induction parity ranks are not integral",
        ));
    }
    let (larger, smaller) = ((total + imbalance) / 2, (total - imbalance) / 2);
    require_degree(larger)?;
    if a < b && !n.is_multiple_of(2) {
        Ok((smaller, larger))
    } else {
        Ok((larger, smaller))
    }
}

pub(super) fn half(value: Integer) -> Result<Integer> {
    let (magnitude, remainder) = value.magnitude.divide_u32(2);
    if remainder != 0 {
        return Err(Error::plain("graded trace division by two is inexact"));
    }
    Ok(Integer {
        magnitude,
        negative: value.negative,
    })
}
