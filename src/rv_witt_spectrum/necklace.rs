use crate::error::{Error, Result};
use crate::model::Model;

use super::{add, Integer};

const TRIAL_LIMIT: usize = 1_000_000;

pub(crate) fn coordinate(
    model: &mut Model,
    index: u128,
    mut ghost: impl FnMut(&mut Model, u128) -> Result<Integer>,
) -> Result<Integer> {
    if index == 0 {
        return Err(Error::plain("primitive period must be positive"));
    }
    // Evaluate the largest argument first so that an arithmetic bound can
    // diagnose a growing sequence before integer factorization.
    let mut sum = ghost(model, index)?;
    if index == 1 {
        return Ok(sum);
    }
    let primes = distinct_primes(model, index)?;
    accumulate(model, &primes, index, true, &mut sum, &mut ghost)?;
    let (magnitude, remainder) = sum.magnitude.divide(model, index)?;
    if remainder != 0 {
        return Err(Error::plain(
            "primitive orbit Möbius sum fails exact division",
        ));
    }
    Ok(Integer {
        magnitude,
        negative: sum.negative,
    })
}

fn distinct_primes(model: &mut Model, mut n: u128) -> Result<Vec<u128>> {
    let mut factors = Vec::new();
    let mut candidate = 2;
    let mut trials = 0;
    while candidate <= n / candidate {
        if trials == TRIAL_LIMIT {
            return Err(Error::plain(
                "primitive period factorization exceeds the trial division limit",
            ));
        }
        model.spend(1)?;
        trials += 1;
        if n.is_multiple_of(candidate) {
            factors.push(candidate);
            while n.is_multiple_of(candidate) {
                model.spend(1)?;
                n /= candidate;
            }
        }
        candidate = if candidate == 2 { 3 } else { candidate + 2 };
    }
    if n > 1 {
        factors.push(n);
    }
    Ok(factors)
}

fn accumulate(
    model: &mut Model,
    primes: &[u128],
    period: u128,
    negative: bool,
    sum: &mut Integer,
    ghost: &mut impl FnMut(&mut Model, u128) -> Result<Integer>,
) -> Result<()> {
    for (i, &prime) in primes.iter().enumerate() {
        model.spend(1)?;
        let reduced = period / prime;
        let value = ghost(model, reduced)?;
        *sum = add(*sum, if negative { value.negate() } else { value })?;
        accumulate(model, &primes[i + 1..], reduced, !negative, sum, ghost)?;
    }
    Ok(())
}
