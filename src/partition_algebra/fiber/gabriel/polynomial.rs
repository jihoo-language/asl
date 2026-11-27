use std::collections::BTreeSet;

use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::Model;

use super::super::{add, multiply, negate, ratio};

type Polynomial = Vec<Rational>;
const ONE: Rational = Rational::integer(1);

pub(super) fn splitting_projector(
    model: &mut Model,
    modulus: &[Rational],
    n: u128,
) -> Result<Option<Polynomial>> {
    if modulus.len() < 3 {
        return Ok(None);
    }
    let mut roots = BTreeSet::new();
    if modulus[0] == Rational::ZERO {
        roots.insert(Rational::ZERO);
    } else {
        for d in 1..=8 {
            for m in -16..=16 {
                roots.insert(ratio(m, d)?);
            }
        }
        if let Ok(n) = i128::try_from(n) {
            for shift in 0..=3 {
                if let Some(a) = n.checked_sub(shift) {
                    for d in 1..=4 {
                        roots.insert(ratio(a, d)?);
                        roots.insert(ratio(-a, d)?);
                    }
                }
            }
        }
    }
    for root in roots {
        // Remove the entire primary power, so the remaining factor is
        // coprime and provides a true CRT projection even on nilpotents.
        let linear = [negate(root)?, ONE];
        let mut complement = modulus.to_vec();
        let mut primary = vec![ONE];
        loop {
            let (quotient, remainder) = division(model, complement.clone(), &linear)?;
            if !remainder.is_empty() {
                break;
            }
            complement = quotient;
            primary = product(model, &primary, &linear)?;
        }
        if primary.len() > 1 && complement.len() > 1 {
            let inverse = inverse_mod(model, &complement, &primary)?;
            let value = product(model, &inverse, &complement)?;
            return Ok(Some(division(model, value, modulus)?.1));
        }
    }
    Ok(None)
}

fn divide(a: Rational, b: Rational) -> Result<Rational> {
    multiply(a, ratio(b.denominator, b.numerator)?)
}

fn trim(value: &mut Polynomial) {
    value.truncate(
        value
            .iter()
            .rposition(|c| *c != Rational::ZERO)
            .map_or(0, |i| i + 1),
    );
}

fn division(
    model: &mut Model,
    mut a: Polynomial,
    b: &[Rational],
) -> Result<(Polynomial, Polynomial)> {
    let leading = *b
        .last()
        .ok_or_else(|| Error::plain("partition polynomial division by zero"))?;
    if a.len() < b.len() {
        return Ok((Vec::new(), a));
    }
    let mut quotient = vec![Rational::ZERO; a.len() - b.len() + 1];
    for (i, value) in quotient.iter_mut().enumerate().rev() {
        model.spend(b.len())?;
        *value = divide(a[i + b.len() - 1], leading)?;
        if *value != Rational::ZERO {
            for (entry, &factor) in a[i..].iter_mut().zip(b) {
                *entry = add(*entry, negate(multiply(*value, factor)?)?)?;
            }
        }
    }
    trim(&mut quotient);
    trim(&mut a);
    Ok((quotient, a))
}

fn product(model: &mut Model, a: &[Rational], b: &[Rational]) -> Result<Polynomial> {
    if a.is_empty() || b.is_empty() {
        return Ok(Vec::new());
    }
    model.spend(a.len() * b.len())?;
    let mut output = vec![Rational::ZERO; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            output[i + j] = add(output[i + j], multiply(x, y)?)?;
        }
    }
    trim(&mut output);
    Ok(output)
}

fn inverse_mod(model: &mut Model, value: &[Rational], modulus: &[Rational]) -> Result<Polynomial> {
    let mut a = modulus.to_vec();
    let mut b = value.to_vec();
    let mut s = Vec::new();
    let mut t = vec![ONE];
    while !b.is_empty() {
        let (q, r) = division(model, a, &b)?;
        let u = product(model, &q, &t)?;
        s.resize(s.len().max(u.len()), Rational::ZERO);
        for (s, u) in s.iter_mut().zip(u) {
            *s = add(*s, negate(u)?)?;
        }
        trim(&mut s);
        (a, b) = (b, r);
        (s, t) = (t, s);
    }
    let [constant] = a.as_slice() else {
        return Err(Error::plain(
            "partition polynomial primary factors are not coprime",
        ));
    };
    s.into_iter().map(|c| divide(c, *constant)).collect()
}
