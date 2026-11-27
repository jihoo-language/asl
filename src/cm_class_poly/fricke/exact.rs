use crate::error::Result;
use crate::model::Model;
use crate::wide_integer::Natural;

use super::super::{Integer, Poly, LIMBS};
use super::certificate;

pub(super) fn integers(coefficients: &[i128]) -> Poly {
    let mut polynomial = Poly(
        coefficients
            .iter()
            .copied()
            .map(Integer::from_i128)
            .collect(),
    );
    polynomial.trim();
    polynomial
}

pub(super) fn zero() -> Poly {
    Poly(Vec::new())
}

pub(super) fn plus(a: &Poly, b: &Poly) -> Result<Poly> {
    let width = a.0.len().max(b.0.len());
    let mut result = Poly(
        (0..width)
            .map(|i| a.coefficient(i).add(b.coefficient(i)))
            .collect::<Result<_>>()?,
    );
    result.trim();
    Ok(result)
}

pub(super) fn scale(a: &Poly, factor: i128) -> Result<Poly> {
    let mut result = Poly(
        a.0.iter()
            .map(|c| c.multiply(factor))
            .collect::<Result<_>>()?,
    );
    result.trim();
    Ok(result)
}

pub(super) fn degree(a: &Poly) -> usize {
    a.0.len().saturating_sub(1)
}

pub(super) fn unit(a: &Poly) -> bool {
    a.0.len() == 1 && a.is_monic()
}

fn derivative(poly: &Poly) -> Result<Poly> {
    Ok(Poly(
        poly.0
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, value)| value.multiply(i as i128))
            .collect::<Result<_>>()?,
    ))
}

fn content_primitive(model: &mut Model, mut poly: Poly) -> Result<Poly> {
    poly.trim();
    let Some(lead) = poly.0.last() else {
        return Ok(poly);
    };
    let mut content = lead.magnitude;
    let negative = lead.negative;
    let one = Natural::<LIMBS>::from_small(1);
    for coefficient in poly.0.iter().rev().skip(1) {
        if content == one {
            break;
        }
        if !coefficient.is_zero() {
            content = content.gcd(model, coefficient.magnitude)?;
        }
    }
    if content != one {
        for coefficient in &mut poly.0 {
            coefficient.magnitude = coefficient
                .magnitude
                .divide_exact(model, content)?
                .ok_or_else(certificate)?;
        }
    }
    if negative {
        poly.0
            .iter_mut()
            .for_each(|coefficient| *coefficient = coefficient.negate());
    }
    Ok(poly)
}

pub(super) fn gcd(model: &mut Model, a: &Poly, b: &Poly) -> Result<Poly> {
    let mut a = content_primitive(model, a.clone())?;
    let mut b = content_primitive(model, b.clone())?;
    while let Some(&right) = b.0.last() {
        if b.0.len() == 1 {
            return Ok(Poly::one());
        }
        while let Some(offset) = a.0.len().checked_sub(b.0.len()) {
            let left = *a.0.last().unwrap();
            let common = left.magnitude.gcd(model, right.magnitude)?;
            let lhs = right
                .magnitude
                .divide_exact(model, common)?
                .ok_or_else(certificate)?;
            let rhs = left
                .magnitude
                .divide_exact(model, common)?
                .ok_or_else(certificate)?;
            let lhs = Integer {
                magnitude: lhs,
                negative: false,
            };
            let rhs = Integer {
                magnitude: rhs,
                negative: left.negative,
            };
            model.spend(LIMBS * (a.0.len() + b.0.len()))?;
            for coefficient in &mut a.0 {
                *coefficient = coefficient.multiply_wide(lhs)?;
            }
            for (coefficient, &subtract) in a.0[offset..].iter_mut().zip(&b.0) {
                *coefficient = coefficient.add(subtract.multiply_wide(rhs)?.negate())?;
            }
            a = content_primitive(model, a)?;
        }
        (a, b) = (b, a);
    }
    // Every caller supplies a monic operand. Gauss's lemma makes the
    // primitive positive gcd monic, so it divides in Z[H].
    if !a.is_monic() {
        return Err(certificate());
    }
    Ok(a)
}

pub(super) fn strata(model: &mut Model, poly: &Poly) -> Result<Vec<(usize, Poly)>> {
    let mut repeated = gcd(model, poly, &derivative(poly)?)?;
    let mut distinct = poly.clone().exact_divide(model, &repeated)?;
    let mut result = Vec::new();
    let mut product = Poly::one();
    let mut exponent = 1;
    while !unit(&distinct) {
        let shared = gcd(model, &distinct, &repeated)?;
        let part = distinct.exact_divide(model, &shared)?;
        if !unit(&part) {
            let factor = part.power(model, exponent)?;
            product = product.multiply(model, &factor)?;
            result.push((exponent, part));
        }
        repeated = repeated.exact_divide(model, &shared)?;
        distinct = shared;
        exponent += 1;
    }
    if !unit(&repeated) || !poly.same(&product) {
        return Err(certificate());
    }
    Ok(result)
}

pub(super) fn radical(model: &mut Model, strata: &[(usize, Poly)]) -> Result<Poly> {
    strata.iter().try_fold(Poly::one(), |product, (_, factor)| {
        product.multiply(model, factor)
    })
}

pub(super) fn show(p: &Poly) -> String {
    format!(
        "[{}]",
        p.0.iter()
            .map(|&c| super::super::decimal(c))
            .collect::<Vec<_>>()
            .join(", ")
    )
}
