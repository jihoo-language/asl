use crate::model::Model;
use crate::rv_qp_quad::Wide;

use super::{arithmetic, certificate, Result};

const LIMBS: usize = 128;
pub(super) type Integer = Wide<LIMBS>;

pub(super) fn power(
    model: &mut Model,
    polynomial: &[Integer],
    exponent: usize,
) -> Result<(Vec<Integer>, Integer)> {
    let degree = polynomial.len() - 1;
    let maximum = degree * exponent;
    let zero = Integer::from_i128(0);
    let mut traces = vec![zero; maximum + 1];
    for k in 1..=maximum {
        let mut sum = zero;
        for j in 1..=degree.min(k - 1) {
            sum = add(sum, multiply(model, polynomial[j], traces[k - j])?)?;
        }
        if k <= degree {
            sum = add(sum, scale(model, polynomial[k], k as i128)?)?;
        }
        traces[k] = sum.negate();
    }
    if exponent == 1 {
        return Ok((polynomial.to_vec(), traces.get(1).copied().unwrap_or(zero)));
    }
    let mut powered = vec![Integer::from_i128(1)];
    for k in 1..=degree {
        let mut sum = zero;
        for j in 1..=k {
            sum = add(sum, multiply(model, powered[k - j], traces[exponent * j])?)?;
        }
        model.spend(LIMBS)?;
        let (magnitude, remainder) = sum.magnitude.divide_u32(k as u32);
        if remainder != 0 {
            return Err(certificate("each powered Newton coefficient is integral"));
        }
        powered.push(Integer {
            magnitude,
            negative: !sum.negative && !magnitude.is_zero(),
        });
    }
    Ok((powered, traces.get(exponent).copied().unwrap_or(zero)))
}

pub(super) fn reciprocity(model: &mut Model, coefficients: &[Integer], q: Integer) -> Result<()> {
    let d = coefficients.len() - 1;
    let genus = d / 2;
    if d != 2 * genus || coefficients[0].signed() != Some(1) {
        return Err(certificate("constant term 1 and even Frobenius degree"));
    }
    let mut power = Integer::from_i128(1);
    for i in (0..=genus).rev() {
        if !equal(
            coefficients[d - i],
            multiply(model, coefficients[i], power)?,
        ) {
            return Err(certificate("exact Frobenius reciprocity"));
        }
        if i != 0 {
            power = multiply(model, power, q)?;
        }
    }
    Ok(())
}

pub(super) fn q_power(model: &mut Model, ell: usize, exponent: usize) -> Result<Integer> {
    let mut value = Integer::from_i128(1);
    for _ in 0..exponent {
        value = scale(model, value, ell as i128)?;
    }
    Ok(value)
}

fn equal(a: Integer, b: Integer) -> bool {
    a.magnitude == b.magnitude && (a.negative == b.negative || a.is_zero())
}

pub(super) fn add(a: Integer, b: Integer) -> Result<Integer> {
    a.add(b).map_err(|_| arithmetic())
}

fn scale(model: &mut Model, a: Integer, b: i128) -> Result<Integer> {
    model.spend(LIMBS + a.magnitude.bits().div_ceil(32))?;
    a.multiply(b).map_err(|_| arithmetic())
}

fn multiply(model: &mut Model, a: Integer, b: Integer) -> Result<Integer> {
    model.spend(LIMBS + a.magnitude.bits().div_ceil(32) * b.magnitude.bits().div_ceil(32))?;
    a.multiply_wide(b).map_err(|_| arithmetic())
}
