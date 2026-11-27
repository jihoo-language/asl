use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::cm::primitive_class_number;
use crate::lee_linear::Rational;
use crate::model::Model;
use crate::siegel::discriminant;

use super::{decimal, multiply, overflow, scale_by_power, Integer};

pub(crate) enum Query {
    Class(i128),
    Jacobi(i128, i128),
    Trace(u128, u128),
}

impl Query {
    pub(crate) fn report(self, model: &mut Model) -> Result<String> {
        match self {
            Self::Class(n) => {
                let value = if n < 0 {
                    0
                } else {
                    twelve_hurwitz(model, n as u128)?
                };
                let value = Rational::new(value, 12)?;
                Ok(format!("rvHurwitzClass[{n}] = {value}"))
            }
            Self::Jacobi(n, r) => {
                let value = Rational::new(jacobi(model, n, r)?, 12)?;
                Ok(format!("rvHurwitzJacobi[{n}, {r}] = {value}"))
            }
            Self::Trace(k, n) => {
                let value = trace(model, k, n)?;
                Ok(format!("rvHeckeTrace[{k}, {n}] = {}", decimal(value)))
            }
        }
    }
}

fn jacobi(model: &mut Model, n: i128, r: i128) -> Result<i128> {
    if n < 0 {
        return Ok(0);
    }
    let d = discriminant(n as u128, 1, r)?;
    if d.negative && !d.is_zero() {
        return Ok(0);
    }
    twelve_hurwitz(model, d.magnitude.small().ok_or_else(class_limit)?)
}

fn twelve_hurwitz(model: &mut Model, n: u128) -> Result<i128> {
    if n == 0 {
        return Ok(-1);
    }
    if !valid_discriminant(n) {
        return Ok(0);
    }
    let mut sum = 0;
    for content in 1..=n.isqrt() {
        model.spend(1).map_err(|_| class_limit())?;
        let square = content * content;
        if !n.is_multiple_of(square) {
            continue;
        }
        let primitive = n / square;
        if valid_discriminant(primitive) {
            let count = primitive_class_number(model, primitive).map_err(|_| class_limit())?;
            let weight = match primitive {
                3 => 4,
                4 => 6,
                _ => 12,
            };
            sum += weight * count as i128;
        }
    }
    Ok(sum)
}

fn valid_discriminant(n: u128) -> bool {
    matches!(n % 4, 0 | 3)
}

fn trace(model: &mut Model, k: u128, n: u128) -> Result<Integer> {
    if k < 2 || !k.is_multiple_of(2) {
        return Err(Error::plain(
            "level-one cusp-form Hecke trace requires even weight at least 2",
        ));
    }
    if n == 0 {
        return Err(Error::plain("Hecke index n must be positive"));
    }
    let four_n = n.checked_mul(4).ok_or_else(class_limit)?;
    let mut sum = Integer::from_i128(0);
    // P_{k−2}(t,n) is even in t. Keep 24 times the trace integral until
    // the final exact division, including the H(0) endpoints.
    for t in 0..=four_n.isqrt() {
        let h = twelve_hurwitz(model, four_n - t * t)?;
        if h != 0 {
            let p = recurrence(model, k - 2, t as i128, n as i128)?;
            let multiplicity = if t == 0 { 1 } else { 2 };
            let term = multiply(model, p, Integer::from_i128(-multiplicity * h))?;
            sum = sum.add(term).map_err(|_| overflow())?;
        }
    }
    for d in 1..=n.isqrt() {
        model.spend(1)?;
        if n.is_multiple_of(d) {
            let other = n / d;
            let multiplicity = if d == other { 1 } else { 2 };
            let factor = scale_by_power(model, Integer::from_i128(-12 * multiplicity), d, k - 1)?;
            sum = sum.add(factor).map_err(|_| overflow())?;
            if k == 2 {
                let divisors = if d == other { d } else { d + other };
                let correction = Integer::from_i128(divisors as i128)
                    .multiply(24)
                    .map_err(|_| overflow())?;
                sum = sum.add(correction).map_err(|_| overflow())?;
            }
        }
    }
    let (magnitude, remainder) = sum.magnitude.divide_u32(24);
    if remainder != 0 {
        return Err(Error::plain("Eichler–Selberg trace division is inexact"));
    }
    Ok(Integer { magnitude, ..sum })
}

fn recurrence(model: &mut Model, degree: u128, t: i128, n: i128) -> Result<Integer> {
    let (mut current, mut previous) = (Integer::from_i128(1), Integer::from_i128(0));
    for bit in (0..u128::BITS - degree.leading_zeros()).rev() {
        // For a=P_j and b=P_{j−1}, P_{2j}=a²−n b² and
        // P_{2j−1}=b(2a−tb). Iterate the binary digits of the degree.
        let a_square = multiply(model, current, current)?;
        let minus_n_b_square = multiply(model, previous, previous)?
            .multiply(-n)
            .map_err(|_| overflow())?;
        let even = a_square.add(minus_n_b_square).map_err(|_| overflow())?;
        let twice_a = current.multiply(2).map_err(|_| overflow())?;
        let minus_t_b = previous.multiply(-t).map_err(|_| overflow())?;
        let before = multiply(
            model,
            previous,
            twice_a.add(minus_t_b).map_err(|_| overflow())?,
        )?;
        (current, previous) = if (degree >> bit) & 1 == 0 {
            (even, before)
        } else {
            let a = even.multiply(t).map_err(|_| overflow())?;
            let b = before.multiply(-n).map_err(|_| overflow())?;
            (a.add(b).map_err(|_| overflow())?, even)
        };
    }
    Ok(current)
}

fn class_limit() -> Error {
    Error::plain("Hurwitz class number enumeration exceeds the computation limit")
}
