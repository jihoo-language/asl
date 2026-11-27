use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::hyperelliptic_jacobian::Divisor;
use crate::hyperelliptic_kummer_surface::hilbert::{
    checked_product, coefficients, intermediate_overflow, power, unsigned, Integer, Magnitude,
    SurfaceKind,
};
use crate::jacobian_level::Basis;
use crate::model::Model;

use super::{require_level_basis, require_translation_order, torsion_fixed};

pub(crate) fn hilbert_betti(model: &mut Model, n: u128, i: u128) -> Result<u128> {
    let Some((degree, _)) = dual_degree(n, i) else {
        return Ok(0);
    };
    if degree == 0 {
        return Ok(1);
    }
    // After the degree-zero generator of length one, every generator has
    // cohomological degree at least its length. Low degrees stabilize.
    let rows = coefficients(
        model,
        n.min(degree),
        Some(degree),
        None,
        SurfaceKind::Abelian,
    )?;
    betti_value(last_coefficient(&rows), degree)
}

pub(crate) fn hilbert_trace(
    model: &mut Model,
    curve: &Curve,
    n: u128,
    r: u128,
    i: u128,
) -> Result<i128> {
    let Some((degree, twist)) = dual_degree(n, i) else {
        return Ok(0);
    };
    let value = if degree == 0 {
        Integer::from_i128(1)
    } else {
        let rows = coefficients(
            model,
            n.min(degree),
            Some(degree),
            Some((curve, r)),
            SurfaceKind::Abelian,
        )?;
        last_coefficient(&rows)
    };
    trace_value(model, curve, r, degree, twist, value)
}

pub(crate) fn betti(model: &mut Model, n: u128, i: u128) -> Result<u128> {
    betti_value(geometric_value(model, n, 1, i)?, i)
}

pub(crate) fn translation_character(model: &mut Model, n: u128, h: u128, i: u128) -> Result<i128> {
    require_translation_order(n, h)?;
    let value = geometric_value(model, n, h, i)?;
    cohomology_sign(value, i).signed().ok_or_else(|| {
        Error::plain("generalized Kummer translation character exceeds signed 128-bit range")
    })
}

fn geometric_value(model: &mut Model, n: u128, h: u128, i: u128) -> Result<Integer> {
    let Some((degree, _)) = dual_degree(n - 1, i) else {
        return Ok(Integer::from_i128(0));
    };
    if degree == 0 {
        return Ok(Integer::from_i128(1));
    }
    let length = n.min(degree);
    let rows = coefficients(model, length, Some(degree), None, SurfaceKind::Abelian)?;
    let mut numerator = rows[length as usize].clone();
    if n <= degree {
        let n = n as usize; // The coefficient table already bounds this length.
        let top = degree as usize;
        for (k, shift, exact) in character_terms(model, n, top, h as usize)? {
            if exact == 0 {
                continue;
            }
            for (coefficient, &source) in numerator[shift..].iter_mut().zip(&rows[k]) {
                let term = checked_product(model, Integer::from_i128(exact), source)?;
                *coefficient = coefficient.add(term).map_err(|_| intermediate_overflow())?;
            }
        }
    }
    quotient(model, numerator, &rows[1])
}

// Length, cohomological shift, and exact-order torsion count for the
// nontrivial divisor terms of the smooth generalized Kummer identity.
pub(crate) fn geometric_torsion_terms(
    model: &mut Model,
    n: usize,
    top: usize,
) -> Result<Vec<(usize, usize, u128)>> {
    character_terms(model, n, top, 1)?
        .into_iter()
        .map(|(k, shift, exact)| {
            let count = u128::try_from(exact)
                .map_err(|_| Error::plain("exact-order torsion count must be nonnegative"))?;
            Ok((k, shift, count))
        })
        .collect()
}

fn character_terms(
    model: &mut Model,
    n: usize,
    top: usize,
    order: usize,
) -> Result<Vec<(usize, usize, i128)>> {
    let mut orders = vec![(1, 1i128)];
    let mut terms = Vec::new();
    for d in 2..=n {
        model.spend(1)?;
        if !n.is_multiple_of(d) {
            continue;
        }
        let k = n / d;
        let shift = 2 * (n - k);
        if shift > top {
            break;
        }
        // Characters of order dividing d sum to d^4 when the translation
        // belongs to d*A[n], which means its order divides n/d.
        let mut exact = if (n / order).is_multiple_of(d) {
            (d as i128)
                .checked_pow(4)
                .ok_or_else(intermediate_overflow)?
        } else {
            0
        };
        // Invert the divisor sums to get the trace on exact-order characters.
        model.spend(orders.len())?;
        for &(e, count) in &orders {
            if d.is_multiple_of(e) {
                exact = exact.checked_sub(count).ok_or_else(intermediate_overflow)?;
            }
        }
        orders.push((d, exact));
        terms.push((k, shift, exact));
    }
    Ok(terms)
}

pub(crate) fn trace(
    model: &mut Model,
    basis: &Basis,
    curve: &Curve,
    r: u128,
    i: u128,
    translation: Option<&Divisor>,
) -> Result<i128> {
    require_level_basis(basis)?;
    basis.require_base_extension(model, curve)?;
    let coordinates = translation
        .map(|t| basis.translation_coordinates(model, t))
        .transpose()?;
    let n = u128::from(basis.prime.get());
    let Some((degree, twist)) = dual_degree(n - 1, i) else {
        return Ok(0);
    };
    let value = match degree {
        0 => Integer::from_i128(1),
        1 => return Ok(0),
        _ => {
            let length = n.min(degree);
            let rows = coefficients(
                model,
                length,
                Some(degree),
                Some((curve, r)),
                SurfaceKind::Abelian,
            )?;
            let value = quotient(model, rows[length as usize].clone(), &rows[1])?;
            if degree == 2 * (n - 1) {
                let power_index = usize::try_from(r).map_err(|_| intermediate_overflow())?;
                let fixed =
                    torsion_fixed(model, basis, curve, power_index, coordinates.as_deref())?;
                let factor = if fixed == 0 {
                    Integer::from_i128(-1)
                } else {
                    nonnegative(fixed - 1)
                };
                let correction = tate(model, curve, r, n - 1, factor)?;
                value.add(correction).map_err(|_| intermediate_overflow())?
            } else {
                value
            }
        }
    };
    trace_value(model, curve, r, degree, twist, value)
}

// The space has dimension 2*length. Mirror high degrees by duality while
// keeping the range check safe even if 4*length is wider than u128.
pub(crate) fn dual_degree(length: u128, i: u128) -> Option<(u128, u128)> {
    if i.div_ceil(4) > length {
        return None;
    }
    let middle = length.saturating_mul(2);
    let twist = i.saturating_sub(middle);
    Some((i.min(middle) - twist, twist))
}

fn last_coefficient(rows: &[Vec<Integer>]) -> Integer {
    *rows.last().unwrap().last().unwrap()
}

fn quotient(
    model: &mut Model,
    mut numerator: Vec<Integer>,
    abelian: &[Integer],
) -> Result<Integer> {
    // The signed cohomology polynomial of A has constant term 1.
    for i in 1..numerator.len() {
        for j in 1..=i.min(4) {
            let term = checked_product(model, abelian[j], numerator[i - j])?;
            numerator[i] = numerator[i]
                .add(term.negate())
                .map_err(|_| intermediate_overflow())?;
        }
    }
    Ok(*numerator.last().unwrap())
}

fn betti_value(value: Integer, degree: u128) -> Result<u128> {
    unsigned(cohomology_sign(value, degree)).ok_or_else(|| {
        Error::plain(
            "Jacobian Hilbert or generalized Kummer Betti number exceeds unsigned 128-bit range",
        )
    })
}

fn trace_value(
    model: &mut Model,
    curve: &Curve,
    r: u128,
    degree: u128,
    twist: u128,
    value: Integer,
) -> Result<i128> {
    let value = tate(model, curve, r, twist, cohomology_sign(value, degree))?;
    value.signed().ok_or_else(|| {
        Error::plain(
            "Jacobian Hilbert or generalized Kummer Frobenius trace exceeds signed 128-bit range",
        )
    })
}

fn cohomology_sign(value: Integer, degree: u128) -> Integer {
    if degree.is_multiple_of(2) {
        value
    } else {
        value.negate()
    }
}

fn nonnegative(value: u128) -> Integer {
    Integer {
        magnitude: Magnitude::from_small(value),
        negative: false,
    }
}

fn tate(model: &mut Model, curve: &Curve, r: u128, twist: u128, value: Integer) -> Result<Integer> {
    if twist == 0 || value.is_zero() {
        return Ok(value);
    }
    let exponent = r.checked_mul(twist).ok_or_else(intermediate_overflow)?;
    let magnitude = power(model, curve.ring.size_of_residue_field()?, exponent)?;
    checked_product(
        model,
        value,
        Integer {
            magnitude,
            negative: false,
        },
    )
}
