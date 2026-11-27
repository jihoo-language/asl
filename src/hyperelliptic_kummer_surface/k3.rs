use std::fmt::Write;

use crate::ec_poly::Field;
use crate::error::{Error, Result};
use crate::finite_field_factor::distinct_degrees;
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::model::Model;
use crate::rv_ore::gcd;
use crate::rv_qp_quad::Wide;
use crate::syntax::Label;

mod artin_tate;

pub(crate) enum Query {
    Zeta,
    Points(u128),
    Picard(u128),
    ArtinTate(u128),
}

pub(crate) struct Command(pub(crate) Label, pub(crate) Query);

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<String> {
        let Self(name, query) = self;
        let curve = model
            .kummer_curve(&name.text)
            .map_err(|error| error.at(name.span))?;
        let q = curve.ring.size_of_residue_field()?;
        let orbits = node_orbits(model, &curve, q)?;
        let weil = curve.weil(model)?;
        Ok(match query {
            Query::Zeta => zeta(&name.text, q, cross_quartic(q, weil)?, &orbits),
            Query::Points(degree) => {
                let (singular, nodes, resolved) = point_counts(model, q, weil, &orbits, degree)?;
                format!("rvKummerK3Points[{}, {degree}] = singular={singular}, nodes={nodes}, resolved={resolved}", name.text)
            }
            Query::Picard(degree) => {
                let (rational, geometric) = picard_ranks(model, q, weil, &orbits, degree)?;
                format!(
                    "rvKummerK3Picard[{}, {degree}] = rank={rational}, geometric rank={geometric}",
                    name.text
                )
            }
            Query::ArtinTate(degree) => {
                artin_tate::report(model, &curve, &name.text, q, weil, &orbits, degree)?
            }
        })
    }
}

pub(crate) fn node_orbits(model: &mut Model, curve: &Curve, q: u64) -> Result<Vec<usize>> {
    let field = Field { ring: &curve.ring };
    let polynomial = curve.polynomial();
    model.spend(polynomial.len() * field.cells())?;
    let mut branch = vec![1]; // The point at infinity.
    distinct_degrees(
        model,
        field,
        Some(q),
        polynomial.to_vec(),
        |_, block, degree| {
            branch.extend(std::iter::repeat_n(degree, (block.len() - 1) / degree));
            Ok(())
        },
    )?;
    model.spend(6 * 6)?;
    let mut orbits = vec![1]; // The zero class of J[2].
    for (i, &length) in branch.iter().enumerate() {
        // Two points in the same cycle have one orbit at each cyclic distance.
        // At the antipodal distance a half-turn preserves the unordered pair.
        orbits.extend(std::iter::repeat_n(length, (length - 1) / 2));
        if length.is_multiple_of(2) {
            orbits.push(length / 2);
        }
        for &other in &branch[..i] {
            let count = gcd(length as u128, other as u128) as usize;
            orbits.extend(std::iter::repeat_n(length * other / count, count));
        }
    }
    orbits.sort_unstable();
    debug_assert_eq!(orbits.iter().sum::<usize>(), 16);
    Ok(orbits)
}

fn coefficient_overflow() -> Error {
    Error::plain("Kummer K3 abelian Frobenius coefficient exceeds signed 128-bit range")
}

pub(super) fn cross_quartic(q: u64, weil: &[i128]) -> Result<[i128; 5]> {
    // The other two eigenvalues of the exterior square are both q.
    let compute = || {
        let q = i128::from(q);
        let q2 = q.checked_mul(q)?;
        let a_squared = weil[1].checked_mul(weil[1])?;
        let b = weil[2];
        let linear = q.checked_mul(2)?.checked_sub(b)?;
        let quadratic = q.checked_mul(a_squared.checked_add(linear)?.checked_sub(b)?)?;
        Some([
            1,
            linear,
            quadratic,
            q2.checked_mul(linear)?,
            q2.checked_mul(q2)?,
        ])
    };
    compute().ok_or_else(coefficient_overflow)
}

fn zeta(name: &str, q: u64, cross: [i128; 5], orbits: &[usize]) -> String {
    let q2 = u128::from(q) * u128::from(q);
    let mut result = format!("rvKummerK3Zeta[{name}] = 1/((1-T)(1-{q2}*T)(1-{q}*T)^2(1");
    for (power, &coefficient) in cross.iter().enumerate().skip(1).filter(|(_, c)| **c != 0) {
        let sign = if coefficient < 0 { '-' } else { '+' };
        write!(result, "{sign}{}*T", coefficient.unsigned_abs()).unwrap();
        if power > 1 {
            write!(result, "^{power}").unwrap();
        }
    }
    result.push(')');
    for orbit in orbits.chunk_by(|a, b| a == b) {
        let length = orbit[0];
        if length == 1 {
            write!(result, "(1-{q}*T)").unwrap();
        } else {
            write!(result, "(1-({q}*T)^{length})").unwrap();
        }
        if orbit.len() > 1 {
            write!(result, "^{}", orbit.len()).unwrap();
        }
    }
    write!(result, "), node orbit lengths={orbits:?}").unwrap();
    result
}

fn point_overflow() -> Error {
    Error::plain("Kummer K3 point count exceeds unsigned 128-bit range")
}

pub(super) fn point_counts(
    model: &mut Model,
    q: u64,
    weil: &[i128],
    orbits: &[usize],
    degree: u128,
) -> Result<(u128, usize, u128)> {
    let exponent = u32::try_from(degree).map_err(|_| point_overflow())?;
    let field_order = u128::from(q)
        .checked_pow(exponent)
        .ok_or_else(point_overflow)?;
    let square = field_order
        .checked_mul(field_order)
        .ok_or_else(point_overflow)?;
    let trace = exterior_trace(model, weil, exponent as usize)?;
    let nodes = fixed_nodes(orbits, degree);
    let singular = square
        .checked_add_signed(1 + trace)
        .ok_or_else(point_overflow)?;
    let resolved = singular
        .checked_add(field_order * nodes as u128)
        .ok_or_else(point_overflow)?;
    Ok((singular, nodes, resolved))
}

fn exterior_trace(model: &mut Model, weil: &[i128], degree: usize) -> Result<i128> {
    Ok(powered_weil_coefficients(model, weil, degree, point_overflow)?.1)
}

pub(crate) fn powered_weil_coefficients(
    model: &mut Model,
    weil: &[i128],
    degree: usize,
    overflow: fn() -> Error,
) -> Result<(i128, i128)> {
    // Both callers bound the exponent by checking the extension size first.
    let traces = h1_power_sums::<8>(model, weil, 2 * degree, overflow)?;
    let first = traces[degree].signed().ok_or_else(overflow)?;
    let second = exterior_power_sum(model, &traces, degree, overflow)?
        .signed()
        .ok_or_else(overflow)?;
    Ok((first, second))
}

pub(super) fn fixed_nodes(orbits: &[usize], degree: u128) -> usize {
    orbits
        .iter()
        .filter(|&&d| degree.is_multiple_of(d as u128))
        .sum()
}

pub(super) fn h1_power_sums<const N: usize>(
    model: &mut Model,
    weil: &[i128],
    degree: usize,
    overflow: fn() -> Error,
) -> Result<Vec<Wide<N>>> {
    let mut traces = vec![Wide::from_i128(0)];
    for exponent in 1..=degree {
        let terms = exponent.min(weil.len() - 1);
        model.spend(terms * N)?;
        let mut trace = Wide::from_i128(0);
        for i in 1..=terms {
            let power = if i == exponent {
                Wide::from_i128(exponent as i128)
            } else {
                traces[exponent - i]
            };
            let term = power.multiply(weil[i]).map_err(|_| overflow())?;
            trace = trace.add(term.negate()).map_err(|_| overflow())?;
        }
        traces.push(trace);
    }
    Ok(traces)
}

pub(super) fn exterior_power_sum<const N: usize>(
    model: &mut Model,
    traces: &[Wide<N>],
    degree: usize,
    overflow: fn() -> Error,
) -> Result<Wide<N>> {
    let sr = traces[degree];
    let s2r = traces[2 * degree];
    let numerator = sr
        .multiply_wide(sr)
        .and_then(|square| square.add(s2r.negate()))
        .map_err(|_| overflow())?;
    let (magnitude, remainder) = numerator.magnitude.divide(model, 2)?;
    debug_assert_eq!(remainder, 0);
    Ok(Wide {
        magnitude,
        negative: numerator.negative,
    })
}

// This list contains all cyclotomic polynomials of degree at most four.
const CYCLOTOMIC: &[(u128, &[i128])] = &[
    (1, &[1, -1]),
    (2, &[1, 1]),
    (3, &[1, 1, 1]),
    (4, &[1, 0, 1]),
    (5, &[1, 1, 1, 1, 1]),
    (6, &[1, -1, 1]),
    (8, &[1, 0, 0, 0, 1]),
    (10, &[1, -1, 1, -1, 1]),
    (12, &[1, 0, -1, 0, 1]),
];

fn picard_ranks(
    model: &mut Model,
    q: u64,
    weil: &[i128],
    orbits: &[usize],
    degree: u128,
) -> Result<(u128, u8)> {
    let exceptional: u128 = orbits.iter().map(|&d| gcd(degree, d as u128)).sum();
    let (rational, geometric) = abelian_picard(model, q, cross_quartic(q, weil)?, degree)?;
    Ok((u128::from(rational) + exceptional, geometric + 16))
}

fn abelian_picard(model: &mut Model, q: u64, cross: [i128; 5], degree: u128) -> Result<(u8, u8)> {
    let mut remaining = cross.to_vec();
    let (mut rational, mut geometric) = (2, 2);
    for &(order, cyclotomic) in CYCLOTOMIC {
        if remaining.len() < cyclotomic.len() {
            continue;
        }
        let q = i128::from(q);
        let factor = cyclotomic
            .iter()
            .enumerate()
            .map(|(i, &coefficient)| {
                q.checked_pow(i as u32)
                    .and_then(|pow| pow.checked_mul(coefficient))
                    .ok_or_else(coefficient_overflow)
            })
            .collect::<Result<Vec<_>>>()?;
        while let Some(quotient) = monic_quotient(model, &remaining, &factor)? {
            let rank = (factor.len() - 1) as u8;
            geometric += rank;
            if degree.is_multiple_of(order) {
                rational += rank;
            }
            remaining = quotient;
        }
    }
    Ok((rational, geometric))
}

// Both polynomials are in descending order, and the denominator is monic.
fn monic_quotient(
    model: &mut Model,
    numerator: &[i128],
    denominator: &[i128],
) -> Result<Option<Vec<i128>>> {
    let Some(degree) = numerator.len().checked_sub(denominator.len()) else {
        return Ok(None);
    };
    model.spend((degree + 1) * denominator.len())?;
    let mut remainder = numerator.to_vec();
    for i in 0..=degree {
        let leading = remainder[i];
        for (j, &coefficient) in denominator.iter().enumerate().skip(1) {
            let term = leading
                .checked_mul(coefficient)
                .ok_or_else(coefficient_overflow)?;
            remainder[i + j] = remainder[i + j]
                .checked_sub(term)
                .ok_or_else(coefficient_overflow)?;
        }
    }
    if remainder[degree + 1..].iter().all(|&c| c == 0) {
        remainder.truncate(degree + 1);
        Ok(Some(remainder))
    } else {
        Ok(None)
    }
}
