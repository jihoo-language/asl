use std::mem::size_of;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::syntax::Label;
use crate::wide_integer::Natural;

use super::k3::{exterior_power_sum, fixed_nodes, h1_power_sums, node_orbits};

const LIMBS: usize = 16;
type Integer = Wide<LIMBS>;
type Magnitude = Natural<LIMBS>;
type Term = (usize, Integer);

pub(crate) enum Command {
    Points(Label, u128, u128),
    Trace(Label, u128, u128, u128),
    Betti(u128, u128),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<String> {
        match self {
            Self::Points(name, n, r) => {
                let curve = model
                    .kummer_curve(&name.text)
                    .map_err(|error| error.at(name.span))?;
                let value = if n == 0 {
                    1
                } else {
                    let value = coefficient(model, n, None, Some((&curve, r)))?;
                    unsigned(value).ok_or_else(|| {
                        Error::plain(
                            "Kummer Hilbert-scheme point count exceeds unsigned 128-bit range",
                        )
                    })?
                };
                Ok(format!(
                    "rvKummerHilbPoints[{}, {n}, {r}] = {value}",
                    name.text
                ))
            }
            Self::Trace(name, n, r, i) => {
                let curve = model
                    .kummer_curve(&name.text)
                    .map_err(|error| error.at(name.span))?;
                let value = trace(model, &curve, n, r, i)?;
                Ok(format!(
                    "rvKummerHilbTrace[{}, {n}, {r}, {i}] = {value}",
                    name.text
                ))
            }
            Self::Betti(n, i) => {
                let value = match graded_index(n, i) {
                    None => 0,
                    Some((_, 0, _)) => 1,
                    Some((length, degree, _)) => {
                        let value = coefficient(model, length, Some(degree), None)?;
                        unsigned(value).ok_or_else(|| {
                            Error::plain(
                                "K3 Hilbert-scheme Betti number exceeds unsigned 128-bit range",
                            )
                        })?
                    }
                };
                Ok(format!("rvKummerHilbBetti[{n}, {i}] = {value}"))
            }
        }
    }
}

fn trace(model: &mut Model, curve: &Curve, n: u128, r: u128, i: u128) -> Result<i128> {
    let Some((length, degree, twist)) = graded_index(n, i) else {
        return Ok(0);
    };
    let mut value = if degree == 0 {
        Integer::from_i128(1)
    } else {
        coefficient(model, length, Some(degree), Some((curve, r)))?
    };
    if twist != 0 && !value.is_zero() {
        let exponent = r.checked_mul(twist).ok_or_else(intermediate_overflow)?;
        let q = curve.ring.size_of_residue_field()?;
        let factor = power(model, q, exponent)?;
        value.magnitude = multiply(value.magnitude, factor)?;
    }
    value.signed().ok_or_else(|| {
        Error::plain("Kummer Hilbert-scheme Frobenius trace exceeds signed 128-bit range")
    })
}

// Return the truncated length, half-degree, and Poincaré-duality Tate twist.
fn graded_index(n: u128, i: u128) -> Option<(u128, u128, u128)> {
    let half = i / 2;
    if !i.is_multiple_of(2) || half.div_ceil(2) > n {
        return None;
    }
    // After removing (1-t)^(-1) from H^0(S), every generator of length m
    // has half-degree at least m/2. In half-degree d the coefficients thus
    // stabilize once the length reaches 2d.
    let degree = if half <= n { half } else { n - (half - n) };
    Some((n.min(degree.saturating_mul(2)), degree, half - degree))
}

fn coefficient(
    model: &mut Model,
    length: u128,
    degree: Option<u128>,
    arithmetic: Option<(&Curve, u128)>,
) -> Result<Integer> {
    let length = usize::try_from(length).map_err(|_| storage_limit())?;
    let degree = degree
        .map(usize::try_from)
        .transpose()
        .map_err(|_| storage_limit())?;
    let top = degree.unwrap_or(0);
    check_storage(model, length, top)?;
    let limit = length.min(degree.unwrap_or(length));
    let surface = match arithmetic {
        None => Surface::Betti,
        Some((curve, r)) => Surface::arithmetic(model, curve, r, limit)?,
    };
    let weights = log_coefficients(model, length, degree, &surface)?;
    let zero = Integer::from_i128(0);
    let mut first = vec![zero; top + 1];
    first[0] = Integer::from_i128(1);
    let mut coefficients = Vec::with_capacity(length + 1);
    coefficients.push(first);
    for n in 1..=length {
        let mut row = vec![zero; top + 1];
        for (j, terms) in weights.iter().enumerate().take(n + 1).skip(1) {
            let previous = &coefficients[n - j];
            for &(shift, weight) in terms {
                let cells = top + 1 - shift;
                model.spend(cells)?;
                for (target, &old) in row[shift..].iter_mut().zip(previous) {
                    if old.is_zero() {
                        continue;
                    }
                    let term = checked_product(model, weight, old)?;
                    *target = target.add(term).map_err(|_| intermediate_overflow())?;
                }
            }
        }
        for cell in &mut row {
            if !cell.is_zero() {
                let (magnitude, remainder) = cell.magnitude.divide(model, n as u128)?;
                if remainder != 0 {
                    return Err(Error::plain(
                        "Kummer Hilbert-scheme Newton division is inexact",
                    ));
                }
                cell.magnitude = magnitude;
            }
        }
        coefficients.push(row);
    }
    Ok(coefficients[length][top])
}

fn log_coefficients(
    model: &mut Model,
    length: usize,
    degree: Option<usize>,
    surface: &Surface,
) -> Result<Vec<Vec<Term>>> {
    let mut result = vec![Vec::new(); length + 1];
    for m in 1..=length {
        for l in 1..=length / m {
            let j = m * l;
            let twist = j - l;
            if degree.is_some_and(|top| twist > top) {
                continue;
            }
            let magnitude = surface
                .power(twist)
                .multiply_u128(m as u128)
                .ok_or_else(intermediate_overflow)?;
            let factor = Integer {
                magnitude,
                negative: false,
            };
            for b in 0..=2 {
                let grade = twist + l * b;
                if degree.is_some_and(|top| grade > top) {
                    break;
                }
                let shift = degree.map_or(0, |_| grade);
                let trace = surface.trace(b, l)?;
                if !trace.is_zero() {
                    result[j].push((shift, checked_product(model, factor, trace)?));
                }
            }
        }
    }
    Ok(result)
}

enum Surface {
    Betti,
    Frobenius {
        powers: Vec<Magnitude>,
        middle: Vec<Integer>,
    },
}

impl Surface {
    fn arithmetic(model: &mut Model, curve: &Curve, r: u128, limit: usize) -> Result<Self> {
        let q = curve.ring.size_of_residue_field()?;
        let q_power = power(model, q, r)?;
        let mut powers = Vec::with_capacity(limit + 1);
        powers.push(Magnitude::from_small(1));
        for l in 1..=limit {
            model.spend(LIMBS)?;
            powers.push(multiply(powers[l - 1], q_power)?);
        }
        // The bounds on the field powers also bound the required H¹ exponents.
        let step = usize::try_from(r).map_err(|_| intermediate_overflow())?;
        let orbits = node_orbits(model, curve, q)?;
        let weil = curve.weil(model)?;
        let traces = h1_power_sums(model, weil, 2 * step * limit, intermediate_overflow)?;
        let mut middle = Vec::with_capacity(limit + 1);
        middle.push(Integer::from_i128(0));
        for (l, &field_order) in powers.iter().enumerate().skip(1) {
            let exterior = exterior_power_sum(model, &traces, step * l, intermediate_overflow)?;
            let nodes = fixed_nodes(&orbits, r * l as u128);
            let magnitude = field_order
                .multiply_u128(nodes as u128)
                .ok_or_else(intermediate_overflow)?;
            let exceptional = Integer {
                magnitude,
                negative: false,
            };
            middle.push(
                exterior
                    .add(exceptional)
                    .map_err(|_| intermediate_overflow())?,
            );
        }
        Ok(Self::Frobenius { powers, middle })
    }

    fn power(&self, exponent: usize) -> Magnitude {
        match self {
            Self::Betti => Magnitude::from_small(1),
            Self::Frobenius { powers, .. } => powers[exponent],
        }
    }

    fn trace(&self, half_degree: usize, l: usize) -> Result<Integer> {
        Ok(match (self, half_degree) {
            (_, 0) | (Self::Betti, 2) => Integer::from_i128(1),
            (Self::Betti, 1) => Integer::from_i128(22),
            (Self::Frobenius { middle, .. }, 1) => middle[l],
            (Self::Frobenius { powers, .. }, 2) => Integer {
                magnitude: multiply(powers[l], powers[l])?,
                negative: false,
            },
            _ => unreachable!("only even surface degrees 0, 2, 4 are used"),
        })
    }
}

fn unsigned(value: Integer) -> Option<u128> {
    if value.negative && !value.is_zero() {
        return None;
    }
    value.magnitude.small()
}

fn multiply(left: Magnitude, right: Magnitude) -> Result<Magnitude> {
    left.checked_product(right)
        .ok_or_else(intermediate_overflow)
}

fn checked_product(model: &mut Model, left: Integer, right: Integer) -> Result<Integer> {
    model.spend(1 + left.magnitude.bits().div_ceil(32) * right.magnitude.bits().div_ceil(32))?;
    left.multiply_wide(right)
        .map_err(|_| intermediate_overflow())
}

fn power(model: &mut Model, q: u64, mut exponent: u128) -> Result<Magnitude> {
    let mut result = Magnitude::from_small(1);
    let mut factor = Magnitude::from_small(u128::from(q));
    while exponent != 0 {
        model.spend(LIMBS)?;
        if exponent & 1 != 0 {
            result = multiply(result, factor)?;
        }
        exponent >>= 1;
        if exponent != 0 {
            factor = multiply(factor, factor)?;
        }
    }
    Ok(result)
}

fn check_storage(model: &mut Model, length: usize, degree: usize) -> Result<()> {
    let entries = length
        .saturating_add(1)
        .saturating_mul(degree.saturating_add(1));
    let words = size_of::<Integer>().div_ceil(size_of::<usize>());
    let rows = length.saturating_add(1).saturating_mul(6);
    let cells = entries.saturating_mul(words).saturating_add(rows);
    if cells > MAX_PIVOT_CELLS {
        return Err(storage_limit());
    }
    model.spend(length)?;
    let terms: usize = (1..=length).map(|m| length / m).sum();
    let words = size_of::<Term>().div_ceil(size_of::<usize>());
    let cells = cells + 3 * terms * words;
    if cells > MAX_PIVOT_CELLS {
        return Err(storage_limit());
    }
    model.spend(cells)
}

fn storage_limit() -> Error {
    Error::plain("Kummer Hilbert-scheme calculation exceeds the coefficient storage limit")
}

fn intermediate_overflow() -> Error {
    Error::plain("Kummer Hilbert-scheme exact intermediate exceeds 512 bits")
}
