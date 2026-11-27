use std::mem::size_of;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::syntax::Label;
use crate::wide_integer::Natural;

use super::k3::{exterior_power_sum, fixed_nodes, h1_power_sums, node_orbits};

pub(crate) mod characteristic;
pub(crate) mod elliptic_genus;

const LIMBS: usize = 16;
pub(crate) type Integer = Wide<LIMBS>;
pub(crate) type Magnitude = Natural<LIMBS>;
type Term = (usize, Integer);

#[derive(Clone, Copy)]
pub(crate) enum SurfaceKind {
    K3,
    Abelian,
}

impl SurfaceKind {
    fn grading_unit(self) -> usize {
        match self {
            Self::K3 => 2,
            Self::Abelian => 1,
        }
    }
}

pub(crate) enum Command {
    Points(Label, u128, u128),
    Trace(Label, u128, u128, u128),
    Betti(u128, u128),
    Characteristic(characteristic::Query),
    EllipticGenus(elliptic_genus::Query),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<String> {
        match self {
            Self::Characteristic(query) => query.report(model),
            Self::EllipticGenus(query) => query.report(model),
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
                let value = betti(model, n, i)?;
                Ok(format!("rvKummerHilbBetti[{n}, {i}] = {value}"))
            }
        }
    }
}

pub(crate) fn betti(model: &mut Model, n: u128, i: u128) -> Result<u128> {
    match graded_index(n, i) {
        None => Ok(0),
        Some((_, 0, _)) => Ok(1),
        Some((length, degree, _)) => {
            let value = coefficient(model, length, Some(degree), None)?;
            unsigned(value).ok_or_else(|| {
                Error::plain("K3 Hilbert-scheme Betti number exceeds unsigned 128-bit range")
            })
        }
    }
}

fn trace(model: &mut Model, curve: &Curve, n: u128, r: u128, i: u128) -> Result<i128> {
    graded_trace(model, curve, n, r, i, |model, length, degree| {
        coefficient(model, length, Some(degree), Some((curve, r)))
    })?
    .signed()
    .ok_or_else(|| {
        Error::plain("Kummer Hilbert-scheme Frobenius trace exceeds signed 128-bit range")
    })
}

pub(super) fn graded_trace(
    model: &mut Model,
    curve: &Curve,
    n: u128,
    r: u128,
    i: u128,
    coefficient: impl FnOnce(&mut Model, u128, u128) -> Result<Integer>,
) -> Result<Integer> {
    let Some((length, degree, twist)) = graded_index(n, i) else {
        return Ok(Integer::from_i128(0));
    };
    let mut value = if degree == 0 {
        Integer::from_i128(1)
    } else {
        coefficient(model, length, degree)?
    };
    if twist != 0 && !value.is_zero() {
        let exponent = r.checked_mul(twist).ok_or_else(intermediate_overflow)?;
        let q = curve.ring.size_of_residue_field()?;
        let factor = power(model, q, exponent)?;
        value.magnitude = multiply(value.magnitude, factor)?;
    }
    Ok(value)
}

// Return the truncated length, half-degree, and Poincaré-duality Tate twist.
pub(super) fn graded_index(n: u128, i: u128) -> Option<(u128, u128, u128)> {
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
    let coefficients = coefficients(model, length, degree, arithmetic, SurfaceKind::K3)?;
    Ok(*coefficients.last().unwrap().last().unwrap())
}

pub(crate) fn coefficients(
    model: &mut Model,
    length: u128,
    degree: Option<u128>,
    arithmetic: Option<(&Curve, u128)>,
    kind: SurfaceKind,
) -> Result<Vec<Vec<Integer>>> {
    coefficients_with_options(model, length, degree, arithmetic, kind, false, usize::MAX)
}

pub(crate) fn inversion_lefschetz(
    model: &mut Model,
    length: u128,
    arithmetic: Option<(&Curve, u128)>,
) -> Result<Vec<Vec<Integer>>> {
    let inversion = arithmetic.is_none_or(|(_, r)| !r.is_multiple_of(2));
    coefficients_with_options(
        model,
        length,
        None,
        arithmetic,
        SurfaceKind::Abelian,
        inversion,
        usize::MAX,
    )
}

pub(crate) fn abelian_symmetric_coefficients(
    model: &mut Model,
    length: u128,
    degree: u128,
    curve: &Curve,
    r: u128,
) -> Result<Vec<Vec<Integer>>> {
    coefficients_with_options(
        model,
        length,
        Some(degree),
        Some((curve, r)),
        SurfaceKind::Abelian,
        false,
        1,
    )
}

pub(super) fn k3_symmetric_coefficients(
    model: &mut Model,
    length: u128,
    degree: u128,
    traces: Option<(&[Magnitude], &[Integer], usize)>,
) -> Result<Vec<Vec<Integer>>> {
    coefficients_from_surface(
        model,
        length,
        Some(degree),
        SurfaceKind::K3,
        false,
        1,
        |model, limit| {
            let Some((fields, middle, step)) = traces else {
                return Ok(Surface::Betti(SurfaceKind::K3));
            };
            model.spend(2 * (limit + 1))?;
            Ok(Surface::Frobenius {
                powers: (0..=limit).map(|i| fields[i * step]).collect(),
                middle: (0..=limit).map(|i| middle[i * step]).collect(),
                first: Vec::new(),
            })
        },
    )
}

pub(crate) fn symmetric_partition_product(
    model: &mut Model,
    parts: &[(u128, u128)],
    rows: &[Vec<Integer>],
) -> Result<Vec<Integer>> {
    let limit = (rows.len() - 1) as u128;
    let width = rows[0].len();
    let mut product = vec![Integer::from_i128(1)];
    for &(_, count) in parts {
        let row = &rows[count.min(limit) as usize];
        product = graded_product(model, &product, row, width)?;
    }
    product.resize(width, Integer::from_i128(0));
    Ok(product)
}

pub(super) fn graded_product(
    model: &mut Model,
    left: &[Integer],
    right: &[Integer],
    width: usize,
) -> Result<Vec<Integer>> {
    if left.is_empty() || right.is_empty() {
        return Ok(Vec::new());
    }
    let width = width.min(left.len() + right.len() - 1);
    model.spend(width)?;
    let mut product = vec![Integer::from_i128(0); width];
    for (j, &left) in left
        .iter()
        .take(width)
        .enumerate()
        .filter(|(_, value)| !value.is_zero())
    {
        model.spend((width - j).min(right.len()))?;
        for (target, &right) in product[j..].iter_mut().zip(right) {
            if !right.is_zero() {
                let term = checked_product(model, left, right)?;
                *target = target.add(term).map_err(|_| intermediate_overflow())?;
            }
        }
    }
    Ok(product)
}

fn coefficients_with_options(
    model: &mut Model,
    length: u128,
    degree: Option<u128>,
    arithmetic: Option<(&Curve, u128)>,
    kind: SurfaceKind,
    inversion: bool,
    last_mode: usize,
) -> Result<Vec<Vec<Integer>>> {
    coefficients_from_surface(
        model,
        length,
        degree,
        kind,
        inversion,
        last_mode,
        |model, limit| match arithmetic {
            None => Ok(Surface::Betti(kind)),
            Some((curve, r)) => Surface::arithmetic(model, curve, r, limit, kind),
        },
    )
}

pub(super) fn semilinear_coefficient(
    model: &mut Model,
    curve: &Curve,
    length: u128,
    degree: Option<u128>,
    r: u128,
    traces: impl FnOnce(&mut Model, &[Magnitude]) -> Result<Vec<Integer>>,
) -> Result<Integer> {
    let rows = coefficients_from_surface(
        model,
        length,
        degree,
        SurfaceKind::K3,
        false,
        usize::MAX,
        |model, limit| {
            let q = curve.ring.size_of_residue_field()?;
            let powers = field_powers(model, q, r, limit)?;
            let middle = traces(model, &powers)?;
            Ok(Surface::Frobenius {
                powers,
                middle,
                first: Vec::new(),
            })
        },
    )?;
    Ok(*rows.last().unwrap().last().unwrap())
}

fn coefficients_from_surface(
    model: &mut Model,
    length: u128,
    degree: Option<u128>,
    kind: SurfaceKind,
    inversion: bool,
    last_mode: usize,
    prepare: impl FnOnce(&mut Model, usize) -> Result<Surface>,
) -> Result<Vec<Vec<Integer>>> {
    let length = usize::try_from(length).map_err(|_| storage_limit())?;
    let degree = degree
        .map(usize::try_from)
        .transpose()
        .map_err(|_| storage_limit())?;
    let top = degree.unwrap_or(0);
    check_storage(model, length, top, 1 + 4 / kind.grading_unit(), last_mode)?;
    let limit = length.min(degree.unwrap_or(length));
    let surface = prepare(model, limit)?;
    let weights = log_coefficients(model, length, degree, &surface, kind, inversion, last_mode)?;
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
    Ok(coefficients)
}

fn log_coefficients(
    model: &mut Model,
    length: usize,
    degree: Option<usize>,
    surface: &Surface,
    kind: SurfaceKind,
    inversion: bool,
    last_mode: usize,
) -> Result<Vec<Vec<Term>>> {
    let mut result = vec![Vec::new(); length + 1];
    let unit = kind.grading_unit();
    for m in 1..=length.min(last_mode) {
        for l in 1..=length / m {
            let j = m * l;
            let twist = j - l;
            let shift = (2 / unit) * twist;
            if degree.is_some_and(|top| shift > top) {
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
            for b in 0..=4 / unit {
                let grade = shift + l * b;
                if degree.is_some_and(|top| grade > top) {
                    break;
                }
                let shift = degree.map_or(0, |_| grade);
                // Odd oscillators give exterior powers. Their logarithmic
                // coefficient has a minus sign at every power l.
                let trace = surface.signed_trace(model, b * unit, l)?;
                // Inversion multiplies odd surface cohomology by −1.
                // At the l-th ghost power this contributes (−1)^l.
                let trace = if inversion && !(b * unit).is_multiple_of(2) && !l.is_multiple_of(2) {
                    trace.negate()
                } else {
                    trace
                };
                if !trace.is_zero() {
                    result[j].push((shift, checked_product(model, factor, trace)?));
                }
            }
        }
    }
    Ok(result)
}

enum Surface {
    Betti(SurfaceKind),
    Frobenius {
        powers: Vec<Magnitude>,
        middle: Vec<Integer>,
        first: Vec<Integer>,
    },
}

impl Surface {
    fn arithmetic(
        model: &mut Model,
        curve: &Curve,
        r: u128,
        limit: usize,
        kind: SurfaceKind,
    ) -> Result<Self> {
        let q = curve.ring.size_of_residue_field()?;
        let powers = field_powers(model, q, r, limit)?;
        // The bounds on the field powers also bound the required H¹ exponents.
        let step = usize::try_from(r).map_err(|_| intermediate_overflow())?;
        let orbits = match kind {
            SurfaceKind::K3 => node_orbits(model, curve, q)?,
            SurfaceKind::Abelian => Vec::new(),
        };
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
        let first = match kind {
            SurfaceKind::Abelian => (0..=limit).map(|l| traces[step * l]).collect(),
            SurfaceKind::K3 => Vec::new(),
        };
        Ok(Self::Frobenius {
            powers,
            middle,
            first,
        })
    }

    fn power(&self, exponent: usize) -> Magnitude {
        match self {
            Self::Betti(_) => Magnitude::from_small(1),
            Self::Frobenius { powers, .. } => powers[exponent],
        }
    }

    fn signed_trace(&self, model: &mut Model, degree: usize, l: usize) -> Result<Integer> {
        let value = match (self, degree) {
            (Self::Betti(SurfaceKind::K3), b) => Integer::from_i128([1, 0, 22, 0, 1][b]),
            (Self::Betti(SurfaceKind::Abelian), b) => Integer::from_i128([1, 4, 6, 4, 1][b]),
            (_, 0) => Integer::from_i128(1),
            (Self::Frobenius { first, .. }, 1) => first[l],
            (Self::Frobenius { middle, .. }, 2) => middle[l],
            (Self::Frobenius { first, powers, .. }, 3) => checked_product(
                model,
                first[l],
                Integer {
                    magnitude: powers[l],
                    negative: false,
                },
            )?,
            (Self::Frobenius { powers, .. }, 4) => Integer {
                magnitude: multiply(powers[l], powers[l])?,
                negative: false,
            },
            _ => unreachable!("a surface has degrees zero through four"),
        };
        Ok(if degree.is_multiple_of(2) {
            value
        } else {
            value.negate()
        })
    }
}

pub(super) fn field_powers(
    model: &mut Model,
    q: u64,
    r: u128,
    limit: usize,
) -> Result<Vec<Magnitude>> {
    let q_power = power(model, q, r)?;
    let mut powers = Vec::with_capacity(limit + 1);
    powers.push(Magnitude::from_small(1));
    for l in 1..=limit {
        model.spend(LIMBS)?;
        powers.push(multiply(powers[l - 1], q_power)?);
    }
    Ok(powers)
}

pub(crate) fn unsigned(value: Integer) -> Option<u128> {
    if value.negative && !value.is_zero() {
        return None;
    }
    value.magnitude.small()
}

fn multiply(left: Magnitude, right: Magnitude) -> Result<Magnitude> {
    left.checked_product(right)
        .ok_or_else(intermediate_overflow)
}

pub(crate) fn checked_product(model: &mut Model, left: Integer, right: Integer) -> Result<Integer> {
    model.spend(1 + left.magnitude.bits().div_ceil(32) * right.magnitude.bits().div_ceil(32))?;
    left.multiply_wide(right)
        .map_err(|_| intermediate_overflow())
}

pub(crate) fn power(model: &mut Model, q: u64, mut exponent: u128) -> Result<Magnitude> {
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

fn check_storage(
    model: &mut Model,
    length: usize,
    degree: usize,
    atoms: usize,
    last_mode: usize,
) -> Result<()> {
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
    let terms: usize = (1..=length.min(last_mode)).map(|m| length / m).sum();
    let words = size_of::<Term>().div_ceil(size_of::<usize>());
    let cells = cells + atoms * terms * words;
    if cells > MAX_PIVOT_CELLS {
        return Err(storage_limit());
    }
    model.spend(cells)
}

fn storage_limit() -> Error {
    Error::plain("Kummer Hilbert-scheme calculation exceeds the coefficient storage limit")
}

pub(crate) fn intermediate_overflow() -> Error {
    Error::plain("Kummer Hilbert-scheme exact intermediate exceeds 512 bits")
}
