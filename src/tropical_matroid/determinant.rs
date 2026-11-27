// The fallback works in Z[sqrt(d)] after clearing each row's denominators.
// Its exact intermediate integer coefficients have a separate scratch bound.
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;
use crate::rv_ore::gcd;
use crate::rv_qp_quad::{self, gauge::Coefficients, Element, Field, Wide};
use crate::wide_integer::Natural;

const LIMBS: usize = 64;
type Integer<const N: usize = LIMBS> = Wide<N>;

#[derive(Clone, Copy)]
struct Entry<const N: usize = LIMBS> {
    a: Integer<N>,
    b: Integer<N>,
}

impl<const N: usize> Entry<N> {
    fn integer(a: i128) -> Self {
        Self {
            a: Integer::from_i128(a),
            b: Integer::from_i128(0),
        }
    }

    fn is_zero(self) -> bool {
        self.a.magnitude.is_zero() && self.b.magnitude.is_zero()
    }

    fn resize<const M: usize>(self) -> Result<Entry<M>> {
        Ok(Entry {
            a: self.a.resize()?,
            b: self.b.resize()?,
        })
    }

    fn same_value(self, other: Self) -> bool {
        let same = |a: Integer<N>, b: Integer<N>| {
            a.magnitude == b.magnitude && (a.is_zero() || a.negative == b.negative)
        };
        same(self.a, other.a) && same(self.b, other.b)
    }

    fn multiply(self, model: &mut Model, other: Self, d: Integer<N>) -> Result<Self> {
        let mut a = multiply(model, self.a, other.a)?;
        if self.b.magnitude.is_zero() && other.b.magnitude.is_zero() {
            return Ok(Self { a, b: self.b });
        }
        let inner = multiply(model, self.b, other.b)?;
        a = a.add(multiply(model, inner, d)?)?;
        let b = multiply(model, self.a, other.b)?.add(multiply(model, self.b, other.a)?)?;
        Ok(Self { a, b })
    }

    fn subtract(self, other: Self) -> Result<Self> {
        Ok(Self {
            a: self.a.add(other.a.negate())?,
            b: self.b.add(other.b.negate())?,
        })
    }

    fn divide_exact(self, model: &mut Model, divisor: Self, d: Integer<N>) -> Result<Self> {
        if self.is_zero() {
            return Ok(self);
        }
        let (numerator, divisor) = if divisor.b.magnitude.is_zero() {
            (self, divisor.a)
        } else {
            let conjugate = Self {
                a: divisor.a,
                b: divisor.b.negate(),
            };
            let norm = divisor.multiply(model, conjugate, d)?;
            debug_assert!(norm.b.magnitude.is_zero());
            (self.multiply(model, conjugate, d)?, norm.a)
        };
        let mut quotient = |value: Integer<N>| -> Result<Integer<N>> {
            value.divide_exact(model, divisor)?.ok_or_else(|| {
                Error::plain("fraction-free tropical elimination encountered an inexact division")
            })
        };
        Ok(Self {
            a: quotient(numerator.a)?,
            b: quotient(numerator.b)?,
        })
    }
}

fn multiply<const N: usize>(
    model: &mut Model,
    lhs: Integer<N>,
    rhs: Integer<N>,
) -> Result<Integer<N>> {
    model.spend(1 + lhs.magnitude.bits().div_ceil(32) * rhs.magnitude.bits().div_ceil(32))?;
    lhs.multiply_wide(rhs)
}

fn scratch(model: &mut Model, rows: usize, width: usize, diagnostic: &str) -> Result<()> {
    let words_per_entry = std::mem::size_of::<Entry>().div_ceil(std::mem::size_of::<usize>());
    let cells = rows * (width * words_per_entry + 3);
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(diagnostic));
    }
    model.spend(cells)
}

fn bareiss_entry<const N: usize>(
    model: &mut Model,
    entry: Entry<N>,
    upper: Entry<N>,
    pivot: Entry<N>,
    factor: Entry<N>,
    previous: Entry<N>,
    d: Integer<N>,
) -> Result<Entry<N>> {
    let numerator = match (!entry.is_zero(), !factor.is_zero() && !upper.is_zero()) {
        (false, false) => return Ok(entry),
        (true, false) => pivot.multiply(model, entry, d)?,
        (false, true) => {
            let cross = factor.multiply(model, upper, d)?;
            Entry {
                a: cross.a.negate(),
                b: cross.b.negate(),
            }
        }
        (true, true) => {
            model.spend(N)?;
            let product = pivot.multiply(model, entry, d)?;
            let cross = factor.multiply(model, upper, d)?;
            product.subtract(cross)?
        }
    };
    // Bareiss's identity makes the quotient exact in the input ring.
    numerator.divide_exact(model, previous, d)
}

fn reduce_pair(model: &mut Model, d: i128, mut matrix: Vec<Vec<Entry>>) -> Result<Option<Entry>> {
    let rank = matrix.len();
    let d = Integer::from_i128(d);
    let mut previous = Entry::integer(1);
    for k in 0..rank {
        model.spend(rank - k)?;
        let Some(pivot_row) = (k..rank).find(|&row| !matrix[row][k].is_zero()) else {
            return Ok(None);
        };
        matrix.swap(k, pivot_row);
        let (head, tail) = matrix.split_at_mut(k + 1);
        let pivot_row = &head[k];
        let pivot = pivot_row[k];
        model.spend(LIMBS)?;
        let identity_scale = pivot.same_value(previous);
        for row in tail {
            model.spend(1)?;
            let factor = row[k];
            let has_factor = !factor.is_zero();
            if identity_scale && !has_factor {
                // The Bareiss numerator scales this row by exactly the divisor.
                continue;
            }
            model.spend(rank - k - 1)?;
            for (entry, upper) in row[k + 1..].iter_mut().zip(&pivot_row[k + 1..]) {
                *entry = bareiss_entry(model, *entry, *upper, pivot, factor, previous, d)?;
            }
        }
        // Row exchanges change the sign only, which is immaterial to valuation.
        previous = pivot;
    }
    Ok(Some(previous))
}

pub(super) fn rational(
    model: &mut Model,
    matrix: Vec<Vec<i128>>,
    prime: Prime,
) -> Result<Option<i64>> {
    scratch(
        model,
        matrix.len(),
        matrix.len(),
        "wide tropical minor exceeds the determinant storage limit",
    )?;
    let matrix = matrix
        .into_iter()
        .map(|row| row.into_iter().map(Entry::integer).collect())
        .collect();
    let Some(det) = reduce_pair(model, 0, matrix)? else {
        return Ok(None);
    };
    Ok(Some(det.a.prime_part(model, prime)?.0))
}

fn integral_row(model: &mut Model, row: Vec<Element>, prime: Prime) -> Result<(Vec<Entry>, i64)> {
    let mut common = Natural::<LIMBS>::from_small(1);
    let mut order = 0;
    for entry in &row {
        if entry.c != 1 {
            model.spend(4 * LIMBS)?;
            let denominator = entry.c as u128;
            let (_, remainder) = common.divide(model, denominator)?;
            let factor = denominator / gcd(denominator, remainder);
            common = common.multiply_u128(factor).ok_or_else(|| {
                Error::plain("common denominator exceeds the wide tropical determinant bound")
            })?;
            order += super::valuation(model, factor, u128::from(prime.get()))?;
        }
    }
    let row = row
        .into_iter()
        .map(|entry| {
            model.spend(8 * LIMBS)?;
            let (factor, residue) = common.divide(model, entry.c as u128)?;
            debug_assert_eq!(residue, 0);
            let coefficient = |x: i128| -> Result<Integer> {
                let magnitude = factor.multiply_u128(x.unsigned_abs()).ok_or_else(|| {
                    Error::plain("cleared coefficient exceeds the wide tropical determinant bound")
                })?;
                Ok(Integer {
                    magnitude,
                    negative: x < 0,
                })
            };
            Ok(Entry {
                a: coefficient(entry.a)?,
                b: coefficient(entry.b)?,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    Ok((row, order))
}

pub(super) fn quadratic(
    model: &mut Model,
    field: Field,
    matrix: Vec<Vec<Element>>,
) -> Result<Option<i64>> {
    scratch(
        model,
        matrix.len(),
        matrix.len(),
        "wide tropical minor exceeds the determinant storage limit",
    )?;
    let prime = field.prime();
    let mut denominator_order = 0;
    let matrix = matrix
        .into_iter()
        .map(|row| {
            let (row, order) = integral_row(model, row, prime)?;
            denominator_order += order;
            Ok(row)
        })
        .collect::<Result<Vec<_>>>()?;
    let Some(det) = reduce_pair(model, field.radicand(), matrix)? else {
        return Ok(None);
    };
    let numerator_order = field
        .wide_local(model, det.a, det.b)?
        .expect("nonzero determinant")
        .0;
    Ok(Some(numerator_order - denominator_order))
}

pub(super) fn explain_limit(error: Error) -> Error {
    if rv_qp_quad::is_coefficient_overflow(&error) {
        Error::plain(format!(
            "wide tropical determinant exceeds its {}-bit intermediate arithmetic bound",
            32 * LIMBS
        ))
    } else {
        error
    }
}

// A product difference of two pairs of retained coefficients can use
// 2*LIMBS+5 limbs. Division by a quadratic pivot may then multiply
// it by the pivot's conjugate before the exact quotient is available.
const DIVISION_LIMBS: usize = 3 * LIMBS + 9;
const NORM_LIMBS: usize = 2 * LIMBS + 5;

fn fraction_free_row(
    model: &mut Model,
    row: &mut [Entry],
    upper: &[Entry],
    j: usize,
    previous: Entry,
    d: Integer,
) -> Result<()> {
    let pivot = upper[j];
    let factor = row[j];
    model.spend(LIMBS)?;
    if factor.is_zero() && pivot.same_value(previous) {
        return Ok(());
    }
    row[j] = Entry::integer(0);
    model.spend(row.len())?;
    for (column, (entry, &upper)) in row.iter_mut().zip(upper).enumerate() {
        if column == j {
            continue;
        }
        let result = bareiss_entry(model, *entry, upper, pivot, factor, previous, d);
        *entry = match result {
            Err(error) if rv_qp_quad::is_coefficient_overflow(&error) => {
                model.spend(7 * DIVISION_LIMBS)?;
                bareiss_entry::<DIVISION_LIMBS>(
                    model,
                    entry.resize()?,
                    upper.resize()?,
                    pivot.resize()?,
                    factor.resize()?,
                    previous.resize()?,
                    d.resize()?,
                )?
                .resize()?
            }
            other => other?,
        };
    }
    Ok(())
}

fn local(model: &mut Model, field: Coefficients, entry: Entry) -> Result<Option<(i64, u32)>> {
    match field {
        Coefficients::Rational(prime) => {
            if entry.a.is_zero() {
                Ok(None)
            } else {
                Ok(Some(entry.a.prime_part(model, prime)?))
            }
        }
        Coefficients::Quadratic(field) => match field.wide_local(model, entry.a, entry.b) {
            Err(error) if rv_qp_quad::is_coefficient_overflow(&error) => {
                model.spend(2 * NORM_LIMBS)?;
                field.wide_local::<NORM_LIMBS>(model, entry.a.resize()?, entry.b.resize()?)
            }
            other => other,
        },
    }
}

fn reduce_weighted_row(
    model: &mut Model,
    row: &[Entry],
    weights: &[i64],
    field: Coefficients,
) -> Result<Option<(usize, Vec<u32>)>> {
    let mut minimum = None;
    let mut locals = Vec::with_capacity(row.len());
    for (j, (&entry, &weight)) in row.iter().zip(weights).enumerate() {
        model.spend(1)?;
        let value = local(model, field, entry)?
            .map(|(order, unit)| (i128::from(order) - i128::from(weight), unit));
        if let Some((order, unit)) = value {
            if minimum.is_none_or(|(best, _, _)| order < best) {
                minimum = Some((order, j, unit));
            }
        }
        locals.push(value);
    }
    let Some((minimum, j, unit)) = minimum else {
        return Ok(None);
    };
    let prime = field.prime();
    let reciprocal = invert(model, prime, unit)?;
    model.spend(locals.len())?;
    let row = locals
        .into_iter()
        .map(|value| match value {
            Some((order, unit)) if order == minimum => prime.mul(unit, reciprocal),
            _ => 0,
        })
        .collect();
    Ok(Some((j, row)))
}

pub(crate) fn weighted_echelon(
    model: &mut Model,
    field: Coefficients,
    rows: Vec<Vec<Element>>,
    weights: &[i64],
) -> Result<Vec<Vec<u32>>> {
    let width = weights.len();
    scratch(
        model,
        width.min(rows.len()),
        width,
        "wide residue echelon exceeds the storage limit",
    )?;
    let d = Integer::from_i128(match field {
        Coefficients::Rational(_) => 0,
        Coefficients::Quadratic(field) => field.radicand(),
    });
    let mut echelon: Vec<(usize, Vec<Entry>)> = Vec::new();
    let mut reduced = Vec::new();
    for row in rows {
        model.spend(1 + row.len())?;
        if row.iter().all(|entry| entry.is_zero()) {
            continue;
        }
        // Clearing each row independently preserves the exact row space.
        let (mut row, _) = integral_row(model, row, field.prime())?;
        let mut previous = Entry::integer(1);
        for (j, upper) in &echelon {
            fraction_free_row(model, &mut row, upper, *j, previous, d)?;
            previous = upper[*j];
        }
        if let Some((j, residue)) = reduce_weighted_row(model, &row, weights, field)? {
            // Scaling each exact weighted row to make this pivot one
            // gives integral rows. Their pivot coordinates are triangular
            // with unit diagonal, so they span the saturated lattice.
            echelon.push((j, row));
            reduced.push(residue);
            if echelon.len() == width {
                break;
            }
        }
    }
    Ok(reduced)
}

pub(crate) fn explain_residue_limit(error: Error) -> Error {
    if rv_qp_quad::is_coefficient_overflow(&error) {
        Error::plain(format!(
            "wide residue echelon exceeds its {}-bit retained coefficient bound",
            32 * LIMBS
        ))
    } else {
        error
    }
}
