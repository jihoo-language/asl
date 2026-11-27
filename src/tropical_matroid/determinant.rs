// The fallback works in Z[sqrt(d)] after clearing each row's denominators.
// Its exact intermediate integer coefficients have a separate scratch bound.
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_ore::gcd;
use crate::rv_qp_quad::{self, Element, Field, Wide};
use crate::wide_integer::Natural;

const LIMBS: usize = 64;
type Integer = Wide<LIMBS>;

#[derive(Clone, Copy)]
struct Entry {
    a: Integer,
    b: Integer,
}

impl Entry {
    fn integer(a: i128) -> Self {
        Self {
            a: Integer::from_i128(a),
            b: Integer::from_i128(0),
        }
    }

    fn is_zero(self) -> bool {
        self.a.magnitude.is_zero() && self.b.magnitude.is_zero()
    }

    fn same_value(self, other: Self) -> bool {
        let same = |a: Integer, b: Integer| {
            a.magnitude == b.magnitude && (a.is_zero() || a.negative == b.negative)
        };
        same(self.a, other.a) && same(self.b, other.b)
    }

    fn multiply(self, model: &mut Model, other: Self, d: Integer) -> Result<Self> {
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

    fn divide_exact(self, model: &mut Model, divisor: Self, d: Integer) -> Result<Self> {
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
        let mut quotient = |value: Integer| -> Result<Integer> {
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

fn multiply(model: &mut Model, lhs: Integer, rhs: Integer) -> Result<Integer> {
    model.spend(1 + lhs.magnitude.bits().div_ceil(32) * rhs.magnitude.bits().div_ceil(32))?;
    lhs.multiply_wide(rhs)
}

fn scratch(model: &mut Model, rank: usize) -> Result<()> {
    let words_per_entry = std::mem::size_of::<Entry>().div_ceil(std::mem::size_of::<usize>());
    let cells = rank * (rank * words_per_entry + 3);
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "wide tropical minor exceeds the determinant storage limit",
        ));
    }
    model.spend(cells)
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
                let numerator = match (!entry.is_zero(), has_factor && !upper.is_zero()) {
                    (false, false) => continue,
                    (true, false) => pivot.multiply(model, *entry, d)?,
                    (false, true) => {
                        let cross = factor.multiply(model, *upper, d)?;
                        Entry {
                            a: cross.a.negate(),
                            b: cross.b.negate(),
                        }
                    }
                    (true, true) => {
                        model.spend(LIMBS)?;
                        let product = pivot.multiply(model, *entry, d)?;
                        let cross = factor.multiply(model, *upper, d)?;
                        product.subtract(cross)?
                    }
                };
                // Bareiss's identity makes this quotient exact in the input ring.
                *entry = numerator.divide_exact(model, previous, d)?;
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
    scratch(model, matrix.len())?;
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
    scratch(model, matrix.len())?;
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
