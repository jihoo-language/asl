use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_wz::Certificate;

use super::{
    cells, check_cells, exponent_overflow, multiply, Accumulator, Exact, Monomial, Period, Ring,
};

type Polynomial = Vec<(Monomial, i128)>;

fn arithmetic() -> Error {
    Error::plain("exact Laurent lattice rank arithmetic exceeds i128")
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn divide(value: i128, divisor: u128) -> i128 {
    if let Ok(divisor) = i128::try_from(divisor) {
        value / divisor
    } else {
        -i128::from(value == i128::MIN)
    }
}

fn primitive(vector: &mut [i128]) {
    let divisor = vector.iter().map(|v| v.unsigned_abs()).fold(0, gcd);
    if divisor > 1 {
        for value in vector {
            *value = divide(*value, divisor);
        }
    }
}

fn validate_columns(model: &mut Model, variables: usize, columns: &[Vec<i128>]) -> Result<()> {
    check_cells(
        2usize
            .saturating_mul(variables)
            .saturating_mul(columns.len()),
    )?;
    let mut basis = BTreeMap::<usize, Vec<i128>>::new();
    for column in columns {
        let mut vector = column.to_vec();
        model.spend(1 + variables)?;
        primitive(&mut vector);
        // Work over the rationals via fraction-free elimination of the
        // integer columns. Remove content to reduce intermediate growth.
        for (&pivot, row) in &basis {
            model.spend(1)?;
            let entry = vector[pivot];
            if entry != 0 {
                let leading = row[pivot];
                let common = gcd(leading.unsigned_abs(), entry.unsigned_abs());
                let (scale, multiple) = (divide(leading, common), divide(entry, common));
                model.spend(2 * variables)?;
                for (x, &y) in vector.iter_mut().zip(row) {
                    *x = x
                        .checked_mul(scale)
                        .and_then(|scaled| {
                            y.checked_mul(multiple).and_then(|v| scaled.checked_sub(v))
                        })
                        .ok_or_else(arithmetic)?;
                }
                primitive(&mut vector);
            }
        }
        model.spend(variables)?;
        let pivot = vector.iter().position(|&entry| entry != 0).ok_or_else(|| {
            Error::plain("Laurent lattice embedding columns are linearly dependent over Q")
        })?;
        basis.insert(pivot, vector);
    }
    Ok(())
}

fn coordinate(
    model: &mut Model,
    variables: usize,
    columns: Option<&[Vec<i128>]>,
    index: usize,
    inverted: bool,
) -> Result<Monomial> {
    model.spend(variables)?;
    let exponents = if let Some(columns) = columns {
        columns[index]
            .iter()
            .map(|&e| {
                if inverted {
                    e.checked_neg().ok_or_else(exponent_overflow)
                } else {
                    Ok(e)
                }
            })
            .collect::<Result<Vec<_>>>()?
    } else {
        let mut exponents = vec![0; variables];
        exponents[index] = if inverted { -1 } else { 1 };
        exponents
    };
    Ok(Monomial {
        power: 0,
        exponents: exponents.into(),
    })
}

fn canonical(
    model: &mut Model,
    variables: usize,
    x_count: usize,
    rank: usize,
    z: i128,
    columns: Option<&[Vec<i128>]>,
) -> Result<Polynomial> {
    let one = Monomial::unit(variables);
    let mut base = vec![(one.clone(), 1)];
    let mut shift = one.clone();
    for i in 0..rank {
        // The first factors are (1+x) and (1+y)/y. Retain the inverse
        // x-product for the second term of the final bracket.
        let exponent = coordinate(model, variables, columns, i, i >= x_count)?;
        base = multiply(model, &base, &[(one.clone(), 1), (exponent, 1)], Exact)?;
        if i < x_count {
            let inverse = coordinate(model, variables, columns, i, true)?;
            shift = shift.product(&inverse)?;
        }
    }
    let mut shifted = Vec::with_capacity(base.len());
    for (term, coefficient) in &base {
        model.spend(variables)?;
        shifted.push((term.product(&shift)?, Exact.mul(model, *coefficient, z)?));
    }
    for i in x_count..rank {
        let exponent = coordinate(model, variables, columns, i, false)?;
        shifted = multiply(model, &shifted, &[(one.clone(), 1), (exponent, 1)], Exact)?;
    }
    let mut expected = Accumulator::default();
    for (term, coefficient) in base.into_iter().chain(shifted) {
        expected.add(model, term, coefficient, Exact)?;
    }
    Ok(expected.finish())
}

pub(crate) fn check(
    model: &mut Model,
    period: &Period,
    certificate: &Certificate,
    columns: Option<&[Vec<i128>]>,
) -> Result<()> {
    if certificate.u == 0 || certificate.z == 0 {
        return Err(Error::plain("WZ summand requires u >= 1 and nonzero z"));
    }
    let x_count = certificate.u - 1;
    let rank = x_count
        .checked_add(certificate.v)
        .ok_or_else(|| Error::plain("hypergeometric Laurent rank is out of range"))?;
    let variables = period.variables;
    if let Some(columns) = columns {
        if columns.len() != rank || columns.iter().any(|column| column.len() != variables) {
            return Err(Error::plain(format!(
                "Laurent lattice embedding requires {rank} columns of length {variables}"
            )));
        }
        if rank > variables {
            return Err(Error::plain(
                "Laurent lattice embedding has more columns than rows",
            ));
        }
        validate_columns(model, variables, columns)?;
    } else if rank != variables {
        return Err(Error::plain(format!(
            "canonical hypergeometric Laurent rank {rank} requires a matching declared dimension or an embedding"
        )));
    }
    let expected = canonical(model, variables, x_count, rank, certificate.z, columns)?;
    model.spend(cells(&expected) + cells(&period.laurent))?;
    if expected != period.laurent {
        return Err(Error::plain(
            "declared Laurent polynomial does not match the canonical hypergeometric Laurent polynomial",
        ));
    }
    Ok(())
}
