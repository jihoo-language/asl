use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::Model;

use super::{add, multiply, ratio, Coordinates};

const PRIMES: [u64; 3] = [1_000_000_007, 1_000_000_009, 998_244_353];

pub(super) struct Kernel {
    // The index is the free column, whose coefficient is one. All other
    // free columns are zero, so the exact rows also order the quotient.
    pub(super) rows: BTreeMap<usize, Coordinates>,
}

pub(super) fn certify(model: &mut Model, gram: &[Vec<i128>]) -> Result<Kernel> {
    for prime in PRIMES {
        let (rref, pivots) = echelon(model, gram, prime)?;
        if let Some(rows) = lift(model, gram, &rref, &pivots, prime)? {
            // The modular pivots certify a nonzero minor of the integer
            // matrix. Each exactly verified rational kernel vector has a
            // distinct free coordinate, giving the matching upper bound.
            return Ok(Kernel { rows });
        }
    }
    Err(Error::plain("could not certify the rational kernel of the partition regular trace form within the reconstruction bound"))
}

fn echelon(
    model: &mut Model,
    gram: &[Vec<i128>],
    prime: u64,
) -> Result<(Vec<Vec<u64>>, Vec<usize>)> {
    let size = gram.len();
    model.spend(size * size)?;
    let mut rows: Vec<Vec<u64>> = gram
        .iter()
        .map(|row| {
            row.iter()
                .map(|&value| value.rem_euclid(i128::from(prime)) as u64)
                .collect()
        })
        .collect();
    let mut pivots = Vec::with_capacity(size);
    for column in 0..size {
        let position = pivots.len();
        model.spend(size - position)?;
        let Some(index) = (position..size).find(|&i| rows[i][column] != 0) else {
            continue;
        };
        rows.swap(position, index);
        let inv = power(rows[position][column], prime - 2, prime);
        model.spend(size - column + 64)?;
        for value in &mut rows[position][column..] {
            *value = *value * inv % prime;
        }
        let pivot = rows[position][column..].to_vec();
        for (i, row) in rows.iter_mut().enumerate() {
            let leading = row[column];
            if i == position || leading == 0 {
                continue;
            }
            model.spend(size - column)?;
            for (value, &factor) in row[column..].iter_mut().zip(&pivot) {
                *value = (*value + prime - leading * factor % prime) % prime;
            }
        }
        pivots.push(column);
    }
    Ok((rows, pivots))
}

fn power(mut base: u64, mut exponent: u64, prime: u64) -> u64 {
    let mut value = 1;
    while exponent != 0 {
        if exponent & 1 == 1 {
            value = value * base % prime;
        }
        base = base * base % prime;
        exponent >>= 1;
    }
    value
}

fn lift(
    model: &mut Model,
    gram: &[Vec<i128>],
    rref: &[Vec<u64>],
    pivots: &[usize],
    prime: u64,
) -> Result<Option<BTreeMap<usize, Coordinates>>> {
    let mut kernel = BTreeMap::new();
    for free in (0..gram.len()).filter(|i| pivots.binary_search(i).is_err()) {
        let mut coordinates = BTreeMap::from([(free, Rational::integer(1))]);
        for (row, &column) in pivots.iter().enumerate() {
            let value = rref[row][free];
            if value != 0 {
                model.spend(32)?;
                let Some((n, d)) = reconstruct(prime - value, prime) else {
                    return Ok(None);
                };
                coordinates.insert(column, ratio(i128::from(n), i128::from(d))?);
            }
        }
        for row in gram {
            model.spend(coordinates.len())?;
            let mut value = Rational::ZERO;
            for (&column, &coefficient) in &coordinates {
                if row[column] != 0 {
                    value = add(
                        value,
                        multiply(Rational::integer(row[column]), coefficient)?,
                    )?;
                }
            }
            if value.numerator != 0 {
                return Ok(None);
            }
        }
        kernel.insert(free, coordinates);
    }
    Ok(Some(kernel))
}

fn reconstruct(residue: u64, prime: u64) -> Option<(i64, i64)> {
    let [prime, residue] = [prime as i64, residue as i64];
    let mut remainders = [prime, residue];
    let mut coefficients = [0, 1];
    let bound = (prime / 2).isqrt();
    while remainders[1] > bound {
        let quotient = remainders[0] / remainders[1];
        remainders = [remainders[1], remainders[0] - quotient * remainders[1]];
        coefficients = [
            coefficients[1],
            coefficients[0] - quotient * coefficients[1],
        ];
    }
    let n = remainders[1];
    let d = coefficients[1];
    (d != 0 && d.abs() <= bound && (n - d * residue).rem_euclid(prime) == 0)
        .then(|| (d.signum() * n, d.abs()))
}
