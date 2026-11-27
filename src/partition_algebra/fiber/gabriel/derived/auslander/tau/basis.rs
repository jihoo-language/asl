use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::Model;

use super::Vector;
use crate::partition_algebra::fiber::{add, multiply, negate, ratio};

pub(super) fn dual(model: &mut Model, g: &[&[i128]]) -> Result<Vec<Vector>> {
    let n = g.len();
    let one = Rational::integer(1);
    let mut rows: Vec<_> = g
        .iter()
        .enumerate()
        .map(|(i, row)| {
            row.iter()
                .copied()
                .map(Rational::integer)
                .chain((0..n).map(|j| Rational::integer(i128::from(i == j))))
                .collect::<Vec<_>>()
        })
        .collect();
    let mut determinant = one;
    for column in 0..n {
        let pivot = (column..n)
            .find(|&i| rows[i][column] != Rational::ZERO)
            .ok_or_else(|| Error::plain("partition chamber g-vectors are singular"))?;
        if pivot != column {
            rows.swap(pivot, column);
            determinant = negate(determinant)?;
        }
        let top = rows[column][column];
        determinant = multiply(determinant, top)?;
        let inverse = ratio(top.denominator, top.numerator)?;
        for value in &mut rows[column] {
            *value = multiply(*value, inverse)?;
        }
        let normalized = rows[column].clone();
        for (i, row) in rows.iter_mut().enumerate() {
            if i == column {
                continue;
            }
            model.spend(2 * n)?;
            let scale = negate(row[column])?;
            for (x, &y) in row.iter_mut().zip(&normalized) {
                *x = add(*x, multiply(y, scale)?)?;
            }
        }
    }
    if determinant.denominator != 1 || determinant.numerator.unsigned_abs() != 1 {
        return Err(Error::plain("partition chamber g-basis is not unimodular"));
    }
    let mut columns = vec![Vec::with_capacity(n); n];
    for row in &rows {
        for (value, column) in row[n..].iter().zip(&mut columns) {
            if value.denominator != 1 {
                return Err(Error::plain("partition c-vector is not integral"));
            }
            column.push(value.numerator);
        }
    }
    for (i, row) in g.iter().enumerate() {
        for (j, column) in columns.iter().enumerate() {
            let dot=row.iter().zip(column).try_fold(0i128,|sum,(&a,&b)| a.checked_mul(b).and_then(|x|sum.checked_add(x)).ok_or_else(|| Error::plain("partition dual g-basis verification arithmetic exceeds signed 128 bits")))?;
            if dot != i128::from(i == j) {
                return Err(Error::plain(
                    "partition c-vectors failed the exact dual-basis identity",
                ));
            }
        }
    }
    Ok(columns)
}
