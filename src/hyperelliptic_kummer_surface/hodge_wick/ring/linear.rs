use crate::error::Result;
use crate::model::Model;

use super::super::rational::{Rational, LIMBS};
use super::super::storage;

pub(super) type Matrix = Vec<Vec<Rational>>;

pub(super) fn matrix(rows: usize, columns: usize) -> Result<Matrix> {
    storage(rows.saturating_mul(columns.saturating_mul(2 * LIMBS + 1).saturating_add(4)))?;
    if rows == 0 {
        return Ok(Vec::new());
    }
    Ok(vec![vec![Rational::zero(); columns]; rows])
}

pub(super) struct Echelon {
    pub(super) rows: Matrix,
    pub(super) pivots: Vec<usize>,
}

pub(super) fn rref(model: &mut Model, mut rows: Matrix, columns: usize) -> Result<Echelon> {
    let mut pivots = Vec::new();
    for column in 0..columns {
        let pivot = pivots.len();
        model.spend(rows.len() - pivot)?;
        let Some(row) = (pivot..rows.len()).find(|&row| !rows[row][column].is_zero()) else {
            continue;
        };
        rows.swap(row, pivot);
        let inverse = rows[pivot][column].inverse()?;
        for entry in &mut rows[pivot][column..] {
            *entry = entry.multiply(model, inverse)?;
        }
        for i in 0..rows.len() {
            if i == pivot || rows[i][column].is_zero() {
                continue;
            }
            let correction = rows[i][column].negate();
            let (target, source) = if i < pivot {
                let (earlier, later) = rows.split_at_mut(pivot);
                (&mut earlier[i], &later[0])
            } else {
                let (earlier, later) = rows.split_at_mut(i);
                (&mut later[0], &earlier[pivot])
            };
            add_multiple(model, &mut target[column..], &source[column..], correction)?;
        }
        pivots.push(column);
        if pivots.len() == rows.len() {
            break;
        }
    }
    rows.truncate(pivots.len());
    Ok(Echelon { rows, pivots })
}

impl Echelon {
    pub(super) fn kernel(&self, columns: usize) -> Result<Matrix> {
        let mut result = matrix(columns - self.pivots.len(), columns)?;
        let free = (0..columns).filter(|i| self.pivots.binary_search(i).is_err());
        for (row, free) in result.iter_mut().zip(free) {
            row[free] = Rational::integer(1);
            for (&pivot, basis) in self.pivots.iter().zip(&self.rows) {
                row[pivot] = basis[free].negate();
            }
        }
        Ok(result)
    }
}

pub(super) fn restrict(model: &mut Model, gram: &Matrix, basis: &Matrix) -> Result<Matrix> {
    let mut product = matrix(basis.len(), gram.len())?;
    for (out, row) in product.iter_mut().zip(basis) {
        for (&factor, source) in row.iter().zip(gram) {
            add_multiple(model, out, source, factor)?;
        }
    }
    let mut result = matrix(basis.len(), basis.len())?;
    for (i, row) in product.iter().enumerate() {
        for j in i..basis.len() {
            model.spend(row.len())?;
            let mut sum = Rational::zero();
            for (&a, &b) in row.iter().zip(&basis[j]) {
                let term = a.multiply(model, b)?;
                sum = sum.add(model, term)?;
            }
            result[i][j] = sum;
            result[j][i] = sum;
        }
    }
    Ok(result)
}

pub(super) fn add_multiple(
    model: &mut Model,
    row: &mut [Rational],
    source: &[Rational],
    factor: Rational,
) -> Result<()> {
    if !factor.is_zero() {
        model.spend(row.len())?;
        for (output, &value) in row.iter_mut().zip(source) {
            let addend = value.multiply(model, factor)?;
            *output = output.add(model, addend)?;
        }
    }
    Ok(())
}

pub(super) fn inertia(model: &mut Model, mut rows: Matrix) -> Result<[usize; 3]> {
    let mut signs = [0; 3];
    let mut pivot = 0;
    while pivot < rows.len() {
        let n = rows.len();
        model.spend((n - pivot).saturating_mul(n - pivot))?;
        if let Some(i) = (pivot..n).find(|&i| !rows[i][i].is_zero()) {
            swap_axis(&mut rows, i, pivot);
            let entry = rows[pivot][pivot];
            let inverse = entry.inverse()?.negate();
            signs[usize::from(entry.negative())] += 1;
            for i in pivot + 1..n {
                let factor = rows[i][pivot].multiply(model, inverse)?;
                let (earlier, later) = rows.split_at_mut(i + 1);
                let row = &mut earlier[i];
                let diagonal = factor.multiply(model, row[pivot])?;
                row[i] = row[i].add(model, diagonal)?;
                for (entry, other) in row[i + 1..].iter_mut().zip(later) {
                    let correction = factor.multiply(model, other[pivot])?;
                    *entry = entry.add(model, correction)?;
                    other[i] = *entry;
                }
            }
            pivot += 1;
        } else if let Some((i, j)) =
            (pivot..n).find_map(|i| (i + 1..n).find(|&j| !rows[i][j].is_zero()).map(|j| (i, j)))
        {
            // The pivot block [[0,b],[b,0]] has one sign of each kind.
            swap_axis(&mut rows, i, pivot);
            swap_axis(&mut rows, j, pivot + 1);
            let inverse = rows[pivot][pivot + 1].inverse()?.negate();
            signs[0] += 1;
            signs[1] += 1;
            for i in pivot + 2..n {
                let a = rows[i][pivot].multiply(model, inverse)?;
                let b = rows[i][pivot + 1].multiply(model, inverse)?;
                let (earlier, later) = rows.split_at_mut(i + 1);
                let row = &mut earlier[i];
                let first = a.multiply(model, row[pivot + 1])?;
                let second = b.multiply(model, row[pivot])?;
                row[i] = row[i].add(model, first)?.add(model, second)?;
                for (entry, other) in row[i + 1..].iter_mut().zip(later) {
                    let first = a.multiply(model, other[pivot + 1])?;
                    let second = b.multiply(model, other[pivot])?;
                    *entry = entry.add(model, first)?.add(model, second)?;
                    other[i] = *entry;
                }
            }
            pivot += 2;
        } else {
            signs[2] = n - pivot;
            break;
        }
    }
    Ok(signs)
}

fn swap_axis(rows: &mut Matrix, i: usize, j: usize) {
    if i != j {
        rows.swap(i, j);
        for row in rows {
            row.swap(i, j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hyperbolic_pivot_preserves_the_residual_sign_and_nullity() {
        let mut model = Model::new();
        let gram = [[0, 2, 1, 0], [2, 0, 3, 0], [1, 3, 0, 0], [0, 0, 0, 0]]
            .into_iter()
            .map(|row| row.into_iter().map(Rational::integer).collect())
            .collect();
        let result = inertia(&mut model, gram).unwrap_or_else(|err| panic!("{err}"));
        assert_eq!(result, [1, 2, 1]);
    }
}
