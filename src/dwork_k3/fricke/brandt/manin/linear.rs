use crate::error::Result;
use crate::lee_linear::Rational as Q;
use crate::model::Model;

use super::{arithmetic, certificate};

pub(super) type Vector = Vec<Q>;
pub(super) type Matrix = Vec<Vector>;

pub(super) const ZERO: Q = Q::ZERO;
pub(super) const ONE: Q = Q::integer(1);

pub(super) fn neg(a: Q) -> Result<Q> {
    Q::new(
        a.numerator.checked_neg().ok_or_else(arithmetic)?,
        a.denominator,
    )
}

pub(super) fn sub(a: Q, b: Q) -> Result<Q> {
    a.add(neg(b)?)
}

pub(super) fn identity(n: usize) -> Matrix {
    (0..n)
        .map(|i| (0..n).map(|j| Q::integer(i128::from(i == j))).collect())
        .collect()
}

pub(super) fn columns(columns: &[Vector], rows: usize) -> Matrix {
    (0..rows)
        .map(|i| columns.iter().map(|v| v[i]).collect())
        .collect()
}

pub(super) fn apply(model: &mut Model, matrix: &Matrix, value: &[Q]) -> Result<Vector> {
    model.spend(matrix.len() * value.len())?;
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .zip(value)
                .try_fold(ZERO, |sum, (&a, &b)| sum.add(a.multiply(b)?))
        })
        .collect()
}

pub(super) fn multiply(model: &mut Model, left: &Matrix, right: &Matrix) -> Result<Matrix> {
    let m = left.len();
    let n = right.first().map_or(0, Vec::len);
    model.spend(m * n * right.len())?;
    let mut out = vec![vec![ZERO; n]; m];
    for (row, l) in out.iter_mut().zip(left) {
        for (r, &a) in right.iter().zip(l) {
            if a.numerator == 0 {
                continue;
            }
            for (c, &b) in row.iter_mut().zip(r) {
                *c = c.add(a.multiply(b)?)?;
            }
        }
    }
    Ok(out)
}

pub(super) fn add_scaled(model: &mut Model, out: &mut [Q], input: &[Q], scale: Q) -> Result<()> {
    model.spend(out.len())?;
    for (a, &b) in out.iter_mut().zip(input) {
        *a = a.add(b.multiply(scale)?)?;
    }
    Ok(())
}

// Reduced row echelon form with pivots restricted to the prescribed prefix.
pub(super) fn rref(
    model: &mut Model,
    matrix: &mut Matrix,
    width: usize,
) -> Result<(Vec<usize>, Q)> {
    model.spend(matrix.len() * width)?;
    let mut pivots = Vec::new();
    let mut determinant = ONE;
    for j in 0..width {
        let i = pivots.len();
        let Some(pivot) = (i..matrix.len()).find(|&k| matrix[k][j].numerator != 0) else {
            continue;
        };
        if i != pivot {
            matrix.swap(i, pivot);
            determinant = neg(determinant)?;
        }
        let a = matrix[i][j];
        if matrix.len() == width {
            determinant = determinant.multiply(a)?;
        }
        let inverse = Q::ratio(a.denominator, a.numerator)?;
        model.spend(matrix[i].len() - j)?;
        for value in &mut matrix[i][j..] {
            *value = value.multiply(inverse)?;
        }
        let row = matrix[i].clone();
        for (k, target) in matrix.iter_mut().enumerate() {
            if k != i && target[j].numerator != 0 {
                let coefficient = neg(target[j])?;
                add_scaled(model, &mut target[j..], &row[j..], coefficient)?;
            }
        }
        pivots.push(j);
    }
    Ok((pivots, determinant))
}

pub(super) fn inverse(model: &mut Model, matrix: &Matrix) -> Result<Option<(Matrix, Q)>> {
    let n = matrix.len();
    let mut augmented: Matrix = matrix
        .iter()
        .zip(identity(n))
        .map(|(row, unit)| row.iter().copied().chain(unit).collect())
        .collect();
    let (pivots, determinant) = rref(model, &mut augmented, n)?;
    if pivots.len() == n {
        Ok(Some((
            augmented.into_iter().map(|row| row[n..].to_vec()).collect(),
            determinant,
        )))
    } else {
        Ok(None)
    }
}

pub(super) fn free_columns(width: usize, pivots: &[usize]) -> Vec<usize> {
    (0..width)
        .filter(|j| pivots.binary_search(j).is_err())
        .collect()
}

pub(super) struct Frame {
    pub(super) vectors: Vec<Vector>,
    pub(super) free: Vec<usize>,
    width: usize,
}

impl Frame {
    pub(super) fn nullspace(
        model: &mut Model,
        width: usize,
        mut equations: Matrix,
    ) -> Result<Self> {
        let (pivots, _) = rref(model, &mut equations, width)?;
        let free = free_columns(width, &pivots);
        let mut vectors = Vec::with_capacity(free.len());
        for &j in &free {
            let mut vector = vec![ZERO; width];
            vector[j] = ONE;
            for (row, &pivot) in pivots.iter().enumerate() {
                vector[pivot] = neg(equations[row][j])?;
            }
            vectors.push(vector);
        }
        Ok(Self {
            vectors,
            free,
            width,
        })
    }

    pub(super) fn coordinates(&self, model: &mut Model, value: &[Q]) -> Result<Option<Vector>> {
        let coordinates = self.free.iter().map(|&j| value[j]).collect::<Vector>();
        let recovered = self.expand(model, &coordinates)?;
        Ok((value == recovered).then_some(coordinates))
    }

    pub(super) fn expand(&self, model: &mut Model, value: &[Q]) -> Result<Vector> {
        let mut out = vec![ZERO; self.width];
        for (v, &c) in self.vectors.iter().zip(value) {
            if c.numerator != 0 {
                add_scaled(model, &mut out, v, c)?;
            }
        }
        Ok(out)
    }

    pub(super) fn restrict(&self, model: &mut Model, operator: &Matrix) -> Result<Matrix> {
        let mut images = Vec::with_capacity(self.vectors.len());
        for basis in &self.vectors {
            let image = apply(model, operator, basis)?;
            images.push(
                self.coordinates(model, &image)?
                    .ok_or_else(|| certificate("cuspidal sector is invariant"))?,
            );
        }
        Ok(columns(&images, images.len()))
    }
}

pub(super) fn characteristic(model: &mut Model, matrix: &Matrix) -> Result<Vector> {
    let n = matrix.len();
    let mut power = identity(n);
    let mut traces = vec![ZERO];
    for _ in 0..n {
        power = multiply(model, &power, matrix)?;
        traces.push(
            power
                .iter()
                .enumerate()
                .try_fold(ZERO, |trace, (j, row)| trace.add(row[j]))?,
        );
    }
    let mut coefficients = vec![ONE];
    for i in 1..=n {
        let sum = (1..=i).try_fold(ZERO, |sum, j| {
            sum.add(coefficients[i - j].multiply(traces[j])?)
        })?;
        coefficients.push(sum.multiply(Q::new(-1, i as i128)?)?);
    }
    Ok(coefficients)
}

pub(super) fn describe(value: &[Q]) -> String {
    let parts = value
        .iter()
        .map(|a| {
            if a.denominator == 1 {
                a.numerator.to_string()
            } else {
                format!("{}/{}", a.numerator, a.denominator)
            }
        })
        .collect::<Vec<_>>();
    format!("[{}]", parts.join(", "))
}

pub(super) fn describe_matrix(matrix: &Matrix) -> String {
    format!(
        "[{}]",
        matrix
            .iter()
            .map(|row| describe(row))
            .collect::<Vec<_>>()
            .join(", ")
    )
}
