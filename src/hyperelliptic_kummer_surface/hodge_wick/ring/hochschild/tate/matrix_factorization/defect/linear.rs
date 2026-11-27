use crate::hw_algebra::CoefficientRing;
use crate::rv_lens::{self, Rational as Scalar};

use super::field::Field;
use super::{Model, Result};

pub(super) type Vector = Vec<Scalar>;

pub(super) struct Rect {
    pub(super) columns: usize,
    pub(super) rows: Vec<Vector>,
}

impl Rect {
    pub(super) fn zero(rows: usize, columns: usize) -> Self {
        Self {
            columns,
            rows: vec![vec![Scalar::zero(); columns]; rows],
        }
    }

    pub(super) fn identity(n: usize) -> Self {
        let mut matrix = Self::zero(n, n);
        for (i, row) in matrix.rows.iter_mut().enumerate() {
            row[i] = Scalar::integer(1);
        }
        matrix
    }

    pub(super) fn compose(&self, model: &mut Model, field: &Field, before: &Self) -> Result<Self> {
        let mut result = Self::zero(self.rows.len(), before.columns);
        for (row, out) in self.rows.iter().zip(&mut result.rows) {
            for (a, b) in row.iter().zip(&before.rows) {
                let minus = rv_lens::negative(model, a)?;
                subtract_multiple(model, field, out, &minus, b, 0)?;
            }
        }
        Ok(result)
    }

    pub(super) fn apply(&self, model: &mut Model, field: &Field, v: &[Scalar]) -> Result<Vector> {
        let mut result = Vec::with_capacity(self.rows.len());
        for row in &self.rows {
            model.spend(self.columns)?;
            let mut coefficient = Scalar::zero();
            for (a, b) in row.iter().zip(v) {
                let term = field.multiply(model, a, b)?;
                coefficient = field.add(model, &coefficient, &term)?;
            }
            result.push(coefficient);
        }
        Ok(result)
    }

    pub(super) fn echelon(&mut self, model: &mut Model, field: &Field) -> Result<Vec<usize>> {
        let mut columns = Vec::new();
        for col in 0..self.columns {
            let rank = columns.len();
            model.spend(1 + self.rows.len() - rank)?;
            let Some(pivot) = (rank..self.rows.len()).find(|&r| !self.rows[r][col].is_zero())
            else {
                continue;
            };
            self.rows.swap(rank, pivot);
            let inverse = field.inverse_unit(model, &self.rows[rank][col])?;
            for entry in &mut self.rows[rank][col..] {
                *entry = field.multiply(model, entry, &inverse)?;
            }
            let pivot = self.rows[rank].clone();
            for (i, row) in self.rows.iter_mut().enumerate() {
                if i != rank {
                    let factor = row[col].clone();
                    subtract_multiple(model, field, row, &factor, &pivot, col)?;
                }
            }
            columns.push(col);
        }
        Ok(columns)
    }

    pub(super) fn kernel(&self, model: &mut Model, pivots: &[usize]) -> Result<Vec<Vector>> {
        let mut basis = Vec::with_capacity(self.columns - pivots.len());
        for j in 0..self.columns {
            if pivots.binary_search(&j).is_ok() {
                continue;
            }
            let mut vector = vec![Scalar::zero(); self.columns];
            vector[j] = Scalar::integer(1);
            for (i, &col) in pivots.iter().enumerate() {
                vector[col] = rv_lens::negative(model, &self.rows[i][j])?;
            }
            basis.push(vector);
        }
        Ok(basis)
    }
}

fn subtract_multiple(
    model: &mut Model,
    field: &Field,
    destination: &mut [Scalar],
    factor: &Scalar,
    source: &[Scalar],
    start: usize,
) -> Result<()> {
    if factor.is_zero() {
        return Ok(());
    }
    model.spend(destination.len() - start)?;
    for (a, b) in destination[start..].iter_mut().zip(&source[start..]) {
        let product = field.multiply(model, factor, b)?;
        *a = field.subtract(model, a, &product)?;
    }
    Ok(())
}

pub(super) fn extend(
    model: &mut Model,
    field: &Field,
    basis: &mut [Option<Vector>],
    mut v: Vector,
) -> Result<bool> {
    for (col, pivot) in basis.iter_mut().enumerate() {
        if v[col].is_zero() {
            continue;
        }
        let factor = v[col].clone();
        if let Some(row) = pivot {
            subtract_multiple(model, field, &mut v, &factor, row, col)?;
        } else {
            let inverse = field.inverse_unit(model, &factor)?;
            for coefficient in &mut v[col..] {
                *coefficient = field.multiply(model, coefficient, &inverse)?;
            }
            *pivot = Some(v);
            return Ok(true);
        }
    }
    Ok(false)
}
