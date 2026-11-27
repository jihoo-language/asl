use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;

use super::{Field, Matrix};

pub(super) const AXES: usize = 4;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(super) struct Monomial(pub(super) [u8; AXES]);

impl Monomial {
    pub(super) const ONE: Self = Self([0; AXES]);

    pub(super) fn linear(axis: usize) -> Self {
        let mut array = [0; AXES];
        array[axis] = 1;
        Self(array)
    }

    pub(super) fn degree(self) -> u8 {
        self.0.iter().sum()
    }

    pub(super) fn product(self, b: Self, cutoff: u8) -> Option<Self> {
        let p = Self(std::array::from_fn(|i| self.0[i] + b.0[i]));
        (p.degree() <= cutoff).then_some(p)
    }
}

pub(super) fn homogeneous(axes: usize, degree: u8) -> Vec<Monomial> {
    fn visit(
        out: &mut Vec<Monomial>,
        powers: &mut [u8; AXES],
        remaining: u8,
        position: usize,
        count: usize,
    ) {
        if position + 1 == count {
            powers[position] = remaining;
            out.push(Monomial(*powers));
        } else {
            for power in 0..=remaining {
                powers[position] = power;
                visit(out, powers, remaining - power, position + 1, count);
            }
        }
    }
    let mut out = Vec::new();
    visit(&mut out, &mut [0; AXES], degree, 0, axes);
    out
}

#[derive(Clone, Eq, PartialEq)]
pub(super) struct Series {
    pub(super) rows: usize,
    pub(super) cols: usize,
    pub(super) terms: BTreeMap<Monomial, Matrix>,
}

impl Series {
    pub(super) fn zero(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            terms: BTreeMap::new(),
        }
    }

    pub(super) fn constant(matrix: Matrix) -> Self {
        let mut out = Self::zero(matrix.rows, matrix.cols);
        if !matrix.is_zero() {
            out.terms.insert(Monomial::ONE, matrix);
        }
        out
    }

    pub(super) fn identity(n: usize) -> Self {
        Self::constant(Matrix::identity(n))
    }

    pub(super) fn add_term(
        &mut self,
        model: &mut Model,
        field: &Field,
        monomial: Monomial,
        matrix: &Matrix,
        scalar: u32,
    ) -> Result<()> {
        if self.rows != matrix.rows || self.cols != matrix.cols {
            return Err(Error::plain(
                "partition formal matrix sum has incompatible dimensions",
            ));
        }
        model.spend(1 + self.rows * self.cols)?;
        let value = self
            .terms
            .entry(monomial)
            .or_insert_with(|| Matrix::zero(self.rows, self.cols));
        value.add_scaled(field, matrix, scalar);
        if value.is_zero() {
            self.terms.remove(&monomial);
        }
        Ok(())
    }

    pub(super) fn add(
        &mut self,
        model: &mut Model,
        field: &Field,
        rhs: &Self,
        scale: u32,
    ) -> Result<()> {
        for (&monomial, matrix) in &rhs.terms {
            self.add_term(model, field, monomial, matrix, scale)?;
        }
        Ok(())
    }

    pub(super) fn compose(
        &self,
        model: &mut Model,
        field: &Field,
        rhs: &Self,
        cutoff: u8,
    ) -> Result<Self> {
        if self.cols != rhs.rows {
            return Err(Error::plain(
                "partition formal matrix product has incompatible dimensions",
            ));
        }
        let mut out = Self::zero(self.rows, rhs.cols);
        model.spend(self.terms.len() * rhs.terms.len() + 1)?;
        for (&left, a) in &self.terms {
            for (&right, b) in &rhs.terms {
                if let Some(power) = left.product(right, cutoff) {
                    let matrix = a.compose(model, field, b)?;
                    out.add_term(model, field, power, &matrix, 1)?;
                }
            }
        }
        Ok(out)
    }
}
