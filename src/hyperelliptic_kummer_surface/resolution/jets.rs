use crate::ec_poly::{Field, Scalar};
use crate::error::Result;
use crate::hw_algebra::{CoefficientRing, Matrix};
use crate::model::Model;

use super::{pivot, Point, QuarticTerm, TangentCone};

pub(super) fn evaluate(
    model: &mut Model,
    field: Field<'_>,
    coordinates: &[Scalar; 4],
    terms: &[QuarticTerm],
) -> Result<Scalar> {
    let mut value = field.zero();
    for (slots, coefficient) in terms {
        let mut term = coefficient.clone();
        for &slot in slots {
            term = field.multiply(model, &term, &coordinates[slot])?;
        }
        value = field.add(model, &value, &term)?;
    }
    Ok(value)
}

pub(super) fn gradient(
    model: &mut Model,
    field: Field<'_>,
    coordinates: &[Scalar; 4],
    terms: &[QuarticTerm],
) -> Result<[Scalar; 4]> {
    let mut result = std::array::from_fn(|_| field.zero());
    for (slots, coefficient) in terms {
        // Replacing one factor at a time retains its multiplicity in the field.
        for (index, &axis) in slots.iter().enumerate() {
            let mut term = coefficient.clone();
            for (other, &slot) in slots.iter().enumerate() {
                if index != other {
                    term = field.multiply(model, &term, &coordinates[slot])?;
                }
            }
            result[axis] = field.add(model, &result[axis], &term)?;
        }
    }
    Ok(result)
}

impl TangentCone {
    pub(super) fn hessian(&self, model: &mut Model, field: Field<'_>) -> Result<Matrix<Scalar>> {
        let mut entries = Vec::with_capacity(9);
        for i in 0..3 {
            for j in 0..3 {
                let coefficient = &self.coefficients[i.min(j)][i.max(j)];
                entries.push(if i == j {
                    // The Hasse quadratic has coefficient H_ii/2 on h_i².
                    field.add(model, coefficient, coefficient)?
                } else {
                    coefficient.clone()
                });
            }
        }
        Matrix::from_entries(3, entries)
    }

    pub(super) fn normal(
        &self,
        model: &mut Model,
        field: Field<'_>,
        base: &Point,
        hessian: &Matrix<Scalar>,
        direction: &[Scalar; 4],
    ) -> Result<[Scalar; 4]> {
        let partial = apply(
            model,
            field,
            hessian,
            self.axes.map(|axis| &direction[axis]),
        )?;
        let mut normal = std::array::from_fn(|_| field.zero());
        for (axis, value) in self.axes.into_iter().zip(partial) {
            normal[axis] = value;
        }
        // At a node H(P)P=0 by the homogeneous Euler identity. Symmetry makes
        // H(P)h annihilate P, whose pivot coordinate is one.
        let sum = dot(model, field, normal.iter().zip(base.coordinates()))?;
        normal[pivot(base)] = field.ring.negate(model, &sum)?;
        Ok(normal)
    }
}

pub(super) fn apply(
    model: &mut Model,
    field: Field<'_>,
    matrix: &Matrix<Scalar>,
    vector: [&Scalar; 3],
) -> Result<[Scalar; 3]> {
    let mut result = std::array::from_fn(|_| field.zero());
    for (row, entry) in matrix.entries().chunks_exact(3).zip(&mut result) {
        *entry = dot(model, field, row.iter().zip(vector))?;
    }
    Ok(result)
}

pub(super) fn dot<'a>(
    model: &mut Model,
    field: Field<'_>,
    terms: impl Iterator<Item = (&'a Scalar, &'a Scalar)>,
) -> Result<Scalar> {
    let mut sum = field.zero();
    for (a, b) in terms {
        let term = field.multiply(model, a, b)?;
        sum = field.add(model, &sum, &term)?;
    }
    Ok(sum)
}
