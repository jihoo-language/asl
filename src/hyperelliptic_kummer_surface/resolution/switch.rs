use crate::ec_poly::{Field, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::hw_closed::Curve;
use crate::model::Model;

use super::{field, pivot, quartic, Point, QuarticTerm, ResolvedPoint, TangentCone};

impl ResolvedPoint {
    pub(crate) fn gauss_switch(&self, model: &mut Model, curve: &Curve) -> Result<Self> {
        let field = field(curve);
        let terms = quartic(model, curve)?;
        let normal = if let Some(direction) = self.direction() {
            let cone = TangentCone::new(model, field, self.base(), &terms)?;
            let hessian = cone.hessian(model, field)?;
            let partial = apply(
                model,
                field,
                &hessian,
                cone.axes.map(|axis| &direction[axis]),
            )?;
            complete_covector(model, field, self.base(), cone.axes, partial)?
        } else {
            gradient(model, field, self.base(), &terms)?
        };
        // The alternating polarity has A²=-I, so it also represents A⁻¹.
        let target_coordinates = polarity(model, field, &normal)?;
        let target = Point::new(model, curve, target_coordinates)?;
        let reciprocal = polarity(model, field, self.base().coordinates())?;
        let incidence = dot(model, field, target.coordinates().iter().zip(&reciprocal))?;
        if !field.is_zero(&incidence) {
            return Err(Error::plain(
                "Kummer switch polarity failed its incidence relation",
            ));
        }
        if !target.is_node(model, curve)? {
            let normal = gradient(model, field, &target, &terms)?;
            if !parallel(model, field, &normal, &reciprocal)? {
                return Err(Error::plain(
                    "Kummer switch failed its reciprocal Gauss relation",
                ));
            }
            return Ok(Self {
                base: target.keep(model)?,
                direction: None,
            });
        }

        // Fix the target direction's pivot coordinate at zero. The restricted
        // Hessian is invertible, with the last equation forced by incidence.
        let cone = TangentCone::new(model, field, &target, &terms)?;
        let hessian = cone.hessian(model, field)?;
        let inverse = hw_algebra::inverse(model, &field, &hessian)?
            .ok_or_else(|| Error::plain("Kummer switch node Hessian is singular"))?;
        let partial = apply(
            model,
            field,
            &inverse,
            cone.axes.map(|axis| &reciprocal[axis]),
        )?;
        let mut direction = std::array::from_fn(|_| field.zero());
        for (axis, value) in cone.axes.into_iter().zip(partial) {
            direction[axis] = value;
        }
        let check = apply(
            model,
            field,
            &hessian,
            cone.axes.map(|axis| &direction[axis]),
        )?;
        if complete_covector(model, field, &target, cone.axes, check)? != reciprocal {
            return Err(Error::plain(
                "Kummer switch failed to recover the exceptional direction",
            ));
        }
        let target = target.keep(model)?;
        Self::exceptional(model, curve, target, direction)
    }
}

impl TangentCone {
    fn hessian(&self, model: &mut Model, field: Field<'_>) -> Result<Matrix<Scalar>> {
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
}

fn polarity(
    model: &mut Model,
    field: Field<'_>,
    [w, s, p, t]: &[Scalar; 4],
) -> Result<[Scalar; 4]> {
    Ok([
        t.clone(),
        field.ring.negate(model, p)?,
        s.clone(),
        field.ring.negate(model, w)?,
    ])
}

fn gradient(
    model: &mut Model,
    field: Field<'_>,
    base: &Point,
    terms: &[QuarticTerm],
) -> Result<[Scalar; 4]> {
    let mut result = std::array::from_fn(|_| field.zero());
    for (slots, coefficient) in terms {
        // Replacing one factor at a time retains its multiplicity in the field.
        for (index, &axis) in slots.iter().enumerate() {
            let mut term = coefficient.clone();
            for (other, &slot) in slots.iter().enumerate() {
                if index != other {
                    term = field.multiply(model, &term, &base.coordinates()[slot])?;
                }
            }
            result[axis] = field.add(model, &result[axis], &term)?;
        }
    }
    Ok(result)
}

fn complete_covector(
    model: &mut Model,
    field: Field<'_>,
    base: &Point,
    axes: [usize; 3],
    partial: [Scalar; 3],
) -> Result<[Scalar; 4]> {
    let mut normal = std::array::from_fn(|_| field.zero());
    for (axis, value) in axes.into_iter().zip(partial) {
        normal[axis] = value;
    }
    // At a node H(P)P=0 by the homogeneous Euler identity. Symmetry makes
    // H(P)h annihilate P, whose pivot coordinate is one.
    let sum = dot(model, field, normal.iter().zip(base.coordinates()))?;
    normal[pivot(base)] = field.ring.negate(model, &sum)?;
    Ok(normal)
}

fn parallel(
    model: &mut Model,
    field: Field<'_>,
    left: &[Scalar; 4],
    right: &[Scalar; 4],
) -> Result<bool> {
    let pivot = right
        .iter()
        .position(|value| !field.is_zero(value))
        .expect("polar covector of a nonzero point");
    if field.is_zero(&left[pivot]) {
        return Ok(false);
    }
    for (a, b) in left.iter().zip(right) {
        if field.multiply(model, a, &right[pivot])? != field.multiply(model, b, &left[pivot])? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn apply(
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

fn dot<'a>(
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
