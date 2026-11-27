use std::rc::Rc;

use crate::ec_poly::{Field, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing};
use crate::hw_closed::Curve;
use crate::model::Model;

use super::jets::{apply, dot, gradient};
use super::{field, quartic, Point, QuarticTerm, ResolvedPoint, TangentCone};

impl ResolvedPoint {
    pub(crate) fn gauss_switch(&self, model: &mut Model, curve: &Curve) -> Result<Self> {
        let terms = quartic(model, curve)?;
        let image = self.gauss_switch_transient(model, curve, &terms)?;
        image.base.account_storage(model)?;
        Ok(image)
    }

    pub(in crate::hyperelliptic_kummer_surface) fn gauss_switch_transient(
        &self,
        model: &mut Model,
        curve: &Curve,
        terms: &[QuarticTerm],
    ) -> Result<Self> {
        let field = field(curve);
        let normal = if let Some(direction) = self.direction() {
            let cone = TangentCone::new(model, field, self.base(), terms)?;
            let hessian = cone.hessian(model, field)?;
            cone.normal(model, field, self.base(), &hessian, direction)?
        } else {
            gradient(model, field, self.base().coordinates(), terms)?
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
            let normal = gradient(model, field, target.coordinates(), terms)?;
            if !parallel(model, field, &normal, &reciprocal)? {
                return Err(Error::plain(
                    "Kummer switch failed its reciprocal Gauss relation",
                ));
            }
            return Ok(Self {
                base: Rc::new(target),
                direction: None,
            });
        }

        // Fix the target direction's pivot coordinate at zero. The restricted
        // Hessian is invertible, with the last equation forced by incidence.
        let cone = TangentCone::new(model, field, &target, terms)?;
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
        if cone.normal(model, field, &target, &hessian, &direction)? != reciprocal {
            return Err(Error::plain(
                "Kummer switch failed to recover the exceptional direction",
            ));
        }
        Self::exceptional_with_quartic(model, curve, Rc::new(target), direction, Some(terms))
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
