use std::rc::Rc;

use crate::ec_poly::{Field, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::model::Model;

use super::jets::{dot, evaluate, gradient};
use super::{
    field, quartic, remove_radial_component, Point, QuarticTerm, ResolvedPoint, TangentCone,
};

impl ResolvedPoint {
    pub(crate) fn node_involution(
        &self,
        model: &mut Model,
        curve: &Curve,
        center: Rc<Point>,
    ) -> Result<Self> {
        if !center.is_node(model, curve)? {
            return Err(Error::plain(
                "rvKummerNodeInvolution center must be a Kummer node",
            ));
        }
        let field = field(curve);
        let terms = quartic(model, curve)?;
        if self.base() == &center {
            return self.project_center_ray(model, curve, field, &terms);
        }

        // The line through the center N and P has equation
        // F(uN+vP)=u v²(uA+vB). The second jet is invariant along N.
        let point = self.base();
        let mut ray = point.coordinates().clone();
        remove_radial_component(model, field, &center, &mut ray)?;
        let cone = TangentCone::new(model, field, &center, &terms)?;
        let a = cone.evaluate(model, field, cone.axes.map(|axis| &ray[axis]))?;

        if let Some(direction) = self.direction() {
            // At another node B=0. The differential of B is N·H(P)h;
            // the differential of A contributes only a radial term at P.
            if field.is_zero(&a) {
                return Err(contained_line());
            }
            let cone = TangentCone::new(model, field, point, &terms)?;
            let hessian = cone.hessian(model, field)?;
            let normal = cone.normal(model, field, point, &hessian, direction)?;
            let db = dot(model, field, center.coordinates().iter().zip(&normal))?;
            let image =
                subtract_scaled(model, field, (center.coordinates(), &db), (direction, &a))?;
            return Self::exceptional(model, curve, point.clone(), image);
        }

        let normal = gradient(model, field, point.coordinates(), &terms)?;
        let b = dot(model, field, center.coordinates().iter().zip(&normal))?;
        if field.is_zero(&a) {
            if field.is_zero(&b) {
                return Err(contained_line());
            }
            return Self::exceptional(model, curve, center, ray);
        }
        let image = subtract_scaled(
            model,
            field,
            (center.coordinates(), &b),
            (point.coordinates(), &a),
        )?;
        let image = Point::new(model, curve, image)?.keep(model)?;
        Self::smooth(model, curve, image)
    }

    fn project_center_ray(
        &self,
        model: &mut Model,
        curve: &Curve,
        field: Field<'_>,
        terms: &[QuarticTerm],
    ) -> Result<Self> {
        let center = self.base();
        let direction = self.direction().ok_or_else(|| {
            Error::plain("Kummer node projection requires a direction over the center")
        })?;
        // For a tangent direction h, F(uN+v h)=v³(uB+vF(h)).
        let normal = gradient(model, field, direction, terms)?;
        let b = dot(model, field, center.coordinates().iter().zip(&normal))?;
        let value = evaluate(model, field, direction, terms)?;
        if field.is_zero(&b) {
            if field.is_zero(&value) {
                return Err(contained_line());
            }
            return Ok(Self {
                base: center.clone(),
                direction: Some(direction.clone()),
            });
        }
        let image = subtract_scaled(
            model,
            field,
            (center.coordinates(), &value),
            (direction, &b),
        )?;
        let image = Point::new(model, curve, image)?.keep(model)?;
        Self::smooth(model, curve, image)
    }
}

fn contained_line() -> Error {
    Error::plain("Kummer node projection encountered a line contained in the quartic")
}

fn subtract_scaled(
    model: &mut Model,
    field: Field<'_>,
    (left, a): (&[Scalar; 4], &Scalar),
    (right, b): (&[Scalar; 4], &Scalar),
) -> Result<[Scalar; 4]> {
    let mut result = std::array::from_fn(|_| field.zero());
    for ((output, left), right) in result.iter_mut().zip(left).zip(right) {
        let left = field.multiply(model, left, a)?;
        let right = field.multiply(model, right, b)?;
        *output = field.subtract(model, &left, &right)?;
    }
    Ok(result)
}
