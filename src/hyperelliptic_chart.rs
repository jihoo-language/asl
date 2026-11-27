use std::rc::Rc;

use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic_jacobian::Divisor;
use crate::model::Model;

#[derive(Clone, Copy)]
pub(crate) enum Direction {
    ToChart,
    ToSource,
}

impl Direction {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::ToChart => "rvJacRechart",
            Self::ToSource => "rvJacUnchart",
        }
    }
}

pub(crate) struct Chart {
    pub(crate) source: Rc<str>,
    pub(crate) branch: Scalar,
    pub(crate) derivative: Scalar,
}

impl Chart {
    pub(crate) fn build(
        model: &mut Model,
        source: Rc<str>,
        curve: &Curve,
        branch: Scalar,
    ) -> Result<(Self, Polynomial)> {
        if curve.ring.precision() != 1 {
            return Err(Error::plain(
                "rvjacchart requires a named residue field of precision one",
            ));
        }
        let field = Field { ring: &curve.ring };
        if !field.is_zero(&poly::evaluate(model, field, curve.polynomial(), &branch)?) {
            return Err(Error::plain("rvjacchart requires a finite branch point"));
        }
        let derivative = poly::derivative(model, field, curve.polynomial())?;
        let derivative = poly::evaluate(model, field, &derivative, &branch)?;
        if field.is_zero(&derivative) {
            return Err(Error::plain("rvjacchart requires a simple branch point"));
        }
        let numerator = poly::trim(vec![derivative.clone(), branch.clone()]);
        let denominator = [field.zero(), field.one()];
        let h = homogeneous_substitution(
            model,
            field,
            curve.polynomial(),
            2 * (curve.grade.genus + 1),
            &numerator,
            &denominator,
        )?;
        let inverse = field.inverse_unit(model, &derivative)?;
        let inverse_square = field.multiply(model, &inverse, &inverse)?;
        let h = poly::scale(model, field, &h, &inverse_square)?;
        model.allocate(4 + branch.len() + derivative.len())?;
        Ok((
            Self {
                source,
                branch,
                derivative,
            },
            h,
        ))
    }

    pub(crate) fn transport(
        &self,
        model: &mut Model,
        direction: Direction,
        target: &Curve,
        divisor: &Divisor,
    ) -> Result<Divisor> {
        let degree = divisor.u.len() - 1;
        if degree == 0 {
            return Ok(Divisor::zero(target));
        }
        let field = Field { ring: &target.ring };
        let zero = field.zero();
        let branch_factor = vec![field.subtract(model, &zero, &self.branch)?, field.one()];
        let zero_factor = vec![zero.clone(), field.one()];
        let inverse = field.inverse_unit(model, &self.derivative)?;
        let (pole, pole_factor, numerator, denominator, y_scale, old_infinity) = match direction {
            Direction::ToChart => (
                &self.branch,
                &branch_factor,
                poly::trim(vec![self.derivative.clone(), self.branch.clone()]),
                &zero_factor,
                inverse,
                &zero,
            ),
            Direction::ToSource => (
                &zero,
                &zero_factor,
                vec![self.derivative.clone()],
                &branch_factor,
                target
                    .ring
                    .power(model, &inverse, target.grade.genus as u64)?,
                &self.branch,
            ),
        };
        // A pole of the coordinate change is a ramification point and can occur
        // only once in a reduced Mumford divisor. It maps to the new infinity.
        let stripped;
        let u = if field.is_zero(&poly::evaluate(model, field, &divisor.u, pole)?) {
            stripped = poly::exact_divide(model, field, &divisor.u, pole_factor)?;
            &stripped
        } else {
            &divisor.u
        };
        let u = homogeneous_substitution(model, field, u, u.len() - 1, &numerator, denominator)?;
        let leading_inverse = field.inverse_unit(model, u.last().expect("nonzero support"))?;
        let u = poly::scale(model, field, &u, &leading_inverse)?;
        let v = homogeneous_substitution(
            model,
            field,
            &divisor.v,
            target.grade.genus + 1,
            &numerator,
            denominator,
        )?;
        let v = poly::scale(model, field, &v, &y_scale)?;
        let mapped = Divisor::new(model, target, u, v)?;
        // Move the base point as well as the finite support. The former infinity
        // is a rational branch point, giving a class of order dividing two.
        if degree.is_multiple_of(2) {
            return Ok(mapped);
        }
        let base_change = Divisor::from_point(model, target, old_infinity, field.zero())?;
        mapped.add(model, target, &base_change)
    }
}

// Return denominator^degree * polynomial(numerator / denominator).
fn homogeneous_substitution(
    model: &mut Model,
    field: Field<'_>,
    polynomial: &[Scalar],
    degree: usize,
    numerator: &[Scalar],
    denominator: &[Scalar],
) -> Result<Polynomial> {
    let Some((leading, rest)) = polynomial.split_last() else {
        return Ok(Vec::new());
    };
    field.check_polynomial(degree.saturating_add(1))?;
    let mut result = vec![leading.clone()];
    let mut denominator_power = vec![field.one()];
    for coefficient in rest.iter().rev() {
        denominator_power = poly::multiply(model, field, &denominator_power, denominator)?;
        result = poly::multiply(model, field, &result, numerator)?;
        let term = poly::scale(model, field, &denominator_power, coefficient)?;
        result = poly::sum(model, field, &result, &term, false)?;
    }
    for _ in rest.len()..degree {
        result = poly::multiply(model, field, &result, denominator)?;
    }
    Ok(result)
}
