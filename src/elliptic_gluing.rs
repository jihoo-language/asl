use std::rc::Rc;

use crate::ec_poly::{self as poly, Field, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic_chart::homogeneous_substitution;
use crate::hyperelliptic_jacobian::Divisor;
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};

pub(crate) struct Command {
    pub(crate) name: Label,
    pub(crate) factors: [Label; 2],
    pub(crate) matrix: [ParsedExpr; 4],
}

#[derive(Clone, Copy)]
pub(crate) enum Side {
    Left,
    Right,
}

impl Side {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Left => "rvEllGlueLeft",
            Self::Right => "rvEllGlueRight",
        }
    }
}

pub(crate) struct Factor {
    pub(crate) name: Rc<str>,
    pub(crate) curve: Rc<Curve>,
}

pub(crate) struct Gluing {
    pub(crate) name: Rc<str>,
    pub(crate) curve: Rc<Curve>,
    factors: [Factor; 2],
    pub(crate) pole: Scalar,
    pub(crate) image_of_infinity: Scalar,
    pub(crate) numerator: Scalar,
    pub(crate) scale: Scalar,
}

fn require_factors([left, right]: &[Factor; 2]) -> Result<()> {
    for factor in [left, right] {
        let f = factor.curve.polynomial();
        if factor.curve.ring.precision() != 1
            || f.len() != 4
            || f.last() != Some(&factor.curve.ring.scalar(1))
        {
            return Err(Error::plain(format!(
                "rvellglue factor '{}' must be a smooth monic cubic over a named precision-one field",
                factor.name
            )));
        }
    }
    if left.curve.grade.ring != right.curve.grade.ring {
        return Err(Error::plain(
            "rvellglue requires both factors over the same named field",
        ));
    }
    Ok(())
}

impl Gluing {
    pub(crate) fn build(
        model: &mut Model,
        name: Rc<str>,
        factors: [Factor; 2],
        matrix: [Scalar; 4],
    ) -> Result<Self> {
        require_factors(&factors)?;
        let [left, right] = &factors;
        let ring = &left.curve.ring;
        let field = Field { ring };
        let f = left.curve.polynomial();
        let g = right.curve.polynomial();
        let [a, b, c, d] = &matrix;
        let ad = field.multiply(model, a, d)?;
        let bc = field.multiply(model, b, c)?;
        let negative_determinant = field.subtract(model, &bc, &ad)?;
        if field.is_zero(&negative_determinant) {
            return Err(Error::plain("rvellglue matrix determinant must be nonzero"));
        }
        if field.is_zero(c) {
            return Err(Error::plain("rvellglue requires a non-affine branch matching. An affine matching gives a geometrically decomposable quotient"));
        }

        let inverse_c = field.inverse_unit(model, c)?;
        let negative_d = ring.negate(model, d)?;
        let pole = field.multiply(model, &negative_d, &inverse_c)?;
        let image_of_infinity = field.multiply(model, a, &inverse_c)?;
        let inverse_c_squared = field.multiply(model, &inverse_c, &inverse_c)?;
        let numerator = field.multiply(model, &negative_determinant, &inverse_c_squared)?;
        let scale = poly::evaluate(model, field, g, &image_of_infinity)?;
        let shifted = homogeneous_substitution(
            model,
            field,
            f,
            3,
            &[pole.clone(), field.one()],
            &[field.one()],
        )?;
        if field.is_zero(&shifted[0]) || field.is_zero(&scale) {
            return Err(Error::plain(
                "rvellglue branch matching requires f(A) and g(B) to be nonzero",
            ));
        }

        // For u=x−A, the branch map is B+m/u. Certify all branch roots
        // simultaneously by clearing the denominator in g(B+m/u).
        let pulled_back = homogeneous_substitution(
            model,
            field,
            g,
            3,
            &[numerator.clone(), image_of_infinity.clone()],
            &[field.zero(), field.one()],
        )?;
        if pulled_back != poly::scale(model, field, &shifted, &scale)? {
            return Err(Error::plain(
                "rvellglue Möbius map does not match the finite branch roots",
            ));
        }
        let sextic = homogeneous_substitution(
            model,
            field,
            &shifted,
            3,
            &[field.zero(), field.zero(), scale.clone()],
            &[field.one()],
        )?;
        let curve = Rc::new(Curve::new(
            model,
            left.curve.grade.ring.clone(),
            ring.clone(),
            sextic,
        )?);
        model.allocate(10 + 4 * field.cells())?;
        Ok(Self {
            name,
            curve,
            factors,
            pole,
            image_of_infinity,
            numerator,
            scale,
        })
    }

    pub(crate) fn factor(&self, side: Side) -> &Factor {
        match side {
            Side::Left => &self.factors[0],
            Side::Right => &self.factors[1],
        }
    }

    pub(crate) fn point(
        &self,
        model: &mut Model,
        side: Side,
        t: &[u64],
        y: &[u64],
    ) -> Result<Divisor> {
        let ring = &self.curve.ring;
        let field = Field { ring };
        let value = poly::evaluate(model, field, self.curve.polynomial(), t)?;
        if ring.multiply(model, y, y)? != value {
            return Err(Error::plain(format!(
                "point does not lie on the glued genus-two curve '{}'",
                self.name
            )));
        }
        let target = &self.factor(side).curve;
        if matches!(side, Side::Right) && poly::zero(t) {
            return Ok(Divisor::zero(target));
        }
        let t_squared = ring.multiply(model, t, t)?;
        let u = field.multiply(model, &self.scale, &t_squared)?;
        let (x, y) = match side {
            Side::Left => (field.add(model, &self.pole, &u)?, y.to_vec()),
            Side::Right => {
                let denominator = ring.multiply(model, &u, t)?;
                let inverse = field.inverse_unit(model, &denominator)?;
                let reciprocal_u = ring.multiply(model, t, &inverse)?;
                let offset = field.multiply(model, &self.numerator, &reciprocal_u)?;
                (
                    field.add(model, &self.image_of_infinity, &offset)?,
                    ring.multiply(model, y, &inverse)?,
                )
            }
        };
        Divisor::from_point(model, target, &x, y)
    }
}
