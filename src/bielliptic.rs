use std::rc::Rc;

use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic_chart::homogeneous_substitution;
use crate::hyperelliptic_jacobian::Divisor;
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};

pub(crate) struct Command {
    pub(crate) source: Label,
    pub(crate) targets: [Label; 2],
    pub(crate) matrix: [ParsedExpr; 4],
}

struct Chart {
    matrix: [Scalar; 4],
    determinant_cubed: Scalar,
}

enum Factor {
    Plus,
    Minus,
}

pub(crate) struct Quotient {
    pub(crate) source: Rc<Curve>,
    pub(crate) target: Rc<Curve>,
    chart: Rc<Chart>,
    factor: Factor,
    scale: Scalar,
}

impl Quotient {
    pub(crate) fn build(
        model: &mut Model,
        source: Rc<Curve>,
        matrix: [Scalar; 4],
    ) -> Result<[Self; 2]> {
        let field = Field { ring: &source.ring };
        let [a, b, c, d] = &matrix;
        let ad = field.multiply(model, a, d)?;
        let bc = field.multiply(model, b, c)?;
        let determinant = field.subtract(model, &ad, &bc)?;
        if field.is_zero(&determinant) {
            return Err(Error::plain("bielliptic chart determinant must be nonzero"));
        }
        let numerator = poly::trim(vec![b.clone(), a.clone()]);
        let denominator = poly::trim(vec![d.clone(), c.clone()]);
        let transformed = homogeneous_substitution(
            model,
            field,
            source.polynomial(),
            6,
            &numerator,
            &denominator,
        )?;
        let invalid_degree =
            || Error::plain("bielliptic chart requires nonzero constant and sextic coefficients");
        let [h0, odd1, h1, odd3, h2, odd5, h3]: [Scalar; 7] =
            transformed.try_into().map_err(|_| invalid_degree())?;
        if field.is_zero(&h0) {
            return Err(invalid_degree());
        }
        if [odd1, odd3, odd5].iter().any(|c| !field.is_zero(c)) {
            return Err(Error::plain("bielliptic chart must give an even sextic"));
        }

        let plus = monic_cubic(model, field, [&h0, &h1, &h2, &h3])?;
        let minus = monic_cubic(model, field, [&h3, &h2, &h1, &h0])?;
        let determinant_cubed = field.ring.power(model, &determinant, 3)?;
        model.allocate(12 + 7 * field.cells())?;
        let chart = Rc::new(Chart {
            matrix,
            determinant_cubed,
        });
        let make = |model: &mut Model, coefficients, factor, scale| {
            let target = Curve::new(
                model,
                source.grade.ring.clone(),
                source.ring.clone(),
                coefficients,
            )?;
            Ok(Self {
                source: source.clone(),
                target: Rc::new(target),
                chart: chart.clone(),
                factor,
                scale,
            })
        };
        Ok([
            make(model, plus, Factor::Plus, h3)?,
            make(model, minus, Factor::Minus, h0)?,
        ])
    }

    pub(crate) fn point(&self, model: &mut Model, x: &[u64], y: &[u64]) -> Result<Divisor> {
        let ring = &self.source.ring;
        let field = Field { ring };
        let f_at_x = poly::evaluate(model, field, self.source.polynomial(), x)?;
        if ring.multiply(model, y, y)? != f_at_x {
            return Err(Error::plain(
                "point does not lie on the bielliptic source curve",
            ));
        }

        // Inverting the Möbius substitution gives t=n/m and Y=det³*y/m³.
        // For the minus factor use 1/t and Y/t³ in place of t and Y.
        let [a, b, c, d] = &self.chart.matrix;
        let dx = ring.multiply(model, d, x)?;
        let cx = ring.multiply(model, c, x)?;
        let n = field.subtract(model, &dx, b)?;
        let m = field.subtract(model, a, &cx)?;
        let (numerator, denominator) = match self.factor {
            Factor::Plus => (n, m),
            Factor::Minus => (m, n),
        };
        if field.is_zero(&denominator) {
            return Ok(Divisor::zero(&self.target));
        }
        let inverse = field.inverse_unit(model, &denominator)?;
        let ratio = field.multiply(model, &numerator, &inverse)?;
        let square = field.multiply(model, &ratio, &ratio)?;
        let image_x = field.multiply(model, &self.scale, &square)?;
        let inverse_cube = ring.power(model, &inverse, 3)?;
        let ordinate = ring.multiply(model, &self.chart.determinant_cubed, y)?;
        let ordinate = field.multiply(model, &ordinate, &inverse_cube)?;
        let image_y = field.multiply(model, &self.scale, &ordinate)?;
        Divisor::from_point(model, &self.target, &image_x, image_y)
    }
}

fn monic_cubic(
    model: &mut Model,
    field: Field<'_>,
    [constant, linear, quadratic, leading]: [&Scalar; 4],
) -> Result<Polynomial> {
    // For the cubic G, use X=leading*u and W=leading*v on v²=G(u).
    let square = field.multiply(model, leading, leading)?;
    Ok(vec![
        field.multiply(model, constant, &square)?,
        field.multiply(model, linear, leading)?,
        quadratic.clone(),
        field.one(),
    ])
}
