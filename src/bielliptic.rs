use std::rc::Rc;

use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic_chart::{homogeneous_substitution, Chart};
use crate::hyperelliptic_jacobian::Divisor;
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod picard;

pub(crate) struct Command {
    pub(crate) source: Label,
    pub(crate) targets: [Label; 2],
    pub(crate) matrix: [ParsedExpr; 4],
}

pub(crate) struct Coordinates {
    numerator: Polynomial,
    denominator: Polynomial,
    offset: Scalar,
    scale: Scalar,
    ordinate_scale: Scalar,
}

pub(crate) struct Quotient {
    pub(crate) source: Rc<Curve>,
    pub(crate) target: Rc<Curve>,
    coordinates: Rc<Coordinates>,
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
        // The inverse substitution gives t=N/M and Y=det³*y/M³.
        let n = poly::trim(vec![field.ring.negate(model, b)?, d.clone()]);
        let m = poly::trim(vec![a.clone(), field.ring.negate(model, c)?]);
        let determinant_cubed = field.ring.power(model, &determinant, 3)?;
        model.allocate(16 + (2 * (n.len() + m.len()) + 6) * field.cells())?;
        let make = |model: &mut Model, coefficients, numerator, denominator, scale| {
            let target = Curve::new(
                model,
                source.grade.ring.clone(),
                source.ring.clone(),
                coefficients,
            )?;
            let ordinate_scale = field.multiply(model, &scale, &determinant_cubed)?;
            Ok(Self {
                source: source.clone(),
                target: Rc::new(target),
                coordinates: Rc::new(Coordinates::new(
                    numerator,
                    denominator,
                    field.zero(),
                    scale,
                    ordinate_scale,
                )),
            })
        };
        // The second quotient uses 1/t and Y/t³.
        Ok([
            make(model, plus, n.clone(), m.clone(), h3)?,
            make(model, minus, m, n, h0)?,
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

        self.coordinates.point(model, &self.target, x, y)
    }
}

impl Coordinates {
    pub(crate) fn new(
        numerator: Polynomial,
        denominator: Polynomial,
        offset: Scalar,
        scale: Scalar,
        ordinate_scale: Scalar,
    ) -> Self {
        Self {
            numerator,
            denominator,
            offset,
            scale,
            ordinate_scale,
        }
    }

    fn rechart(&self, model: &mut Model, field: Field<'_>, chart: &Chart) -> Result<Self> {
        // x=r+d₀/z and y=d₀*w/z³. Multiply the two linear forms by z.
        let numerator = [chart.derivative.clone(), chart.branch.clone()];
        let denominator = [field.zero(), field.one()];
        let substitute = |model: &mut Model, form: &Polynomial| {
            homogeneous_substitution(model, field, form, 1, &numerator, &denominator)
        };
        let numerator = substitute(model, &self.numerator)?;
        let denominator = substitute(model, &self.denominator)?;
        let ordinate_scale = field.multiply(model, &self.ordinate_scale, &chart.derivative)?;
        model.allocate(8 + (numerator.len() + denominator.len() + 3) * field.cells())?;
        Ok(Self::new(
            numerator,
            denominator,
            self.offset.clone(),
            self.scale.clone(),
            ordinate_scale,
        ))
    }

    fn point(&self, model: &mut Model, target: &Curve, x: &[u64], y: &[u64]) -> Result<Divisor> {
        let field = Field { ring: &target.ring };
        let numerator = poly::evaluate(model, field, &self.numerator, x)?;
        let denominator = poly::evaluate(model, field, &self.denominator, x)?;
        self.image(model, target, &numerator, &denominator, y)
    }

    fn infinity(&self, model: &mut Model, target: &Curve) -> Result<Divisor> {
        let zero = target.ring.scalar(0);
        // At genus-two odd-degree infinity, y/z³ vanishes in the local ring.
        self.image(
            model,
            target,
            self.numerator.get(1).unwrap_or(&zero),
            self.denominator.get(1).unwrap_or(&zero),
            &zero,
        )
    }

    fn image(
        &self,
        model: &mut Model,
        target: &Curve,
        numerator: &Scalar,
        denominator: &Scalar,
        y: &[u64],
    ) -> Result<Divisor> {
        let ring = &target.ring;
        let field = Field { ring };
        if field.is_zero(denominator) {
            return Ok(Divisor::zero(target));
        }
        let inverse = field.inverse_unit(model, denominator)?;
        let ratio = field.multiply(model, numerator, &inverse)?;
        let square = field.multiply(model, &ratio, &ratio)?;
        let image_x = field.multiply(model, &self.scale, &square)?;
        let image_x = field.add(model, &self.offset, &image_x)?;
        let inverse_cube = ring.power(model, &inverse, 3)?;
        let ordinate = ring.multiply(model, &self.ordinate_scale, y)?;
        let image_y = field.multiply(model, &ordinate, &inverse_cube)?;
        Divisor::from_point(model, target, &image_x, image_y)
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
