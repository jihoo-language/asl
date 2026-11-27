use std::rc::Rc;

use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic_chart::Chart;
use crate::hyperelliptic_jacobian::{quadratic_trace, Divisor};
use crate::model::Model;
use crate::syntax::Label;

use super::Coordinates;

#[derive(Clone, Copy)]
pub(crate) enum Direction {
    Push,
    Pull,
}

impl Direction {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Push => "rvBiellipticPush",
            Self::Pull => "rvBiellipticPull",
        }
    }
}

pub(crate) struct Map {
    pub(crate) source: Rc<Curve>,
    pub(crate) target: Rc<Curve>,
    coordinates: Rc<Coordinates>,
}

impl Map {
    pub(crate) fn resolve(model: &mut Model, source: &Label, target: &Label) -> Result<Self> {
        let quotient = model.bielliptic_quotient(target)?;
        let curve = model
            .jacobian_curve(&source.text)
            .map_err(|error| error.at(source.span))?;
        if Rc::ptr_eq(&curve, &quotient.source) {
            return Ok(Self {
                source: curve,
                target: quotient.target.clone(),
                coordinates: quotient.coordinates.clone(),
            });
        }
        let wrong_source = || {
            Error::new(
                source.span,
                format!(
                    "RVJac({}) must be on the named bielliptic source of '{}' or a direct branch chart of that source",
                    source.text, target.text,
                ),
            )
        };
        let chart = model.jacobian_chart(source).map_err(|_| wrong_source())?;
        let original = model.hww_curve(&Label {
            text: chart.source.to_string(),
            span: source.span,
        })?;
        if !Rc::ptr_eq(&original, &quotient.source) {
            return Err(wrong_source());
        }
        Self::from_chart(
            model,
            curve,
            quotient.target.clone(),
            &quotient.coordinates,
            &chart,
        )
    }

    pub(crate) fn from_chart(
        model: &mut Model,
        source: Rc<Curve>,
        target: Rc<Curve>,
        coordinates: &Coordinates,
        chart: &Chart,
    ) -> Result<Self> {
        let field = Field { ring: &source.ring };
        let coordinates = Rc::new(coordinates.rechart(model, field, chart)?);
        Ok(Self {
            source,
            target,
            coordinates,
        })
    }

    pub(crate) fn push(&self, model: &mut Model, divisor: &Divisor) -> Result<Divisor> {
        let field = Field {
            ring: &self.target.ring,
        };
        let degree = divisor.u.len() - 1;
        if degree == 0 {
            return Ok(Divisor::zero(&self.target));
        }
        // All multiplicities above the pole contribute the elliptic infinity.
        let denominator = &self.coordinates.denominator;
        let support = if denominator.len() == 2 {
            let square = poly::multiply(model, field, denominator, denominator)?;
            let poles = poly::monic_gcd(model, field, divisor.u.clone(), square)?;
            poly::exact_divide(model, field, &divisor.u, &poles)?
        } else {
            divisor.u.clone()
        };
        let image = match support.as_slice() {
            [_] => Divisor::zero(&self.target),
            [constant, _] => {
                let x = field.ring.negate(model, constant)?;
                let y = poly::evaluate(model, field, &divisor.v, &x)?;
                self.coordinates.point(model, &self.target, &x, &y)?
            }
            [constant, linear, _] => {
                // Write u(x)=(x-center)²-k, using -1/2=(p-1)/2.
                let center = field
                    .ring
                    .scale(model, linear, (field.ring.prime() - 1) / 2)?;
                let square = field.multiply(model, &center, &center)?;
                let k = field.subtract(model, &square, constant)?;
                if field.is_zero(&k) {
                    let y = poly::evaluate(model, field, &divisor.v, &center)?;
                    let point = self.coordinates.point(model, &self.target, &center, &y)?;
                    point.add(model, &self.target, &point)?
                } else {
                    self.trace_image(model, field, &support, &divisor.v, &center, &k)?
                }
            }
            _ => unreachable!("reduced genus-two support"),
        };
        // Subtract deg(u) times the source base point's image, which is killed by 2.
        if degree.is_multiple_of(2) {
            return Ok(image);
        }
        let base = self.coordinates.infinity(model, &self.target)?;
        image.add(model, &self.target, &base)
    }

    fn trace_image(
        &self,
        model: &mut Model,
        field: Field<'_>,
        u: &[Scalar],
        v: &[Scalar],
        center: &Scalar,
        k: &Scalar,
    ) -> Result<Divisor> {
        let coordinates = &self.coordinates;
        let inverse = poly::inverse_mod(model, field, &coordinates.denominator, u)?;
        let ratio = poly::multiply_mod(model, field, &coordinates.numerator, &inverse, u)?;
        let x = poly::multiply_mod(model, field, &ratio, &ratio, u)?;
        let x = poly::scale(model, field, &x, &coordinates.scale)?;
        let inverse_cube = poly::power_mod(model, field, &inverse, 3, u)?;
        let y = poly::multiply_mod(model, field, v, &inverse_cube, u)?;
        let y = poly::scale(model, field, &y, &coordinates.ordinate_scale)?;

        // Evaluate the constant part at x=center and put the linear part on ω.
        let split = |model: &mut Model, coordinate: Polynomial| -> Result<[Scalar; 2]> {
            let real = poly::evaluate(model, field, &coordinate, center)?;
            let imaginary = coordinate.get(1).cloned().unwrap_or_else(|| field.zero());
            Ok([real, imaginary])
        };
        let [x0, x1] = split(model, x)?;
        let x0 = field.add(model, &x0, &coordinates.offset)?;
        let [y0, y1] = split(model, y)?;
        let u = [
            vec![field.ring.negate(model, &x0)?, field.one()],
            poly::trim(vec![field.ring.negate(model, &x1)?]),
        ];
        quadratic_trace::trace(model, &self.target, k, u, [vec![y0], vec![y1]])
    }

    pub(crate) fn pull(&self, model: &mut Model, divisor: &Divisor) -> Result<Divisor> {
        if divisor.u.len() == 1 {
            return Ok(Divisor::zero(&self.source));
        }
        let field = Field {
            ring: &self.source.ring,
        };
        let x = field.ring.negate(model, &divisor.u[0])?;
        let y = divisor.coefficient(&self.target, true, 0);
        let coordinates = &self.coordinates;
        let x = field.subtract(model, &x, &coordinates.offset)?;
        let n = &coordinates.numerator;
        let m = &coordinates.denominator;
        let n_squared = poly::multiply(model, field, n, n)?;
        let m_squared = poly::multiply(model, field, m, m)?;
        let left = poly::scale(model, field, &n_squared, &coordinates.scale)?;
        let right = poly::scale(model, field, &m_squared, &x)?;
        let support = poly::sum(model, field, &left, &right, true)?;
        let support = poly::monic(model, field, support)?;
        let inverse = field.inverse_unit(model, &coordinates.ordinate_scale)?;
        let scale = field.multiply(model, &y, &inverse)?;
        let m_cubed = poly::multiply(model, field, &m_squared, m)?;
        let ordinate = poly::scale(model, field, &m_cubed, &scale)?;
        // The degree drop omits the source infinity, and the pole fiber is
        // linearly equivalent to twice that infinity.
        Divisor::new(model, &self.source, support, ordinate)
    }
}
