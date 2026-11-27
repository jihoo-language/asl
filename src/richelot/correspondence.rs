use std::rc::Rc;

use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic_chart::{homogeneous_substitution, Chart};
use crate::hyperelliptic_jacobian::{self, Divisor};
use crate::model::Model;
use crate::syntax::Label;

use super::Isogeny;

pub(crate) struct Correspondence {
    pub(crate) source: Rc<Curve>,
    pub(crate) target: Rc<Curve>,
    chart: Option<Rc<Chart>>,
    isogeny: Rc<Isogeny>,
    linear: usize,
}

impl Correspondence {
    pub(crate) fn resolve(model: &Model, target: &Label) -> Result<Self> {
        let curve = model
            .jacobian_curve(&target.text)
            .map_err(|error| error.at(target.span))?;
        let chart = model.jacobian_chart(target).ok();
        let isogeny = if let Some(chart) = &chart {
            model.richelot(&Label {
                text: chart.source.to_string(),
                span: target.span,
            })?
        } else {
            model.richelot(target)?
        };
        let source = model.hww_curve(&Label {
            text: isogeny.source.clone(),
            span: target.span,
        })?;
        if source.polynomial().len() != 6 || hyperelliptic_jacobian::require_curve(&source).is_err()
        {
            return Err(Error::new(
                target.span,
                "Richelot correspondence requires a monic quintic source",
            ));
        }
        let linear = isogeny
            .source_factors
            .iter()
            .position(|factor| factor.len() == 2)
            .expect("a quintic splitting has one linear factor");
        Ok(Self {
            source,
            target: curve,
            chart,
            isogeny,
            linear,
        })
    }

    pub(crate) fn point(&self, model: &mut Model, x: &[u64], y: &[u64]) -> Result<Divisor> {
        let field = Field {
            ring: &self.source.ring,
        };
        let f_at_x = poly::evaluate(model, field, self.source.polynomial(), x)?;
        if field.ring.multiply(model, y, y)? != f_at_x {
            return Err(Error::plain(
                "point does not lie on the Richelot source curve",
            ));
        }

        // Rotate the splitting so the last pair contains the source infinity.
        let i = (self.linear + 1) % 3;
        let j = (self.linear + 2) % 3;
        let g1 = poly::evaluate(model, field, &self.isogeny.source_factors[i], x)?;
        let g2 = poly::evaluate(model, field, &self.isogeny.source_factors[j], x)?;
        let stored = &self.isogeny.target_factors;
        if poly::zero(y) {
            let support = if field.is_zero(&g1) {
                &stored[j]
            } else if field.is_zero(&g2) {
                &stored[i]
            } else {
                return Ok(Divisor::zero(&self.target));
            };
            return self.represent(model, field, support, &[]);
        }

        let raw_first = poly::scale(model, field, &stored[0], &self.isogeny.determinant)?;
        let brackets = [raw_first.as_slice(), &stored[1], &stored[2]];
        let left = poly::scale(model, field, brackets[i], &g1)?;
        let right = poly::scale(model, field, brackets[j], &g2)?;
        let support = poly::sum(model, field, &left, &right, false)?;
        // On the fiber, the cofactor identity gives w² = f_target(z) for
        // w = G1(x) B1(z) (x-z) / y.
        let ordinate = poly::multiply(model, field, &left, &[x.to_vec(), field.signed(-1)])?;
        let y_inverse = field.ring.inverse(model, y)?;
        let ordinate = poly::scale(model, field, &ordinate, &y_inverse)?;
        self.represent(model, field, &support, &ordinate)
    }

    fn represent(
        &self,
        model: &mut Model,
        field: Field<'_>,
        support: &[Scalar],
        ordinate: &[Scalar],
    ) -> Result<Divisor> {
        let (u, v): (Polynomial, Polynomial) = if let Some(chart) = &self.chart {
            let numerator = poly::trim(vec![chart.derivative.clone(), chart.branch.clone()]);
            let denominator = [field.zero(), field.one()];
            // The fiber over the source infinity is a hyperelliptic fiber.
            // Transform both projective points before discarding any at the
            // new infinity, which is the base point in the target chart.
            let u = homogeneous_substitution(model, field, support, 2, &numerator, &denominator)?;
            let v = homogeneous_substitution(model, field, ordinate, 3, &numerator, &denominator)?;
            let inverse = field.inverse_unit(model, &chart.derivative)?;
            (u, poly::scale(model, field, &v, &inverse)?)
        } else {
            (support.to_vec(), ordinate.to_vec())
        };
        let leading = u
            .last()
            .ok_or_else(|| Error::plain("Richelot fiber support vanished"))?;
        let inverse = field.inverse_unit(model, leading)?;
        let u = poly::scale(model, field, &u, &inverse)?;
        Divisor::new(model, &self.target, u, v)
    }
}
