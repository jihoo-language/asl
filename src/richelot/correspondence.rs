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

mod push;

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

    pub(crate) fn twisted(
        &self,
        model: &mut Model,
        source: Rc<Curve>,
        k: &Scalar,
        inverse_powers: &[Scalar; 6],
    ) -> Result<(Self, Scalar)> {
        let field = Field { ring: &source.ring };
        let scale = |model: &mut Model, factor: &Polynomial| {
            poly::scale_homogeneous(model, field, factor, inverse_powers, factor.len() - 1)
        };
        let [a, b, c] = &self.isogeny.source_factors;
        let source_factors = [scale(model, a)?, scale(model, b)?, scale(model, c)?];
        let (determinant, target_factors) = super::transform(model, field, &source_factors)?
            .ok_or_else(|| {
                Error::plain("Richelot twist has zero determinant despite the nonzero determinant of the original splitting")
            })?;
        let polynomial = super::product(model, field, &target_factors)?;
        let uncharted = Curve::new(
            model,
            source.grade.ring.clone(),
            source.ring.clone(),
            polynomial,
        )?;
        let (chart, target, target_scale) = if let Some(chart) = &self.chart {
            // The branch is at r/k and its derivative is d/k^4. Thus the chart
            // coordinate becomes t/k^3 and its monic twist is by k^3.
            let branch = field.multiply(model, &chart.branch, &inverse_powers[1])?;
            let (chart, polynomial) =
                Chart::build(model, chart.source.clone(), &uncharted, branch)?;
            let target = Curve::new(
                model,
                source.grade.ring.clone(),
                source.ring.clone(),
                polynomial,
            )?;
            (Some(Rc::new(chart)), target, field.ring.power(model, k, 3)?)
        } else {
            (None, uncharted, k.clone())
        };
        let isogeny = Isogeny {
            source: self.isogeny.source.clone(),
            source_factors,
            target_factors,
            determinant,
        };
        Ok((
            Self {
                source,
                target: Rc::new(target),
                chart,
                isogeny: Rc::new(isogeny),
                linear: self.linear,
            },
            target_scale,
        ))
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

        let [i, j] = self.quadratic_indices();
        let g1 = poly::evaluate(model, field, &self.isogeny.source_factors[i], x)?;
        let g2 = poly::evaluate(model, field, &self.isogeny.source_factors[j], x)?;
        if poly::zero(y) {
            return if field.is_zero(&g1) {
                self.branch_pair(model, field, j)
            } else if field.is_zero(&g2) {
                self.branch_pair(model, field, i)
            } else {
                Ok(Divisor::zero(&self.target))
            };
        }
        let y_inverse = field.ring.inverse(model, y)?;
        let ratio = field.multiply(model, &g1, &y_inverse)?;
        let x_ratio = field.ring.multiply(model, x, &ratio)?;
        let (u, v) = self.fiber(model, field, [&g1, &g2], [&x_ratio, &ratio])?;
        self.represent(model, field, u, v)
    }

    pub(crate) fn source_name(&self) -> &str {
        &self.isogeny.source
    }

    fn quadratic_indices(&self) -> [usize; 2] {
        // Rotate the splitting so the last pair contains the source infinity.
        [(self.linear + 1) % 3, (self.linear + 2) % 3]
    }

    fn branch_pair(&self, model: &mut Model, field: Field<'_>, index: usize) -> Result<Divisor> {
        let support = self.isogeny.target_factors[index].clone();
        let (u, v) = self.chart_coordinates(model, field, support, Vec::new())?;
        self.represent(model, field, u, v)
    }

    // One base-field component of L = G1(x)B1 + G2(x)B2 and
    // w = B1 · ((x G1(x)/y) − (G1(x)/y) z).
    fn fiber(
        &self,
        model: &mut Model,
        field: Field<'_>,
        [g1, g2]: [&Scalar; 2],
        [x_ratio, ratio]: [&Scalar; 2],
    ) -> Result<(Polynomial, Polynomial)> {
        let [i, j] = self.quadratic_indices();
        let stored = &self.isogeny.target_factors;
        let raw_first = poly::scale(model, field, &stored[0], &self.isogeny.determinant)?;
        let brackets = [raw_first.as_slice(), &stored[1], &stored[2]];
        let left = poly::scale(model, field, brackets[i], g1)?;
        let right = poly::scale(model, field, brackets[j], g2)?;
        let support = poly::sum(model, field, &left, &right, false)?;
        // On the fiber, the cofactor identity gives w² = f_target(z) for
        // w = G1(x) B1(z) (x-z) / y.
        let negated_ratio = field.ring.negate(model, ratio)?;
        let ordinate =
            poly::multiply(model, field, brackets[i], &[x_ratio.clone(), negated_ratio])?;
        self.chart_coordinates(model, field, support, ordinate)
    }

    fn chart_coordinates(
        &self,
        model: &mut Model,
        field: Field<'_>,
        support: Polynomial,
        ordinate: Polynomial,
    ) -> Result<(Polynomial, Polynomial)> {
        if let Some(chart) = &self.chart {
            let numerator = poly::trim(vec![chart.derivative.clone(), chart.branch.clone()]);
            let denominator = [field.zero(), field.one()];
            // The fiber over the source infinity is a hyperelliptic fiber.
            // Transform both projective points before discarding any at the
            // new infinity, which is the base point in the target chart.
            let u = homogeneous_substitution(model, field, &support, 2, &numerator, &denominator)?;
            let v = homogeneous_substitution(model, field, &ordinate, 3, &numerator, &denominator)?;
            let inverse = field.inverse_unit(model, &chart.derivative)?;
            Ok((u, poly::scale(model, field, &v, &inverse)?))
        } else {
            Ok((support, ordinate))
        }
    }

    fn represent(
        &self,
        model: &mut Model,
        field: Field<'_>,
        u: Polynomial,
        v: Polynomial,
    ) -> Result<Divisor> {
        let leading = u
            .last()
            .ok_or_else(|| Error::plain("Richelot fiber support vanished"))?;
        let inverse = field.inverse_unit(model, leading)?;
        let u = poly::scale(model, field, &u, &inverse)?;
        Divisor::new(model, &self.target, u, v)
    }
}
