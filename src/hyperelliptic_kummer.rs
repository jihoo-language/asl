use std::collections::BTreeMap;

use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::Result;
use crate::finite_field_factor::{distinct_degrees, quadratic_sign};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic_jacobian::Divisor;
use crate::model::Model;

pub(crate) struct Support {
    pub(crate) polynomial: Polynomial,
    pub(crate) counts: BTreeMap<usize, usize>,
}

impl Support {
    pub(crate) fn is_trivial(&self) -> bool {
        self.polynomial.len() == 1
    }

    pub(crate) fn into_coefficient(self, curve: &Curve, index: usize) -> Scalar {
        self.polynomial
            .into_iter()
            .nth(index)
            .unwrap_or_else(|| curve.ring.scalar(0))
    }

    fn collect(
        &mut self,
        model: &mut Model,
        field: Field<'_>,
        q: Option<u64>,
        remaining: Polynomial,
        unit: &[Scalar],
    ) -> Result<()> {
        distinct_degrees(model, field, q, remaining, |model, block, degree| {
            self.collect_degree(model, field, q, &block, degree, unit)
        })
    }

    fn collect_degree(
        &mut self,
        model: &mut Model,
        field: Field<'_>,
        q: Option<u64>,
        block: &[Scalar],
        degree: usize,
        unit: &[Scalar],
    ) -> Result<()> {
        let sign = quadratic_sign(model, field, q, block, degree, unit)?;
        let sign_plus_one = poly::sum(model, field, &sign, &[field.one()], false)?;
        model.spend(block.len() * field.cells())?;
        let nonsquares = poly::monic_gcd(model, field, block.to_vec(), sign_plus_one)?;
        if nonsquares.len() > 1 {
            model.spend(1)?;
            *self.counts.entry(degree).or_default() += (nonsquares.len() - 1) / degree;
            self.polynomial = poly::multiply(model, field, &self.polynomial, &nonsquares)?;
        }
        Ok(())
    }
}

pub(crate) fn support(model: &mut Model, curve: &Curve, divisor: &Divisor) -> Result<Support> {
    let field = Field { ring: &curve.ring };
    let mut result = Support {
        polynomial: vec![field.one()],
        counts: BTreeMap::new(),
    };
    let degree = divisor.u.len() - 1;
    if degree == 0 {
        return Ok(result);
    }
    let f = curve.polynomial();
    field.check_polynomial(f.len())?;
    model.spend((f.len() + divisor.u.len()) * field.cells())?;
    let common = poly::monic_gcd(model, field, f.to_vec(), divisor.u.clone())?;
    let other = poly::exact_divide(model, field, f, &common)?;
    let q = field.ring.size_of_residue_field().ok();

    let even_degree = degree.is_multiple_of(2);
    if common.len() > 1 {
        // At a shared branch factor, replace the vanishing x−theta by its local
        // unit 1/f'(theta). Cancelling the common factors gives u0/f0 there,
        // with sign (-1)^(deg(u)−1). In square classes 1/f0 equals f0.
        let u0 = poly::exact_divide(model, field, &divisor.u, &common)?;
        let mut unit = poly::multiply_mod(model, field, &u0, &other, &common)?;
        if even_degree {
            unit = poly::scale(model, field, &unit, &field.signed(-1))?;
        }
        result.collect(model, field, q, common, &unit)?;
    }
    // Away from the common branch factors the descent map is (-1)^deg(u) u.
    if even_degree {
        result.collect(model, field, q, other, &divisor.u)?;
    } else {
        let unit = poly::scale(model, field, &divisor.u, &field.signed(-1))?;
        result.collect(model, field, q, other, &unit)?;
    }
    Ok(result)
}
