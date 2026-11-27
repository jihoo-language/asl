use std::collections::BTreeMap;

use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::Result;
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
        mut remaining: Polynomial,
        unit: &[Scalar],
    ) -> Result<()> {
        let x = field.polynomial(&[0, 1]);
        let mut x_power = x.clone();
        let mut degree = 1;
        while degree <= (remaining.len() - 1) / 2 {
            x_power = frobenius(model, field, q, x_power, &remaining)?;
            let difference = poly::sum(model, field, &x_power, &x, true)?;
            model.spend(remaining.len() * field.cells())?;
            let block = monic_gcd(model, field, remaining.clone(), difference)?;
            if block.len() > 1 {
                self.collect_degree(model, field, q, &block, degree, unit)?;
                remaining = poly::exact_divide(model, field, &remaining, &block)?;
                x_power = poly::remainder(model, field, x_power, &remaining)?;
            }
            degree += 1;
        }
        // All remaining factors have degree at least `degree`, so a nonconstant
        // remainder here must be irreducible.
        if remaining.len() > 1 {
            self.collect_degree(model, field, q, &remaining, remaining.len() - 1, unit)?;
        }
        Ok(())
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
        let nonsquares = monic_gcd(model, field, block.to_vec(), sign_plus_one)?;
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
    let common = monic_gcd(model, field, f.to_vec(), divisor.u.clone())?;
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

fn monic_gcd(
    model: &mut Model,
    field: Field<'_>,
    left: Polynomial,
    right: Polynomial,
) -> Result<Polynomial> {
    let gcd = poly::gcd(model, field, left, right)?;
    let leading = gcd.last().expect("one input is a monic branch polynomial");
    if *leading == field.one() {
        Ok(gcd)
    } else {
        let inverse = field.inverse_unit(model, leading)?;
        poly::scale(model, field, &gcd, &inverse)
    }
}

fn frobenius(
    model: &mut Model,
    field: Field<'_>,
    q: Option<u64>,
    mut value: Polynomial,
    modulus: &[Scalar],
) -> Result<Polynomial> {
    if let Some(q) = q {
        return poly::power_mod(model, field, &value, q, modulus);
    }
    for _ in 0..field.ring.degree() {
        value = poly::power_mod(model, field, &value, field.ring.prime(), modulus)?;
    }
    Ok(value)
}

fn quadratic_sign(
    model: &mut Model,
    field: Field<'_>,
    q: Option<u64>,
    block: &[Scalar],
    degree: usize,
    unit: &[Scalar],
) -> Result<Polynomial> {
    if let Some(size) = q.and_then(|q| q.checked_pow(degree as u32)) {
        return poly::power_mod(model, field, unit, (size - 1) / 2, block);
    }
    // The base-p digits of (p^(e*d)−1)/2 all equal (p−1)/2.
    let p = field.ring.prime();
    let digit_power = poly::power_mod(model, field, unit, (p - 1) / 2, block)?;
    model.spend(digit_power.len() * field.cells())?;
    let mut sign = digit_power.clone();
    for _ in 1..field.ring.degree() * degree {
        let raised = poly::power_mod(model, field, &sign, p, block)?;
        sign = poly::multiply_mod(model, field, &raised, &digit_power, block)?;
    }
    Ok(sign)
}
