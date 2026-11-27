use super::{field, require_curve, Divisor};
use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::model::Model;

/// Trace a Mumford pair over K[ω]/(ω²−k) to the Jacobian over K.
/// Each pair of polynomials gives its real and imaginary components in that order.
pub(crate) fn trace(
    model: &mut Model,
    curve: &Curve,
    k: &Scalar,
    u: [Polynomial; 2],
    v: [Polynomial; 2],
) -> Result<Divisor> {
    require_curve(curve)?;
    let field = field(curve);
    if field.is_zero(k) {
        return Err(Error::plain("rvJacQuadTrace radicand must be nonzero"));
    }
    let [a, b] = u.map(poly::trim);
    if a.last() != Some(&field.one()) || b.len() >= a.len() {
        return Err(Error::plain(
            "rvJacQuadTrace relative Mumford U must be nonzero and monic",
        ));
    }
    if a.len() == 1 {
        return Ok(Divisor::zero(curve));
    }
    let [p, q] = relative_remainder(model, field, k, v.map(poly::trim), &a, &b)?;
    let p_squared = poly::multiply(model, field, &p, &p)?;
    let q_squared = poly::multiply(model, field, &q, &q)?;
    let k_q_squared = poly::scale(model, field, &q_squared, k)?;
    let square = poly::sum(model, field, &p_squared, &k_q_squared, false)?;
    let real_norm = poly::sum(model, field, curve.polynomial(), &square, true)?;
    let pq = poly::multiply(model, field, &p, &q)?;
    let imaginary_norm = poly::scale(model, field, &pq, &field.signed(-2))?;
    let norm = relative_remainder(model, field, k, [real_norm, imaginary_norm], &a, &b)?;
    if norm.iter().any(|part| !part.is_empty()) {
        return Err(Error::plain("rvJacQuadTrace U must divide f−V²"));
    }

    model.spend((a.len() + b.len()) * field.cells())?;
    let common = poly::monic_gcd(model, field, a.clone(), b.clone())?;
    let mut shared = Divisor::zero(curve);
    if common.len() > 1 {
        // On shared support PQ vanishes modulo G. If P vanishes and Q is a unit,
        // the two lifts are opposite to full multiplicity. If both vanish, the
        // place is ramified and simple. Where P is a unit, Q vanishes and they agree.
        model.spend((common.len() + p.len()) * field.cells())?;
        let opposite = poly::monic_gcd(model, field, common.clone(), p.clone())?;
        if opposite.len() < common.len() {
            let agreeing = poly::exact_divide(model, field, &common, &opposite)?;
            model.spend(p.len() * field.cells())?;
            let interpolation = poly::remainder(model, field, p.clone(), &agreeing)?;
            shared = Divisor::new(model, curve, agreeing, interpolation)?;
            shared = shared.add(model, curve, &shared)?;
        }
    }
    if a.len() == common.len() {
        return Ok(shared);
    }

    let (a, b) = if common.len() > 1 {
        (
            poly::exact_divide(model, field, &a, &common)?,
            poly::exact_divide(model, field, &b, &common)?,
        )
    } else {
        (a, b)
    };
    // The remaining support and its conjugate are coprime. Their union is cut
    // out by M = A²−kB², and ω = −A/B on U = 0 since B is a unit modulo M.
    let a_squared = poly::multiply(model, field, &a, &a)?;
    let b_squared = poly::multiply(model, field, &b, &b)?;
    let k_b_squared = poly::scale(model, field, &b_squared, k)?;
    let modulus = poly::sum(model, field, &a_squared, &k_b_squared, true)?;
    let inverse_b = poly::inverse_mod(model, field, &b, &modulus)?;
    let a_over_b = poly::multiply_mod(model, field, &a, &inverse_b, &modulus)?;
    let twist = poly::multiply_mod(model, field, &q, &a_over_b, &modulus)?;
    let interpolation = poly::sum(model, field, &p, &twist, true)?;
    let residual = Divisor::new(model, curve, modulus, interpolation)?;
    residual.add(model, curve, &shared)
}

// Monic long division over K[ω]/(ω²−k) requires no inversions in that algebra.
fn relative_remainder(
    model: &mut Model,
    field: Field<'_>,
    k: &Scalar,
    [mut real, mut imaginary]: [Polynomial; 2],
    a: &[Scalar],
    b: &[Scalar],
) -> Result<[Polynomial; 2]> {
    let degree = a.len() - 1;
    let width = real.len().max(imaginary.len());
    if width <= degree {
        return Ok([real, imaginary]);
    }
    field.check_polynomial(width.saturating_mul(2))?;
    model.spend(width.saturating_mul(2))?;
    let zero = field.zero();
    real.resize(width, zero.clone());
    imaginary.resize(width, zero.clone());
    for index in (degree..width).rev() {
        let real_top = real.pop().expect("dividend has a coefficient at index");
        let imaginary_top = imaginary
            .pop()
            .expect("dividend has a coefficient at index");
        if field.is_zero(&real_top) && field.is_zero(&imaginary_top) {
            continue;
        }
        let k_imaginary_top = field.multiply(model, k, &imaginary_top)?;
        let order = index - degree;
        model.spend(degree)?;
        for (offset, a_coefficient) in a[..degree].iter().enumerate() {
            let b_coefficient = b.get(offset).unwrap_or(&zero);
            let real_a = field.multiply(model, &real_top, a_coefficient)?;
            let real_b = field.multiply(model, &k_imaginary_top, b_coefficient)?;
            let imaginary_a = field.multiply(model, &imaginary_top, a_coefficient)?;
            let imaginary_b = field.multiply(model, &real_top, b_coefficient)?;
            let real_product = field.add(model, &real_a, &real_b)?;
            let imaginary_product = field.add(model, &imaginary_a, &imaginary_b)?;
            let j = order + offset;
            real[j] = field.subtract(model, &real[j], &real_product)?;
            imaginary[j] = field.subtract(model, &imaginary[j], &imaginary_product)?;
        }
    }
    Ok([poly::trim(real), poly::trim(imaginary)])
}
