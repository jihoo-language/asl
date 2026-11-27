use std::fmt::Write;

use crate::ec_poly;
use crate::error::{Error, Result};
use crate::hw_closed::{Curve, CurveOp};
use crate::model::{Model, Value};
use crate::rv_ore::gcd;

use super::{picard_ranks, powered_weil_coefficients};

pub(super) fn report(
    model: &mut Model,
    curve: &Curve,
    name: &str,
    q: u64,
    weil: &[i128],
    orbits: &[usize],
    degree: u128,
) -> Result<String> {
    let exponent = u32::try_from(degree).map_err(|_| overflow())?;
    let field_order = i128::from(q).checked_pow(exponent).ok_or_else(overflow)?;
    let (a, b) = powered_weil_coefficients(model, weil, exponent as usize, overflow)?;
    let abelian = abelian_factor(field_order, a, b)?;
    let product = abelian
        .checked_mul(exceptional_factor(model, orbits, degree)?)
        .ok_or_else(overflow)?;
    let (rank, geometric_rank) = picard_ranks(model, q, weil, orbits, degree)?;
    let sign = if rank.is_multiple_of(2) { "-" } else { "" };
    let mut result = format!(
        "rvKummerK3ArtinTate[{name}, {degree}] = field order={field_order}, rank={rank}, geometric rank={geometric_rank}, finite NS discriminant times Brauer order={sign}{product}, finite NS discriminant square class=[{sign}{product}]"
    );
    if geometric_rank == 22 {
        let sigma = artin_invariant(model, curve)?;
        let discriminant = u128::from(curve.ring.prime())
            .checked_pow(2 * sigma)
            .ok_or_else(overflow)?;
        write!(
            result,
            ", Artin invariant={sigma}, geometric NS discriminant=-{discriminant}"
        )
        .unwrap();
        if rank == 22 {
            let brauer_order = brauer_order(product, discriminant)?;
            write!(result, ", Brauer order={brauer_order}").unwrap();
        }
    }
    Ok(result)
}

fn overflow() -> Error {
    Error::plain("Kummer K3 Artin–Tate computation exceeds the checked 128-bit arithmetic range")
}

fn abelian_factor(field_order: i128, a: i128, b: i128) -> Result<u128> {
    // The two unavoidable eigenvalues Q are already removed. The remaining
    // reciprocal quartic has multiplicity zero, two, or four at T=1/Q.
    let compute = || {
        let four_q = field_order.checked_mul(4)?;
        let delta = a
            .checked_mul(a)?
            .checked_sub(b.checked_mul(4)?)?
            .checked_add(four_q.checked_mul(2)?)?;
        if delta != 0 {
            return Some(delta);
        }
        let u = a / 2;
        let remainder = four_q.checked_sub(u.checked_mul(u)?)?;
        Some(if remainder != 0 {
            remainder
        } else {
            field_order
        })
    };
    let factor = compute().ok_or_else(overflow)?;
    if factor <= 0 {
        return Err(Error::plain(
            "Kummer K3 Artin–Tate abelian factor must be positive",
        ));
    }
    Ok(factor as u128)
}

fn exceptional_factor(model: &mut Model, orbits: &[usize], degree: u128) -> Result<u128> {
    let mut product = 1u128;
    for &length in orbits {
        let length = length as u128;
        let cycles = gcd(length, degree);
        model.spend(1 + cycles as usize)?;
        let factor = (length / cycles)
            .checked_pow(cycles as u32)
            .ok_or_else(overflow)?;
        product = product.checked_mul(factor).ok_or_else(overflow)?;
    }
    Ok(product)
}

fn artin_invariant(model: &mut Model, curve: &Curve) -> Result<u32> {
    let Value::HWWMatrix(cartier) = curve.evaluate(model, CurveOp::Alpha(1))? else {
        unreachable!("the first Hasse–Witt matrix is the Cartier–Manin matrix up to transpose");
    };
    model.spend(cartier.entries().len() * curve.ring.degree())?;
    Ok(
        if cartier
            .entries()
            .iter()
            .all(|coefficient| ec_poly::zero(coefficient))
        {
            1
        } else {
            2
        },
    )
}

fn brauer_order(product: u128, discriminant: u128) -> Result<u128> {
    if product.is_multiple_of(discriminant) {
        let order = product / discriminant;
        if order > 0 && order.isqrt().pow(2) == order {
            return Ok(order);
        }
    }
    Err(Error::plain(
        "Kummer K3 rank-22 Artin–Tate product does not give a positive square Brauer order",
    ))
}
