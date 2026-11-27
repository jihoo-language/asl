use crate::error::{Error, Result};
use crate::model::Model;

use super::super::Magnitude;
use super::inverse::{self, shifted_degree, Chamber};
use super::{eta_product, multiply, overflow, Cutoff, Integer};

const ZERO: Integer = Integer {
    magnitude: Magnitude::ZERO,
    negative: false,
};

pub(super) fn polar(
    model: &mut Model,
    m: i128,
    n: i128,
    r: i128,
    chamber: Chamber,
) -> Result<Integer> {
    let index = positive_index(m)?;
    let Some(degree) = shifted_degree(n) else {
        return Ok(ZERO);
    };
    let on_side = match chamber {
        Chamber::Inside => r > 0,
        Chamber::Outside => r < 0,
    };
    let radius = r.unsigned_abs();
    let mut total = if on_side {
        let t = Integer {
            magnitude: Magnitude::from_small(radius),
            negative: false,
        };
        inverse::wall_jump(model, m, n, t)?
    } else {
        ZERO
    };

    // For s ≠ 0 the sign of t=r−2ms must agree with s. After using
    // y-inversion this gives 1 ≤ |s| and 2m|s| < |r|. The q shift
    // |s|(|r|−m|s|) strictly increases on this interval.
    let bound = radius.saturating_sub(1) / (2 * index);
    for s in 1..=bound {
        model.spend(1)?;
        let offset = index * s;
        let Some(shift) = s.checked_mul(radius - offset) else {
            break;
        };
        if shift > degree {
            break;
        }
        let n = unshift(degree - shift);
        let t = Integer {
            magnitude: Magnitude::from_small(radius - 2 * offset),
            negative: false,
        };
        let term = inverse::wall_jump(model, m, n, t)?;
        total = total.add(term).map_err(|_| overflow())?;
    }
    Ok(total)
}

pub(super) fn finite(model: &mut Model, m: i128, n: i128, r: i128) -> Result<Integer> {
    let Some((n, r)) = reduce(model, m, n, r)? else {
        return Ok(ZERO);
    };
    // In 0 ≤ r ≤ m all nonzero Appell–Lerch translates vanish. The
    // s=0 term also vanishes in the outside expansion at this r.
    finite_row(model, m, n, r)
}

pub(super) fn normalized(model: &mut Model, m: i128, n: i128, r: i128) -> Result<Integer> {
    let Some((n, r)) = reduce(model, m, n, r)? else {
        return Ok(ZERO);
    };
    if n < 0 {
        return Ok(ZERO);
    }
    let cutoff = Cutoff::new(model, m as u128 + 1, n as u128)?;
    let eta = eta_product(model, cutoff.q, 24)?;
    let mut sum = ZERO;
    // Δ=q∏(1−qʲ)^24. Descending finite rows reuse the inverse-product
    // rectangle populated by the first query.
    for (j, tau) in eta.into_iter().enumerate() {
        if !tau.is_zero() {
            let coefficient = finite_row(model, m, n - j as i128 - 1, r)?;
            let term = multiply(model, coefficient, tau)?;
            sum = sum.add(term).map_err(|_| overflow())?;
        }
    }
    Ok(sum)
}

fn finite_row(model: &mut Model, m: i128, n: i128, r: i128) -> Result<Integer> {
    inverse::coefficient(model, m, n, Integer::from_i128(r), Chamber::Outside)
}

fn positive_index(m: i128) -> Result<u128> {
    if m <= 0 {
        return Err(Error::plain("Igusa mock Jacobi index m must be positive"));
    }
    Ok(m as u128)
}

fn reduce(model: &mut Model, m: i128, n: i128, r: i128) -> Result<Option<(i128, i128)>> {
    let m = positive_index(m)?;
    let Some(degree) = shifted_degree(n) else {
        return Ok(None);
    };
    // Evenness in r lets us use the nonnegative representative of the
    // centered residue. Preserve the discriminant 4mn−r² during reduction.
    let radius = r.unsigned_abs();
    let remainder = radius % (2 * m);
    let r = remainder.min(2 * m - remainder);
    let difference = Magnitude::from_small(radius - r)
        .multiply_u128(radius + r)
        .ok_or_else(overflow)?;
    let (divided, remainder) = difference.divide_u32(4);
    debug_assert_eq!(remainder, 0);
    let (shift, remainder) = divided.divide(model, m)?;
    debug_assert_eq!(remainder, 0);
    let reduced_n = shift.small().and_then(|shift| degree.checked_sub(shift));
    Ok(reduced_n.map(|n| (unshift(n), r as i128)))
}

fn unshift(degree: u128) -> i128 {
    if degree == 0 {
        -1
    } else {
        (degree - 1) as i128
    }
}
