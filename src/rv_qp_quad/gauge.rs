use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;

use super::{normalize_wide, Element, Field, Wide};

// A bracket numerator is a sum of signed terms with up to five i128
// factors and occupies at most 638 bits. Its norm fits 1408 bits.
const LIMBS: usize = 44;
type Integer = Wide<LIMBS>;

#[derive(Clone, Copy)]
pub(crate) enum Coefficients {
    Rational(Prime),
    Quadratic(Field),
}

impl Coefficients {
    pub(crate) fn prime(self) -> Prime {
        match self {
            Self::Rational(p) => p,
            Self::Quadratic(field) => field.prime(),
        }
    }

    fn radicand(self) -> i128 {
        match self {
            Self::Rational(_) => 0,
            Self::Quadratic(field) => field.d,
        }
    }

    pub(crate) fn order(self, model: &mut Model, x: Element) -> Result<Option<i64>> {
        match self {
            Self::Rational(p) => x.rational_valuation(model, p),
            Self::Quadratic(field) => field.valuation(model, x),
        }
    }

    pub(crate) fn angle(self, model: &mut Model, x: Element, order: i64) -> Result<u32> {
        match self {
            Self::Rational(p) => x.rational_angle(model, p),
            Self::Quadratic(field) => field.angular_component(model, x, order),
        }
    }

    pub(crate) fn bracket(
        self,
        model: &mut Model,
        x: [Element; 2],
        y: [Element; 2],
    ) -> Result<Option<(i64, u32)>> {
        model.spend(1)?;
        let d = self.radicand();
        if let Some((a, b, c)) = bracket_small(d, x, y) {
            let value = Element::new(a, b, c)?;
            let Some(v) = self.order(model, value)? else {
                return Ok(None);
            };
            return Ok(Some((v, self.angle(model, value, v)?)));
        }
        let (real, imaginary, denominator) = bracket_wide(model, d, x, y)?;
        if real.is_zero() && imaginary.is_zero() {
            return Ok(None);
        }
        let p = self.prime();
        let (v, unit) = match self {
            Self::Rational(_) => real.prime_part(model, p)?,
            Self::Quadratic(field) => field
                .wide_local(model, real, imaginary)?
                .ok_or_else(|| Error::plain("nonzero bracket lost its quadratic numerator"))?,
        };
        let (e, denom_unit) = denominator.prime_part(model, p)?;
        let denom_inverse = invert(model, p, denom_unit)?;
        Ok(Some((v - e, p.mul(unit, denom_inverse))))
    }

    pub(crate) fn ratio(self, model: &mut Model, x: Element, y: Element) -> Result<Element> {
        if y.is_zero() {
            return Err(Error::plain("residue chart division by zero"));
        }
        model.spend(1)?;
        if x == y {
            return Ok(Element::ONE);
        }
        if y == Element::ONE || x.is_zero() {
            return Ok(x);
        }
        let d = self.radicand();
        if let Some((a, b, c)) = ratio_small(d, x, y) {
            return Element::new(a, b, c);
        }
        // Multiply by the conjugate before normalizing the quotient.
        // The inverse on its own may have an oversized denominator.
        model.spend(6 * 16 + 15 * LIMBS)?;
        let real = Integer::product(x.a, y.a)
            .add(Integer::product(x.b, y.b).multiply(d)?.negate())?
            .multiply(y.c)?;
        let imaginary = Integer::product(x.b, y.a)
            .add(Integer::product(x.a, y.b).negate())?
            .multiply(y.c)?;
        let denominator = Integer::product(y.a, y.a)
            .add(Integer::product(y.b, y.b).multiply(d)?.negate())?
            .multiply(x.c)?;
        normalize_wide(model, real, imaginary, denominator)
    }

    pub(crate) fn subtract_product(
        self,
        model: &mut Model,
        x: Element,
        a: Element,
        y: Element,
    ) -> Result<Element> {
        let d = self.radicand();
        let left = [x, y];
        let right = [a, Element::ONE];
        model.spend(1)?;
        if let Some((a, b, c)) = bracket_small(d, left, right) {
            return Element::new(a, b, c);
        }
        let (a, b, c) = bracket_wide(model, d, left, right)?;
        normalize_wide(model, a, b, c)
    }
}

fn product_small(d: i128, x: Element, y: Element) -> Option<(i128, i128, i128)> {
    let a =
        x.a.checked_mul(y.a)?
            .checked_add(d.checked_mul(x.b.checked_mul(y.b)?)?)?;
    let b = x.a.checked_mul(y.b)?.checked_add(x.b.checked_mul(y.a)?)?;
    let c = x.c.checked_mul(y.c)?;
    Some((a, b, c))
}

fn bracket_small(
    d: i128,
    [a, b]: [Element; 2],
    [c, e]: [Element; 2],
) -> Option<(i128, i128, i128)> {
    let (x, y, m) = product_small(d, a, e)?;
    let (u, v, n) = product_small(d, c, b)?;
    let (left, right) = {
        let g = super::gcd(m as u128, n as u128) as i128;
        (n / g, m / g)
    };
    Some((
        x.checked_mul(left)?.checked_sub(u.checked_mul(right)?)?,
        y.checked_mul(left)?.checked_sub(v.checked_mul(right)?)?,
        m.checked_mul(left)?,
    ))
}

fn product_wide(d: i128, x: Element, y: Element) -> Result<(Integer, Integer)> {
    let a = Integer::product(x.a, y.a).add(Integer::product(x.b, y.b).multiply(d)?)?;
    let b = Integer::product(x.a, y.b).add(Integer::product(x.b, y.a))?;
    Ok((a, b))
}

fn bracket_wide(
    model: &mut Model,
    d: i128,
    [a, b]: [Element; 2],
    [c, e]: [Element; 2],
) -> Result<(Integer, Integer, Integer)> {
    model.spend(9 * 16 + 28 * LIMBS)?;
    let (x, y) = product_wide(d, a, e)?;
    let (u, v) = product_wide(d, c, b)?;
    // Sum the signed products before asking for a field coordinate,
    // valuation, or residue. Individual products need not fit i128.
    let real = x
        .multiply(c.c)?
        .multiply(b.c)?
        .add(u.multiply(a.c)?.multiply(e.c)?.negate())?;
    let imaginary = y
        .multiply(c.c)?
        .multiply(b.c)?
        .add(v.multiply(a.c)?.multiply(e.c)?.negate())?;
    let denominator = Integer::product(a.c, e.c).multiply(c.c)?.multiply(b.c)?;
    Ok((real, imaginary, denominator))
}

fn ratio_small(d: i128, x: Element, y: Element) -> Option<(i128, i128, i128)> {
    let n =
        y.a.checked_mul(y.a)?
            .checked_sub(d.checked_mul(y.b.checked_mul(y.b)?)?)?;
    let a =
        x.a.checked_mul(y.a)?
            .checked_sub(d.checked_mul(x.b.checked_mul(y.b)?)?)?
            .checked_mul(y.c)?;
    let b =
        x.b.checked_mul(y.a)?
            .checked_sub(x.a.checked_mul(y.b)?)?
            .checked_mul(y.c)?;
    let c = x.c.checked_mul(n)?;
    if c < 0 {
        Some((a.checked_neg()?, b.checked_neg()?, c.checked_neg()?))
    } else if c == 0 {
        None
    } else {
        Some((a, b, c))
    }
}
