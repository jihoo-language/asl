use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;
use crate::rv_ore::gcd;
use crate::wide_integer::Natural;

pub(crate) mod gauge;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Element {
    pub(crate) a: i128,
    pub(crate) b: i128,
    pub(crate) c: i128,
}

impl Element {
    pub(crate) const ZERO: Self = Self { a: 0, b: 0, c: 1 };
    pub(crate) const ONE: Self = Self { a: 1, b: 0, c: 1 };

    pub(crate) fn new(a: i128, b: i128, c: i128) -> Result<Self> {
        if c <= 0 {
            return Err(Error::plain("quadratic field denominator must be positive"));
        }
        let content = gcd(gcd(a.unsigned_abs(), b.unsigned_abs()), c as u128) as i128;
        Ok(Self {
            a: a / content,
            b: b / content,
            c: c / content,
        })
    }

    pub(crate) fn is_zero(self) -> bool {
        self.a == 0 && self.b == 0
    }

    pub(crate) fn rational_inverse(self, model: &mut Model) -> Result<Self> {
        debug_assert_eq!(self.b, 0);
        if self.a == 0 {
            return Err(Error::plain("zero has no inverse"));
        }
        model.spend(1)?;
        if self.a > 0 {
            Self::new(self.c, 0, self.a)
        } else {
            Self::new(-self.c, 0, negate(self.a)?)
        }
    }

    pub(crate) fn rational_valuation(self, model: &mut Model, prime: Prime) -> Result<Option<i64>> {
        debug_assert_eq!(self.b, 0);
        if self.a == 0 {
            return Ok(None);
        }
        let p = u128::from(prime.get());
        Ok(Some(
            valuation(model, self.a.unsigned_abs(), p)? - valuation(model, self.c as u128, p)?,
        ))
    }

    pub(crate) fn rational_angle(self, model: &mut Model, prime: Prime) -> Result<u32> {
        debug_assert!(self.b == 0 && self.a != 0);
        let p = i128::from(prime.get());
        let mut num = self.a;
        let mut den = self.c;
        while num % p == 0 {
            model.spend(1)?;
            num /= p;
        }
        while den % p == 0 {
            model.spend(1)?;
            den /= p;
        }
        let inverse = invert(model, prime, den.rem_euclid(p) as u32)?;
        Ok(prime.mul(num.rem_euclid(p) as u32, inverse))
    }

    pub(crate) fn is_prime_multiple_of(self, other: Self, prime: Prime) -> bool {
        let p = i128::from(prime.get());
        let common = gcd(p as u128, other.c as u128) as i128;
        self.c == other.c / common
            && other.a.checked_mul(p / common) == Some(self.a)
            && other.b.checked_mul(p / common) == Some(self.b)
    }
}

#[derive(Clone, Copy)]
pub(crate) struct Field {
    prime: Prime,
    d: i128,
    root: u32,
}

impl Field {
    pub(crate) fn new(model: &mut Model, p: usize, d: i128, root: i128) -> Result<Self> {
        if p < 3 || p.is_multiple_of(2) {
            return Err(Error::plain(
                "quadratic p-adic embedding requires an odd prime",
            ));
        }
        let prime = Prime::new(model, p)?;
        if d >= 0 && (d as u128).isqrt().pow(2) == d as u128 {
            return Err(Error::plain("quadratic field radicand must be nonsquare"));
        }
        let modulus = i128::from(prime.get());
        let r = root.rem_euclid(modulus) as u32;
        if r == 0 || prime.mul(r, r) != d.rem_euclid(modulus) as u32 {
            return Err(Error::plain(
                "quadratic embedding requires a nonzero root r with r^2 = d mod p",
            ));
        }
        Ok(Self { prime, d, root: r })
    }

    pub(crate) fn prime(self) -> Prime {
        self.prime
    }

    pub(crate) fn radicand(self) -> i128 {
        self.d
    }

    pub(crate) fn add(self, model: &mut Model, x: Element, y: Element) -> Result<Element> {
        affine_combination(model, x, y, false)
    }

    pub(crate) fn negate(self, model: &mut Model, x: Element) -> Result<Element> {
        model.spend(1)?;
        Element::new(negate(x.a)?, negate(x.b)?, x.c)
    }

    pub(crate) fn subtract(self, model: &mut Model, x: Element, y: Element) -> Result<Element> {
        affine_combination(model, x, y, true)
    }

    pub(crate) fn conjugate(self, model: &mut Model, x: Element) -> Result<Element> {
        model.spend(1)?;
        Element::new(x.a, negate(x.b)?, x.c)
    }

    pub(crate) fn multiply(self, model: &mut Model, x: Element, y: Element) -> Result<Element> {
        product(model, self.d, x, y)
    }

    pub(crate) fn inverse(self, model: &mut Model, x: Element) -> Result<Element> {
        if x.is_zero() {
            return Err(Error::plain("zero has no inverse in the quadratic field"));
        }
        if x.b == 0 {
            return x.rational_inverse(model);
        }
        model.spend(1)?;
        let norm = self.numerator_norm(model, x.a, x.b)?;
        if let Some(n) = norm.signed() {
            let small = || -> Option<_> {
                let cross = gcd(x.c as u128, n.unsigned_abs()) as i128;
                let mut unit = x.c / cross;
                let mut denominator = n / cross;
                if denominator < 0 {
                    denominator = denominator.checked_neg()?;
                    unit = -unit;
                }
                Some((x.a.checked_mul(unit)?, x.b.checked_mul(-unit)?, denominator))
            };
            if let Some((a, b, c)) = small() {
                return Element::new(a, b, c);
            }
        }
        model.spend(2 * 16)?;
        normalize_wide(
            model,
            QuadWide::product(x.a, x.c),
            QuadWide::product(x.b, x.c).negate(),
            norm,
        )
    }

    pub(crate) fn power(
        self,
        model: &mut Model,
        mut x: Element,
        mut exponent: usize,
    ) -> Result<Element> {
        let mut value = Element::ONE;
        while exponent != 0 {
            if !exponent.is_multiple_of(2) {
                value = self.multiply(model, value, x)?;
            }
            exponent /= 2;
            if exponent != 0 {
                x = self.multiply(model, x, x)?;
            }
        }
        Ok(value)
    }

    fn numerator_norm(self, model: &mut Model, a: i128, b: i128) -> Result<QuadWide> {
        let small = || -> Option<i128> {
            a.checked_mul(a)?
                .checked_sub(self.d.checked_mul(b.checked_mul(b)?)?)
        };
        if let Some(norm) = small() {
            return Ok(QuadWide::from_i128(norm));
        }
        model.spend(2 * 16 + 4 * 8 + 2 * QUADRATIC_LIMBS)?;
        QuadWide::product(a, a).add(QuadWide::product(b, b).multiply(self.d)?.negate())
    }

    pub(crate) fn valuation(self, model: &mut Model, x: Element) -> Result<Option<i64>> {
        if x.is_zero() {
            return Ok(None);
        }
        let p = u128::from(self.prime.get());
        let (a, b, common) = strip_numerator(model, x, p as i128)?;
        let (chosen, _) = self.local_numerator(a, b);
        let norm_order = if chosen == 0 {
            // One of these residues is a unit after common powers were
            // stripped, so the norm gives the full chosen branch order.
            self.numerator_norm(model, a, b)?
                .prime_part(model, self.prime)?
                .0
        } else {
            0
        };
        Ok(Some(
            common + norm_order - valuation(model, x.c as u128, p)?,
        ))
    }

    // The caller chooses a limb count that also accommodates the norm
    // of a canceled numerator. No i128 field coefficient is required.
    pub(crate) fn wide_local<const N: usize>(
        self,
        model: &mut Model,
        mut a: Wide<N>,
        mut b: Wide<N>,
    ) -> Result<Option<(i64, u32)>> {
        if a.is_zero() && b.is_zero() {
            return Ok(None);
        }
        let p = self.prime;
        let mut order = 0;
        let (residue_a, residue_b) = loop {
            model.spend(2 * N)?;
            let (qa, residue_a) = a.magnitude.divide_u32(p.get());
            let (qb, residue_b) = b.magnitude.divide_u32(p.get());
            if residue_a != 0 || residue_b != 0 {
                break (
                    if a.negative {
                        p.negate(residue_a)
                    } else {
                        residue_a
                    },
                    if b.negative {
                        p.negate(residue_b)
                    } else {
                        residue_b
                    },
                );
            }
            a.magnitude = qa;
            b.magnitude = qb;
            order += 1;
        };
        let b_root = p.mul(residue_b, self.root);
        let chosen = p.add(residue_a, b_root);
        if chosen != 0 {
            return Ok(Some((order, chosen)));
        }
        let an = a.magnitude.bits().div_ceil(32);
        let bn = b.magnitude.bits().div_ceil(32);
        model.spend(an * an + bn * bn + 8 * bn + 3 * N)?;
        let norm = a
            .multiply_wide(a)?
            .add(b.multiply_wide(b)?.multiply(self.d)?.negate())?;
        let (extra, unit) = norm.prime_part(model, p)?;
        let inverse = invert(model, p, p.add(residue_a, p.negate(b_root)))?;
        Ok(Some((order + extra, p.mul(unit, inverse))))
    }

    fn local_numerator(self, a: i128, b: i128) -> (u32, u32) {
        let p = i128::from(self.prime.get());
        let a = a.rem_euclid(p) as u32;
        let b_root = self.prime.mul(b.rem_euclid(p) as u32, self.root);
        (
            self.prime.add(a, b_root),
            self.prime.add(a, self.prime.negate(b_root)),
        )
    }

    pub(crate) fn residue(self, model: &mut Model, x: Element, precision: usize) -> Result<u128> {
        if precision == 0 {
            return Err(Error::plain("quadratic residue precision must be positive"));
        }
        let p = u128::from(self.prime.get());
        let modulus = prime_power(p, precision)?;
        if x.is_zero() {
            return Ok(0);
        }
        let (a, b, common) = strip_numerator(model, x, p as i128)?;
        let denominator_order = valuation(model, x.c as u128, p)?;
        let primitive = Element { a, b, c: 1 };
        let (chosen, _) = self.local_numerator(a, b);
        let (norm_order, norm_unit) = if chosen == 0 {
            let (order, unit) = self
                .numerator_norm(model, a, b)?
                .peel_prime(model, self.prime)?;
            (order, Some(unit))
        } else {
            (0, None)
        };
        let value_order =
            usize::try_from(common + norm_order - denominator_order).map_err(|_| {
                Error::plain("quadratic residue requires an integral element at the chosen prime")
            })?;
        if value_order >= precision {
            return Ok(0);
        }
        let unit_precision = precision - value_order;
        let unit_modulus = prime_power(p, unit_precision)?;
        let denominator = x.c as u128 / prime_power(p, denominator_order as usize)?;
        let (unit_numerator, unit_denominator) = if let Some(norm) = norm_unit {
            // The primitive numerator is N divided by its unit conjugate.
            // After removing the norm's p-power, only the output unit
            // precision is needed in the conjugate and denominator.
            let conjugate =
                self.numerator_residue(model, primitive, unit_precision, unit_modulus, true)?;
            (
                norm.modulus(model, unit_modulus)?,
                mul_mod(model, denominator % unit_modulus, conjugate, unit_modulus)?,
            )
        } else {
            (
                self.numerator_residue(model, primitive, unit_precision, unit_modulus, false)?,
                denominator % unit_modulus,
            )
        };
        let inverse = self.lift_inverse(model, unit_denominator, unit_precision)?;
        Ok((modulus / unit_modulus) * mul_mod(model, unit_numerator, inverse, unit_modulus)?)
    }

    pub(crate) fn angular_component(
        self,
        model: &mut Model,
        x: Element,
        order: i64,
    ) -> Result<u32> {
        if x.b == 0 {
            return x.rational_angle(model, self.prime);
        }
        let p = u128::from(self.prime.get());
        let (a, b, common) = strip_numerator(model, x, p as i128)?;
        let (chosen, conjugate) = self.local_numerator(a, b);
        let denominator_order = valuation(model, x.c as u128, p)?;
        let (norm_order, coefficient) = if chosen == 0 {
            let norm = self.numerator_norm(model, a, b)?;
            let (norm_order, unit) = norm.prime_part(model, self.prime)?;
            let inverse_conjugate = invert(model, self.prime, conjugate)?;
            (norm_order, self.prime.mul(unit, inverse_conjugate))
        } else {
            (0, chosen)
        };
        if common + norm_order - denominator_order != order {
            return Err(Error::plain(
                "quadratic angular component needs the valuation of its input",
            ));
        }
        let denominator_unit = x.c as u128 / prime_power(p, denominator_order as usize)?;
        let inverse = invert(model, self.prime, (denominator_unit % p) as u32)?;
        Ok(self.prime.mul(coefficient, inverse))
    }

    fn numerator_residue(
        self,
        model: &mut Model,
        x: Element,
        precision: usize,
        modulus: u128,
        conjugate: bool,
    ) -> Result<u128> {
        let a = reduce(x.a, modulus);
        if x.b == 0 {
            return Ok(a);
        }
        let b = reduce(x.b, modulus);
        let root = self.lift_root(model, precision)?;
        let term = mul_mod(model, b, root, modulus)?;
        Ok(if conjugate {
            subtract_mod(a, term, modulus)
        } else {
            (a + term) % modulus
        })
    }

    fn lift_root(self, model: &mut Model, precision: usize) -> Result<u128> {
        let p = u128::from(self.prime.get());
        let twice_root = self.prime.add(self.root, self.root);
        let slope_inverse = u128::from(invert(model, self.prime, twice_root)?);
        let mut root = u128::from(self.root);
        let mut modulus = p;
        for _ in 1..precision {
            model.spend(1)?;
            let next = modulus * p;
            let square = mul_mod(model, root, root, next)?;
            let error = subtract_mod(reduce(self.d, next), square, next);
            let digit = (error / modulus * slope_inverse) % p;
            root += digit * modulus;
            modulus = next;
        }
        Ok(root)
    }

    fn lift_inverse(self, model: &mut Model, unit: u128, precision: usize) -> Result<u128> {
        let p = u128::from(self.prime.get());
        let initial = u128::from(invert(model, self.prime, (unit % p) as u32)?);
        let mut inverse = initial;
        let mut modulus = p;
        for _ in 1..precision {
            model.spend(1)?;
            let next = modulus * p;
            let error = subtract_mod(1, mul_mod(model, unit % next, inverse, next)?, next);
            let digit = (error / modulus * initial) % p;
            inverse += digit * modulus;
            modulus = next;
        }
        Ok(inverse)
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Negate,
    Conjugate,
    Inverse,
    Power(usize),
}

impl Unary {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Negate => "neg",
            Self::Conjugate => "rvQpQuadConj",
            Self::Inverse => "rvQpQuadInv",
            Self::Power(_) => "rvQpQuadPow",
        }
    }

    pub(crate) fn apply(self, model: &mut Model, field: Field, input: Element) -> Result<Element> {
        match self {
            Self::Negate => field.negate(model, input),
            Self::Conjugate => field.conjugate(model, input),
            Self::Inverse => field.inverse(model, input),
            Self::Power(k) => field.power(model, input, k),
        }
    }
}

pub(crate) fn persist(model: &mut Model, value: Element) -> Result<Value> {
    model.allocate(6)?;
    Ok(Value::RVQpQuad(value))
}

fn product(model: &mut Model, d: i128, x: Element, y: Element) -> Result<Element> {
    model.spend(1)?;
    let g = gcd(gcd(x.a.unsigned_abs(), x.b.unsigned_abs()), y.c as u128) as i128;
    let h = gcd(gcd(y.a.unsigned_abs(), y.b.unsigned_abs()), x.c as u128) as i128;
    let (a, b, e, f) = (x.a / g, x.b / g, y.a / h, y.b / h);
    let small = || -> Option<_> {
        let real = a
            .checked_mul(e)?
            .checked_add(d.checked_mul(b.checked_mul(f)?)?)?;
        let imaginary = a.checked_mul(f)?.checked_add(b.checked_mul(e)?)?;
        Some((real, imaginary, (x.c / h).checked_mul(y.c / g)?))
    };
    if let Some((a, b, c)) = small() {
        return Element::new(a, b, c);
    }
    model.spend(5 * 16 + 4 * 8 + 4 * QUADRATIC_LIMBS)?;
    normalize_wide(
        model,
        QuadWide::product(a, e).add(QuadWide::product(b, f).multiply(d)?)?,
        QuadWide::product(a, f).add(QuadWide::product(b, e))?,
        QuadWide::product(x.c / h, y.c / g),
    )
}

fn affine_combination(model: &mut Model, x: Element, y: Element, minus: bool) -> Result<Element> {
    model.spend(1)?;
    let common = gcd(x.c as u128, y.c as u128) as i128;
    let (sx, sy) = (y.c / common, x.c / common);
    let small = || -> Option<_> {
        let combine = if minus {
            i128::checked_sub
        } else {
            i128::checked_add
        };
        Some((
            combine(x.a.checked_mul(sx)?, y.a.checked_mul(sy)?)?,
            combine(x.b.checked_mul(sx)?, y.b.checked_mul(sy)?)?,
            x.c.checked_mul(sx)?,
        ))
    };
    if let Some((a, b, c)) = small() {
        return Element::new(a, b, c);
    }
    model.spend(5 * 16 + 4 * QUADRATIC_LIMBS)?;
    let combine = |a, b| {
        let right = QuadWide::product(b, sy);
        QuadWide::product(a, sx).add(if minus { right.negate() } else { right })
    };
    normalize_wide(
        model,
        combine(x.a, y.a)?,
        combine(x.b, y.b)?,
        QuadWide::product(x.c, sx),
    )
}

// Three i128 factors have magnitude at most 2^381. The signed sum
// defining a norm or real product coordinate uses at most 382 bits.
const QUADRATIC_LIMBS: usize = 12;
type QuadWide = Wide<QUADRATIC_LIMBS>;

#[derive(Clone, Copy)]
pub(crate) struct Wide<const N: usize> {
    pub(crate) magnitude: Natural<N>,
    pub(crate) negative: bool,
}

impl<const N: usize> Wide<N> {
    pub(crate) fn from_i128(value: i128) -> Self {
        Self {
            magnitude: Natural::from_small(value.unsigned_abs()),
            negative: value < 0,
        }
    }

    pub(crate) fn product(a: i128, b: i128) -> Self {
        Self {
            magnitude: Natural::product(a, b),
            negative: (a < 0) != (b < 0),
        }
    }

    pub(crate) fn resize<const M: usize>(self) -> Result<Wide<M>> {
        Ok(Wide {
            magnitude: self.magnitude.resize().ok_or_else(arithmetic)?,
            negative: self.negative,
        })
    }

    pub(crate) fn multiply(self, coefficient: i128) -> Result<Self> {
        let magnitude = self
            .magnitude
            .multiply_u128(coefficient.unsigned_abs())
            .ok_or_else(arithmetic)?;
        Ok(Self {
            magnitude,
            negative: self.negative != (coefficient < 0),
        })
    }

    pub(crate) fn is_zero(self) -> bool {
        self.magnitude.is_zero()
    }

    pub(crate) fn multiply_wide(self, rhs: Self) -> Result<Self> {
        let magnitude = self
            .magnitude
            .checked_product(rhs.magnitude)
            .ok_or_else(arithmetic)?;
        Ok(Self {
            magnitude,
            negative: self.negative != rhs.negative,
        })
    }

    pub(crate) fn divide_exact(self, model: &mut Model, rhs: Self) -> Result<Option<Self>> {
        if rhs.is_zero() {
            return Err(Error::plain("exact wide division by zero"));
        }
        Ok(self
            .magnitude
            .divide_exact(model, rhs.magnitude)?
            .map(|magnitude| Self {
                magnitude,
                negative: self.negative != rhs.negative,
            }))
    }

    pub(crate) fn negate(self) -> Self {
        Self {
            negative: !self.negative,
            ..self
        }
    }

    pub(crate) fn add(mut self, mut other: Self) -> Result<Self> {
        if self.negative == other.negative {
            self.magnitude
                .add_multiple(other.magnitude, 1)
                .ok_or_else(arithmetic)?;
        } else {
            if self.magnitude.cmp(other.magnitude).is_lt() {
                (self, other) = (other, self);
            }
            self.magnitude = self.magnitude.difference(other.magnitude);
        }
        Ok(self)
    }

    pub(crate) fn signed(self) -> Option<i128> {
        self.magnitude.signed(self.negative)
    }

    pub(crate) fn prime_part(self, model: &mut Model, p: Prime) -> Result<(i64, u32)> {
        let (order, unit) = self.peel_prime(model, p)?;
        let residue = unit.magnitude.residue(p.get());
        Ok((
            order,
            if unit.negative {
                p.negate(residue)
            } else {
                residue
            },
        ))
    }

    fn peel_prime(mut self, model: &mut Model, p: Prime) -> Result<(i64, Self)> {
        if self.magnitude.is_zero() {
            return Err(Error::plain("nonzero quadratic numerator has zero norm"));
        }
        let mut v = 0;
        if let Some(mut small) = self.magnitude.small() {
            let p = u128::from(p.get());
            while small.is_multiple_of(p) {
                model.spend(1)?;
                small /= p;
                v += 1;
            }
            self.magnitude = Natural::from_small(small);
            return Ok((v, self));
        }
        loop {
            model.spend(N)?;
            let (quotient, residue) = self.magnitude.divide_u32(p.get());
            if residue != 0 {
                return Ok((v, self));
            }
            self.magnitude = quotient;
            v += 1;
        }
    }

    fn modulus(self, model: &mut Model, modulus: u128) -> Result<u128> {
        let (_, residue) = self.magnitude.divide(model, modulus)?;
        Ok(if self.negative && residue != 0 {
            modulus - residue
        } else {
            residue
        })
    }
}

fn normalize_wide<const N: usize>(
    model: &mut Model,
    a: Wide<N>,
    b: Wide<N>,
    c: Wide<N>,
) -> Result<Element> {
    if c.magnitude.is_zero() {
        return Err(Error::plain("quadratic field denominator must be nonzero"));
    }
    let common = a
        .magnitude
        .gcd(model, b.magnitude)?
        .gcd(model, c.magnitude)?;
    let common = Wide {
        magnitude: common,
        negative: false,
    };
    let mut coefficient = |x: Wide<N>, negative| {
        x.divide_exact(model, common)?
            .and_then(|value| value.magnitude.signed(negative))
            .ok_or_else(arithmetic)
    };
    let real = coefficient(a, a.negative != c.negative)?;
    let imaginary = coefficient(b, b.negative != c.negative)?;
    let denominator = coefficient(c, false)?;
    Ok(Element {
        a: real,
        b: imaginary,
        c: denominator,
    })
}

fn strip_numerator(model: &mut Model, x: Element, p: i128) -> Result<(i128, i128, i64)> {
    debug_assert!(!x.is_zero());
    let (mut a, mut b) = (x.a, x.b);
    let mut order = 0;
    while a % p == 0 && b % p == 0 {
        model.spend(1)?;
        a /= p;
        b /= p;
        order += 1;
    }
    Ok((a, b, order))
}

fn valuation(model: &mut Model, mut nonzero: u128, p: u128) -> Result<i64> {
    debug_assert_ne!(nonzero, 0);
    let mut value = 0;
    while nonzero.is_multiple_of(p) {
        model.spend(1)?;
        nonzero /= p;
        value += 1;
    }
    Ok(value)
}

fn prime_power(prime: u128, exponent: usize) -> Result<u128> {
    u32::try_from(exponent)
        .ok()
        .and_then(|exp| prime.checked_pow(exp))
        .filter(|&modulus| modulus <= i128::MAX as u128)
        .ok_or_else(modulus_error)
}

fn reduce(value: i128, modulus: u128) -> u128 {
    value.rem_euclid(modulus as i128) as u128
}

fn subtract_mod(left: u128, right: u128, modulus: u128) -> u128 {
    if left >= right {
        left - right
    } else {
        modulus - (right - left)
    }
}

pub(crate) fn mul_mod(model: &mut Model, mut a: u128, mut b: u128, modulus: u128) -> Result<u128> {
    model.spend(1)?;
    if let Some(value) = a.checked_mul(b) {
        return Ok(value % modulus);
    }
    // Both residues are smaller than an i128-sized modulus, so addition
    // of two residues is exact in u128 even when their product is wider.
    let mut result = 0;
    while b != 0 {
        model.spend(1)?;
        if !b.is_multiple_of(2) {
            result = (result + a) % modulus;
        }
        a = (a + a) % modulus;
        b /= 2;
    }
    Ok(result)
}

fn negate(a: i128) -> Result<i128> {
    a.checked_neg().ok_or_else(arithmetic)
}

const COEFFICIENT_RANGE: &str = "quadratic exact integer arithmetic is out of range";

pub(crate) fn is_coefficient_overflow(error: &Error) -> bool {
    error.message == COEFFICIENT_RANGE
}

fn arithmetic() -> Error {
    Error::plain(COEFFICIENT_RANGE)
}

fn modulus_error() -> Error {
    Error::plain("quadratic residue modulus at the required precision is out of range")
}
