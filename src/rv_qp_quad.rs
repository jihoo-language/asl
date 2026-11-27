use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;
use crate::rv_ore::gcd;

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
        model.spend(1)?;
        let g = gcd(gcd(x.a.unsigned_abs(), x.b.unsigned_abs()), y.c as u128) as i128;
        let h = gcd(gcd(y.a.unsigned_abs(), y.b.unsigned_abs()), x.c as u128) as i128;
        let (a, b, e, f) = (x.a / g, x.b / g, y.a / h, y.b / h);
        let real = add(multiply(a, e)?, multiply(self.d, multiply(b, f)?)?)?;
        let imaginary = add(multiply(a, f)?, multiply(b, e)?)?;
        Element::new(real, imaginary, multiply(x.c / h, y.c / g)?)
    }

    pub(crate) fn inverse(self, model: &mut Model, x: Element) -> Result<Element> {
        if x.is_zero() {
            return Err(Error::plain("zero has no inverse in the quadratic field"));
        }
        model.spend(1)?;
        if x.b == 0 {
            return if x.a > 0 {
                Element::new(x.c, 0, x.a)
            } else {
                Element::new(-x.c, 0, negate(x.a)?)
            };
        }
        let norm = self.numerator_norm(x.a, x.b)?;
        let cross = gcd(x.c as u128, norm.unsigned_abs()) as i128;
        let mut unit = x.c / cross;
        let mut denominator = norm / cross;
        if denominator < 0 {
            denominator = negate(denominator)?;
            unit = -unit;
        }
        Element::new(multiply(x.a, unit)?, multiply(x.b, -unit)?, denominator)
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

    fn numerator_norm(self, a: i128, b: i128) -> Result<i128> {
        subtract(multiply(a, a)?, multiply(self.d, multiply(b, b)?)?)
    }

    pub(crate) fn valuation(self, model: &mut Model, x: Element) -> Result<Option<i64>> {
        if x.is_zero() {
            return Ok(None);
        }
        let p = i128::from(self.prime.get());
        let (mut a, mut b) = (x.a, x.b);
        let mut value = 0;
        while a % p == 0 && b % p == 0 {
            model.spend(1)?;
            a /= p;
            b /= p;
            value += 1;
        }
        let reduced_a = a.rem_euclid(p) as u32;
        let reduced_b = b.rem_euclid(p) as u32;
        if self
            .prime
            .add(reduced_a, self.prime.mul(reduced_b, self.root))
            == 0
        {
            // After stripping the common p-powers, cancellation at the
            // chosen prime means the conjugate is a unit. The norm gives
            // exactly the eventual cancellation depth in the chosen branch.
            let norm = self.numerator_norm(a, b)?;
            value += valuation(model, norm.unsigned_abs(), p as u128)?;
        }
        value -= valuation(model, x.c as u128, p as u128)?;
        Ok(Some(value))
    }

    pub(crate) fn residue(self, model: &mut Model, x: Element, precision: usize) -> Result<u128> {
        if precision == 0 {
            return Err(Error::plain("quadratic residue precision must be positive"));
        }
        let p = u128::from(self.prime.get());
        let denom_order = valuation(model, x.c as u128, p)? as usize;
        let root_precision = precision
            .checked_add(denom_order)
            .ok_or_else(modulus_error)?;
        let modulus = prime_power(p, precision)?;
        let enlarged = prime_power(p, root_precision)?;
        let divisor = enlarged / modulus;
        let root = self.lift_root(model, root_precision)?;
        let a = reduce(x.a, enlarged);
        let b = reduce(x.b, enlarged);
        let reduced = (a + mul_mod(model, b, root, enlarged)?) % enlarged;
        if !reduced.is_multiple_of(divisor) {
            return Err(Error::plain(
                "quadratic residue requires an integral element at the chosen prime",
            ));
        }
        let denominator_unit = x.c as u128 / divisor;
        let inverse = self.lift_inverse(model, denominator_unit, precision)?;
        mul_mod(model, reduced / divisor, inverse, modulus)
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

fn affine_combination(model: &mut Model, x: Element, y: Element, minus: bool) -> Result<Element> {
    model.spend(1)?;
    let common = gcd(x.c as u128, y.c as u128) as i128;
    let (sx, sy) = (y.c / common, x.c / common);
    let sign = if minus { subtract } else { add };
    Element::new(
        sign(multiply(x.a, sx)?, multiply(y.a, sy)?)?,
        sign(multiply(x.b, sx)?, multiply(y.b, sy)?)?,
        multiply(x.c, sx)?,
    )
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

fn mul_mod(model: &mut Model, mut a: u128, mut b: u128, modulus: u128) -> Result<u128> {
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

fn add(a: i128, b: i128) -> Result<i128> {
    a.checked_add(b).ok_or_else(arithmetic)
}

fn subtract(a: i128, b: i128) -> Result<i128> {
    a.checked_sub(b).ok_or_else(arithmetic)
}

fn multiply(a: i128, b: i128) -> Result<i128> {
    a.checked_mul(b).ok_or_else(arithmetic)
}

fn negate(a: i128) -> Result<i128> {
    a.checked_neg().ok_or_else(arithmetic)
}

fn arithmetic() -> Error {
    Error::plain("quadratic exact integer arithmetic is out of range")
}

fn modulus_error() -> Error {
    Error::plain("quadratic residue modulus at the required precision is out of range")
}
