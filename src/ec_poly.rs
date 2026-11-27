use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing};
use crate::model::Model;
use crate::rv_unram::Ring;

pub(crate) type Scalar = Vec<u64>;
pub(crate) type Polynomial = Vec<Scalar>;

#[derive(Clone, Copy)]
pub(crate) struct Field<'a> {
    pub(crate) ring: &'a Ring,
}

impl CoefficientRing for Field<'_> {
    type Element = Scalar;

    fn cells(&self) -> usize {
        3 + self.ring.degree()
    }

    fn zero(&self) -> Scalar {
        self.ring.scalar(0)
    }

    fn one(&self) -> Scalar {
        self.ring.scalar(1)
    }

    fn is_zero(&self, value: &Scalar) -> bool {
        zero(value)
    }

    fn is_unit(&self, value: &Scalar) -> bool {
        !zero(value)
    }

    fn add(&self, model: &mut Model, a: &Scalar, b: &Scalar) -> Result<Scalar> {
        self.ring.add(model, a, b)
    }

    fn subtract(&self, model: &mut Model, a: &Scalar, b: &Scalar) -> Result<Scalar> {
        let modulus = self.ring.modulus();
        model.spend(self.ring.degree())?;
        Ok(a.iter()
            .zip(b)
            .map(|(&x, &y)| crate::rv_unram::sub_mod(x, y, modulus))
            .collect())
    }

    fn multiply(&self, model: &mut Model, a: &Scalar, b: &Scalar) -> Result<Scalar> {
        // A coefficient in the base ring needs only scalar multiplication,
        // including when the entire polynomial has been extended from it.
        model.spend(a.len().saturating_sub(1))?;
        if zero(&a[1..]) {
            return self.ring.scale(model, b, a[0]);
        }
        model.spend(b.len().saturating_sub(1))?;
        if zero(&b[1..]) {
            return self.ring.scale(model, a, b[0]);
        }
        self.ring.multiply(model, a, b)
    }

    fn inverse_unit(&self, model: &mut Model, a: &Scalar) -> Result<Scalar> {
        self.ring.inverse(model, a)
    }
}

impl Field<'_> {
    pub(crate) fn signed(self, c: i128) -> Scalar {
        self.ring
            .scalar(c.rem_euclid(i128::from(self.ring.prime())) as u64)
    }

    pub(crate) fn polynomial(self, constants: &[i128]) -> Polynomial {
        trim(constants.iter().map(|&c| self.signed(c)).collect())
    }

    pub(crate) fn check_polynomial(self, len: usize) -> Result<()> {
        if len.saturating_mul(self.cells()) > crate::complex::MAX_PIVOT_CELLS {
            return Err(Error::plain(
                "finite-field polynomial exceeds the storage limit",
            ));
        }
        Ok(())
    }
}

pub(crate) fn zero(value: &[u64]) -> bool {
    value.iter().all(|&c| c == 0)
}

pub(crate) fn trim(mut values: Polynomial) -> Polynomial {
    values.truncate(values.iter().rposition(|c| !zero(c)).map_or(0, |i| i + 1));
    values
}

pub(crate) fn scale(
    model: &mut Model,
    field: Field<'_>,
    a: &[Scalar],
    value: &Scalar,
) -> Result<Polynomial> {
    if zero(value) {
        return Ok(Vec::new());
    }
    model.spend(a.len())?;
    a.iter().map(|c| field.multiply(model, c, value)).collect()
}

pub(crate) fn sum(
    model: &mut Model,
    field: Field<'_>,
    a: &[Scalar],
    b: &[Scalar],
    negate: bool,
) -> Result<Polynomial> {
    let width = a.len().max(b.len());
    field.check_polynomial(width)?;
    model.spend(width)?;
    let zero = field.zero();
    let values = (0..width)
        .map(|i| {
            let x = a.get(i).unwrap_or(&zero);
            let y = b.get(i).unwrap_or(&zero);
            if negate {
                field.subtract(model, x, y)
            } else {
                field.add(model, x, y)
            }
        })
        .collect::<Result<_>>()?;
    Ok(trim(values))
}

pub(crate) fn multiply(
    model: &mut Model,
    field: Field<'_>,
    a: &[Scalar],
    b: &[Scalar],
) -> Result<Polynomial> {
    hw_algebra::multiply_polynomials(model, &field, a, b).map(trim)
}

pub(crate) fn derivative(model: &mut Model, field: Field<'_>, a: &[Scalar]) -> Result<Polynomial> {
    model.spend(a.len())?;
    let result = a
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| field.ring.scale(model, c, i as u64 % field.ring.prime()))
        .collect::<Result<_>>()?;
    Ok(trim(result))
}

pub(crate) fn evaluate(
    model: &mut Model,
    field: Field<'_>,
    a: &[Scalar],
    x: &[u64],
) -> Result<Scalar> {
    let mut value = field.zero();
    for c in a.iter().rev() {
        let product = field.ring.multiply(model, &value, x)?;
        value = field.add(model, &product, c)?;
    }
    Ok(value)
}

pub(crate) fn remainder(
    model: &mut Model,
    field: Field<'_>,
    numerator: Polynomial,
    denominator: &[Scalar],
) -> Result<Polynomial> {
    division(model, field, numerator, denominator, |_, _| {})
}

pub(crate) fn divide(
    model: &mut Model,
    field: Field<'_>,
    numerator: Polynomial,
    denominator: &[Scalar],
) -> Result<(Polynomial, Polynomial)> {
    let size = numerator
        .len()
        .saturating_add(1)
        .saturating_sub(denominator.len());
    field.check_polynomial(size)?;
    model.spend(size)?;
    let mut quotient = vec![field.zero(); size];
    let remainder = division(model, field, numerator, denominator, |order, c| {
        quotient[order] = c
    })?;
    Ok((trim(quotient), remainder))
}

fn division(
    model: &mut Model,
    field: Field<'_>,
    mut numerator: Polynomial,
    denominator: &[Scalar],
    mut quotient: impl FnMut(usize, Scalar),
) -> Result<Polynomial> {
    let leading = denominator
        .last()
        .ok_or_else(|| Error::plain("polynomial division by zero"))?;
    if numerator.len() < denominator.len() {
        return Ok(numerator);
    }
    let one = field.one();
    let inverse = if leading == &one {
        one
    } else {
        field.inverse_unit(model, leading)?
    };
    let width = denominator.len() - 1;
    while numerator.len() >= denominator.len() {
        let top = numerator.pop().expect("nonempty dividend");
        model.spend(1 + width)?;
        if zero(&top) {
            continue;
        }
        let order = numerator.len() - width;
        let factor = field.multiply(model, &top, &inverse)?;
        for (dst, c) in numerator[order..].iter_mut().zip(denominator) {
            if zero(c) {
                continue;
            }
            let product = field.multiply(model, &factor, c)?;
            *dst = field.subtract(model, dst, &product)?;
        }
        quotient(order, factor);
        while numerator.last().is_some_and(|c| zero(c)) {
            numerator.pop();
            model.spend(1)?;
        }
    }
    Ok(numerator)
}

pub(crate) fn exact_divide(
    model: &mut Model,
    field: Field<'_>,
    a: &[Scalar],
    b: &[Scalar],
) -> Result<Polynomial> {
    model.spend(a.len())?;
    let (quotient, remainder) = divide(model, field, a.to_vec(), b)?;
    if !remainder.is_empty() {
        return Err(Error::plain(
            "expected an exact polynomial division over the residue field",
        ));
    }
    Ok(quotient)
}

pub(crate) fn multiply_mod(
    model: &mut Model,
    field: Field<'_>,
    a: &[Scalar],
    b: &[Scalar],
    modulus: &[Scalar],
) -> Result<Polynomial> {
    let product = multiply(model, field, a, b)?;
    remainder(model, field, product, modulus)
}

pub(crate) fn power_mod(
    model: &mut Model,
    field: Field<'_>,
    base: &[Scalar],
    mut exponent: u64,
    modulus: &[Scalar],
) -> Result<Polynomial> {
    let mut result = remainder(model, field, vec![field.one()], modulus)?;
    let mut power = remainder(model, field, base.to_vec(), modulus)?;
    while exponent != 0 {
        model.spend(1)?;
        if exponent & 1 != 0 {
            result = multiply_mod(model, field, &result, &power, modulus)?;
        }
        exponent >>= 1;
        if exponent != 0 {
            power = multiply_mod(model, field, &power, &power, modulus)?;
        }
    }
    Ok(result)
}

pub(crate) fn gcd(
    model: &mut Model,
    field: Field<'_>,
    mut a: Polynomial,
    mut b: Polynomial,
) -> Result<Polynomial> {
    while !b.is_empty() {
        let c = remainder(model, field, a, &b)?;
        (a, b) = (b, c);
    }
    Ok(a)
}

pub(crate) fn resultant(
    model: &mut Model,
    field: Field<'_>,
    mut a: Polynomial,
    b: &[Scalar],
) -> Result<Scalar> {
    if a.is_empty() || b.is_empty() {
        return Ok(field.zero());
    }
    model.spend(b.len() * field.cells())?;
    let mut b = b.to_vec();
    let one = field.one();
    let mut result = field.one();
    loop {
        let m = a.len() - 1;
        let n = b.len() - 1;
        let leading = &b[n];
        if n == 0 {
            if leading == &one {
                return Ok(result);
            }
            let factor = field.ring.power(model, leading, m as u64)?;
            return field.multiply(model, &result, &factor);
        }
        let remainder = remainder(model, field, a, &b)?;
        if remainder.is_empty() {
            return Ok(field.zero());
        }
        // Res(a,b) = (−1)^(mn) lc(b)^(m−deg r) Res(b,r) for r = a mod b.
        if leading != &one {
            let factor = field
                .ring
                .power(model, leading, (m + 1 - remainder.len()) as u64)?;
            result = field.multiply(model, &result, &factor)?;
        }
        if m % 2 == 1 && n % 2 == 1 {
            result = field.ring.negate(model, &result)?;
        }
        (a, b) = (b, remainder);
    }
}

pub(crate) fn coprime(
    model: &mut Model,
    field: Field<'_>,
    a: Polynomial,
    b: Polynomial,
) -> Result<bool> {
    Ok(gcd(model, field, a, b)?.len() == 1)
}

// The returned gcd is monic and equals left·a+right·b.
pub(crate) fn bezout(
    model: &mut Model,
    field: Field<'_>,
    a: &[Scalar],
    b: &[Scalar],
) -> Result<(Polynomial, Polynomial, Polynomial)> {
    model.spend(a.len() + b.len())?;
    let (mut r0, mut r1) = (a.to_vec(), b.to_vec());
    let (mut a0, mut a1) = (vec![field.one()], Vec::new());
    let (mut b0, mut b1) = (Vec::new(), vec![field.one()]);
    while !r1.is_empty() {
        let (quotient, remainder) = divide(model, field, r0, &r1)?;
        let next_a = multiply(model, field, &quotient, &a1)?;
        let next_a = sum(model, field, &a0, &next_a, true)?;
        let next_b = multiply(model, field, &quotient, &b1)?;
        let next_b = sum(model, field, &b0, &next_b, true)?;
        (r0, r1) = (r1, remainder);
        (a0, a1) = (a1, next_a);
        (b0, b1) = (b1, next_b);
    }
    let leading = r0
        .last()
        .ok_or_else(|| Error::plain("Bézout gcd of two zero polynomials"))?;
    if *leading == field.one() {
        return Ok((r0, a0, b0));
    }
    let inverse = field.inverse_unit(model, leading)?;
    Ok((
        scale(model, field, &r0, &inverse)?,
        scale(model, field, &a0, &inverse)?,
        scale(model, field, &b0, &inverse)?,
    ))
}

pub(crate) fn inverse_mod(
    model: &mut Model,
    field: Field<'_>,
    a: &[Scalar],
    modulus: &[Scalar],
) -> Result<Polynomial> {
    model.spend(a.len() + modulus.len())?;
    let mut r0 = modulus.to_vec();
    let mut r1 = remainder(model, field, a.to_vec(), modulus)?;
    let mut v0 = Vec::new();
    let mut v1 = vec![field.one()];
    while !r1.is_empty() {
        let (quotient, remainder) = divide(model, field, r0, &r1)?;
        let change = multiply_mod(model, field, &quotient, &v1, modulus)?;
        let value = sum(model, field, &v0, &change, true)?;
        (r0, r1) = (r1, remainder);
        (v0, v1) = (v1, value);
    }
    if r0.len() != 1 {
        return Err(Error::plain(
            "rational denominator is not a unit in the finite kernel algebra",
        ));
    }
    let inv = field.inverse_unit(model, &r0[0])?;
    scale(model, field, &v0, &inv)
}

pub(crate) fn power_sums(
    model: &mut Model,
    field: Field<'_>,
    monic: &[Scalar],
) -> Result<[Scalar; 4]> {
    let n = monic.len() - 1;
    let mut sums = [
        field.signed(n as i128),
        field.zero(),
        field.zero(),
        field.zero(),
    ];
    for i in 1..=3 {
        let mut sum = field.zero();
        for j in 1..=n.min(i - 1) {
            let term = field.multiply(model, &monic[n - j], &sums[i - j])?;
            sum = field.add(model, &sum, &term)?;
        }
        if i <= n {
            let last = field
                .ring
                .scale(model, &monic[n - i], i as u64 % field.ring.prime())?;
            sum = field.add(model, &sum, &last)?;
        }
        sums[i] = field.ring.negate(model, &sum)?;
    }
    Ok(sums)
}
