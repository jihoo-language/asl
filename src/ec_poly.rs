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
        let b = self.ring.negate(model, b)?;
        self.ring.add(model, a, &b)
    }

    fn multiply(&self, model: &mut Model, a: &Scalar, b: &Scalar) -> Result<Scalar> {
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
            return Err(Error::plain("isogeny polynomial exceeds the storage limit"));
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
    mut numerator: Polynomial,
    denominator: &[Scalar],
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
    while numerator.len() >= denominator.len() {
        let top = numerator.pop().expect("nonempty dividend");
        let width = denominator.len() - 1;
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
        while numerator.last().is_some_and(|c| zero(c)) {
            numerator.pop();
            model.spend(1)?;
        }
    }
    Ok(numerator)
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

pub(crate) fn coprime(
    model: &mut Model,
    field: Field<'_>,
    mut a: Polynomial,
    mut b: Polynomial,
) -> Result<bool> {
    while !b.is_empty() {
        let c = remainder(model, field, a, &b)?;
        (a, b) = (b, c);
    }
    Ok(a.len() == 1)
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
