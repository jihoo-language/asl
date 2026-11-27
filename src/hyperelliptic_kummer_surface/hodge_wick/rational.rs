use std::cmp::Ordering;
use std::fmt::Write;

use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::wide_integer::Natural;

pub(super) const LIMBS: usize = 32;
pub(super) type Integer = Wide<LIMBS>;
type Magnitude = Natural<LIMBS>;

pub(super) fn arithmetic() -> Error {
    Error::plain("Hodge Wick arithmetic exceeds the checked 1024-bit exact range")
}

#[derive(Clone, Copy)]
pub(crate) struct Rational {
    numerator: Integer,
    denominator: Magnitude,
}

impl PartialEq for Rational {
    fn eq(&self, other: &Self) -> bool {
        self.numerator.negative == other.numerator.negative
            && self.numerator.magnitude == other.numerator.magnitude
            && self.denominator == other.denominator
    }
}

impl Eq for Rational {}

impl Rational {
    pub(super) fn zero() -> Self {
        Self::integer(0)
    }

    pub(super) fn integer(value: i128) -> Self {
        Self::wide(Integer::from_i128(value))
    }

    pub(super) fn wide(mut numerator: Integer) -> Self {
        if numerator.is_zero() {
            numerator.negative = false;
        }
        Self {
            numerator,
            denominator: Magnitude::from_small(1),
        }
    }

    pub(super) fn input(model: &mut Model, value: (i128, i128)) -> Result<Self> {
        if value.1 <= 0 {
            return Err(Error::plain(
                "Hodge Wick rational denominator must be positive",
            ));
        }
        Self::normalize(
            model,
            Integer::from_i128(value.0),
            Magnitude::from_small(value.1 as u128),
        )
    }

    pub(super) fn is_zero(self) -> bool {
        self.numerator.is_zero()
    }

    pub(super) fn negative(self) -> bool {
        self.numerator.negative
    }

    pub(super) fn negate(self) -> Self {
        if self.is_zero() {
            self
        } else {
            Self {
                numerator: self.numerator.negate(),
                ..self
            }
        }
    }

    pub(super) fn inverse(self) -> Result<Self> {
        if self.is_zero() {
            return Err(Error::plain(
                "division by zero in Hodge Wick rational arithmetic",
            ));
        }
        Ok(Self {
            numerator: Integer {
                magnitude: self.denominator,
                negative: self.numerator.negative,
            },
            denominator: self.numerator.magnitude,
        })
    }

    pub(super) fn add(self, model: &mut Model, rhs: Self) -> Result<Self> {
        if self.is_zero() {
            return Ok(rhs);
        }
        if rhs.is_zero() {
            return Ok(self);
        }
        model.spend(LIMBS)?;
        if self.denominator == rhs.denominator {
            return Self::normalize(
                model,
                self.numerator
                    .add(rhs.numerator)
                    .map_err(|_| arithmetic())?,
                self.denominator,
            );
        }
        let common = gcd(model, self.denominator, rhs.denominator)?;
        let left = divide(model, self.denominator, common)?;
        let right = divide(model, rhs.denominator, common)?;
        let sum = product(self.numerator, right)?
            .add(product(rhs.numerator, left)?)
            .map_err(|_| arithmetic())?;
        if sum.is_zero() {
            return Ok(Self::zero());
        }
        let cancel = gcd(model, sum.magnitude, common)?;
        let denominator = left
            .checked_product(divide(model, rhs.denominator, cancel)?)
            .ok_or_else(arithmetic)?;
        Ok(Self {
            numerator: Integer {
                magnitude: divide(model, sum.magnitude, cancel)?,
                negative: sum.negative,
            },
            denominator,
        })
    }

    pub(super) fn multiply(self, model: &mut Model, rhs: Self) -> Result<Self> {
        if self.is_zero() || rhs.is_zero() {
            return Ok(Self::zero());
        }
        model.spend(LIMBS)?;
        let a = gcd(model, self.numerator.magnitude, rhs.denominator)?;
        let b = gcd(model, rhs.numerator.magnitude, self.denominator)?;
        let numerator = divide(model, self.numerator.magnitude, a)?
            .checked_product(divide(model, rhs.numerator.magnitude, b)?)
            .ok_or_else(arithmetic)?;
        let denominator = divide(model, self.denominator, b)?
            .checked_product(divide(model, rhs.denominator, a)?)
            .ok_or_else(arithmetic)?;
        Ok(Self {
            numerator: Integer {
                magnitude: numerator,
                negative: self.negative() != rhs.negative(),
            },
            denominator,
        })
    }

    fn normalize(model: &mut Model, numerator: Integer, denominator: Magnitude) -> Result<Self> {
        if numerator.is_zero() {
            return Ok(Self::zero());
        }
        let common = gcd(model, numerator.magnitude, denominator)?;
        Ok(Self {
            numerator: Integer {
                magnitude: divide(model, numerator.magnitude, common)?,
                negative: numerator.negative,
            },
            denominator: divide(model, denominator, common)?,
        })
    }

    pub(super) fn compare(self, rhs: Self) -> Ordering {
        self.numerator
            .negative
            .cmp(&rhs.numerator.negative)
            .then_with(|| self.numerator.magnitude.cmp(rhs.numerator.magnitude))
            .then_with(|| self.denominator.cmp(rhs.denominator))
    }

    pub(super) fn small_parts(self) -> Option<(i128, u128)> {
        let numerator = self.numerator.magnitude.small()?;
        let numerator = if self.negative() && numerator == 1 << 127 {
            i128::MIN
        } else {
            let numerator = i128::try_from(numerator).ok()?;
            if self.negative() {
                -numerator
            } else {
                numerator
            }
        };
        Some((numerator, self.denominator.small()?))
    }

    pub(super) fn text(self, model: &mut Model) -> Result<String> {
        let mut text = String::new();
        if self.negative() {
            text.push('-');
        }
        text.push_str(&decimal(model, self.numerator.magnitude)?);
        if self.denominator.small() != Some(1) {
            text.push('/');
            text.push_str(&decimal(model, self.denominator)?);
        }
        Ok(text)
    }
}

fn gcd(model: &mut Model, a: Magnitude, b: Magnitude) -> Result<Magnitude> {
    if a.small() == Some(1) || b.small() == Some(1) {
        return Ok(Magnitude::from_small(1));
    }
    a.gcd(model, b)
}

fn divide(model: &mut Model, value: Magnitude, divisor: Magnitude) -> Result<Magnitude> {
    if divisor.small() == Some(1) {
        return Ok(value);
    }
    value.divide_exact(model, divisor)?.ok_or_else(arithmetic)
}

fn product(value: Integer, factor: Magnitude) -> Result<Integer> {
    Ok(Integer {
        magnitude: value
            .magnitude
            .checked_product(factor)
            .ok_or_else(arithmetic)?,
        negative: value.negative,
    })
}

fn decimal(model: &mut Model, mut value: Magnitude) -> Result<String> {
    if let Some(small) = value.small() {
        return Ok(small.to_string());
    }
    let mut digits = Vec::new();
    while !value.is_zero() {
        model.spend(LIMBS)?;
        let (quotient, remainder) = value.divide_u32(1_000_000_000);
        digits.push(remainder);
        value = quotient;
    }
    let mut text = digits.pop().unwrap().to_string();
    while let Some(digit) = digits.pop() {
        write!(text, "{digit:09}").unwrap();
    }
    Ok(text)
}

pub(super) struct Q;

impl CoefficientRing for Q {
    type Element = Rational;

    fn cells(&self) -> usize {
        2 * LIMBS + 1
    }
    fn zero(&self) -> Rational {
        Rational::zero()
    }
    fn one(&self) -> Rational {
        Rational::integer(1)
    }
    fn is_zero(&self, a: &Rational) -> bool {
        a.is_zero()
    }
    fn is_unit(&self, a: &Rational) -> bool {
        !a.is_zero()
    }
    fn add(&self, model: &mut Model, a: &Rational, b: &Rational) -> Result<Rational> {
        a.add(model, *b)
    }
    fn subtract(&self, model: &mut Model, a: &Rational, b: &Rational) -> Result<Rational> {
        a.add(model, b.negate())
    }
    fn multiply(&self, model: &mut Model, a: &Rational, b: &Rational) -> Result<Rational> {
        a.multiply(model, *b)
    }
    fn inverse_unit(&self, _: &mut Model, a: &Rational) -> Result<Rational> {
        a.inverse()
    }
}
