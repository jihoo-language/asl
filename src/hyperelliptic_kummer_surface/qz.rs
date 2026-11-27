use std::cmp::Ordering;

use super::rational_sparse::{decimal, Integer, Magnitude, WORK_LIMBS};
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::wide_integer::Natural;

const LIMBS: usize = 8;

#[derive(Clone, Copy, Eq, PartialEq)]
pub(crate) struct Scalar {
    numerator: Natural<LIMBS>,
    denominator: Natural<LIMBS>,
}

impl Scalar {
    pub(crate) fn zero() -> Self {
        Self {
            numerator: Natural::ZERO,
            denominator: Natural::from_small(1),
        }
    }

    pub(crate) fn literal(model: &mut Model, numerator: i128, denominator: u128) -> Result<Self> {
        Self::from_wide(
            model,
            Integer::from_i128(numerator),
            Magnitude::from_small(denominator),
        )
    }

    pub(super) fn from_wide(
        model: &mut Model,
        numerator: Integer,
        denominator: Magnitude,
    ) -> Result<Self> {
        if denominator.is_zero() {
            return Err(Error::plain("Q/Z denominator must be positive"));
        }
        let (_, residue) = numerator.magnitude.div_rem(model, denominator)?;
        if residue.is_zero() {
            return Ok(Self::zero());
        }
        let residue = if numerator.negative {
            denominator.difference(residue)
        } else {
            residue
        };
        let common = residue.gcd(model, denominator)?;
        Ok(Self {
            numerator: narrow_quotient(model, residue, common)?,
            denominator: narrow_quotient(model, denominator, common)?,
        })
    }

    pub(crate) fn add(self, model: &mut Model, other: Self) -> Result<Self> {
        let a = self
            .denominator
            .resize::<WORK_LIMBS>()
            .ok_or_else(arithmetic)?;
        let b = other
            .denominator
            .resize::<WORK_LIMBS>()
            .ok_or_else(arithmetic)?;
        let common = a.gcd(model, b)?;
        let a_reduced = a.divide_exact(model, common)?.ok_or_else(arithmetic)?;
        let b_reduced = b.divide_exact(model, common)?.ok_or_else(arithmetic)?;
        let denominator = a.checked_product(b_reduced).ok_or_else(arithmetic)?;
        let numerator = self
            .scaled_numerator(b_reduced)?
            .add(other.scaled_numerator(a_reduced)?)?;
        Self::from_wide(model, numerator, denominator)
    }

    pub(crate) fn negate(self) -> Self {
        if self.numerator.is_zero() {
            return self;
        }
        Self {
            numerator: self.denominator.difference(self.numerator),
            denominator: self.denominator,
        }
    }

    pub(crate) fn scale(self, model: &mut Model, scalar: i128) -> Result<Self> {
        let numerator = self.scaled_numerator(Magnitude::from_small(scalar.unsigned_abs()))?;
        let numerator = if scalar < 0 {
            numerator.negate()
        } else {
            numerator
        };
        Self::from_wide(
            model,
            numerator,
            self.denominator.resize().ok_or_else(arithmetic)?,
        )
    }

    pub(crate) fn compare(self, model: &mut Model, other: Self) -> Result<Ordering> {
        model.spend(WORK_LIMBS)?;
        let left = self.scaled_numerator(other.denominator.resize().ok_or_else(arithmetic)?)?;
        let right = other.scaled_numerator(self.denominator.resize().ok_or_else(arithmetic)?)?;
        Ok(left.magnitude.cmp(right.magnitude))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(2 * LIMBS)?;
        Ok(Value::QZ(self))
    }

    pub(crate) fn root_order(self) -> Option<u128> {
        self.denominator.small()
    }

    pub(crate) fn root_index(self, order: usize) -> Option<usize> {
        let denominator = self.denominator.small()?;
        let order = order as u128;
        if !order.is_multiple_of(denominator) {
            return None;
        }
        Some((self.numerator.small()? * (order / denominator)) as usize)
    }

    pub(crate) fn order(self) -> String {
        decimal(self.denominator)
    }

    pub(super) fn order_matches(self, bound: u128) -> bool {
        self.denominator.small() == Some(bound)
    }

    pub(crate) fn text(self) -> String {
        if self.numerator.is_zero() {
            "0".into()
        } else {
            format!("{}/{}", decimal(self.numerator), decimal(self.denominator))
        }
    }

    fn scaled_numerator(self, factor: Magnitude) -> Result<Integer> {
        let numerator: Magnitude = self.numerator.resize().ok_or_else(arithmetic)?;
        let magnitude = numerator.checked_product(factor).ok_or_else(arithmetic)?;
        Ok(Integer {
            magnitude,
            negative: false,
        })
    }
}

fn narrow_quotient(
    model: &mut Model,
    value: Magnitude,
    divisor: Magnitude,
) -> Result<Natural<LIMBS>> {
    value
        .divide_exact(model, divisor)?
        .and_then(Natural::resize)
        .ok_or_else(arithmetic)
}

fn arithmetic() -> Error {
    Error::plain("Q/Z scalar exceeds the supported exact coefficient range")
}
