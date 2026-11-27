use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::syntax::ParsedExpr;

use super::{
    degree_limit, lift, narrow, one, product, Polynomial, SchurShape, Spectrum, MAX_DEGREE,
};

mod plethysm;

pub(crate) enum Expression {
    Seed(SchurShape),
    Scale(i128, Box<ParsedExpr>),
    Plethysm(Box<[ParsedExpr; 2]>),
    Evaluate(Box<[ParsedExpr; 2]>),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Character {
    weight: u128,
    terms: BTreeMap<SchurShape, i128>,
}

impl Character {
    pub(crate) fn seed(shape: SchurShape) -> Self {
        Self {
            weight: shape.weight(),
            terms: BTreeMap::from([(shape, 1)]),
        }
    }

    pub(crate) fn zero(weight: u128) -> Self {
        Self {
            weight,
            terms: BTreeMap::new(),
        }
    }

    fn constant(value: i128) -> Result<Self> {
        let mut zero = Self::zero(0);
        if value != 0 {
            zero.terms.insert(SchurShape::new(Vec::new())?, value);
        }
        Ok(zero)
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(
            4 + self
                .terms
                .keys()
                .map(|shape| 6 + 2 * shape.rows().len())
                .sum::<usize>(),
        )?;
        Ok(Value::RVSchurCharacter(Rc::new(self)))
    }

    pub(crate) fn render(&self) -> String {
        let terms: Vec<_> = self
            .terms
            .iter()
            .rev()
            .map(|(shape, &coefficient)| format!("({:?},{coefficient})", shape.rows()))
            .collect();
        format!("RVSchurCharacter({}){{{}}}", self.weight, terms.join(", "))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        let rows = self
            .terms
            .keys()
            .zip(other.terms.keys())
            .map(|(a, b)| a.rows().len().min(b.rows().len()) + 1)
            .sum();
        model.spend(rows)?;
        Ok(self.cmp(other))
    }

    pub(crate) fn coefficient(&self, shape: &SchurShape) -> String {
        let coefficient = self.terms.get(shape).copied().unwrap_or(0);
        format!("rvSchurCharacterCoeff[{:?}] = {coefficient}", shape.rows())
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        debug_assert_eq!(self.weight, rhs.weight);
        model.spend(self.terms.len() + rhs.terms.len())?;
        let mut terms = self.terms.clone();
        for (shape, &coefficient) in &rhs.terms {
            let sum = terms
                .get(shape)
                .copied()
                .unwrap_or(0)
                .checked_add(coefficient)
                .ok_or_else(overflow)?;
            if sum == 0 {
                terms.remove(shape);
            } else {
                terms.insert(shape.clone(), sum);
            }
        }
        Ok(Self {
            weight: self.weight,
            terms,
        })
    }

    pub(crate) fn scale(&self, model: &mut Model, factor: i128) -> Result<Self> {
        if factor == 0 {
            return Ok(Self::zero(self.weight));
        }
        model.spend(self.terms.len())?;
        let terms = self
            .terms
            .iter()
            .map(|(shape, &value)| {
                Ok((
                    shape.clone(),
                    value.checked_mul(factor).ok_or_else(overflow)?,
                ))
            })
            .collect::<Result<_>>()?;
        Ok(Self {
            weight: self.weight,
            terms,
        })
    }

    pub(crate) fn evaluation(&self, model: &mut Model, input: &Spectrum) -> Result<Spectrum> {
        let mut sum = Spectrum::zero();
        for (shape, &coefficient) in &self.terms {
            let term = input.schur(model, shape)?;
            let value = multiple(model, &term, coefficient)?;
            sum = sum.add(model, &value)?;
        }
        Ok(sum)
    }
}

pub(crate) fn output_weight(outer: u128, inner: u128) -> Result<u128> {
    outer
        .checked_mul(inner)
        .ok_or_else(|| Error::plain("Schur plethysm weight exceeds unsigned 128-bit range"))
}

fn multiple(model: &mut Model, input: &Spectrum, coefficient: i128) -> Result<Spectrum> {
    let exponent = coefficient.unsigned_abs();
    let mut evaluate = |factor: &[i128]| -> Result<Box<[i128]>> {
        if factor.len() == 1 || exponent == 1 {
            return Ok(factor.into());
        }
        if exponent > (MAX_DEGREE / (factor.len() - 1)) as u128 {
            return Err(degree_limit());
        }
        let mut base = lift(factor);
        let mut exponent = exponent;
        let mut result: Polynomial = vec![one()];
        while exponent != 0 {
            if !exponent.is_multiple_of(2) {
                result = product(model, &result, &base)?;
            }
            exponent /= 2;
            if exponent != 0 {
                base = product(model, &base, &base)?;
            }
        }
        narrow(result)
    };
    // P and Q are already relatively prime, and powers preserve that fact.
    let value = Spectrum {
        numerator: evaluate(&input.numerator)?,
        denominator: evaluate(&input.denominator)?,
    };
    Ok(if coefficient < 0 {
        value.negate()
    } else {
        value
    })
}

fn overflow() -> Error {
    Error::plain("Schur character coefficient exceeds signed 128-bit range")
}
