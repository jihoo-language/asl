use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::syntax::ParsedExpr;

use super::{
    degree_limit, lift, narrow, one, product, Polynomial, SchurShape, Spectrum, MAX_DEGREE,
};

mod kronecker;
mod plethysm;
mod stable;
mod table;
pub(crate) use stable::StableCharacter;

pub(crate) enum Expression {
    Seed(SchurShape),
    Scale(i128, Box<ParsedExpr>),
    Plethysm(Box<[ParsedExpr; 2]>),
    Evaluate(Box<[ParsedExpr; 2]>),
    Kronecker(Box<[ParsedExpr; 2]>),
    StableSeed(SchurShape),
    StableScale(i128, Box<ParsedExpr>),
    StableProduct(Box<[ParsedExpr; 2]>),
    Specialize(u128, Box<ParsedExpr>),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Character {
    weight: u128,
    terms: Terms,
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
        store_terms(model, &self.terms)?;
        Ok(Value::RVSchurCharacter(Rc::new(self)))
    }

    pub(crate) fn render(&self) -> String {
        format!(
            "RVSchurCharacter({}){{{}}}",
            self.weight,
            render_terms(&self.terms)
        )
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        charge_compare(model, &self.terms, &other.terms)?;
        Ok(self.cmp(other))
    }

    pub(crate) fn coefficient(&self, shape: &SchurShape) -> String {
        coefficient_report("rvSchurCharacterCoeff", &self.terms, shape)
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        debug_assert_eq!(self.weight, rhs.weight);
        Ok(Self {
            weight: self.weight,
            terms: add_terms(model, &self.terms, &rhs.terms)?,
        })
    }

    pub(crate) fn scale(&self, model: &mut Model, factor: i128) -> Result<Self> {
        Ok(Self {
            weight: self.weight,
            terms: scale_terms(model, &self.terms, factor)?,
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

type Terms = BTreeMap<SchurShape, i128>;

fn store_terms(model: &mut Model, terms: &Terms) -> Result<()> {
    model.allocate(
        4 + terms
            .keys()
            .map(|shape| 6 + 2 * shape.rows().len())
            .sum::<usize>(),
    )
}

fn render_terms(terms: &Terms) -> String {
    terms
        .iter()
        .rev()
        .map(|(shape, &coefficient)| format!("({:?},{coefficient})", shape.rows()))
        .collect::<Vec<_>>()
        .join(", ")
}

fn charge_compare(model: &mut Model, a: &Terms, b: &Terms) -> Result<()> {
    let rows = a
        .keys()
        .zip(b.keys())
        .map(|(a, b)| a.rows().len().min(b.rows().len()) + 1)
        .sum();
    model.spend(rows)
}

fn coefficient_report(name: &str, terms: &Terms, shape: &SchurShape) -> String {
    let coefficient = terms.get(shape).copied().unwrap_or(0);
    format!("{name}[{:?}] = {coefficient}", shape.rows())
}

fn add_terms(model: &mut Model, left: &Terms, right: &Terms) -> Result<Terms> {
    model.spend(left.len() + right.len())?;
    let mut terms = left.clone();
    for (shape, &coefficient) in right {
        add_term(&mut terms, shape.clone(), coefficient)?;
    }
    Ok(terms)
}

fn add_term(terms: &mut Terms, shape: SchurShape, coefficient: i128) -> Result<()> {
    let sum = terms
        .get(&shape)
        .copied()
        .unwrap_or(0)
        .checked_add(coefficient)
        .ok_or_else(overflow)?;
    if sum == 0 {
        terms.remove(&shape);
    } else {
        terms.insert(shape, sum);
    }
    Ok(())
}

fn scale_terms(model: &mut Model, terms: &Terms, factor: i128) -> Result<Terms> {
    if factor == 0 {
        return Ok(Terms::new());
    }
    model.spend(terms.len())?;
    terms
        .iter()
        .map(|(shape, &value)| {
            Ok((
                shape.clone(),
                value.checked_mul(factor).ok_or_else(overflow)?,
            ))
        })
        .collect()
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
