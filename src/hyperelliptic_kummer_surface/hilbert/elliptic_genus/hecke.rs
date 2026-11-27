use std::cmp::Ordering;
use std::collections::{btree_map::Entry, BTreeMap};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_ore::gcd;
use crate::syntax::ParsedExpr;

pub(crate) use super::decimal;
use super::{add, hurwitz, multiply, Integer, Magnitude, INTEGER_WORDS, TERM_LIMIT};

#[derive(Clone, Copy)]
pub(crate) enum Pair {
    Multiply,
    SameAction,
}

pub(crate) enum Report {
    Trace,
    Characteristic,
    CharacteristicCoefficient(u128),
}

impl Report {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            Self::Trace => "rvHeckeActionTrace",
            Self::Characteristic => "rvHeckeChar",
            Self::CharacteristicCoefficient(_) => "rvHeckeCharCoeff",
        }
    }

    pub(crate) fn write(
        self,
        model: &mut Model,
        weight: u128,
        operator: &Operator,
    ) -> Result<String> {
        let name = self.name();
        match self {
            Self::Trace => Ok(format!(
                "{name} = {}",
                decimal(operator.trace(model, weight)?)
            )),
            Self::Characteristic => {
                let coefficients = operator.characteristic(model, weight, dimension(weight))?;
                let coefficients: Vec<_> = coefficients.into_iter().map(decimal).collect();
                Ok(format!("{name} = [{}]", coefficients.join(", ")))
            }
            Self::CharacteristicCoefficient(j) => Ok(format!(
                "{name}[{j}] = {}",
                decimal(operator.characteristic_coefficient(model, weight, j)?)
            )),
        }
    }
}

pub(crate) enum Expression {
    Basis(u128, u128),
    Scale(i128, Box<ParsedExpr>),
    SameAction(Box<[ParsedExpr; 2]>),
}

pub(super) type Support = BTreeMap<u128, Integer>;

#[derive(Clone)]
pub(crate) struct Operator(Support);

impl Operator {
    pub(crate) fn zero() -> Self {
        Self(Support::new())
    }

    pub(crate) fn basis(index: u128) -> Result<Self> {
        if index == 0 {
            return Err(Error::plain("Hecke index n must be positive"));
        }
        Ok(Self(Support::from([(index, integer(1))])))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + self.0.len() * (INTEGER_WORDS + 8))?;
        Ok(Value::ModularHecke(Rc::new(self)))
    }

    pub(crate) fn terms(&self) -> impl Iterator<Item = (u128, Integer)> + '_ {
        self.0
            .iter()
            .map(|(&index, &coefficient)| (index, coefficient))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(1 + INTEGER_WORDS * self.0.len().min(rhs.0.len()))?;
        for ((n, a), (m, b)) in self.terms().zip(rhs.terms()) {
            let ordering = n.cmp(&m).then_with(|| compare_integers(a, b));
            if !ordering.is_eq() {
                return Ok(ordering);
            }
        }
        Ok(self.0.len().cmp(&rhs.0.len()))
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        self.sum(model, rhs, false)
    }

    fn sum(&self, model: &mut Model, rhs: &Self, subtract: bool) -> Result<Self> {
        model.spend(1 + INTEGER_WORDS * self.0.len())?;
        let mut result = self.0.clone();
        for (n, a) in rhs.terms() {
            add(model, &mut result, n, if subtract { a.negate() } else { a })?;
        }
        Ok(Self(result))
    }

    pub(crate) fn scale(&self, model: &mut Model, scalar: Integer) -> Result<Self> {
        if scalar.is_zero() {
            return Ok(Self::zero());
        }
        let result: Result<_> = self
            .terms()
            .map(|(n, a)| Ok((n, multiply(model, a, scalar)?)))
            .collect();
        Ok(Self(result?))
    }

    pub(super) fn scalar(&self) -> Option<Integer> {
        match self.0.len() {
            0 => Some(integer(0)),
            1 => self.0.get(&1).copied(),
            _ => None,
        }
    }

    pub(crate) fn multiply(&self, model: &mut Model, weight: u128, rhs: &Self) -> Result<Self> {
        if let Some(scalar) = self.scalar() {
            return rhs.scale(model, scalar);
        }
        if let Some(scalar) = rhs.scalar() {
            return self.scale(model, scalar);
        }
        let mut result = Support::new();
        for (n, a) in self.terms() {
            for (m, b) in rhs.terms() {
                let product = n.checked_mul(m).ok_or_else(|| {
                    Error::plain(
                        "level-one Hecke operator index product exceeds unsigned 128-bit range",
                    )
                })?;
                let coefficient = multiply(model, a, b)?;
                let common = gcd(n, m);
                for small in 1..=common.isqrt() {
                    model.spend(1)?;
                    if common.is_multiple_of(small) {
                        let large = common / small;
                        add_product(model, &mut result, weight, product, small, coefficient)?;
                        if large != small {
                            add_product(model, &mut result, weight, product, large, coefficient)?;
                        }
                    }
                }
            }
        }
        Ok(Self(result))
    }

    pub(crate) fn same_action(&self, model: &mut Model, weight: u128, rhs: &Self) -> Result<bool> {
        if dimension(weight) == 0 {
            return Ok(true);
        }
        let difference = self.sum(model, rhs, true)?;
        if difference.0.is_empty() {
            return Ok(true);
        }
        if difference.scalar().is_some() {
            return Ok(false);
        }
        let square = difference.multiply(model, weight, &difference)?;
        Ok(square.trace(model, weight)?.is_zero())
    }

    fn trace(&self, model: &mut Model, weight: u128) -> Result<Integer> {
        self.trace_cached(model, weight, &mut Support::new())
    }

    pub(super) fn trace_cached(
        &self,
        model: &mut Model,
        weight: u128,
        traces: &mut Support,
    ) -> Result<Integer> {
        let dimension = dimension(weight);
        let mut sum = integer(0);
        if dimension == 0 {
            return Ok(sum);
        }
        for (n, a) in self.terms() {
            let full = traces.len() >= TERM_LIMIT;
            let trace = if n == 1 {
                unsigned(dimension)
            } else {
                match traces.entry(n) {
                    Entry::Occupied(entry) => *entry.get(),
                    Entry::Vacant(entry) => {
                        if full {
                            return Err(capacity());
                        }
                        *entry.insert(hurwitz::trace(model, weight, n)?)
                    }
                }
            };
            sum = sum
                .add(multiply(model, a, trace)?)
                .map_err(|_| overflow())?;
        }
        Ok(sum)
    }

    fn characteristic_coefficient(
        &self,
        model: &mut Model,
        weight: u128,
        j: u128,
    ) -> Result<Integer> {
        let degree = dimension(weight);
        if j > degree {
            return Ok(integer(0));
        }
        if j == 0 {
            return Ok(integer(1));
        }
        if let Some(scalar) = self.scalar() {
            return scalar_coefficient(model, degree, j, scalar);
        }
        Ok(*self.characteristic(model, weight, j)?.last().unwrap())
    }

    pub(super) fn characteristic(
        &self,
        model: &mut Model,
        weight: u128,
        j: u128,
    ) -> Result<Vec<Integer>> {
        let count = prefix_capacity(j)?;
        let mut coefficients = Vec::with_capacity(count + 1);
        coefficients.push(integer(1));
        if count == 0 {
            return Ok(coefficients);
        }
        if let Some(scalar) = self.scalar() {
            for i in 1..=j {
                let next = if scalar.is_zero() {
                    model.spend(1)?;
                    integer(0)
                } else {
                    let next =
                        binomial_step(model, *coefficients.last().unwrap(), dimension(weight), i)?;
                    multiply(model, next, scalar.negate())?
                };
                coefficients.push(next);
            }
            return Ok(coefficients);
        }

        let mut basis_traces = Support::new();
        let mut traces = Vec::with_capacity(count);
        let mut power = self.clone();
        for i in 1..=count {
            traces.push(power.trace_cached(model, weight, &mut basis_traces)?);
            let mut sum = integer(0);
            for (&coefficient, &trace) in coefficients.iter().zip(traces.iter().rev()) {
                sum = sum
                    .add(multiply(model, coefficient, trace)?)
                    .map_err(|_| overflow())?;
            }
            coefficients.push(divide(model, sum.negate(), i as u128)?);
            if i < count {
                power = power.multiply(model, weight, self)?;
            }
        }
        Ok(coefficients)
    }
}

pub(crate) fn validate_weight(weight: u128) -> Result<()> {
    if weight < 2 || !weight.is_multiple_of(2) {
        return Err(Error::plain(
            "RVHecke requires even level-one weight at least 2",
        ));
    }
    Ok(())
}

pub(super) fn dimension(weight: u128) -> u128 {
    // S_k = Delta M_{k-12}, with M_* = C[E_4,E_6].
    (weight / 12).saturating_sub(u128::from(weight % 12 == 2))
}

pub(crate) fn integer(value: i128) -> Integer {
    Integer::from_i128(value)
}

pub(super) fn compare_integers(a: Integer, b: Integer) -> Ordering {
    match (a.negative, b.negative) {
        (true, true) => b.magnitude.cmp(a.magnitude),
        (false, false) => a.magnitude.cmp(b.magnitude),
        (false, true) => Ordering::Greater,
        (true, false) => Ordering::Less,
    }
}

fn unsigned(value: u128) -> Integer {
    Integer {
        magnitude: Magnitude::from_small(value),
        negative: false,
    }
}

fn add_product(
    model: &mut Model,
    result: &mut Support,
    weight: u128,
    product: u128,
    divisor: u128,
    coefficient: Integer,
) -> Result<()> {
    let index = (product / divisor) / divisor;
    let term = super::scale_by_power(model, coefficient, divisor, weight - 1)?;
    add(model, result, index, term)
}

fn scalar_coefficient(
    model: &mut Model,
    degree: u128,
    j: u128,
    scalar: Integer,
) -> Result<Integer> {
    if scalar.is_zero() {
        return Ok(integer(0));
    }
    let negative = !scalar.negative && j % 2 == 1;
    let mut binomial = integer(1);
    for i in 1..=j.min(degree - j) {
        binomial = binomial_step(model, binomial, degree, i)?;
    }
    let mut power = Integer {
        negative: false,
        ..scalar
    };
    let mut exponent = j;
    while exponent > 0 {
        if exponent % 2 == 1 {
            binomial = multiply(model, binomial, power)?;
        }
        exponent /= 2;
        if exponent > 0 {
            power = multiply(model, power, power)?;
        }
    }
    Ok(Integer {
        negative,
        ..binomial
    })
}

fn binomial_step(model: &mut Model, current: Integer, n: u128, i: u128) -> Result<Integer> {
    let numerator = n - i + 1;
    let common = gcd(numerator, i);
    let a = divide(model, current, i / common)?;
    multiply(model, a, unsigned(numerator / common))
}

fn divide(model: &mut Model, value: Integer, divisor: u128) -> Result<Integer> {
    let (magnitude, remainder) = value.magnitude.divide(model, divisor)?;
    if remainder != 0 {
        return Err(Error::plain(
            "Hecke characteristic coefficient division is inexact",
        ));
    }
    Ok(Integer { magnitude, ..value })
}

fn prefix_capacity(j: u128) -> Result<usize> {
    usize::try_from(j)
        .ok()
        .filter(|&count| count < TERM_LIMIT / 2)
        .ok_or_else(capacity)
}

fn capacity() -> Error {
    Error::plain("level-one Hecke calculation exceeds the coefficient storage limit")
}

fn overflow() -> Error {
    Error::plain("level-one Hecke exact intermediate exceeds 512 bits")
}
