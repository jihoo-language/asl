use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::rc::Rc;

use crate::elliptic_modular::decimal;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_qp_quad::Wide;
use crate::rv_witt_spectrum::norm_power;
use crate::syntax::ParsedExpr;
use crate::wide_integer::Natural;

mod diagram;
mod dimension;
pub(crate) mod karoubi;
pub(crate) mod local;
pub(crate) mod square;
use diagram::Diagram;
pub(crate) use dimension::Report;

const MAX_RANK: usize = 32;
const MAX_TERMS: usize = 100_000;
type Integer = Wide<16>;
type Monomial = (Diagram, u128);
type Polynomial = BTreeMap<u128, Integer>;

pub(crate) enum Expression {
    Basis(usize, Vec<u128>, bool),
    Identity(usize),
    Unary(Unary, Box<ParsedExpr>),
    Binary(Binary, Box<[ParsedExpr; 2]>),
    Localize(Box<ParsedExpr>),
    Divide(local::Denominator, Box<ParsedExpr>),
    SquareUnit(usize, usize),
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Scale(i128, u128),
    Reflect,
}

#[derive(Clone, Copy)]
pub(crate) enum Binary {
    Compose,
    Tensor,
}

pub(crate) enum Query {
    Coefficient(Vec<u128>, u128),
    Entry(u128, Vec<u128>, Vec<u128>),
    Trace,
    TraceAt(u128),
    ZeroAt(u128),
}

pub(crate) fn require_rank(rank: u128) -> Result<usize> {
    if rank > MAX_RANK as u128 {
        return Err(Error::plain(
            "partition diagram rank exceeds the limit of 32",
        ));
    }
    Ok(rank as usize)
}

impl Binary {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Compose => "rvPartitionCompose",
            Self::Tensor => "rvPartitionTensor",
        }
    }

    pub(crate) fn rank(self, a: usize, b: usize) -> Result<usize> {
        match self {
            Self::Compose if a == b => Ok(a),
            Self::Compose => Err(Error::plain("partition composition requires equal ranks")),
            Self::Tensor => require_rank((a + b) as u128),
        }
    }
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Algebra {
    rank: usize,
    terms: BTreeMap<Monomial, i128>,
}

impl Algebra {
    pub(crate) fn zero(rank: usize) -> Self {
        Self {
            rank,
            terms: BTreeMap::new(),
        }
    }

    pub(crate) fn basis(
        model: &mut Model,
        rank: usize,
        labels: Vec<u128>,
        orbit: bool,
    ) -> Result<Self> {
        model.spend(labels.len())?;
        let diagram = Diagram::literal(rank, labels)?;
        let terms = if orbit {
            diagram.orbit(model)?
        } else {
            BTreeMap::from([((diagram, 0), 1)])
        };
        Ok(Self { rank, terms })
    }

    pub(crate) fn identity(rank: usize) -> Self {
        let diagram = Diagram::identity(rank);
        Self {
            rank,
            terms: BTreeMap::from([((diagram, 0), 1)]),
        }
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(4 + (8 + 2 * self.rank) * self.terms.len())?;
        Ok(Value::RVPartitionAlgebra(Rc::new(self)))
    }

    pub(crate) fn render(&self) -> String {
        let mut out = format!("RVPartitionAlgebra({}){{", self.rank);
        for (i, ((diagram, degree), coefficient)) in self.terms.iter().enumerate() {
            if i != 0 {
                out.push_str(", ");
            }
            write!(out, "({:?},{degree},{coefficient})", diagram.labels()).unwrap();
        }
        out.push('}');
        out
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend((1 + 2 * self.rank) * self.terms.len().min(rhs.terms.len()))?;
        Ok(self.cmp(rhs))
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        debug_assert_eq!(self.rank, rhs.rank);
        model.spend((1 + 2 * self.rank) * (self.terms.len() + rhs.terms.len()))?;
        let mut terms = self.terms.clone();
        for (key, &value) in &rhs.terms {
            let sum = terms
                .get(key)
                .copied()
                .unwrap_or(0)
                .checked_add(value)
                .ok_or_else(coefficient_overflow)?;
            if sum == 0 {
                terms.remove(key);
            } else {
                terms.insert(key.clone(), sum);
                require_terms(terms.len())?;
            }
        }
        Ok(Self {
            rank: self.rank,
            terms,
        })
    }

    pub(crate) fn unary(&self, model: &mut Model, operation: Unary) -> Result<Self> {
        if matches!(operation, Unary::Scale(0, _)) {
            return Ok(Self::zero(self.rank));
        }
        model.spend((1 + 2 * self.rank) * self.terms.len())?;
        let terms = self
            .terms
            .iter()
            .map(|((diagram, degree), &coefficient)| match operation {
                Unary::Reflect => Ok(((diagram.reflect(self.rank), *degree), coefficient)),
                Unary::Scale(scale, shift) => {
                    let degree = add_degrees(*degree, shift)?;
                    let value = coefficient
                        .checked_mul(scale)
                        .ok_or_else(coefficient_overflow)?;
                    Ok(((diagram.clone(), degree), value))
                }
            })
            .collect::<Result<_>>()?;
        Ok(Self {
            rank: self.rank,
            terms,
        })
    }

    pub(crate) fn binary(&self, model: &mut Model, operation: Binary, rhs: &Self) -> Result<Self> {
        let rank = operation.rank(self.rank, rhs.rank)?;
        let mut evaluated = BTreeMap::new();
        for ((left, i), &a) in &self.terms {
            for ((right, j), &b) in &rhs.terms {
                model.spend(4 + 4 * (self.rank + rhs.rank))?;
                let (diagram, loops) = match operation {
                    Binary::Compose => left.compose(right, rank),
                    Binary::Tensor => (left.tensor(self.rank, right, rhs.rank), 0),
                };
                let degree = add_degrees(add_degrees(*i, *j)?, loops as u128)?;
                accumulate(&mut evaluated, (diagram, degree), Integer::product(a, b))?;
            }
        }
        let terms = evaluated
            .into_iter()
            .map(|(key, value)| Ok((key, value.signed().ok_or_else(coefficient_overflow)?)))
            .collect::<Result<_>>()?;
        Ok(Self { rank, terms })
    }

    pub(crate) fn query(&self, model: &mut Model, query: Query) -> Result<String> {
        match query {
            Query::Coefficient(labels, degree) => {
                model.spend(labels.len())?;
                let diagram = Diagram::literal(self.rank, labels)?;
                let key = (diagram, degree);
                let coefficient = self.terms.get(&key).copied().unwrap_or(0);
                Ok(format!(
                    "rvPartitionCoefficient[{:?}, {degree}] = {coefficient}",
                    key.0.labels()
                ))
            }
            Query::Trace => {
                let polynomial = self.trace_polynomial(model)?;
                Ok(format!(
                    "rvPartitionTrace = {}",
                    polynomial_report(polynomial)
                ))
            }
            Query::TraceAt(n) => {
                let polynomial = self.trace_polynomial(model)?;
                let value = evaluate_polynomial(model, polynomial, n)?;
                Ok(format!("rvPartitionTraceAt[{n}] = {}", decimal(value)))
            }
            Query::ZeroAt(n) => Ok(format!(
                "rvPartitionZeroAt[{n}] = {}",
                self.zero_at(model, n)?
            )),
            Query::Entry(dimension, upper, lower) => {
                let value = self.entry(model, dimension, upper, lower)?;
                Ok(format!(
                    "rvPartitionEntry[{dimension}] = {}",
                    decimal(value)
                ))
            }
        }
    }

    fn trace_polynomial(&self, model: &mut Model) -> Result<Polynomial> {
        let mut polynomial = Polynomial::new();
        for ((diagram, degree), &value) in &self.terms {
            model.spend(1 + 3 * self.rank)?;
            let degree = add_degrees(*degree, diagram.closure_components(self.rank) as u128)?;
            accumulate(&mut polynomial, degree, Integer::from_i128(value))?;
        }
        Ok(polynomial)
    }

    fn entry(
        &self,
        model: &mut Model,
        n: u128,
        mut upper: Vec<u128>,
        lower: Vec<u128>,
    ) -> Result<Integer> {
        if upper.len() != self.rank || lower.len() != self.rank {
            return Err(Error::plain(
                "partition matrix index has the wrong tensor rank",
            ));
        }
        upper.extend(lower);
        if upper.iter().any(|&index| index >= n) {
            return Err(Error::plain("partition matrix index lies outside 0..N"));
        }
        self.evaluate_entry(model, n, &upper)
    }

    fn evaluate_entry(&self, model: &mut Model, n: u128, indices: &[u128]) -> Result<Integer> {
        let mut polynomial = Polynomial::new();
        for ((diagram, degree), &value) in &self.terms {
            model.spend(1 + 2 * self.rank)?;
            if diagram.supports(indices) {
                accumulate(&mut polynomial, *degree, Integer::from_i128(value))?;
            }
        }
        evaluate_polynomial(model, polynomial, n)
    }

    fn zero_at(&self, model: &mut Model, n: u128) -> Result<bool> {
        if self.terms.is_empty() || (n == 0 && self.rank != 0) {
            return Ok(true);
        }
        if self.rank > 4 {
            return Err(Error::plain("partition matrix-zero test exceeds rank 4"));
        }
        diagram::all_patterns(model, 2 * self.rank, n, &mut |model, indices| {
            Ok(self.evaluate_entry(model, n, indices)?.is_zero())
        })
    }
}

fn polynomial_report(polynomial: Polynomial) -> String {
    let terms: Vec<_> = polynomial
        .into_iter()
        .map(|(degree, value)| format!("({degree},{})", decimal(value)))
        .collect();
    format!("[{}]", terms.join(", "))
}

fn evaluate_polynomial(
    model: &mut Model,
    polynomial: Polynomial,
    dimension: u128,
) -> Result<Integer> {
    let mut descending = polynomial.into_iter().rev();
    let Some((mut previous, mut value)) = descending.next() else {
        return Ok(Integer::from_i128(0));
    };
    let base = natural(dimension);
    for (degree, coefficient) in descending {
        if value.is_zero() {
            value = coefficient;
        } else {
            let place = norm_power(model, base, previous - degree)?;
            value = value
                .multiply_wide(place)
                .and_then(|value| value.add(coefficient))
                .map_err(|_| intermediate_overflow())?;
        }
        previous = degree;
    }
    if value.is_zero() {
        return Ok(value);
    }
    let place = norm_power(model, base, previous)?;
    value
        .multiply_wide(place)
        .map_err(|_| intermediate_overflow())
}

fn natural(value: u128) -> Integer {
    Integer {
        magnitude: Natural::from_small(value),
        negative: false,
    }
}

fn accumulate<K: Ord>(map: &mut BTreeMap<K, Integer>, key: K, coefficient: Integer) -> Result<()> {
    let value = map
        .get(&key)
        .copied()
        .unwrap_or_else(|| Integer::from_i128(0));
    let value = value
        .add(coefficient)
        .map_err(|_| intermediate_overflow())?;
    if value.is_zero() {
        map.remove(&key);
    } else {
        map.insert(key, value);
        require_terms(map.len())?;
    }
    Ok(())
}

fn require_terms(count: usize) -> Result<()> {
    if count > MAX_TERMS {
        return Err(Error::plain(
            "partition diagram algebra exceeds the term storage limit",
        ));
    }
    Ok(())
}

fn add_degrees(a: u128, b: u128) -> Result<u128> {
    a.checked_add(b)
        .ok_or_else(|| Error::plain("partition diagram t degree exceeds unsigned 128-bit range"))
}

fn coefficient_overflow() -> Error {
    Error::plain("partition diagram coefficient exceeds signed 128-bit range")
}

fn intermediate_overflow() -> Error {
    Error::plain("partition diagram exact intermediate exceeds 512 bits")
}
