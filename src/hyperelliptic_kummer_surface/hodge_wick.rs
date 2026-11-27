use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::hw_algebra::{self, Matrix};
use crate::model::{Model, Value, MAX_SHOWN_BYTES};
use crate::syntax::{Label, ParsedExpr};

use super::nakajima::Grade;

mod field;
mod moments;
mod rational;
pub(crate) mod ring;
use field::Etale;
use rational::{Rational, LIMBS, Q};

pub(crate) type Entry = (i128, i128);

pub(crate) enum Command {
    Declare(Label, usize, Vec<i128>, Vec<Vec<Entry>>),
    Trace(Label, Vec<Entry>),
    Pair(bool, Box<[ParsedExpr; 2]>),
    Integral(bool, Box<ParsedExpr>),
    Ring(Box<ring::Command>),
}

pub(crate) enum Expression {
    Unit(Label, u128),
    Divisor(Label, u128, usize),
    Quadratic(Label, u128, Vec<Entry>),
    Inverse(Label, u128),
    Scale(bool, Entry, Box<ParsedExpr>),
    Cup(bool, Box<[ParsedExpr; 2]>),
    Cast(bool, Box<ParsedExpr>),
    Basis(Label, i128, i128, usize),
    Primitive(i128, usize, Box<ParsedExpr>),
}

pub(crate) struct Form {
    gram: Matrix<Rational>,
    inverse: Matrix<Rational>,
    field: Etale,
    multiplicity: usize,
    dimension: usize,
    slices: RefCell<BTreeMap<(i128, i128), Rc<ring::Slice>>>,
    koszul: ring::koszul::Cache,
    resolutions: ring::resolution::Cache,
    hochschild: ring::hochschild::Cache,
}

impl Form {
    pub(crate) fn new(
        model: &mut Model,
        m: usize,
        polynomial: Vec<i128>,
        gram: Vec<Vec<Entry>>,
    ) -> Result<Self> {
        if polynomial.len() < 2 || polynomial.last() != Some(&1) {
            return Err(Error::plain(
                "Hodge Wick étale polynomial must be monic of positive degree",
            ));
        }
        let rho = gram.len();
        let dimension = m
            .checked_mul(polynomial.len() - 1)
            .and_then(|t| t.checked_add(rho));
        let Some(dimension @ (22 | 23)) = dimension.filter(|_| m > 0) else {
            return Err(Error::plain(
                "Hodge Wick requires positive étale rank m and total dimension 22 or 23",
            ));
        };
        if gram.iter().any(|row| row.len() != rho) {
            return Err(Error::plain(
                "Hodge Wick divisor Gram matrix must be square",
            ));
        }
        let entries = gram
            .into_iter()
            .flatten()
            .map(|entry| Rational::input(model, entry))
            .collect::<Result<_>>()?;
        let gram = Matrix::from_entries(rho, entries)?;
        for i in 0..rho {
            for j in 0..i {
                if gram.entry(i, j) != gram.entry(j, i) {
                    return Err(Error::plain(
                        "Hodge Wick divisor Gram matrix must be symmetric",
                    ));
                }
            }
        }
        let inverse = hw_algebra::inverse(model, &Q, &gram)?
            .ok_or_else(|| Error::plain("Hodge Wick divisor Gram matrix is singular"))?;
        let field = Etale::new(model, polynomial)?;
        model.allocate(2 * rho * rho * (2 * LIMBS + 1) + 4 * (field.degree() + 1))?;
        Ok(Self {
            gram,
            inverse,
            field,
            multiplicity: m,
            dimension,
            slices: RefCell::default(),
            koszul: RefCell::default(),
            resolutions: RefCell::default(),
            hochschild: RefCell::default(),
        })
    }

    pub(crate) fn validate_grade(&self, grade: &Grade) -> Result<()> {
        if grade.length < 1 {
            return Err(Error::plain("Hodge Wick Hilbert length must be positive"));
        }
        let expected = if grade.length == 1 { 22 } else { 23 };
        if self.dimension != expected {
            return Err(Error::plain(format!(
                "Hodge Wick length {} requires total Beauville rank {expected}",
                grade.length
            )));
        }
        if grade.degree < 0 || grade.degree % 2 != 0 {
            return Err(Error::plain(
                "Hodge Wick cohomological degree must be a nonnegative even integer",
            ));
        }
        Ok(())
    }

    pub(crate) fn divisor(&self, index: usize) -> Result<Tensor> {
        if index >= self.gram.size() {
            return Err(Error::plain(
                "Hodge Wick divisor index is outside the marked Gram matrix",
            ));
        }
        Ok(Tensor::variable(index as u8))
    }

    pub(crate) fn quadratic(&self, model: &mut Model, entries: &[Entry]) -> Result<Tensor> {
        let mut tensor = Tensor::default();
        for (i, entry) in self.field.element(model, entries)?.into_iter().enumerate() {
            tensor.insert(
                model,
                vec![Power {
                    index: (self.gram.size() + i) as u8,
                    exponent: 1,
                }],
                entry,
            )?;
        }
        Ok(tensor)
    }

    pub(crate) fn inverse(&self, model: &mut Model) -> Result<Tensor> {
        let mut tensor = Tensor::variable(self.gram.size() as u8);
        for i in 0..self.gram.size() {
            for j in i..self.gram.size() {
                let mut entry = *self.inverse.entry(i, j);
                let mut powers = vec![Power {
                    index: i as u8,
                    exponent: 1,
                }];
                if i == j {
                    powers[0].exponent = 2;
                } else {
                    powers.push(Power {
                        index: j as u8,
                        exponent: 1,
                    });
                    entry = entry.add(model, entry)?;
                }
                tensor.insert(model, powers, entry)?;
            }
        }
        Ok(tensor)
    }

    pub(crate) fn trace(&self, model: &mut Model, entries: &[Entry]) -> Result<String> {
        self.field.trace(model, entries)?.text(model)
    }

    pub(crate) fn integral(&self, model: &mut Model, tensor: &Tensor) -> Result<String> {
        moments::integral(model, self, tensor)?.text(model)
    }
}

#[derive(Clone, Copy, Eq, PartialEq, Ord, PartialOrd)]
struct Power {
    index: u8,
    exponent: u128,
}

#[derive(Clone, Default)]
pub(crate) struct Tensor {
    terms: BTreeMap<Vec<Power>, Rational>,
    cells: usize,
}

impl Tensor {
    pub(crate) fn unit() -> Self {
        Self::monomial(Vec::new())
    }

    fn variable(index: u8) -> Self {
        Self::monomial(vec![Power { index, exponent: 1 }])
    }

    fn monomial(powers: Vec<Power>) -> Self {
        Self {
            cells: term_cells(&powers),
            terms: BTreeMap::from([(powers, Rational::integer(1))]),
        }
    }

    fn insert(&mut self, model: &mut Model, powers: Vec<Power>, value: Rational) -> Result<()> {
        if value.is_zero() {
            return Ok(());
        }
        model.spend(1 + powers.len())?;
        if let Some(previous) = self.terms.get_mut(&powers) {
            *previous = previous.add(model, value)?;
            if previous.is_zero() {
                self.terms.remove(&powers);
                self.cells -= term_cells(&powers);
            }
        } else {
            let cells = self.cells.saturating_add(term_cells(&powers));
            storage(cells)?;
            self.terms.insert(powers, value);
            self.cells = cells;
        }
        Ok(())
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(self.cells + 1)?;
        Ok(Value::HodgeTensor(Rc::new(self)))
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let mut sum = Self::default();
        for (powers, &coefficient) in self.terms.iter().chain(&rhs.terms) {
            sum.insert(model, powers.clone(), coefficient)?;
        }
        Ok(sum)
    }

    pub(crate) fn scale(&self, model: &mut Model, scalar: Entry) -> Result<Self> {
        let scalar = Rational::input(model, scalar)?;
        let mut result = Self::default();
        if !scalar.is_zero() {
            for (powers, &coefficient) in &self.terms {
                let value = coefficient.multiply(model, scalar)?;
                result.insert(model, powers.clone(), value)?;
            }
        }
        Ok(result)
    }

    pub(crate) fn cup(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        model
            .spend(self.terms.len().saturating_mul(rhs.terms.len()))
            .map_err(|_| work_limit())?;
        let mut result = Self::default();
        for (left, &a) in &self.terms {
            for (right, &b) in &rhs.terms {
                let coefficient = a.multiply(model, b)?;
                result.insert(model, multiply_powers(left, right)?, coefficient)?;
            }
        }
        Ok(result)
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        for ((left, a), (right, b)) in self.terms.iter().zip(&rhs.terms) {
            model.spend(1 + left.len() + right.len())?;
            let order = left.cmp(right).then_with(|| a.compare(*b));
            if !order.is_eq() {
                return Ok(order);
            }
        }
        Ok(self.terms.len().cmp(&rhs.terms.len()))
    }

    pub(crate) fn text(&self, model: &mut Model, form: &Form) -> Result<String> {
        let mut out = "{".to_owned();
        for (position, (powers, coefficient)) in self.terms.iter().enumerate() {
            if position > 0 {
                out.push_str(", ");
            }
            out.push_str(&coefficient.text(model)?);
            for power in powers {
                let i = usize::from(power.index);
                if i < form.gram.size() {
                    write!(out, "*D[{i}]").unwrap();
                } else {
                    write!(out, "*Q[{}]", i - form.gram.size()).unwrap();
                }
                if power.exponent != 1 {
                    write!(out, "^{}", power.exponent).unwrap();
                }
            }
            if out.len() >= MAX_SHOWN_BYTES {
                return Err(Error::plain("Hodge Wick tensor exceeds the report limit"));
            }
        }
        out.push('}');
        Ok(out)
    }
}

fn multiply_powers(left: &[Power], right: &[Power]) -> Result<Vec<Power>> {
    let mut left = left.iter().peekable();
    let mut right = right.iter().peekable();
    let mut result = Vec::new();
    while let (Some(a), Some(b)) = (left.peek(), right.peek()) {
        match a.index.cmp(&b.index) {
            Ordering::Less => result.push(*left.next().unwrap()),
            Ordering::Greater => result.push(*right.next().unwrap()),
            Ordering::Equal => {
                let a = left.next().unwrap();
                let b = right.next().unwrap();
                let exponent = a.exponent.checked_add(b.exponent).ok_or_else(work_limit)?;
                result.push(Power {
                    index: a.index,
                    exponent,
                });
            }
        }
    }
    result.extend(left.chain(right));
    Ok(result)
}

fn term_cells(powers: &[Power]) -> usize {
    2 * LIMBS + 8 + 5 * powers.len()
}

fn storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        Err(work_limit())
    } else {
        Ok(())
    }
}

fn work_limit() -> Error {
    Error::plain("Hodge Wick computation exceeds the work or storage limit")
}
