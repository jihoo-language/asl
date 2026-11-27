use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value, MAX_SHOWN_BYTES};
use crate::syntax::{Label, ParsedExpr};

use super::super::{moments, storage, work_limit, Entry, Form, Rational, Tensor, LIMBS};
use super::{linear, Grade as RingGrade, Monomial};

mod algebra;
mod complex;
mod operations;
use algebra::Algebra;
use complex::Complex;

pub(crate) type Cache = RefCell<BTreeMap<i128, Rc<Algebra>>>;

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Diff,
    Class,
    Representative,
    Bv,
    FromTangent,
}

#[derive(Clone, Copy)]
pub(crate) enum Bilinear {
    Cup,
    Bracket,
    Insert,
}

pub(crate) enum Command {
    Dims(Label, i128, i128, i128),
    Coefficient(usize, Vec<usize>, Box<ParsedExpr>),
}

pub(crate) enum Expression {
    Term(Label, i128, i128, i128, usize, Vec<usize>, Entry),
    Basis(Label, i128, i128, i128, usize),
    Scale(Entry, Box<ParsedExpr>),
    Unary(Unary, Box<ParsedExpr>),
    Bilinear(Bilinear, Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) basis: Rc<str>,
    pub(crate) length: i128,
    pub(crate) arity: i128,
    pub(crate) shift: i128,
}

impl Grade {
    pub(crate) fn ring(&self, weight: i128) -> RingGrade {
        RingGrade {
            basis: self.basis.clone(),
            length: self.length,
            degree: 2 * weight,
        }
    }

    pub(crate) fn with_arity_offset(&self, offset: i128) -> Result<Self> {
        Ok(Self {
            arity: self.arity.checked_add(offset).ok_or_else(grading_range)?,
            ..self.clone()
        })
    }

    pub(crate) fn product_grade(&self, other: &Self, op: Bilinear) -> Result<Self> {
        if self.basis != other.basis || self.length != other.length {
            return Err(Error::plain(
                "Hochschild operands must have the same marking and Hilbert length",
            ));
        }
        let offset = match op {
            Bilinear::Cup => 0,
            Bilinear::Bracket | Bilinear::Insert => -1,
        };
        Ok(Self {
            arity: self
                .arity
                .checked_add(other.arity)
                .and_then(|a| a.checked_add(offset))
                .ok_or_else(grading_range)?,
            shift: self
                .shift
                .checked_add(other.shift)
                .ok_or_else(grading_range)?,
            ..self.clone()
        })
    }

    pub(crate) fn possible(&self) -> bool {
        if self.arity < 0 {
            return false;
        }
        let top = self.length.saturating_mul(2);
        let smallest = self.arity;
        let largest = self.arity.saturating_mul(top);
        self.shift >= -largest && self.shift <= top - smallest
    }
}

fn grading_range() -> Error {
    Error::plain("Verbitsky Hochschild grading is outside the signed integer range")
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
struct Key {
    word: Vec<usize>,
    output: usize,
}

#[derive(Default)]
pub(crate) struct Cochain {
    entries: BTreeMap<Key, Rational>,
    cells: usize,
}

impl Cochain {
    fn insert(&mut self, model: &mut Model, key: Key, coefficient: Rational) -> Result<()> {
        use std::collections::btree_map::Entry;
        if coefficient.is_zero() {
            return Ok(());
        }
        model.spend(1 + key.word.len())?;
        let cells = key.word.len() + 2 * LIMBS + 10;
        match self.entries.entry(key) {
            Entry::Vacant(entry) => {
                storage(self.cells + cells)?;
                entry.insert(coefficient);
                self.cells += cells;
            }
            Entry::Occupied(mut entry) => {
                let sum = entry.get().add(model, coefficient)?;
                if sum.is_zero() {
                    entry.remove();
                    self.cells -= cells;
                } else {
                    entry.insert(sum);
                }
            }
        }
        Ok(())
    }

    fn empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(self.cells + 1)?;
        Ok(Value::VerbitskyHochschild(Rc::new(self)))
    }

    pub(crate) fn text(&self, model: &mut Model) -> Result<String> {
        let mut report = "[".to_owned();
        for (i, (key, coefficient)) in self.entries.iter().enumerate() {
            if i > 0 {
                report.push_str(", ");
            }
            write!(
                report,
                "({:?}, {}, {})",
                key.word,
                key.output,
                coefficient.text(model)?
            )
            .unwrap();
            if report.len() >= MAX_SHOWN_BYTES {
                return Err(Error::plain(
                    "Verbitsky Hochschild report exceeds the output limit",
                ));
            }
        }
        report.push(']');
        Ok(report)
    }

    pub(crate) fn add(&self, model: &mut Model, other: &Self) -> Result<Self> {
        let mut result = Self::default();
        for (key, &coefficient) in self.entries.iter().chain(&other.entries) {
            result.insert(model, key.clone(), coefficient)?;
        }
        Ok(result)
    }

    pub(crate) fn scale(&self, model: &mut Model, scalar: Entry) -> Result<Self> {
        let scalar = Rational::input(model, scalar)?;
        let mut result = Self::default();
        if !scalar.is_zero() {
            for (key, &coefficient) in &self.entries {
                let value = coefficient.multiply(model, scalar)?;
                result.insert(model, key.clone(), value)?;
            }
        }
        Ok(result)
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        for ((key, value), (other_key, other_value)) in self.entries.iter().zip(&other.entries) {
            model.spend(1 + key.word.len().min(other_key.word.len()))?;
            let order = key.cmp(other_key);
            if !order.is_eq() {
                return Ok(order);
            }
            let order = value.compare(*other_value);
            if !order.is_eq() {
                return Ok(order);
            }
        }
        Ok(self.entries.len().cmp(&other.entries.len()))
    }
}

impl Form {
    pub(crate) fn hochschild_validate(&self, grade: &Grade) -> Result<()> {
        self.validate_grade(&grade.ring(0))
    }

    pub(crate) fn hochschild_dimensions(
        &self,
        model: &mut Model,
        grade: &Grade,
    ) -> Result<[usize; 4]> {
        self.hochschild_validate(grade)?;
        if !grade.possible() {
            return Ok([0; 4]);
        }
        let algebra = self.hochschild_algebra(model, grade)?;
        let complex = Complex::new(model, &algebra, grade)?;
        Ok(complex.dimensions())
    }

    pub(crate) fn hochschild_term(
        &self,
        model: &mut Model,
        grade: &Grade,
        word: Vec<usize>,
        output: usize,
        coefficient: Entry,
    ) -> Result<Cochain> {
        let algebra = self.hochschild_algebra(model, grade)?;
        let key = Key { word, output };
        algebra.validate_key(grade, &key)?;
        let mut cochain = Cochain::default();
        let coefficient = Rational::input(model, coefficient)?;
        cochain.insert(model, key, coefficient)?;
        Ok(cochain)
    }

    pub(crate) fn hochschild_coefficient(
        &self,
        model: &mut Model,
        grade: &Grade,
        value: &Cochain,
        word: Vec<usize>,
        output: usize,
    ) -> Result<String> {
        let algebra = self.hochschild_algebra(model, grade)?;
        let key = Key { word, output };
        algebra.validate_key(grade, &key)?;
        value
            .entries
            .get(&key)
            .copied()
            .unwrap_or_else(Rational::zero)
            .text(model)
    }

    pub(crate) fn hochschild_basis(
        &self,
        model: &mut Model,
        grade: &Grade,
        index: usize,
    ) -> Result<Cochain> {
        self.hochschild_validate(grade)?;
        if !grade.possible() {
            return Err(basis_error());
        }
        let algebra = self.hochschild_algebra(model, grade)?;
        Complex::new(model, &algebra, grade)?.basis(model, index)
    }

    pub(crate) fn hochschild_class(
        &self,
        model: &mut Model,
        grade: &Grade,
        value: &Cochain,
    ) -> Result<Cochain> {
        if value.empty() {
            return Ok(Cochain::default());
        }
        let algebra = self.hochschild_algebra(model, grade)?;
        Complex::new(model, &algebra, grade)?.project(model, value)
    }
}

fn basis_error() -> Error {
    Error::plain("Verbitsky Hochschild cohomology basis index out of range")
}
