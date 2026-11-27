use std::cmp::Ordering;

use super::{
    linear, storage, work_limit, Algebra, BTreeMap, Cochain, Entry, Error, Form, Label, Model,
    ParsedExpr, Rational, Rc, Result, Value, Write, LIMBS, MAX_SHOWN_BYTES,
};

mod complex;
mod maps;
mod mixed;
use complex::Complex;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) enum Flavor {
    Hochschild,
    Cyclic,
    Negative,
    Periodic,
}

impl Flavor {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Hochschild => "hh",
            Self::Cyclic => "hc",
            Self::Negative => "hn",
            Self::Periodic => "hp",
        }
    }

    fn allows(self, exponent: i128) -> bool {
        match self {
            Self::Hochschild => exponent == 0,
            Self::Cyclic => exponent <= 0,
            Self::Negative => exponent >= 0,
            Self::Periodic => true,
        }
    }
}

pub(crate) struct SyntaxGrade(
    pub(crate) Label,
    pub(crate) i128,
    pub(crate) Flavor,
    pub(crate) i128,
    pub(crate) i128,
);

impl From<SyntaxGrade> for Grade {
    fn from(SyntaxGrade(name, length, flavor, degree, weight): SyntaxGrade) -> Self {
        Self {
            basis: name.text.into(),
            length,
            flavor,
            degree,
            weight,
        }
    }
}

pub(crate) enum Command {
    Dims(SyntaxGrade),
    Sbi(Label, i128, i128, i128),
    Pair(Box<[ParsedExpr; 2]>),
}

pub(crate) enum Expression {
    Term(SyntaxGrade, i128, Vec<usize>, Entry),
    Basis(SyntaxGrade, usize),
    Scale(Entry, Box<ParsedExpr>),
    Apply(Operation, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Diff,
    Connes,
    Class,
    Representative,
    Primitive,
    S,
    I,
    Boundary,
    U,
    Eval,
    Periodic,
    Lift,
}

impl Operation {
    pub(crate) fn named(name: &str) -> Option<Self> {
        if !name.starts_with("rvVerbitsky") {
            return None;
        }
        [
            Self::Diff,
            Self::Connes,
            Self::Class,
            Self::Representative,
            Self::Primitive,
            Self::S,
            Self::I,
            Self::Boundary,
            Self::U,
            Self::Eval,
            Self::Periodic,
            Self::Lift,
        ]
        .into_iter()
        .find(|op| op.name() == name)
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Diff => "rvVerbitskyCyclicDiff",
            Self::Connes => "rvVerbitskyChainConnes",
            Self::Class => "rvVerbitskyCyclicClass",
            Self::Representative => "rvVerbitskyCyclicRepresentative",
            Self::Primitive => "rvVerbitskyCyclicPrimitive",
            Self::S => "rvVerbitskyCyclicS",
            Self::I => "rvVerbitskyCyclicI",
            Self::Boundary => "rvVerbitskyCyclicB",
            Self::U => "rvVerbitskyNegativeU",
            Self::Eval => "rvVerbitskyNegativeEval",
            Self::Periodic => "rvVerbitskyNegativePeriodic",
            Self::Lift => "rvVerbitskyNegativeLift",
        }
    }
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) basis: Rc<str>,
    pub(crate) length: i128,
    pub(crate) flavor: Flavor,
    pub(crate) degree: i128,
    pub(crate) weight: i128,
}

impl Grade {
    fn algebra_grade(&self) -> super::Grade {
        super::Grade {
            basis: self.basis.clone(),
            length: self.length,
            arity: 0,
            shift: 0,
        }
    }

    pub(crate) fn possible(&self) -> bool {
        self.weight >= 0
            && match self.flavor {
                Flavor::Hochschild => (0..=self.weight).contains(&self.degree),
                Flavor::Cyclic => self.degree >= 0,
                Flavor::Negative => self.degree <= self.weight,
                Flavor::Periodic => true,
            }
            && (self.weight > 0 || self.degree.rem_euclid(2) == 0)
    }

    fn with_offset(&self, offset: i128) -> Result<Self> {
        Ok(Self {
            degree: self.degree.checked_add(offset).ok_or_else(grading_error)?,
            ..self.clone()
        })
    }

    fn with_flavor(&self, flavor: Flavor) -> Self {
        Self {
            flavor,
            ..self.clone()
        }
    }

    fn power(&self, arity: usize) -> Option<i128> {
        let arity = arity as i128;
        if arity.rem_euclid(2) != self.degree.rem_euclid(2) {
            return None;
        }
        arity.div_euclid(2).checked_sub(self.degree.div_euclid(2))
    }
}

fn grading_error() -> Error {
    Error::plain("Verbitsky cyclic grading is outside the signed integer range")
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
struct Key {
    power: i128,
    word: Vec<usize>,
}

#[derive(Default)]
pub(crate) struct Chain {
    entries: BTreeMap<Key, Rational>,
    cells: usize,
}

impl Chain {
    fn insert(&mut self, model: &mut Model, key: Key, coefficient: Rational) -> Result<()> {
        use std::collections::btree_map::Entry;
        if coefficient.is_zero() {
            return Ok(());
        }
        model.spend(key.word.len() + 1)?;
        let cells = key.word.len() + 2 * LIMBS + 13;
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
        Ok(Value::VerbitskyCyclic(Rc::new(self)))
    }

    pub(crate) fn text(&self, model: &mut Model) -> Result<String> {
        let mut report = "[".to_owned();
        for (i, (key, coefficient)) in self.entries.iter().enumerate() {
            if i > 0 {
                report.push_str(", ");
            }
            write!(
                report,
                "({}, {:?}, {})",
                key.power,
                key.word,
                coefficient.text(model)?
            )
            .unwrap();
            if report.len() >= MAX_SHOWN_BYTES {
                return Err(Error::plain(
                    "Verbitsky cyclic report exceeds the output limit",
                ));
            }
        }
        report.push(']');
        Ok(report)
    }

    pub(crate) fn add(&self, model: &mut Model, other: &Self) -> Result<Self> {
        let mut result = Self::default();
        for (key, &value) in self.entries.iter().chain(&other.entries) {
            result.insert(model, key.clone(), value)?;
        }
        Ok(result)
    }

    pub(crate) fn scale(&self, model: &mut Model, scalar: Entry) -> Result<Self> {
        let scalar = Rational::input(model, scalar)?;
        let mut result = Self::default();
        if !scalar.is_zero() {
            for (key, &value) in &self.entries {
                let value = value.multiply(model, scalar)?;
                result.insert(model, key.clone(), value)?;
            }
        }
        Ok(result)
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        for ((key, value), (k, v)) in self.entries.iter().zip(&other.entries) {
            model.spend(1 + key.word.len().min(k.word.len()))?;
            let order = key.cmp(k).then_with(|| value.compare(*v));
            if !order.is_eq() {
                return Ok(order);
            }
        }
        Ok(self.entries.len().cmp(&other.entries.len()))
    }
}

impl Form {
    pub(crate) fn cyclic_validate(&self, grade: &Grade) -> Result<()> {
        self.hochschild_validate(&grade.algebra_grade())
    }

    fn cyclic_algebra(&self, model: &mut Model, grade: &Grade) -> Result<Rc<Algebra>> {
        self.hochschild_algebra(model, &grade.algebra_grade())
    }

    pub(crate) fn cyclic_dimensions(&self, model: &mut Model, grade: &Grade) -> Result<[usize; 4]> {
        self.cyclic_validate(grade)?;
        if !grade.possible() {
            return Ok([0; 4]);
        }
        let algebra = self.cyclic_algebra(model, grade)?;
        Ok(Complex::new(model, &algebra, grade)?.dimensions())
    }

    pub(crate) fn cyclic_term(
        &self,
        model: &mut Model,
        grade: &Grade,
        power: i128,
        word: Vec<usize>,
        coefficient: Entry,
    ) -> Result<Chain> {
        let Some(arity) = word.len().checked_sub(1) else {
            return Err(Error::plain(
                "Verbitsky chain requires a leading ring factor",
            ));
        };
        if word[1..].contains(&0) {
            return Err(Error::plain(
                "normalized chain bars must have positive weight",
            ));
        }
        if !grade.flavor.allows(power) || grade.power(arity) != Some(power) {
            return Err(Error::plain(
                "Verbitsky cyclic term has the wrong homological degree or u exponent",
            ));
        }
        let algebra = self.cyclic_algebra(model, grade)?;
        let weight = word.iter().try_fold(0i128, |sum, &i| {
            algebra
                .weights
                .get(i)
                .map(|&w| sum + w as i128)
                .ok_or_else(|| {
                    Error::plain("Verbitsky cyclic term index is outside its ring basis")
                })
        })?;
        if weight != grade.weight {
            return Err(Error::plain(
                "Verbitsky cyclic term has the wrong internal weight",
            ));
        }
        let mut result = Chain::default();
        let coefficient = Rational::input(model, coefficient)?;
        result.insert(model, Key { power, word }, coefficient)?;
        Ok(result)
    }

    pub(crate) fn cyclic_basis(
        &self,
        model: &mut Model,
        grade: &Grade,
        index: usize,
    ) -> Result<Chain> {
        if !grade.possible() {
            return Err(basis_error());
        }
        let algebra = self.cyclic_algebra(model, grade)?;
        Complex::new(model, &algebra, grade)?.basis(model, index)
    }

    pub(crate) fn cyclic_pair(
        &self,
        model: &mut Model,
        grade: &Grade,
        cochain: &Cochain,
        chain: &Chain,
    ) -> Result<String> {
        let algebra = self.cyclic_algebra(model, grade)?;
        let mut sum = Rational::zero();
        for (key, &coefficient) in &chain.entries {
            let start = super::Key {
                word: key.word[1..].to_vec(),
                output: 0,
            };
            let mut trace = Rational::zero();
            for (entry, &value) in cochain
                .entries
                .range(start.clone()..)
                .take_while(|(entry, _)| entry.word == start.word)
            {
                let product = algebra.product(key.word[0], entry.output);
                if let Some(&(_, top)) = product
                    .iter()
                    .find(|&&(index, _)| index == algebra.top_index())
                {
                    let term = top.multiply(model, value)?;
                    trace = trace.add(model, term)?;
                }
            }
            let term = trace.multiply(model, coefficient)?;
            sum = sum.add(model, term)?;
        }
        sum.multiply(model, algebra.trace_factor())?.text(model)
    }
}

fn basis_error() -> Error {
    Error::plain("Verbitsky cyclic homology basis index out of range")
}
