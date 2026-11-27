use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value, MAX_SHOWN_BYTES};
use crate::syntax::{Label, ParsedExpr};

use super::super::rational::Rational;
use super::super::{multiply_powers, storage, work_limit, Entry, Form, Power, Tensor};
use super::{koszul, Grade as RingGrade, Monomial};

mod builder;
mod pfaffian;
pub(crate) use builder::Resolution;

pub(crate) type Cache = RefCell<BTreeMap<i128, Rc<Resolution>>>;

pub(crate) enum Command {
    Shifts(Label, i128),
    Entry(Label, i128, usize, usize, usize),
    Pfaffian(Label, i128),
}

pub(crate) enum Expression {
    Basis(Label, i128, usize, usize, Vec<u128>),
    Times(Vec<u128>, Box<ParsedExpr>),
    Differential(Box<ParsedExpr>),
    Augment(Box<ParsedExpr>),
    Scale(Entry, Box<ParsedExpr>),
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) basis: Rc<str>,
    pub(crate) length: i128,
    pub(crate) homological: i128,
    pub(crate) weight: i128,
}

impl Grade {
    pub(crate) fn differential(&self) -> Result<Self> {
        Ok(Self {
            homological: self
                .homological
                .checked_sub(1)
                .ok_or_else(koszul::grading_range)?,
            ..self.clone()
        })
    }

    pub(crate) fn ring_grade(&self) -> Result<RingGrade> {
        if self.homological != 0 {
            return Err(Error::plain(
                "Verbitsky augmentation requires a homological degree zero vector",
            ));
        }
        let degree = self
            .weight
            .checked_mul(2)
            .filter(|&d| d >= 0)
            .ok_or_else(koszul::grading_range)?;
        Ok(RingGrade {
            basis: self.basis.clone(),
            length: self.length,
            degree,
        })
    }

    fn koszul_grade(&self) -> koszul::Grade {
        koszul::Grade {
            basis: self.basis.clone(),
            length: self.length,
            homological: self.homological,
            weight: self.weight,
        }
    }
}

#[derive(Default)]
pub(crate) struct Vector {
    entries: BTreeMap<usize, Tensor>,
    cells: usize,
}

impl Vector {
    fn insert(
        &mut self,
        model: &mut Model,
        index: usize,
        powers: Monomial,
        coefficient: Rational,
    ) -> Result<()> {
        if coefficient.is_zero() {
            return Ok(());
        }
        let prior = self.entries.get(&index).map_or(0, |term| term.cells + 8);
        let remaining = self.cells - prior;
        let target = self.entries.entry(index).or_default();
        target.insert(model, powers, coefficient)?;
        self.cells = if target.terms.is_empty() {
            self.entries.remove(&index);
            remaining
        } else {
            remaining + target.cells + 8
        };
        storage(self.cells)
    }

    fn add_tensor(&mut self, model: &mut Model, index: usize, value: &Tensor) -> Result<()> {
        for (powers, &coefficient) in &value.terms {
            self.insert(model, index, powers.clone(), coefficient)?;
        }
        Ok(())
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(self.cells + 1)?;
        Ok(Value::VerbitskyFree(Rc::new(self)))
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let mut sum = Self::default();
        for (&i, value) in self.entries.iter().chain(&rhs.entries) {
            sum.add_tensor(model, i, value)?;
        }
        Ok(sum)
    }

    pub(crate) fn scale(&self, model: &mut Model, scalar: Entry) -> Result<Self> {
        let scale = Rational::input(model, scalar)?;
        let mut result = Self::default();
        if !scale.is_zero() {
            for (&index, polynomial) in &self.entries {
                for (powers, &coefficient) in &polynomial.terms {
                    let value = coefficient.multiply(model, scale)?;
                    result.insert(model, index, powers.clone(), value)?;
                }
            }
        }
        Ok(result)
    }

    pub(crate) fn times(&self, model: &mut Model, form: &Form, powers: &[u128]) -> Result<Self> {
        let (monomial, _) = form.free_monomial(powers)?;
        let mut result = Self::default();
        for (&index, polynomial) in &self.entries {
            for (powers, &coefficient) in &polynomial.terms {
                result.insert(
                    model,
                    index,
                    multiply_powers(powers, &monomial)?,
                    coefficient,
                )?;
            }
        }
        Ok(result)
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        for ((a, left), (b, right)) in self.entries.iter().zip(&rhs.entries) {
            model.spend(1)?;
            let order = a.cmp(b);
            if !order.is_eq() {
                return Ok(order);
            }
            let order = left.compare(model, right)?;
            if !order.is_eq() {
                return Ok(order);
            }
        }
        Ok(self.entries.len().cmp(&rhs.entries.len()))
    }

    pub(crate) fn text(&self, model: &mut Model, form: &Form) -> Result<String> {
        let mut text = "[".to_owned();
        let mut separator = "";
        for (&index, polynomial) in &self.entries {
            write!(
                text,
                "{separator}({index}, {})",
                polynomial.text(model, form)?
            )
            .unwrap();
            check_report(&text)?;
            separator = ", ";
        }
        text.push(']');
        Ok(text)
    }
}

pub(super) fn check_report(text: &str) -> Result<()> {
    if text.len() >= MAX_SHOWN_BYTES {
        Err(Error::plain(
            "Verbitsky free resolution exceeds the report limit",
        ))
    } else {
        Ok(())
    }
}

impl Form {
    fn variable_count(&self) -> usize {
        self.gram.size() + self.field.degree()
    }

    pub(crate) fn validate_free(&self, grade: &Grade) -> Result<()> {
        self.validate_koszul(&grade.koszul_grade())
    }

    pub(crate) fn supports_free(&self, grade: &Grade) -> bool {
        grade.weight >= 0
            && grade.homological >= 0
            && grade.homological <= self.variable_count() as i128
    }

    fn polynomial_weight(&self, powers: &[Power]) -> Result<u128> {
        powers
            .iter()
            .try_fold(0u128, |weight, power| {
                let degree = 1 + u128::from(usize::from(power.index) >= self.gram.size());
                weight.checked_add(degree.checked_mul(power.exponent)?)
            })
            .ok_or_else(koszul::grading_range)
    }

    fn free_monomial(&self, dense: &[u128]) -> Result<(Monomial, u128)> {
        if dense.len() > self.variable_count() {
            return Err(Error::plain(
                "Verbitsky free polynomial has too many variable exponents",
            ));
        }
        let powers = dense
            .iter()
            .enumerate()
            .filter_map(|(i, &exponent)| {
                (exponent != 0).then_some(Power {
                    index: i as u8,
                    exponent,
                })
            })
            .collect::<Vec<_>>();
        let weight = self.polynomial_weight(&powers)?;
        Ok((powers, weight))
    }

    pub(crate) fn free_times_grade(&self, grade: &Grade, dense: &[u128]) -> Result<Grade> {
        let (_, weight) = self.free_monomial(dense)?;
        let weight = i128::try_from(weight)
            .ok()
            .and_then(|w| grade.weight.checked_add(w))
            .ok_or_else(koszul::grading_range)?;
        Ok(Grade {
            weight,
            ..grade.clone()
        })
    }

    pub(crate) fn free_basis(
        &self,
        model: &mut Model,
        grade: &Grade,
        index: usize,
        dense: &[u128],
    ) -> Result<(Grade, Vector)> {
        let (powers, weight) = self.free_monomial(dense)?;
        let p = usize::try_from(grade.homological).map_err(|_| invalid_index())?;
        if p > self.variable_count() || (p == 0 && index != 0) {
            return Err(invalid_index());
        }
        let shift = if p == 0 {
            0
        } else {
            self.resolution(model, grade)?.shift(p, index)?
        };
        let weight = i128::try_from(weight)
            .ok()
            .and_then(|w| w.checked_add(shift as i128))
            .ok_or_else(koszul::grading_range)?;
        let mut vector = Vector::default();
        vector.insert(model, index, powers, Rational::integer(1))?;
        Ok((
            Grade {
                weight,
                ..grade.clone()
            },
            vector,
        ))
    }

    pub(crate) fn free_differential(
        &self,
        model: &mut Model,
        grade: &Grade,
        x: &Vector,
    ) -> Result<Vector> {
        if grade.homological <= 0 || x.entries.is_empty() {
            return Ok(Vector::default());
        }
        self.resolution(model, grade)?
            .differential(model, grade.homological as usize, x)
    }

    pub(crate) fn free_augment(
        &self,
        model: &mut Model,
        grade: &Grade,
        x: &Vector,
    ) -> Result<Tensor> {
        let ring = grade.ring_grade()?;
        if let Some(polynomial) = x.entries.get(&0) {
            self.project(model, &ring, polynomial)
        } else {
            Ok(Tensor::default())
        }
    }

    pub(crate) fn free_shifts(&self, model: &mut Model, grade: &Grade) -> Result<String> {
        let resolution = self.resolution(model, grade)?;
        let shifts = resolution.shifts();
        let ranks: Vec<_> = shifts.iter().map(|degree| degree.len()).collect();
        let text = format!("ranks={ranks:?} shifts={shifts:?}");
        check_report(&text)?;
        Ok(text)
    }

    pub(crate) fn free_entry(
        &self,
        model: &mut Model,
        grade: &Grade,
        p: usize,
        i: usize,
        j: usize,
    ) -> Result<String> {
        self.resolution(model, grade)?
            .entry(p, i, j)?
            .unwrap_or(&Tensor::default())
            .text(model, self)
    }

    pub(crate) fn free_pfaffian(&self, model: &mut Model, grade: &Grade) -> Result<String> {
        if self.variable_count() != 3 {
            return Err(Error::plain(
                "Buchsbaum–Eisenbud presentation requires codimension three",
            ));
        }
        let resolution = self.resolution(model, grade)?;
        pfaffian::presentation(model, self, &resolution)
    }
}

fn invalid_index() -> Error {
    Error::plain("Verbitsky free module index is outside the minimal resolution")
}
