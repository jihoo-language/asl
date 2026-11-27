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
use super::Grade as RingGrade;

mod complex;
pub(crate) use complex::Homology;

pub(crate) type Cache = RefCell<BTreeMap<(i128, i128, i128), Rc<Homology>>>;

pub(crate) enum Command {
    Betti(Label, i128, i128, i128),
    Table(Label, i128),
    Relation(Box<ParsedExpr>),
}

pub(crate) enum Expression {
    Term(Vec<usize>, Box<ParsedExpr>),
    Differential(Box<ParsedExpr>),
    Cast(bool, Box<ParsedExpr>),
    Basis(Label, i128, i128, i128, usize),
    Product(bool, Box<[ParsedExpr; 2]>),
    Scale(bool, Entry, Box<ParsedExpr>),
}

#[derive(Clone, Hash, Eq, PartialEq)]
pub(crate) struct Grade {
    pub(crate) basis: Rc<str>,
    pub(crate) length: i128,
    pub(crate) homological: i128,
    pub(crate) weight: i128,
}

impl Grade {
    fn ring_grade(&self, weight: u128) -> Result<RingGrade> {
        let degree = weight
            .checked_mul(2)
            .and_then(|d| i128::try_from(d).ok())
            .ok_or_else(grading_range)?;
        Ok(RingGrade {
            basis: self.basis.clone(),
            length: self.length,
            degree,
        })
    }

    pub(crate) fn differential(&self) -> Result<Self> {
        let homological = self.homological.checked_sub(1).ok_or_else(grading_range)?;
        Ok(Self {
            homological,
            ..self.clone()
        })
    }

    pub(crate) fn product(&self, right: &Self) -> Result<Self> {
        if self.basis != right.basis || self.length != right.length {
            return Err(Error::plain(
                "Koszul product requires the same marking and Hilbert length",
            ));
        }
        let homological = self
            .homological
            .checked_add(right.homological)
            .ok_or_else(grading_range)?;
        let weight = self
            .weight
            .checked_add(right.weight)
            .ok_or_else(grading_range)?;
        Ok(Self {
            homological,
            weight,
            ..self.clone()
        })
    }
}

pub(crate) fn grading_range() -> Error {
    Error::plain("Verbitsky Koszul grade exceeds the supported integer range")
}

fn indexing() -> Error {
    Error::plain("inconsistent Verbitsky Koszul basis index")
}

#[derive(Default)]
pub(crate) struct Chain {
    terms: BTreeMap<u32, Tensor>,
    cells: usize,
}

impl Chain {
    pub(crate) fn term(model: &mut Model, mask: u32, coefficient: &Tensor) -> Result<Self> {
        let mut output = Self::default();
        output.insert(model, mask, coefficient, false)?;
        Ok(output)
    }

    fn insert(&mut self, model: &mut Model, mask: u32, source: &Tensor, minus: bool) -> Result<()> {
        if source.terms.is_empty() {
            return Ok(());
        }
        let before = self.terms.get(&mask).map_or(0, |value| value.cells + 8);
        let remaining = self.cells - before;
        let target = self.terms.entry(mask).or_default();
        for (powers, &coefficient) in &source.terms {
            target.insert(
                model,
                powers.clone(),
                if minus {
                    coefficient.negate()
                } else {
                    coefficient
                },
            )?;
            storage(remaining + target.cells + 8)?;
        }
        self.cells = if target.terms.is_empty() {
            self.terms.remove(&mask);
            remaining
        } else {
            remaining + target.cells + 8
        };
        storage(self.cells)?;
        Ok(())
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(self.cells + 1)?;
        Ok(Value::VerbitskyKoszul(Rc::new(self)))
    }

    pub(crate) fn add(&self, model: &mut Model, right: &Self) -> Result<Self> {
        let mut sum = Self::default();
        for (&mask, value) in self.terms.iter().chain(&right.terms) {
            sum.insert(model, mask, value, false)?;
        }
        Ok(sum)
    }

    pub(crate) fn scale(&self, model: &mut Model, scalar: Entry) -> Result<Self> {
        let scalar = Rational::input(model, scalar)?;
        let mut output = Self::default();
        if !scalar.is_zero() {
            for (&mask, value) in &self.terms {
                let mut tensor = Tensor::default();
                for (powers, &coefficient) in &value.terms {
                    let value = coefficient.multiply(model, scalar)?;
                    tensor.insert(model, powers.clone(), value)?;
                }
                output.insert(model, mask, &tensor, false)?;
            }
        }
        Ok(output)
    }

    pub(crate) fn compare(&self, model: &mut Model, right: &Self) -> Result<Ordering> {
        for ((a, left), (b, right)) in self.terms.iter().zip(&right.terms) {
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
        Ok(self.terms.len().cmp(&right.terms.len()))
    }

    pub(crate) fn text(&self, model: &mut Model, form: &Form) -> Result<String> {
        let mut output = "[".to_owned();
        let mut separator = "";
        for (&mask, value) in &self.terms {
            let indices: Vec<_> = bits(mask).collect();
            write!(
                output,
                "{separator}({indices:?}, {})",
                value.text(model, form)?
            )
            .unwrap();
            if output.len() >= MAX_SHOWN_BYTES {
                return Err(Error::plain(
                    "Verbitsky Koszul representative exceeds the report limit",
                ));
            }
            separator = ", ";
        }
        output.push(']');
        Ok(output)
    }
}

impl Form {
    fn variables(&self) -> usize {
        self.gram.size() + self.field.degree()
    }

    pub(crate) fn supports_koszul(&self, grade: &Grade) -> bool {
        if grade.homological < 0
            || grade.homological > self.variables() as i128
            || grade.weight < grade.homological
        {
            return false;
        }
        let p = grade.homological as usize;
        let excess = (grade.weight - grade.homological) as u128;
        let least = (p.saturating_sub(self.gram.size()) as u128)
            .max(excess.saturating_sub(2 * grade.length as u128));
        let greatest = (p.min(self.field.degree()) as u128).min(excess);
        least <= greatest && (self.gram.size() != 0 || (excess - least).is_multiple_of(2))
    }

    pub(crate) fn validate_koszul(&self, grade: &Grade) -> Result<()> {
        self.validate_grade(&grade.ring_grade(0)?)
    }

    fn mask_weight(&self, mask: u32) -> u128 {
        u128::from(mask.count_ones() + (mask >> self.gram.size()).count_ones())
    }

    pub(crate) fn term_grade(&self, source: &RingGrade, indices: &[usize]) -> Result<(u32, Grade)> {
        if indices.last().is_some_and(|&i| i >= self.variables())
            || !indices.windows(2).all(|pair| pair[0] < pair[1])
        {
            return Err(Error::plain(
                "Koszul generator indices must be strictly increasing and in range",
            ));
        }
        let mask = indices.iter().fold(0, |mask, &i| mask | 1u32 << i);
        let weight = source.degree / 2 + self.mask_weight(mask) as i128;
        Ok((
            mask,
            Grade {
                basis: source.basis.clone(),
                length: source.length,
                homological: indices.len() as i128,
                weight,
            },
        ))
    }

    fn coefficient_grade(&self, chain: &Grade, mask: u32) -> Result<Option<RingGrade>> {
        let Ok(weight) = u128::try_from(chain.weight) else {
            return Ok(None);
        };
        let Some(weight) = weight
            .checked_sub(self.mask_weight(mask))
            .filter(|&w| w <= 2 * chain.length as u128)
        else {
            return Ok(None);
        };
        chain.ring_grade(weight).map(Some)
    }

    pub(crate) fn koszul_diff(
        &self,
        model: &mut Model,
        grade: &Grade,
        chain: &Chain,
    ) -> Result<Chain> {
        let mut result = Chain::default();
        for (&mask, value) in &chain.terms {
            for (position, i) in bits(mask).enumerate() {
                let next = mask ^ (1u32 << i);
                let Some(coefficient_grade) = self.coefficient_grade(grade, next)? else {
                    continue;
                };
                let raw = value.cup(model, &Tensor::variable(i as u8))?;
                let value = self.project(model, &coefficient_grade, &raw)?;
                result.insert(model, next, &value, !position.is_multiple_of(2))?;
            }
        }
        Ok(result)
    }

    pub(crate) fn koszul_product(
        &self,
        model: &mut Model,
        grade: &Grade,
        left: &Chain,
        right: &Chain,
    ) -> Result<Chain> {
        let mut product = Chain::default();
        model
            .spend(left.terms.len().saturating_mul(right.terms.len()))
            .map_err(|_| work_limit())?;
        for (&a, x) in &left.terms {
            for (&b, y) in &right.terms {
                if a & b != 0 {
                    continue;
                }
                let mask = a | b;
                let Some(coefficient_grade) = self.coefficient_grade(grade, mask)? else {
                    continue;
                };
                // Every generator of the right factor that precedes a
                // generator of the left factor contributes one exterior swap.
                let swaps: u32 = bits(a).map(|i| (b & ((1u32 << i) - 1)).count_ones()).sum();
                let raw = x.cup(model, y)?;
                let value = self.project(model, &coefficient_grade, &raw)?;
                product.insert(model, mask, &value, !swaps.is_multiple_of(2))?;
            }
        }
        Ok(product)
    }

    pub(crate) fn minimal_relation(
        &self,
        model: &mut Model,
        grade: &Grade,
        cycle: &Chain,
    ) -> Result<String> {
        self.lift_minimal_relation(model, grade, cycle)?
            .text(model, self)
    }

    pub(super) fn lift_minimal_relation(
        &self,
        model: &mut Model,
        grade: &Grade,
        cycle: &Chain,
    ) -> Result<Tensor> {
        if grade.homological != 1 {
            return Err(Error::plain(
                "Verbitsky minimal relations require homological degree one",
            ));
        }
        let mut polynomial = Tensor::default();
        for (&mask, value) in &cycle.terms {
            let generator = Power {
                index: mask.trailing_zeros() as u8,
                exponent: 1,
            };
            for (powers, &coefficient) in &value.terms {
                let term = multiply_powers(powers, &[generator])?;
                polynomial.insert(model, term, coefficient)?;
            }
        }
        let quotient_grade = grade.ring_grade(grade.weight.max(0) as u128)?;
        if !self
            .project(model, &quotient_grade, &polynomial)?
            .terms
            .is_empty()
        {
            return Err(Error::plain(
                "Verbitsky minimal relation has nonzero image in the quotient",
            ));
        }
        Ok(polynomial)
    }

    pub(crate) fn tor_table(&self, model: &mut Model, grade: &Grade) -> Result<String> {
        let variables = self.variables();
        let max_weight =
            (2 * grade.length as u128).saturating_add((variables + self.field.degree()) as u128);
        let max_weight = usize::try_from(max_weight).map_err(|_| work_limit())?;
        storage(
            (variables + 3)
                .saturating_mul(max_weight.saturating_add(1))
                .saturating_mul(2),
        )?;
        model
            .spend((variables + 1) * (max_weight + 1))
            .map_err(|_| work_limit())?;
        let mut table = vec![vec![0usize; max_weight + 1]; variables + 1];
        for (p, row) in table.iter_mut().enumerate() {
            for (w, value) in row.iter_mut().enumerate() {
                *value = self.tor_betti(
                    model,
                    &Grade {
                        homological: p as i128,
                        weight: w as i128,
                        ..grade.clone()
                    },
                )?;
            }
        }
        let mut coefficients = vec![0i128; max_weight + 1];
        for (w, x) in coefficients
            .iter_mut()
            .take(2 * grade.length as usize + 1)
            .enumerate()
        {
            *x = self.rank(&grade.ring_grade(w as u128)?)? as i128;
        }
        for i in 0..variables {
            let a = 1 + usize::from(i >= self.gram.size());
            for w in (a..=max_weight).rev() {
                coefficients[w] -= coefficients[w - a];
            }
        }
        let mut euler = vec![0i128; max_weight + 1];
        let mut result = "betti=[".to_owned();
        let mut separator = "";
        for (p, row) in table.iter().enumerate() {
            for (w, &betti) in row.iter().enumerate() {
                if betti != table[variables - p][max_weight - w] {
                    return Err(Error::plain(
                        "Verbitsky Koszul Betti numbers violate Gorenstein symmetry",
                    ));
                }
                let sign = if p.is_multiple_of(2) { 1 } else { -1 };
                euler[w] += sign * betti as i128;
                if betti != 0 {
                    write!(result, "{separator}({p}, {w}, {betti})").unwrap();
                    if result.len() >= MAX_SHOWN_BYTES {
                        return Err(Error::plain(
                            "Verbitsky Koszul table exceeds the report limit",
                        ));
                    }
                    separator = ", ";
                }
            }
        }
        if euler != coefficients {
            return Err(Error::plain(
                "Verbitsky Koszul Betti numbers disagree with the weighted Hilbert numerator",
            ));
        }
        let numerator = coefficients
            .into_iter()
            .enumerate()
            .filter(|&(_, x)| x != 0)
            .collect::<Vec<_>>();
        write!(result, "] numerator={numerator:?}").unwrap();
        if result.len() >= MAX_SHOWN_BYTES {
            return Err(Error::plain(
                "Verbitsky Koszul table exceeds the report limit",
            ));
        }
        Ok(result)
    }
}

fn bits(mut mask: u32) -> impl Iterator<Item = usize> {
    std::iter::from_fn(move || {
        if mask == 0 {
            return None;
        }
        let index = mask.trailing_zeros() as usize;
        mask &= mask - 1;
        Some(index)
    })
}
