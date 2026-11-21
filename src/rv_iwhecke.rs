use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_iwahori;
use crate::rv_iwkl::{self, Window};

const MAX_CELLS: usize = 300_000;

pub(crate) type Monomial = (usize, i128);
pub(crate) type Polynomial = Rc<[Monomial]>;

#[derive(PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Term {
    pub(crate) window: Window,
    pub(crate) coefficient: Polynomial,
}

pub(crate) type Support = Rc<[Term]>;

fn overflow() -> Error {
    Error::plain("affine Hecke arithmetic is out of range")
}

fn storage_limit() -> Error {
    Error::plain("affine Hecke calculation exceeds the storage limit")
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(storage_limit());
    }
    Ok(())
}

fn window_cells(rank: usize) -> Result<usize> {
    rank.checked_mul(2)
        .and_then(|cells| cells.checked_add(1))
        .ok_or_else(storage_limit)
}

pub(crate) fn persist(model: &mut Model, support: Vec<Term>) -> Result<Value> {
    let mut cells = 0;
    for term in &support {
        cells += window_cells(term.window.len())? + 3 * term.coefficient.len();
        check_storage(cells)?;
    }
    model.allocate(cells)?;
    Ok(Value::RVIwHecke(support.into()))
}

pub(crate) fn zero() -> Value {
    Value::RVIwHecke([].into())
}

pub(crate) fn standard(model: &mut Model, window: Vec<i128>) -> Result<Value> {
    check_storage(
        window_cells(window.len())?
            .checked_add(3)
            .ok_or_else(storage_limit)?,
    )?;
    rv_iwahori::validate_window(model, &window)?;
    persist(
        model,
        vec![Term {
            window: window.into(),
            coefficient: [(0, 1)].into(),
        }],
    )
}

pub(crate) fn scalar(model: &mut Model, rank: usize, degree: usize) -> Result<Value> {
    check_storage(
        window_cells(rank)?
            .checked_add(3)
            .ok_or_else(storage_limit)?,
    )?;
    model.spend(rank)?;
    let window: Window = (0..rank).map(|index| index as i128).collect();
    persist(
        model,
        vec![Term {
            window,
            coefficient: [(degree, 1)].into(),
        }],
    )
}

pub(crate) fn canonical(model: &mut Model, window: &[i128]) -> Result<Value> {
    let expansion = rv_iwkl::expansion(model, window)?;
    let support = expansion
        .into_iter()
        .map(|(window, coefficient)| {
            let coefficient = coefficient
                .into_iter()
                .enumerate()
                .filter(|&(_, scalar)| scalar != 0)
                .collect::<Vec<_>>();
            Term {
                window,
                coefficient: coefficient.into(),
            }
        })
        .collect();
    persist(model, support)
}

pub(crate) fn compare(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Ordering> {
    for (left, right) in left.iter().zip(right) {
        model.spend(
            1 + left.window.len().min(right.window.len())
                + left.coefficient.len().min(right.coefficient.len()),
        )?;
        let comparison = left.cmp(right);
        if comparison != Ordering::Equal {
            return Ok(comparison);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

pub(crate) fn add(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Value> {
    let mut sum = Accumulator::default();
    for term in left.iter().chain(right) {
        sum.add(model, &term.window, &term.coefficient, 0, 1)?;
    }
    persist(model, sum.finish())
}

pub(crate) fn neg(model: &mut Model, terms: &[Term]) -> Result<Value> {
    let mut result = Vec::with_capacity(terms.len());
    for term in terms {
        model.spend(1 + term.coefficient.len())?;
        let coefficient = term
            .coefficient
            .iter()
            .map(|&(power, coefficient)| {
                Ok((power, coefficient.checked_neg().ok_or_else(overflow)?))
            })
            .collect::<Result<Vec<_>>>()?;
        result.push(Term {
            window: term.window.clone(),
            coefficient: coefficient.into(),
        });
    }
    persist(model, result)
}

#[derive(Default)]
pub(crate) struct Accumulator {
    terms: BTreeMap<Window, BTreeMap<usize, i128>>,
    cells: usize,
}

impl Accumulator {
    pub(crate) fn add(
        &mut self,
        model: &mut Model,
        window: &Window,
        polynomial: &[Monomial],
        shift: usize,
        scalar: i128,
    ) -> Result<()> {
        let size = window_cells(window.len())?;
        model.spend((1 + self.terms.len().max(1).ilog2() as usize) * size + polynomial.len())?;
        if !self.terms.contains_key(window) {
            self.cells += size;
            check_storage(self.cells)?;
        }
        let output = self.terms.entry(window.clone()).or_default();
        for &(degree, coefficient) in polynomial {
            let degree = degree.checked_add(shift).ok_or_else(overflow)?;
            let coefficient = coefficient.checked_mul(scalar).ok_or_else(overflow)?;
            if let Some(previous) = output.get_mut(&degree) {
                *previous = previous.checked_add(coefficient).ok_or_else(overflow)?;
                if *previous == 0 {
                    output.remove(&degree);
                    self.cells -= 3;
                }
            } else {
                self.cells += 3;
                check_storage(self.cells)?;
                output.insert(degree, coefficient);
            }
        }
        if output.is_empty() {
            self.terms.remove(window);
            self.cells -= size;
        }
        Ok(())
    }

    pub(crate) fn finish(self) -> Vec<Term> {
        self.terms
            .into_iter()
            .map(|(window, coefficient)| Term {
                window,
                coefficient: coefficient.into_iter().collect::<Vec<_>>().into(),
            })
            .collect()
    }
}

fn add_three(first: i128, second: i128, third: i128) -> Result<i128> {
    first
        .checked_add(second)
        .and_then(|sum| sum.checked_add(third))
        .or_else(|| {
            second
                .checked_add(third)
                .and_then(|sum| sum.checked_add(first))
        })
        .ok_or_else(overflow)
}

fn rotate(model: &mut Model, terms: &[Term], component: i128) -> Result<Vec<Term>> {
    let mut output = Vec::with_capacity(terms.len());
    for term in terms {
        let size = term.window.len();
        let rank = size as i128;
        model.spend(size)?;
        let mut window = Vec::with_capacity(size);
        if component == 0 {
            output.push(Term {
                window: term.window.clone(),
                coefficient: term.coefficient.clone(),
            });
            continue;
        }
        let residue = component.rem_euclid(rank);
        for index in 0..size {
            let argument = index as i128 + residue;
            let offset = argument.div_euclid(rank) * rank - residue;
            let entry = term.window[argument.rem_euclid(rank) as usize];
            window.push(add_three(component, offset, entry)?);
        }
        output.push(Term {
            window: window.into(),
            coefficient: term.coefficient.clone(),
        });
    }
    Ok(output)
}

fn simple_product(model: &mut Model, terms: &[Term], simple: usize) -> Result<Vec<Term>> {
    let mut product = Accumulator::default();
    for term in terms {
        let reflected = rv_iwkl::right_simple(model, &term.window, simple)?;
        if rv_iwkl::descent(&term.window, simple)? {
            product.add(model, &reflected, &term.coefficient, 1, 1)?;
            product.add(model, &term.window, &term.coefficient, 1, 1)?;
            product.add(model, &term.window, &term.coefficient, 0, -1)?;
        } else {
            product.add(model, &reflected, &term.coefficient, 0, 1)?;
        }
    }
    Ok(product.finish())
}

pub(crate) fn multiply(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Value> {
    let support = product(model, left, right)?;
    persist(model, support)
}

pub(crate) fn product(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Vec<Term>> {
    if left.is_empty() || right.is_empty() {
        return Ok(Vec::new());
    }
    let mut product = Accumulator::default();
    for term in right {
        let component = rv_iwahori::component(&term.window)?;
        let mut window = rv_iwkl::translate(model, &term.window, component)?;
        let mut word = Vec::new();
        while let Some(simple) = rv_iwkl::wall(model, &window)? {
            check_storage(1 + word.len())?;
            word.push(simple);
            window = rv_iwkl::right_simple(model, &window, simple)?;
        }
        let mut multiplied = rotate(model, left, component)?;
        for simple in word.into_iter().rev() {
            multiplied = simple_product(model, &multiplied, simple)?;
        }
        for image in multiplied {
            for &(shift, scalar) in term.coefficient.iter() {
                product.add(model, &image.window, &image.coefficient, shift, scalar)?;
            }
        }
    }
    Ok(product.finish())
}

pub(crate) fn coefficient<'a>(
    model: &mut Model,
    terms: &'a [Term],
    window: &[i128],
) -> Result<&'a [Monomial]> {
    model.spend(1 + window.len() * (1 + terms.len().max(1).ilog2() as usize))?;
    let result = terms.binary_search_by(|term| term.window.as_ref().cmp(window));
    Ok(result.map_or(&[][..], |index| terms[index].coefficient.as_ref()))
}

pub(crate) fn coefficients(model: &mut Model, coefficient: &[Monomial]) -> Result<Vec<i128>> {
    let size = match coefficient.last() {
        Some(&(degree, _)) => degree.checked_add(1).ok_or_else(overflow)?,
        None => 1,
    };
    check_storage(size.checked_mul(2).ok_or_else(storage_limit)?)?;
    model.spend(size)?;
    let mut result = vec![0; size];
    for &(power, scalar) in coefficient {
        result[power] = scalar;
    }
    Ok(result)
}

fn power(model: &mut Model, mut scalar: i128, mut degree: usize) -> Result<i128> {
    let mut result = 1i128;
    while degree != 0 {
        model.spend(1)?;
        if degree & 1 != 0 {
            result = result.checked_mul(scalar).ok_or_else(overflow)?;
        }
        degree >>= 1;
        if degree != 0 {
            scalar = scalar.checked_mul(scalar).ok_or_else(overflow)?;
        }
    }
    Ok(result)
}

pub(crate) fn evaluate(model: &mut Model, coefficient: &[Monomial], at: i128) -> Result<i128> {
    let mut result = 0i128;
    let mut degree = coefficient.last().map_or(0, |&(degree, _)| degree);
    for &(power_degree, scalar) in coefficient.iter().rev() {
        let factor = power(model, at, degree - power_degree)?;
        result = result
            .checked_mul(factor)
            .and_then(|result| result.checked_add(scalar))
            .ok_or_else(overflow)?;
        degree = power_degree;
    }
    result
        .checked_mul(power(model, at, degree)?)
        .ok_or_else(overflow)
}
