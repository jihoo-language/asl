use std::collections::BTreeMap;
use std::fmt::Write;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;

use super::super::linear;
use super::builder::{Axis, Catalog};
use super::{check_report, multiply_powers, storage, Form, Monomial, Rational, Resolution, Tensor};

struct Unknown {
    i: usize,
    j: usize,
    powers: Monomial,
}

pub(super) fn presentation(
    model: &mut Model,
    form: &Form,
    resolution: &Resolution,
) -> Result<String> {
    let shifts = resolution.shifts();
    let generators = shifts[1];
    let top = shifts[3];
    let size = generators.len();
    if top.len() != 1 || size < 3 || size.is_multiple_of(2) {
        return Err(Error::plain(
            "Verbitsky resolution has inconsistent codimension-three Gorenstein shifts",
        ));
    }
    let bound = top[0];
    let ideal = (0..size)
        .map(|index| {
            resolution
                .entry(1, 0, index)?
                .ok_or_else(|| Error::plain("missing Verbitsky ideal generator"))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut catalog = Catalog::default();
    let axis = Axis::new(model, form, &mut catalog, generators, bound)?;
    let mut unknowns = Vec::new();
    let mut cells = 0;
    for (i, &a) in generators.iter().enumerate() {
        for (j, &b) in generators.iter().enumerate().skip(i + 1) {
            // Constants are excluded by minimality, even for a matching
            // pair of degree shifts in the middle of the resolution.
            let Some(weight) = bound.checked_sub(a + b).filter(|&weight| weight > 0) else {
                continue;
            };
            for powers in catalog.at(model, form, weight)?.iter() {
                cells += 10 + 5 * powers.len();
                storage(cells)?;
                unknowns.push(Unknown {
                    i,
                    j,
                    powers: powers.clone(),
                });
            }
        }
    }
    let mut equations = linear::matrix(axis.len(), unknowns.len())?;
    for (column, entry) in unknowns.iter().enumerate() {
        for (row, generator, minus) in [
            (entry.i, ideal[entry.j], false),
            (entry.j, ideal[entry.i], true),
        ] {
            for (term, &coefficient) in &generator.terms {
                let powers = multiply_powers(term, &entry.powers)?;
                equations[axis.locate(row, &powers)?][column] = if minus {
                    coefficient.negate()
                } else {
                    coefficient
                };
            }
        }
    }
    let solution = linear::rref(model, equations, unknowns.len())?.kernel(unknowns.len())?;
    for row in &solution {
        if let Some(report) = certify(model, form, &ideal, &unknowns, row)? {
            return Ok(report);
        }
    }

    // If each basis element is singular, bounded polynomial nonvanishing
    // search uses the grid 0..m for a form of degree m in each parameter.
    let mut digits = vec![0; solution.len()];
    let middle = (size - 1) / 2;
    while advance(&mut digits, middle) {
        model.spend(1 + solution.len())?;
        if digits.iter().filter(|&&digit| digit != 0).count() < 2 {
            continue;
        }
        let mut values = vec![Rational::zero(); unknowns.len()];
        for (&digit, row) in digits.iter().zip(&solution) {
            linear::add_multiple(model, &mut values, row, Rational::integer(digit as i128))?;
        }
        if let Some(report) = certify(model, form, &ideal, &unknowns, &values)? {
            return Ok(report);
        }
    }
    Err(Error::plain(
        "no full-rank Buchsbaum–Eisenbud matrix was found",
    ))
}

fn advance(digits: &mut [usize], bound: usize) -> bool {
    for digit in digits {
        if *digit < bound {
            *digit += 1;
            return true;
        }
        *digit = 0;
    }
    false
}

fn certify(
    model: &mut Model,
    form: &Form,
    ideal: &[&Tensor],
    unknowns: &[Unknown],
    values: &[Rational],
) -> Result<Option<String>> {
    let mut upper: BTreeMap<(usize, usize), Tensor> = BTreeMap::new();
    let mut cells = 0;
    for (entry, &coefficient) in unknowns.iter().zip(values) {
        if !coefficient.is_zero() {
            let polynomial = upper.entry((entry.i, entry.j)).or_default();
            let before = polynomial.cells;
            polynomial.insert(model, entry.powers.clone(), coefficient)?;
            cells += polynomial.cells - before;
            storage(cells)?;
        }
    }
    let mut pfaffians = Pfaffians {
        upper: &upper,
        memo: BTreeMap::new(),
        cells: 0,
    };
    let first = pfaffians.pfaffian(model, &(1..ideal.len()).collect::<Vec<_>>())?;
    if first.terms.is_empty() {
        return Ok(None);
    }
    let (powers, &coefficient) = ideal[0].terms.first_key_value().unwrap();
    let Some(&pfaffian) = first.terms.get(powers) else {
        return Err(certification_error());
    };
    let scalar = pfaffian.multiply(model, coefficient.inverse()?)?;
    for (i, generator) in ideal.iter().enumerate() {
        let complement: Vec<_> = (0..ideal.len()).filter(|&j| j != i).collect();
        let pfaffian = pfaffians.pfaffian(model, &complement)?;
        let factor = if i.is_multiple_of(2) {
            scalar
        } else {
            scalar.negate()
        };
        if pfaffian.terms.len() != generator.terms.len() {
            return Err(certification_error());
        }
        for ((u, &x), (v, &y)) in pfaffian.terms.iter().zip(&generator.terms) {
            if u != v || x != y.multiply(model, factor)? {
                return Err(certification_error());
            }
        }
    }

    let mut out = format!("scale={} matrix=[", scalar.text(model)?);
    let mut row_separator = "";
    for i in 0..ideal.len() {
        let mut sum = Tensor::default();
        write!(out, "{row_separator}[").unwrap();
        let mut column_separator = "";
        for (j, &generator) in ideal.iter().enumerate() {
            let entry = upper.get(&(i.min(j), i.max(j)));
            if let Some(entry) = entry {
                let negative = i > j;
                let product = entry.cup(model, generator)?;
                add_signed(model, &mut sum, &product, negative)?;
                let polynomial = if negative {
                    entry.scale(model, (-1, 1))?
                } else {
                    entry.clone()
                };
                write!(out, "{column_separator}{}", polynomial.text(model, form)?).unwrap();
            } else {
                write!(out, "{column_separator}{{}}").unwrap();
            }
            column_separator = ", ";
            check_report(&out)?;
        }
        if !sum.terms.is_empty() {
            return Err(certification_error());
        }
        out.push(']');
        row_separator = ", ";
    }
    out.push(']');
    check_report(&out)?;
    Ok(Some(out))
}

struct Pfaffians<'a> {
    upper: &'a BTreeMap<(usize, usize), Tensor>,
    memo: BTreeMap<Vec<usize>, Rc<Tensor>>,
    cells: usize,
}

impl Pfaffians<'_> {
    fn pfaffian(&mut self, model: &mut Model, indices: &[usize]) -> Result<Rc<Tensor>> {
        if indices.is_empty() {
            return Ok(Rc::new(Tensor::unit()));
        }
        if let Some(result) = self.memo.get(indices) {
            return Ok(Rc::clone(result));
        }
        model.spend(indices.len())?;
        let mut sum = Tensor::default();
        for (position, &index) in indices.iter().enumerate().skip(1) {
            let Some(entry) = self.upper.get(&(indices[0], index)) else {
                continue;
            };
            let remaining = indices
                .iter()
                .enumerate()
                .skip(1)
                .filter_map(|(j, &index)| (j != position).then_some(index))
                .collect::<Vec<_>>();
            let minor = self.pfaffian(model, &remaining)?;
            let product = entry.cup(model, &minor)?;
            add_signed(model, &mut sum, &product, position.is_multiple_of(2))?;
        }
        self.cells += 10 + sum.cells + 2 * indices.len();
        storage(self.cells)?;
        let sum = Rc::new(sum);
        self.memo.insert(indices.to_vec(), Rc::clone(&sum));
        Ok(sum)
    }
}

fn add_signed(model: &mut Model, target: &mut Tensor, input: &Tensor, minus: bool) -> Result<()> {
    for (powers, &value) in &input.terms {
        target.insert(
            model,
            powers.clone(),
            if minus { value.negate() } else { value },
        )?;
    }
    Ok(())
}

fn certification_error() -> Error {
    Error::plain("Buchsbaum–Eisenbud matrix failed exact polynomial certification")
}
