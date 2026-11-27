use std::cmp::Ordering;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};

const MAX_ORDER: usize = 4096;
pub(crate) type Polynomial = Rc<[i128]>;

pub(crate) fn check_order(order: usize) -> Result<()> {
    if !(1..=MAX_ORDER).contains(&order) {
        return Err(Error::plain(format!(
            "cyclotomic order must be between 1 and {MAX_ORDER}"
        )));
    }
    Ok(())
}

fn arithmetic() -> Error {
    Error::plain("cyclotomic integer arithmetic is out of range")
}

fn trim(coefficients: &mut Vec<i128>) {
    coefficients.truncate(
        coefficients
            .iter()
            .rposition(|&c| c != 0)
            .map_or(0, |i| i + 1),
    );
}

// Both exact division in the factorization z^m-1 and reduction in the
// quotient have the same elimination step, since every divisor is monic.
fn divide(
    model: &mut Model,
    dividend: &mut Vec<i128>,
    divisor: &[i128],
    mut quotient: impl FnMut(usize, i128),
) -> Result<()> {
    let degree = divisor.len() - 1;
    for index in (degree..dividend.len()).rev() {
        let factor = dividend[index];
        model.spend(1)?;
        if factor == 0 {
            continue;
        }
        let start = index - degree;
        quotient(start, factor);
        model.spend(degree)?;
        for (entry, &c) in dividend[start..index].iter_mut().zip(divisor) {
            *entry = factor
                .checked_mul(c)
                .and_then(|term| entry.checked_sub(term))
                .ok_or_else(arithmetic)?;
        }
        dividend[index] = 0;
    }
    dividend.truncate(degree);
    trim(dividend);
    Ok(())
}

pub(crate) fn modulus(model: &mut Model, order: usize) -> Result<Polynomial> {
    check_order(order)?;
    model.spend(order)?;
    let mut factors = Vec::<(usize, Vec<i128>)>::new();
    for d in (1..=order).filter(|d| order.is_multiple_of(*d)) {
        let mut power = vec![0; d + 1];
        power[0] = -1;
        power[d] = 1;
        for (e, divisor) in &factors {
            if !d.is_multiple_of(*e) {
                continue;
            }
            let mut quotient = vec![0; power.len() - divisor.len() + 1];
            divide(model, &mut power, divisor, |i, c| quotient[i] = c)?;
            if !power.is_empty() {
                return Err(Error::plain("cyclotomic factor division failed"));
            }
            power = quotient;
        }
        factors.push((d, power));
    }
    Ok(factors.pop().unwrap().1.into())
}

pub(crate) fn persist(model: &mut Model, mut coefficients: Vec<i128>) -> Result<Value> {
    trim(&mut coefficients);
    model.allocate(2 * coefficients.len())?;
    Ok(Value::RVCyclo(coefficients.into()))
}

pub(crate) fn zero() -> Value {
    Value::RVCyclo([].into())
}

pub(crate) fn finish(model: &mut Model, order: usize, coefficients: Vec<i128>) -> Result<Value> {
    let coefficients = reduced(model, order, coefficients)?;
    persist(model, coefficients)
}

pub(crate) fn reduced(
    model: &mut Model,
    order: usize,
    mut coefficients: Vec<i128>,
) -> Result<Vec<i128>> {
    let modulus = model.cyclotomic(order)?;
    divide(model, &mut coefficients, &modulus, |_, _| {})?;
    Ok(coefficients)
}

pub(crate) fn root(model: &mut Model, order: usize, power: i128) -> Result<Value> {
    let power = power.rem_euclid(order as i128) as usize;
    let mut coefficients = vec![0; power + 1];
    coefficients[power] = 1;
    finish(model, order, coefficients)
}

pub(crate) fn compare(model: &mut Model, a: &Polynomial, b: &Polynomial) -> Result<Ordering> {
    if Rc::ptr_eq(a, b) {
        return Ok(Ordering::Equal);
    }
    model.spend(a.len().min(b.len()))?;
    Ok(a.cmp(b))
}

pub(crate) fn add(model: &mut Model, a: &[i128], b: &[i128]) -> Result<Value> {
    let mut sum = a.to_vec();
    accumulate(model, &mut sum, b, false)?;
    persist(model, sum)
}

pub(crate) fn accumulate(
    model: &mut Model,
    target: &mut Vec<i128>,
    coefficients: &[i128],
    subtract: bool,
) -> Result<()> {
    model.spend(target.len().max(coefficients.len()))?;
    target.resize(target.len().max(coefficients.len()), 0);
    for (entry, &coefficient) in target.iter_mut().zip(coefficients) {
        *entry = if subtract {
            entry.checked_sub(coefficient)
        } else {
            entry.checked_add(coefficient)
        }
        .ok_or_else(arithmetic)?;
    }
    trim(target);
    Ok(())
}

pub(crate) fn scale(model: &mut Model, coefficients: &[i128], factor: i128) -> Result<Value> {
    if factor == 0 {
        return Ok(zero());
    }
    model.spend(coefficients.len())?;
    let sum = coefficients
        .iter()
        .map(|entry| entry.checked_mul(factor).ok_or_else(arithmetic))
        .collect::<Result<Vec<_>>>()?;
    persist(model, sum)
}

pub(crate) fn multiply(model: &mut Model, order: usize, a: &[i128], b: &[i128]) -> Result<Value> {
    let coefficients = product_coefficients(model, order, a, b)?;
    persist(model, coefficients)
}

pub(crate) fn product_coefficients(
    model: &mut Model,
    order: usize,
    a: &[i128],
    b: &[i128],
) -> Result<Vec<i128>> {
    if a.is_empty() || b.is_empty() {
        return Ok(Vec::new());
    }
    model.spend(a.len() + b.len())?;
    let a: Vec<_> = a
        .iter()
        .copied()
        .enumerate()
        .filter(|(_, c)| *c != 0)
        .collect();
    let b: Vec<_> = b
        .iter()
        .copied()
        .enumerate()
        .filter(|(_, c)| *c != 0)
        .collect();
    model.spend(a.len() * b.len())?;
    let mut coefficients = vec![0i128; order];
    for &(i, x) in &a {
        for &(j, y) in &b {
            let entry = &mut coefficients[(i + j) % order];
            *entry = x
                .checked_mul(y)
                .and_then(|term| entry.checked_add(term))
                .ok_or_else(arithmetic)?;
        }
    }
    reduced(model, order, coefficients)
}

pub(crate) fn shifted(
    model: &mut Model,
    order: usize,
    source: &[i128],
    exponent: i128,
) -> Result<Vec<i128>> {
    if source.is_empty() {
        return Ok(Vec::new());
    }
    model.spend(source.len())?;
    let power = exponent.rem_euclid(order as i128) as usize;
    if power == 0 {
        return Ok(source.to_vec());
    }
    let mut result = vec![0; order];
    for (i, &coefficient) in source.iter().enumerate() {
        result[(i + power) % order] = coefficient;
    }
    reduced(model, order, result)
}

pub(crate) fn embed(
    model: &mut Model,
    source: usize,
    target: usize,
    coefficients: &[i128],
) -> Result<Value> {
    if !target.is_multiple_of(source) {
        return Err(Error::plain(
            "cyclotomic embedding requires the source order to divide the target order",
        ));
    }
    let step = target / source;
    model.spend(coefficients.len())?;
    let mut result = vec![0; target];
    for (i, &coefficient) in coefficients.iter().enumerate() {
        result[step * i] = coefficient;
    }
    finish(model, target, result)
}

pub(crate) fn conjugate(model: &mut Model, order: usize, coefficients: &[i128]) -> Result<Value> {
    let coefficients = substitute_root(model, order, coefficients, order - 1)?;
    persist(model, coefficients)
}

pub(crate) fn substitute_root(
    model: &mut Model,
    order: usize,
    coefficients: &[i128],
    power: usize,
) -> Result<Vec<i128>> {
    let power = power % order;
    let mut changed = vec![0i128; order];
    model.spend(order + coefficients.len())?;
    for (index, &coefficient) in coefficients.iter().enumerate() {
        let value = &mut changed[(power * index) % order];
        *value = value.checked_add(coefficient).ok_or_else(arithmetic)?;
    }
    reduced(model, order, changed)
}
