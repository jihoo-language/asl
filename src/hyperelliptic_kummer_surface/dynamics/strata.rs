use std::mem::size_of;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::model::Model;
use crate::syntax::Label;

use super::super::{
    hilbert::{
        checked_product, field_powers, graded_index, intermediate_overflow,
        k3_symmetric_coefficients, power, symmetric_partition_product, Integer, Magnitude,
    },
    stratum::Partition,
};
use super::{
    semilinear::{middle_traces, node_permutation, word_and_base},
    Word,
};

pub(super) fn report(
    model: &mut Model,
    word_name: &Label,
    base: &Label,
    shape_name: &Label,
    r: u128,
    i: u128,
    coefficient: Option<u128>,
) -> Result<String> {
    let (word, curve) = word_and_base(model, word_name, base)?;
    let shape = model.hilbert_stratum(shape_name)?;
    word.basis
        .require_base_extension(model, &curve)
        .map_err(|error| error.at(base.span))?;
    let value = if coefficient == Some(0) {
        Integer::from_i128(1)
    } else if let Some((degree, twist)) = stratum_grade(&shape, i) {
        let action = Action {
            word: &word,
            curve: &curve,
            base,
            r,
        };
        evaluate(model, &action, &shape, degree, twist, coefficient)?
    } else {
        Integer::from_i128(0)
    };
    let value = value.signed().ok_or_else(|| {
        Error::plain("Hilbert–Chow stratum trace or Euler factor exceeds signed 128-bit range")
    })?;
    let (word, base, shape) = (&word_name.text, &base.text, &shape_name.text);
    Ok(match coefficient {
        Some(j) => {
            format!("rvK3StratumEulerFactor[{word}, {base}, {shape}, {r}, {i}, {j}] = {value}")
        }
        None => format!("rvK3StratumFrobTrace[{word}, {base}, {shape}, {r}, {i}] = {value}"),
    })
}

fn stratum_grade(shape: &Partition, i: u128) -> Option<(u128, u128)> {
    let half = (i / 2).checked_sub(shape.shift())?;
    let unshifted = 2 * half + i % 2;
    let (_, degree, twist) = graded_index(shape.length, unshifted)?;
    Some((degree, twist + shape.shift()))
}

struct Action<'a> {
    word: &'a Word,
    curve: &'a Curve,
    base: &'a Label,
    r: u128,
}

fn evaluate(
    model: &mut Model,
    action: &Action<'_>,
    shape: &Partition,
    degree: u128,
    twist: u128,
    coefficient: Option<u128>,
) -> Result<Integer> {
    let steps = coefficient.unwrap_or(1);
    if degree == 0 {
        if steps > 1 {
            return Ok(Integer::from_i128(0));
        }
        let value = tate(model, action, 1, twist, Integer::from_i128(1))?;
        return Ok(if coefficient.is_some() {
            value.negate()
        } else {
            value
        });
    }
    let length = shape
        .parts()
        .iter()
        .map(|&(_, c)| c)
        .max()
        .unwrap_or(0)
        .min(degree);
    if steps > 1 {
        let rank = partition_value(model, shape, length, degree, None)?;
        if rank.magnitude.cmp(Magnitude::from_small(steps)).is_lt() {
            return Ok(Integer::from_i128(0));
        }
    }
    let steps = usize::try_from(steps).map_err(|_| storage_limit())?;
    let words = size_of::<Integer>().div_ceil(size_of::<usize>());
    if steps.saturating_add(1).saturating_mul(2 * words) > MAX_PIVOT_CELLS {
        return Err(storage_limit());
    }
    let length_usize = usize::try_from(length).map_err(|_| storage_limit())?;
    let through = length_usize
        .checked_mul(steps)
        .ok_or_else(intermediate_overflow)?;
    let last = action
        .r
        .checked_mul(through as u128)
        .ok_or_else(intermediate_overflow)?;
    let q = action.curve.ring.size_of_residue_field()?;
    // Bound the series before allocating it. The same marked matrix yields
    // all Ψ^(r*s*k) traces used by the symmetric and determinant recurrences.
    power(model, q, last)?;
    let fields = field_powers(model, q, action.r, through)?;
    let permutation = node_permutation(model, action.word, action.curve, action.base)?;
    let middle = middle_traces(
        model,
        action.word,
        action.curve,
        &permutation,
        action.r,
        &fields,
        intermediate_overflow,
    )?;
    model.spend(steps)?;
    let mut traces = Vec::with_capacity(steps);
    for step in 1..=steps {
        let value = partition_value(model, shape, length, degree, Some((&fields, &middle, step)))?;
        traces.push(tate(model, action, step, twist, value)?);
    }
    if coefficient.is_some() {
        newton_coefficient(model, &traces)
    } else {
        Ok(traces[0])
    }
}

fn partition_value(
    model: &mut Model,
    shape: &Partition,
    length: u128,
    degree: u128,
    traces: Option<(&[Magnitude], &[Integer], usize)>,
) -> Result<Integer> {
    let rows = k3_symmetric_coefficients(model, length, degree, traces)?;
    let product = symmetric_partition_product(model, shape.parts(), &rows)?;
    Ok(*product.last().unwrap())
}

fn tate(
    model: &mut Model,
    action: &Action<'_>,
    step: usize,
    twist: u128,
    value: Integer,
) -> Result<Integer> {
    if twist == 0 || value.is_zero() {
        return Ok(value);
    }
    let exponent = action
        .r
        .checked_mul(step as u128)
        .and_then(|r| r.checked_mul(twist))
        .ok_or_else(intermediate_overflow)?;
    let q = action.curve.ring.size_of_residue_field()?;
    let magnitude = power(model, q, exponent)?;
    checked_product(
        model,
        value,
        Integer {
            magnitude,
            negative: false,
        },
    )
}

fn newton_coefficient(model: &mut Model, traces: &[Integer]) -> Result<Integer> {
    model.spend(traces.len() + 1)?;
    let mut coefficients = Vec::with_capacity(traces.len() + 1);
    coefficients.push(Integer::from_i128(1));
    for j in 1..=traces.len() {
        let mut sum = Integer::from_i128(0);
        for (&trace, &coefficient) in traces.iter().zip(coefficients.iter().rev()) {
            let term = checked_product(model, trace, coefficient)?;
            sum = sum.add(term).map_err(|_| intermediate_overflow())?;
        }
        let (magnitude, remainder) = sum.magnitude.divide(model, j as u128)?;
        if remainder != 0 {
            return Err(Error::plain(
                "Hilbert–Chow Euler factor Newton division is inexact",
            ));
        }
        coefficients.push(Integer {
            magnitude,
            negative: !sum.negative,
        });
    }
    Ok(*coefficients.last().unwrap())
}

fn storage_limit() -> Error {
    Error::plain("Hilbert–Chow Euler factor exceeds the coefficient storage limit")
}
