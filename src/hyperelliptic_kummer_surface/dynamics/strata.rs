use std::mem::size_of;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::model::Model;
use crate::rv_witt_spectrum::{self, determinant_from_ghosts, Bridge, Spectrum};
use crate::syntax::Label;

use super::super::{
    hilbert::{
        checked_product, field_powers, graded_index, intermediate_overflow,
        k3_symmetric_coefficients, power, symmetric_partition_product, Integer, Magnitude,
    },
    specht::Decoration,
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
    let shape = model.hilbert_stratum(shape_name)?;
    report_summand(
        model,
        word_name,
        base,
        shape_name,
        Summand::Hilbert(&shape),
        Index { r, i, coefficient },
    )
}

pub(super) fn specht_report(
    model: &mut Model,
    word_name: &Label,
    base: &Label,
    shape_name: &Label,
    r: u128,
    i: u128,
    coefficient: Option<u128>,
) -> Result<String> {
    let shape = model.specht_stratum(shape_name)?;
    report_summand(
        model,
        word_name,
        base,
        shape_name,
        Summand::Specht(&shape),
        Index { r, i, coefficient },
    )
}

pub(crate) fn spectrum(model: &mut Model, bridge: &Bridge) -> Result<Spectrum> {
    let (word, curve) = word_and_base(model, &bridge.word, &bridge.base)?;
    word.basis
        .require_base_extension(model, &curve)
        .map_err(|error| error.at(bridge.base.span))?;
    let shape;
    let decoration;
    let summand = if bridge.specht {
        decoration = model.specht_stratum(&bridge.stratum)?;
        Summand::Specht(&decoration)
    } else {
        shape = model.hilbert_stratum(&bridge.stratum)?;
        Summand::Hilbert(&shape)
    };
    let Some((degree, twist)) = summand.grade(model, bridge.i)? else {
        return Ok(Spectrum::zero());
    };
    let action = Action {
        word: &word,
        curve: &curve,
        base: &bridge.base,
        r: bridge.r,
    };
    if degree == 0 {
        let eigenvalue = tate(model, &action, 1, twist, Integer::from_i128(1))?;
        return Spectrum::effective(vec![Integer::from_i128(1), eigenvalue.negate()]);
    }
    let length = summand.symmetric_length(degree);
    let rank = partition_value(model, &summand, length, degree, None)?;
    let rank = rank
        .magnitude
        .small()
        .and_then(|n| usize::try_from(n).ok())
        .unwrap_or(usize::MAX);
    rv_witt_spectrum::require_degree(rank)?;
    if rank == 0 {
        return Ok(Spectrum::zero());
    }
    let traces = power_traces(model, &action, &summand, length, degree, twist, rank)?;
    Spectrum::effective(determinant_from_ghosts(model, &traces)?)
}

struct Index {
    r: u128,
    i: u128,
    coefficient: Option<u128>,
}

fn report_summand(
    model: &mut Model,
    word_name: &Label,
    base: &Label,
    shape_name: &Label,
    summand: Summand<'_>,
    index: Index,
) -> Result<String> {
    let Index { r, i, coefficient } = index;
    let (word, curve) = word_and_base(model, word_name, base)?;
    word.basis
        .require_base_extension(model, &curve)
        .map_err(|error| error.at(base.span))?;
    let value = if coefficient == Some(0) {
        Integer::from_i128(1)
    } else if let Some((degree, twist)) = summand.grade(model, i)? {
        let action = Action {
            word: &word,
            curve: &curve,
            base,
            r,
        };
        evaluate(model, &action, &summand, degree, twist, coefficient)?
    } else {
        Integer::from_i128(0)
    };
    let value = value.signed().ok_or_else(|| {
        Error::plain("Hilbert–Chow stratum trace or Euler factor exceeds signed 128-bit range")
    })?;
    let (word, base, shape) = (&word_name.text, &base.text, &shape_name.text);
    let operation = match (&summand, coefficient.is_some()) {
        (Summand::Hilbert(_), false) => "rvK3StratumFrobTrace",
        (Summand::Hilbert(_), true) => "rvK3StratumEulerFactor",
        (Summand::Specht(_), false) => "rvK3SpechtFrobTrace",
        (Summand::Specht(_), true) => "rvK3SpechtEulerFactor",
    };
    Ok(match coefficient {
        Some(j) => format!("{operation}[{word}, {base}, {shape}, {r}, {i}, {j}] = {value}"),
        None => format!("{operation}[{word}, {base}, {shape}, {r}, {i}] = {value}"),
    })
}

enum Summand<'a> {
    Hilbert(&'a Partition),
    Specht(&'a Decoration),
}

impl Summand<'_> {
    fn shape(&self) -> &Partition {
        match self {
            Self::Hilbert(shape) => shape,
            Self::Specht(decoration) => &decoration.shape,
        }
    }

    fn grade(&self, model: &mut Model, i: u128) -> Result<Option<(u128, u128)>> {
        let shift = match self {
            Self::Hilbert(shape) => shape.shift(),
            Self::Specht(_) => 0,
        };
        let Some(half) = (i / 2).checked_sub(shift) else {
            return Ok(None);
        };
        let Some((_, degree, twist)) = graded_index(self.shape().length, 2 * half + i % 2) else {
            return Ok(None);
        };
        if let Self::Specht(decoration) = self {
            if !decoration.supports_degree(model, degree)? {
                return Ok(None);
            }
        }
        Ok(Some((degree, twist + shift)))
    }

    fn symmetric_length(&self, degree: u128) -> u128 {
        self.shape()
            .parts()
            .iter()
            .map(|&(_, count)| count)
            .max()
            .unwrap_or(0)
            .min(degree)
    }
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
    summand: &Summand<'_>,
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
    let length = summand.symmetric_length(degree);
    if steps > 1 {
        let rank = partition_value(model, summand, length, degree, None)?;
        if rank.magnitude.cmp(Magnitude::from_small(steps)).is_lt() {
            return Ok(Integer::from_i128(0));
        }
    }
    let steps = usize::try_from(steps).map_err(|_| storage_limit())?;
    let traces = power_traces(model, action, summand, length, degree, twist, steps)?;
    if coefficient.is_some() {
        Ok(*determinant_from_ghosts(model, &traces)?.last().unwrap())
    } else {
        Ok(traces[0])
    }
}

fn power_traces(
    model: &mut Model,
    action: &Action<'_>,
    summand: &Summand<'_>,
    length: u128,
    degree: u128,
    twist: u128,
    steps: usize,
) -> Result<Vec<Integer>> {
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
        let value = partition_value(
            model,
            summand,
            length,
            degree,
            Some((&fields, &middle, step)),
        )?;
        traces.push(tate(model, action, step, twist, value)?);
    }
    Ok(traces)
}

fn partition_value(
    model: &mut Model,
    summand: &Summand<'_>,
    length: u128,
    degree: u128,
    traces: Option<(&[Magnitude], &[Integer], usize)>,
) -> Result<Integer> {
    let rows = k3_symmetric_coefficients(model, length, degree, traces)?;
    match summand {
        Summand::Hilbert(shape) => {
            let product = symmetric_partition_product(model, shape.parts(), &rows)?;
            Ok(*product.last().unwrap())
        }
        Summand::Specht(decoration) => decoration.coefficient(model, &rows),
    }
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

fn storage_limit() -> Error {
    Error::plain("Hilbert–Chow Euler factor exceeds the coefficient storage limit")
}
