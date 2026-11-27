use std::cmp::Ordering;

use crate::error::{Error, Result};
use crate::model::Model;

use super::orthogonal::DIMENSION;
use super::reduced::{grade, Grade};
use super::{arithmetic, decimal, guard, work_limit, Magnitude, WEYL_LIMBS};

struct Factors {
    rho: usize,
    degree: usize,
    dimension: usize,
}

impl Factors {
    fn new(rho: u128, degree: u128) -> Result<Self> {
        if !(1..=(DIMENSION - 2) as u128).contains(&rho) {
            return Err(Error::plain("Hodge polarization rank must be in 1..=21"));
        }
        let transcendental = DIMENSION as u128 - rho;
        if degree == 0 || !transcendental.is_multiple_of(degree) {
            return Err(Error::plain(
                "real multiplication field degree must be a positive divisor of 23 minus the polarization rank",
            ));
        }
        let dimension = transcendental / degree;
        if dimension < 3 {
            return Err(Error::plain(
                "real multiplication requires at least three transcendental directions per real place",
            ));
        }
        Ok(Self {
            rho: rho as usize,
            degree: degree as usize,
            dimension: dimension as usize,
        })
    }
}

pub(super) fn dimension(
    model: &mut Model,
    n: u128,
    d: i128,
    rho: u128,
    degree: u128,
    primitive: bool,
) -> Result<String> {
    let group = Factors::new(rho, degree)?;
    let count = match grade(n, d, primitive)? {
        Grade::Absent => Magnitude::default(),
        Grade::Axis => Magnitude::from_small(1),
        Grade::Interior {
            charge,
            weight_bound,
        } => {
            let cap = series_cap(weight_bound)?;
            let orthogonal = orthogonal_series(model, group.dimension as u128, cap)?;
            let mut neutral = orthogonal.clone();
            for _ in 1..group.degree {
                neutral = convolution(model, &neutral, &orthogonal)?;
            }
            for size in 1..=cap {
                for _ in 1..group.rho {
                    inverse_factor(model, &mut neutral, size)?;
                }
            }
            charged_coefficient(model, &mut neutral, charge, primitive)?
        }
    };
    let name = if primitive {
        "rvHilbRealMultiplicationPrimitiveHodge"
    } else {
        "rvHilbRealMultiplicationHodge"
    };
    Ok(format!(
        "{name}[{n}, {d}, {rho}, {degree}] = {}",
        decimal(model, count)?
    ))
}

pub(super) fn oscillator(model: &mut Model, m: u128, k: u128) -> Result<String> {
    if m < 2 {
        return Err(Error::plain(
            "orthogonal oscillator dimension must be at least two",
        ));
    }
    let cap = series_cap(k)?;
    let b = orthogonal_series(model, m, cap)?;
    Ok(format!(
        "rvHilbSOOscillator[{m}, {k}] = {}",
        decimal(model, b[cap])?
    ))
}

fn series_cap(cap: u128) -> Result<usize> {
    let cap = usize::try_from(cap).map_err(|_| work_limit())?;
    guard(cap.saturating_add(1).saturating_mul(5 * WEYL_LIMBS))?;
    Ok(cap)
}

// Cauchy decomposes Sym(V ⊗ Q^infinity) into Schur functors. An SO(m)
// invariant occurs once when all m padded rows are even (O(m)-trivial)
// or all odd (the determinant representation of O(m)).
fn orthogonal_series(model: &mut Model, m: u128, cap: usize) -> Result<Vec<Magnitude>> {
    let rank = usize::try_from(m).unwrap_or(cap).min(cap);
    let odd_energy = if m <= cap as u128 {
        rank * (rank + 1) / 2
    } else {
        cap + 1
    };
    let mut result = vec![Magnitude::default(); cap + 1];
    let mut scratch = result.clone();
    let mut odd = if odd_energy <= cap {
        vec![1; rank]
    } else {
        Vec::new()
    };
    even_shapes(
        model,
        rank,
        cap,
        0,
        &mut Vec::new(),
        &mut |model, shape, energy| {
            add_shape(model, &mut result[energy..], &mut scratch, shape)?;
            if energy + odd_energy <= cap {
                odd.fill(1);
                for (out, &part) in odd.iter_mut().zip(shape) {
                    *out += part;
                }
                add_shape(
                    model,
                    &mut result[energy + odd_energy..],
                    &mut scratch,
                    &odd,
                )?;
            }
            Ok(())
        },
    )?;
    Ok(result)
}

fn even_shapes(
    model: &mut Model,
    rank: usize,
    budget: usize,
    energy: usize,
    shape: &mut Vec<usize>,
    emit: &mut impl FnMut(&mut Model, &[usize], usize) -> Result<()>,
) -> Result<()> {
    model.spend(shape.len() + 1).map_err(|_| work_limit())?;
    emit(model, shape, energy)?;
    if shape.len() == rank {
        return Ok(());
    }
    let row = shape.len() + 1;
    let bound = budget / row;
    let last = shape.last().copied().unwrap_or(bound).min(bound);
    for part in (2..=last).step_by(2) {
        let weight = row * part;
        shape.push(part);
        even_shapes(model, rank, budget - weight, energy + weight, shape, emit)?;
        shape.pop();
    }
    Ok(())
}

fn add_shape(
    model: &mut Model,
    destination: &mut [Magnitude],
    scratch: &mut [Magnitude],
    shape: &[usize],
) -> Result<()> {
    let scratch = &mut scratch[..destination.len()];
    model
        .spend((WEYL_LIMBS + 1) * scratch.len())
        .map_err(|_| work_limit())?;
    scratch.fill(Magnitude::default());
    scratch[0] = Magnitude::from_small(1);
    for (row, &width) in shape.iter().enumerate() {
        let below = &shape[row + 1..];
        for column in 0..width {
            let hook = width - column + below.partition_point(|&len| len > column);
            inverse_factor(model, scratch, hook)?;
        }
    }
    for (target, &value) in destination.iter_mut().zip(scratch.iter()) {
        target.add_multiple(value, 1).ok_or_else(arithmetic)?;
    }
    Ok(())
}

fn inverse_factor(model: &mut Model, series: &mut [Magnitude], shift: usize) -> Result<()> {
    model
        .spend(WEYL_LIMBS * series.len().saturating_sub(shift))
        .map_err(|_| work_limit())?;
    for index in shift..series.len() {
        let lower = series[index - shift];
        series[index]
            .add_multiple(lower, 1)
            .ok_or_else(arithmetic)?;
    }
    Ok(())
}

fn convolution(
    model: &mut Model,
    left: &[Magnitude],
    right: &[Magnitude],
) -> Result<Vec<Magnitude>> {
    // Charge for the full iteration bound before multiplying exact coefficients.
    model
        .spend(
            WEYL_LIMBS
                .saturating_mul(left.len())
                .saturating_mul(left.len().saturating_add(1))
                / 2,
        )
        .map_err(|_| work_limit())?;
    let mut output = vec![Magnitude::default(); left.len()];
    for (i, &x) in left.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (value, &y) in output[i..].iter_mut().zip(right) {
            if !y.is_zero() {
                value
                    .add_multiple(x.checked_product(y).ok_or_else(arithmetic)?, 1)
                    .ok_or_else(arithmetic)?;
            }
        }
    }
    Ok(output)
}

// After removing Q^a, the coefficient of z^a in the two charged
// oscillators is sum_k Q^(2k)/((Q;Q)_(a+k)(Q;Q)_k), where (Q;Q)_s
// is the product of 1-Q^j for 1 <= j <= s. Increase k incrementally.
fn charged_coefficient(
    model: &mut Model,
    series: &mut [Magnitude],
    charge: u128,
    primitive: bool,
) -> Result<Magnitude> {
    let cap = series.len() - 1;
    for size in 1..=charge.min(cap as u128) as usize {
        inverse_factor(model, series, size)?;
    }
    let mut current = Magnitude::default();
    let mut previous = Magnitude::default();
    for k in 0..=cap / 2 {
        let index = cap - 2 * k;
        let series = &mut series[..=index];
        let shift = charge.saturating_add(k as u128);
        if k > 0 {
            inverse_factor(model, series, k)?;
            inverse_factor(model, series, shift.min(index as u128 + 1) as usize)?;
        }
        current
            .add_multiple(series[index], 1)
            .ok_or_else(arithmetic)?;
        if primitive && index > 0 {
            // The lower degree has excess one less and charge a+1.
            // Its summand adds one more inverse factor of length a+k+1.
            let extra = shift.saturating_add(1).min(index as u128) as usize;
            model
                .spend(WEYL_LIMBS * ((index - 1) / extra + 1))
                .map_err(|_| work_limit())?;
            for offset in (0..index).step_by(extra) {
                previous
                    .add_multiple(series[index - 1 - offset], 1)
                    .ok_or_else(arithmetic)?;
            }
        }
    }
    if current.cmp(previous) == Ordering::Less {
        return Err(Error::plain(
            "primitive real-multiplication Hodge dimension is negative",
        ));
    }
    Ok(current.difference(previous))
}
