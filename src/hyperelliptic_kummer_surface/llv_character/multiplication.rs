use std::cmp::Ordering;

use crate::error::{Error, Result};
use crate::model::Model;

use super::orthogonal::DIMENSION;
use super::reduced::{grade, Grade};
use super::{arithmetic, decimal, guard, work_limit, Magnitude, WEYL_LIMBS};

#[derive(Clone, Copy)]
pub(super) enum Group {
    Orthogonal,
    Unitary,
}

struct Factors {
    rho: usize,
    degree: usize,
    dimension: usize,
}

impl Factors {
    fn new(group: Group, rho: u128, degree: u128) -> Result<Self> {
        if !(1..=(DIMENSION - 2) as u128).contains(&rho) {
            return Err(Error::plain("Hodge polarization rank must be in 1..=21"));
        }
        let transcendental = DIMENSION as u128 - rho;
        let dimension = match group {
            Group::Orthogonal => {
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
                dimension
            }
            Group::Unitary => {
                if degree == 0
                    || degree > transcendental / 2
                    || !transcendental.is_multiple_of(2 * degree)
                {
                    return Err(Error::plain(
                        "CM field degree 2e must be a positive divisor of 23 minus the polarization rank",
                    ));
                }
                transcendental / (2 * degree)
            }
        };
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
    group: Group,
    primitive: bool,
) -> Result<String> {
    let factors = Factors::new(group, rho, degree)?;
    let count = match grade(n, d, primitive)? {
        Grade::Absent => Magnitude::default(),
        Grade::Axis => Magnitude::from_small(1),
        Grade::Interior {
            charge,
            weight_bound,
        } => {
            let cap = series_cap(weight_bound)?;
            let factor = match group {
                Group::Orthogonal => orthogonal_series(model, factors.dimension as u128, cap)?,
                Group::Unitary => unitary_series(model, factors.dimension as u128, cap, 0)?,
            };
            let mut neutral = factor.clone();
            for _ in 1..factors.degree {
                neutral = convolution(model, &neutral, &factor)?;
            }
            for size in 1..=cap {
                for _ in 1..factors.rho {
                    inverse_factor(model, &mut neutral, size)?;
                }
            }
            charged_coefficient(model, &mut neutral, charge, primitive)?
        }
    };
    let name = match (group, primitive) {
        (Group::Orthogonal, false) => "rvHilbRealMultiplicationHodge",
        (Group::Orthogonal, true) => "rvHilbRealMultiplicationPrimitiveHodge",
        (Group::Unitary, false) => "rvHilbComplexMultiplicationHodge",
        (Group::Unitary, true) => "rvHilbComplexMultiplicationPrimitiveHodge",
    };
    Ok(format!(
        "{name}[{n}, {d}, {rho}, {degree}] = {}",
        decimal(model, count)?
    ))
}

pub(super) fn so_oscillator(model: &mut Model, m: u128, k: u128) -> Result<String> {
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

pub(super) fn gl_oscillator(model: &mut Model, m: u128, k: u128, s: i128) -> Result<String> {
    if m == 0 {
        return Err(Error::plain(
            "general linear oscillator dimension must be positive",
        ));
    }
    let shift = s.unsigned_abs();
    let excess = rectangle_energy(m, shift).and_then(|energy| k.checked_sub(energy));
    let count = match excess {
        None => Magnitude::default(),
        Some(0) => Magnitude::from_small(1),
        Some(excess) => {
            let cap = series_cap(excess)?;
            unitary_series(model, m, cap, shift)?[cap]
        }
    };
    Ok(format!(
        "rvHilbGLOscillator[{m}, {k}, {s}] = {}",
        decimal(model, count)?
    ))
}

fn rectangle_energy(rank: u128, shift: u128) -> Option<u128> {
    if shift == 0 {
        return Some(0);
    }
    let triangular = if rank.is_multiple_of(2) {
        (rank / 2).checked_mul(rank + 1)?
    } else {
        rank.checked_mul(rank / 2 + 1)?
    };
    triangular.checked_mul(shift)
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
    let mut result = vec![Magnitude::default(); cap + 1];
    let mut scratch = result.clone();
    let odd_energy = rectangle_energy(m, 1)
        .filter(|&energy| energy <= cap as u128)
        .map(|energy| energy as usize);
    partitions(
        model,
        rank,
        2,
        cap,
        0,
        &mut Vec::new(),
        &mut |model, shape, energy| {
            add_schur_term(
                model,
                &mut result[energy..],
                &mut scratch,
                &[SchurShape::plain(shape)],
            )?;
            if let Some(offset) = odd_energy.filter(|&odd| energy + odd <= cap) {
                add_schur_term(
                    model,
                    &mut result[energy + offset..],
                    &mut scratch,
                    &[SchurShape {
                        shape,
                        rank: m,
                        shift: 1,
                    }],
                )?;
            }
            Ok(())
        },
    )?;
    Ok(result)
}

// Opposite Schur factors pair under GL(m). In determinant sector s
// the two diagrams differ by the rectangle (|s|^m), whose energy was
// removed before choosing the truncation cap.
fn unitary_series(model: &mut Model, m: u128, cap: usize, shift: u128) -> Result<Vec<Magnitude>> {
    let rank = usize::try_from(m).unwrap_or(cap).min(cap);
    let mut result = vec![Magnitude::default(); cap + 1];
    let mut scratch = result.clone();
    partitions(
        model,
        rank,
        1,
        cap / 2,
        0,
        &mut Vec::new(),
        &mut |model, shape, energy| {
            add_schur_term(
                model,
                &mut result[2 * energy..],
                &mut scratch,
                &[
                    SchurShape::plain(shape),
                    SchurShape {
                        shape,
                        rank: m,
                        shift,
                    },
                ],
            )
        },
    )?;
    Ok(result)
}

fn partitions(
    model: &mut Model,
    rank: usize,
    step: usize,
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
    for part in (step..=last).step_by(step) {
        let weight = row * part;
        shape.push(part);
        partitions(
            model,
            rank,
            step,
            budget - weight,
            energy + weight,
            shape,
            emit,
        )?;
        shape.pop();
    }
    Ok(())
}

struct SchurShape<'a> {
    shape: &'a [usize],
    rank: u128,
    shift: u128,
}

impl<'a> SchurShape<'a> {
    fn plain(shape: &'a [usize]) -> Self {
        Self {
            shape,
            rank: shape.len() as u128,
            shift: 0,
        }
    }

    fn denominator(&self, model: &mut Model, series: &mut [Magnitude]) -> Result<()> {
        let cap = (series.len() - 1) as u128;
        for (row, &extra) in self.shape.iter().enumerate() {
            let width = self
                .shift
                .checked_add(extra as u128)
                .ok_or_else(arithmetic)?;
            let lower = &self.shape[row + 1..];
            // A factor with hook greater than cap is 1 in the truncated
            // series. Only the rightmost cap columns can contribute.
            for arm in 1..=width.min(cap) {
                let column = width - arm;
                let leg = if column < self.shift {
                    self.rank - row as u128 - 1
                } else {
                    let column = (column - self.shift) as usize;
                    lower.partition_point(|&width| width > column) as u128
                };
                let hook = arm.saturating_add(leg);
                if hook <= cap {
                    inverse_factor(model, series, hook as usize)?;
                }
            }
        }
        if self.shift > 0 {
            // In the rectangular tail only the bottom cap rows contribute.
            for leg in 0..(self.rank - self.shape.len() as u128).min(cap) {
                for hook in (leg + 1)..=cap.min(leg.saturating_add(self.shift)) {
                    inverse_factor(model, series, hook as usize)?;
                }
            }
        }
        Ok(())
    }
}

fn add_schur_term(
    model: &mut Model,
    destination: &mut [Magnitude],
    scratch: &mut [Magnitude],
    shapes: &[SchurShape<'_>],
) -> Result<()> {
    let scratch = &mut scratch[..destination.len()];
    model
        .spend((WEYL_LIMBS + 1) * scratch.len())
        .map_err(|_| work_limit())?;
    scratch.fill(Magnitude::default());
    scratch[0] = Magnitude::from_small(1);
    for shape in shapes {
        shape.denominator(model, scratch)?;
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
            "primitive Hodge invariant dimension is negative",
        ));
    }
    Ok(current.difference(previous))
}
