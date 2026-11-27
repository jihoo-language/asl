use std::collections::BTreeMap;
use std::fmt::Write;

use crate::error::{Error, Result};
use crate::model::{Model, MAX_SHOWN_BYTES};

use super::weyl::{self, RootSystem};
use super::{arithmetic, decimal, guard, report_limit, work_limit};
use super::{Magnitude, Oscillators};

const RANK: usize = 11;

enum Grade {
    Absent,
    Axis,
    Interior { charge: u128, weight_bound: u128 },
}

fn grade(n: u128, d: i128, primitive: bool) -> Result<Grade> {
    if n == 1 {
        return Err(Error::plain(
            "these LLV reports exclude length one, whose LLV type is D12",
        ));
    }
    if d < 0 || d % 2 != 0 {
        return Ok(Grade::Absent);
    }
    let half = d as u128 / 2;
    let charge = n.abs_diff(half);
    if charge > n || (primitive && half > n) {
        return Ok(Grade::Absent);
    }
    let weight_bound = n - charge;
    if weight_bound == 0 {
        Ok(Grade::Axis)
    } else {
        Ok(Grade::Interior {
            charge,
            weight_bound,
        })
    }
}

pub(super) fn multiplicity(
    model: &mut Model,
    n: u128,
    d: i128,
    mut shape: Vec<u128>,
) -> Result<String> {
    trim(&mut shape);
    if shape.len() > RANK || !shape.windows(2).all(|pair| pair[0] >= pair[1]) {
        return Err(Error::plain(
            "type B11 highest weight must be weakly decreasing with at most eleven nonzero coordinates",
        ));
    }
    let count = match grade(n, d, false)? {
        Grade::Axis if shape.is_empty() => 1,
        Grade::Interior {
            charge,
            weight_bound,
        } => match support_excess(&shape, weight_bound) {
            None => 0,
            Some(0) => 1,
            Some(_) => {
                Oscillators::new(model, n)?.extracted(model, &shape, RANK, Some(charge as usize))?
            }
        },
        _ => 0,
    };
    Ok(format!(
        "rvHilbReducedLLVMultiplicity[{n}, {d}]{shape:?} = {count}"
    ))
}

pub(super) fn primitive_multiplicity(
    model: &mut Model,
    n: u128,
    d: i128,
    mut shape: Vec<i128>,
) -> Result<String> {
    trim(&mut shape);
    if shape.len() > RANK
        || shape.iter().take(RANK - 1).any(|&part| part < 0)
        || !shape
            .windows(2)
            .all(|pair| pair[0] >= 0 && pair[0] as u128 >= pair[1].unsigned_abs())
    {
        return Err(Error::plain(
            "type D11 highest weight requires eleven or fewer coordinates with nu1 >= ... >= nu10 >= |nu11|",
        ));
    }
    // The two possible signs at the last coordinate have the same branching
    // multiplicity from B11. The sign is still retained in displayed weights.
    let absolute: Vec<_> = shape.iter().map(|part| part.unsigned_abs()).collect();
    let count = match grade(n, d, true)? {
        Grade::Axis if shape.is_empty() => 1,
        Grade::Interior {
            charge,
            weight_bound,
        } if support_excess(&absolute, weight_bound).is_some() => {
            let oscillators = Oscillators::new(model, n)?;
            let charge = charge as usize;
            let mut lower = [0; RANK];
            lower[..absolute.len()].copy_from_slice(&absolute);
            let mut suffix = [0; RANK + 1];
            for i in (0..RANK).rev() {
                suffix[i] = (i as u128 + 1) * lower[i] + suffix[i + 1];
            }
            let mut count = 0i128;
            interlacers(
                model,
                &lower,
                &suffix,
                weight_bound,
                &mut Vec::new(),
                &mut |model, upper| {
                    let difference = degree_difference(&oscillators, model, upper, charge, true)?;
                    count = count.checked_add(difference).ok_or_else(arithmetic)?;
                    Ok(())
                },
            )?;
            nonnegative(count)?
        }
        _ => 0,
    };
    Ok(format!(
        "rvHilbPrimitiveLLVMultiplicity[{n}, {d}]{shape:?} = {count}"
    ))
}

fn trim<T: Default + PartialEq>(shape: &mut Vec<T>) {
    while shape.last() == Some(&T::default()) {
        shape.pop();
    }
}

// Restore one neutral oscillator to view the eleven pairs as the standard
// SO(23) representation. In the symmetric powers of its length-m copies,
// a harmonic Young column of height h needs h distinct modes by column
// antisymmetry. Its least cost is 1 + ... + h. Orthogonal contractions add
// pairs of tensor factors, and the determinant-twisted orthogonal lift has
// a taller first column. The resulting minimum is sum i*mu_i (one based).
// At equality the unique semistandard filling uses the row number as mode.
// The neutral correction and normalized degree oscillator have constant 1.
fn support_excess(shape: &[u128], budget: u128) -> Option<u128> {
    shape
        .iter()
        .enumerate()
        .try_fold(budget, |left, (i, &part)| {
            left.checked_sub(part.checked_mul(i as u128 + 1)?)
        })
}

fn coefficient(
    oscillators: &Oscillators,
    model: &mut Model,
    shape: &[u128],
    charge: usize,
) -> Result<i128> {
    match oscillators
        .n
        .checked_sub(charge)
        .and_then(|bound| support_excess(shape, bound as u128))
    {
        None => Ok(0),
        Some(0) => Ok(1),
        Some(_) => oscillators.extracted(model, shape, RANK, Some(charge)),
    }
}

fn degree_difference(
    oscillators: &Oscillators,
    model: &mut Model,
    shape: &[u128],
    charge: usize,
    primitive: bool,
) -> Result<i128> {
    let count = coefficient(oscillators, model, shape, charge)?;
    if primitive {
        // On the primitive half, lowering the cohomological degree by 2
        // raises its absolute degree charge by 1.
        count
            .checked_sub(coefficient(oscillators, model, shape, charge + 1)?)
            .ok_or_else(arithmetic)
    } else {
        Ok(count)
    }
}

fn nonnegative(count: i128) -> Result<i128> {
    if count < 0 {
        Err(Error::plain(
            "primitive LLV branching produced a negative multiplicity",
        ))
    } else {
        Ok(count)
    }
}

pub(super) fn decomposition(
    model: &mut Model,
    n: u128,
    d: i128,
    primitive: bool,
) -> Result<String> {
    let name = if primitive {
        "rvHilbPrimitiveLLVDecomposition"
    } else {
        "rvHilbReducedLLVDecomposition"
    };
    let mut report = format!("{name}[{n}, {d}] = components=[");
    let mut components = BTreeMap::new();
    match grade(n, d, primitive)? {
        Grade::Absent => {}
        Grade::Axis => {
            components.insert(Vec::new(), 1);
        }
        Grade::Interior {
            charge,
            weight_bound,
        } => {
            let oscillators = Oscillators::new(model, n)?;
            let mut emit = |model: &mut Model, shape: &[u128]| {
                let count =
                    degree_difference(&oscillators, model, shape, charge as usize, primitive)?;
                if count == 0 {
                    return Ok(());
                }
                if primitive {
                    branch(model, shape, &mut |shape| {
                        let entry = components.entry(shape).or_insert(0i128);
                        *entry = entry.checked_add(count).ok_or_else(arithmetic)?;
                        guard(components.len().saturating_mul(4 * RANK + 8))
                    })
                } else {
                    components.insert(shape.iter().map(|&part| part as i128).collect(), count);
                    guard(components.len().saturating_mul(4 * RANK + 8))
                }
            };
            graded_shapes(model, weight_bound, &mut Vec::new(), &mut emit)?;
        }
    }
    let roots = if primitive {
        RootSystem::D
    } else {
        RootSystem::B
    };
    let mut total = Magnitude::default();
    let mut first = true;
    for (shape, count) in components.into_iter().rev() {
        if nonnegative(count)? == 0 {
            continue;
        }
        let absolute: Vec<_> = shape.iter().map(|part| part.unsigned_abs()).collect();
        let dimension = weyl::dimension_for(model, &absolute, RANK, roots)?;
        total
            .add_multiple(
                dimension
                    .multiply_u128(count as u128)
                    .ok_or_else(arithmetic)?,
                1,
            )
            .ok_or_else(arithmetic)?;
        let delimiter = if first { "" } else { ", " };
        first = false;
        write!(
            report,
            "{delimiter}({shape:?}, {count}, {})",
            decimal(model, dimension)?
        )
        .unwrap();
        if report.len() > MAX_SHOWN_BYTES {
            return Err(report_limit());
        }
    }
    write!(report, "] total={}", decimal(model, total)?).unwrap();
    if report.len() > MAX_SHOWN_BYTES {
        return Err(report_limit());
    }
    Ok(report)
}

fn graded_shapes(
    model: &mut Model,
    budget: u128,
    shape: &mut Vec<u128>,
    emit: &mut impl FnMut(&mut Model, &[u128]) -> Result<()>,
) -> Result<()> {
    emit(model, shape)?;
    if shape.len() == RANK {
        return Ok(());
    }
    let row = shape.len() as u128 + 1;
    let bound = budget / row;
    let max_part = shape.last().copied().unwrap_or(bound).min(bound);
    for part in 1..=max_part {
        model.spend(1)?;
        shape.push(part);
        graded_shapes(model, budget - row * part, shape, emit)?;
        shape.pop();
    }
    Ok(())
}

fn branch(
    model: &mut Model,
    shape: &[u128],
    emit: &mut impl FnMut(Vec<i128>) -> Result<()>,
) -> Result<()> {
    let mut upper = [0; RANK];
    for (part, &value) in upper.iter_mut().zip(shape) {
        *part = value as i128;
    }
    let mut limits = [(0, 0); RANK];
    for i in 0..RANK - 1 {
        limits[i] = (upper[i + 1], upper[i]);
    }
    limits[RANK - 1] = (-upper[RANK - 1], upper[RANK - 1]);
    let branches = limits
        .iter()
        .try_fold(1usize, |size, &(low, high)| {
            size.checked_mul(usize::try_from(high - low + 1).ok()?)
        })
        .ok_or_else(work_limit)?;
    guard(branches.saturating_mul(4 * RANK + 8))?;
    model.spend(branches * RANK)?;
    branch_at(&limits, &mut [0; RANK], 0, emit)
}

fn branch_at(
    limits: &[(i128, i128); RANK],
    shape: &mut [i128; RANK],
    i: usize,
    emit: &mut impl FnMut(Vec<i128>) -> Result<()>,
) -> Result<()> {
    if i == RANK {
        let mut canonical = shape.to_vec();
        trim(&mut canonical);
        return emit(canonical);
    }
    let (low, high) = limits[i];
    for value in low..=high {
        shape[i] = value;
        branch_at(limits, shape, i + 1, emit)?;
    }
    Ok(())
}

fn interlacers(
    model: &mut Model,
    lower: &[u128; RANK],
    suffix: &[u128; RANK + 1],
    budget: u128,
    upper: &mut Vec<u128>,
    emit: &mut impl FnMut(&mut Model, &[u128]) -> Result<()>,
) -> Result<()> {
    let i = upper.len();
    if i == RANK {
        return emit(model, upper);
    }
    let row = i as u128 + 1;
    let available = (budget - suffix[i + 1]) / row;
    let high = if i == 0 {
        available
    } else {
        available.min(lower[i - 1])
    };
    for part in lower[i]..=high {
        model.spend(1)?;
        upper.push(part);
        interlacers(model, lower, suffix, budget - row * part, upper, emit)?;
        upper.pop();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_harmonic_support_shortcut_matches_unpruned_weyl_extraction() {
        use super::super::partitions;

        let mut model = Model::new();
        let oscillators = Oscillators::new(&mut model, 8).unwrap_or_else(|err| panic!("{err}"));
        for charge in 0..=oscillators.n {
            for weight in 0..=oscillators.n - charge {
                partitions(
                    &mut model,
                    weight,
                    weight,
                    &mut Vec::new(),
                    &mut |model, shape| {
                        let exact = oscillators.extracted(model, shape, RANK, Some(charge))?;
                        assert_eq!(coefficient(&oscillators, model, shape, charge)?, exact);
                        Ok(())
                    },
                )
                .unwrap_or_else(|err| panic!("{err}"));
            }
        }
    }

    #[test]
    fn the_middle_exterior_power_branches_into_both_signed_constituents() {
        let mut model = Model::new();
        let mut branches = Vec::new();
        branch(&mut model, &[1; RANK], &mut |shape| {
            branches.push(shape);
            Ok(())
        })
        .unwrap_or_else(|err| panic!("{err}"));
        let mut negative = vec![1; RANK];
        negative[RANK - 1] = -1;
        assert_eq!(branches, [negative, vec![1; RANK - 1], vec![1; RANK]]);

        let branch_dimensions: Vec<_> = branches
            .iter()
            .map(|shape| {
                let absolute: Vec<_> = shape.iter().map(|part| part.unsigned_abs()).collect();
                weyl::dimension_for(&mut model, &absolute, RANK, RootSystem::D)
                    .unwrap_or_else(|err| panic!("{err}"))
                    .small()
                    .unwrap()
            })
            .collect();
        assert_eq!(branch_dimensions, [352_716, 646_646, 352_716]);
        let unrestricted = weyl::dimension_for(&mut model, &[1; RANK], RANK, RootSystem::B)
            .unwrap_or_else(|err| panic!("{err}"))
            .small()
            .unwrap();
        assert_eq!(unrestricted, branch_dimensions.iter().sum::<u128>());
    }
}
