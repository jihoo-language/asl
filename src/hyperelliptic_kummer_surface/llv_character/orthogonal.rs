use crate::error::Result;
use crate::model::Model;

use super::work_limit;

pub(super) const DIMENSION: usize = 23;
pub(super) const RANK: usize = DIMENSION / 2;

// The SO(N) -> SO(N-1) interlacings have the same nonnegative intervals
// in absolute coordinates. The only two-valued coordinate is the last
// one for an even-dimensional child, including the SO(2) weight.
pub(super) fn branch(
    model: &mut Model,
    dimension: usize,
    shape: &[u128],
    max_nonzero: usize,
    emit: &mut impl FnMut(&mut Model, &[u128], u32) -> Result<()>,
) -> Result<()> {
    let rank = (dimension - 1) / 2;
    let signed = !dimension.is_multiple_of(2);
    let mut limits = [(0, 0); RANK];
    let mut branches = 1usize;
    for (i, range) in limits[..rank].iter_mut().enumerate() {
        let low = shape.get(i + 1).copied().unwrap_or(0);
        let high = if i < max_nonzero {
            shape.get(i).copied().unwrap_or(0)
        } else {
            0
        };
        if low > high {
            return Ok(());
        }
        let width = usize::try_from(high - low)
            .ok()
            .and_then(|count| count.checked_add(1))
            .ok_or_else(work_limit)?;
        branches = branches.checked_mul(width).ok_or_else(work_limit)?;
        *range = (low, high);
    }
    model
        .spend(branches.saturating_mul(rank))
        .map_err(|_| work_limit())?;
    branch_at(model, &limits[..rank], &mut [0; RANK], 0, signed, emit)
}

fn branch_at(
    model: &mut Model,
    limits: &[(u128, u128)],
    shape: &mut [u128; RANK],
    i: usize,
    signed: bool,
    emit: &mut impl FnMut(&mut Model, &[u128], u32) -> Result<()>,
) -> Result<()> {
    let Some(&(low, high)) = limits.get(i) else {
        let length = shape[..i]
            .iter()
            .rposition(|&part| part != 0)
            .map_or(0, |index| index + 1);
        let signs = if signed && length == i { 2 } else { 1 };
        return emit(model, &shape[..length], signs);
    };
    for value in low..=high {
        shape[i] = value;
        branch_at(model, limits, shape, i + 1, signed, emit)?;
    }
    Ok(())
}
