use std::collections::VecDeque;

use crate::complex::{for_each_boundary, xor_sorted, OrderComplex, Simplex, MAX_PIVOT_CELLS};
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};

pub(crate) fn betti(model: &mut Model, ty: TypeId) -> Result<Vec<usize>> {
    let mut complex = OrderComplex::new(model, ty)?;
    let mut current = complex.vertices();
    let mut lower_rank = 0;
    let mut ranks = Vec::new();
    while !current.is_empty() {
        let next = complex.extend(model, &current)?;
        let upper_rank = boundary_rank(model, &current, &next)?;
        ranks.push(current.len() - lower_rank - upper_rank);
        lower_rank = upper_rank;
        current = next;
    }
    while ranks.last() == Some(&0) {
        ranks.pop();
    }
    Ok(ranks)
}

fn boundary_rank(model: &mut Model, rows: &[Simplex], columns: &[Simplex]) -> Result<usize> {
    if columns.is_empty() {
        return Ok(0);
    }
    let mut pivots = vec![None::<Vec<usize>>; rows.len()];
    let mut pivot_cells = 0;
    let mut rank = 0;
    for_each_boundary(model, rows, columns, |model, _, mut column| {
        column.sort_unstable();
        while let Some(&leading) = column.last() {
            if let Some(pivot) = &pivots[leading] {
                model.spend(column.len() + pivot.len())?;
                column = xor_sorted(&column, pivot);
            } else {
                pivot_cells += column.len();
                if pivot_cells > MAX_PIVOT_CELLS {
                    return Err(Error::plain("homology reduction exceeds the storage limit"));
                }
                pivots[leading] = Some(column);
                rank += 1;
                break;
            }
        }
        Ok(())
    })?;
    Ok(rank)
}

pub(crate) fn homotopy(
    model: &mut Model,
    ty: TypeId,
    left: &Value,
    right: &Value,
) -> Result<Option<usize>> {
    if model.equal(left, right)? {
        return Ok(Some(0));
    }
    if comparable(model, ty, left, right)? {
        return Ok(Some(1));
    }

    let maps = model.elements(ty)?;
    let start = model.index(ty, left)?;
    let goal = model.index(ty, right)?;
    let mut visited = vec![false; maps.len()];
    let mut queue = VecDeque::from([(start, 0)]);
    visited[start] = true;
    while let Some((current, distance)) = queue.pop_front() {
        for (next, seen) in visited.iter_mut().enumerate() {
            model.spend(1)?;
            if !*seen && comparable(model, ty, &maps[current], &maps[next])? {
                if next == goal {
                    return Ok(Some(distance + 1));
                }
                *seen = true;
                queue.push_back((next, distance + 1));
            }
        }
    }
    Ok(None)
}

fn comparable(model: &mut Model, ty: TypeId, left: &Value, right: &Value) -> Result<bool> {
    Ok(model.leq(ty, left, right)? || model.leq(ty, right, left)?)
}
