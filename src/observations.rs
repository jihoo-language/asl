use std::collections::{HashMap, VecDeque};

use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};

const MAX_SIMPLICES: usize = 50_000;
const MAX_SIMPLEX_VERTICES: usize = 300_000;
const MAX_PIVOT_CELLS: usize = 300_000;

type Simplex = Vec<usize>;

pub(crate) fn betti(model: &mut Model, ty: TypeId) -> Result<Vec<usize>> {
    let elements = model.elements(ty)?;
    if elements.is_empty() {
        return Ok(Vec::new());
    }
    let mut upper = vec![Vec::new(); elements.len()];
    let mut edges = 0;
    for (left, x) in elements.iter().enumerate() {
        for (right, y) in elements.iter().enumerate().skip(left + 1) {
            let comparison = if model.leq(ty, x, y)? {
                Some((left, right))
            } else if model.leq(ty, y, x)? {
                Some((right, left))
            } else {
                None
            };
            if let Some((smaller, larger)) = comparison {
                edges += 1;
                if edges + elements.len() > MAX_SIMPLICES {
                    return Err(Error::plain("order complex exceeds the simplex limit"));
                }
                upper[smaller].push(larger);
            }
        }
    }

    let mut current: Vec<Simplex> = (0..elements.len()).map(|index| vec![index]).collect();
    let mut simplices = current.len();
    let mut simplex_vertices = simplices;
    let mut lower_rank = 0;
    let mut ranks = Vec::new();
    loop {
        let mut next = Vec::new();
        for simplex in &current {
            for &larger in &upper[*simplex.last().unwrap()] {
                simplices += 1;
                simplex_vertices += simplex.len() + 1;
                if simplices > MAX_SIMPLICES || simplex_vertices > MAX_SIMPLEX_VERTICES {
                    return Err(Error::plain("order complex exceeds the simplex limit"));
                }
                model.spend(simplex.len() + 1)?;
                let mut extended = simplex.clone();
                extended.push(larger);
                next.push(extended);
            }
        }
        let upper_rank = boundary_rank(model, &current, &next)?;
        ranks.push(current.len() - lower_rank - upper_rank);
        if next.is_empty() {
            break;
        }
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
    let indices: HashMap<&[usize], usize> = rows
        .iter()
        .enumerate()
        .map(|(index, row)| (row.as_slice(), index))
        .collect();
    let mut pivots = vec![None::<Vec<usize>>; rows.len()];
    let mut pivot_cells = 0;
    let mut rank = 0;
    let mut face = Vec::new();
    for simplex in columns {
        let mut column = Vec::with_capacity(simplex.len());
        for vertex in 0..simplex.len() {
            model.spend(simplex.len())?;
            face.clear();
            face.extend_from_slice(&simplex[..vertex]);
            face.extend_from_slice(&simplex[vertex + 1..]);
            column.push(indices[face.as_slice()]);
        }
        column.sort_unstable();

        // Over F2, adding columns is the symmetric difference of row indices.
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
    }
    Ok(rank)
}

fn xor_sorted(mut left: &[usize], mut right: &[usize]) -> Vec<usize> {
    let mut result = Vec::with_capacity(left.len() + right.len());
    while let (Some((&x, xs)), Some((&y, ys))) = (left.split_first(), right.split_first()) {
        match x.cmp(&y) {
            std::cmp::Ordering::Less => {
                result.push(x);
                left = xs;
            }
            std::cmp::Ordering::Greater => {
                result.push(y);
                right = ys;
            }
            std::cmp::Ordering::Equal => {
                left = xs;
                right = ys;
            }
        }
    }
    result.extend_from_slice(left);
    result.extend_from_slice(right);
    result
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
