use crate::error::Result;
use crate::model::Model;

use super::super::{arithmetic, guard};

pub(super) fn horizontal(outer: &[usize], inner: &[usize]) -> bool {
    outer
        .iter()
        .skip(1)
        .enumerate()
        .all(|(row, &part)| inner.get(row).copied().unwrap_or(0) >= part)
}

// The skew diagram is read from right to left, one row at a time.
// The two neighboring assignments are already known when a cell is visited.
pub(super) fn coefficient(
    model: &mut Model,
    outer: &[usize],
    inner: &[usize],
    content: &[usize],
) -> Result<u32> {
    model.spend(1)?;
    if inner.is_empty() {
        return Ok(u32::from(outer == content));
    }
    if content.is_empty() {
        return Ok(u32::from(outer == inner));
    }
    if content.len() == 1 {
        return Ok(u32::from(horizontal(outer, inner)));
    }
    if inner.len() == 1 {
        return Ok(u32::from(horizontal(outer, content)));
    }

    let boxes = outer.iter().sum::<usize>();
    guard(16 * boxes)?;
    model.spend(boxes)?;
    let mut neighbors = Vec::new();
    let mut upper = Vec::new();
    for (i, &width) in outer.iter().enumerate() {
        let start = inner.get(i).copied().unwrap_or(0);
        let mut row = vec![None; width];
        for j in (start..width).rev() {
            let right = row.get(j + 1).copied().flatten();
            let above = upper.get(j).copied().flatten();
            row[j] = Some(neighbors.len());
            neighbors.push((right, above));
        }
        upper = row;
    }
    count_tableaux(model, &neighbors, content)
}

fn count_tableaux(
    model: &mut Model,
    neighbors: &[(Option<usize>, Option<usize>)],
    content: &[usize],
) -> Result<u32> {
    let mut assigned = vec![0usize; neighbors.len()];
    let mut next = assigned.clone();
    let mut used = vec![0; content.len()];
    let mut count = 0u32;
    let mut depth = 0;
    loop {
        let (right, above) = neighbors[depth];
        let low = above.map_or(0, |i| assigned[i] + 1).max(next[depth]);
        let high = right.map_or(content.len() - 1, |i| assigned[i]);
        model.spend(high.saturating_add(2).saturating_sub(low).max(1))?;
        let choice = (low..=high).find(|&value| {
            used[value] < content[value] && (value == 0 || used[value] < used[value - 1])
        });
        if let Some(value) = choice {
            next[depth] = value + 1;
            if depth + 1 == neighbors.len() {
                count = count.checked_add(1).ok_or_else(arithmetic)?;
            } else {
                assigned[depth] = value;
                used[value] += 1;
                depth += 1;
                next[depth] = 0;
            }
        } else if depth == 0 {
            return Ok(count);
        } else {
            depth -= 1;
            used[assigned[depth]] -= 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn elementary_pieri_and_nontrivial_lattice_tableaux_have_exact_multiplicities() {
        let mut model = Model::new();
        let mut lr = |outer: &[usize], inner: &[usize], content: &[usize]| {
            coefficient(&mut model, outer, inner, content).unwrap_or_else(|err| panic!("{err}"))
        };
        assert_eq!(lr(&[3, 2, 1], &[2, 1], &[2, 1]), 2);
        assert_eq!(lr(&[4, 1, 1], &[2, 1], &[2, 1]), 1);
        assert_eq!(lr(&[2, 2, 2], &[2, 1], &[2, 1]), 1);
        assert_eq!(lr(&[3, 3], &[2, 1], &[2, 1]), 1);
        assert_eq!(lr(&[3, 2, 1], &[2, 2], &[2]), 1);
        assert_eq!(lr(&[2, 2, 2], &[2, 2, 1], &[1]), 1);
    }
}
