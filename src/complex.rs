use std::cmp::Ordering;
use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::model::{Model, TypeId};

const MAX_SIMPLICES: usize = 50_000;
const MAX_SIMPLEX_VERTICES: usize = 300_000;
pub(crate) const MAX_PIVOT_CELLS: usize = 300_000;

pub(crate) type Simplex = Vec<usize>;

pub(crate) struct OrderComplex {
    upper: Vec<Vec<usize>>,
    simplices: usize,
    simplex_vertices: usize,
}

impl OrderComplex {
    pub(crate) fn new(model: &mut Model, ty: TypeId) -> Result<Self> {
        let elements = model.elements(ty)?;
        let mut complex = Self {
            upper: vec![Vec::new(); elements.len()],
            simplices: elements.len(),
            simplex_vertices: elements.len(),
        };
        if model.intrinsically_discrete(ty) {
            return Ok(complex);
        }
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
                    complex.upper[smaller].push(larger);
                }
            }
        }
        Ok(complex)
    }

    pub(crate) fn vertices(&self) -> Vec<Simplex> {
        (0..self.upper.len()).map(|index| vec![index]).collect()
    }

    pub(crate) fn extend(&mut self, model: &mut Model, layer: &[Simplex]) -> Result<Vec<Simplex>> {
        let mut next = Vec::new();
        for simplex in layer {
            for &larger in &self.upper[*simplex.last().unwrap()] {
                self.simplices += 1;
                self.simplex_vertices += simplex.len() + 1;
                if self.simplices > MAX_SIMPLICES || self.simplex_vertices > MAX_SIMPLEX_VERTICES {
                    return Err(Error::plain("order complex exceeds the simplex limit"));
                }
                model.spend(simplex.len() + 1)?;
                let mut extended = simplex.clone();
                extended.push(larger);
                next.push(extended);
            }
        }
        Ok(next)
    }
}

pub(crate) fn for_each_boundary(
    model: &mut Model,
    rows: &[Simplex],
    columns: &[Simplex],
    mut use_boundary: impl FnMut(&mut Model, usize, Vec<usize>) -> Result<()>,
) -> Result<()> {
    if columns.is_empty() {
        return Ok(());
    }
    let indices: HashMap<&[usize], usize> = rows
        .iter()
        .enumerate()
        .map(|(index, row)| (row.as_slice(), index))
        .collect();
    let mut face = Vec::new();
    for (index, simplex) in columns.iter().enumerate() {
        let mut boundary = Vec::with_capacity(simplex.len());
        for vertex in 0..simplex.len() {
            model.spend(simplex.len())?;
            face.clear();
            face.extend_from_slice(&simplex[..vertex]);
            face.extend_from_slice(&simplex[vertex + 1..]);
            boundary.push(indices[face.as_slice()]);
        }
        use_boundary(model, index, boundary)?;
    }
    Ok(())
}

// Sparse F2 vectors are sorted lists of coordinates with coefficient one.
pub(crate) fn xor_sorted<T: Ord + Clone>(mut left: &[T], mut right: &[T]) -> Vec<T> {
    let mut result = Vec::with_capacity(left.len() + right.len());
    while let (Some((x, xs)), Some((y, ys))) = (left.split_first(), right.split_first()) {
        match x.cmp(y) {
            Ordering::Less => {
                result.push(x.clone());
                left = xs;
            }
            Ordering::Greater => {
                result.push(y.clone());
                right = ys;
            }
            Ordering::Equal => {
                left = xs;
                right = ys;
            }
        }
    }
    result.extend_from_slice(left);
    result.extend_from_slice(right);
    result
}
