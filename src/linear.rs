use std::collections::BTreeMap;

use crate::complex::{xor_sorted, MAX_PIVOT_CELLS};
use crate::error::{Error, Result};
use crate::model::Model;

// Supports are sorted, duplicate-free lists of nonzero coordinates over F2.
// Each row has its pivot in the last coordinate of its support.
#[derive(Default)]
pub(crate) struct Basis {
    pivots: BTreeMap<usize, Row>,
    cells: usize,
}

struct Row {
    image: Box<[usize]>,
    witness: Box<[usize]>,
}

impl Basis {
    pub(crate) fn cells(&self) -> usize {
        self.cells
    }

    pub(crate) fn contains_pivot(&self, coordinate: usize) -> bool {
        self.pivots.contains_key(&coordinate)
    }

    pub(crate) fn rank(&self) -> usize {
        self.pivots.len()
    }

    pub(crate) fn vectors(&self) -> impl Iterator<Item = &[usize]> {
        self.pivots.values().map(|row| row.image.as_ref())
    }

    pub(crate) fn vector(&self, coordinate: usize) -> Option<&[usize]> {
        self.pivots.get(&coordinate).map(|row| row.image.as_ref())
    }

    pub(crate) fn insert(&mut self, model: &mut Model, image: Vec<usize>) -> Result<bool> {
        Ok(self
            .insert_with_witness(model, image, Vec::new())?
            .is_none())
    }

    // A dependent column returns its accumulated witness for the zero image.
    pub(crate) fn insert_with_witness(
        &mut self,
        model: &mut Model,
        image: Vec<usize>,
        witness: Vec<usize>,
    ) -> Result<Option<Vec<usize>>> {
        self.insert_traced(model, image, witness, |_, _, _| Ok(()))
    }

    pub(crate) fn insert_traced(
        &mut self,
        model: &mut Model,
        mut image: Vec<usize>,
        mut witness: Vec<usize>,
        trace: impl FnOnce(&mut Model, &[usize], &[usize]) -> Result<()>,
    ) -> Result<Option<Vec<usize>>> {
        model.spend(1 + image.len() + witness.len())?;
        while let Some(&leading) = image.last() {
            if let Some(row) = self.pivots.get(&leading) {
                image = xor(model, &image, &row.image)?;
                witness = xor(model, &witness, &row.witness)?;
            } else {
                let cells = self.cells + 1 + image.len() + witness.len();
                if cells > MAX_PIVOT_CELLS {
                    return Err(Error::plain("linear reduction exceeds the storage limit"));
                }
                trace(model, &image, &witness)?;
                self.cells = cells;
                self.pivots.insert(
                    leading,
                    Row {
                        image: image.into_boxed_slice(),
                        witness: witness.into_boxed_slice(),
                    },
                );
                return Ok(None);
            }
        }
        trace(model, &image, &witness)?;
        Ok(Some(witness))
    }

    pub(crate) fn reduce(&self, model: &mut Model, image: Vec<usize>) -> Result<Vec<usize>> {
        self.remainder(model, image, None)
    }

    pub(crate) fn reduce_with_witness(
        &self,
        model: &mut Model,
        image: Vec<usize>,
        mut witness: Vec<usize>,
    ) -> Result<(Vec<usize>, Vec<usize>)> {
        let image = self.remainder(model, image, Some(&mut witness))?;
        Ok((image, witness))
    }

    fn remainder(
        &self,
        model: &mut Model,
        mut image: Vec<usize>,
        mut witness: Option<&mut Vec<usize>>,
    ) -> Result<Vec<usize>> {
        let mut remainder = Vec::new();
        while let Some(&leading) = image.last() {
            if let Some(row) = self.pivots.get(&leading) {
                image = xor(model, &image, &row.image)?;
                if let Some(witness) = witness.as_mut() {
                    **witness = xor(model, witness, &row.witness)?;
                }
            } else {
                model.spend(1)?;
                remainder.push(image.pop().unwrap());
            }
        }
        remainder.reverse();
        Ok(remainder)
    }
}

#[derive(Default)]
pub(crate) struct QuotientBasis {
    boundaries: Basis,
    cycles: Basis,
    generators: Vec<usize>,
}

impl QuotientBasis {
    pub(crate) fn build(
        model: &mut Model,
        boundaries: Basis,
        size: usize,
        mut differential: impl FnMut(&mut Model, usize) -> Result<Vec<usize>>,
    ) -> Result<Self> {
        let mut images = Basis::default();
        let mut cycles = Basis::default();
        let mut generators = Vec::new();
        for coordinate in 0..size {
            if boundaries.contains_pivot(coordinate) {
                continue;
            }
            // Nonpivot coordinates complement the boundaries. A dependence
            // among their images gives a cycle with this coordinate as pivot.
            let image = differential(model, coordinate)?;
            if let Some(cycle) = images.insert_with_witness(model, image, vec![coordinate])? {
                let dependence =
                    cycles.insert_with_witness(model, cycle, vec![generators.len()])?;
                debug_assert!(dependence.is_none());
                generators.push(coordinate);
            }
        }
        model.allocate(boundaries.cells() + cycles.cells() + generators.len())?;
        Ok(Self {
            boundaries,
            cycles,
            generators,
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.generators.len()
    }

    // A cycle span in the complement of the boundary pivots gives the
    // same ordered generators as build. Each resulting representative
    // has zero entries at the other cycle pivots.
    pub(crate) fn from_cycle_span(
        model: &mut Model,
        boundaries: Basis,
        span: Basis,
    ) -> Result<Self> {
        let mut cycles = Basis::default();
        let mut generators = Vec::new();
        for row in span.vectors() {
            let row = cycles.reduce(model, row.to_vec())?;
            let pivot = *row.last().unwrap();
            cycles.insert_with_witness(model, row, vec![generators.len()])?;
            generators.push(pivot);
        }
        model.allocate(boundaries.cells() + cycles.cells() + generators.len())?;
        Ok(Self {
            boundaries,
            cycles,
            generators,
        })
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        support: Vec<usize>,
    ) -> Result<Option<Vec<usize>>> {
        let support = self.boundaries.reduce(model, support)?;
        let (remainder, coordinates) =
            self.cycles
                .reduce_with_witness(model, support, Vec::new())?;
        Ok(remainder.is_empty().then_some(coordinates))
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<usize>> {
        let mut support = Vec::new();
        for &index in class {
            let pivot = self.generators[index];
            support = xor(model, &support, self.cycles.vector(pivot).unwrap())?;
        }
        model.spend(support.len())?;
        Ok(support)
    }

    pub(crate) fn is_boundary(&self, model: &mut Model, support: Vec<usize>) -> Result<bool> {
        Ok(self.boundaries.reduce(model, support)?.is_empty())
    }
}

pub(crate) fn xor(model: &mut Model, left: &[usize], right: &[usize]) -> Result<Vec<usize>> {
    model.spend(left.len() + right.len())?;
    Ok(xor_sorted(left, right))
}
