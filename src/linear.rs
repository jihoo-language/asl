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
        mut image: Vec<usize>,
        mut witness: Vec<usize>,
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

pub(crate) fn xor(model: &mut Model, left: &[usize], right: &[usize]) -> Result<Vec<usize>> {
    model.spend(left.len() + right.len())?;
    Ok(xor_sorted(left, right))
}
