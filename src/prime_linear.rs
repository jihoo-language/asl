use std::cmp::Ordering;
use std::collections::BTreeMap;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;

pub(crate) type Coordinate = (usize, u32);

// Vectors contain only nonzero coefficients and are sorted by coordinate.
// The last coordinate of every stored row is a pivot with coefficient one.
#[derive(Default)]
pub(crate) struct Basis {
    rows: BTreeMap<usize, Row>,
    cells: usize,
}

struct Row {
    image: Box<[Coordinate]>,
    witness: Box<[Coordinate]>,
}

impl Basis {
    pub(crate) fn cells(&self) -> usize {
        self.cells
    }

    pub(crate) fn insert(
        &mut self,
        model: &mut Model,
        prime: Prime,
        image: Vec<Coordinate>,
    ) -> Result<()> {
        self.insert_with_witness(model, prime, image, Vec::new())?;
        Ok(())
    }

    // An independent column is normalized before storage. A dependent
    // column returns the accumulated witness whose image is zero.
    pub(crate) fn insert_with_witness(
        &mut self,
        model: &mut Model,
        prime: Prime,
        mut image: Vec<Coordinate>,
        mut witness: Vec<Coordinate>,
    ) -> Result<Option<Vec<Coordinate>>> {
        model.spend(1 + image.len() + witness.len())?;
        while let Some(&(leading, coefficient)) = image.last() {
            if let Some(row) = self.rows.get(&leading) {
                let factor = prime.negate(coefficient);
                image = add_scaled(model, prime, &image, &row.image, factor)?;
                witness = add_scaled(model, prime, &witness, &row.witness, factor)?;
            } else {
                let cells = self.cells + 1 + 2 * (image.len() + witness.len());
                check_storage(cells)?;
                if coefficient != 1 {
                    let inverse = invert(model, prime, coefficient)?;
                    model.spend(image.len() + witness.len())?;
                    for (_, value) in image.iter_mut().chain(&mut witness) {
                        *value = prime.mul(*value, inverse);
                    }
                }
                self.cells = cells;
                self.rows.insert(
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

    pub(crate) fn reduce(
        &self,
        model: &mut Model,
        prime: Prime,
        image: Vec<Coordinate>,
    ) -> Result<Vec<Coordinate>> {
        self.remainder(model, prime, image, None)
    }

    pub(crate) fn reduce_with_witness(
        &self,
        model: &mut Model,
        prime: Prime,
        image: Vec<Coordinate>,
    ) -> Result<(Vec<Coordinate>, Vec<Coordinate>)> {
        let mut witness = Vec::new();
        let remainder = self.remainder(model, prime, image, Some(&mut witness))?;
        Ok((remainder, witness))
    }

    fn remainder(
        &self,
        model: &mut Model,
        prime: Prime,
        mut image: Vec<Coordinate>,
        mut witness: Option<&mut Vec<Coordinate>>,
    ) -> Result<Vec<Coordinate>> {
        let mut remainder = Vec::new();
        while let Some((leading, coefficient)) = image.pop() {
            if let Some(row) = self.rows.get(&leading) {
                image.push((leading, coefficient));
                image = add_scaled(model, prime, &image, &row.image, prime.negate(coefficient))?;
                if let Some(witness) = witness.as_mut() {
                    // The original vector equals the remainder plus the
                    // images of these positively accumulated coefficients.
                    **witness = add_scaled(model, prime, witness, &row.witness, coefficient)?;
                }
            } else {
                model.spend(1)?;
                remainder.push((leading, coefficient));
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
        prime: Prime,
        boundaries: Basis,
        size: usize,
        mut differential: impl FnMut(&mut Model, usize) -> Result<Vec<Coordinate>>,
    ) -> Result<Self> {
        let mut images = Basis::default();
        let mut cycles = Basis::default();
        let mut generators = Vec::new();
        for coordinate in 0..size {
            if boundaries.rows.contains_key(&coordinate) {
                continue;
            }
            // Nonpivot coordinates form a complement to the boundaries.
            // A dependence gives a cycle led by the new coordinate.
            let image = differential(model, coordinate)?;
            if let Some(cycle) =
                images.insert_with_witness(model, prime, image, vec![(coordinate, 1)])?
            {
                let result =
                    cycles.insert_with_witness(model, prime, cycle, vec![(generators.len(), 1)])?;
                debug_assert!(result.is_none());
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

    pub(crate) fn project(
        &self,
        model: &mut Model,
        prime: Prime,
        support: Vec<Coordinate>,
    ) -> Result<Option<Vec<Coordinate>>> {
        let support = self.boundaries.reduce(model, prime, support)?;
        let (remainder, coordinates) = self.cycles.reduce_with_witness(model, prime, support)?;
        Ok(remainder.is_empty().then_some(coordinates))
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        prime: Prime,
        class: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let mut support = Vec::new();
        for &(index, coefficient) in class {
            let pivot = self.generators[index];
            let row = &self.cycles.rows[&pivot];
            support = add_scaled(model, prime, &support, &row.image, coefficient)?;
        }
        model.spend(support.len())?;
        Ok(support)
    }

    pub(crate) fn is_boundary(
        &self,
        model: &mut Model,
        prime: Prime,
        support: Vec<Coordinate>,
    ) -> Result<bool> {
        Ok(self.boundaries.reduce(model, prime, support)?.is_empty())
    }
}

pub(crate) fn add_scaled(
    model: &mut Model,
    prime: Prime,
    left: &[Coordinate],
    right: &[Coordinate],
    factor: u32,
) -> Result<Vec<Coordinate>> {
    if factor == 0 {
        model.spend(left.len())?;
        return Ok(left.to_vec());
    }
    let size = left.len() + right.len();
    model.spend(size)?;
    let mut output = Vec::with_capacity(size.min(MAX_PIVOT_CELLS / 2));
    let (mut i, mut j) = (0, 0);
    while let (Some(&(a, x)), Some(&(b, y))) = (left.get(i), right.get(j)) {
        match a.cmp(&b) {
            Ordering::Less => {
                output.push((a, x));
                i += 1;
            }
            Ordering::Equal => {
                let coefficient = prime.add(x, prime.mul(factor, y));
                if coefficient != 0 {
                    output.push((a, coefficient));
                }
                i += 1;
                j += 1;
            }
            Ordering::Greater => {
                output.push((b, prime.mul(factor, y)));
                j += 1;
            }
        }
    }
    output.extend_from_slice(&left[i..]);
    output.extend(
        right[j..]
            .iter()
            .map(|&(index, value)| (index, prime.mul(value, factor))),
    );
    check_storage(2 * output.len())?;
    Ok(output)
}

// Fermat inversion uses at most two field multiplications per exponent bit.
fn invert(model: &mut Model, prime: Prime, mut value: u32) -> Result<u32> {
    let mut power = prime.get() - 2;
    let mut inverse = 1;
    while power != 0 {
        model.spend(2)?;
        if power & 1 != 0 {
            inverse = prime.mul(inverse, value);
        }
        value = prime.mul(value, value);
        power >>= 1;
    }
    Ok(inverse)
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "odd-prime linear reduction exceeds the storage limit",
        ));
    }
    Ok(())
}
