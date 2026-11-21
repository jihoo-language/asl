use std::cmp::Ordering;
use std::collections::BTreeMap;

use crate::binary_field::BinaryField;
use crate::error::{Error, Result};
use crate::model::Model;

const MAX_CELLS: usize = 300_000;

pub(crate) type Entry = (usize, usize);
pub(crate) type Vector = Vec<Entry>;

#[derive(Default)]
pub(crate) struct Basis {
    pivots: BTreeMap<usize, Pivot>,
    cells: usize,
}

struct Pivot {
    image: Box<[Entry]>,
    witness: Box<[Entry]>,
}

impl Basis {
    pub(crate) fn cells(&self) -> usize {
        self.cells
    }

    pub(crate) fn rank(&self) -> usize {
        self.pivots.len()
    }

    pub(crate) fn insert(
        &mut self,
        model: &mut Model,
        field: BinaryField,
        image: Vector,
    ) -> Result<()> {
        self.insert_with_witness(model, field, image, Vec::new())?;
        Ok(())
    }

    // Each stored pivot has coefficient one. A dependent column returns its
    // tracked linear combination of source coordinates.
    pub(crate) fn insert_with_witness(
        &mut self,
        model: &mut Model,
        field: BinaryField,
        mut image: Vector,
        mut witness: Vector,
    ) -> Result<Option<Vector>> {
        model.spend(1 + image.len() + witness.len())?;
        while let Some(&(leading, coefficient)) = image.last() {
            if let Some(pivot) = self.pivots.get(&leading) {
                image = add_multiple(model, field, &image, &pivot.image, coefficient)?;
                witness = add_multiple(model, field, &witness, &pivot.witness, coefficient)?;
            } else {
                let cells = self.cells + 1 + 2 * (image.len() + witness.len());
                check_storage(cells)?;
                if coefficient != 1 {
                    let inverse = field.inverse(model, coefficient)?;
                    model.spend(field.degree() * (image.len() + witness.len()))?;
                    for (_, value) in image.iter_mut().chain(&mut witness) {
                        *value = field.multiply(*value, inverse);
                    }
                }
                self.cells = cells;
                self.pivots.insert(
                    leading,
                    Pivot {
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
        field: BinaryField,
        mut image: Vector,
        mut witness: Option<&mut Vector>,
    ) -> Result<Vector> {
        let mut remainder = Vec::new();
        while let Some((leading, coefficient)) = image.pop() {
            if let Some(pivot) = self.pivots.get(&leading) {
                image.push((leading, coefficient));
                image = add_multiple(model, field, &image, &pivot.image, coefficient)?;
                if let Some(witness) = &mut witness {
                    **witness = add_multiple(model, field, witness, &pivot.witness, coefficient)?;
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

pub(crate) struct Quotient {
    boundaries: Basis,
    cycles: Basis,
    generators: Vec<usize>,
}

impl Quotient {
    pub(crate) fn build(
        model: &mut Model,
        field: BinaryField,
        boundaries: Basis,
        size: usize,
        mut differential: impl FnMut(&mut Model, usize) -> Result<Vector>,
    ) -> Result<Self> {
        let mut images = Basis::default();
        let mut cycles = Basis::default();
        let mut generators = Vec::new();
        // Skipping the leading boundary coordinates chooses a complement to
        // the incoming image. The kernel on this complement lifts homology.
        for coordinate in 0..size {
            if boundaries.pivots.contains_key(&coordinate) {
                continue;
            }
            let image = differential(model, coordinate)?;
            if let Some(cycle) =
                images.insert_with_witness(model, field, image, vec![(coordinate, 1)])?
            {
                let reduced =
                    cycles.insert_with_witness(model, field, cycle, vec![(generators.len(), 1)])?;
                debug_assert!(reduced.is_none());
                generators.push(coordinate);
            }
        }
        model.allocate(boundaries.cells + cycles.cells + generators.len())?;
        Ok(Self {
            boundaries,
            cycles,
            generators,
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.generators.len()
    }

    pub(crate) fn representative(&self, generator: usize) -> &[Entry] {
        &self.cycles.pivots[&self.generators[generator]].image
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        field: BinaryField,
        vector: Vector,
    ) -> Result<Vector> {
        let reduced = self.boundaries.reduce(model, field, vector, None)?;
        let mut coordinates = Vec::new();
        let remainder = self
            .cycles
            .reduce(model, field, reduced, Some(&mut coordinates))?;
        if !remainder.is_empty() {
            return Err(Error::plain("specialized image is not a Margolis cycle"));
        }
        Ok(coordinates)
    }
}

pub(crate) fn add_multiple(
    model: &mut Model,
    field: BinaryField,
    mut left: &[Entry],
    mut right: &[Entry],
    factor: usize,
) -> Result<Vector> {
    if factor == 0 {
        model.spend(left.len())?;
        return Ok(left.to_vec());
    }
    model.spend(left.len() + right.len() * if factor == 1 { 1 } else { field.degree() })?;
    let mut result = Vec::new();
    let scale = |value| {
        if factor == 1 {
            value
        } else {
            field.multiply(factor, value)
        }
    };
    while let (Some((&(a, x), rest_left)), Some((&(b, y), rest_right))) =
        (left.split_first(), right.split_first())
    {
        match a.cmp(&b) {
            Ordering::Less => {
                result.push((a, x));
                left = rest_left;
            }
            Ordering::Equal => {
                let coefficient = x ^ scale(y);
                if coefficient != 0 {
                    result.push((a, coefficient));
                }
                left = rest_left;
                right = rest_right;
            }
            Ordering::Greater => {
                result.push((b, scale(y)));
                right = rest_right;
            }
        }
    }
    result.extend_from_slice(left);
    result.extend(right.iter().map(|&(row, value)| (row, scale(value))));
    check_storage(2 * result.len())?;
    Ok(result)
}

pub(crate) fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain(
            "binary-field reduction exceeds the storage limit",
        ));
    }
    Ok(())
}
