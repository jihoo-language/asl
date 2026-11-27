use std::collections::{BTreeMap, VecDeque};
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::{add_scaled, Basis, Coordinate};

struct FlagStep {
    space: Basis,
    image_rank: usize,
}

#[derive(Default)]
struct Flag {
    steps: BTreeMap<usize, FlagStep>,
    pending: VecDeque<usize>,
    cells: usize,
}

impl Flag {
    fn include(
        &mut self,
        model: &mut Model,
        prime: Prime,
        space: Basis,
        retained: usize,
    ) -> Result<()> {
        let rank = space.rank();
        if let Some(existing) = self.steps.get(&rank) {
            return contained(model, prime, &space, &existing.space);
        }
        self.cells += 6 + space.cells();
        check_storage(retained + self.cells)?;
        self.steps.insert(
            rank,
            FlagStep {
                space,
                image_rank: 0,
            },
        );
        self.pending.push_back(rank);
        Ok(())
    }
}

type Column = Rc<[Coordinate]>;

pub(crate) fn final_type(
    model: &mut Model,
    prime: Prime,
    dimension: usize,
    retained: impl Fn() -> usize,
    mut frobenius: impl FnMut(&mut Model, usize) -> Result<Column>,
    mut verschiebung: impl FnMut(&mut Model, usize) -> Result<Column>,
) -> Result<Vec<usize>> {
    check_storage(dimension.saturating_mul(4))?;
    let mut flag = Flag::default();
    flag.include(model, prime, Basis::default(), retained())?;
    let mut ambient = Basis::default();
    for i in 0..dimension {
        ambient.insert(model, prime, vec![(i, 1)])?;
        check_storage(retained() + flag.cells + ambient.cells())?;
    }
    flag.include(model, prime, ambient, retained())?;

    while let Some(rank) = flag.pending.pop_front() {
        let step = flag.steps.get_mut(&rank).unwrap();
        let mut image = Basis::default();
        let mut quotient = Basis::default();
        for vector in step.space.vectors(prime) {
            let mut column = Vec::new();
            for (index, value) in vector {
                let mapped = verschiebung(model, index)?;
                column = add_scaled(model, prime, &column, &mapped, value)?;
            }
            image.insert(model, prime, column)?;
            check_storage(retained() + flag.cells + image.cells() + quotient.cells())?;
        }
        for vector in step.space.vectors(prime) {
            quotient.insert(model, prime, vector)?;
            check_storage(retained() + flag.cells + image.cells() + quotient.cells())?;
        }
        let mut preimage = Basis::default();
        // Dependencies of F-columns modulo the target give its preimage.
        for i in 0..dimension {
            let column = frobenius(model, i)?;
            if let Some(vector) =
                quotient.insert_with_witness(model, prime, column.to_vec(), vec![(i, 1)])?
            {
                preimage.insert(model, prime, vector)?;
            }
            check_storage(
                retained() + flag.cells + image.cells() + quotient.cells() + preimage.cells(),
            )?;
        }
        drop(quotient);
        step.image_rank = image.rank();
        flag.include(model, prime, image, retained() + preimage.cells())?;
        flag.include(model, prime, preimage, retained())?;
    }

    let genus = dimension / 2;
    let mut previous = &flag.steps[&0];
    let mut previous_rank = 0;
    let mut sequence = Vec::with_capacity(genus);
    for (&rank, step) in flag.steps.iter().skip(1) {
        contained(model, prime, &previous.space, &step.space)?;
        let width = rank - previous_rank;
        let rise = step.image_rank.checked_sub(previous.image_rank);
        if rise != Some(0) && rise != Some(width) {
            return Err(Error::plain("EO flag has invalid Verschiebung ranks"));
        }
        for index in previous_rank + 1..=rank.min(genus) {
            let increment = if rise == Some(width) {
                index - previous_rank
            } else {
                0
            };
            sequence.push(previous.image_rank + increment);
        }
        previous = step;
        previous_rank = rank;
    }
    Ok(sequence)
}

fn contained(model: &mut Model, prime: Prime, lower: &Basis, upper: &Basis) -> Result<()> {
    for vector in lower.vectors(prime) {
        if !upper.reduce(model, prime, vector)?.is_empty() {
            return Err(Error::plain("canonical EO subspaces are not nested"));
        }
    }
    Ok(())
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain("canonical EO flag exceeds the storage limit"));
    }
    Ok(())
}
