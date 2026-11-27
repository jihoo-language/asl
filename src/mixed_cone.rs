use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::mixed::{Complex, Differential, Grade, Id};
use crate::model::Model;
use crate::prime_linear::{add_scaled, Coordinate};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Include,
    Lift,
    Project,
}

pub(crate) struct Cone {
    source: Id,
    target: Id,
    target_len: usize,
    map: Rc<[Coordinate]>,
}

impl Cone {
    pub(crate) fn build(
        model: &mut Model,
        map_id: Id,
        map: Rc<[Coordinate]>,
    ) -> Result<(Self, Complex)> {
        let hom = model.mixed_hom(map_id)?;
        let maps = model.mixed_complex(map_id);
        if !maps.apply(model, Differential::Total, &map)?.is_empty() {
            return Err(Error::plain("mixcone requires a closed degree-zero map"));
        }
        let source = model.mixed_complex(hom.source);
        let target = model.mixed_complex(hom.target);
        let dimension = source.len() + target.len();
        let mut cells = 8 * dimension + 2 * map.len();
        check_storage(cells)?;
        model.spend(dimension + map.len())?;

        let mut degrees = Vec::with_capacity(dimension);
        degrees.extend((0..target.len()).map(|index| target.raw_degree(index)));
        for index in 0..source.len() {
            let degree = source.raw_degree(index).checked_add(1).ok_or_else(|| {
                Error::plain("suspension degree in mixed mapping cone is out of range")
            })?;
            degrees.push(degree);
        }

        let mut columns = Vec::with_capacity(dimension);
        for index in 0..target.len() {
            let column = target.total_column(model, index)?;
            cells += 2 * column.len();
            check_storage(cells)?;
            columns.push(column);
        }
        let mut map_entries = map.iter().map(|&(i, c)| (hom.split(i), c)).peekable();
        for index in 0..source.len() {
            let mut column = Vec::new();
            while let Some(((_, target_index), coefficient)) =
                map_entries.next_if(|&((source_index, _), _)| source_index == index)
            {
                column.push((target_index, coefficient));
            }
            let negative = source.total_column(model, index)?;
            cells += 2 * negative.len();
            check_storage(cells)?;
            column.extend(
                negative
                    .into_iter()
                    .map(|(i, c)| (target.len() + i, source.prime.negate(c))),
            );
            columns.push(column);
        }

        let complex = Complex::from_total(model, source.prime, degrees, columns)?;
        let cone = Self {
            source: hom.source,
            target: hom.target,
            target_len: target.len(),
            map,
        };
        Ok((cone, complex))
    }

    pub(crate) fn map_grade(&self, cone: Id, argument: Grade, op: Operation) -> Result<Grade> {
        let (source, target, offset) = match op {
            Operation::Include => (self.target, cone, 0),
            Operation::Lift => (self.source, cone, 1),
            Operation::Project => (cone, self.source, -1),
        };
        if argument.complex != source {
            return Err(Error::plain(
                "mixed mapping cone operation has the wrong endpoint",
            ));
        }
        Grade {
            complex: target,
            b: argument.b,
        }
        .shift(offset)
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        op: Operation,
        input: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        model.spend(input.len())?;
        let output = match op {
            Operation::Include => input.to_vec(),
            Operation::Lift => input
                .iter()
                .map(|&(i, c)| (self.target_len + i, c))
                .collect(),
            Operation::Project => {
                let start = input.partition_point(|&(i, _)| i < self.target_len);
                input[start..]
                    .iter()
                    .map(|&(i, c)| (i - self.target_len, c))
                    .collect()
            }
        };
        Ok(output)
    }

    pub(crate) fn check_shear(&self, original: &Self, model: &Model, map_id: Id) -> Result<()> {
        let hom = model.mixed_hom(map_id)?;
        if self.source != original.source
            || self.target != original.target
            || self.source != hom.source
            || self.target != hom.target
        {
            return Err(Error::plain(
                "mixed mapping cone shear has incompatible endpoints",
            ));
        }
        Ok(())
    }

    pub(crate) fn shear(
        &self,
        model: &mut Model,
        original: &Self,
        map_id: Id,
        homotopy: &[Coordinate],
        input: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let maps = model.mixed_complex(map_id);
        let hom = model.mixed_hom(map_id)?;
        let boundary = maps.apply(model, Differential::Total, homotopy)?;
        if add_scaled(model, maps.prime, &original.map, &boundary, 1)? != *self.map {
            return Err(Error::plain(
                "mixed cone maps must differ by the differential of the given homotopy",
            ));
        }
        let projected = original.apply(model, Operation::Project, input)?;
        let correction = hom.evaluate(model, homotopy, &projected)?;
        add_scaled(model, maps.prime, input, &correction, maps.prime.negate(1))
    }
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "mixed mapping cone exceeds the matrix storage limit",
        ));
    }
    Ok(())
}
