use std::collections::btree_map::Entry;
use std::collections::BTreeMap;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::mixed::{Complex, Differential, Id};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::Coordinate;

pub(crate) struct Hom {
    pub(crate) source: Id,
    pub(crate) target: Id,
    rows: usize,
    columns: usize,
    prime: Prime,
}

impl Hom {
    pub(crate) fn build(model: &mut Model, source: Id, target: Id) -> Result<(Self, Complex)> {
        let domain = model.mixed_complex(source);
        let range = model.mixed_complex(target);
        if domain.prime != range.prime {
            return Err(Error::plain(
                "mixhom endpoints must have the same prime field",
            ));
        }
        let hom = Self {
            source,
            target,
            rows: range.len(),
            columns: domain.len(),
            prime: domain.prime,
        };
        let dimension = hom.columns.saturating_mul(hom.rows);
        check_storage(6usize.saturating_mul(dimension))?;
        model.spend(dimension)?;
        let mut degrees = Vec::with_capacity(dimension);
        for x in 0..hom.columns {
            for y in 0..hom.rows {
                let degree = i128::from(range.raw_degree(y)) - i128::from(domain.raw_degree(x));
                degrees.push(
                    i64::try_from(degree)
                        .map_err(|_| Error::plain("mixed Hom basis degree is out of range"))?,
                );
            }
        }
        let mut cells = 6 * dimension;
        let down = hom.matrices(
            model,
            &domain,
            &range,
            Differential::Down,
            &degrees,
            &mut cells,
        )?;
        let up = hom.matrices(
            model,
            &domain,
            &range,
            Differential::Up,
            &degrees,
            &mut cells,
        )?;
        let complex = Complex::from_columns(model, hom.prime, degrees, down, up)?;
        Ok((hom, complex))
    }

    fn matrices(
        &self,
        model: &mut Model,
        source: &Complex,
        target: &Complex,
        operation: Differential,
        degrees: &[i64],
        cells: &mut usize,
    ) -> Result<Vec<Vec<Coordinate>>> {
        if degrees.is_empty() {
            return Ok(Vec::new());
        }
        let mut source_rows = vec![Vec::new(); self.columns];
        let mut row_cells = self.columns;
        for x in 0..self.columns {
            let column = source.column(operation, x);
            model.spend(1 + column.len())?;
            row_cells += 2 * column.len();
            check_storage(*cells + row_cells)?;
            for &(z, coefficient) in column {
                source_rows[z].push((x, coefficient));
            }
        }
        let mut columns = Vec::with_capacity(degrees.len());
        for (x, row) in source_rows.iter().enumerate() {
            for y in 0..self.rows {
                let target_column = target.column(operation, y);
                let count = row.len() + target_column.len();
                *cells += 2 * count;
                check_storage(*cells + row_cells)?;
                model.spend(1 + count * (1 + count.max(1).ilog2() as usize))?;
                let mut column = Vec::with_capacity(count);
                column.extend(
                    target_column
                        .iter()
                        .map(|&(z, coefficient)| (self.index(x, z), coefficient)),
                );
                let even = degrees[self.index(x, y)] % 2 == 0;
                column.extend(row.iter().map(|&(z, coefficient)| {
                    (
                        self.index(z, y),
                        if even {
                            self.prime.negate(coefficient)
                        } else {
                            coefficient
                        },
                    )
                }));
                column.sort_unstable();
                columns.push(column);
            }
        }
        Ok(columns)
    }

    pub(crate) fn index(&self, source: usize, target: usize) -> usize {
        source * self.rows + target
    }

    fn split(&self, index: usize) -> (usize, usize) {
        (index / self.rows, index % self.rows)
    }

    pub(crate) fn identity(&self, model: &mut Model) -> Result<Vec<Coordinate>> {
        if self.source != self.target {
            return Err(Error::plain("mixId requires a mixed endomorphism complex"));
        }
        model.spend(self.columns)?;
        Ok((0..self.columns).map(|i| (self.index(i, i), 1)).collect())
    }

    pub(crate) fn evaluate(
        &self,
        model: &mut Model,
        map: &[Coordinate],
        argument: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        if argument.is_empty() {
            return Ok(Vec::new());
        }
        let comparisons = 1 + argument.len().ilog2() as usize;
        model.spend(map.len() * comparisons)?;
        let mut result = Sum::new(self.prime);
        for &(index, coefficient) in map {
            let (source, target) = self.split(index);
            if let Ok(i) = argument.binary_search_by_key(&source, |&(x, _)| x) {
                result.add(model, target, coefficient, argument[i].1)?;
            }
        }
        Ok(result.finish())
    }

    pub(crate) fn compose(
        &self,
        model: &mut Model,
        outer_hom: &Self,
        inner_hom: &Self,
        outer: &[Coordinate],
        inner: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        if outer.is_empty() || inner.is_empty() {
            return Ok(Vec::new());
        }
        let mut result = Sum::new(self.prime);
        let comparisons = 1 + outer.len().ilog2() as usize;
        model.spend(inner.len() * comparisons * 2)?;
        for &(index, coefficient) in inner {
            let (source, intermediate) = inner_hom.split(index);
            let first = outer.partition_point(|&(i, _)| outer_hom.split(i).0 < intermediate);
            let last = outer.partition_point(|&(i, _)| outer_hom.split(i).0 <= intermediate);
            for &(i, outer_coefficient) in &outer[first..last] {
                let (_, target) = outer_hom.split(i);
                result.add(
                    model,
                    self.index(source, target),
                    outer_coefficient,
                    coefficient,
                )?;
            }
        }
        Ok(result.finish())
    }
}

struct Sum {
    prime: Prime,
    entries: BTreeMap<usize, u32>,
}

impl Sum {
    fn new(prime: Prime) -> Self {
        Self {
            prime,
            entries: BTreeMap::new(),
        }
    }

    fn add(&mut self, model: &mut Model, index: usize, a: u32, b: u32) -> Result<()> {
        model.spend(1 + self.entries.len().max(1).ilog2() as usize)?;
        let scalar = self.prime.mul(a, b);
        let cells = 3 * (self.entries.len() + 1);
        match self.entries.entry(index) {
            Entry::Occupied(mut entry) => {
                let coefficient = self.prime.add(*entry.get(), scalar);
                if coefficient == 0 {
                    entry.remove();
                } else {
                    *entry.get_mut() = coefficient;
                }
            }
            Entry::Vacant(entry) => {
                check_storage(cells)?;
                entry.insert(scalar);
            }
        }
        Ok(())
    }

    fn finish(self) -> Vec<Coordinate> {
        self.entries.into_iter().collect()
    }
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain("mixed Hom complex exceeds the storage limit"));
    }
    Ok(())
}
