use std::cmp::Reverse;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::mixed::{Complex, Differential};
use crate::model::Model;
use crate::prime_linear::{add_scaled, invert, Basis, Coordinate};

pub(crate) struct Torsion {
    pub(crate) length: u64,
    pub(crate) null_degree: i64,
    pub(crate) witness: Box<[Coordinate]>,
}

pub(crate) struct Bar {
    pub(crate) birth: i64,
    pub(crate) cycle: Box<[Coordinate]>,
    pub(crate) torsion: Option<Torsion>,
}

pub(crate) struct Bars {
    bars: Vec<Bar>,
}

impl Bars {
    pub(crate) fn build(model: &mut Model, complex: &Complex) -> Result<Self> {
        let size = complex.len();
        let mut order: Vec<usize> = (0..size).collect();
        model.spend(size * (1 + size.max(1).ilog2() as usize))?;
        order.sort_unstable_by_key(|&i| (Reverse(complex.raw_degree(i)), i));
        let mut positions = vec![0; size];
        for (i, &generator) in order.iter().enumerate() {
            positions[generator] = i;
        }

        let prime = complex.prime;
        let mut rows = Basis::default();
        let mut hits = vec![false; size];
        let mut candidates = Vec::new();
        let mut bars = Vec::new();
        let mut retained = 4 * size;
        check_storage(retained)?;
        for (column, &generator) in order.iter().enumerate() {
            let degree = complex.raw_degree(generator);
            let mut image = complex.apply(model, Differential::Total, &[(generator, 1)])?;
            reindex(model, &mut image, &positions)?;
            let (mut image, used) = rows.reduce_with_witness(model, prime, image)?;
            let mut witness = add_scaled(model, prime, &[(column, 1)], &used, prime.negate(1))?;
            let Some(&(pivot, coefficient)) = image.last() else {
                retained += 1 + 2 * witness.len();
                check_storage(retained + rows.cells())?;
                candidates.push((column, witness));
                continue;
            };

            let birth = complex.raw_degree(order[pivot]);
            let delay = i128::from(birth) - i128::from(degree) + 1;
            if delay < 0 || delay % 2 != 0 {
                return Err(Error::plain("mixed bar pivot violates homogeneity"));
            }
            let length = u64::try_from(delay / 2)
                .map_err(|_| Error::plain("mixed bar length is out of range"))?;
            hits[pivot] = true;
            if coefficient != 1 {
                let factor = invert(model, prime, coefficient)?;
                model.spend(image.len() + witness.len())?;
                for (_, entry) in image.iter_mut().chain(&mut witness) {
                    *entry = prime.mul(*entry, factor);
                }
            }
            if length == 0 {
                rows.insert_with_witness(model, prime, image, witness)?;
            } else {
                retained += 1 + 2 * (image.len() + witness.len());
                check_storage(retained + rows.cells())?;
                rows.insert_with_witness(model, prime, image.clone(), witness.clone())?;
                reindex(model, &mut image, &order)?;
                reindex(model, &mut witness, &order)?;
                bars.push(Bar {
                    birth,
                    cycle: image.into(),
                    torsion: Some(Torsion {
                        length,
                        null_degree: degree,
                        witness: witness.into(),
                    }),
                });
            }
            check_storage(retained + rows.cells())?;
        }

        // Incoming pivots replace the corresponding zero-column cycles
        // in the triangular kernel basis. The other cycles are free.
        for (position, mut cycle) in candidates {
            if !hits[position] {
                let birth = complex.raw_degree(order[position]);
                reindex(model, &mut cycle, &order)?;
                bars.push(Bar {
                    birth,
                    cycle: cycle.into(),
                    torsion: None,
                });
            }
        }

        model.spend(bars.len() * (1 + bars.len().max(1).ilog2() as usize))?;
        bars.sort_by_key(|bar| {
            (
                bar.torsion.is_some(),
                bar.birth,
                bar.torsion.as_ref().map(|t| t.length),
            )
        });
        model.allocate(
            bars.iter()
                .map(|bar| {
                    4 + 2 * bar.cycle.len()
                        + bar.torsion.as_ref().map_or(0, |t| 3 + 2 * t.witness.len())
                })
                .sum(),
        )?;
        Ok(Self { bars })
    }

    pub(crate) fn bars(&self) -> &[Bar] {
        &self.bars
    }

    pub(crate) fn bar(&self, i: usize) -> Result<&Bar> {
        self.bars
            .get(i)
            .ok_or_else(|| Error::plain("mixed bar index is out of range"))
    }

    pub(crate) fn rank(&self, model: &mut Model, source: i64, target: i64) -> Result<usize> {
        let source = i128::from(source);
        let target = i128::from(target);
        if target > source || (source - target) % 2 != 0 {
            return Err(Error::plain(
                "mixed bar rank requires a nonnegative power of u",
            ));
        }
        model.spend(self.bars.len())?;
        Ok(self
            .bars
            .iter()
            .filter(|bar| {
                let birth = i128::from(bar.birth);
                birth >= source
                    && (birth - source) % 2 == 0
                    && bar
                        .torsion
                        .as_ref()
                        .is_none_or(|t| (birth - target) / 2 < i128::from(t.length))
            })
            .count())
    }

    pub(crate) fn periodic_rank(&self, model: &mut Model, b: i64) -> Result<usize> {
        model.spend(self.bars.len())?;
        Ok(self
            .bars
            .iter()
            .filter(|bar| bar.torsion.is_none() && bar.birth.rem_euclid(2) == b.rem_euclid(2))
            .count())
    }
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "mixed barcode reduction exceeds the storage limit",
        ));
    }
    Ok(())
}

fn reindex(model: &mut Model, entries: &mut [Coordinate], positions: &[usize]) -> Result<()> {
    model.spend(entries.len() * (1 + entries.len().max(1).ilog2() as usize))?;
    for (index, _) in entries.iter_mut() {
        *index = positions[*index];
    }
    entries.sort_unstable();
    Ok(())
}
