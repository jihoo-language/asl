use crate::error::{Error, Result};
use crate::model::Model;

use super::Field;

struct Pivot {
    image: Vec<u32>,
    preimage: Vec<u32>,
}

pub(super) struct Solver {
    columns: Vec<Vec<u32>>,
    pivots: Vec<Option<Pivot>>,
    pub(super) rank: usize,
}

impl Solver {
    pub(super) fn new(
        model: &mut Model,
        field: &Field,
        width: usize,
        columns: Vec<Vec<u32>>,
    ) -> Result<Self> {
        let mut pivots: Vec<Option<Pivot>> = (0..width).map(|_| None).collect();
        let mut rank = 0;
        for (i, original) in columns.iter().enumerate() {
            if original.len() != width {
                return Err(Error::plain(
                    "partition formal linear map has inconsistent column lengths",
                ));
            }
            let mut image = original.clone();
            let mut preimage = vec![0; columns.len()];
            preimage[i] = 1;
            for (j, entry) in pivots.iter().enumerate() {
                if let Some(row) = entry {
                    let scale = image[j];
                    if scale != 0 {
                        model.spend(image.len() + preimage.len())?;
                        subtract(field, &mut image, &row.image, scale);
                        subtract(field, &mut preimage, &row.preimage, scale);
                    }
                }
            }
            if let Some(pivot) = image.iter().position(|&x| x != 0) {
                let inv = field.inverse(model, image[pivot])?;
                for x in image.iter_mut().chain(&mut preimage) {
                    *x = field.mul(inv, *x);
                }
                pivots[pivot] = Some(Pivot { image, preimage });
                rank += 1;
            }
        }
        Ok(Self {
            columns,
            pivots,
            rank,
        })
    }

    pub(super) fn solve(
        &self,
        model: &mut Model,
        field: &Field,
        rhs: &[u32],
    ) -> Result<(Vec<u32>, Vec<u32>)> {
        debug_assert_eq!(rhs.len(), self.pivots.len());
        let mut residual = rhs.to_vec();
        let mut solution = vec![0; self.columns.len()];
        for (i, row) in self.pivots.iter().enumerate() {
            if let Some(row) = row {
                let scale = residual[i];
                if scale == 0 {
                    continue;
                }
                model.spend(residual.len() + solution.len())?;
                subtract(field, &mut residual, &row.image, scale);
                for (x, &v) in solution.iter_mut().zip(&row.preimage) {
                    *x = field.add(*x, field.mul(scale, v));
                }
            }
        }
        Ok((solution, residual))
    }

    pub(super) fn obstruction(
        &self,
        model: &mut Model,
        field: &Field,
        rhs: &[u32],
        residual: &[u32],
    ) -> Result<(Vec<u32>, u32)> {
        let pivot = residual.iter().position(|&x| x != 0).ok_or_else(|| {
            Error::plain("partition formal obstruction requires a nonzero residual")
        })?;
        let mut row = Vec::with_capacity(rhs.len());
        for i in 0..rhs.len() {
            let mut unit = vec![0; rhs.len()];
            unit[i] = 1;
            let (_, reduced) = self.solve(model, field, &unit)?;
            row.push(reduced[pivot]);
        }
        let dot = |rhs: &[u32]| {
            row.iter()
                .zip(rhs)
                .fold(0, |x, (&a, &b)| field.add(x, field.mul(a, b)))
        };
        let value = dot(rhs);
        if self.columns.iter().any(|c| dot(c) != 0) || value == 0 || value != residual[pivot] {
            return Err(Error::plain(
                "partition formal obstruction functional failed its certificate",
            ));
        }
        Ok((row, value))
    }
}

fn subtract(field: &Field, target: &mut [u32], source: &[u32], scale: u32) {
    for (a, &b) in target.iter_mut().zip(source) {
        *a = field.sub(*a, field.mul(b, scale));
    }
}
