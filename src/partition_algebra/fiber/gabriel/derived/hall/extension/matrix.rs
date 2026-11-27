use crate::error::{Error, Result};
use crate::model::Model;

use super::super::linear::Span;
use super::super::{Field, Hall, Shape};

#[derive(Clone)]
pub(super) struct Matrix {
    pub(super) rows: usize,
    pub(super) cols: usize,
    entries: Vec<u32>,
}

impl Matrix {
    pub(super) fn new(rows: usize, cols: usize, entries: Vec<u32>) -> Self {
        debug_assert_eq!(entries.len(), rows * cols);
        Self {
            rows,
            cols,
            entries,
        }
    }

    pub(super) fn zero(rows: usize, cols: usize) -> Self {
        Self::new(rows, cols, vec![0; rows * cols])
    }

    pub(super) fn get(&self, i: usize, j: usize) -> u32 {
        self.entries[i * self.cols + j]
    }

    fn set(&mut self, i: usize, j: usize, a: u32) {
        self.entries[i * self.cols + j] = a;
    }

    pub(super) fn identity(n: usize) -> Self {
        let mut one = Self::zero(n, n);
        for i in 0..n {
            one.set(i, i, 1);
        }
        one
    }

    pub(super) fn compose(&self, model: &mut Model, field: &Field, rhs: &Self) -> Result<Self> {
        if self.cols != rhs.rows {
            return Err(Error::plain(
                "partition extension path matrices have incompatible dimensions",
            ));
        }
        model.spend(1 + self.rows * rhs.cols * (self.cols + 1))?;
        let mut product = Self::zero(self.rows, rhs.cols);
        for i in 0..self.rows {
            for k in 0..self.cols {
                let x = self.get(i, k);
                if x == 0 {
                    continue;
                }
                for j in 0..rhs.cols {
                    let y = rhs.get(k, j);
                    product.set(i, j, field.add(product.get(i, j), field.mul(x, y)));
                }
            }
        }
        Ok(product)
    }

    pub(super) fn is_zero(&self) -> bool {
        self.entries.iter().all(|&x| x == 0)
    }

    pub(super) fn rank(&self, model: &mut Model, field: &Field) -> Result<usize> {
        let mut span = Span::new(self.cols);
        if self.cols == 0 {
            return Ok(0);
        }
        for row in self.entries.chunks(self.cols) {
            span.insert(model, field, row.to_vec())?;
        }
        Ok(span.rank)
    }

    pub(super) fn upper_block(a: &Self, b: &Self, values: &[u32]) -> Self {
        debug_assert_eq!(values.len(), a.rows * b.cols);
        let mut out = Self::zero(a.rows + b.rows, a.cols + b.cols);
        for i in 0..a.rows {
            for j in 0..a.cols {
                out.set(i, j, a.get(i, j));
            }
            for j in 0..b.cols {
                out.set(i, a.cols + j, values[i * b.cols + j]);
            }
        }
        for i in 0..b.rows {
            for j in 0..b.cols {
                out.set(i + a.rows, j + a.cols, b.get(i, j));
            }
        }
        out
    }
}

pub(super) struct Representation {
    pub(super) dims: Vec<usize>,
    pub(super) arrows: Vec<Option<Matrix>>,
}

impl Representation {
    pub(super) fn new(hall: &Hall, shape: &Shape) -> Self {
        let mut positions = vec![Vec::new(); hall.vertices];
        let mut copy = 0;
        for &(i, l, m) in &shape.0 {
            for _ in 0..m {
                for (layer, &v) in hall.path_heads[usize::from(i)][..usize::from(l)]
                    .iter()
                    .enumerate()
                {
                    positions[v].push((copy, layer));
                }
                copy += 1;
            }
        }
        let dims: Vec<_> = positions.iter().map(Vec::len).collect();
        let arrows = hall
            .successor
            .iter()
            .enumerate()
            .map(|(i, &next)| {
                let j = next?;
                let mut matrix = Matrix::zero(dims[j], dims[i]);
                for (column, &(copy, layer)) in positions[i].iter().enumerate() {
                    if let Some(row) = positions[j]
                        .iter()
                        .position(|&(c, k)| c == copy && k == layer + 1)
                    {
                        matrix.set(row, column, 1);
                    }
                }
                Some(matrix)
            })
            .collect();
        Self { dims, arrows }
    }

    pub(super) fn path(
        &self,
        model: &mut Model,
        hall: &Hall,
        head: usize,
        length: usize,
    ) -> Result<Option<(usize, Matrix)>> {
        let mut at = head;
        let mut path = Matrix::identity(self.dims[head]);
        for _ in 0..length {
            let Some(next) = hall.successor[at] else {
                return Ok(None);
            };
            path = self.arrows[at]
                .as_ref()
                .unwrap()
                .compose(model, &hall.field, &path)?;
            at = next;
        }
        Ok(Some((at, path)))
    }
}
