use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::prime_linear::{Basis, Coordinate, QuotientBasis};
use crate::sheaf::{self, Matrix};

use super::{linear, Disk};

pub(crate) struct Pair {
    pub(crate) a: Rc<Matrix>,
    pub(crate) b: Rc<Matrix>,
}

impl Pair {
    pub(crate) fn dual(&self, model: &mut Model) -> Result<Self> {
        Ok(Self {
            a: Rc::new(linear::transpose(model, &self.a)?),
            b: Rc::new(linear::transpose(model, &self.b)?),
        })
    }

    pub(crate) fn compose(
        &self,
        model: &mut Model,
        prime: crate::odd_cobar::Prime,
        before: &Self,
    ) -> Result<Self> {
        Ok(Self {
            a: Rc::new(Matrix::compose(model, prime, &self.a, &before.a)?),
            b: Rc::new(Matrix::compose(model, prime, &self.b, &before.b)?),
        })
    }
}

pub(crate) struct Hom {
    pub(crate) source: Rc<Disk>,
    pub(crate) target: Rc<Disk>,
    kernel: QuotientBasis,
    nearby_cells: usize,
}

impl Hom {
    pub(crate) fn build(model: &mut Model, source: Rc<Disk>, target: Rc<Disk>) -> Result<Self> {
        let prime = source.prime;
        if target.prime != prime {
            return Err(Error::plain(
                "disk Hom requires the same prime field on both sides",
            ));
        }
        let (n, m) = (source.nearby, source.vanishing);
        let (nn, mm) = (target.nearby, target.vanishing);
        let nearby_cells = n.saturating_mul(nn);
        let size = nearby_cells.saturating_add(m.saturating_mul(mm));
        let second = n.saturating_mul(mm);
        sheaf::check_cells(
            size.saturating_add(second)
                .saturating_add(m.saturating_mul(nn)),
        )?;
        let u_transpose = linear::transpose(model, &source.u)?;
        let v_transpose = linear::transpose(model, &source.v)?;
        let kernel = QuotientBasis::build(model, prime, Basis::default(), size, |model, index| {
            // The two residuals are u_target A - B u_source and
            // v_target B - A v_source. Both are stored column by column.
            let mut output = Vec::new();
            if index < nearby_cells {
                let (col, row) = (index / nn, index % nn);
                let forward = target.u.column(row);
                let backward = v_transpose.column(col);
                model.spend(forward.len() + backward.len())?;
                output.extend(forward.iter().map(|&(i, x)| (col * mm + i, x)));
                output.extend(
                    backward
                        .iter()
                        .map(|&(j, x)| (second + j * nn + row, prime.negate(x))),
                );
            } else {
                let offset = index - nearby_cells;
                let (col, row) = (offset / mm, offset % mm);
                let backward = u_transpose.column(col);
                let forward = target.v.column(row);
                model.spend(backward.len() + forward.len())?;
                output.extend(
                    backward
                        .iter()
                        .map(|&(j, x)| (j * mm + row, prime.negate(x))),
                );
                output.extend(forward.iter().map(|&(i, x)| (second + col * nn + i, x)));
            }
            Ok(output)
        })?;
        Ok(Self {
            source,
            target,
            kernel,
            nearby_cells,
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.kernel.rank()
    }

    pub(crate) fn literal(&self, model: &mut Model, a: Vec<i64>, b: Vec<i64>) -> Result<Value> {
        let p = self.source.prime;
        let a = Rc::new(linear::literal(
            model,
            p,
            self.target.nearby,
            self.source.nearby,
            a,
        )?);
        let b = Rc::new(linear::literal(
            model,
            p,
            self.target.vanishing,
            self.source.vanishing,
            b,
        )?);
        self.encode(model, &Pair { a, b })
    }

    pub(crate) fn identity(&self, model: &mut Model) -> Result<Value> {
        let a = Rc::new(Matrix::identity(model, self.source.nearby)?);
        let b = Rc::new(Matrix::identity(model, self.source.vanishing)?);
        self.encode(model, &Pair { a, b })
    }

    pub(crate) fn encode(&self, model: &mut Model, pair: &Pair) -> Result<Value> {
        let mut coordinates = Vec::new();
        coordinates.extend(pair.a.columns().flat_map(|(col, entries)| {
            entries
                .iter()
                .map(move |&(i, x)| (col * self.target.nearby + i, x))
        }));
        coordinates.extend(pair.b.columns().flat_map(|(col, entries)| {
            entries
                .iter()
                .map(move |&(i, x)| (self.nearby_cells + col * self.target.vanishing + i, x))
        }));
        let value = self
            .kernel
            .project(model, self.source.prime, coordinates)?
            .ok_or_else(|| Error::plain("disk morphism must intertwine both u and v"))?;
        crate::odd_ext::persist(model, value)
    }

    pub(crate) fn decode(&self, model: &mut Model, value: &[Coordinate]) -> Result<Pair> {
        let coordinates = self
            .kernel
            .representative(model, self.source.prime, value)?;
        let mid = coordinates.partition_point(|&(i, _)| i < self.nearby_cells);
        Ok(Pair {
            a: Rc::new(reshape(&coordinates[..mid], self.target.nearby, 0)?),
            b: Rc::new(reshape(
                &coordinates[mid..],
                self.target.vanishing,
                self.nearby_cells,
            )?),
        })
    }

    pub(crate) fn report(&self, model: &mut Model, value: &[Coordinate]) -> Result<String> {
        let pair = self.decode(model, value)?;
        let a = linear::report(model, &pair.a, self.target.nearby, self.source.nearby)?;
        let b = linear::report(model, &pair.b, self.target.vanishing, self.source.vanishing)?;
        Ok(format!("{a} / {b}"))
    }
}

fn reshape(support: &[Coordinate], height: usize, offset: usize) -> Result<Matrix> {
    let mut columns = Vec::new();
    for entries in support.chunk_by(|x, y| (x.0 - offset) / height == (y.0 - offset) / height) {
        let j = (entries[0].0 - offset) / height;
        let col: Box<[_]> = entries
            .iter()
            .map(|&(i, x)| ((i - offset) % height, x))
            .collect();
        columns.push((j, col));
    }
    Matrix::from_columns(columns)
}
