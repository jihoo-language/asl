use std::cmp::Ordering;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::super::Matrix;
use super::{Field, Spec};

mod distribution;
mod fourier;
mod primary;
use distribution::Distribution;

const TOTAL_DIM_LIMIT: usize = 64;

pub(crate) struct Request {
    pub(crate) frame: Label,
    pub(crate) method: Label,
    pub(crate) dimensions: Vec<u128>,
}

pub(crate) enum Expression {
    Literal(Label, Vec<u128>, Vec<Vec<u128>>),
    Triangular(Box<[ParsedExpr; 2]>, Vec<Vec<u128>>),
    Character(u128, Box<ParsedExpr>),
    Fourier(Request, u128),
}

pub(crate) enum Report {
    Trace(Box<ParsedExpr>),
    Critical(Box<ParsedExpr>),
    Distribution(Request),
    Inversion(Request),
}

pub(crate) struct Representation {
    pub(crate) frame: Rc<str>,
    dimensions: Rc<[usize]>,
    arrows: Vec<Matrix>,
}

impl Representation {
    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(
            4 + self.dimensions.len()
                + self
                    .arrows
                    .iter()
                    .map(|a| 3 + a.entries().len())
                    .sum::<usize>(),
        )?;
        Ok(Value::RVPartitionQuiver(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        model.spend(
            1 + self.dimensions.len()
                + self.arrows.iter().map(|a| a.entries().len()).sum::<usize>(),
        )?;
        Ok(self
            .frame
            .cmp(&other.frame)
            .then_with(|| self.dimensions.cmp(&other.dimensions))
            .then_with(|| {
                self.arrows
                    .iter()
                    .map(Matrix::entries)
                    .cmp(other.arrows.iter().map(Matrix::entries))
            }))
    }

    pub(crate) fn render(&self) -> String {
        let arrows: Vec<_> = self.arrows.iter().map(Matrix::entries).collect();
        format!(
            "RVPartitionQuiver({},{:?}){arrows:?}",
            self.frame, self.dimensions
        )
    }
}

impl Spec {
    pub(crate) fn ambient_dimensions(&self, dims: Vec<u128>) -> Result<Rc<[usize]>> {
        if dims.len() != self.hall.vertices {
            return Err(Error::plain(
                "partition quiver dimension vector requires one dimension per vertex",
            ));
        }
        let dims = dims
            .into_iter()
            .map(|d| {
                usize::try_from(d)
                    .ok()
                    .filter(|&n| n <= TOTAL_DIM_LIMIT)
                    .ok_or_else(dim_limit)
            })
            .collect::<Result<Vec<_>>>()?;
        if dims.iter().sum::<usize>() > TOTAL_DIM_LIMIT {
            return Err(dim_limit());
        }
        Ok(dims.into())
    }

    pub(crate) fn ambient_characteristic(&self) -> usize {
        self.hall.field.characteristic() as usize
    }

    fn potential_terms(&self) -> Result<&super::Terms> {
        self.potential.as_ref().ok_or_else(super::no_potential)
    }

    fn matrix_sizes(&self, dims: &[usize]) -> Vec<(usize, usize)> {
        self.symbols[..self.arrows]
            .iter()
            .map(|arrow| (dims[arrow.target], dims[arrow.source]))
            .collect()
    }

    pub(crate) fn ambient_literal(
        &self,
        model: &mut Model,
        dimensions: Rc<[usize]>,
        entries: Vec<Vec<u128>>,
    ) -> Result<Representation> {
        let sizes = self.matrix_sizes(&dimensions);
        let arrows = matrix_codes(model, &self.hall.field, &sizes, entries)?;
        Ok(Representation {
            frame: self.name.clone(),
            dimensions,
            arrows,
        })
    }

    pub(crate) fn ambient_triangular(
        &self,
        model: &mut Model,
        x: &Representation,
        y: &Representation,
        mixed: Vec<Vec<u128>>,
    ) -> Result<Representation> {
        let dims = x
            .dimensions
            .iter()
            .zip(y.dimensions.iter())
            .map(|(a, b)| (a + b) as u128)
            .collect();
        let dimensions = self.ambient_dimensions(dims)?;
        let sizes: Vec<_> = x
            .arrows
            .iter()
            .zip(&y.arrows)
            .map(|(x, y)| (x.rows, y.cols))
            .collect();
        let mixed = matrix_codes(model, &self.hall.field, &sizes, mixed)?;
        let arrows = x
            .arrows
            .iter()
            .zip(&y.arrows)
            .zip(&mixed)
            .map(|((a, b), block)| Matrix::upper_block(a, b, block.entries()))
            .collect();
        let result = Representation {
            frame: self.name.clone(),
            dimensions,
            arrows,
        };
        if self.potential.is_some() {
            let field = &self.hall.field;
            let sum = field.add(self.ambient_trace(model, x)?, self.ambient_trace(model, y)?);
            if self.ambient_trace(model, &result)? != sum {
                return Err(Error::plain(
                    "partition quiver potential trace failed block triangular additivity",
                ));
            }
        }
        Ok(result)
    }

    fn ambient_trace(&self, model: &mut Model, rep: &Representation) -> Result<u32> {
        self.representation_trace(model, &rep.dimensions, |letter| &rep.arrows[letter])
    }

    pub(crate) fn ambient_report(
        &self,
        model: &mut Model,
        rep: &Representation,
        critical: bool,
    ) -> Result<String> {
        if critical {
            for relation in &self.gradients {
                if !self
                    .path_matrix(model, &rep.dimensions, relation, |letter| {
                        &rep.arrows[letter]
                    })?
                    .is_zero()
                {
                    return Ok(format!("rvPartitionQuiverCritical[{}] = false", self.name));
                }
            }
            return Ok(format!("rvPartitionQuiverCritical[{}] = true", self.name));
        }
        Ok(format!(
            "rvPartitionQuiverTrace[{}] = {}",
            self.name,
            self.ambient_trace(model, rep)?
        ))
    }

    pub(crate) fn ambient_distribution_report(
        &self,
        model: &mut Model,
        dims: Rc<[usize]>,
        method: &Label,
        inversion: bool,
    ) -> Result<String> {
        let counts = self.ambient_distribution(model, dims, method)?;
        if inversion {
            return counts.inversion(model, &self.hall.field, &self.name);
        }
        Ok(counts.render(&self.name))
    }

    pub(crate) fn ambient_fourier(
        &self,
        model: &mut Model,
        dims: Rc<[usize]>,
        method: &Label,
        lambda: u128,
    ) -> Result<Value> {
        let lambda = field_code(&self.hall.field, lambda)?;
        let counts = self.ambient_distribution(model, dims, method)?;
        counts.fourier(model, &self.hall.field, lambda)
    }

    pub(crate) fn ambient_character(
        &self,
        model: &mut Model,
        rep: &Representation,
        lambda: u128,
    ) -> Result<Value> {
        let lambda = field_code(&self.hall.field, lambda)?;
        let field = &self.hall.field;
        let value = field.mul(lambda, self.ambient_trace(model, rep)?);
        let exponent = fourier::absolute_trace(model, field, value)?;
        crate::rv_cyclotomic::root(model, self.ambient_characteristic(), i128::from(exponent))
    }

    fn ambient_distribution(
        &self,
        model: &mut Model,
        dims: Rc<[usize]>,
        method: &Label,
    ) -> Result<Distribution> {
        let weights = self.potential_terms()?;
        let field = &self.hall.field;
        let shape = self.matrix_sizes(&dims);
        let entries = shape.iter().map(|(r, s)| r * s).sum();
        let (total, group) = distribution::mass(field, &dims, entries)?;
        let structured = self.components.iter().all(|(vertices, length, _)| {
            matches!((vertices.as_slice(), *length), ([_], 2) | ([_, _], 3))
        });
        let framed = match method.text.as_str() {
            "auto" if !structured => self.direct_counts(model, dims, &shape, total)?,
            "direct" => self.direct_counts(model, dims, &shape, total)?,
            "auto" | "structure" => {
                let mut counts = distribution::unit(field);
                for (vertices, length, _) in &self.components {
                    let start = vertices[0];
                    let coefficient = weights
                        .iter()
                        .find(|(w, _)| w.start == start)
                        .map(|(_, &c)| c)
                        .ok_or_else(|| {
                            Error::plain("partition quiver potential is missing a cyclic component")
                        })?;
                    let one = match (vertices.as_slice(), *length) {
                        ([v], 2) => primary::counts(model, field, dims[*v], coefficient)?,
                        ([u,v], 3) => distribution::two_cycle(model, field, [dims[*u],dims[*v]], coefficient)?,
                        _ => return Err(Error::plain("structured partition quiver counting requires a cubic loop or a quartic two-cycle")),
                    };
                    counts = distribution::convolve(model, field, &counts, &one, 1)?;
                }
                counts
            }
            _ => {
                return Err(Error::new(
                    method.span,
                    "partition quiver count method must be auto, structure, or direct",
                ))
            }
        };
        Distribution::new(framed, total, group)
    }

    fn direct_counts(
        &self,
        model: &mut Model,
        dims: Rc<[usize]>,
        sizes: &[(usize, usize)],
        total: u128,
    ) -> Result<Vec<u128>> {
        const DIRECT_LIMIT: u128 = 500_000;
        if total > DIRECT_LIMIT {
            return Err(Error::plain(
                "direct partition quiver count exceeds 500,000 matrix assignments",
            ));
        }
        let q = self.hall.field.order();
        let mut codes = vec![0; sizes.iter().map(|(r, s)| r * s).sum()];
        model.spend(total as usize * (codes.len() + 1))?;
        let mut counts = vec![0; q as usize];
        for _ in 0..total {
            let mut offset = 0;
            let arrows = sizes
                .iter()
                .map(|&(r, s)| {
                    let end = offset + r * s;
                    let matrix = Matrix::new(r, s, codes[offset..end].to_vec());
                    offset = end;
                    matrix
                })
                .collect();
            let rep = Representation {
                frame: self.name.clone(),
                dimensions: dims.clone(),
                arrows,
            };
            let value = self.ambient_trace(model, &rep)?;
            counts[value as usize] += 1;
            for digit in &mut codes {
                *digit += 1;
                if *digit != q {
                    break;
                }
                *digit = 0;
            }
        }
        Ok(counts)
    }
}

fn matrix_codes(
    model: &mut Model,
    field: &Field,
    sizes: &[(usize, usize)],
    entries: Vec<Vec<u128>>,
) -> Result<Vec<Matrix>> {
    if sizes.len() != entries.len() {
        return Err(Error::plain(
            "partition quiver requires one matrix entry list per original arrow",
        ));
    }
    sizes
        .iter()
        .zip(entries)
        .map(|(&(r, s), codes)| {
            if r * s != codes.len() {
                return Err(Error::plain(
                    "partition quiver arrow entry count does not match its matrix dimensions",
                ));
            }
            model.spend(1 + codes.len())?;
            let entries = codes
                .into_iter()
                .map(|c| field_code(field, c))
                .collect::<Result<Vec<_>>>()?;
            Ok(Matrix::new(r, s, entries))
        })
        .collect()
}

fn dim_limit() -> Error {
    Error::plain("partition quiver total dimension exceeds 64")
}

fn field_code(field: &Field, code: u128) -> Result<u32> {
    u32::try_from(code)
        .ok()
        .filter(|&x| x < field.order())
        .ok_or_else(|| Error::plain("partition quiver entry is outside the finite field codes"))
}
