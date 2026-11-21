use std::collections::btree_map::Entry;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::binary_field::BinaryField;
use crate::comodule::Generator;
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::rank_jets::{Matrix, Path};

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct PathId(pub(crate) usize);

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct Diagram {
    pub(crate) source: TypeId,
    pub(crate) target: TypeId,
    pub(crate) path: PathId,
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct Grade {
    pub(crate) diagram: Diagram,
    pub(crate) parity: usize,
}

impl Grade {
    pub(crate) fn differential(self) -> Self {
        Self {
            parity: self.parity ^ 1,
            ..self
        }
    }

    pub(crate) fn compose(outer: Self, inner: Self) -> Result<Self> {
        if outer.diagram.path != inner.diagram.path {
            return Err(Error::plain("ecomp requires the same Margolis arc"));
        }
        if outer.diagram.source != inner.diagram.target {
            return Err(Error::plain("ecomp requires composable JetHom endpoints"));
        }
        Ok(Self {
            diagram: Diagram {
                source: inner.diagram.source,
                target: outer.diagram.target,
                path: outer.diagram.path,
            },
            parity: outer.parity ^ inner.parity,
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct Term {
    pub(crate) power: usize,
    pub(crate) source: Generator,
    pub(crate) target: Generator,
}

impl Term {
    pub(crate) fn parity(self) -> usize {
        (self.source.degree ^ self.target.degree) & 1
    }
}

pub(crate) type WeightedTerm = (Term, usize);
pub(crate) type Support = Rc<[WeightedTerm]>;

const TERM_CELLS: usize = 6;

pub(crate) struct Complex {
    matrix: Matrix,
    sources: Vec<Generator>,
    targets: Vec<Generator>,
}

impl Complex {
    pub(crate) fn build(model: &mut Model, diagram: Diagram) -> Result<Self> {
        let path = model.jet_path(diagram.path);
        let matrix = Matrix::build(model, diagram.source, diagram.target, &path)?;
        let cells = 2 * (matrix.source_size() + matrix.target_size());
        check_storage(matrix.cells() + cells)?;
        model.spend(cells)?;
        model.allocate(matrix.cells() + cells)?;
        let sources = matrix.source().generators().collect();
        let targets = matrix.target().generators().collect();
        Ok(Self {
            matrix,
            sources,
            targets,
        })
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.matrix.is_zero()
    }

    pub(crate) fn differential(
        &self,
        model: &mut Model,
        path: &Path,
        terms: &[WeightedTerm],
    ) -> Result<Vec<WeightedTerm>> {
        let mut sum = Sum::default();
        for &(term, factor) in terms {
            model.spend(path.field.degree())?;
            for (image, coefficient) in self.column(model, path, term)? {
                sum.add(model, image, path.field.multiply(factor, coefficient))?;
            }
        }
        Ok(sum.finish())
    }

    pub(crate) fn column(
        &self,
        model: &mut Model,
        path: &Path,
        term: Term,
    ) -> Result<Vec<WeightedTerm>> {
        let source = self
            .matrix
            .source()
            .coordinate(term.source)
            .ok_or_else(outside_basis)?;
        let target = self
            .matrix
            .target()
            .coordinate(term.target)
            .ok_or_else(outside_basis)?;
        self.matrix
            .column(model, path, term.power, source, target)
            .map(|image| {
                image
                    .into_iter()
                    .map(|entry| {
                        let term = Term {
                            power: entry.power,
                            source: self.sources[entry.source],
                            target: self.targets[entry.target],
                        };
                        (term, entry.coefficient)
                    })
                    .collect()
            })
    }

    pub(crate) fn basis(
        &self,
        model: &mut Model,
        parity: usize,
        powers: std::ops::Range<usize>,
        limit: usize,
    ) -> Result<Vec<Term>> {
        if powers.is_empty() {
            return Ok(Vec::new());
        }
        let [x_even, x_odd] = parity_ranks(model, &self.sources)?;
        let [y_even, y_odd] = parity_ranks(model, &self.targets)?;
        let paired = if parity == 0 {
            x_even
                .saturating_mul(y_even)
                .saturating_add(x_odd.saturating_mul(y_odd))
        } else {
            x_even
                .saturating_mul(y_odd)
                .saturating_add(x_odd.saturating_mul(y_even))
        };
        let size = paired
            .checked_mul(powers.len())
            .filter(|&size| size <= limit)
            .ok_or_else(|| Error::plain("JetHom basis exceeds the finite basis limit"))?;
        if size == 0 {
            return Ok(Vec::new());
        }
        model.spend(size + self.sources.len() + self.targets.len())?;
        let mut targets = [Vec::new(), Vec::new()];
        for &target in &self.targets {
            targets[target.degree & 1].push(target);
        }
        let mut basis = Vec::with_capacity(size);
        for power in powers {
            model.spend(self.sources.len())?;
            for &source in &self.sources {
                basis.extend(
                    targets[(source.degree & 1) ^ parity]
                        .iter()
                        .map(|&target| Term {
                            power,
                            source,
                            target,
                        }),
                );
            }
        }
        Ok(basis)
    }
}

fn parity_ranks(model: &mut Model, generators: &[Generator]) -> Result<[usize; 2]> {
    model.spend(generators.len())?;
    let even = generators
        .iter()
        .filter(|generator| generator.degree & 1 == 0)
        .count();
    Ok([even, generators.len() - even])
}

fn outside_basis() -> Error {
    Error::plain("term is outside its JetHom basis")
}

#[derive(Default)]
struct Sum(BTreeMap<Term, usize>);

impl Sum {
    fn add(&mut self, model: &mut Model, term: Term, coefficient: usize) -> Result<()> {
        if coefficient == 0 {
            return Ok(());
        }
        model.spend(1 + self.0.len().max(1).ilog2() as usize)?;
        let size = self.0.len();
        match self.0.entry(term) {
            Entry::Occupied(mut entry) => {
                *entry.get_mut() ^= coefficient;
                if *entry.get() == 0 {
                    entry.remove();
                }
            }
            Entry::Vacant(entry) => {
                check_storage(TERM_CELLS * (size + 1))?;
                entry.insert(coefficient);
            }
        }
        Ok(())
    }

    fn finish(self) -> Vec<WeightedTerm> {
        self.0.into_iter().collect()
    }
}

pub(crate) fn zero() -> Value {
    Value::JetHom([].into())
}

pub(crate) fn persist(model: &mut Model, terms: Vec<WeightedTerm>) -> Result<Value> {
    model.allocate(TERM_CELLS * terms.len())?;
    Ok(Value::JetHom(terms.into()))
}

pub(crate) fn add(model: &mut Model, a: &[WeightedTerm], b: &[WeightedTerm]) -> Result<Value> {
    let mut sum = Sum::default();
    for &(term, coefficient) in a.iter().chain(b) {
        sum.add(model, term, coefficient)?;
    }
    persist(model, sum.finish())
}

pub(crate) fn scale(
    model: &mut Model,
    field: BinaryField,
    terms: &[WeightedTerm],
    scalar: usize,
) -> Result<Value> {
    if scalar == 0 {
        return Ok(zero());
    }
    model.spend(field.degree() * terms.len())?;
    persist(
        model,
        terms
            .iter()
            .map(|&(term, coefficient)| (term, field.multiply(coefficient, scalar)))
            .collect(),
    )
}

pub(crate) fn compose(
    model: &mut Model,
    path: &Path,
    outer: &[WeightedTerm],
    inner: &[WeightedTerm],
) -> Result<Vec<WeightedTerm>> {
    model.spend(outer.len().saturating_mul(inner.len()))?;
    let mut sum = Sum::default();
    for &(g, b) in outer {
        for &(f, a) in inner {
            if f.target == g.source && g.power < path.order - f.power {
                model.spend(path.field.degree())?;
                let term = Term {
                    source: f.source,
                    target: g.target,
                    power: f.power + g.power,
                };
                sum.add(model, term, path.field.multiply(a, b))?;
            }
        }
    }
    Ok(sum.finish())
}

pub(crate) fn elements(model: &mut Model, grade: Grade, limit: usize) -> Result<Vec<Value>> {
    let path = model.jet_path(grade.diagram.path);
    let order = path.field.order();
    let max_rank = limit.ilog2() as usize / path.field.degree();
    let basis =
        model
            .jet_complex(grade.diagram)?
            .basis(model, grade.parity, 0..path.order, max_rank)?;
    let count = order.pow(basis.len() as u32);
    model.allocate(count + TERM_CELLS * basis.len() * (count - count / order))?;
    model.spend(count * (1 + basis.len()))?;
    Ok((0..count)
        .map(|mut digits| {
            Value::JetHom(
                basis
                    .iter()
                    .filter_map(|&term| {
                        let coefficient = digits % order;
                        digits /= order;
                        (coefficient != 0).then_some((term, coefficient))
                    })
                    .collect(),
            )
        })
        .collect())
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain("JetHom calculation exceeds the storage limit"));
    }
    Ok(())
}
