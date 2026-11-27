use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::syntax::ParsedExpr;

use super::local::LocalAlgebra;
use super::vacillating::{self, Path, Tail, Tower};
use super::{diagram, Algebra, Binary::Compose, Diagram, Unary::Reflect};

pub(crate) struct Key {
    pub(crate) rank: usize,
    pub(super) tail: Tail,
}

impl Key {
    pub(crate) fn new(rank: usize, tail: Vec<u128>) -> Result<Self> {
        vacillating::require_literal_rank(rank)?;
        let tail = vacillating::checked_tail(rank, tail)
            .filter(|tail| tail.iter().map(|&row| usize::from(row)).sum::<usize>() <= rank)
            .ok_or_else(|| {
                Error::plain("invalid Deligne terminal tail for the given tensor rank")
            })?;
        Ok(Self { rank, tail })
    }
}

pub(crate) enum Report {
    Paths(Key),
    Matrix(Key, Option<u128>, Box<ParsedExpr>),
    Overlaps(Key, Option<u128>),
}

pub(crate) struct Block {
    rank: usize,
    tail: Tail,
    paths: Vec<Path>,
    dimension: LocalAlgebra,
    bridges: RefCell<BTreeMap<usize, Rc<Bridge>>>,
    rotated: RefCell<BTreeMap<usize, LocalAlgebra>>,
    overlaps: RefCell<BTreeMap<(usize, usize), LocalAlgebra>>,
}

struct Bridge {
    forward: LocalAlgebra,
    inverse: LocalAlgebra,
}

impl Block {
    pub(super) fn new(model: &mut Model, tower: &Tower, key: &Key) -> Result<Self> {
        let paths = vacillating::terminal_paths(model, key.rank, &key.tail)?;
        let pivot = tower.projector(model, &paths[0])?;
        let dimension = pivot.trace_scalar(model)?;
        if dimension.is_zero() {
            return Err(Error::plain("Deligne block has zero generic dimension"));
        }
        dimension.charge(model)?;
        model.allocate(paths.len() * (2 * key.rank + 1) * (1 + key.rank))?;
        let bridge = Bridge {
            forward: pivot.clone(),
            inverse: pivot,
        };
        Ok(Self {
            rank: key.rank,
            tail: key.tail.clone(),
            paths,
            dimension,
            bridges: RefCell::new(BTreeMap::from([(0, Rc::new(bridge))])),
            rotated: RefCell::default(),
            overlaps: RefCell::default(),
        })
    }

    fn require_finite(&self, model: &mut Model, tower: &Tower, n: u128) -> Result<()> {
        if self.dimension.scalar_display(model, Some(n))? == "0" {
            return Err(Error::plain(format!(
                "Deligne terminal has dimension zero at t={n}"
            )));
        }
        for path in &self.paths {
            tower.projector(model, path)?.require_regular_at(n)?;
        }
        Ok(())
    }

    fn within(&self, indices: [u128; 2]) -> Result<[usize; 2]> {
        if indices.iter().any(|&i| i >= self.paths.len() as u128) {
            return Err(Error::plain("Deligne block path index is out of range"));
        }
        Ok(indices.map(|i| i as usize))
    }

    pub(crate) fn path_report(&self) -> String {
        let paths: Vec<_> = self.paths.iter().map(Path::vertices).collect();
        format!(
            "rvDeligneBlockPaths[{}, {:?}] = {paths:?}",
            self.rank, self.tail
        )
    }

    fn bridge(&self, model: &mut Model, tower: &Tower, i: usize) -> Result<Rc<Bridge>> {
        if let Some(bridge) = self.bridges.borrow().get(&i) {
            return Ok(bridge.clone());
        }
        let pivot = tower.projector(model, &self.paths[0])?;
        let target = tower.projector(model, &self.paths[i])?;
        let weight: usize = self.tail.iter().map(|&row| usize::from(row)).sum();
        let mut result = None;
        diagram::all_patterns(
            model,
            2 * self.rank,
            (2 * self.rank) as u128,
            &mut |model, pattern| {
                let diagram = Diagram::literal(self.rank, pattern.to_vec())?;
                // A diagram with fewer propagating blocks factors through a
                // smaller tensor power, where this terminal simple is absent.
                if diagram.propagating(self.rank) < weight {
                    return Ok(true);
                }
                let basis = LocalAlgebra::embed(Algebra {
                    rank: self.rank,
                    terms: BTreeMap::from([((diagram, 0), 1)]),
                })?;
                let on_pivot = basis.binary(model, Compose, &pivot)?;
                let u = target.binary(model, Compose, &on_pivot)?;
                if u.is_zero() {
                    return Ok(true);
                }
                result = Some(u);
                Ok(false)
            },
        )?;
        let forward = result
            .ok_or_else(|| Error::plain("no Deligne path bridge exists in the diagram basis"))?;
        let reflected = forward.unary(model, Reflect)?;
        let norm = reflected
            .trace_pair(model, &forward)?
            .divide_by_scalar(model, &self.dimension)?;
        let inverse = reflected.divide_by_scalar(model, &norm)?;
        if inverse.binary(model, Compose, &forward)? != pivot
            || forward.binary(model, Compose, &inverse)? != target
        {
            return Err(Error::plain(
                "Deligne path bridge failed its scalar corner identities",
            ));
        }
        forward.charge(model)?;
        inverse.charge(model)?;
        let bridge = Rc::new(Bridge { forward, inverse });
        self.bridges.borrow_mut().insert(i, bridge.clone());
        Ok(bridge)
    }

    pub(crate) fn unit(
        &self,
        model: &mut Model,
        tower: &Tower,
        indices: [u128; 2],
    ) -> Result<LocalAlgebra> {
        let [i, j] = self.within(indices)?;
        if i == j {
            return tower.projector(model, &self.paths[i]);
        }
        let target = self.bridge(model, tower, i)?;
        let source = self.bridge(model, tower, j)?;
        target.forward.binary(model, Compose, &source.inverse)
    }

    pub(crate) fn entry(
        &self,
        model: &mut Model,
        tower: &Tower,
        indices: [u128; 2],
        input: &LocalAlgebra,
    ) -> Result<LocalAlgebra> {
        let [i, j] = self.within(indices)?;
        if input.is_zero() {
            return Ok(LocalAlgebra::zero(0));
        }
        if i == j {
            let e = tower.projector(model, &self.paths[i])?;
            return e
                .trace_pair(model, input)?
                .divide_by_scalar(model, &self.dimension);
        }
        let right = self.bridge(model, tower, j)?;
        let column = input.binary(model, Compose, &right.forward)?;
        self.coefficient(model, tower, i, &column)
    }

    fn coefficient(
        &self,
        model: &mut Model,
        tower: &Tower,
        i: usize,
        column: &LocalAlgebra,
    ) -> Result<LocalAlgebra> {
        if column.is_zero() {
            return Ok(LocalAlgebra::zero(0));
        }
        let left = self.bridge(model, tower, i)?;
        left.inverse
            .trace_pair(model, column)?
            .divide_by_scalar(model, &self.dimension)
    }

    pub(crate) fn matrix(
        &self,
        model: &mut Model,
        tower: &Tower,
        input: &LocalAlgebra,
        n: Option<u128>,
    ) -> Result<String> {
        if let Some(n) = n {
            self.require_finite(model, tower, n)?;
            input.require_regular_at(n)?;
            if !input.is_zero() {
                for i in 0..self.paths.len() {
                    let bridge = self.bridge(model, tower, i)?;
                    bridge.forward.require_regular_at(n)?;
                    bridge.inverse.require_regular_at(n)?;
                }
            }
        }
        let size = self.paths.len();
        let mut matrix: Vec<Vec<String>> = (0..size).map(|_| Vec::with_capacity(size)).collect();
        for j in 0..size {
            let column = if input.is_zero() {
                LocalAlgebra::zero(self.rank)
            } else {
                let bridge = self.bridge(model, tower, j)?;
                input.binary(model, Compose, &bridge.forward)?
            };
            for (i, row) in matrix.iter_mut().enumerate() {
                row.push(
                    self.coefficient(model, tower, i, &column)?
                        .scalar_display(model, n)?,
                );
            }
        }
        Ok(self.matrix_report("rvDeligneBlock", &matrix, n))
    }

    pub(crate) fn overlap_entry(
        &self,
        model: &mut Model,
        tower: &Tower,
        indices: [u128; 2],
    ) -> Result<LocalAlgebra> {
        let [i, j] = self.within(indices)?;
        self.overlap(model, tower, i, j)
    }

    fn overlap(
        &self,
        model: &mut Model,
        tower: &Tower,
        i: usize,
        j: usize,
    ) -> Result<LocalAlgebra> {
        if let Some(overlap) = self.overlaps.borrow().get(&(i, j)) {
            return Ok(overlap.clone());
        }
        let e = tower.projector(model, &self.paths[i])?;
        let rotated = self.rotated.borrow().get(&j).cloned();
        let f = if let Some(rotated) = rotated {
            rotated
        } else {
            let rotated = tower.projector(model, &self.paths[j])?.rotate(model)?;
            rotated.charge(model)?;
            self.rotated.borrow_mut().insert(j, rotated.clone());
            rotated
        };
        let value = e
            .trace_pair(model, &f)?
            .divide_by_scalar(model, &self.dimension)?;
        value.charge(model)?;
        self.overlaps.borrow_mut().insert((i, j), value.clone());
        Ok(value)
    }

    pub(crate) fn overlap_matrix(
        &self,
        model: &mut Model,
        tower: &Tower,
        n: Option<u128>,
    ) -> Result<String> {
        if let Some(n) = n {
            self.require_finite(model, tower, n)?;
        }
        let size = self.paths.len();
        let mut matrix = Vec::with_capacity(size);
        let mut columns = vec![LocalAlgebra::zero(0); size];
        let one = LocalAlgebra::identity(0);
        for i in 0..size {
            let mut row = Vec::with_capacity(size);
            let mut sum = LocalAlgebra::zero(0);
            for (j, column) in columns.iter_mut().enumerate() {
                let value = self.overlap(model, tower, i, j)?;
                sum = sum.add(model, &value)?;
                *column = column.add(model, &value)?;
                row.push(value.scalar_display(model, n)?);
            }
            if sum != one {
                return Err(Error::plain("Deligne recoupling row does not sum to one"));
            }
            matrix.push(row);
        }
        if columns.iter().any(|col| *col != one) {
            return Err(Error::plain(
                "Deligne recoupling column does not sum to one",
            ));
        }
        Ok(self.matrix_report("rvDeligneRacah", &matrix, n))
    }

    fn matrix_report(&self, name: &str, matrix: &[Vec<String>], n: Option<u128>) -> String {
        let rows: Vec<_> = matrix
            .iter()
            .map(|row| format!("[{}]", row.join(", ")))
            .collect();
        let dimension = n.map_or_else(String::new, |n| format!(", {n}"));
        let suffix = if n.is_some() { "At" } else { "" };
        format!(
            "{name}{suffix}[{}, {:?}{dimension}] = [{}]",
            self.rank,
            self.tail,
            rows.join(", ")
        )
    }
}
