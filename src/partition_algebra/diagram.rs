use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;

use super::dimension::dimensions;
use super::{coefficient_overflow, natural, Monomial, MAX_RANK, MAX_TERMS};

const MAX_VERTICES: usize = 2 * MAX_RANK;
const NONE: u8 = u8::MAX;

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(super) struct Diagram(Box<[u8]>);

impl Diagram {
    pub(super) fn literal(rank: usize, labels: Vec<u128>) -> Result<Self> {
        if labels.len() != 2 * rank {
            return Err(Error::plain(
                "partition diagram must label exactly 2r vertices",
            ));
        }
        let mut blocks = BTreeMap::new();
        let labels = labels
            .into_iter()
            .map(|label| {
                let next = blocks.len() as u8;
                *blocks.entry(label).or_insert(next)
            })
            .collect();
        Ok(Self(labels))
    }

    pub(super) fn identity(rank: usize) -> Self {
        let row = 0..rank as u8;
        Self(row.clone().chain(row).collect())
    }

    pub(super) fn labels(&self) -> &[u8] {
        &self.0
    }

    fn blocks(&self) -> usize {
        self.0
            .iter()
            .max()
            .map_or(0, |&label| usize::from(label) + 1)
    }

    pub(super) fn reflect(&self, rank: usize) -> Self {
        Self::canonical(self.0[rank..].iter().chain(&self.0[..rank]).copied())
    }

    pub(super) fn tensor(&self, left: usize, rhs: &Self, right: usize) -> Self {
        let offset = self.blocks() as u8;
        Self::canonical(
            self.0[..left]
                .iter()
                .copied()
                .chain(rhs.0[..right].iter().map(|&label| label + offset))
                .chain(self.0[left..].iter().copied())
                .chain(rhs.0[right..].iter().map(|&label| label + offset)),
        )
    }

    pub(super) fn compose(&self, rhs: &Self, rank: usize) -> (Self, usize) {
        let mut graph = Connectivity::new(3 * rank);
        graph.install(self, 0);
        graph.install(rhs, rank);
        let mut visible = [false; 3 * MAX_RANK];
        let exterior = (0..rank).chain(2 * rank..3 * rank);
        let diagram = Self::canonical(exterior.map(|vertex| {
            let root = graph.root(vertex);
            visible[root] = true;
            root as u8
        }));
        let mut internal = 0;
        for vertex in rank..2 * rank {
            let root = graph.root(vertex);
            if !visible[root] {
                visible[root] = true;
                internal += 1;
            }
        }
        (diagram, internal)
    }

    pub(super) fn closure_components(&self, rank: usize) -> usize {
        let mut graph = Connectivity::new(2 * rank);
        graph.install(self, 0);
        for i in 0..rank {
            graph.join(i, i + rank);
        }
        (0..2 * rank)
            .filter(|&i| graph.parents[i] as usize == i)
            .count()
    }

    pub(super) fn supports(&self, indices: &[u128]) -> bool {
        let mut seen = [None; MAX_VERTICES];
        for (&block, &index) in self.0.iter().zip(indices) {
            let slot = &mut seen[block as usize];
            if *slot.get_or_insert(index) != index {
                return false;
            }
        }
        true
    }

    pub(super) fn orbit(&self, model: &mut Model) -> Result<BTreeMap<Monomial, i128>> {
        let blocks = self.blocks();
        let (bell, _) = dimensions(model, blocks, 0)?;
        if bell
            .magnitude
            .cmp(natural(MAX_TERMS as u128).magnitude)
            .is_gt()
        {
            return Err(Error::plain(
                "partition orbit coarsenings exceed the term storage limit",
            ));
        }
        let mut result = BTreeMap::new();
        self.coarsen(model, blocks, &mut Vec::new(), &mut Vec::new(), &mut result)?;
        Ok(result)
    }

    fn coarsen(
        &self,
        model: &mut Model,
        blocks: usize,
        grouped: &mut Vec<u8>,
        sizes: &mut Vec<usize>,
        result: &mut BTreeMap<Monomial, i128>,
    ) -> Result<()> {
        if grouped.len() == blocks {
            model.spend(1 + self.0.len() + blocks)?;
            let mut coefficient = 1i128;
            for &size in sizes.iter() {
                for factor in 1..size {
                    coefficient = coefficient
                        .checked_mul(-(factor as i128))
                        .ok_or_else(coefficient_overflow)?;
                }
            }
            let diagram = Self(
                self.0
                    .iter()
                    .map(|&label| grouped[label as usize])
                    .collect(),
            );
            result.insert((diagram, 0), coefficient);
            return Ok(());
        }
        let existing = sizes.len();
        for target in 0..=existing {
            model.spend(1)?;
            if target == existing {
                sizes.push(0);
            }
            sizes[target] += 1;
            grouped.push(target as u8);
            self.coarsen(model, blocks, grouped, sizes, result)?;
            grouped.pop();
            sizes[target] -= 1;
            if target == existing {
                sizes.pop();
            }
        }
        Ok(())
    }

    fn canonical(labels: impl Iterator<Item = u8>) -> Self {
        let mut seen = [NONE; 3 * MAX_RANK];
        let mut next = 0;
        let labels = labels
            .map(|label| {
                let slot = &mut seen[label as usize];
                if *slot == NONE {
                    *slot = next;
                    next += 1;
                }
                *slot
            })
            .collect();
        Self(labels)
    }
}

struct Connectivity {
    parents: [u8; 3 * MAX_RANK],
}

impl Connectivity {
    fn new(size: usize) -> Self {
        let mut parents = [0; 3 * MAX_RANK];
        for (i, slot) in parents.iter_mut().take(size).enumerate() {
            *slot = i as u8;
        }
        Self { parents }
    }

    fn root(&mut self, mut vertex: usize) -> usize {
        while self.parents[vertex] as usize != vertex {
            let parent = self.parents[vertex] as usize;
            self.parents[vertex] = self.parents[parent];
            vertex = self.parents[vertex] as usize;
        }
        vertex
    }

    fn join(&mut self, a: usize, b: usize) {
        let (a, b) = (self.root(a), self.root(b));
        self.parents[a.max(b)] = a.min(b) as u8;
    }

    fn install(&mut self, diagram: &Diagram, offset: usize) {
        let mut first = [None; MAX_VERTICES];
        for (vertex, &label) in diagram.0.iter().enumerate() {
            let vertex = vertex + offset;
            let slot = &mut first[label as usize];
            self.join(*slot.get_or_insert(vertex), vertex);
        }
    }
}

// S_N-orbits of pairs of tensors are the realizable restricted growth
// words. Short-circuit as soon as one orbit has a nonzero entry.
pub(super) fn all_patterns(
    model: &mut Model,
    vertices: usize,
    dimension: u128,
    test: &mut impl FnMut(&mut Model, &[u128]) -> Result<bool>,
) -> Result<bool> {
    fn visit(
        model: &mut Model,
        vertices: usize,
        dimension: u128,
        used: u128,
        labels: &mut Vec<u128>,
        test: &mut impl FnMut(&mut Model, &[u128]) -> Result<bool>,
    ) -> Result<bool> {
        if labels.len() == vertices {
            return test(model, labels);
        }
        for label in 0..(used + 1).min(dimension) {
            model.spend(1)?;
            labels.push(label);
            if !visit(
                model,
                vertices,
                dimension,
                used.max(label + 1),
                labels,
                test,
            )? {
                return Ok(false);
            }
            labels.pop();
        }
        Ok(true)
    }
    visit(model, vertices, dimension, 0, &mut Vec::new(), test)
}
