use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::integer_smith::{self, Matrix, Track};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::stable_marks::{self, Tree};

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) prime: Prime,
    pub(crate) marks: usize,
}

impl Grade {
    pub(crate) fn new(prime: Prime, marks: usize) -> Result<Self> {
        if prime.get() == 2 {
            return Err(Error::plain(
                "tame hyperelliptic cover requires odd residue characteristic",
            ));
        }
        if marks < 4 || !marks.is_multiple_of(2) {
            return Err(Error::plain(
                "hyperelliptic cover requires an even number of at least four branch marks",
            ));
        }
        stable_marks::check_count(marks)?;
        Ok(Self { prime, marks })
    }

    fn genus(self) -> usize {
        self.marks / 2 - 1
    }
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
struct Vertex {
    base: usize,
    genus: usize,
    degree: u8,
    deck: usize,
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
struct Edge {
    from: usize,
    to: usize,
    base: usize,
    length: u128,
    dilation: u8,
    deck: usize,
}

impl Edge {
    fn fraction(self) -> (u128, u8) {
        if self.dilation == 2 && self.length.is_multiple_of(2) {
            (self.length / 2, 1)
        } else {
            (self.length, self.dilation)
        }
    }
}

#[derive(Clone, Copy)]
struct Fiber {
    first: usize,
    size: usize,
}

impl Fiber {
    fn at(self, sheet: usize) -> usize {
        self.first + sheet % self.size
    }
}

type Cycle = Vec<(usize, i32)>;

pub(crate) struct Cover {
    grade: Grade,
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    branches: Vec<usize>,
    cycles: Vec<Cycle>,
}

struct Base {
    splits: Vec<([usize; 2], u128, bool)>,
    branches: Vec<usize>,
    directions: Vec<Vec<usize>>,
}

impl Base {
    fn read(model: &mut Model, tree: &Tree) -> Result<Self> {
        let mut locations = BTreeMap::new();
        for (v, directions) in tree.vertices().iter().enumerate() {
            for direction in directions {
                model.spend(1 + direction.len() * (1 + locations.len().max(1).ilog2() as usize))?;
                if locations.insert(direction.as_slice(), v).is_some() {
                    return Err(Error::plain("stable branch tree repeats a direction"));
                }
            }
        }

        let mut directions = vec![Vec::new(); tree.vertices().len()];
        let mut splits = Vec::new();
        for (i, (a, b, length)) in tree.splits().iter().enumerate() {
            model.spend(1 + tree.marks)?;
            let Some((&first, &second)) =
                locations.get(a.as_slice()).zip(locations.get(b.as_slice()))
            else {
                return Err(Error::plain("stable branch split is missing an endpoint"));
            };
            if first == second || *length == 0 {
                return Err(Error::plain(
                    "invalid finite edge in the stable branch tree",
                ));
            }
            directions[first].push(i);
            directions[second].push(i);
            let endpoints = [first.min(second), first.max(second)];
            splits.push((endpoints, *length, !a.len().is_multiple_of(2)));
        }

        let mut branches = Vec::with_capacity(tree.marks);
        for mark in 0..tree.marks {
            model.spend(1 + locations.len().max(1).ilog2() as usize)?;
            let Some(&vertex) = locations.get(&[mark][..]) else {
                return Err(Error::plain("stable branch mark has no incident vertex"));
            };
            directions[vertex].push(splits.len() + mark);
            branches.push(vertex);
        }
        if locations.len() != 2 * splits.len() + branches.len() {
            return Err(Error::plain(
                "stable branch tree direction count disagrees with incidences",
            ));
        }
        Ok(Self {
            splits,
            branches,
            directions,
        })
    }
}

impl Cover {
    pub(crate) fn build(model: &mut Model, grade: Grade, tree: &Tree) -> Result<Rc<Self>> {
        debug_assert_eq!(grade.marks, tree.marks);
        let base = Base::read(model, tree)?;
        let mut vertices = Vec::new();
        let mut fibers = Vec::new();
        for (v, directions) in tree.vertices().iter().enumerate() {
            model.spend(1 + directions.len())?;
            let odd = directions
                .iter()
                .filter(|d| !d.len().is_multiple_of(2))
                .count();
            if !odd.is_multiple_of(2) {
                return Err(Error::plain(
                    "hyperelliptic vertex has an odd number of ramified directions",
                ));
            }
            let first = vertices.len();
            let size = if odd == 0 { 2 } else { 1 };
            fibers.push(Fiber { first, size });
            for sheet in 0..size {
                vertices.push(Vertex {
                    base: v,
                    genus: odd.saturating_sub(2) / 2,
                    degree: (2 / size) as u8,
                    deck: first + (size - 1 - sheet),
                });
            }
        }
        let mut edges = Vec::new();
        for (i, &([from, to], length, odd)) in base.splits.iter().enumerate() {
            let count = if odd { 1 } else { 2 };
            let first = edges.len();
            for sheet in 0..count {
                model.spend(1)?;
                edges.push(Edge {
                    from: fibers[from].at(sheet),
                    to: fibers[to].at(sheet),
                    base: i,
                    length,
                    dilation: (2 / count) as u8,
                    deck: first + (count - 1 - sheet),
                });
            }
        }
        model.spend(1 + base.branches.len())?;
        let branches = base.branches.iter().map(|&v| fibers[v].first).collect();
        let mut result = Self {
            grade,
            vertices,
            edges,
            branches,
            cycles: Vec::new(),
        };
        result.certify_map(model, &base)?;
        result.cycles = fundamental_cycles(model, result.vertices.len(), &result.edges)?;
        if result.cycles.len() + result.vertices.iter().map(|v| v.genus).sum::<usize>()
            != grade.genus()
        {
            return Err(Error::plain(
                "hyperelliptic graph genus disagrees with the branch count",
            ));
        }
        let cells = result.cells();
        integer_smith::storage(cells)?;
        model.allocate(cells)?;
        Ok(Rc::new(result))
    }

    fn certify_map(&self, model: &mut Model, base: &Base) -> Result<()> {
        let mut degrees = vec![0usize; base.directions.len()];
        let mut local = vec![BTreeMap::<usize, usize>::new(); self.vertices.len()];
        for (i, v) in self.vertices.iter().enumerate() {
            model.spend(1)?;
            let deck = self.vertices[v.deck];
            if deck.deck != i
                || deck.base != v.base
                || deck.genus != v.genus
                || deck.degree != v.degree
            {
                return Err(Error::plain(
                    "invalid hyperelliptic deck action on vertices",
                ));
            }
            degrees[v.base] += usize::from(v.degree);
        }
        for (i, e) in self.edges.iter().enumerate() {
            model.spend(
                1 + local[e.from].len().max(1).ilog2() as usize
                    + local[e.to].len().max(1).ilog2() as usize,
            )?;
            let [a, b] = base.splits[e.base].0;
            let expected = e.dilation == if base.splits[e.base].2 { 2 } else { 1 };
            if !expected || self.vertices[e.from].base != a || self.vertices[e.to].base != b {
                return Err(Error::plain(
                    "hyperelliptic edge has the wrong branch-tree image",
                ));
            }
            let partner = self.edges[e.deck];
            if partner.deck != i
                || partner.from != self.vertices[e.from].deck
                || partner.to != self.vertices[e.to].deck
                || partner.base != e.base
                || partner.length != e.length
                || partner.dilation != e.dilation
            {
                return Err(Error::plain("invalid hyperelliptic deck action on edges"));
            }
            *local[e.from].entry(e.base).or_insert(0) += usize::from(e.dilation);
            *local[e.to].entry(e.base).or_insert(0) += usize::from(e.dilation);
        }
        for (i, &v) in self.branches.iter().enumerate() {
            model.spend(1 + local[v].len().max(1).ilog2() as usize)?;
            if self.vertices[v].base != base.branches[i] || self.vertices[v].deck != v {
                return Err(Error::plain("invalid hyperelliptic ramified mark"));
            }
            *local[v].entry(base.splits.len() + i).or_insert(0) += 2;
        }
        if degrees.iter().any(|&d| d != 2) {
            return Err(Error::plain(
                "hyperelliptic graph has the wrong degree over a base vertex",
            ));
        }
        for (vertex, dilation) in self.vertices.iter().zip(&local) {
            let directions = &base.directions[vertex.base];
            model.spend(1 + directions.len() * (1 + dilation.len().max(1).ilog2() as usize))?;
            if dilation.len() != directions.len()
                || directions
                    .iter()
                    .any(|d| dilation.get(d) != Some(&usize::from(vertex.degree)))
            {
                return Err(Error::plain(
                    "hyperelliptic graph map fails harmonicity in a branch direction",
                ));
            }
        }
        Ok(())
    }

    fn cells(&self) -> usize {
        8 + 5 * self.vertices.len()
            + 8 * self.edges.len()
            + self.branches.len()
            + self
                .cycles
                .iter()
                .map(|cycle| 3 + 2 * cycle.len())
                .sum::<usize>()
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(self.cells().min(rhs.cells()))?;
        Ok((&self.vertices, &self.edges, &self.branches).cmp(&(
            &rhs.vertices,
            &rhs.edges,
            &rhs.branches,
        )))
    }

    fn metric(&self, model: &mut Model, multiplier: usize) -> Result<Vec<i128>> {
        if multiplier == 0 {
            return Err(Error::plain(
                "hyperelliptic ramification multiplier must be positive",
            ));
        }
        self.edges
            .iter()
            .enumerate()
            .map(|(i, &edge)| {
                model.spend(1)?;
                let multiplier = multiplier as u128;
                let (length, denominator) = edge.fraction();
                let d = u128::from(denominator);
                if !multiplier.is_multiple_of(d) {
                    return Err(Error::plain(format!(
                    "hyperelliptic edge {i} has half-integral length at multiplier {multiplier}",
                )));
                }
                length
                    .checked_mul(multiplier / d)
                    .and_then(|length| i128::try_from(length).ok())
                    .ok_or_else(metric_overflow)
            })
            .collect()
    }

    fn monodromy(&self, model: &mut Model, multiplier: usize) -> Result<Vec<Vec<i128>>> {
        let lengths = self.metric(model, multiplier)?;
        let rank = self.cycles.len();
        integer_smith::storage(rank * (3 + 2 * rank))?;
        let mut pairing = vec![vec![0; rank]; rank];
        for (i, a) in self.cycles.iter().enumerate() {
            for (j, b) in self.cycles.iter().enumerate().skip(i) {
                let value = period_inner(model, a, b, &lengths)?;
                pairing[i][j] = value;
                pairing[j][i] = value;
            }
        }
        Ok(pairing)
    }

    fn component_group(&self, model: &mut Model, multiplier: usize) -> Result<String> {
        let pairing = self.monodromy(model, multiplier)?;
        let rank = pairing.len();
        let mut matrix = Matrix::new(model, rank, rank)?;
        for (i, row) in pairing.into_iter().enumerate() {
            for (j, value) in row.into_iter().enumerate().filter(|&(_, x)| x != 0) {
                matrix.write(model, i, j, value)?;
            }
        }
        let smith = integer_smith::reduce(model, matrix, Track::None, Track::None)?;
        if smith.factors.len() != rank {
            return Err(Error::plain(
                "integral hyperelliptic monodromy pairing is degenerate",
            ));
        }
        let exponent = smith.factors.last().copied().unwrap_or(1);
        let mut order = 1i128;
        let mut elementary_divisors = Vec::new();
        for value in smith.factors {
            model.spend(1)?;
            order = order.checked_mul(value).ok_or_else(|| {
                Error::plain("hyperelliptic component-group order exceeds signed i128")
            })?;
            if value > 1 {
                elementary_divisors.push(value);
            }
        }
        Ok(format!(
            "rvHyperComponent = factors={elementary_divisors:?}, exponent={exponent}, order={order}",
        ))
    }

    fn graph(&self, model: &mut Model) -> Result<String> {
        model.spend(self.cells())?;
        let vertices: Vec<_> = self
            .vertices
            .iter()
            .map(|v| (v.base, v.genus, v.degree, v.deck))
            .collect();
        let edges: Vec<_> = self
            .edges
            .iter()
            .map(|e| (e.from, e.to, e.base, e.fraction(), e.dilation, e.deck))
            .collect();
        let branches: Vec<_> = self
            .branches
            .iter()
            .enumerate()
            .map(|(i, &v)| (i, v, 2))
            .collect();
        Ok(format!(
            "rvHyperGraph = vertices={vertices:?}, edges={edges:?}, branches={branches:?}"
        ))
    }

    pub(crate) fn summary(&self, model: &mut Model) -> Result<String> {
        model.spend(1)?;
        Ok(format!(
            "RVHyperCover({}, {}){{genus={}, toric={}}}",
            self.grade.prime.get(),
            self.grade.marks,
            self.grade.genus(),
            self.cycles.len()
        ))
    }
}

fn fundamental_cycles(model: &mut Model, count: usize, edges: &[Edge]) -> Result<Vec<Cycle>> {
    if count == 0 {
        return Err(Error::plain("hyperelliptic graph has no vertices"));
    }
    let mut incident = vec![Vec::new(); count];
    for (i, e) in edges.iter().enumerate() {
        model.spend(1)?;
        incident[e.from].push((i, e.to));
        incident[e.to].push((i, e.from));
    }
    let mut parent = vec![None; count];
    let mut reached = vec![false; count];
    let mut in_tree = vec![false; edges.len()];
    let mut queue = vec![0];
    reached[0] = true;
    let mut i = 0;
    while i < queue.len() {
        let v = queue[i];
        model.spend(1 + incident[v].len())?;
        for &(edge, w) in &incident[v] {
            if !reached[w] {
                reached[w] = true;
                parent[w] = Some((v, edge));
                in_tree[edge] = true;
                queue.push(w);
            }
        }
        i += 1;
    }
    if queue.len() != count {
        return Err(Error::plain("hyperelliptic graph is disconnected"));
    }
    let mut basis = Vec::new();
    let mut cells = 0;
    for (edge, e) in edges.iter().enumerate().filter(|(i, _)| !in_tree[*i]) {
        let mut cycle = BTreeMap::from([(edge, 1)]);
        root_path(model, &parent, edges, &mut cycle, e.to, 1)?;
        root_path(model, &parent, edges, &mut cycle, e.from, -1)?;
        model.spend(count + cycle.len())?;
        let row: Cycle = cycle
            .into_iter()
            .filter(|&(_, coefficient)| coefficient != 0)
            .collect();
        let mut boundary = vec![0i32; count];
        for &(i, c) in &row {
            boundary[edges[i].from] -= c;
            boundary[edges[i].to] += c;
        }
        if boundary.iter().any(|&c| c != 0) {
            return Err(Error::plain(
                "constructed hyperelliptic period cycle has nonzero boundary",
            ));
        }
        cells += 3 + 2 * row.len();
        integer_smith::storage(cells)?;
        basis.push(row);
    }
    Ok(basis)
}

fn root_path(
    model: &mut Model,
    parents: &[Option<(usize, usize)>],
    edges: &[Edge],
    cycle: &mut BTreeMap<usize, i32>,
    mut child: usize,
    sign: i32,
) -> Result<()> {
    while let Some((parent, edge)) = parents[child] {
        model.spend(1 + cycle.len().max(1).ilog2() as usize)?;
        let orientation = if edges[edge].from == child {
            sign
        } else {
            -sign
        };
        *cycle.entry(edge).or_insert(0) += orientation;
        child = parent;
    }
    Ok(())
}

fn period_inner(model: &mut Model, a: &Cycle, b: &Cycle, lengths: &[i128]) -> Result<i128> {
    model.spend(1 + a.len() + b.len())?;
    let (mut i, mut j) = (0, 0);
    let mut sum = 0i128;
    while i < a.len() && j < b.len() {
        let (k, x) = a[i];
        let (l, y) = b[j];
        match k.cmp(&l) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                let length = lengths[k]
                    .checked_mul(i128::from(x * y))
                    .ok_or_else(metric_overflow)?;
                sum = sum.checked_add(length).ok_or_else(metric_overflow)?;
                i += 1;
                j += 1;
            }
        }
    }
    Ok(sum)
}

fn metric_overflow() -> Error {
    Error::plain("hyperelliptic edge length or monodromy pairing exceeds signed i128")
}

pub(crate) enum Query {
    Graph,
    Harmonic,
    Genus,
    Cycles,
    Metric(usize),
    Monodromy(usize),
    Component(usize),
}

impl Query {
    pub(crate) fn report(self, model: &mut Model, cover: &Cover) -> Result<String> {
        match self {
            Self::Graph => cover.graph(model),
            Self::Harmonic => Ok("rvHyperHarmonic = true".to_owned()),
            Self::Genus => {
                model.spend(cover.vertices.len())?;
                let genera: Vec<_> = cover.vertices.iter().map(|v| v.genus).collect();
                Ok(format!(
                    "rvHyperGenus = {}, toric={}, components={genera:?}",
                    cover.grade.genus(),
                    cover.cycles.len()
                ))
            }
            Self::Cycles => {
                model.spend(cover.cycles.iter().map(|row| 1 + row.len()).sum())?;
                Ok(format!("rvHyperCycles = {:?}", cover.cycles))
            }
            Self::Metric(multiplier) => Ok(format!(
                "rvHyperMetric = {:?}",
                cover.metric(model, multiplier)?
            )),
            Self::Monodromy(multiplier) => Ok(format!(
                "rvHyperMonodromy = {:?}",
                cover.monodromy(model, multiplier)?
            )),
            Self::Component(multiplier) => cover.component_group(model, multiplier),
        }
    }
}
