use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::integer_smith::{self as integral, Matrix, Track, Vector};
use crate::model::{Model, Value};
use crate::rv_lens::{self, Rational};
use crate::rv_tl;
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Command {
    Braid(Label, usize, Vec<i128>),
    Group(Label, Grade),
    Euler(Label),
}

pub(crate) enum Query {
    Term(Label, Vec<usize>, Vec<usize>, i128),
    Generator(Label, Grade, usize, Generator),
    EulerAt(Label, usize),
    Operation(Operation, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Generator {
    Basis,
    Class,
    Null,
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Differential,
    Project,
    Representative,
    Scale(i128),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Differential => "khD",
            Self::Project => "khClass",
            Self::Representative => "khRep",
            Self::Scale(_) => "khScale",
        }
    }
}

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) enum Form {
    Chain,
    Class,
}

impl Form {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Chain => "KhC",
            Self::Class => "KhH",
        }
    }
}

type State = u128;
pub(crate) type Chain = Rc<[(Cell, i128)]>;
pub(crate) type Class = Rc<[(usize, i128)]>;

#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Cell {
    pub(crate) state: State,
    pub(crate) xs: Rc<[usize]>,
}

impl Cell {
    pub(crate) fn cells(&self) -> usize {
        7 + self.xs.len()
    }
}

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) i: i64,
    pub(crate) j: i64,
}

impl Grade {
    pub(crate) fn shift(self, delta: i64) -> Result<Self> {
        let i = self
            .i
            .checked_add(delta)
            .ok_or_else(|| Error::plain("Khovanov homological degree is out of range"))?;
        Ok(Self { i, j: self.j })
    }
}

#[derive(Clone, Copy)]
struct Crossing {
    negative: bool,
    before: [usize; 2],
    after: [usize; 2],
}

pub(crate) struct Braid {
    crossings: Vec<Crossing>,
    ends: Vec<usize>,
    negative: usize,
    groups: RefCell<HashMap<Grade, Rc<Group>>>,
}

impl Braid {
    pub(crate) fn new(model: &mut Model, strands: usize, word: &[i128]) -> Result<Self> {
        if word.len() > State::BITS as usize {
            return Err(Error::plain("Khovanov braid has more than 128 crossings"));
        }
        let vertices = strands
            .checked_add(2 * word.len())
            .ok_or_else(integral::arithmetic)?;
        integral::storage(vertices)?;
        model.spend(vertices)?;
        let mut ends: Vec<_> = (0..strands).collect();
        let mut crossings = Vec::with_capacity(word.len());
        let mut negative = 0;
        // Like rvTLBraid, the input word acts from right to left. Index
        // resolutions by positions in the original word.
        for (r, &letter) in word.iter().enumerate().rev() {
            let index = letter.unsigned_abs();
            if index == 0 || index >= strands as u128 {
                return Err(Error::plain(
                    "Khovanov braid letter must have index between 1 and strands minus one",
                ));
            }
            let left = index as usize - 1;
            let right = left + 1;
            let before = [ends[left], ends[right]];
            let first = strands + 2 * (word.len() - 1 - r);
            let after = [first, first + 1];
            ends[left..=right].copy_from_slice(&after);
            negative += usize::from(letter < 0);
            crossings.push(Crossing {
                negative: letter < 0,
                before,
                after,
            });
        }
        crossings.reverse();
        model.allocate(ends.len() + 5 * crossings.len() + 1)?;
        Ok(Self {
            crossings,
            ends,
            negative,
            groups: RefCell::default(),
        })
    }

    pub(crate) fn group(&self, model: &mut Model, grade: Grade) -> Result<Rc<Group>> {
        if let Some(group) = self.groups.borrow().get(&grade) {
            return Ok(group.clone());
        }
        let group = Rc::new(Group::build(model, self, grade)?);
        self.groups.borrow_mut().insert(grade, group.clone());
        Ok(group)
    }

    pub(crate) fn invariants(&self, model: &mut Model, grade: Grade) -> Result<Invariants> {
        if let Some(group) = self.groups.borrow().get(&grade) {
            model.spend(group.invariants.torsion.len())?;
            return Ok(group.invariants.clone());
        }
        let mut context = Context::new(self);
        let basis = context.basis(model, grade)?;
        let chains = basis.len();
        let outgoing = incident(model, &mut context, &basis, false)?;
        let outgoing = incidence_smith(model, chains, outgoing, Track::None)?;
        let incoming = incident(model, &mut context, &basis, true)?;
        // The reversed incidence rows give the transpose of the incoming
        // differential, with the same Smith factors. The image of the
        // outgoing differential is free, so its kernel splits integrally.
        let incoming = incidence_smith(model, chains, incoming, Track::None)?;
        let free = chains
            .checked_sub(outgoing.factors.len() + incoming.factors.len())
            .ok_or_else(|| Error::plain("Khovanov differential does not square to zero"))?;
        let torsion = incoming.factors.into_iter().filter(|&m| m > 1).collect();
        Ok(Invariants {
            chains,
            free,
            torsion,
        })
    }

    pub(crate) fn basis_term(
        &self,
        model: &mut Model,
        grade: Grade,
        index: usize,
    ) -> Result<Chain> {
        let cell = Context::new(self)
            .basis(model, grade)?
            .into_iter()
            .nth(index)
            .ok_or_else(|| Error::plain("Khovanov basis index is out of range"))?;
        persist(model, vec![(cell, 1)])
    }

    pub(crate) fn chain_elements(&self, model: &mut Model, grade: Grade) -> Result<Vec<Value>> {
        if !Context::new(self).basis(model, grade)?.is_empty() {
            return Err(Error::plain(
                "cannot enumerate a nonzero free Khovanov chain group",
            ));
        }
        model.allocate(1)?;
        Ok(vec![Value::KhChain(Rc::default())])
    }

    pub(crate) fn possible(&self, grade: Grade) -> bool {
        self.grading(grade).is_some()
    }

    pub(crate) fn components(&self, model: &mut Model) -> Result<usize> {
        model.spend(4 * self.vertices())?;
        integral::storage(2 * self.vertices())?;
        let mut components = Disjoint::new(self.vertices());
        for crossing in &self.crossings {
            components.join(crossing.before[0], crossing.after[1]);
            components.join(crossing.before[1], crossing.after[0]);
        }
        for (i, &end) in self.ends.iter().enumerate() {
            components.join(i, end);
        }
        Ok((0..self.vertices())
            .filter(|&v| components.root(v) == v)
            .count())
    }

    pub(crate) fn lee_possible(&self, i: i64, f: i64) -> bool {
        let ones = i128::from(i) + self.negative as i128;
        (0..=self.crossings.len() as i128).contains(&ones)
            && f <= self.shift(ones as usize) + self.vertices() as i64
    }

    pub(crate) fn lee_lower_bound(&self) -> i64 {
        self.shift(0) - self.vertices() as i64
    }

    fn vertices(&self) -> usize {
        self.ends.len() + 2 * self.crossings.len()
    }

    fn shift(&self, ones: usize) -> i64 {
        ones as i64 + self.crossings.len() as i64 - 3 * self.negative as i64
    }

    fn grade(&self, cell: &Cell, circles: usize) -> Grade {
        let ones = cell.state.count_ones() as usize;
        Grade {
            i: ones as i64 - self.negative as i64,
            j: self.shift(ones) + circles as i64 - 2 * cell.xs.len() as i64,
        }
    }

    fn grading(&self, grade: Grade) -> Option<(usize, i64)> {
        let ones = i128::from(grade.i) + self.negative as i128;
        if ones < 0 || ones > self.crossings.len() as i128 {
            return None;
        }
        // Each smoothing change switches the parity of the circle count.
        // Starting from the all-vertical state gives j ≡ strands + crossings.
        if grade.j.rem_euclid(2) as usize != (self.ends.len() + self.crossings.len()) % 2 {
            return None;
        }
        let ones = ones as usize;
        let positive = self.crossings.len() - self.negative;
        let horizontal = ones + self.negative - 2 * ones.saturating_sub(positive);
        let bound = (self.ends.len() + horizontal) as u128;
        let labels = i128::from(grade.j) - i128::from(self.shift(ones));
        (labels.unsigned_abs() <= bound).then_some((ones, labels as i64))
    }

    pub(crate) fn term(
        &self,
        model: &mut Model,
        ones: &[usize],
        xs: Vec<usize>,
        coefficient: i128,
    ) -> Result<(Grade, Chain)> {
        model.spend(ones.len() + xs.len())?;
        let mut state = 0;
        for &i in ones {
            if i >= self.crossings.len() {
                return Err(Error::plain("Khovanov crossing index is out of range"));
            }
            let bit = 1u128 << i;
            if state & bit != 0 {
                return Err(Error::plain("repeated Khovanov crossing index"));
            }
            state |= bit;
        }
        let resolution = self.resolve(model, state)?;
        if xs.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(Error::plain(
                "Khovanov circle indices must be strictly increasing",
            ));
        }
        if xs.last().is_some_and(|&i| i >= resolution.anchors.len()) {
            return Err(Error::plain("Khovanov circle index is out of range"));
        }
        let cell = Cell {
            state,
            xs: xs.into(),
        };
        let grade = self.grade(&cell, resolution.anchors.len());
        let chain = if coefficient == 0 {
            Rc::default()
        } else {
            persist(model, vec![(cell, coefficient)])?
        };
        Ok((grade, chain))
    }

    fn resolve(&self, model: &mut Model, state: State) -> Result<Resolution> {
        integral::storage(5 * self.vertices())?;
        model.spend(4 * self.vertices())?;
        let mut components = Disjoint::new(self.vertices());
        for (i, edge) in self.crossings.iter().enumerate() {
            let horizontal = (state & (1 << i) != 0) != edge.negative;
            if horizontal {
                components.join(edge.before[0], edge.before[1]);
                components.join(edge.after[0], edge.after[1]);
            } else {
                components.join(edge.before[0], edge.after[0]);
                components.join(edge.before[1], edge.after[1]);
            }
        }
        for (i, &end) in self.ends.iter().enumerate() {
            components.join(i, end);
        }
        let mut root_label = vec![usize::MAX; self.vertices()];
        let mut labels = Vec::with_capacity(self.vertices());
        let mut anchors = Vec::new();
        for vertex in 0..self.vertices() {
            let root = components.root(vertex);
            if root_label[root] == usize::MAX {
                root_label[root] = anchors.len();
                anchors.push(vertex);
            }
            labels.push(root_label[root]);
        }
        Ok(Resolution { labels, anchors })
    }

    pub(crate) fn differential(&self, model: &mut Model, input: &Chain) -> Result<Chain> {
        self.integral_differential(model, input, Saddle::Khovanov)
    }

    pub(crate) fn lee_perturbation(&self, model: &mut Model, input: &Chain) -> Result<Chain> {
        self.integral_differential(model, input, Saddle::Perturbation)
    }

    fn integral_differential(
        &self,
        model: &mut Model,
        input: &Chain,
        saddle: Saddle,
    ) -> Result<Chain> {
        let mut context = Context::new(self);
        let mut output = ChainSum::default();
        for (cell, coefficient) in input.iter() {
            context.edges_for(model, cell, false, saddle, |model, target, sign| {
                output.add(
                    model,
                    target,
                    coefficient
                        .checked_mul(sign)
                        .ok_or_else(integral::arithmetic)?,
                )
            })?;
        }
        persist(model, output.terms.into_iter().collect())
    }

    pub(crate) fn euler(&self, model: &mut Model) -> Result<Vec<(i64, i128)>> {
        let state_count = 1usize
            .checked_shl(self.crossings.len() as u32)
            .ok_or_else(|| Error::plain("Khovanov cube exceeds the work limit"))?;
        model.spend(state_count)?;
        let mut euler = BTreeMap::<i64, i128>::new();
        for state in 0..state_count {
            let resolution = self.resolve(model, state as u128)?;
            let ones = state.count_ones() as usize;
            let circles = resolution.anchors.len();
            let j = self.shift(ones) + circles as i64;
            let sign = if ones.abs_diff(self.negative).is_multiple_of(2) {
                1i128
            } else {
                -1i128
            };
            let mut coefficient = sign;
            for k in 0..=circles {
                model.spend(1 + euler.len().max(1).ilog2() as usize)?;
                let e = j - 2 * k as i64;
                let slot = euler.entry(e).or_default();
                *slot = slot
                    .checked_add(coefficient)
                    .ok_or_else(integral::arithmetic)?;
                if *slot == 0 {
                    euler.remove(&e);
                }
                integral::storage(6 * euler.len())?;
                if k != circles {
                    let numerator = (circles - k) as u128;
                    let denominator = (k + 1) as u128;
                    let common = rv_lens::gcd(numerator, denominator);
                    coefficient = (coefficient / (denominator / common) as i128)
                        .checked_mul((numerator / common) as i128)
                        .ok_or_else(integral::arithmetic)?;
                }
            }
        }
        Ok(euler.into_iter().collect())
    }

    pub(crate) fn euler_at(&self, model: &mut Model, level: usize) -> Result<(usize, Value)> {
        let order = rv_tl::Grade::checked(level, self.ends.len(), 0)?.order;
        let mut sum = Rational::zero();
        for (degree, mut coefficient) in self.euler(model)? {
            if !degree.unsigned_abs().is_multiple_of(2) {
                coefficient = coefficient.checked_neg().ok_or_else(integral::arithmetic)?;
            }
            let term = rv_lens::shifted(
                model,
                order,
                &Rational::integer(coefficient),
                -2 * i128::from(degree),
            )?;
            sum = rv_lens::sum(model, &sum, &term)?;
        }
        Ok((order, rv_lens::rational_value(model, sum)?))
    }
}

struct Disjoint {
    parents: Vec<usize>,
    ranks: Vec<u8>,
}

impl Disjoint {
    fn new(size: usize) -> Self {
        Self {
            parents: (0..size).collect(),
            ranks: vec![0; size],
        }
    }

    fn root(&mut self, mut i: usize) -> usize {
        while i != self.parents[i] {
            self.parents[i] = self.parents[self.parents[i]];
            i = self.parents[i];
        }
        i
    }

    fn join(&mut self, a: usize, b: usize) {
        let (mut a, mut b) = (self.root(a), self.root(b));
        if a != b {
            if self.ranks[a] < self.ranks[b] {
                std::mem::swap(&mut a, &mut b);
            }
            self.parents[b] = a;
            self.ranks[a] += u8::from(self.ranks[a] == self.ranks[b]);
        }
    }
}

struct Resolution {
    labels: Vec<usize>,
    anchors: Vec<usize>,
}

impl Resolution {
    fn involved(&self, edge: Crossing) -> Vec<usize> {
        let mut output: Vec<_> = edge
            .before
            .into_iter()
            .chain(edge.after)
            .map(|i| self.labels[i])
            .collect();
        output.sort_unstable();
        output.dedup();
        output
    }

    fn cells(&self) -> usize {
        8 + self.labels.len() + self.anchors.len()
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Saddle {
    Khovanov,
    Perturbation,
    Lee,
}

pub(crate) struct Context<'a> {
    braid: &'a Braid,
    cache: HashMap<State, Rc<Resolution>>,
    cache_cells: usize,
}

impl<'a> Context<'a> {
    pub(crate) fn new(braid: &'a Braid) -> Self {
        Self {
            braid,
            cache: HashMap::new(),
            cache_cells: 0,
        }
    }

    fn resolution(&mut self, model: &mut Model, state: State) -> Result<Rc<Resolution>> {
        model.spend(1)?;
        if let Some(found) = self.cache.get(&state) {
            return Ok(found.clone());
        }
        let entry = Rc::new(self.braid.resolve(model, state)?);
        // The cache is optional. A sparse differential can stream additional
        // resolutions after its working cache fills.
        if integral::storage(self.cache_cells + entry.cells()).is_ok() {
            self.cache_cells += entry.cells();
            self.cache.insert(state, entry.clone());
        }
        Ok(entry)
    }

    fn edges(
        &mut self,
        model: &mut Model,
        cell: &Cell,
        reverse: bool,
        visit: impl FnMut(&mut Model, Cell, i128) -> Result<()>,
    ) -> Result<()> {
        self.edges_for(model, cell, reverse, Saddle::Khovanov, visit)
    }

    pub(crate) fn lee_edges(
        &mut self,
        model: &mut Model,
        cell: &Cell,
        reverse: bool,
        visit: impl FnMut(&mut Model, Cell, i128) -> Result<()>,
    ) -> Result<()> {
        self.edges_for(model, cell, reverse, Saddle::Lee, visit)
    }

    fn edges_for(
        &mut self,
        model: &mut Model,
        cell: &Cell,
        reverse: bool,
        saddle: Saddle,
        mut visit: impl FnMut(&mut Model, Cell, i128) -> Result<()>,
    ) -> Result<()> {
        let old = self.resolution(model, cell.state)?;
        for (i, &edge) in self.braid.crossings.iter().enumerate() {
            model.spend(1)?;
            let bit = 1u128 << i;
            if (cell.state & bit != 0) != reverse {
                continue;
            }
            let state = cell.state ^ bit;
            let new = self.resolution(model, state)?;
            let a = old.involved(edge);
            let b = new.involved(edge);
            let mut untouched = Vec::with_capacity(cell.xs.len());
            let mut involved_x = 0;
            model.spend(cell.xs.len())?;
            for &x in cell.xs.iter() {
                if a.contains(&x) {
                    involved_x += 1;
                } else {
                    untouched.push(new.labels[old.anchors[x]]);
                }
            }
            let sign = if (cell.state & (bit - 1)).count_ones().is_multiple_of(2) {
                1
            } else {
                -1
            };
            let mut output = |extra: &[usize]| {
                model.spend(1 + untouched.len() + extra.len())?;
                let mut xs = untouched.clone();
                xs.extend(extra);
                model.spend(xs.len() * (1 + xs.len().max(1).ilog2() as usize))?;
                xs.sort_unstable();
                visit(
                    model,
                    Cell {
                        state,
                        xs: xs.into(),
                    },
                    sign,
                )
            };
            if saddle != Saddle::Perturbation {
                // Reverse traversal uses the coordinate transpose of the
                // Frobenius maps. The cube sign is unchanged by toggling bit i.
                match (a.as_slice(), b.as_slice(), involved_x, reverse) {
                    ([_, _], [_], 0, false) | ([_, _], [_], 1, true) => output(&[])?,
                    ([_, _], [x], 1, false) | ([_, _], [x], 2, true) => output(&[*x])?,
                    ([_, _], [_], 2, false) | ([_, _], [_], 0, true) => {}
                    ([_], [x, y], 0, false) | ([_], [x, y], 1, true) => {
                        output(&[*x])?;
                        output(&[*y])?;
                    }
                    ([_], [x, y], 1, false) => output(&[*x, *y])?,
                    ([_], [_, _], 0, true) => output(&[])?,
                    _ => {
                        return Err(Error::plain(
                            "Khovanov saddle did not merge or split circles",
                        ))
                    }
                }
            }
            if saddle != Saddle::Khovanov {
                // Lee contributes m(X,X)=1 and the 1 tensor 1 term in
                // Delta(X). Traversing backward takes their transposes.
                match (a.as_slice(), b.as_slice(), involved_x, reverse) {
                    ([_, _], [_], 2, false) | ([_], [_, _], 1, false) => output(&[])?,
                    ([_, _], [x], 0, true) => output(&[*x])?,
                    ([_], [x, y], 0, true) => output(&[*x, *y])?,
                    _ => {}
                }
            }
        }
        Ok(())
    }

    fn basis(&mut self, model: &mut Model, grade: Grade) -> Result<Vec<Cell>> {
        let Some((ones, label_degree)) = self.braid.grading(grade) else {
            return Ok(Vec::new());
        };
        let mut cells = 0;
        let mut basis = Vec::new();
        self.visit_resolutions(model, ones, |model, state, circles| {
            let difference = circles as i64 - label_degree;
            if difference < 0 || difference as usize > 2 * circles || difference % 2 != 0 {
                return Ok(());
            }
            choose(model, circles, difference as usize / 2, |_, xs| {
                let cell = Cell {
                    state,
                    xs: xs.into(),
                };
                cells += cell.cells();
                integral::storage(cells)?;
                basis.push(cell);
                Ok(())
            })
        })?;
        model.spend(cells * (1 + basis.len().max(1).ilog2() as usize))?;
        basis.sort_unstable();
        Ok(basis)
    }

    pub(crate) fn lee_grade(&mut self, model: &mut Model, cell: &Cell) -> Result<Grade> {
        let circles = self.resolution(model, cell.state)?.anchors.len();
        Ok(self.braid.grade(cell, circles))
    }

    pub(crate) fn lee_basis(&mut self, model: &mut Model, i: i64) -> Result<Vec<(Cell, i64)>> {
        let ones = i128::from(i) + self.braid.negative as i128;
        if ones < 0 || ones > self.braid.crossings.len() as i128 {
            return Ok(Vec::new());
        }
        let ones = ones as usize;
        let shift = self.braid.shift(ones);
        let mut cells = 0;
        let mut basis = Vec::new();
        self.visit_resolutions(model, ones, |model, state, circles| {
            for count in 0..=circles {
                let j = shift + circles as i64 - 2 * count as i64;
                choose(model, circles, count, |_, xs| {
                    let cell = Cell {
                        state,
                        xs: xs.into(),
                    };
                    cells += 1 + cell.cells();
                    integral::storage(cells)?;
                    basis.push((cell, j));
                    Ok(())
                })?;
            }
            Ok(())
        })?;
        model.spend(cells * (1 + basis.len().max(1).ilog2() as usize))?;
        basis.sort_unstable_by(|(a, j), (b, k)| k.cmp(j).then_with(|| a.cmp(b)));
        Ok(basis)
    }

    fn visit_resolutions(
        &mut self,
        model: &mut Model,
        ones: usize,
        mut visit: impl FnMut(&mut Model, State, usize) -> Result<()>,
    ) -> Result<()> {
        choose(model, self.braid.crossings.len(), ones, |model, ones| {
            let state = ones.iter().fold(0u128, |acc, i| acc | (1u128 << i));
            let resolution = self.resolution(model, state)?;
            visit(model, state, resolution.anchors.len())
        })
    }
}

fn choose(
    model: &mut Model,
    n: usize,
    k: usize,
    mut visit: impl FnMut(&mut Model, &[usize]) -> Result<()>,
) -> Result<()> {
    model.spend(1 + k)?;
    let mut indices: Vec<_> = (0..k).collect();
    loop {
        model.spend(1 + k)?;
        visit(model, &indices)?;
        let Some(next) = (0..k).rev().find(|&i| indices[i] < n - k + i) else {
            return Ok(());
        };
        indices[next] += 1;
        let start = indices[next];
        for (offset, index) in indices[next + 1..].iter_mut().enumerate() {
            *index = start + offset + 1;
        }
    }
}

#[derive(Default)]
struct ChainSum {
    terms: BTreeMap<Cell, i128>,
    cells: usize,
}

impl ChainSum {
    fn add(&mut self, model: &mut Model, cell: Cell, coefficient: i128) -> Result<()> {
        if coefficient == 0 {
            return Ok(());
        }
        model.spend(cell.cells() * (1 + self.terms.len().max(1).ilog2() as usize))?;
        use std::collections::btree_map::Entry;
        match self.terms.entry(cell) {
            Entry::Vacant(entry) => {
                self.cells += entry.key().cells();
                integral::storage(self.cells)?;
                entry.insert(coefficient);
            }
            Entry::Occupied(mut entry) => {
                let total = entry
                    .get()
                    .checked_add(coefficient)
                    .ok_or_else(integral::arithmetic)?;
                if total == 0 {
                    self.cells -= entry.key().cells();
                    entry.remove();
                } else {
                    *entry.get_mut() = total;
                }
            }
        }
        Ok(())
    }
}

fn persist(model: &mut Model, terms: Vec<(Cell, i128)>) -> Result<Chain> {
    let cells = terms.iter().map(|(cell, _)| cell.cells()).sum::<usize>();
    integral::storage(cells)?;
    model.allocate(cells)?;
    Ok(terms.into())
}

pub(crate) fn add_chains(model: &mut Model, a: &Chain, b: &Chain) -> Result<Chain> {
    let mut sum = ChainSum::default();
    for (cell, coefficient) in a.iter().chain(b.iter()) {
        sum.add(model, cell.clone(), *coefficient)?;
    }
    persist(model, sum.terms.into_iter().collect())
}

pub(crate) fn scale_chain(model: &mut Model, input: &Chain, factor: i128) -> Result<Chain> {
    if factor == 0 {
        return Ok(Rc::default());
    }
    model.spend(input.len())?;
    let output = input
        .iter()
        .map(|(cell, value)| {
            let value = value.checked_mul(factor).ok_or_else(integral::arithmetic)?;
            Ok((cell.clone(), value))
        })
        .collect::<Result<_>>()?;
    persist(model, output)
}

#[derive(Clone)]
pub(crate) struct Invariants {
    pub(crate) chains: usize,
    pub(crate) free: usize,
    pub(crate) torsion: Vec<i128>,
}

pub(crate) struct Group {
    basis: Vec<Cell>,
    previous: Vec<Cell>,
    outgoing: integral::Smith,
    incoming: integral::Smith,
    indices: Vec<usize>,
    invariants: Invariants,
}

impl Group {
    fn build(model: &mut Model, braid: &Braid, grade: Grade) -> Result<Self> {
        let mut context = Context::new(braid);
        let basis = context.basis(model, grade)?;
        let after = incident(model, &mut context, &basis, false)?;
        let outgoing = incidence_smith(model, basis.len(), after, Track::Both)?;
        let rank = outgoing.factors.len();
        let kernel = basis.len() - rank;

        let before = incident(model, &mut context, &basis, true)?;
        let mut previous = Vec::with_capacity(before.len());
        let incoming = if before.is_empty() {
            integral::Smith::empty()
        } else {
            let mut incoming = Matrix::new(model, kernel, before.len())?;
            for (column, (source, incidences)) in before.into_iter().enumerate() {
                previous.push(source);
                // The domain basis of the outgoing Smith split has its kernel in
                // the trailing coordinates. Express the actual incoming column
                // in these coordinates before doing Smith in that kernel.
                if rank == 0 {
                    for (row, value) in incidences {
                        incoming.write(model, row, column, value)?;
                    }
                } else {
                    for row in 0..basis.len() {
                        let value = dot(model, outgoing.columns.coordinate(row), &incidences)?;
                        if row >= rank {
                            if value != 0 {
                                incoming.write(model, row - rank, column, value)?;
                            }
                        } else if value != 0 {
                            return Err(Error::plain(
                                "Khovanov differential does not square to zero",
                            ));
                        }
                    }
                }
            }
            integral::reduce(model, incoming, Track::Both, Track::Basis)?
        };
        let free = kernel - incoming.factors.len();
        let mut torsion = Vec::new();
        let mut indices: Vec<_> = (incoming.factors.len()..kernel).collect();
        for (i, &factor) in incoming.factors.iter().enumerate() {
            if factor > 1 {
                torsion.push(factor);
                indices.push(i);
            }
        }
        let cells = basis
            .iter()
            .chain(previous.iter())
            .map(Cell::cells)
            .sum::<usize>()
            + incoming.rows.cells()
            + incoming.columns.cells()
            + outgoing.columns.cells()
            + 2 * (outgoing.factors.len() + incoming.factors.len() + torsion.len())
            + indices.len()
            + 1;
        integral::storage(cells)?;
        model.allocate(cells)?;
        let invariants = Invariants {
            chains: basis.len(),
            free,
            torsion,
        };
        Ok(Self {
            basis,
            previous,
            outgoing,
            incoming,
            indices,
            invariants,
        })
    }

    pub(crate) fn elements(&self, model: &mut Model, limit: usize) -> Result<Vec<Value>> {
        if self.invariants.free != 0 {
            return Err(Error::plain(
                "cannot enumerate a Khovanov group with a free summand",
            ));
        }
        let mut size = 1usize;
        for &factor in &self.invariants.torsion {
            let rad = usize::try_from(factor).unwrap_or(usize::MAX);
            size = size.saturating_mul(rad);
            if size > limit {
                return Err(Error::plain(
                    "Khovanov torsion group exceeds the finite element limit",
                ));
            }
        }
        let rank = self.invariants.torsion.len();
        model.allocate(size * (1 + 3 * rank))?;
        model.spend(size * (1 + rank) * (1 + size.ilog2() as usize))?;
        let mut values: Vec<Class> = Vec::with_capacity(size);
        for mut code in 0..size {
            let mut coordinates = Vec::new();
            for (i, &factor) in self.invariants.torsion.iter().enumerate() {
                let radix = factor as usize;
                let digit = code % radix;
                if digit != 0 {
                    coordinates.push((i, digit as i128));
                }
                code /= radix;
            }
            values.push(coordinates.into());
        }
        values.sort_unstable();
        Ok(values.into_iter().map(Value::KhClass).collect())
    }

    pub(crate) fn order(&self, generator: usize) -> Result<Option<i128>> {
        let coordinate = self.index(generator)?;
        Ok(self.incoming.factors.get(coordinate).copied())
    }

    fn index(&self, generator: usize) -> Result<usize> {
        self.indices
            .get(generator)
            .copied()
            .ok_or_else(|| Error::plain("Khovanov generator index is out of range"))
    }

    pub(crate) fn generator(&self, model: &mut Model, generator: usize) -> Result<Class> {
        self.index(generator)?;
        persist_class(model, [(generator, 1)].into())
    }

    pub(crate) fn project(&self, model: &mut Model, chain: &Chain) -> Result<Class> {
        if chain.is_empty() {
            return Ok(Rc::default());
        }
        let mut raw = Vector::new();
        for (cell, value) in chain.iter() {
            model.spend(cell.cells() * (1 + self.basis.len().max(1).ilog2() as usize))?;
            let index = self
                .basis
                .binary_search(cell)
                .map_err(|_| Error::plain("Khovanov chain term has the wrong bidegree"))?;
            raw.insert(index, *value);
        }
        let rank = self.outgoing.factors.len();
        let kernel = if rank == 0 {
            raw
        } else {
            let mut kernel = Vector::new();
            for row in 0..self.basis.len() {
                let value = dot(model, self.outgoing.columns.coordinate(row), &raw)?;
                if value == 0 {
                    continue;
                }
                if row < rank {
                    return Err(Error::plain("Khovanov chain is not a cycle"));
                }
                kernel.insert(row - rank, value);
            }
            kernel
        };
        let mut result = Vector::new();
        let incoming_rank = self.incoming.factors.len();
        for (generator, &coordinate) in self.indices.iter().enumerate() {
            model.spend(1)?;
            let mut value = if incoming_rank == 0 {
                kernel.get(&coordinate).copied().unwrap_or(0)
            } else {
                dot(model, self.incoming.rows.coordinate(coordinate), &kernel)?
            };
            if let Some(&order) = self.incoming.factors.get(coordinate) {
                value = value.rem_euclid(order);
            }
            if value != 0 {
                result.insert(generator, value);
            }
        }
        persist_class(model, result)
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &Class) -> Result<Chain> {
        model.spend(class.len() * (1 + class.len().max(1).ilog2() as usize))?;
        let mut kernel = Vector::new();
        for &(generator, factor) in class.iter() {
            let i = self.index(generator)?;
            if self.incoming.factors.is_empty() {
                kernel.insert(i, factor);
            } else {
                integral::add_scaled(model, &mut kernel, self.incoming.rows.basis(i), factor)?;
            }
        }
        let mut raw = Vector::new();
        let rank = self.outgoing.factors.len();
        if rank == 0 {
            raw = kernel;
        } else {
            for (i, factor) in kernel {
                integral::add_scaled(
                    model,
                    &mut raw,
                    self.outgoing.columns.basis(i + rank),
                    factor,
                )?;
            }
        }
        raw_chain(model, &self.basis, &raw)
    }

    pub(crate) fn null(&self, model: &mut Model, generator: usize) -> Result<Chain> {
        let i = self.index(generator)?;
        if self.incoming.factors.get(i).is_none() {
            return Err(Error::plain(
                "a free Khovanov generator has no torsion null",
            ));
        }
        raw_chain(model, &self.previous, self.incoming.columns.basis(i))
    }

    pub(crate) fn add(&self, model: &mut Model, a: &Class, b: &Class) -> Result<Class> {
        let mut result = Vector::new();
        model.spend((a.len() + b.len()) * (1 + (a.len() + b.len()).max(1).ilog2() as usize))?;
        result.extend(a.iter().copied());
        for &(i, coefficient) in b.iter() {
            let first = result.get(&i).copied().unwrap_or(0);
            let value = if let Some(order) = self.order(i)? {
                ((first as u128 + coefficient as u128) % order as u128) as i128
            } else {
                first
                    .checked_add(coefficient)
                    .ok_or_else(integral::arithmetic)?
            };
            if value == 0 {
                result.remove(&i);
            } else {
                result.insert(i, value);
            }
        }
        persist_class(model, result)
    }

    pub(crate) fn scale(&self, model: &mut Model, a: &Class, factor: i128) -> Result<Class> {
        if factor == 0 {
            return Ok(Rc::default());
        }
        model.spend(a.len())?;
        let mut result = Vector::new();
        for &(i, coefficient) in a.iter() {
            let value = if let Some(order) = self.order(i)? {
                modular_product(model, coefficient, factor, order)?
            } else {
                coefficient
                    .checked_mul(factor)
                    .ok_or_else(integral::arithmetic)?
            };
            if value != 0 {
                result.insert(i, value);
            }
        }
        persist_class(model, result)
    }
}

fn dot(model: &mut Model, left: &Vector, right: &Vector) -> Result<i128> {
    let (small, large) = if left.len() <= right.len() {
        (left, right)
    } else {
        (right, left)
    };
    model.spend(1 + small.len() * (1 + large.len().max(1).ilog2() as usize))?;
    let mut sum = 0i128;
    for (index, value) in small {
        if let Some(&other) = large.get(index) {
            let product = value.checked_mul(other).ok_or_else(integral::arithmetic)?;
            sum = sum.checked_add(product).ok_or_else(integral::arithmetic)?;
        }
    }
    Ok(sum)
}

fn modular_product(model: &mut Model, a: i128, b: i128, modulus: i128) -> Result<i128> {
    let (mut a, mut b) = (a as u128, b.rem_euclid(modulus) as u128);
    let modulus = modulus as u128;
    if let Some(product) = a.checked_mul(b) {
        return Ok((product % modulus) as i128);
    }
    model.spend(u128::BITS as usize)?;
    let mut result = 0;
    while b != 0 {
        if b % 2 != 0 {
            result = (result + a) % modulus;
        }
        a = (2 * a) % modulus;
        b /= 2;
    }
    Ok(result as i128)
}

fn incidence_smith(
    model: &mut Model,
    columns: usize,
    adjacency: BTreeMap<Cell, Vector>,
    track: Track,
) -> Result<integral::Smith> {
    if adjacency.is_empty() {
        return Ok(integral::Smith::empty());
    }
    let mut matrix = Matrix::new(model, adjacency.len(), columns)?;
    for (row, entries) in adjacency.into_values().enumerate() {
        for (column, value) in entries {
            matrix.write(model, row, column, value)?;
        }
    }
    integral::reduce(model, matrix, Track::None, track)
}

fn incident(
    model: &mut Model,
    context: &mut Context<'_>,
    basis: &[Cell],
    reverse: bool,
) -> Result<BTreeMap<Cell, Vector>> {
    let mut cells = 0;
    let mut incidences = BTreeMap::<Cell, Vector>::new();
    for (index, cell) in basis.iter().enumerate() {
        context.edges(model, cell, reverse, |model, opposite, sign| {
            model.spend(opposite.cells() * (1 + incidences.len().max(1).ilog2() as usize))?;
            let width = opposite.cells();
            let column = incidences.entry(opposite).or_insert_with(|| {
                cells += width + 3;
                Vector::new()
            });
            cells += 6;
            integral::storage(cells)?;
            if column.insert(index, sign).is_some() {
                return Err(Error::plain(
                    "Khovanov saddle produced a duplicate incidence",
                ));
            }
            Ok(())
        })?;
    }
    Ok(incidences)
}

fn raw_chain(model: &mut Model, basis: &[Cell], vector: &Vector) -> Result<Chain> {
    model.spend(vector.len())?;
    let terms = vector
        .iter()
        .map(|(&i, &value)| (basis[i].clone(), value))
        .collect();
    persist(model, terms)
}

fn persist_class(model: &mut Model, vector: Vector) -> Result<Class> {
    integral::storage(3 * vector.len())?;
    model.allocate(3 * vector.len())?;
    Ok(vector.into_iter().collect())
}

pub(crate) fn compare_chains(model: &mut Model, a: &Chain, b: &Chain) -> Result<Ordering> {
    if Rc::ptr_eq(a, b) {
        return Ok(Ordering::Equal);
    }
    for ((g, c), (h, d)) in a.iter().zip(b.iter()) {
        model.spend(g.cells().min(h.cells()))?;
        let ordering = g.cmp(h).then_with(|| c.cmp(d));
        if !ordering.is_eq() {
            return Ok(ordering);
        }
    }
    Ok(a.len().cmp(&b.len()))
}

pub(crate) fn compare_classes(model: &mut Model, a: &Class, b: &Class) -> Result<Ordering> {
    if Rc::ptr_eq(a, b) {
        return Ok(Ordering::Equal);
    }
    model.spend(a.len().min(b.len()))?;
    Ok(a.cmp(b))
}
