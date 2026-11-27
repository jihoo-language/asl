use std::cmp::Ordering;
use std::collections::{btree_map::Entry, BTreeMap, BTreeSet};
use std::fmt;

use super::{persist, require_cycle, Chain, Terms};
use crate::error::{Error, Result};
use crate::integer_smith::storage;
use crate::khovanov::{Braid, Cell, Context};
use crate::lee_linear::{arithmetic, Echelon, Rational, Vector};
use crate::model::Model;
use crate::rv_lens::gcd;
use crate::syntax::{Label, ParsedExpr};

#[derive(Clone, Copy, Eq, PartialEq)]
pub(crate) struct Parameter {
    p: usize,
    q: usize,
}

impl Parameter {
    const ZERO: Self = Self { p: 0, q: 1 };
    const TWO: Self = Self { p: 2, q: 1 };

    pub(crate) fn new(p: usize, q: usize) -> Result<Self> {
        if q == 0 || p as u128 > 2 * q as u128 {
            return Err(Error::plain(
                "annular Rasmussen parameter requires 0 <= p/q <= 2 with q positive",
            ));
        }
        let divisor = gcd(p as u128, q as u128) as usize;
        Ok(Self {
            p: p / divisor,
            q: q / divisor,
        })
    }

    fn between(self, next: Self) -> Result<Self> {
        Self::new(
            self.p.checked_add(next.p).ok_or_else(arithmetic)?,
            self.q.checked_add(next.q).ok_or_else(arithmetic)?,
        )
    }

    pub(crate) fn pair(self) -> (usize, usize) {
        (self.p, self.q)
    }
}

impl Ord for Parameter {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.p as u128 * other.q as u128).cmp(&(other.p as u128 * self.q as u128))
    }
}

impl PartialOrd for Parameter {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for Parameter {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.q == 1 {
            write!(f, "{}", self.p)
        } else {
            write!(f, "{}/{}", self.p, self.q)
        }
    }
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Degree {
    pub(crate) j: i64,
    pub(crate) k: i64,
}

impl Degree {
    fn numerator(self, t: Parameter) -> i128 {
        i128::from(self.j) * t.q as i128 - i128::from(self.k) * t.p as i128
    }

    fn at(self, t: Parameter) -> Result<Rational> {
        Rational::new(self.numerator(t), t.q as i128)
    }

    fn intersection(self, other: Self) -> Result<Option<Parameter>> {
        let p = i128::from(self.j) - i128::from(other.j);
        let q = i128::from(self.k) - i128::from(other.k);
        if q == 0 || p.signum() != q.signum() || p.unsigned_abs() >= 2 * q.unsigned_abs() {
            return Ok(None);
        }
        let p = usize::try_from(p.unsigned_abs()).map_err(|_| arithmetic())?;
        let q = usize::try_from(q.unsigned_abs()).map_err(|_| arithmetic())?;
        Parameter::new(p, q).map(Some)
    }
}

pub(crate) type Piece = (Parameter, Parameter, Degree);

pub(crate) enum Command {
    Value(Label, Parameter),
    Profile(Label),
    Weight(Parameter, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Witness {
    Representative,
    Primitive,
    Dual,
}

pub(crate) enum Query {
    Orientation(Label),
    Witness(Label, Parameter, Witness),
    DualD(Box<ParsedExpr>),
    Pair(Box<[ParsedExpr; 2]>),
}

fn orientation_terms(model: &mut Model, braid: &Braid) -> Result<Chain> {
    let (_, seed) = braid.psi_cell(model)?;
    let strands = seed.xs.len();
    let count = 1usize
        .checked_shl(strands as u32)
        .ok_or_else(|| Error::plain("Lee orientation cycle exceeds the computation limit"))?;
    storage(
        count
            .saturating_mul(9)
            .saturating_add((count / 2).saturating_mul(strands)),
    )?;
    model.spend(count * (1 + strands))?;
    let mut terms = Vec::with_capacity(count);
    for label in 0..count {
        let xs: Vec<_> = (0..strands).filter(|&i| label & (1 << i) != 0).collect();
        // The oriented circles are nested in order of the input strands.
        // Expand (1+X) on the even positions and (1-X) on the odd ones.
        let parity = xs.iter().filter(|&&i| i % 2 == 1).count();
        let sign = if parity.is_multiple_of(2) { 1 } else { -1 };
        terms.push((
            Cell {
                state: seed.state,
                xs: xs.into(),
            },
            Rational::integer(sign),
        ));
    }
    model.spend(count * (1 + strands) * (1 + count.ilog2() as usize))?;
    terms.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let terms: Chain = terms.into();
    require_cycle(model, &mut Context::new(braid), &terms)?;
    Ok(terms)
}

pub(crate) fn orientation(model: &mut Model, braid: &Braid) -> Result<Chain> {
    let terms = orientation_terms(model, braid)?;
    persist(model, terms.to_vec())
}

#[derive(Default)]
struct Rows {
    ids: BTreeMap<Cell, usize>,
    cells: Vec<(Cell, Degree)>,
    width: usize,
}

impl Rows {
    fn intern(
        &mut self,
        model: &mut Model,
        context: &mut Context<'_>,
        cell: Cell,
    ) -> Result<usize> {
        model.spend(cell.cells() * (1 + self.ids.len().max(1).ilog2() as usize))?;
        match self.ids.entry(cell) {
            Entry::Occupied(found) => Ok(*found.get()),
            Entry::Vacant(entry) => {
                let cell = entry.key();
                self.width += 2 * cell.cells() + 7;
                storage(self.width)?;
                let j = context.lee_grade(model, cell)?.j;
                let k = context.annular_k(model, cell)?;
                let position = self.cells.len();
                self.cells.push((cell.clone(), Degree { j, k }));
                entry.insert(position);
                Ok(position)
            }
        }
    }
}

struct Problem {
    target: Vec<(Cell, Degree)>,
    sources: Vec<(Cell, i64)>,
    incoming: Vec<Vector>,
    oriented: Vector,
}

impl Problem {
    fn build(model: &mut Model, braid: &Braid) -> Result<Self> {
        let orientation = orientation_terms(model, braid)?;
        let mut context = Context::with_annular(braid, true);
        let mut rows = Rows::default();
        let mut oriented = Vec::with_capacity(orientation.len());
        for (cell, coefficient) in orientation.iter() {
            let i = rows.intern(model, &mut context, cell.clone())?;
            oriented.push((i, coefficient.numerator));
        }

        let sources = context.lee_basis(model, -1)?;
        let incoming_base = sources
            .iter()
            .map(|(cell, _)| cell.cells() + 6)
            .sum::<usize>();
        let mut incoming_width = incoming_base + 4 * oriented.len();
        storage(incoming_width + rows.width)?;
        let mut incoming = Vec::with_capacity(sources.len());
        for (cell, j) in &sources {
            let degree = Degree {
                j: *j,
                k: context.annular_k(model, cell)?,
            };
            let mut terms = Vec::new();
            context.lee_edges(model, cell, false, |model, adjacent, sign| {
                model.spend(1)?;
                terms.push((adjacent, sign));
                Ok(())
            })?;
            let mut column = Vec::with_capacity(terms.len());
            for (adjacent, sign) in terms {
                let i = rows.intern(model, &mut context, adjacent)?;
                let target = rows.cells[i].1;
                if target.j < degree.j
                    || target.numerator(Parameter::TWO) < degree.numerator(Parameter::TWO)
                {
                    return Err(Error::plain(
                        "Lee differential violates the annular Rasmussen filtration",
                    ));
                }
                column.push((i, sign));
            }
            incoming_width += 4 * column.len();
            storage(incoming_width + rows.width)?;
            incoming.push(column);
        }
        Ok(Self {
            target: rows.cells,
            sources,
            incoming,
            oriented,
        })
    }

    fn reduction(&self, model: &mut Model, t: Parameter, track: bool) -> Result<Reduction> {
        let mut positions: Vec<_> = (0..self.target.len()).collect();
        let weights: Vec<_> = self
            .target
            .iter()
            .map(|(_, degree)| degree.numerator(t))
            .collect();
        let cost = 1 + positions.len().max(1).ilog2() as usize;
        model.spend(positions.len() * cost)?;
        // The Lee Echelon uses the largest coordinate as pivot. Place
        // least-filtered terms last, where elimination tests them first.
        positions.sort_unstable_by(|&i, &j| {
            weights[j]
                .cmp(&weights[i])
                .then_with(|| self.target[i].0.cmp(&self.target[j].0))
        });
        let mut index = vec![0; positions.len()];
        for (i, &j) in positions.iter().enumerate() {
            index[j] = i;
        }
        let remap = |model: &mut Model, column: &Vector| -> Result<Vector> {
            model.spend(column.len() * (1 + column.len().max(1).ilog2() as usize))?;
            let mut column: Vector = column.iter().map(|&(i, sign)| (index[i], sign)).collect();
            column.sort_unstable_by_key(|&(i, _)| i);
            debug_assert!(column.windows(2).all(|pair| pair[0].0 != pair[1].0));
            Ok(column)
        };

        let mut boundary = Echelon::default();
        for (source, column) in self.incoming.iter().enumerate() {
            if !column.is_empty() {
                let column = remap(model, column)?;
                let witness = if track {
                    vec![(source + 1, 1)]
                } else {
                    Vec::new()
                };
                boundary.insert(model, column, witness)?;
            }
        }
        let oriented = remap(model, &self.oriented)?;
        let seed = if track { vec![(0, 1)] } else { Vec::new() };
        let (image, witness) = boundary.reduced(model, oriented, seed)?;
        let &(last, coefficient) = image
            .last()
            .ok_or_else(|| Error::plain("canonical Lee orientation cycle is a boundary"))?;
        Ok(Reduction {
            boundary,
            positions,
            image,
            witness,
            last,
            coefficient,
        })
    }

    fn best_degree(&self, model: &mut Model, t: Parameter) -> Result<Degree> {
        let solved = self.reduction(model, t, false)?;
        Ok(self.target[solved.positions[solved.last]].1)
    }
}

struct Reduction {
    boundary: Echelon,
    positions: Vec<usize>,
    image: Vector,
    witness: Vector,
    last: usize,
    coefficient: i128,
}

impl Reduction {
    fn orientation_scale(&self) -> i128 {
        let (i, coefficient) = self.witness[0];
        debug_assert_eq!(i, 0);
        coefficient
    }

    fn dual(&self, model: &mut Model, problem: &Problem) -> Result<Chain> {
        let mut values = vec![Rational::ZERO; self.positions.len()];
        // This nonpivot coordinate detects the reduced orientation.
        // All earlier (hence all strictly higher weight) entries vanish.
        values[self.last] = Rational::ratio(self.orientation_scale(), self.coefficient)?;
        let mut width = 0;
        let mut terms = Vec::new();
        for (i, &position) in self.positions.iter().enumerate().skip(self.last) {
            model.spend(1)?;
            if let Some(row) = self.boundary.pivot(i) {
                let (_, pivot) = row[row.len() - 1];
                model.spend(row.len())?;
                let mut dot = Rational::ZERO;
                for &(k, coefficient) in &row[..row.len() - 1] {
                    dot = dot.add(values[k].multiply(Rational::integer(coefficient))?)?;
                }
                values[i] = dot.multiply(Rational::ratio(-1, pivot)?)?;
            }
            if values[i] != Rational::ZERO {
                let cell = &problem.target[position].0;
                width += cell.cells() + 2;
                storage(width)?;
                terms.push((cell.clone(), values[i]));
            }
        }
        sorted(model, terms)
    }
}

fn sorted(model: &mut Model, mut terms: Vec<(Cell, Rational)>) -> Result<Chain> {
    let cells = terms.iter().map(|(cell, _)| cell.cells()).sum::<usize>();
    model.spend(cells * (1 + terms.len().max(1).ilog2() as usize))?;
    terms.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    persist(model, terms)
}

pub(crate) fn witness(
    model: &mut Model,
    braid: &Braid,
    t: Parameter,
    kind: Witness,
) -> Result<Chain> {
    let problem = Problem::build(model, braid)?;
    let solved = problem.reduction(model, t, true)?;
    let denominator = solved.orientation_scale();
    let terms = match kind {
        Witness::Representative => solved
            .image
            .iter()
            .map(|&(i, coefficient)| {
                let cell = problem.target[solved.positions[i]].0.clone();
                Ok((cell, Rational::ratio(coefficient, denominator)?))
            })
            .collect::<Result<_>>()?,
        Witness::Primitive => solved
            .witness
            .iter()
            .skip(1)
            .map(|&(i, coefficient)| {
                Ok((
                    problem.sources[i - 1].0.clone(),
                    Rational::ratio(coefficient, denominator)?,
                ))
            })
            .collect::<Result<_>>()?,
        Witness::Dual => return solved.dual(model, &problem),
    };
    sorted(model, terms)
}

pub(crate) fn value(model: &mut Model, braid: &Braid, t: Parameter) -> Result<Rational> {
    let problem = Problem::build(model, braid)?;
    problem.best_degree(model, t)?.at(t)
}

pub(crate) fn profile(model: &mut Model, braid: &Braid) -> Result<Vec<Piece>> {
    let problem = Problem::build(model, braid)?;
    let mut degrees = BTreeSet::new();
    for (_, degree) in &problem.target {
        model.spend(1 + degrees.len().max(1).ilog2() as usize)?;
        degrees.insert(*degree);
    }
    let forms: Vec<_> = degrees.into_iter().collect();
    let mut cuts = BTreeSet::from([Parameter::ZERO, Parameter::TWO]);
    for (i, &left) in forms.iter().enumerate() {
        model.spend(forms.len() - i - 1)?;
        for &right in &forms[i + 1..] {
            if let Some(t) = left.intersection(right)? {
                model.spend(1 + cuts.len().ilog2() as usize)?;
                cuts.insert(t);
                storage(4 * cuts.len())?;
            }
        }
    }
    let cuts: Vec<_> = cuts.into_iter().collect();
    let mut pieces: Vec<Piece> = Vec::new();
    for endpoints in cuts.windows(2) {
        let (left, right) = (endpoints[0], endpoints[1]);
        // The comparison of every pair of weights is constant between
        // successive crossings of their lines, so the pivot is fixed.
        let active = problem.best_degree(model, left.between(right)?)?;
        if let Some((_, end, previous)) = pieces.last_mut() {
            if *previous == active {
                *end = right;
                continue;
            }
            debug_assert_eq!(previous.numerator(left), active.numerator(left));
        }
        pieces.push((left, right, active));
    }
    Ok(pieces)
}

pub(crate) fn dual_d(model: &mut Model, braid: &Braid, support: &Chain) -> Result<Chain> {
    let mut context = Context::new(braid);
    let mut terms = Terms::default();
    for (cell, coefficient) in support.iter() {
        context.lee_edges(model, cell, true, |model, target, sign| {
            terms.insert(
                model,
                target,
                coefficient.multiply(Rational::integer(sign))?,
            )
        })?;
    }
    persist(model, terms.values.into_iter().collect())
}

pub(crate) fn pair(model: &mut Model, dual: &Chain, chain: &Chain) -> Result<Rational> {
    let (mut a, mut b) = (dual.iter().peekable(), chain.iter().peekable());
    let mut result = Rational::ZERO;
    while let (Some((x, c)), Some((y, d))) = (a.peek(), b.peek()) {
        model.spend(x.cells().min(y.cells()))?;
        let order = x.cmp(y);
        if order.is_eq() {
            result = result.add(c.multiply(*d)?)?;
        }
        if !order.is_gt() {
            a.next();
        }
        if !order.is_lt() {
            b.next();
        }
    }
    Ok(result)
}

pub(crate) fn weight(
    model: &mut Model,
    braid: &Braid,
    support: &Chain,
    t: Parameter,
    dual: bool,
) -> Result<Option<Rational>> {
    let mut extrema = None;
    let mut context = Context::with_annular(braid, true);
    for (cell, _) in support.iter() {
        let j = context.lee_grade(model, cell)?.j;
        let k = context.annular_k(model, cell)?;
        let value = Degree { j, k }.numerator(t);
        extrema = Some(extrema.map_or(value, |current: i128| {
            if dual {
                current.max(value)
            } else {
                current.min(value)
            }
        }));
    }
    extrema
        .map(|weight| Rational::new(weight, t.q as i128))
        .transpose()
}
