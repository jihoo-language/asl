use std::collections::{btree_map::Entry, BTreeMap};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::integer_smith::storage;
use crate::khovanov::{Braid, Cell, Chain as IntegralChain, Context};
use crate::lee_linear::{checked_multiply, Echelon, Rational, Vector};
use crate::model::Model;
use crate::rv_lens::gcd;

pub(crate) mod dt;

pub(crate) type Chain = Rc<[(Cell, Rational)]>;
type Basis = Vec<(Cell, i64)>;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) i: i64,
    pub(crate) f: i64,
}

impl Grade {
    pub(crate) fn shift(self, direction: i64) -> Result<Self> {
        let i = self
            .i
            .checked_add(direction)
            .ok_or_else(|| Error::plain("Lee homological degree is out of range"))?;
        Ok(Self { i, ..self })
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Lift,
    Perturbation,
    Differential,
    Filter(i64),
    Scale(Rational),
    Null,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Lift => "leeLift",
            Self::Perturbation => "leePhi",
            Self::Differential => "leeD",
            Self::Filter(_) => "leeFilter",
            Self::Scale(_) => "leeScale",
            Self::Null => "leeNull",
        }
    }
}

pub(crate) fn lift(model: &mut Model, chain: &IntegralChain) -> Result<Chain> {
    model.spend(chain.len())?;
    persist(
        model,
        chain
            .iter()
            .map(|(cell, c)| (cell.clone(), Rational::integer(*c)))
            .collect(),
    )
}

pub(crate) fn scale(model: &mut Model, chain: &Chain, scalar: Rational) -> Result<Chain> {
    if scalar == Rational::ZERO {
        return Ok(Rc::default());
    }
    model.spend(chain.len())?;
    let terms = chain
        .iter()
        .map(|(cell, c)| Ok((cell.clone(), c.multiply(scalar)?)))
        .collect::<Result<_>>()?;
    persist(model, terms)
}

pub(crate) fn add(model: &mut Model, left: &Chain, right: &Chain) -> Result<Chain> {
    let mut terms = Terms::default();
    for (cell, c) in left.iter().chain(right.iter()) {
        terms.insert(model, cell.clone(), *c)?;
    }
    persist(model, terms.values.into_iter().collect())
}

pub(crate) fn filter(model: &mut Model, braid: &Braid, chain: &Chain, f: i64) -> Result<Chain> {
    let mut context = Context::new(braid);
    for (cell, _) in chain.iter() {
        let j = context.lee_grade(model, cell)?.j;
        if j < f {
            return Err(Error::plain(format!(
                "Lee cochain has a quantum {j} term below filtration {f}"
            )));
        }
    }
    Ok(chain.clone())
}

pub(crate) fn differential(model: &mut Model, braid: &Braid, chain: &Chain) -> Result<Chain> {
    let mut context = Context::new(braid);
    let terms = differential_terms(model, &mut context, chain)?;
    persist(model, terms.values.into_iter().collect())
}

fn differential_terms(
    model: &mut Model,
    context: &mut Context<'_>,
    chain: &Chain,
) -> Result<Terms> {
    let mut terms = Terms::default();
    for (cell, c) in chain.iter() {
        context.lee_edges(model, cell, false, |model, target, sign| {
            terms.insert(model, target, c.multiply(Rational::integer(sign))?)
        })?;
    }
    Ok(terms)
}

fn require_cycle(model: &mut Model, context: &mut Context<'_>, chain: &Chain) -> Result<()> {
    if differential_terms(model, context, chain)?.values.is_empty() {
        Ok(())
    } else {
        Err(Error::plain("Lee cochain is not a cycle"))
    }
}

#[derive(Default)]
struct Terms {
    values: BTreeMap<Cell, Rational>,
    cells: usize,
}

impl Terms {
    fn insert(&mut self, model: &mut Model, cell: Cell, coefficient: Rational) -> Result<()> {
        let cells = cell.cells() + 2;
        model.spend(cells * (1 + self.values.len().max(1).ilog2() as usize))?;
        match self.values.entry(cell) {
            Entry::Vacant(entry) => {
                self.cells += cells;
                storage(self.cells)?;
                entry.insert(coefficient);
            }
            Entry::Occupied(mut entry) => {
                let sum = entry.get().add(coefficient)?;
                if sum == Rational::ZERO {
                    self.cells -= cells;
                    entry.remove();
                } else {
                    *entry.get_mut() = sum;
                }
            }
        }
        Ok(())
    }
}

fn persist(model: &mut Model, value: Vec<(Cell, Rational)>) -> Result<Chain> {
    let cells = value.iter().map(|(cell, _)| cell.cells() + 2).sum();
    storage(cells)?;
    model.allocate(cells)?;
    Ok(value.into())
}

fn to_integer(
    model: &mut Model,
    coordinate: &BTreeMap<&Cell, usize>,
    chain: &Chain,
) -> Result<(Vector, i128)> {
    let mut denominator = 1i128;
    for (_, c) in chain.iter() {
        model.spend(1)?;
        let common = gcd(denominator as u128, c.denominator as u128) as i128;
        denominator = checked_multiply(denominator / common, c.denominator)?;
    }
    let mut terms = Vec::with_capacity(chain.len());
    for (cell, c) in chain.iter() {
        model.spend(cell.cells() * (1 + coordinate.len().max(1).ilog2() as usize))?;
        let &id = coordinate
            .get(cell)
            .ok_or_else(|| Error::plain("Lee cochain has the wrong homological degree"))?;
        terms.push((
            id,
            checked_multiply(c.numerator, denominator / c.denominator)?,
        ));
    }
    model.spend(terms.len() * (1 + terms.len().max(1).ilog2() as usize))?;
    terms.sort_unstable_by_key(|x| x.0);
    Ok((terms, denominator))
}

fn coordinates<'a>(model: &mut Model, basis: &'a Basis) -> Result<BTreeMap<&'a Cell, usize>> {
    model.spend(basis.len() * (1 + basis.len().max(1).ilog2() as usize))?;
    storage(5 * basis.len())?;
    Ok(basis
        .iter()
        .enumerate()
        .map(|(i, (cell, _))| (cell, i))
        .collect())
}

#[derive(Default)]
struct NeighborIds {
    ids: BTreeMap<Cell, usize>,
    cells: usize,
}

impl NeighborIds {
    fn column(
        &mut self,
        model: &mut Model,
        context: &mut Context<'_>,
        cell: &Cell,
        reverse: bool,
    ) -> Result<Vector> {
        let mut column = Vec::new();
        context.lee_edges(model, cell, reverse, |model, adjacent, sign| {
            let size = self.ids.len();
            model.spend(adjacent.cells() * (1 + size.max(1).ilog2() as usize))?;
            let id = match self.ids.entry(adjacent) {
                Entry::Occupied(entry) => *entry.get(),
                Entry::Vacant(entry) => {
                    self.cells += entry.key().cells() + 2;
                    storage(self.cells)?;
                    *entry.insert(size)
                }
            };
            storage(4 * (column.len() + 1))?;
            column.push((id, sign));
            Ok(())
        })?;
        model.spend(column.len() * (1 + column.len().max(1).ilog2() as usize))?;
        column.sort_unstable_by_key(|x| x.0);
        debug_assert!(column.windows(2).all(|row| row[0].0 < row[1].0));
        Ok(column)
    }
}

struct Cut {
    j: i64,
    start: usize,
    end: usize,
    cycles: usize,
    projected: usize,
}

fn filtration_ranks(model: &mut Model, braid: &Braid, i: i64) -> Result<Vec<(i64, usize)>> {
    let mut context = Context::new(braid);
    let basis = context.lee_basis(model, i)?;
    let mut cuts = Vec::new();
    let mut start = 0;
    for same in basis.chunk_by(|a, b| a.1 == b.1) {
        let end = start + same.len();
        cuts.push(Cut {
            j: same[0].1,
            start,
            end,
            cycles: 0,
            projected: 0,
        });
        start = end;
    }
    {
        let mut outgoing = Echelon::default();
        let mut targets = NeighborIds::default();
        for cut in &mut cuts {
            for (cell, _) in &basis[cut.start..cut.end] {
                let column = targets.column(model, &mut context, cell, false)?;
                outgoing.insert(model, column, Vec::new())?;
            }
            cut.cycles = cut.end - outgoing.rank();
        }
    }
    // The rank of the projection of im d_in onto q < j is the rank
    // of the same rows of d_in. Adding transpose columns in ascending
    // quantum degree gives every cut without rebuilding the matrix.
    let mut incoming = Echelon::default();
    let mut sources = NeighborIds::default();
    for cut in cuts.iter_mut().rev() {
        cut.projected = incoming.rank();
        for (cell, _) in &basis[cut.start..cut.end] {
            let transpose = sources.column(model, &mut context, cell, true)?;
            incoming.insert(model, transpose, Vec::new())?;
        }
    }
    let full = incoming.rank();
    cuts.into_iter()
        .map(|cut| {
            let image = (cut.cycles + cut.projected)
                .checked_sub(full)
                .ok_or_else(|| Error::plain("Lee differential does not square to zero"))?;
            Ok((cut.j, image))
        })
        .collect()
}

pub(crate) fn rank(model: &mut Model, braid: &Braid, i: i64, f: i64) -> Result<usize> {
    if !braid.lee_possible(i, f) {
        return Ok(0);
    }
    Ok(filtration_ranks(model, braid, i)?
        .into_iter()
        .take_while(|&(j, _)| j >= f)
        .last()
        .map_or(0, |(_, rank)| rank))
}

pub(crate) fn levels(model: &mut Model, braid: &Braid, i: i64) -> Result<Vec<i64>> {
    let mut previous = 0;
    let mut levels = Vec::new();
    for (j, rank) in filtration_ranks(model, braid, i)? {
        let jump = rank
            .checked_sub(previous)
            .ok_or_else(|| Error::plain("Lee filtration ranks are not nested"))?;
        model.spend(1 + jump)?;
        levels.extend(std::iter::repeat_n(j, jump));
        previous = rank;
    }
    levels.reverse();
    Ok(levels)
}

pub(crate) fn rasmussen(model: &mut Model, braid: &Braid) -> Result<i64> {
    if braid.components(model)? != 1 {
        return Err(Error::plain("Rasmussen s is defined here for knots only"));
    }
    let levels = levels(model, braid, 0)?;
    match levels.as_slice() {
        &[lo, hi] if hi - lo == 2 => Ok((lo + hi) / 2),
        _ => Err(Error::plain(
            "knot Lee cohomology must have two surviving levels two apart",
        )),
    }
}

struct Boundaries {
    span: Echelon,
    sources: Vec<Cell>,
}

fn boundaries(
    model: &mut Model,
    context: &mut Context<'_>,
    basis: &Basis,
    tracked: bool,
) -> Result<Boundaries> {
    let mut columns = BTreeMap::<Cell, Vector>::new();
    let mut cells = 0;
    for (i, (cell, _)) in basis.iter().enumerate() {
        context.lee_edges(model, cell, true, |model, source, sign| {
            let size = columns.len();
            model.spend(source.cells() * (1 + size.max(1).ilog2() as usize))?;
            let column = match columns.entry(source) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => {
                    cells += entry.key().cells() + 5;
                    entry.insert(Vec::new())
                }
            };
            cells += 4;
            storage(cells)?;
            column.push((i, sign));
            Ok(())
        })?;
    }
    let mut span = Echelon::default();
    let mut sources = Vec::new();
    for (cell, column) in columns {
        debug_assert!(column.windows(2).all(|row| row[0].0 < row[1].0));
        let witness = if tracked {
            let coordinate = sources.len();
            sources.push(cell);
            vec![(coordinate, 1)]
        } else {
            Vec::new()
        };
        span.insert(model, column, witness)?;
    }
    Ok(Boundaries { span, sources })
}

fn scan_representatives(
    model: &mut Model,
    context: &mut Context<'_>,
    basis: &Basis,
    span: &mut Echelon,
    cutoff: i64,
    mut visit: impl FnMut(&mut Model, i64, &Vector, &Echelon) -> Result<bool>,
) -> Result<()> {
    let mut outgoing = Echelon::default();
    let mut targets = NeighborIds::default();
    for (i, (cell, j)) in basis
        .iter()
        .enumerate()
        .take_while(|(_, (_, j))| *j >= cutoff)
    {
        let column = targets.column(model, context, cell, false)?;
        let Some(cycle) = outgoing.insert(model, column, vec![(i, 1)])? else {
            continue;
        };
        // The kernel witness is already an actual closed cochain in
        // this filtration. Add it to the span of all boundaries and
        // earlier cycles before accepting it as a new quotient class.
        if span.insert(model, cycle.clone(), Vec::new())?.is_none()
            && visit(model, *j, &cycle, span)?
        {
            break;
        }
    }
    Ok(())
}

pub(crate) fn representative(
    model: &mut Model,
    braid: &Braid,
    i: i64,
    f: i64,
    k: usize,
) -> Result<Chain> {
    let mut context = Context::new(braid);
    let basis = context.lee_basis(model, i)?;
    let mut span = boundaries(model, &mut context, &basis, false)?.span;
    let mut found = None;
    let mut index = 0;
    scan_representatives(
        model,
        &mut context,
        &basis,
        &mut span,
        f,
        |model, j, cycle, _| {
            if j == f {
                if index == k {
                    let mut chain: Vec<_> = cycle
                        .iter()
                        .map(|&(row, c)| (basis[row].0.clone(), Rational::integer(c)))
                        .collect();
                    model.spend(chain.len() * (1 + chain.len().max(1).ilog2() as usize))?;
                    chain.sort_unstable_by(|a, b| a.0.cmp(&b.0));
                    found = Some(persist(model, chain)?);
                    return Ok(true);
                }
                index += 1;
            }
            Ok(false)
        },
    )?;
    found.ok_or_else(|| Error::plain("Lee generator index is out of range at this filtration jump"))
}

pub(crate) fn class_filtration(
    model: &mut Model,
    braid: &Braid,
    i: i64,
    chain: &Chain,
) -> Result<Option<i64>> {
    if chain.is_empty() {
        return Ok(None);
    }
    let mut context = Context::new(braid);
    require_cycle(model, &mut context, chain)?;
    let basis = context.lee_basis(model, i)?;
    let coord = coordinates(model, &basis)?;
    let (x, _) = to_integer(model, &coord, chain)?;
    let mut span = boundaries(model, &mut context, &basis, false)?.span;
    if span.spans(model, &x)? {
        return Ok(None);
    }
    let mut found = None;
    scan_representatives(
        model,
        &mut context,
        &basis,
        &mut span,
        i64::MIN,
        |model, j, _, span| {
            if span.spans(model, &x)? {
                found = Some(j);
                return Ok(true);
            }
            Ok(false)
        },
    )?;
    found
        .map(Some)
        .ok_or_else(|| Error::plain("closed Lee cochain is outside its cohomology span"))
}

pub(crate) fn independent(
    model: &mut Model,
    braid: &Braid,
    i: i64,
    chains: &[Chain],
) -> Result<bool> {
    let mut context = Context::new(braid);
    for chain in chains {
        require_cycle(model, &mut context, chain)?;
    }
    let basis = context.lee_basis(model, i)?;
    let coord = coordinates(model, &basis)?;
    let mut span = boundaries(model, &mut context, &basis, false)?.span;
    for chain in chains {
        let (x, _) = to_integer(model, &coord, chain)?;
        if span.insert(model, x, Vec::new())?.is_some() {
            return Ok(false);
        }
    }
    Ok(true)
}

pub(crate) fn null(model: &mut Model, braid: &Braid, i: i64, chain: &Chain) -> Result<Chain> {
    if chain.is_empty() {
        return Ok(Rc::default());
    }
    let mut context = Context::new(braid);
    require_cycle(model, &mut context, chain)?;
    let basis = context.lee_basis(model, i)?;
    let coord = coordinates(model, &basis)?;
    let (x, denominator) = to_integer(model, &coord, chain)?;
    let Boundaries { span, sources } = boundaries(model, &mut context, &basis, true)?;
    // Keep the coefficient of the input as the final tracked coordinate.
    // A zero image gives input_scale * x + d(witness) = 0.
    let sentinel = sources.len();
    let relation = span
        .relation(model, x, vec![(sentinel, 1)])?
        .ok_or_else(|| Error::plain("Lee cycle is not a boundary"))?;
    let Some((&(key, input_scale), witness)) = relation.split_last() else {
        return Err(Error::plain(
            "Lee boundary reduction lost its input coefficient",
        ));
    };
    if key != sentinel {
        return Err(Error::plain(
            "Lee boundary reduction lost its input coefficient",
        ));
    }
    let mut result = Vec::with_capacity(witness.len());
    for &(source, coefficient) in witness {
        model.spend(1)?;
        let coefficient =
            Rational::ratio(coefficient, input_scale)?.multiply(Rational::new(-1, denominator)?)?;
        result.push((sources[source].clone(), coefficient));
    }
    persist(model, result)
}

pub(crate) enum Command {
    Rank(crate::syntax::Label, Grade),
    Levels(crate::syntax::Label, i64),
    Rasmussen(crate::syntax::Label),
    Filtration(Box<crate::syntax::ParsedExpr>),
    Independent(Vec<crate::syntax::ParsedExpr>),
}

pub(crate) enum Query {
    Representative(crate::syntax::Label, Grade, usize),
    Operation(Operation, Box<crate::syntax::ParsedExpr>),
}
