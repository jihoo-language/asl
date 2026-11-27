use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::elliptic_modular::decimal;
use crate::error::{Error, Result};
use crate::model::Model;

use super::dimension::dimensions;
use super::local::{Denominator, Factor, LocalAlgebra};
use super::recoupling::{Block, Key};
use super::{
    accumulate, coefficient_overflow, diagram, evaluate_polynomial, intermediate_overflow, natural,
    Algebra, Binary, Diagram, Integer, Unary,
};

const LITERAL_RANK: usize = 3;
const COUNT_RANK: usize = 16;
pub(super) type Tail = Vec<u8>;
type Eigenvalue = [i128; 3];

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Path(Vec<Tail>);

impl Path {
    pub(crate) fn new(rank: usize, tails: Vec<Vec<u128>>) -> Result<Self> {
        require_literal_rank(rank)?;
        if tails.len() != 2 * rank + 1 || !tails[0].is_empty() {
            return Err(Error::plain(
                "Deligne path must start at the empty tail and have 2r+1 vertices",
            ));
        }
        let mut out = vec![Vec::new()];
        for (i, tail) in tails.into_iter().enumerate().skip(1) {
            let tail = checked_tail(rank, tail).ok_or_else(|| {
                Error::plain(format!("invalid Deligne path tail at doubled level {i}"))
            })?;
            if !children(out.last().unwrap(), i).contains(&tail) {
                return Err(Error::plain(format!(
                    "invalid Deligne path branch at doubled level {i}"
                )));
            }
            out.push(tail);
        }
        Ok(Self(out))
    }

    pub(crate) fn rank(&self) -> usize {
        (self.0.len() - 1) / 2
    }

    pub(super) fn vertices(&self) -> &[Tail] {
        &self.0
    }

    fn eigenvalues_at(&self, model: &mut Model, n: u128) -> Result<String> {
        let mut values = Vec::new();
        for (level, tail) in self.0.iter().enumerate() {
            let poly = eigenvalue(tail, level)
                .into_iter()
                .enumerate()
                .filter(|(_, v)| *v != 0)
                .map(|(i, v)| (i as u128, Integer::from_i128(v)))
                .collect();
            let mut value = evaluate_polynomial(model, poly, n)?;
            let (magnitude, remainder) = value.magnitude.divide_u32(2);
            if remainder != 0 {
                return Err(Error::plain(
                    "Deligne transposition eigenvalue is nonintegral",
                ));
            }
            value.magnitude = magnitude;
            values.push(decimal(value));
        }
        Ok(format!(
            "rvDelignePathEigenvalues[{}, {n}] = [{}]",
            self.rank(),
            values.join(", ")
        ))
    }
}

pub(super) fn checked_tail(rank: usize, tail: Vec<u128>) -> Option<Tail> {
    (tail.len() <= rank
        && tail.iter().all(|&row| row > 0 && row <= rank as u128)
        && tail.windows(2).all(|w| w[0] >= w[1]))
    .then(|| tail.into_iter().map(|row| row as u8).collect())
}

pub(super) fn terminal_paths(model: &mut Model, rank: usize, tail: &[u8]) -> Result<Vec<Path>> {
    let mut states = vec![vec![vec![]]];
    for level in 1..=2 * rank {
        let mut next = Vec::new();
        for path in &states {
            let children = children(path.last().unwrap(), level);
            model.spend(path.len() * children.len())?;
            for child in children {
                let mut extended = path.clone();
                extended.push(child);
                next.push(extended);
            }
        }
        states = next;
    }
    states.retain(|path| path.last().is_some_and(|row| row == tail));
    states.sort_unstable();
    if states.is_empty() {
        return Err(Error::plain("Deligne terminal has no vacillating paths"));
    }
    Ok(states.into_iter().map(Path).collect())
}

pub(crate) enum Report {
    Counts(usize),
    Eigenvalues(Path, u128),
}

impl Report {
    pub(crate) fn execute(self, model: &mut Model) -> Result<String> {
        match self {
            Self::Eigenvalues(path, n) => path.eigenvalues_at(model, n),
            Self::Counts(rank) => path_counts(model, rank),
        }
    }
}

pub(crate) fn require_literal_rank(rank: usize) -> Result<()> {
    if rank > LITERAL_RANK {
        return Err(Error::plain(
            "literal Deligne path diagrams are limited to rank 3",
        ));
    }
    Ok(())
}

pub(crate) struct Tower {
    rank: usize,
    operators: Vec<LocalAlgebra>,
    projectors: RefCell<BTreeMap<Vec<Tail>, LocalAlgebra>>,
    blocks: RefCell<BTreeMap<Tail, Rc<Block>>>,
}

impl Tower {
    pub(crate) fn new(model: &mut Model, rank: usize) -> Result<Self> {
        require_literal_rank(rank)?;
        let mut operators = Vec::new();
        for j in 0..=2 * rank {
            let current = operator(model, rank, j)?;
            if current.unary(model, Unary::Reflect)? != current {
                return Err(Error::plain(
                    "Deligne transposition tower operator is not self-adjoint",
                ));
            }
            for old in &operators {
                let left = current.binary(model, Binary::Compose, old)?;
                let right = old.binary(model, Binary::Compose, &current)?;
                if left != right {
                    return Err(Error::plain(
                        "Deligne transposition tower operators do not commute",
                    ));
                }
            }
            current.charge(model)?;
            operators.push(current);
        }
        if operators[0] != scalar(model, rank, eigenvalue(&[], 0))? {
            return Err(Error::plain(
                "Deligne tower base transposition operator is not scalar",
            ));
        }
        Ok(Self {
            rank,
            operators,
            projectors: RefCell::default(),
            blocks: RefCell::default(),
        })
    }

    pub(crate) fn block(&self, model: &mut Model, key: &Key) -> Result<Rc<Block>> {
        debug_assert_eq!(self.rank, key.rank);
        if let Some(block) = self.blocks.borrow().get(&key.tail) {
            return Ok(block.clone());
        }
        let block = Rc::new(Block::new(model, self, key)?);
        self.blocks
            .borrow_mut()
            .insert(key.tail.clone(), block.clone());
        Ok(block)
    }

    pub(crate) fn operator(&self, level: usize) -> Result<LocalAlgebra> {
        self.operators
            .get(level)
            .cloned()
            .ok_or_else(|| Error::plain("Deligne tower level exceeds twice the rank"))
    }

    pub(crate) fn projector(&self, model: &mut Model, path: &Path) -> Result<LocalAlgebra> {
        let mut current = LocalAlgebra::identity(self.rank);
        for level in 1..path.0.len() {
            let prefix = &path.0[..=level];
            let cached = self.projectors.borrow().get(prefix).cloned();
            if let Some(cached) = cached {
                current = cached;
                continue;
            }
            current = self.step(model, current, &path.0[level - 1], &path.0[level], level)?;
            current.charge(model)?;
            self.projectors
                .borrow_mut()
                .insert(prefix.to_vec(), current.clone());
        }
        Ok(current)
    }

    fn step(
        &self,
        model: &mut Model,
        mut value: LocalAlgebra,
        parent: &[u8],
        target: &[u8],
        level: usize,
    ) -> Result<LocalAlgebra> {
        let wanted = eigenvalue(target, level);
        let op = &self.operators[level];
        for sibling in children(parent, level)
            .into_iter()
            .filter(|tail| tail != target)
        {
            let eig = eigenvalue(&sibling, level);
            let factor = subtract_scalar(model, op, self.rank, eig)?;
            let composed = value.binary(model, Binary::Compose, &factor)?;
            value = divide_difference(model, composed, wanted, eig)?;
        }
        let difference = subtract_scalar(model, op, self.rank, wanted)?;
        if !value.binary(model, Binary::Compose, &difference)?.is_zero() || value.is_zero() {
            return Err(Error::plain(format!(
                "Deligne path eigenprojector failed at doubled level {level}"
            )));
        }
        Ok(value)
    }
}

fn operator(model: &mut Model, rank: usize, level: usize) -> Result<LocalAlgebra> {
    let mut terms = BTreeMap::new();
    diagram::all_patterns(
        model,
        2 * rank,
        (2 * rank) as u128,
        &mut |model, pattern| {
            let Some(value) = orbit_value(rank, level, pattern) else {
                return Ok(true);
            };
            let diagram = Diagram::literal(rank, pattern.to_vec())?;
            for ((basis, _), mu) in diagram.orbit(model)? {
                for (degree, &coefficient) in value.iter().enumerate().filter(|(_, c)| **c != 0) {
                    model.spend(1)?;
                    accumulate(
                        &mut terms,
                        (basis.clone(), degree as u128),
                        Integer::product(mu, coefficient),
                    )?;
                }
            }
            Ok(true)
        },
    )?;
    let terms = terms
        .into_iter()
        .map(|(key, value)| Ok((key, value.signed().ok_or_else(coefficient_overflow)?)))
        .collect::<Result<_>>()?;
    LocalAlgebra::embed(Algebra { rank, terms })?.divide(model, Denominator::new([1, 0, 0, 0])?)
}

// Twice the matrix entry, as a polynomial in t. The fixed suffix is
// compared even though the transpositions only act on the active slots.
fn orbit_value(rank: usize, level: usize, pattern: &[u128]) -> Option<Eigenvalue> {
    let active = level / 2;
    let (upper, lower) = pattern.split_at(rank);
    if upper[active..] != lower[active..] {
        return None;
    }
    let avoid = if level.is_multiple_of(2) {
        None
    } else {
        Some(lower[active])
    };
    if let Some((&a, &b)) = lower[..active]
        .iter()
        .zip(&upper[..active])
        .find(|(a, b)| a != b)
    {
        if avoid == Some(a) || avoid == Some(b) {
            return None;
        }
        let correct = lower[..active]
            .iter()
            .zip(&upper[..active])
            .all(|(&from, &to)| {
                to == if from == a {
                    b
                } else if from == b {
                    a
                } else {
                    from
                }
            });
        return correct.then_some([2, 0, 0]);
    }
    let used = lower[..active]
        .iter()
        .copied()
        .chain(avoid)
        .fold(0u16, |bits, letter| bits | (1 << letter));
    Some(binomial(used.count_ones() as i128))
}

fn binomial(offset: i128) -> Eigenvalue {
    [offset * (offset + 1), -(2 * offset + 1), 1]
}

fn eigenvalue(tail: &[u8], level: usize) -> Eigenvalue {
    let weight: i128 = tail.iter().map(|&row| i128::from(row)).sum();
    let mut polynomial = binomial(weight + i128::from(!level.is_multiple_of(2)));
    let contents: i128 = tail
        .iter()
        .enumerate()
        .map(|(i, &row)| {
            let row = i128::from(row);
            row * (row - 1) / 2 - i as i128 * row
        })
        .sum();
    polynomial[0] += 2 * (contents - weight);
    polynomial
}

fn scalar(model: &mut Model, rank: usize, polynomial: Eigenvalue) -> Result<LocalAlgebra> {
    let diagram = Diagram::identity(rank);
    let terms = polynomial
        .into_iter()
        .enumerate()
        .filter(|(_, v)| *v != 0)
        .map(|(i, v)| ((diagram.clone(), i as u128), v))
        .collect();
    LocalAlgebra::embed(Algebra { rank, terms })?.divide(model, Denominator::new([1, 0, 0, 0])?)
}

fn subtract_scalar(
    model: &mut Model,
    op: &LocalAlgebra,
    rank: usize,
    eigenvalue: Eigenvalue,
) -> Result<LocalAlgebra> {
    let negative = scalar(model, rank, eigenvalue.map(|x| -x))?;
    op.add(model, &negative)
}

fn divide_difference(
    model: &mut Model,
    value: LocalAlgebra,
    a: Eigenvalue,
    b: Eigenvalue,
) -> Result<LocalAlgebra> {
    let raw = std::array::from_fn::<_, 3, _>(|i| a[i] - b[i]);
    if raw[2] != 0 || raw[..2].iter().any(|x| x % 2 != 0) {
        return Err(Error::plain(
            "invalid Deligne sibling eigenvalue difference",
        ));
    }
    let [constant, slope, _] = raw.map(|x| x / 2);
    let (sign, factor) = match (slope, constant) {
        (0, c @ (-1 | 1)) => (c, None),
        (0, c @ (-3 | -2 | 2 | 3)) => (c.signum(), Some(Factor::Prime(c.unsigned_abs() as u8))),
        (m @ (-1 | 1), c) => {
            let root = -c * m;
            let root =
                u8::try_from(root).map_err(|_| Error::plain("invalid Deligne sibling root"))?;
            (m, Some(Factor::Root(root)))
        }
        _ => {
            return Err(Error::plain(
                "Deligne sibling eigenvalues are not separated by the localization",
            ))
        }
    };
    let value = if sign < 0 {
        value.unary(model, Unary::Scale(-1, 0))?
    } else {
        value
    };
    if let Some(factor) = factor {
        value.divide(model, Denominator::factor(factor)?)
    } else {
        Ok(value)
    }
}

fn children(tail: &[u8], level: usize) -> Vec<Tail> {
    let mut states = vec![tail.to_vec()];
    if level.is_multiple_of(2) {
        for i in 0..tail.len() {
            if i == 0 || tail[i - 1] > tail[i] {
                let mut child = tail.to_vec();
                child[i] += 1;
                states.push(child);
            }
        }
        let mut child = tail.to_vec();
        child.push(1);
        states.push(child);
    } else {
        for (i, &row) in tail.iter().enumerate() {
            if tail.get(i + 1).is_none_or(|&next| next < row) {
                let mut child = tail.to_vec();
                child[i] -= 1;
                if row == 1 {
                    child.pop();
                }
                states.push(child);
            }
        }
    }
    states
}

fn path_counts(model: &mut Model, rank: usize) -> Result<String> {
    if rank > COUNT_RANK {
        return Err(Error::plain(
            "Deligne path multiplicities are limited to rank 16",
        ));
    }
    let mut states = BTreeMap::from([(Vec::new(), 1u128)]);
    for level in 1..=2 * rank {
        let mut next = BTreeMap::new();
        for (tail, count) in states {
            model.spend(1 + tail.len().pow(2))?;
            for child in children(&tail, level) {
                let slot = next.entry(child).or_insert(0u128);
                *slot = slot.checked_add(count).ok_or_else(intermediate_overflow)?;
            }
        }
        states = next;
    }
    let mut squares = natural(0);
    let mut total = 0u128;
    for &count in states.values() {
        total = total.checked_add(count).ok_or_else(intermediate_overflow)?;
        let value = natural(count);
        squares = value
            .multiply_wide(value)
            .and_then(|q| squares.add(q))
            .map_err(|_| intermediate_overflow())?;
    }
    let (generic, _) = dimensions(model, 2 * rank, 0)?;
    if generic.magnitude != squares.magnitude {
        return Err(Error::plain(
            "Deligne path multiplicities do not exhaust the partition algebra",
        ));
    }
    let terms = states
        .into_iter()
        .map(|(tail, count)| format!("({tail:?},{count})"))
        .collect::<Vec<_>>();
    Ok(format!(
        "rvDelignePathCounts[{rank}] = paths={total}, dimension={}, multiplicities=[{}]",
        decimal(squares),
        terms.join(", ")
    ))
}
