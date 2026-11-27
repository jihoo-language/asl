use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_hive::{self, Boundary, MAX_RANK};

const MAX_CELLS: usize = 300_000;
pub(crate) type Shape = Rc<[usize]>;
pub(crate) type Support = Rc<[(Schubert, i128)]>;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) rank: usize,
    pub(crate) ambient: usize,
}

impl Grade {
    pub(crate) fn checked(rank: usize, ambient: usize) -> Result<Self> {
        if rank == 0 || rank >= ambient {
            return Err(Error::plain("quantum Grassmannian requires 0 < k < n"));
        }
        let grade = Self { rank, ambient };
        if rank.min(grade.width()) > MAX_RANK {
            return Err(Error::plain(format!(
                "quantum Grassmannian requires min(k,n-k) to be at most {MAX_RANK}"
            )));
        }
        Ok(grade)
    }

    pub(crate) fn complement(self) -> Self {
        Self {
            rank: self.width(),
            ..self
        }
    }

    pub(crate) fn width(self) -> usize {
        self.ambient - self.rank
    }

    pub(crate) fn shape(self, model: &mut Model, mut parts: Vec<usize>) -> Result<Shape> {
        if parts.len() > self.rank || parts.iter().any(|&part| part > self.width()) {
            return Err(Error::plain(
                "Schubert partition lies outside the Grassmannian rectangle",
            ));
        }
        model.spend(parts.len())?;
        if parts.windows(2).any(|pair| pair[0] < pair[1]) {
            return Err(Error::plain("Schubert partition must be weakly decreasing"));
        }
        parts.truncate(parts.partition_point(|&part| part != 0));
        storage(cells(parts.len()))?;
        Ok(parts.into())
    }
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Schubert {
    pub(crate) quantum: usize,
    pub(crate) shape: Shape,
}

fn arithmetic() -> Error {
    Error::plain("quantum Grassmannian integer arithmetic is out of range")
}

fn cells(rank: usize) -> usize {
    rank + 4
}

fn storage(count: usize) -> Result<()> {
    if count > MAX_CELLS {
        return Err(Error::plain(
            "quantum Grassmannian calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

pub(crate) fn persist(model: &mut Model, terms: Vec<(Schubert, i128)>) -> Result<Value> {
    let mut count = 0;
    for (term, _) in &terms {
        count += cells(term.shape.len());
        storage(count)?;
    }
    model.allocate(count)?;
    Ok(Value::RVQGr(terms.into()))
}

pub(crate) fn zero() -> Value {
    Value::RVQGr([].into())
}

pub(crate) fn basis(model: &mut Model, grade: Grade, parts: Vec<usize>) -> Result<Value> {
    let shape = grade.shape(model, parts)?;
    persist(model, vec![(Schubert { quantum: 0, shape }, 1)])
}

pub(crate) fn compare(model: &mut Model, left: &Support, right: &Support) -> Result<Ordering> {
    if Rc::ptr_eq(left, right) {
        return Ok(Ordering::Equal);
    }
    for (a, b) in left.iter().zip(right.iter()) {
        model.spend(cells(a.0.shape.len()))?;
        let order = a.cmp(b);
        if order != Ordering::Equal {
            return Ok(order);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

#[derive(Default)]
pub(crate) struct Sum {
    pub(crate) terms: BTreeMap<Schubert, i128>,
    pub(crate) cells: usize,
}

impl Sum {
    pub(crate) fn insert(
        &mut self,
        model: &mut Model,
        term: Schubert,
        value: i128,
        kept: usize,
    ) -> Result<()> {
        let width = cells(term.shape.len());
        model.spend(width * (1 + self.terms.len().max(1).ilog2() as usize))?;
        if value == 0 {
            return Ok(());
        }
        let coefficient = self
            .terms
            .get(&term)
            .copied()
            .unwrap_or(0)
            .checked_add(value)
            .ok_or_else(arithmetic)?;
        if coefficient == 0 {
            self.terms.remove(&term);
            self.cells -= width;
        } else {
            if !self.terms.contains_key(&term) {
                storage(kept + self.cells + width)?;
                self.cells += width;
            }
            self.terms.insert(term, coefficient);
        }
        Ok(())
    }
}

pub(crate) fn add(model: &mut Model, a: &Support, b: &Support) -> Result<Value> {
    let mut sum = Sum::default();
    for (term, coefficient) in a.iter().chain(b.iter()) {
        sum.insert(model, term.clone(), *coefficient, 0)?;
    }
    persist(model, sum.terms.into_iter().collect())
}

pub(crate) fn scale(model: &mut Model, a: &Support, factor: i128) -> Result<Value> {
    if factor == 0 {
        return Ok(zero());
    }
    model.spend(a.len())?;
    let terms = a
        .iter()
        .map(|(term, c)| Ok((term.clone(), c.checked_mul(factor).ok_or_else(arithmetic)?)))
        .collect::<Result<Vec<_>>>()?;
    persist(model, terms)
}

pub(crate) fn shift(model: &mut Model, a: &Support, power: usize) -> Result<Value> {
    model.spend(a.len())?;
    let terms = a
        .iter()
        .map(|(term, c)| {
            let quantum = term.quantum.checked_add(power).ok_or_else(arithmetic)?;
            Ok((
                Schubert {
                    quantum,
                    shape: term.shape.clone(),
                },
                *c,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    persist(model, terms)
}

pub(crate) fn classical(model: &mut Model, a: &Support) -> Result<Value> {
    let length = a.partition_point(|(term, _)| term.quantum == 0);
    model.spend(length)?;
    persist(model, a[..length].to_vec())
}

pub(crate) fn transpose_shape(model: &mut Model, shape: &[usize], kept: usize) -> Result<Shape> {
    let rows = shape.first().copied().unwrap_or(0);
    if rows > MAX_CELLS {
        return Err(Error::plain(
            "transposed Schubert partition exceeds the storage limit",
        ));
    }
    storage(kept + cells(rows))?;
    model.spend(rows * (1 + shape.len().max(1).ilog2() as usize))?;
    Ok((1..=rows)
        .map(|column| shape.partition_point(|&row| row >= column))
        .collect())
}

pub(crate) fn transpose(model: &mut Model, a: &Support) -> Result<Value> {
    let mut sum = Sum::default();
    for (term, c) in a.iter() {
        let shape = transpose_shape(model, &term.shape, sum.cells)?;
        sum.insert(
            model,
            Schubert {
                quantum: term.quantum,
                shape,
            },
            *c,
            0,
        )?;
    }
    persist(model, sum.terms.into_iter().collect())
}

fn rim_hook(model: &mut Model, grade: Grade, lambda: &[i128]) -> Result<Option<(Schubert, i128)>> {
    model.spend(grade.rank * grade.rank)?;
    // A classical output fits twice the rectangle width, so working in i128
    // also handles a width near usize::MAX without overflowing a beta number.
    let mut quantum = 0usize;
    let ambient = grade.ambient as i128;
    let mut residues = Vec::with_capacity(grade.rank);
    for (i, &part) in lambda.iter().enumerate() {
        let beta = part + (grade.rank - 1 - i) as i128;
        quantum += (beta / ambient) as usize;
        residues.push(((beta % ambient) as usize, i));
    }
    residues.sort_unstable_by(|a, b| b.cmp(a));
    if residues.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Ok(None);
    }
    let mut parity = ((grade.rank - 1) * quantum) % 2;
    for (i, &(_, original)) in residues.iter().enumerate() {
        for &(_, next) in &residues[i + 1..] {
            parity ^= usize::from(original > next);
        }
    }
    let shape = residues
        .into_iter()
        .enumerate()
        .map(|(i, (residue, _))| residue - (grade.rank - 1 - i))
        .take_while(|&part| part != 0)
        .collect();
    let sign = if parity == 0 { 1 } else { -1 };
    Ok(Some((Schubert { quantum, shape }, sign)))
}

fn basis_product(
    model: &mut Model,
    grade: Grade,
    left: &Shape,
    right: &Shape,
    kept: usize,
) -> Result<Sum> {
    let mut sum = Sum::default();
    if left.is_empty() || right.is_empty() {
        let term = Schubert {
            quantum: 0,
            shape: left.max(right).clone(),
        };
        sum.insert(model, term, 1, kept)?;
        return Ok(sum);
    }
    if grade.rank > grade.width() {
        let left = transpose_shape(model, left, kept)?;
        let right = transpose_shape(model, right, kept + cells(left.len()))?;
        for (term, c) in basis_product(model, grade.complement(), &left, &right, kept)?.terms {
            let shape = transpose_shape(model, &term.shape, kept + sum.cells)?;
            let term = Schubert { shape, ..term };
            sum.insert(model, term, c, kept)?;
        }
        return Ok(sum);
    }
    if grade.rank == 1 {
        // In the projective-space ring, h^n=q.
        let degree = left[0] as u128 + right[0] as u128;
        let ambient = grade.ambient as u128;
        let part = (degree % ambient) as usize;
        let shape = if part == 0 { [].into() } else { [part].into() };
        let term = Schubert {
            quantum: (degree / ambient) as usize,
            shape,
        };
        sum.insert(model, term, 1, kept)?;
        return Ok(sum);
    }
    let left = padded(left, grade.rank);
    let right = padded(right, grade.rank);
    let a: Rc<[i128]> = left.iter().rev().copied().collect();
    let b: Rc<[i128]> = right.iter().rev().copied().collect();
    rv_hive::tensor_candidates(model, &left, &right, |model, shape| {
        let Some((term, sign)) = rim_hook(model, grade, shape)? else {
            return Ok(());
        };
        let c = shape.iter().rev().copied().collect();
        let boundary = Boundary([a.clone(), b.clone(), c]);
        let coefficient =
            i128::try_from(rv_hive::count(model, &boundary)?).map_err(|_| arithmetic())?;
        sum.insert(model, term, sign * coefficient, kept)
    })?;
    if sum.terms.values().any(|&coefficient| coefficient < 0) {
        return Err(Error::plain(
            "quantum rim-hook reduction produced a negative basis multiplicity",
        ));
    }
    Ok(sum)
}

fn padded(shape: &[usize], rank: usize) -> Vec<i128> {
    let mut result: Vec<_> = shape.iter().map(|&part| part as i128).collect();
    result.resize(rank, 0);
    result
}

pub(crate) fn multiply(model: &mut Model, grade: Grade, a: &Support, b: &Support) -> Result<Value> {
    let mut sum = Sum::default();
    for (x, c) in a.iter() {
        for (y, d) in b.iter() {
            model.spend(cells(x.shape.len()) + cells(y.shape.len()))?;
            let power = x.quantum.checked_add(y.quantum).ok_or_else(arithmetic)?;
            let factor = c.checked_mul(*d).ok_or_else(arithmetic)?;
            let product = basis_product(model, grade, &x.shape, &y.shape, sum.cells)?;
            for (mut term, coefficient) in product.terms {
                term.quantum = term.quantum.checked_add(power).ok_or_else(arithmetic)?;
                let coefficient = factor.checked_mul(coefficient).ok_or_else(arithmetic)?;
                sum.insert(model, term, coefficient, 0)?;
            }
        }
    }
    persist(model, sum.terms.into_iter().collect())
}

pub(crate) fn coefficient(
    model: &mut Model,
    grade: Grade,
    a: &Support,
    quantum: usize,
    parts: Vec<usize>,
) -> Result<i128> {
    let key = Schubert {
        quantum,
        shape: grade.shape(model, parts)?,
    };
    model.spend(cells(key.shape.len()) * (1 + a.len().max(1).ilog2() as usize))?;
    Ok(a.binary_search_by(|(term, _)| term.cmp(&key))
        .map(|i| a[i].1)
        .unwrap_or(0))
}

pub(crate) fn gw(
    model: &mut Model,
    mut grade: Grade,
    quantum: usize,
    shapes: [Vec<usize>; 3],
) -> Result<i128> {
    let [a, b, c] = shapes.map(|parts| grade.shape(model, parts));
    let (mut a, mut b, mut c) = (a?, b?, c?);
    let dimension = grade.rank as u128 * grade.width() as u128;
    let size: u128 = a
        .iter()
        .chain(b.iter())
        .chain(c.iter())
        .map(|&part| part as u128)
        .sum();
    if dimension.checked_add(quantum as u128 * grade.ambient as u128) != Some(size) {
        return Ok(0);
    }
    if grade.rank > grade.width() {
        a = transpose_shape(model, &a, 0)?;
        b = transpose_shape(model, &b, 0)?;
        c = transpose_shape(model, &c, 0)?;
        grade = grade.complement();
    }
    let dual: Shape = (0..grade.rank)
        .rev()
        .map(|i| grade.width() - c.get(i).copied().unwrap_or(0))
        .take_while(|&part| part != 0)
        .collect();
    let key = Schubert {
        quantum,
        shape: dual,
    };
    Ok(basis_product(model, grade, &a, &b, 0)?
        .terms
        .get(&key)
        .copied()
        .unwrap_or(0))
}

// Quantum Pieri for the hyperplane. In beta coordinates the eigenvalue
// equation moves one beta number up by one modulo n. On wraparound its
// reordering sign cancels the sign (-1)^(k-1) attached to q.
pub(crate) fn hyperplane(
    model: &mut Model,
    grade: Grade,
    shape: &[usize],
    kept: usize,
) -> Result<Vec<Schubert>> {
    if grade.rank > grade.width() {
        let shape = transpose_shape(model, shape, kept)?;
        let mut terms = hyperplane(model, grade.complement(), &shape, kept)?;
        let mut cells_kept = kept;
        for term in &mut terms {
            term.shape = transpose_shape(model, &term.shape, cells_kept)?;
            cells_kept += cells(term.shape.len());
        }
        return Ok(terms);
    }
    model.spend(grade.rank * grade.rank)?;
    let beta: Vec<_> = (0..grade.rank)
        .map(|i| shape.get(i).copied().unwrap_or(0) + grade.rank - 1 - i)
        .collect();
    let mut terms = Vec::new();
    for (i, &height) in beta.iter().enumerate() {
        let wrap = height == grade.ambient - 1;
        let next = if wrap { 0 } else { height + 1 };
        if (wrap && beta.last() == Some(&0)) || (i > 0 && beta[i - 1] == next) {
            continue;
        }
        let mut changed = beta.clone();
        changed[i] = next;
        if wrap {
            changed.rotate_left(1);
        }
        let shape = beta_shape(&changed);
        terms.push(Schubert {
            quantum: usize::from(wrap),
            shape,
        });
    }
    Ok(terms)
}

pub(crate) fn beta_shape(descending: &[usize]) -> Shape {
    descending
        .iter()
        .enumerate()
        .map(|(i, &height)| height - (descending.len() - 1 - i))
        .take_while(|&part| part != 0)
        .collect()
}
