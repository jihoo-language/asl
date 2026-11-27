use std::cmp::Ordering;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_phodge::valuation;

pub(crate) mod determinant;

const MAX_GROUND: usize = 64;

pub(crate) type Basis = Vec<usize>;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) prime: u32,
    pub(crate) ground: usize,
    pub(crate) rank: usize,
}

impl Grade {
    pub(crate) fn new(prime: u32, ground: usize, rank: usize) -> Result<Self> {
        if ground > MAX_GROUND {
            return Err(Error::plain(
                "tropical Grassmannian has a ground-set limit of 64",
            ));
        }
        if rank > ground {
            return Err(Error::plain("valuated matroid rank exceeds its ground set"));
        }
        Ok(Self {
            prime,
            ground,
            rank,
        })
    }

    pub(crate) fn dual(self) -> Self {
        Self {
            rank: self.ground - self.rank,
            ..self
        }
    }
}

pub(crate) struct Matroid {
    pub(crate) grade: Grade,
    values: Vec<(Basis, u64)>,
}

impl Matroid {
    pub(crate) fn build(
        model: &mut Model,
        grade: Grade,
        mut valuation: impl FnMut(&mut Model, &[usize]) -> Result<Option<i64>>,
    ) -> Result<Rc<Self>> {
        check_candidates(grade)?;
        let mut values = Vec::new();
        let mut basis: Vec<_> = (0..grade.rank).collect();
        let mut minimum = i64::MAX;
        loop {
            model.spend(1 + basis.len())?;
            if let Some(v) = valuation(model, &basis)? {
                minimum = minimum.min(v);
                check_cells(values.len() + 1, grade.rank)?;
                values.push((basis.clone(), v));
            }
            let pivot = (0..basis.len())
                .rev()
                .find(|&i| basis[i] < grade.ground - grade.rank + i);
            let Some(pivot) = pivot else {
                break;
            };
            let start = basis[pivot] + 1;
            for (i, entry) in basis[pivot..].iter_mut().enumerate() {
                *entry = start + i;
            }
        }
        if values.is_empty() {
            return Err(Error::plain("tropical Plücker data have no finite basis"));
        }
        let values = values
            .into_iter()
            .map(|(basis, v)| (basis, v.abs_diff(minimum)))
            .collect();
        Self::persist(model, grade, values)
    }

    fn persist(model: &mut Model, grade: Grade, values: Vec<(Basis, u64)>) -> Result<Rc<Self>> {
        let cells = check_cells(values.len(), grade.rank)?;
        model.allocate(cells)?;
        Ok(Rc::new(Self { grade, values }))
    }

    pub(crate) fn entries(&self) -> &[(Basis, u64)] {
        &self.values
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        for (a, b) in self.values.iter().zip(&other.values) {
            model.spend(self.grade.rank + 2)?;
            let order = a.cmp(b);
            if order != Ordering::Equal {
                return Ok(order);
            }
        }
        Ok(self.values.len().cmp(&other.values.len()))
    }

    fn lookup(&self, model: &mut Model, basis: &[usize]) -> Result<Option<u64>> {
        let search = 1 + self.values.len().ilog2() as usize;
        model.spend(search * (1 + basis.len()))?;
        Ok(self
            .values
            .binary_search_by(|(key, _)| key.as_slice().cmp(basis))
            .ok()
            .map(|i| self.values[i].1))
    }

    fn canonical_subset(&self, mut set: Basis, size: usize) -> Result<Basis> {
        if set.len() != size {
            return Err(Error::plain(format!(
                "expected a subset of size {size}, found {}",
                set.len()
            )));
        }
        set.sort_unstable();
        if set.last().is_some_and(|&i| i >= self.grade.ground) {
            return Err(Error::plain(
                "tropical subset contains an element outside its ground set",
            ));
        }
        if set.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(Error::plain("tropical subset contains a repeated element"));
        }
        Ok(set)
    }

    pub(crate) fn valuation(&self, model: &mut Model, basis: Basis) -> Result<Option<u64>> {
        let basis = self.canonical_subset(basis, self.grade.rank)?;
        self.lookup(model, &basis)
    }

    pub(crate) fn initial<'a>(
        &'a self,
        model: &mut Model,
        weights: &[i64],
    ) -> Result<(i128, Vec<&'a Basis>)> {
        if weights.len() != self.grade.ground {
            return Err(Error::plain(format!(
                "initial matroid needs {} weights, found {}",
                self.grade.ground,
                weights.len()
            )));
        }
        let mut cost = i128::MAX;
        let mut bases = Vec::new();
        for (key, v) in &self.values {
            model.spend(1 + key.len())?;
            let weighted =
                i128::from(*v) - key.iter().map(|&i| i128::from(weights[i])).sum::<i128>();
            if weighted < cost {
                bases.clear();
                cost = weighted;
            }
            if weighted == cost {
                bases.push(key);
            }
        }
        Ok((cost, bases))
    }

    pub(crate) fn dual(&self, model: &mut Model) -> Result<Rc<Self>> {
        let grade = self.grade.dual();
        check_cells(self.values.len(), grade.rank)?;
        let sorted_work = 1 + self.values.len().ilog2() as usize;
        model.spend(self.values.len() * (grade.ground + sorted_work * (1 + grade.rank)))?;
        let mut values = Vec::with_capacity(self.values.len());
        for (basis, val) in &self.values {
            let mut next = basis.iter().copied().peekable();
            let complement = (0..grade.ground)
                .filter(|point| {
                    if next.peek() == Some(point) {
                        next.next();
                        false
                    } else {
                        true
                    }
                })
                .collect();
            values.push((complement, *val));
        }
        values.sort_unstable();
        Self::persist(model, grade, values)
    }

    pub(crate) fn exchange(
        &self,
        model: &mut Model,
        left: Basis,
        right: Basis,
        removed: usize,
    ) -> Result<(usize, u128, u128)> {
        let left = self.canonical_subset(left, self.grade.rank)?;
        let right = self.canonical_subset(right, self.grade.rank)?;
        let Some(a) = self.lookup(model, &left)? else {
            return Err(Error::plain(
                "valuated exchange requires a finite left basis",
            ));
        };
        let Some(b) = self.lookup(model, &right)? else {
            return Err(Error::plain(
                "valuated exchange requires a finite right basis",
            ));
        };
        if left.binary_search(&removed).is_err() || right.binary_search(&removed).is_ok() {
            return Err(Error::plain(
                "exchange element must be in the left basis outside the right",
            ));
        }
        let old = u128::from(a) + u128::from(b);
        for &added in &right {
            model.spend(2 * self.grade.rank + 1)?;
            if left.binary_search(&added).is_ok() {
                continue;
            }
            let changed_left = exchange_one(&left, removed, added);
            let changed_right = exchange_one(&right, added, removed);
            if let (Some(u), Some(v)) = (
                self.lookup(model, &changed_left)?,
                self.lookup(model, &changed_right)?,
            ) {
                let changed = u128::from(u) + u128::from(v);
                if changed <= old {
                    return Ok((added, old, changed));
                }
            }
        }
        Err(Error::plain("valuated basis exchange has no witness"))
    }

    pub(crate) fn plucker_tie(
        &self,
        model: &mut Model,
        left: Basis,
        right: Basis,
    ) -> Result<(Option<u128>, Vec<usize>)> {
        if self.grade.rank == 0 || self.grade.rank == self.grade.ground {
            return Err(Error::plain(
                "tropical Plücker relations require rank strictly between zero and ground size",
            ));
        }
        let left = self.canonical_subset(left, self.grade.rank - 1)?;
        let right = self.canonical_subset(right, self.grade.rank + 1)?;
        let mut cost = None;
        let mut tied = Vec::new();
        for &added in &right {
            model.spend(2 * self.grade.rank + 1)?;
            if left.binary_search(&added).is_ok() {
                continue;
            }
            let u = insert_sorted(&left, added);
            let v: Basis = right.iter().copied().filter(|&i| i != added).collect();
            let a = self.lookup(model, &u)?;
            let b = self.lookup(model, &v)?;
            let current = a.zip(b).map(|(a, b)| u128::from(a) + u128::from(b));
            if current.is_some_and(|v| cost.is_none_or(|c| v < c)) {
                cost = current;
                tied.clear();
            }
            if current == cost {
                tied.push(added);
            }
        }
        if cost.is_some() && tied.len() < 2 {
            return Err(Error::plain("tropical Plücker minimum has no tie"));
        }
        Ok((cost, tied))
    }
}

fn insert_sorted(set: &[usize], added: usize) -> Basis {
    let split = set.partition_point(|&i| i < added);
    set[..split]
        .iter()
        .copied()
        .chain([added])
        .chain(set[split..].iter().copied())
        .collect()
}

fn exchange_one(basis: &[usize], removed: usize, added: usize) -> Basis {
    let mut changed: Basis = basis.iter().copied().filter(|&i| i != removed).collect();
    let place = changed.partition_point(|&i| i < added);
    changed.insert(place, added);
    changed
}

fn check_candidates(grade: Grade) -> Result<()> {
    let mut count = 1;
    for i in 0..grade.rank.min(grade.ground - grade.rank) {
        count = count * (grade.ground - i) / (i + 1);
        if count > MAX_PIVOT_CELLS {
            return Err(Error::plain(
                "tropical Grassmannian exceeds its Plücker-minor work limit",
            ));
        }
    }
    Ok(())
}

fn check_cells(count: usize, rank: usize) -> Result<usize> {
    let cells = 5 + count.saturating_mul(5 + rank);
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "tropical Plücker valuations exceed their storage limit",
        ));
    }
    Ok(cells)
}

pub(crate) fn rational(
    model: &mut Model,
    grade: Grade,
    mut rows: Vec<Vec<i128>>,
) -> Result<Rc<Matroid>> {
    check_rows(grade, &rows)?;
    // Clearing the scalar content of each row gives the same normalized vector.
    model.spend(grade.rank * grade.ground)?;
    for row in &mut rows {
        let content = row
            .iter()
            .fold(0, |d, x| crate::rv_ore::gcd(d, x.unsigned_abs()));
        if content > 1 {
            for x in row.iter_mut().filter(|x| **x != 0) {
                *x = crate::rv_phodge::divide_content(*x, content);
            }
        }
    }
    Matroid::build(model, grade, |model, basis| {
        let matrix = minor(model, &rows, basis)?;
        match integer_determinant(model, matrix)? {
            IntegerDeterminant::Exact(Some(det)) => Ok(Some(valuation(
                model,
                det.unsigned_abs(),
                u128::from(grade.prime),
            )?)),
            IntegerDeterminant::Exact(None) => Ok(None),
            IntegerDeterminant::Overflow => {
                let matrix = minor(model, &rows, basis)?;
                determinant::rational(model, matrix, Prime::validated(grade.prime))
                    .map_err(determinant::explain_limit)
            }
        }
    })
}

pub(crate) fn quadratic(
    model: &mut Model,
    grade: Grade,
    field: crate::rv_qp_quad::Field,
    rows: &[Vec<crate::rv_qp_quad::Element>],
) -> Result<Rc<Matroid>> {
    check_rows(grade, rows)?;
    Matroid::build(model, grade, |model, basis| {
        let matrix = minor(model, rows, basis)?;
        match quadratic_determinant_valuation(model, field, matrix) {
            Err(error) if crate::rv_qp_quad::is_coefficient_overflow(&error) => {
                let matrix = minor(model, rows, basis)?;
                determinant::quadratic(model, field, matrix).map_err(determinant::explain_limit)
            }
            result => result,
        }
    })
}

pub(crate) fn check_rows<T>(grade: Grade, rows: &[Vec<T>]) -> Result<()> {
    if rows.len() != grade.rank || rows.iter().any(|row| row.len() != grade.ground) {
        return Err(Error::plain(format!(
            "tropical row matrix must have {} rows of length {}",
            grade.rank, grade.ground,
        )));
    }
    Ok(())
}

fn minor<T: Copy>(model: &mut Model, rows: &[Vec<T>], basis: &[usize]) -> Result<Vec<Vec<T>>> {
    model.spend(basis.len() * basis.len())?;
    Ok(rows
        .iter()
        .map(|row| basis.iter().map(|&column| row[column]).collect())
        .collect())
}

enum IntegerDeterminant {
    Exact(Option<i128>),
    Overflow,
}

fn integer_determinant(
    model: &mut Model,
    mut matrix: Vec<Vec<i128>>,
) -> Result<IntegerDeterminant> {
    let rank = matrix.len();
    let mut previous = 1;
    for k in 0..rank {
        model.spend(rank - k)?;
        let Some(pivot_row) = (k..rank).find(|&row| matrix[row][k] != 0) else {
            return Ok(IntegerDeterminant::Exact(None));
        };
        matrix.swap(k, pivot_row);
        let pivot = matrix[k][k];
        let (above, below) = matrix.split_at_mut(k + 1);
        let pivot_row = &above[k];
        for row in below {
            model.spend(rank - k)?;
            let scale = row[k];
            for j in k + 1..rank {
                let Some(updated) = row[j]
                    .checked_mul(pivot)
                    .and_then(|lhs| {
                        scale
                            .checked_mul(pivot_row[j])
                            .and_then(|rhs| lhs.checked_sub(rhs))
                    })
                    .and_then(|numerator| numerator.checked_div(previous))
                else {
                    return Ok(IntegerDeterminant::Overflow);
                };
                row[j] = updated;
            }
        }
        // The valuation ignores the sign of row swaps.
        previous = pivot;
    }
    Ok(IntegerDeterminant::Exact(Some(previous)))
}

fn quadratic_determinant_valuation(
    model: &mut Model,
    field: crate::rv_qp_quad::Field,
    mut matrix: Vec<Vec<crate::rv_qp_quad::Element>>,
) -> Result<Option<i64>> {
    let rank = matrix.len();
    let mut result = 0;
    for k in 0..rank {
        model.spend(rank - k)?;
        let Some(pivot_row) = (k..rank).find(|&row| !matrix[row][k].is_zero()) else {
            return Ok(None);
        };
        matrix.swap(k, pivot_row);
        let pivot = matrix[k][k];
        // Exact elimination gives the minor as a signed product of pivots.
        result += field.valuation(model, pivot)?.expect("nonzero pivot");
        let (above, below) = matrix.split_at_mut(k + 1);
        let pivot_row = &above[k];
        let mut inverse = None;
        for row in below {
            model.spend(1)?;
            if row[k].is_zero() {
                continue;
            }
            let inv = match inverse {
                Some(value) => value,
                None => {
                    let value = field.inverse(model, pivot)?;
                    inverse = Some(value);
                    value
                }
            };
            let scale = field.multiply(model, row[k], inv)?;
            for j in k + 1..rank {
                let term = field.multiply(model, scale, pivot_row[j])?;
                row[j] = field.subtract(model, row[j], term)?;
            }
        }
    }
    Ok(Some(result))
}

pub(crate) enum Query {
    Plucker(Basis),
    Support,
    Initial(Vec<i64>),
    Exchange(Basis, Basis, usize),
    Tie(Basis, Basis),
}

impl Query {
    pub(crate) fn report(self, model: &mut Model, matroid: &Matroid) -> Result<String> {
        Ok(match self {
            Self::Plucker(basis) => {
                let value = matroid.valuation(model, basis)?;
                format!("rvValPlucker = {}", format_valuation(value))
            }
            Self::Support => {
                model.spend(matroid.values.len() * (matroid.grade.rank + 1))?;
                let bases: Vec<_> = matroid.values.iter().map(|(basis, _)| basis).collect();
                format!("rvValSupport = {bases:?}")
            }
            Self::Initial(weights) => {
                let (cost, bases) = matroid.initial(model, &weights)?;
                format!("rvValInitial = ({cost}, {bases:?})")
            }
            Self::Exchange(left, right, removed) => {
                let (added, old, new) = matroid.exchange(model, left, right, removed)?;
                format!("rvValExchange = ({added}, {old}, {new})")
            }
            Self::Tie(left, right) => {
                let (cost, ties) = matroid.plucker_tie(model, left, right)?;
                format!("rvValTie = ({}, {ties:?})", format_valuation(cost))
            }
        })
    }
}

fn format_valuation<T: ToString>(valuation: Option<T>) -> String {
    valuation.map_or_else(|| "infinity".to_owned(), |v| v.to_string())
}
