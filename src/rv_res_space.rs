use std::cmp::Ordering;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;
use crate::rv_qp_quad::gauge::Coefficients as FieldArithmetic;
use crate::rv_qp_quad::{Element, Field};
use crate::tropical_matroid;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) prime: Prime,
    pub(crate) width: usize,
}

impl Grade {
    pub(crate) fn new(prime: Prime, width: usize) -> Result<Self> {
        // Use the same coordinate bound as the tropical row spaces.
        tropical_matroid::Grade::new(prime.get(), width, 0)?;
        Ok(Self { prime, width })
    }

    pub(crate) fn count(self) -> usize {
        // Gaussian binomial recurrence, including the zero space.
        let mut counts = vec![0usize; self.width + 1];
        counts[0] = 1;
        for n in 1..=self.width {
            for k in (1..=n).rev() {
                let power = (self.prime.get() as usize).saturating_pow((n - k) as u32);
                counts[k] = counts[k].saturating_add(power.saturating_mul(counts[k - 1]));
            }
        }
        counts.into_iter().fold(0usize, usize::saturating_add)
    }
}

pub(crate) struct Space {
    pub(crate) grade: Grade,
    pub(crate) rows: Vec<Vec<u32>>,
}

impl Space {
    pub(crate) fn zero(grade: Grade) -> Self {
        Self {
            grade,
            rows: Vec::new(),
        }
    }

    pub(crate) fn entries(model: &mut Model, grade: Grade, rows: Vec<Vec<i128>>) -> Result<Self> {
        validate_rows(grade, &rows)?;
        model.spend(rows.len() * (grade.width + 1))?;
        let rows = rows
            .into_iter()
            .map(|row| row.into_iter().map(|x| residue(grade.prime, x)).collect())
            .collect();
        Self::new(model, grade, rows)
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(4 + self.rows.len() * (3 + self.grade.width))?;
        Ok(Value::RVResSpace(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        model.spend(self.rows.len().min(other.rows.len()) * (1 + self.grade.width))?;
        Ok(self.rows.cmp(&other.rows))
    }

    fn new(model: &mut Model, grade: Grade, mut rows: Vec<Vec<u32>>) -> Result<Self> {
        let mut i = 0;
        for j in 0..grade.width {
            model.spend(1 + rows.len() - i)?;
            let Some(pivot) = (i..rows.len()).find(|&k| rows[k][j] != 0) else {
                continue;
            };
            rows.swap(i, pivot);
            let c = invert(model, grade.prime, rows[i][j])?;
            model.spend(grade.width - j)?;
            for entry in &mut rows[i][j..] {
                *entry = grade.prime.mul(c, *entry);
            }
            let (above, rest) = rows.split_at_mut(i);
            let (pivot, below) = rest.split_first_mut().expect("a pivot row exists");
            for row in above.iter_mut().chain(below.iter_mut()) {
                let c = row[j];
                if c == 0 {
                    continue;
                }
                model.spend(grade.width - j)?;
                subtract_scaled(grade.prime, &mut row[j..], &pivot[j..], c);
            }
            i += 1;
        }
        rows.truncate(i);
        Ok(Self { grade, rows })
    }

    fn contains(&self, model: &mut Model, mut vector: Vec<u32>) -> Result<bool> {
        for row in &self.rows {
            model.spend(1 + self.grade.width)?;
            let j = pivot(row);
            let c = vector[j];
            if c != 0 {
                subtract_scaled(self.grade.prime, &mut vector[j..], &row[j..], c);
            }
        }
        Ok(vector.iter().all(|&x| x == 0))
    }

    pub(crate) fn contains_vector(&self, model: &mut Model, vector: Vec<i128>) -> Result<bool> {
        if vector.len() != self.grade.width {
            return Err(Error::plain(format!(
                "residue vector must have {} coordinates",
                self.grade.width
            )));
        }
        model.spend(vector.len())?;
        self.contains(
            model,
            vector
                .into_iter()
                .map(|x| residue(self.grade.prime, x))
                .collect(),
        )
    }

    pub(crate) fn is_subspace(&self, model: &mut Model, other: &Self) -> Result<bool> {
        model.spend(self.rows.len() * (self.grade.width + 1))?;
        for row in &self.rows {
            if !other.contains(model, row.clone())? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub(crate) fn sum(&self, model: &mut Model, other: &Self) -> Result<Self> {
        model.spend((self.rows.len() + other.rows.len()) * (self.grade.width + 1))?;
        Self::new(
            model,
            self.grade,
            self.rows.iter().chain(&other.rows).cloned().collect(),
        )
    }

    pub(crate) fn orthogonal(&self, model: &mut Model) -> Result<Self> {
        let width = self.grade.width;
        model.spend(width * (1 + self.rows.len()))?;
        let pivots: Vec<_> = self.rows.iter().map(|row| pivot(row)).collect();
        let mut rows = Vec::with_capacity(width - pivots.len());
        for j in 0..width {
            if pivots.binary_search(&j).is_err() {
                let mut vector = vec![0; width];
                vector[j] = 1;
                for (&i, row) in pivots.iter().zip(&self.rows) {
                    vector[i] = self.grade.prime.negate(row[j]);
                }
                rows.push(vector);
            }
        }
        Self::new(model, self.grade, rows)
    }

    pub(crate) fn intersection(&self, model: &mut Model, other: &Self) -> Result<Self> {
        let left = self.orthogonal(model)?;
        let right = other.orthogonal(model)?;
        left.sum(model, &right)?.orthogonal(model)
    }

    pub(crate) fn support(&self, model: &mut Model) -> Result<String> {
        let rank = self.rows.len();
        let grade = tropical_matroid::Grade::new(self.grade.prime.get(), self.grade.width, rank)?;
        let matroid = tropical_matroid::Matroid::build(model, grade, |model, basis| {
            model.spend(rank * (rank + 1))?;
            let rows = self
                .rows
                .iter()
                .map(|row| basis.iter().map(|&j| row[j]).collect())
                .collect();
            let minor = Self::new(
                model,
                Grade {
                    width: rank,
                    ..self.grade
                },
                rows,
            )?;
            Ok((minor.rows.len() == rank).then_some(0))
        })?;
        let support: Vec<_> = matroid.entries().iter().map(|(basis, _)| basis).collect();
        Ok(format!("rvResSupport = {support:?}"))
    }
}

// Called only after the interpreter checks the finite type's enumeration limit.
pub(crate) fn enumerate(model: &mut Model, grade: Grade, count: usize) -> Result<Vec<Value>> {
    let mut spaces = Vec::with_capacity(count);
    for mask in 0..(1u128 << grade.width) {
        let pivots: Vec<_> = (0..grade.width).filter(|&j| (mask >> j) & 1 != 0).collect();
        model.spend(1 + grade.width * (1 + pivots.len()))?;
        let free: Vec<_> = pivots
            .iter()
            .enumerate()
            .flat_map(|(i, &j)| {
                (j + 1..grade.width)
                    .filter(move |&k| (mask >> k) & 1 == 0)
                    .map(move |k| (i, k))
            })
            .collect();
        let choices = (grade.prime.get() as usize).pow(free.len() as u32);
        let mut rows = vec![vec![0; grade.width]; pivots.len()];
        for (i, &j) in pivots.iter().enumerate() {
            rows[i][j] = 1;
        }
        for choice in 0..choices {
            model.spend(1 + rows.len() * (1 + grade.width))?;
            let mut digits = choice;
            for &(i, j) in &free {
                let p = grade.prime.get() as usize;
                rows[i][j] = (digits % p) as u32;
                digits /= p;
            }
            spaces.push(
                Space {
                    grade,
                    rows: rows.clone(),
                }
                .persist(model)?,
            );
        }
    }
    debug_assert_eq!(spaces.len(), count);
    model.spend(count * (1 + count.ilog2() as usize) * (1 + grade.width * grade.width))?;
    spaces.sort_by(|x, y| {
        let (Value::RVResSpace(x), Value::RVResSpace(y)) = (x, y) else {
            unreachable!()
        };
        x.rows.cmp(&y.rows)
    });
    Ok(spaces)
}

fn subtract_scaled(prime: Prime, target: &mut [u32], row: &[u32], scale: u32) {
    let scale = prime.negate(scale);
    for (a, &b) in target.iter_mut().zip(row) {
        *a = prime.add(*a, prime.mul(scale, b));
    }
}

fn pivot(row: &[u32]) -> usize {
    row.iter()
        .position(|&x| x != 0)
        .expect("RREF stores only nonzero rows")
}

fn residue(prime: Prime, x: i128) -> u32 {
    x.rem_euclid(i128::from(prime.get())) as u32
}

pub(crate) fn validate_rows<T>(grade: Grade, rows: &[Vec<T>]) -> Result<()> {
    if rows.len().saturating_mul(1 + grade.width.saturating_mul(8)) > MAX_PIVOT_CELLS {
        return Err(Error::plain("residue row matrix exceeds the storage limit"));
    }
    if rows.iter().any(|row| row.len() != grade.width) {
        return Err(Error::plain(format!(
            "residue row matrix must have rows of length {}",
            grade.width
        )));
    }
    Ok(())
}

pub(crate) fn rational(
    model: &mut Model,
    grade: Grade,
    weights: &[i64],
    rows: Vec<Vec<(i128, i128)>>,
) -> Result<Space> {
    validate_rows(grade, &rows)?;
    let rows = rows
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|(a, c)| Element::new(a, 0, c))
                .collect()
        })
        .collect::<Result<Vec<_>>>()?;
    specialize(
        model,
        grade,
        weights,
        FieldArithmetic::Rational(grade.prime),
        rows,
    )
}

pub(crate) fn quadratic(
    model: &mut Model,
    grade: Grade,
    weights: &[i64],
    field: Field,
    rows: Vec<Vec<Element>>,
) -> Result<Space> {
    validate_rows(grade, &rows)?;
    specialize(
        model,
        grade,
        weights,
        FieldArithmetic::Quadratic(field),
        rows,
    )
}

fn specialize(
    model: &mut Model,
    grade: Grade,
    weights: &[i64],
    field: FieldArithmetic,
    mut rows: Vec<Vec<Element>>,
) -> Result<Space> {
    let width = grade.width;
    if weights.len() != width {
        return Err(Error::plain(format!(
            "residue specialization needs {width} weights, found {}",
            weights.len()
        )));
    }
    if let [row] = rows.as_slice() {
        return specialize_line(model, grade, weights, field, row);
    }
    let mut pivots = Vec::new();
    for i in 0..rows.len().min(width) {
        // These are the valuations in the matrix with column j scaled by
        // p^(-w_j). We never actually construct the powers of p.
        let Some((_, j, pivot_row)) = weighted_pivot(model, field, &rows, weights, i)? else {
            break;
        };
        rows.swap(i, pivot_row);
        normalize_pivot(model, field, &mut rows[i], j)?;
        let (above, rest) = rows.split_at_mut(i);
        let (pivot, below) = rest.split_first_mut().expect("a pivot row exists");
        for row in above.iter_mut().chain(below.iter_mut()) {
            eliminate_pivot(model, field, row, pivot, j)?;
        }
        pivots.push(j);
    }

    // In the chosen identity chart, every entry after scaling and pivot
    // normalization is integral. The pivot columns stay an identity.
    let mut reduced = Vec::with_capacity(pivots.len());
    for (row, pivot) in rows.into_iter().zip(pivots) {
        model.spend(1 + width)?;
        let mut result = vec![0; width];
        for (j, x) in row.into_iter().enumerate() {
            if let Some(order) = field.order(model, x)? {
                let adjusted =
                    i128::from(order) + i128::from(weights[pivot]) - i128::from(weights[j]);
                if adjusted < 0 {
                    return Err(Error::plain(
                        "weighted residue chart has a nonintegral entry",
                    ));
                }
                if adjusted == 0 {
                    result[j] = field.angle(model, x, order)?;
                }
            }
        }
        reduced.push(result);
    }
    Space::new(model, grade, reduced)
}

// In a one-row chart the reduced entries are angular ratios. There is
// no need to represent the exact coefficient ratios in the number field.
fn specialize_line(
    model: &mut Model,
    grade: Grade,
    weights: &[i64],
    field: FieldArithmetic,
    row: &[Element],
) -> Result<Space> {
    let mut locals = Vec::with_capacity(row.len());
    let mut minimum = None;
    for (j, &entry) in row.iter().enumerate() {
        model.spend(1)?;
        let local = field
            .order(model, entry)?
            .map(|order| (i128::from(order) - i128::from(weights[j]), order));
        if let Some((weight, order)) = local {
            if minimum.is_none_or(|(best, _, _)| weight < best) {
                minimum = Some((weight, j, order));
            }
        }
        locals.push(local);
    }
    let Some((weight, pivot, order)) = minimum else {
        return Ok(Space::zero(grade));
    };
    let angular = field.angle(model, row[pivot], order)?;
    let inverse = invert(model, grade.prime, angular)?;
    let mut result = vec![0; row.len()];
    for (j, local) in locals.into_iter().enumerate() {
        model.spend(1)?;
        if let Some((local_weight, order)) = local {
            if local_weight == weight {
                result[j] = grade.prime.mul(field.angle(model, row[j], order)?, inverse);
            }
        }
    }
    Space::new(model, grade, vec![result])
}

fn weighted_pivot(
    model: &mut Model,
    field: FieldArithmetic,
    rows: &[Vec<Element>],
    weights: &[i64],
    first: usize,
) -> Result<Option<(i128, usize, usize)>> {
    let mut best = None;
    for (i, row) in rows.iter().enumerate().skip(first) {
        model.spend(1 + weights.len())?;
        for (j, (&x, &w)) in row.iter().zip(weights).enumerate() {
            if let Some(order) = field.order(model, x)? {
                let candidate = (i128::from(order) - i128::from(w), j, i);
                if best.is_none_or(|old| candidate < old) {
                    best = Some(candidate);
                }
            }
        }
    }
    Ok(best)
}

fn normalize_pivot(
    model: &mut Model,
    field: FieldArithmetic,
    row: &mut [Element],
    j: usize,
) -> Result<()> {
    model.spend(1 + row.len())?;
    let pivot = row[j];
    row[j] = Element::ZERO;
    if pivot != Element::ONE && row.iter().any(|x| !x.is_zero()) {
        for x in row.iter_mut().filter(|x| !x.is_zero()) {
            *x = field.ratio(model, *x, pivot)?;
        }
    }
    row[j] = Element::ONE;
    Ok(())
}

fn eliminate_pivot(
    model: &mut Model,
    field: FieldArithmetic,
    row: &mut [Element],
    pivot: &[Element],
    j: usize,
) -> Result<()> {
    model.spend(1)?;
    if row[j].is_zero() {
        return Ok(());
    }
    model.spend(row.len())?;
    let scale = row[j];
    row[j] = Element::ZERO;
    for (k, (x, &y)) in row.iter_mut().zip(pivot).enumerate() {
        if k != j && !y.is_zero() {
            *x = field.subtract_product(model, *x, scale, y)?;
        }
    }
    Ok(())
}

pub(crate) enum Query {
    Dimension,
    Contains(Vec<i128>),
    Support,
}

impl Query {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            Self::Dimension => "rvResDim",
            Self::Contains(_) => "rvResContains",
            Self::Support => "rvResSupport",
        }
    }

    pub(crate) fn report(self, model: &mut Model, space: &Space) -> Result<String> {
        match self {
            Self::Dimension => Ok(format!("rvResDim = {}", space.rows.len())),
            Self::Contains(vector) => Ok(format!(
                "rvResContains = {}",
                space.contains_vector(model, vector)?
            )),
            Self::Support => space.support(model),
        }
    }
}
