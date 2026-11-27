use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_hive::{self, Boundary};

const MAX_CELLS: usize = 300_000;

pub(crate) type Weight = Rc<[i128]>;
pub(crate) type Support = Rc<[(Weight, i128)]>;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) rank: usize,
    pub(crate) level: usize,
}

impl Grade {
    pub(crate) fn checked(rank: usize, level: usize) -> Result<Self> {
        if !(2..=rv_hive::MAX_RANK).contains(&rank) {
            return Err(Error::plain(format!(
                "fusion rank must be between 2 and {}",
                rv_hive::MAX_RANK
            )));
        }
        Ok(Self { rank, level })
    }

    pub(crate) fn normalize(self, model: &mut Model, mut weight: Vec<i128>) -> Result<Weight> {
        if weight.len() != self.rank {
            return Err(Error::plain("fusion weight must match the ring rank"));
        }
        model.spend(self.rank)?;
        if weight.windows(2).any(|pair| pair[0] > pair[1]) {
            return Err(Error::plain("fusion weight must be weakly increasing"));
        }
        let center = weight[0];
        for entry in &mut weight {
            *entry = entry.checked_sub(center).ok_or_else(arithmetic)?;
        }
        if weight[self.rank - 1] > self.level as i128 {
            return Err(Error::plain("fusion weight lies outside its level alcove"));
        }
        Ok(weight.into())
    }
}

fn arithmetic() -> Error {
    Error::plain("fusion integer arithmetic is out of range")
}

fn size(rank: usize) -> usize {
    2 * rank + 3
}

fn storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain("fusion calculation exceeds the storage limit"));
    }
    Ok(())
}

fn persist(model: &mut Model, terms: Vec<(Weight, i128)>, rank: usize) -> Result<Value> {
    let cells = size(rank) * terms.len();
    storage(cells)?;
    model.allocate(cells)?;
    Ok(Value::RVFusion(terms.into()))
}

pub(crate) fn zero() -> Value {
    Value::RVFusion([].into())
}

pub(crate) fn basis(model: &mut Model, grade: Grade, weight: Vec<i128>) -> Result<Value> {
    let weight = grade.normalize(model, weight)?;
    persist(model, vec![(weight, 1)], grade.rank)
}

pub(crate) fn compare(model: &mut Model, left: &Support, right: &Support) -> Result<Ordering> {
    if Rc::ptr_eq(left, right) {
        return Ok(Ordering::Equal);
    }
    for (a, b) in left.iter().zip(right.iter()) {
        model.spend(1 + a.0.len())?;
        let order = a.cmp(b);
        if order != Ordering::Equal {
            return Ok(order);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

#[derive(Default)]
struct Sum(BTreeMap<Weight, i128>);

impl Sum {
    fn add(&mut self, model: &mut Model, weight: &Weight, factor: i128, kept: usize) -> Result<()> {
        let cells = size(weight.len());
        model.spend(cells * (1 + self.0.len().max(1).ilog2() as usize))?;
        if factor == 0 {
            return Ok(());
        }
        let value = self
            .0
            .get(weight)
            .copied()
            .unwrap_or(0)
            .checked_add(factor)
            .ok_or_else(arithmetic)?;
        if value == 0 {
            self.0.remove(weight);
        } else {
            if !self.0.contains_key(weight) {
                storage(kept + cells * (1 + self.0.len()))?;
            }
            self.0.insert(weight.clone(), value);
        }
        Ok(())
    }
}

pub(crate) fn add(model: &mut Model, grade: Grade, a: &Support, b: &Support) -> Result<Value> {
    let mut output = Sum::default();
    for (weight, coefficient) in a.iter().chain(b.iter()) {
        output.add(model, weight, *coefficient, 0)?;
    }
    persist(model, output.0.into_iter().collect(), grade.rank)
}

pub(crate) fn scale(model: &mut Model, grade: Grade, support: &Support, c: i128) -> Result<Value> {
    if c == 0 {
        return Ok(zero());
    }
    model.spend(support.len())?;
    let terms = support
        .iter()
        .map(|(weight, coefficient)| {
            Ok((
                weight.clone(),
                coefficient.checked_mul(c).ok_or_else(arithmetic)?,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    persist(model, terms, grade.rank)
}

pub(crate) fn dual_weight(weight: &[i128]) -> Weight {
    let width = *weight.last().unwrap();
    weight.iter().rev().map(|&entry| width - entry).collect()
}

pub(crate) fn dual(model: &mut Model, grade: Grade, support: &Support) -> Result<Value> {
    model.spend(size(grade.rank) * support.len() * (1 + support.len().max(1).ilog2() as usize))?;
    let mut terms: Vec<_> = support
        .iter()
        .map(|(weight, coefficient)| (dual_weight(weight), *coefficient))
        .collect();
    terms.sort_unstable();
    persist(model, terms, grade.rank)
}

pub(crate) fn coefficient(
    model: &mut Model,
    grade: Grade,
    support: &Support,
    weight: Vec<i128>,
) -> Result<i128> {
    let weight = grade.normalize(model, weight)?;
    model.spend(size(grade.rank) * (1 + support.len().max(1).ilog2() as usize))?;
    Ok(support
        .binary_search_by(|(entry, _)| entry.cmp(&weight))
        .map(|index| support[index].1)
        .unwrap_or(0))
}

// The inputs below are normalized alcove weights, so their entries are at most
// usize::MAX and their rank is at most MAX_RANK. Sums/products of coordinates
// in the partition and folding calculations fit i128. Ring coefficients do
// not have such a bound and are checked separately.
struct ClassicalOutputs {
    minimum: Vec<i128>,
    highest_prefix: Vec<i128>,
    least_tail: Vec<i128>,
}

impl ClassicalOutputs {
    fn new(a: &[i128], b: &[i128]) -> Self {
        let rank = a.len();
        let weights = a.iter().rev().zip(b.iter().rev());
        let minimum: Vec<_> = weights.clone().map(|(&x, &y)| x.max(y)).collect();
        let mut highest_prefix = vec![0];
        for (&x, &y) in weights {
            highest_prefix.push(highest_prefix.last().unwrap() + x + y);
        }
        let mut least_tail = vec![0; rank + 1];
        for index in (0..rank).rev() {
            least_tail[index] = least_tail[index + 1] + minimum[index];
        }
        Self {
            minimum,
            highest_prefix,
            least_tail,
        }
    }

    fn visit(
        &self,
        model: &mut Model,
        mut solution: impl FnMut(&mut Model, &[i128]) -> Result<()>,
    ) -> Result<()> {
        let total = *self.highest_prefix.last().unwrap();
        let mut weight = Vec::with_capacity(self.minimum.len());
        self.next(model, &mut weight, total, total, &mut solution)
    }

    fn next(
        &self,
        model: &mut Model,
        weight: &mut Vec<i128>,
        total: i128,
        remaining: i128,
        solution: &mut impl FnMut(&mut Model, &[i128]) -> Result<()>,
    ) -> Result<()> {
        model.spend(1)?;
        let index = weight.len();
        let slots = (self.minimum.len() - index) as i128;
        if slots == 0 {
            return solution(model, weight);
        }
        // Nonzero LR coefficients require componentwise containment of both
        // inputs and dominance by their entrywise sum, in decreasing order.
        let lower = self.minimum[index].max((remaining + slots - 1) / slots);
        let upper = weight
            .last()
            .copied()
            .unwrap_or(total)
            .min(remaining - self.least_tail[index + 1])
            .min(self.highest_prefix[index + 1] - (total - remaining));
        for entry in lower..=upper {
            weight.push(entry);
            self.next(model, weight, total, remaining - entry, solution)?;
            weight.pop();
        }
        Ok(())
    }
}

fn fold(model: &mut Model, grade: Grade, highest: &[i128]) -> Result<Option<(Weight, i128)>> {
    let n = grade.rank as i128;
    let period = n + grade.level as i128;
    let mut quotients = 0;
    let mut residues = Vec::with_capacity(grade.rank);
    model.spend(size(grade.rank) + grade.rank * grade.rank)?;
    for (index, &entry) in highest.iter().enumerate() {
        let shifted = entry + n - 1 - index as i128;
        quotients += shifted.div_euclid(period);
        residues.push((shifted.rem_euclid(period), index));
    }
    residues.sort_unstable();
    if residues.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Ok(None);
    }
    // Preserve the coordinate sum by raising the smallest h residues by one
    // period and all residues by q periods, where quotients = q*n + h.
    // The common q-period shift cancels after determinant normalization.
    let h = quotients.rem_euclid(n) as usize;
    for (value, _) in &mut residues[..h] {
        *value += period;
    }
    residues.rotate_left(h);
    let mut sign = 1;
    for (index, &(_, original)) in residues.iter().rev().enumerate() {
        for &(_, later) in residues.iter().rev().skip(index + 1) {
            if original > later {
                sign = -sign;
            }
        }
    }
    let center = residues[0].0;
    let weight = residues
        .into_iter()
        .enumerate()
        .map(|(index, (entry, _))| entry - center - index as i128)
        .collect();
    Ok(Some((weight, sign)))
}

fn basis_product(
    model: &mut Model,
    grade: Grade,
    left: &Weight,
    right: &Weight,
    kept: usize,
) -> Result<Sum> {
    let mut result = Sum::default();
    if left[grade.rank - 1] == 0 || right[grade.rank - 1] == 0 {
        result.add(model, left.max(right), 1, kept)?;
        return Ok(result);
    }
    ClassicalOutputs::new(left, right).visit(model, |model, highest| {
        let Some((weight, sign)) = fold(model, grade, highest)? else {
            return Ok(());
        };
        let highest: Weight = highest.iter().rev().copied().collect();
        let boundary = Boundary([left.clone(), right.clone(), highest]);
        let multiplicity =
            i128::try_from(rv_hive::count(model, &boundary)?).map_err(|_| arithmetic())?;
        result.add(model, &weight, sign * multiplicity, kept)
    })?;
    if result.0.values().any(|&coefficient| coefficient < 0) {
        return Err(Error::plain(
            "Kac–Walton fold produced a negative basis multiplicity",
        ));
    }
    Ok(result)
}

pub(crate) fn multiply_support(
    model: &mut Model,
    grade: Grade,
    a: &Support,
    b: &Support,
) -> Result<Support> {
    let mut output = Sum::default();
    for (left, c) in a.iter() {
        for (right, d) in b.iter() {
            model.spend(grade.rank)?;
            let factor = c.checked_mul(*d).ok_or_else(arithmetic)?;
            let product =
                basis_product(model, grade, left, right, size(grade.rank) * output.0.len())?;
            for (weight, multiplicity) in product.0 {
                let contribution = factor.checked_mul(multiplicity).ok_or_else(arithmetic)?;
                output.add(model, &weight, contribution, 0)?;
            }
        }
    }
    Ok(output.0.into_iter().collect())
}

pub(crate) fn multiply(model: &mut Model, grade: Grade, a: &Support, b: &Support) -> Result<Value> {
    let product = multiply_support(model, grade, a, b)?;
    model.allocate(size(grade.rank) * product.len())?;
    Ok(Value::RVFusion(product))
}

fn unit(grade: Grade) -> Support {
    [(vec![0; grade.rank].into(), 1)].into()
}

pub(crate) fn visit_weights(
    model: &mut Model,
    grade: Grade,
    limit: usize,
    limit_error: &'static str,
    mut visit: impl FnMut(&mut Model, &[i128]) -> Result<()>,
) -> Result<()> {
    let mut labels = 1i128;
    for i in 1..grade.rank {
        labels = labels * (grade.level as i128 + i as i128) / i as i128;
        if labels > limit as i128 {
            return Err(Error::plain(limit_error));
        }
    }
    model.spend(labels as usize * grade.rank)?;
    let mut coordinates = vec![0i128; grade.rank];
    let level = grade.level as i128;
    loop {
        visit(model, &coordinates)?;
        let Some(index) = (1..grade.rank).rev().find(|&i| coordinates[i] < level) else {
            break;
        };
        let next = coordinates[index] + 1;
        coordinates[index..].fill(next);
    }
    Ok(())
}

pub(crate) fn handle_support(model: &mut Model, grade: Grade) -> Result<Support> {
    let mut result = Sum::default();
    visit_weights(
        model,
        grade,
        4096,
        "fusion handle exceeds the limit of 4096 integrable weights",
        |model, coordinates| {
            let conjugate = dual_weight(coordinates);
            let weight: Weight = coordinates.into();
            let product = basis_product(
                model,
                grade,
                &weight,
                &conjugate,
                size(grade.rank) * result.0.len(),
            )?;
            for (weight, multiplicity) in product.0 {
                result.add(model, &weight, multiplicity, 0)?;
            }
            Ok(())
        },
    )?;
    Ok(result.0.into_iter().collect())
}

pub(crate) fn handle(model: &mut Model, grade: Grade) -> Result<Value> {
    let handle = handle_support(model, grade)?;
    model.allocate(size(grade.rank) * handle.len())?;
    Ok(Value::RVFusion(handle))
}

pub(crate) fn puncture_product(
    model: &mut Model,
    grade: Grade,
    punctures: Vec<Weight>,
) -> Result<Support> {
    let mut product = unit(grade);
    for weight in punctures {
        let puncture: Support = [(weight, 1)].into();
        product = multiply_support(model, grade, &product, &puncture)?;
    }
    Ok(product)
}

pub(crate) fn basis_product_support(
    model: &mut Model,
    grade: Grade,
    a: &Weight,
    b: &Weight,
) -> Result<Support> {
    Ok(basis_product(model, grade, a, b, 0)?
        .0
        .into_iter()
        .collect())
}

pub(crate) fn conformal_dimension(
    model: &mut Model,
    grade: Grade,
    mut genus: usize,
    punctures: Vec<Vec<i128>>,
) -> Result<i128> {
    let punctures = punctures
        .into_iter()
        .map(|weight| grade.normalize(model, weight))
        .collect::<Result<_>>()?;
    let mut product = puncture_product(model, grade, punctures)?;
    if genus > 0 {
        let mut handle = handle_support(model, grade)?;
        loop {
            if genus & 1 == 1 {
                product = multiply_support(model, grade, &product, &handle)?;
            }
            genus >>= 1;
            if genus == 0 {
                break;
            }
            handle = multiply_support(model, grade, &handle, &handle)?;
        }
    }
    Ok(product
        .first()
        .filter(|(weight, _)| weight[grade.rank - 1] == 0)
        .map_or(0, |(_, coefficient)| *coefficient))
}
