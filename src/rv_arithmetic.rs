use std::cmp::Ordering;
use std::collections::{btree_map::Entry, BTreeMap};

use crate::binary_field::BinaryField;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId};
use crate::rank_variety::{Endpoints, ProjectiveSpace};

#[derive(Clone, Copy)]
pub(crate) enum Series {
    Closed,
    Zeta,
}

impl Series {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Closed => "rvclosed",
            Self::Zeta => "rvzeta",
        }
    }
}

pub(crate) fn orbit(model: &mut Model, modulus: usize, point: &[usize]) -> Result<Vec<Vec<usize>>> {
    let field = BinaryField::new(model, modulus)?;
    model.spend(1 + point.len())?;
    field.check_point(point, point.len(), "rvorbit")?;
    let pivot = point
        .iter()
        .copied()
        .find(|&coefficient| coefficient != 0)
        .ok_or_else(|| Error::plain("the zero vector does not define a projective point"))?;
    let mut current = point.to_vec();
    if pivot != 1 {
        let inverse = field.inverse(model, pivot)?;
        model.spend(field.degree() * current.len())?;
        for coefficient in &mut current {
            *coefficient = field.multiply(*coefficient, inverse);
        }
    }

    let mut result = Vec::with_capacity(field.degree());
    for _ in 0..field.degree() {
        model.allocate(1 + current.len())?;
        result.push(current.clone());
        square_point(model, field, &mut current)?;
        if current == result[0] {
            // Preserve the Frobenius order and choose the lexicographic start.
            model.spend(result.len() * current.len())?;
            let start = result
                .iter()
                .enumerate()
                .min_by_key(|(_, point)| *point)
                .unwrap()
                .0;
            result.rotate_left(start);
            return Ok(result);
        }
    }
    Err(Error::plain("Frobenius orbit did not close in its field"))
}

type Counts = BTreeMap<usize, Vec<u128>>;

pub(crate) fn evaluate(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    selected: &[usize],
    max_degree: usize,
    series: Series,
) -> Result<Counts> {
    model.check_q_endpoints(source, target, selected)?;
    if selected.is_empty() || max_degree == 0 {
        return Ok(BTreeMap::new());
    }
    if max_degree >= usize::BITS as usize {
        return Err(Error::plain("binary extension degree is out of range"));
    }

    let endpoints = Endpoints::build(model, source, target, selected)?;
    let mut scratch = vec![0; selected.len()];
    let mut counts = BTreeMap::<usize, Vec<u128>>::new();
    for degree in 1..=max_degree {
        let field = BinaryField::of_degree(model, degree)?;
        let space = ProjectiveSpace::new(field, selected.len(), None)?;
        space.visit(model, |model, point| {
            if exact_degree_representative(model, field, point, &mut scratch)? {
                let dimension = endpoints.dimension(model, field, point)?;
                let row = match counts.entry(dimension) {
                    Entry::Occupied(entry) => entry.into_mut(),
                    Entry::Vacant(entry) => {
                        model.allocate(1 + 2 * max_degree)?;
                        entry.insert(vec![0; max_degree])
                    }
                };
                row[degree - 1] += 1;
            }
            Ok(())
        })?;
    }
    match series {
        Series::Closed => Ok(counts),
        Series::Zeta => counts
            .into_iter()
            .map(|(dimension, closed)| zeta(model, &closed).map(|row| (dimension, row)))
            .collect(),
    }
}

fn exact_degree_representative(
    model: &mut Model,
    field: BinaryField,
    point: &[usize],
    current: &mut [usize],
) -> Result<bool> {
    current.copy_from_slice(point);
    for period in 1..=field.degree() {
        square_point(model, field, current)?;
        model.spend(point.len())?;
        match Ord::cmp(current, point) {
            Ordering::Less => return Ok(false),
            Ordering::Equal => return Ok(period == field.degree()),
            Ordering::Greater => {}
        }
    }
    Err(Error::plain("Frobenius orbit did not close in its field"))
}

fn square_point(model: &mut Model, field: BinaryField, point: &mut [usize]) -> Result<()> {
    model.spend(1 + field.degree() * point.len())?;
    for coefficient in point {
        *coefficient = field.multiply(*coefficient, *coefficient);
    }
    Ok(())
}

fn zeta(model: &mut Model, closed: &[u128]) -> Result<Vec<u128>> {
    let max_degree = closed.len();
    let mut rational = vec![0u128; max_degree + 1];
    for (index, &count) in closed.iter().enumerate() {
        let degree = index + 1;
        let points = count
            .checked_mul(degree as u128)
            .ok_or_else(coefficient_overflow)?;
        // A degree-d orbit contributes all d geometric points over F_(2^r)
        // precisely when d divides r.
        for extension in (degree..=max_degree).step_by(degree) {
            model.spend(1)?;
            rational[extension] = rational[extension]
                .checked_add(points)
                .ok_or_else(coefficient_overflow)?;
        }
    }

    let mut coefficients = Vec::with_capacity(max_degree + 1);
    coefficients.push(1u128);
    for n in 1..=max_degree {
        model.spend(n)?;
        let mut numerator = 0u128;
        for r in 1..=n {
            let contribution = rational[r]
                .checked_mul(coefficients[n - r])
                .ok_or_else(coefficient_overflow)?;
            numerator = numerator
                .checked_add(contribution)
                .ok_or_else(coefficient_overflow)?;
        }
        if !numerator.is_multiple_of(n as u128) {
            return Err(Error::plain(
                "Hasse–Weil recurrence did not give an integer coefficient",
            ));
        }
        coefficients.push(numerator / n as u128);
    }
    model.allocate(2 * coefficients.len())?;
    Ok(coefficients)
}

fn coefficient_overflow() -> Error {
    Error::plain("Hasse–Weil coefficient is out of range")
}
