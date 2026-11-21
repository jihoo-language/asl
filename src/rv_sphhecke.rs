use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_iwhecke::{self, Accumulator, Monomial, Polynomial};

const MAX_CELLS: usize = 300_000;

#[derive(PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Term {
    pub(crate) coweight: Rc<[i128]>,
    pub(crate) coefficient: Polynomial,
}

pub(crate) type Support = Rc<[Term]>;

fn arithmetic() -> Error {
    Error::plain("spherical Hecke arithmetic is out of range")
}

fn storage() -> Error {
    Error::plain("spherical Hecke calculation exceeds the storage limit")
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(storage());
    }
    Ok(())
}

fn cells(rank: usize, coefficients: usize) -> Result<usize> {
    rank.checked_mul(2)
        .and_then(|size| {
            coefficients
                .checked_mul(3)
                .and_then(|extra| size.checked_add(extra))
        })
        .and_then(|size| size.checked_add(1))
        .ok_or_else(storage)
}

fn persist(model: &mut Model, support: Vec<Term>) -> Result<Value> {
    let mut total = 0;
    for term in &support {
        total += cells(term.coweight.len(), term.coefficient.len())?;
        check_storage(total)?;
    }
    model.allocate(total)?;
    Ok(Value::RVSphHecke(support.into()))
}

pub(crate) fn zero() -> Value {
    Value::RVSphHecke([].into())
}

pub(crate) fn validate(model: &mut Model, coweight: &[i128]) -> Result<()> {
    model.spend(coweight.len())?;
    if coweight.windows(2).any(|pair| pair[0] > pair[1]) {
        return Err(Error::plain(
            "the spherical coweight must be weakly increasing",
        ));
    }
    Ok(())
}

pub(crate) fn basis(model: &mut Model, coweight: Vec<i128>) -> Result<Value> {
    validate(model, &coweight)?;
    persist(
        model,
        vec![Term {
            coweight: coweight.into(),
            coefficient: [(0, 1)].into(),
        }],
    )
}

pub(crate) fn scalar(model: &mut Model, rank: usize, degree: usize) -> Result<Value> {
    check_storage(cells(rank, 1)?)?;
    model.spend(rank)?;
    persist(
        model,
        vec![Term {
            coweight: vec![0; rank].into(),
            coefficient: [(degree, 1)].into(),
        }],
    )
}

pub(crate) fn compare(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Ordering> {
    for (a, b) in left.iter().zip(right) {
        model.spend(
            1 + a.coweight.len().min(b.coweight.len())
                + a.coefficient.len().min(b.coefficient.len()),
        )?;
        let order = a.cmp(b);
        if order != Ordering::Equal {
            return Ok(order);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

pub(crate) fn add(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Value> {
    let mut sum = Accumulator::default();
    for term in left.iter().chain(right) {
        sum.add(model, &term.coweight, &term.coefficient, 0, 1)?;
    }
    persist(
        model,
        sum.finish()
            .into_iter()
            .map(|term| Term {
                coweight: term.window,
                coefficient: term.coefficient,
            })
            .collect(),
    )
}

pub(crate) fn neg(model: &mut Model, terms: &[Term]) -> Result<Value> {
    let mut result = Vec::with_capacity(terms.len());
    for term in terms {
        model.spend(1 + term.coefficient.len())?;
        let coefficient = term
            .coefficient
            .iter()
            .map(|&(power, coefficient)| {
                Ok((power, coefficient.checked_neg().ok_or_else(arithmetic)?))
            })
            .collect::<Result<Vec<_>>>()?;
        result.push(Term {
            coweight: term.coweight.clone(),
            coefficient: coefficient.into(),
        });
    }
    persist(model, result)
}

pub(crate) fn coefficient<'a>(
    model: &mut Model,
    terms: &'a [Term],
    coweight: &[i128],
) -> Result<&'a [Monomial]> {
    model.spend(1 + coweight.len() * (1 + terms.len().max(1).ilog2() as usize))?;
    let result = terms.binary_search_by(|term| term.coweight.as_ref().cmp(coweight));
    Ok(result.map_or(&[][..], |index| terms[index].coefficient.as_ref()))
}

fn factorial(rank: usize) -> Result<usize> {
    let mut value = 1usize;
    for i in 2..=rank {
        value = value.checked_mul(i).ok_or_else(storage)?;
        check_storage(value)?;
    }
    Ok(value)
}

fn orbit_size(coweight: &[i128], permutations: usize) -> Result<usize> {
    let mut distinct = permutations;
    for repeated in coweight.chunk_by(|a, b| a == b) {
        distinct /= factorial(repeated.len())?;
    }
    permutations.checked_mul(distinct).ok_or_else(storage)
}

fn next_permutation(values: &mut [i128]) -> bool {
    let Some(pivot) = values.windows(2).rposition(|pair| pair[0] < pair[1]) else {
        return false;
    };
    let successor = values
        .iter()
        .rposition(|value| *value > values[pivot])
        .unwrap();
    values.swap(pivot, successor);
    values[pivot + 1..].reverse();
    true
}

fn expand(model: &mut Model, rank: usize, support: &[Term]) -> Result<Vec<rv_iwhecke::Term>> {
    if support.is_empty() {
        return Ok(Vec::new());
    }
    let permutations = factorial(rank)?;
    let mut expanded = Vec::new();
    let mut workspace = 0usize;
    for term in support {
        let count = orbit_size(&term.coweight, permutations)?;
        let required = cells(rank, term.coefficient.len())?
            .checked_mul(count)
            .ok_or_else(storage)?;
        workspace = workspace.checked_add(required).ok_or_else(storage)?;
        check_storage(workspace)?;
        let mut exponents = term.coweight.to_vec();
        let mut residues = Vec::with_capacity(rank);
        loop {
            residues.clear();
            residues.extend((0..rank).map(|row| row as i128));
            let base = exponents
                .iter()
                .map(|&exponent| exponent.checked_mul(rank as i128).ok_or_else(arithmetic))
                .collect::<Result<Vec<_>>>()?;
            loop {
                model.spend(1 + rank)?;
                let window = residues
                    .iter()
                    .zip(&base)
                    .map(|(&row, &shift)| row.checked_add(shift).ok_or_else(arithmetic))
                    .collect::<Result<Vec<_>>>()?;
                expanded.push(rv_iwhecke::Term {
                    window: window.into(),
                    coefficient: term.coefficient.clone(),
                });
                if !next_permutation(&mut residues) {
                    break;
                }
            }
            if !next_permutation(&mut exponents) {
                break;
            }
        }
    }
    model.spend(expanded.len() * (1 + rank) * (1 + expanded.len().max(1).ilog2() as usize))?;
    expanded.sort_unstable_by(|a, b| a.window.cmp(&b.window));
    Ok(expanded)
}

pub(crate) fn embed(model: &mut Model, rank: usize, support: &[Term]) -> Result<Value> {
    let expanded = expand(model, rank, support)?;
    rv_iwhecke::persist(model, expanded)
}

fn poincare(model: &mut Model, rank: usize) -> Result<Vec<i128>> {
    let mut polynomial = vec![1i128];
    for factor in 2..=rank {
        let mut product = vec![0i128; polynomial.len() + factor - 1];
        model.spend(polynomial.len() * factor)?;
        for (power, coefficient) in polynomial.into_iter().enumerate() {
            for entry in &mut product[power..power + factor] {
                *entry = entry.checked_add(coefficient).ok_or_else(arithmetic)?;
            }
        }
        polynomial = product;
    }
    Ok(polynomial)
}

fn divide(model: &mut Model, polynomial: &[Monomial], divisor: &[i128]) -> Result<Polynomial> {
    let degree = divisor.len() - 1;
    let mut remainder: BTreeMap<usize, i128> = polynomial.iter().copied().collect();
    let mut quotient = Vec::new();
    while let Some((&leading, &coefficient)) = remainder.last_key_value() {
        if leading < degree {
            return Err(Error::plain(
                "spherical Hecke normalization is not an exact polynomial quotient",
            ));
        }
        let shift = leading - degree;
        model.spend(1 + divisor.len() * (1 + remainder.len().ilog2() as usize))?;
        for (degree, &value) in divisor.iter().enumerate() {
            let power = shift + degree;
            let value = coefficient.checked_mul(value).ok_or_else(arithmetic)?;
            let updated = remainder
                .get(&power)
                .copied()
                .unwrap_or(0)
                .checked_sub(value)
                .ok_or_else(arithmetic)?;
            if updated == 0 {
                remainder.remove(&power);
            } else {
                remainder.insert(power, updated);
            }
        }
        check_storage(3 * (1 + quotient.len() + remainder.len()))?;
        quotient.push((shift, coefficient));
    }
    quotient.reverse();
    Ok(quotient.into())
}

pub(crate) fn multiply(
    model: &mut Model,
    rank: usize,
    left: &[Term],
    right: &[Term],
) -> Result<Value> {
    if left.is_empty() || right.is_empty() {
        return Ok(zero());
    }
    let left = expand(model, rank, left)?;
    let right = expand(model, rank, right)?;
    let product = rv_iwhecke::product(model, &left, &right)?;
    let permutations = factorial(rank)?;
    let divisor = poincare(model, rank)?;
    let mut orbits = BTreeMap::<Vec<i128>, (Polynomial, usize)>::new();
    for term in product {
        model.spend(1 + rank * (1 + rank.max(1).ilog2() as usize) + term.coefficient.len())?;
        let mut coweight: Vec<_> = term
            .window
            .iter()
            .map(|&value| value.div_euclid(rank as i128))
            .collect();
        coweight.sort_unstable();
        let (coefficient, found) = orbits
            .entry(coweight)
            .or_insert((term.coefficient.clone(), 0));
        if *coefficient != term.coefficient {
            return Err(Error::plain(
                "spherical Hecke product has unequal coefficients on a double coset",
            ));
        }
        *found += 1;
    }
    let mut normalized = Vec::with_capacity(orbits.len());
    for (coweight, (coefficient, found)) in orbits {
        if found != orbit_size(&coweight, permutations)? {
            return Err(Error::plain(
                "spherical Hecke product has an incomplete double coset",
            ));
        }
        let coefficient = divide(model, &coefficient, &divisor)?;
        normalized.push(Term {
            coweight: coweight.into(),
            coefficient,
        });
    }
    persist(model, normalized)
}
