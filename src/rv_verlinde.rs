use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_cyclotomic;
use crate::rv_fusion::{Grade, Support};
use crate::rv_satake;

pub(crate) fn order(grade: Grade) -> Result<usize> {
    let order = grade
        .level
        .checked_add(grade.rank)
        .and_then(|value| value.checked_mul(grade.rank))
        .ok_or_else(|| Error::plain("Verlinde torus order is out of range"))?;
    rv_cyclotomic::check_order(order)?;
    Ok(order)
}

pub(crate) fn torus(model: &mut Model, grade: Grade, label: Vec<i128>) -> Result<Vec<i128>> {
    let label = grade.normalize(model, label)?;
    // The determinant shift produces the SU(n) torus with product one.
    let shifted: Vec<_> = label
        .iter()
        .rev()
        .enumerate()
        .map(|(i, &entry)| entry + (grade.rank - 1 - i) as i128)
        .collect();
    let total: i128 = shifted.iter().sum();
    Ok(shifted
        .into_iter()
        .map(|entry| grade.rank as i128 * entry - total)
        .collect())
}

pub(crate) fn character(
    model: &mut Model,
    order: usize,
    torus: &[i128],
    terms: &Support,
) -> Result<Value> {
    let mut coefficients = vec![0i128; order];
    let arithmetic = || Error::plain("Verlinde character integer arithmetic is out of range");
    for (weight, scalar) in terms.iter() {
        for term in rv_satake::schur_terms(model, weight, false)? {
            model.spend(term.exponents.len())?;
            let power: i128 = term.exponents.iter().zip(torus).map(|(&a, &b)| a * b).sum();
            let entry = &mut coefficients[power.rem_euclid(order as i128) as usize];
            *entry = scalar
                .checked_mul(term.coefficient)
                .and_then(|value| entry.checked_add(value))
                .ok_or_else(arithmetic)?;
        }
    }
    rv_cyclotomic::finish(model, order, coefficients)
}
