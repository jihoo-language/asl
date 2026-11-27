use std::cmp::Ordering;
use std::collections::{btree_map::Entry, BTreeMap};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};

const MAX_TERMS: usize = 60_000;
pub(crate) type Support = Rc<[(Powers, i128)]>;

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Powers {
    pub(crate) q: usize,
    pub(crate) theta: usize,
}

pub(crate) fn arithmetic() -> Error {
    Error::plain("quantum differential integer arithmetic is out of range")
}

pub(crate) fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

#[derive(Default)]
pub(crate) struct Sum(BTreeMap<Powers, i128>);

impl Sum {
    pub(crate) fn add(&mut self, model: &mut Model, powers: Powers, value: i128) -> Result<()> {
        let width = self.0.len();
        model.spend(1 + width.max(1).ilog2() as usize)?;
        if value == 0 {
            return Ok(());
        }
        match self.0.entry(powers) {
            Entry::Occupied(mut slot) => {
                let coefficient = value.checked_add(*slot.get()).ok_or_else(arithmetic)?;
                if coefficient == 0 {
                    slot.remove();
                } else {
                    *slot.get_mut() = coefficient;
                }
            }
            Entry::Vacant(slot) => {
                if width == MAX_TERMS {
                    return Err(Error::plain("Ore operator exceeds the term storage limit"));
                }
                slot.insert(value);
            }
        }
        Ok(())
    }

    pub(crate) fn finish(self, model: &mut Model) -> Result<Value> {
        model.allocate(4 * self.0.len())?;
        Ok(Value::RVOre(self.0.into_iter().collect()))
    }
}

pub(crate) fn zero() -> Value {
    Value::RVOre([].into())
}

pub(crate) fn term(model: &mut Model, q: usize, theta: usize, coefficient: i128) -> Result<Value> {
    let mut operator = Sum::default();
    operator.add(model, Powers { q, theta }, coefficient)?;
    operator.finish(model)
}

pub(crate) fn add(model: &mut Model, left: &Support, right: &Support) -> Result<Value> {
    let mut result = Sum::default();
    for &(powers, coefficient) in left.iter().chain(right.iter()) {
        result.add(model, powers, coefficient)?;
    }
    result.finish(model)
}

pub(crate) fn negate(model: &mut Model, terms: &Support) -> Result<Value> {
    let mut result = Sum::default();
    for &(powers, coefficient) in terms.iter() {
        result.add(
            model,
            powers,
            coefficient.checked_neg().ok_or_else(arithmetic)?,
        )?;
    }
    result.finish(model)
}

pub(crate) fn compare(model: &mut Model, left: &Support, right: &Support) -> Result<Ordering> {
    if !Rc::ptr_eq(left, right) {
        model.spend(left.len().min(right.len()))?;
    }
    Ok(left.cmp(right))
}

pub(crate) fn coefficient(
    model: &mut Model,
    terms: &Support,
    q: usize,
    theta: usize,
) -> Result<i128> {
    model.spend(1 + terms.len().max(1).ilog2() as usize)?;
    Ok(terms
        .binary_search_by_key(&Powers { q, theta }, |&(powers, _)| powers)
        .map(|index| terms[index].1)
        .unwrap_or(0))
}

pub(crate) fn multiply(model: &mut Model, left: &Support, right: &Support) -> Result<Value> {
    let mut result = Sum::default();
    for &(a, c) in left.iter() {
        for &(b, d) in right.iter() {
            model.spend(1)?;
            let q = a.q.checked_add(b.q).ok_or_else(arithmetic)?;
            let highest = a.theta.checked_add(b.theta).ok_or_else(arithmetic)?;
            let scalar = c.checked_mul(d).ok_or_else(arithmetic)?;
            if b.q == 0 || a.theta == 0 {
                result.add(model, Powers { q, theta: highest }, scalar)?;
                continue;
            }
            model.spend(a.theta)?;
            let mut binomial = 1i128;
            let mut shift = scalar;
            for dropped in 0..=a.theta {
                let theta = highest - dropped;
                let value = shift.checked_mul(binomial).ok_or_else(arithmetic)?;
                result.add(model, Powers { q, theta }, value)?;
                if dropped != a.theta {
                    // Update C(j,r) by cancelling the denominator before multiplying.
                    let numerator = a.theta - dropped;
                    let denominator = dropped + 1;
                    let common = gcd(numerator as u128, denominator as u128) as usize;
                    binomial = (binomial / (denominator / common) as i128)
                        .checked_mul((numerator / common) as i128)
                        .ok_or_else(arithmetic)?;
                    shift = shift.checked_mul(b.q as i128).ok_or_else(arithmetic)?;
                }
            }
        }
    }
    result.finish(model)
}
