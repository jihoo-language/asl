use std::cmp::Ordering;
use std::rc::Rc;

use super::{Integer, K3Vector};
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod og6;

const RANK: usize = 8;
const MIDDLE: usize = 6;
const WEDGE_SIGNS: [i128; MIDDLE] = [-1, 1, -1, -1, 1, -1];

pub(crate) enum Expression {
    Literal(Label, i128, i128, [i128; MIDDLE]),
    Basis(Label, usize),
    Kummer(Box<ParsedExpr>),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Vector {
    // Rank, six wedge coordinates, and last Mukai coordinate.
    coordinates: [i128; RANK],
}

impl Vector {
    pub(crate) fn new(r: i128, wedge: [i128; MIDDLE], s: i128) -> Self {
        let mut result = Self::zero();
        result.coordinates[0] = r;
        result.coordinates[1..=MIDDLE].copy_from_slice(&wedge);
        result.coordinates[RANK - 1] = s;
        result
    }

    pub(crate) fn zero() -> Self {
        Self {
            coordinates: [0; RANK],
        }
    }

    pub(crate) fn basis(index: usize) -> Result<Self> {
        let mut result = Self::zero();
        let coordinate = result
            .coordinates
            .get_mut(index)
            .ok_or_else(|| Error::plain("rvAbMukaiBasis index must be between 0 and 7"))?;
        *coordinate = 1;
        Ok(result)
    }

    pub(crate) fn coordinates(&self) -> (i128, &[i128], i128) {
        (
            self.coordinates[0],
            &self.coordinates[1..=MIDDLE],
            self.coordinates[RANK - 1],
        )
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(2 * RANK + 1)?;
        Ok(Value::AbMukai(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        model.spend(RANK)?;
        Ok(self.cmp(other))
    }

    pub(crate) fn add(&self, model: &mut Model, other: &Self) -> Result<Self> {
        model.spend(RANK)?;
        let mut result = Self::zero();
        for (out, (&x, &y)) in result
            .coordinates
            .iter_mut()
            .zip(self.coordinates.iter().zip(&other.coordinates))
        {
            *out = x.checked_add(y).ok_or_else(coordinate_range)?;
        }
        Ok(result)
    }

    pub(crate) fn negate(&self, model: &mut Model) -> Result<Self> {
        model.spend(RANK)?;
        let mut result = Self::zero();
        for (out, &x) in result.coordinates.iter_mut().zip(&self.coordinates) {
            *out = x.checked_neg().ok_or_else(coordinate_range)?;
        }
        Ok(result)
    }

    fn wide_pair(&self, model: &mut Model, other: &Self) -> Result<Integer> {
        model.spend(30 * RANK)?;
        let (r, wedge, s) = self.coordinates();
        let (rr, other_wedge, ss) = other.coordinates();
        let mut pairing = Integer::product(r, ss)
            .negate()
            .add(Integer::product(s, rr).negate())?;
        for ((&x, &y), &sign) in wedge.iter().zip(other_wedge.iter().rev()).zip(&WEDGE_SIGNS) {
            let term = Integer::product(x, y);
            pairing = pairing.add(if sign == 1 { term } else { term.negate() })?;
        }
        Ok(pairing)
    }

    pub(crate) fn pair(&self, model: &mut Model, other: &Self) -> Result<i128> {
        self.wide_pair(model, other)?
            .signed()
            .ok_or_else(|| Error::plain("abelian Mukai pairing exceeds the signed 128-bit range"))
    }

    pub(crate) fn kummer_middle(&self, model: &mut Model) -> Result<K3Vector> {
        model.spend(MIDDLE)?;
        let mut wedge = [0; MIDDLE];
        for (out, &x) in wedge.iter_mut().zip(self.coordinates().1) {
            *out = x.checked_mul(2).ok_or_else(|| {
                Error::plain("Kummer wedge coordinate exceeds the signed 128-bit range")
            })?;
        }
        K3Vector::new(wedge, [0; 16])
    }
}

fn coordinate_range() -> Error {
    Error::plain("abelian Mukai coordinate exceeds the signed 128-bit range")
}
