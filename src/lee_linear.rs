use std::collections::{btree_map::Entry, BTreeMap};

use crate::error::{Error, Result};
use crate::integer_smith::storage;
use crate::model::Model;
use crate::rv_lens::gcd;

pub(crate) type Vector = Vec<(usize, i128)>;

pub(crate) fn arithmetic() -> Error {
    Error::plain("Lee exact rational arithmetic exceeds the signed 128-bit coefficient range")
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Rational {
    pub(crate) numerator: i128,
    pub(crate) denominator: i128,
}

impl Rational {
    pub(crate) const ZERO: Self = Self::integer(0);

    pub(crate) const fn integer(numerator: i128) -> Self {
        Self {
            numerator,
            denominator: 1,
        }
    }

    pub(crate) fn new(numerator: i128, denominator: i128) -> Result<Self> {
        if denominator <= 0 {
            return Err(Error::plain("Lee rational denominator must be positive"));
        }
        Self::ratio(numerator, denominator)
    }

    pub(crate) fn ratio(numerator: i128, denominator: i128) -> Result<Self> {
        if denominator == 0 {
            return Err(Error::plain("Lee rational denominator cannot vanish"));
        }
        let divisor = gcd(numerator.unsigned_abs(), denominator.unsigned_abs());
        let mut numerator = divide(numerator, divisor);
        let mut denominator = divide(denominator, divisor);
        if denominator < 0 {
            numerator = numerator.checked_neg().ok_or_else(arithmetic)?;
            denominator = denominator.checked_neg().ok_or_else(arithmetic)?;
        }
        Ok(Self {
            numerator,
            denominator,
        })
    }

    pub(crate) fn add(self, other: Self) -> Result<Self> {
        let common = gcd(self.denominator as u128, other.denominator as u128) as i128;
        let a = self.denominator / common;
        let b = other.denominator / common;
        let numerator = checked_multiply(self.numerator, b)?
            .checked_add(checked_multiply(other.numerator, a)?)
            .ok_or_else(arithmetic)?;
        let reduction = gcd(numerator.unsigned_abs(), common as u128) as i128;
        let denominator = checked_multiply(a, other.denominator / reduction)?;
        Self::ratio(numerator / reduction, denominator)
    }

    pub(crate) fn multiply(self, other: Self) -> Result<Self> {
        let a = gcd(self.numerator.unsigned_abs(), other.denominator as u128) as i128;
        let b = gcd(other.numerator.unsigned_abs(), self.denominator as u128) as i128;
        let numerator = checked_multiply(self.numerator / a, other.numerator / b)?;
        let denominator = checked_multiply(self.denominator / b, other.denominator / a)?;
        Ok(Self {
            numerator,
            denominator,
        })
    }
}

pub(crate) fn checked_multiply(a: i128, b: i128) -> Result<i128> {
    a.checked_mul(b).ok_or_else(arithmetic)
}

fn divide(value: i128, divisor: u128) -> i128 {
    if let Ok(divisor) = i128::try_from(divisor) {
        return value / divisor;
    }
    // The only possible nonzero multiple of 2^127 in an i128 is MIN.
    if value == 0 {
        0
    } else {
        -1
    }
}

struct Row {
    image: Vector,
    witness: Vector,
}

impl Row {
    fn cells(&self) -> usize {
        8 + 4 * (self.image.len() + self.witness.len())
    }

    fn primitive(&mut self, model: &mut Model) -> Result<()> {
        model.spend(1 + self.image.len() + self.witness.len())?;
        let divisor = self
            .image
            .iter()
            .chain(&self.witness)
            .fold(0, |d, (_, v)| gcd(d, v.unsigned_abs()));
        if divisor > 1 {
            for (_, v) in self.image.iter_mut().chain(&mut self.witness) {
                *v = divide(*v, divisor);
            }
        }
        Ok(())
    }
}

#[derive(Default)]
pub(crate) struct Echelon {
    rows: BTreeMap<usize, Row>,
    cells: usize,
}

impl Echelon {
    pub(crate) fn rank(&self) -> usize {
        self.rows.len()
    }

    // An inserted image is independent unless an exact dependence is
    // returned. Witnesses are transformed by the very same row operations.
    pub(crate) fn insert(
        &mut self,
        model: &mut Model,
        image: Vector,
        witness: Vector,
    ) -> Result<Option<Vector>> {
        let row = self.reduce(model, Row { image, witness })?;
        let Some(&(pivot, _)) = row.image.last() else {
            return Ok(Some(row.witness));
        };
        self.cells += row.cells();
        storage(self.cells)?;
        let Entry::Vacant(entry) = self.rows.entry(pivot) else {
            unreachable!()
        };
        entry.insert(row);
        Ok(None)
    }

    pub(crate) fn spans(&self, model: &mut Model, image: &Vector) -> Result<bool> {
        // Do not track coordinates on a membership-only query.
        let row = self.reduce(
            model,
            Row {
                image: image.clone(),
                witness: Vec::new(),
            },
        )?;
        Ok(row.image.is_empty())
    }

    pub(crate) fn relation(
        &self,
        model: &mut Model,
        image: Vector,
        witness: Vector,
    ) -> Result<Option<Vector>> {
        let row = self.reduce(model, Row { image, witness })?;
        Ok(row.image.is_empty().then_some(row.witness))
    }

    pub(crate) fn reduced(
        &self,
        model: &mut Model,
        image: Vector,
        witness: Vector,
    ) -> Result<(Vector, Vector)> {
        let row = self.reduce(model, Row { image, witness })?;
        Ok((row.image, row.witness))
    }

    pub(crate) fn pivot(&self, index: usize) -> Option<&Vector> {
        self.rows.get(&index).map(|row| &row.image)
    }

    // The greatest coordinate is the pivot. Stopping at a new pivot is
    // sufficient for exact span tests as all earlier rows are triangular.
    fn reduce(&self, model: &mut Model, mut row: Row) -> Result<Row> {
        storage(self.cells + row.cells())?;
        let track = !row.witness.is_empty();
        row.primitive(model)?;
        while let Some(&(pivot, coefficient)) = row.image.last() {
            model.spend(1 + self.rows.len().max(1).ilog2() as usize)?;
            let Some(known) = self.rows.get(&pivot) else {
                break;
            };
            let leading = known.image.last().expect("nonzero pivot row").1;
            let divisor = gcd(leading.unsigned_abs(), coefficient.unsigned_abs());
            let p = divide(leading, divisor);
            let q = divide(coefficient, divisor);
            storage(self.cells + 2 * row.cells() + known.cells())?;
            row.image = combine(model, &row.image, p, &known.image, q)?;
            if track {
                row.witness = combine(model, &row.witness, p, &known.witness, q)?;
            }
            row.primitive(model)?;
        }
        Ok(row)
    }
}

// Signed exact arithmetic over Q by primitive integer representatives.
fn combine(model: &mut Model, left: &Vector, a: i128, right: &Vector, b: i128) -> Result<Vector> {
    model.spend(1 + left.len() + right.len())?;
    let mut result = Vec::with_capacity(left.len() + right.len());
    let mut left = left.iter().peekable();
    let mut right = right.iter().peekable();
    while left.peek().is_some() || right.peek().is_some() {
        let key = match (left.peek(), right.peek()) {
            (Some(x), Some(y)) => x.0.min(y.0),
            (Some(x), None) => x.0,
            (None, Some(y)) => y.0,
            _ => unreachable!(),
        };
        let x = left.next_if(|x| x.0 == key).map_or(0, |x| x.1);
        let y = right.next_if(|y| y.0 == key).map_or(0, |y| y.1);
        let c = checked_multiply(a, x)?
            .checked_sub(checked_multiply(b, y)?)
            .ok_or_else(arithmetic)?;
        if c != 0 {
            result.push((key, c));
        }
    }
    Ok(result)
}

impl std::fmt::Display for Rational {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.denominator == 1 {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "{}/{}", self.numerator, self.denominator)
        }
    }
}
