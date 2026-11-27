use std::collections::BTreeMap;

use crate::error::Result;
use crate::model::Model;

use super::{arithmetic, guard_cells, Denominator, Integer, LIMBS};

pub(super) type Vector = Vec<(usize, Integer)>;

pub(super) struct Column {
    pub(super) image: Vector,
    pub(super) witness: Vector,
}

impl Column {
    fn cells(&self) -> usize {
        8 + (LIMBS + 3) * (self.image.len() + self.witness.len())
    }

    fn primitive(&mut self, model: &mut Model) -> Result<()> {
        let one = Denominator::from_small(1);
        let mut common = Denominator::default();
        for (_, entry) in self.image.iter().chain(&self.witness) {
            common = common.gcd(model, entry.magnitude)?;
            if common == one {
                return Ok(());
            }
        }
        if !common.is_zero() {
            model.spend(self.image.len() + self.witness.len())?;
            for (_, entry) in self.image.iter_mut().chain(&mut self.witness) {
                *entry = divide(model, *entry, common)?;
            }
        }
        Ok(())
    }
}

#[derive(Default)]
pub(super) struct Echelon {
    pivots: BTreeMap<usize, Column>,
    cells: usize,
}

impl Echelon {
    pub(super) fn insert(&mut self, model: &mut Model, column: Column) -> Result<bool> {
        let column = self.reduce(model, column)?;
        let Some(&(pivot, _)) = column.image.last() else {
            return Ok(false);
        };
        self.cells += column.cells();
        guard_cells(self.cells)?;
        self.pivots.insert(pivot, column);
        Ok(true)
    }

    pub(super) fn solve(&self, model: &mut Model, column: Column) -> Result<Column> {
        self.reduce(model, column)
    }

    fn reduce(&self, model: &mut Model, mut column: Column) -> Result<Column> {
        guard_cells(self.cells + column.cells())?;
        column.primitive(model)?;
        while let Some(&(coordinate, right)) = column.image.last() {
            let Some(known) = self.pivots.get(&coordinate) else {
                break;
            };
            guard_cells(self.cells + 2 * column.cells() + known.cells())?;
            let left = known.image.last().expect("nonzero pivot column").1;
            let common = left.magnitude.gcd(model, right.magnitude)?;
            let a = divide(model, left, common)?;
            let b = divide(model, right, common)?.negate();
            column.image = combine(model, &column.image, a, &known.image, b)?;
            column.witness = combine(model, &column.witness, a, &known.witness, b)?;
            column.primitive(model)?;
        }
        Ok(column)
    }
}

pub(super) fn divide(model: &mut Model, value: Integer, divisor: Denominator) -> Result<Integer> {
    Ok(Integer {
        magnitude: super::super::divide(model, value.magnitude, divisor)?,
        negative: value.negative,
    })
}

fn combine(
    model: &mut Model,
    left: &Vector,
    a: Integer,
    right: &Vector,
    b: Integer,
) -> Result<Vector> {
    model.spend(1 + (left.len() + right.len()) * LIMBS)?;
    let mut result = Vec::with_capacity(left.len() + right.len());
    let mut left = left.iter().peekable();
    let mut right = right.iter().peekable();
    while left.peek().is_some() || right.peek().is_some() {
        let coordinate = match (left.peek(), right.peek()) {
            (Some(x), Some(y)) => x.0.min(y.0),
            (Some(x), None) | (None, Some(x)) => x.0,
            (None, None) => unreachable!(),
        };
        let x = left.next_if(|entry| entry.0 == coordinate);
        let y = right.next_if(|entry| entry.0 == coordinate);
        let weighted = |entry: Option<&(usize, Integer)>, factor| match entry {
            Some((_, value)) => value.multiply_wide(factor).map_err(|_| arithmetic()),
            None => Ok(Integer::from_i128(0)),
        };
        let coefficient = weighted(x, a)?
            .add(weighted(y, b)?)
            .map_err(|_| arithmetic())?;
        if !coefficient.is_zero() {
            result.push((coordinate, coefficient));
        }
    }
    Ok(result)
}
