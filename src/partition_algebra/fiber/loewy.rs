use crate::error::{Error, Result};
use crate::lee_linear::{Echelon, Vector};
use crate::model::Model;
use crate::rv_lens::gcd;

use super::elimination::Kernel;
use super::{arithmetic, checked_mul, Coordinates, Fiber, FULL_LOEWY_RANK};

pub(super) struct Powers {
    // Entry zero describes J^2. Include the first zero layer, if present.
    layers: Vec<Layer>,
}

#[derive(Default)]
struct Layer {
    basis: Vec<Coordinates>,
    factors: Vec<(usize, usize)>,
    span: Echelon,
}

impl Powers {
    pub(super) fn new(model: &mut Model, fiber: &Fiber, radical: &Kernel) -> Result<Self> {
        let mut layers: Vec<Layer> = Vec::new();
        if radical.rows.is_empty() {
            return Ok(Self { layers });
        }
        if fiber.parameter.rank > FULL_LOEWY_RANK {
            return Err(Error::plain(
                "higher Loewy powers of a nonzero radical are limited to rank 2",
            ));
        }
        let first: Vec<_> = radical.rows.values().cloned().collect();
        let mut previous = &first;
        loop {
            let mut next = Layer::default();
            for (i, a) in previous.iter().enumerate() {
                for (j, b) in first.iter().enumerate() {
                    let product = fiber.multiply(model, a, b)?;
                    let integral = integral_row(&product)?;
                    if next.span.insert(model, integral, Vec::new())?.is_none() {
                        model.allocate(5 + 6 * product.len())?;
                        next.basis.push(product);
                        next.factors.push((i, j));
                    }
                }
            }
            let finished = next.basis.is_empty();
            if next.basis.len() >= previous.len() {
                return Err(Error::plain("partition radical power did not decrease"));
            }
            layers.push(next);
            if finished {
                return Ok(Self { layers });
            }
            previous = &layers.last().unwrap().basis;
        }
    }

    fn layer(&self, power: u128) -> Option<&Layer> {
        power
            .checked_sub(2)
            .and_then(|index| usize::try_from(index).ok())
            .and_then(|i| self.layers.get(i))
    }

    pub(super) fn dimensions(&self, dimension: usize, radical: usize) -> Vec<usize> {
        [dimension, radical]
            .into_iter()
            .chain(self.layers.iter().map(|layer| layer.basis.len()))
            .collect()
    }

    pub(super) fn basis(&self, power: u128) -> &[Coordinates] {
        self.layer(power)
            .map_or(&[], |layer| layer.basis.as_slice())
    }

    pub(super) fn contains(
        &self,
        model: &mut Model,
        power: u128,
        value: &Coordinates,
    ) -> Result<bool> {
        let Some(layer) = self.layer(power) else {
            return Ok(value.is_empty());
        };
        layer.span.spans(model, &integral_row(value)?)
    }

    pub(super) fn witness(&self, power: u128, index: u128) -> Result<(usize, usize)> {
        if power < 2 {
            return Err(Error::plain(
                "partition Loewy product witnesses require power at least 2",
            ));
        }
        let index = usize::try_from(index).ok();
        self.layer(power)
            .and_then(|layer| layer.factors.get(index?))
            .copied()
            .ok_or_else(|| Error::plain("partition Loewy product witness index is out of range"))
    }
}

pub(super) fn integral_row(row: &Coordinates) -> Result<Vector> {
    Ok(clear_denominators(row)?.0)
}

pub(super) fn clear_denominators(row: &Coordinates) -> Result<(Vector, i128)> {
    let common = row.values().try_fold(1i128, |common, coefficient| {
        let factor = gcd(common as u128, coefficient.denominator as u128) as i128;
        checked_mul(common / factor, coefficient.denominator)
    })?;
    let vector = row
        .iter()
        .map(|(&index, coefficient)| {
            let factor = common
                .checked_div(coefficient.denominator)
                .ok_or_else(arithmetic)?;
            Ok((index, checked_mul(coefficient.numerator, factor)?))
        })
        .collect::<Result<_>>()?;
    Ok((vector, common))
}
