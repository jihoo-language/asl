use crate::error::Result;
use crate::model::Model;

use super::{
    coordinate_range, grade_range, insert_factor, length_weight, pair_weight, remove_factor,
    Accumulator, Factor, Gram, Integer, Mode, State, COLORS, LIMBS, POINT,
};

pub(super) type ColorPair = (u8, u8, i128);

impl State {
    pub(crate) fn boundary(&self, model: &mut Model) -> Result<Self> {
        if self.terms.is_empty() {
            return Ok(Self::zero());
        }
        let mut split_pairs = 0usize;
        let mut split_unit = false;
        for monomial in self.terms.keys() {
            model.spend(1 + monomial.len())?;
            for factor in monomial {
                let mode = factor.mode;
                let pairs = usize::try_from(mode.length / 2).unwrap_or(usize::MAX);
                split_pairs = split_pairs.saturating_add(pairs);
                split_unit |= pairs != 0 && mode.color == 0;
            }
        }
        // Each pair of positive lengths has at least one nonzero diagonal term.
        model.spend(split_pairs)?;
        let unit_diagonal = if split_unit {
            let inverse = model.nakajima_inverse_gram()?;
            diagonal_of_unit(model, &inverse)?
        } else {
            Vec::new()
        };
        let gram = model.nakajima_gram()?;
        let mut output = Accumulator::default();
        for (monomial, &coefficient) in &self.terms {
            join(model, &mut output, &gram, monomial, coefficient)?;
            split(model, &mut output, &unit_diagonal, monomial, coefficient)?;
        }
        output.finish()
    }
}

pub(super) fn diagonal_of_unit(model: &mut Model, inverse: &Gram) -> Result<Vec<ColorPair>> {
    model.spend(COLORS * COLORS)?;
    let mut terms = Vec::new();
    for left in 0..COLORS as u8 {
        for right in 0..COLORS as u8 {
            let weight = inverse.entry(left, right);
            if weight != 0 {
                terms.push((left, right, weight));
            }
        }
    }
    Ok(terms)
}

pub(super) fn surface_product(gram: &Gram, left: u8, right: u8) -> (u8, i128) {
    match (left, right) {
        (0, color) | (color, 0) => (color, 1),
        (POINT, _) | (_, POINT) => (POINT, 0),
        (left, right) => (POINT, gram.entry(left, right)),
    }
}

fn join(
    model: &mut Model,
    output: &mut Accumulator,
    gram: &Gram,
    monomial: &[Factor],
    coefficient: i128,
) -> Result<()> {
    for (i, left) in monomial.iter().enumerate() {
        for (j, right) in monomial.iter().enumerate().skip(i) {
            model.spend(1)?;
            let (color, product) = surface_product(gram, left.mode.color, right.mode.color);
            if product == 0 || (i == j && left.exponent == 1) {
                continue;
            }
            model.spend(1 + monomial.len() + LIMBS)?;
            let contribution = pair_weight(left, right, i == j)?
                .multiply(product)
                .and_then(|value| value.multiply(coefficient))
                .map_err(|_| coordinate_range())?;
            let length = left
                .mode
                .length
                .checked_add(right.mode.length)
                .ok_or_else(grade_range)?;
            let mut key = monomial.to_vec();
            remove_factor(&mut key, j);
            remove_factor(&mut key, i);
            insert_factor(&mut key, Mode { length, color })?;
            output.add(key, contribution)?;
        }
    }
    Ok(())
}

fn split(
    model: &mut Model,
    output: &mut Accumulator,
    unit_diagonal: &[ColorPair],
    monomial: &[Factor],
    coefficient: i128,
) -> Result<()> {
    for (position, factor) in monomial.iter().enumerate() {
        let Mode { length, color } = factor.mode;
        if length == 1 {
            continue;
        }
        let point_terms = [(color, POINT, 1), (POINT, color, 1)];
        let diagonal = match color {
            0 => unit_diagonal,
            POINT => &point_terms[..1],
            _ => &point_terms[..],
        };
        let multiplicity = length_weight(factor, factor.exponent)?;
        let single = Integer::product(coefficient, -multiplicity);
        let half = Integer::product(coefficient, -(multiplicity / 2));
        for left in 1..=length / 2 {
            let right = length - left;
            for &(a, b, weight) in diagonal {
                model.spend(1)?;
                if left == right && a > b {
                    continue;
                }
                // A pair unchanged by swapping the two newborn factors still
                // carries the factor 1/2. Equal lengths make it integral.
                let contribution = if left == right && a == b {
                    half
                } else {
                    single
                }
                .multiply(weight)
                .map_err(|_| coordinate_range())?;
                model.spend(1 + monomial.len() + LIMBS)?;
                let mut key = monomial.to_vec();
                remove_factor(&mut key, position);
                insert_factor(
                    &mut key,
                    Mode {
                        length: left,
                        color: a,
                    },
                )?;
                insert_factor(
                    &mut key,
                    Mode {
                        length: right,
                        color: b,
                    },
                )?;
                output.add(key, contribution)?;
            }
        }
    }
    Ok(())
}
