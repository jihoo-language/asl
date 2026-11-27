use crate::error::{Error, Result};
use crate::hw_algebra::{self, Matrix};
use crate::hyperelliptic::Residues;
use crate::model::Model;
use crate::odd_cobar::Prime;

use super::{
    coordinate_range, grade_range, insert_factor, remove_factor, Accumulator, Factor, Gram,
    Integer, Mode, State, COLORS, LIMBS,
};

impl Gram {
    pub(crate) fn inverse(&self, model: &mut Model) -> Result<Self> {
        // This prime admits a compact centered lift for the fixed marked basis.
        // Verify the lift over Z, without relying on a bound for its entries.
        let ring = Residues::new(Prime::validated(2_147_483_647), 1)?;
        let matrix = Matrix::from_entries(
            COLORS,
            self.0
                .iter()
                .flatten()
                .map(|&entry| ring.coefficient(entry))
                .collect(),
        )?;
        let inverse = hw_algebra::inverse(model, &ring, &matrix)?.ok_or_else(|| {
            Error::plain("Nakajima Gram matrix is singular modulo the lifting prime")
        })?;
        let center = |x: u64| {
            if x > ring.modulus / 2 {
                i128::from(x) - i128::from(ring.modulus)
            } else {
                i128::from(x)
            }
        };
        let entries: [[i128; COLORS]; COLORS] =
            std::array::from_fn(|i| std::array::from_fn(|j| center(*inverse.entry(i, j))));
        model.spend(COLORS.pow(3))?;
        for (i, row) in self.0.iter().enumerate() {
            if entries[i][i] % 2 != 0 {
                return Err(Error::plain(
                    "Nakajima inverse Gram matrix has an odd diagonal",
                ));
            }
            for j in 0..COLORS {
                // The original Gram uses a fixed 22-dimensional basis of doubled
                // coordinates in {0,1,2}; the lifted entries have at most 31 bits.
                let product: i128 = row
                    .iter()
                    .zip(&entries)
                    .map(|(&value, column)| value * column[j])
                    .sum();
                if product != i128::from(i == j) {
                    return Err(Error::plain(
                        "Nakajima integral inverse Gram certificate failed",
                    ));
                }
            }
        }
        model.allocate(COLORS * COLORS)?;
        Ok(Self(entries))
    }
}

impl State {
    pub(crate) fn virasoro(&self, model: &mut Model, index: i128, length: i128) -> Result<Self> {
        if self.terms.is_empty() || index > length {
            return Ok(Self::zero());
        }
        if index == 0 {
            return self.scale(model, length);
        }
        let added = if index < -1 {
            let added = index.checked_neg().ok_or_else(grade_range)?;
            // Every pair of lengths has nonzero inverse-Gram entries. Charge a
            // lower bound before processing an action with new creation modes.
            let pairs = usize::try_from(added / 2).unwrap_or(usize::MAX);
            model.spend(pairs.saturating_mul(self.terms.len()))?;
            Some(added)
        } else {
            None
        };
        let mut output = Accumulator::default();
        for (monomial, &coefficient) in &self.terms {
            move_modes(model, &mut output, monomial, coefficient, index)?;
        }
        if index > 1 {
            let gram = model.nakajima_gram()?;
            for (monomial, &coefficient) in &self.terms {
                contract_pairs(model, &mut output, &gram, monomial, coefficient, index)?;
            }
        } else if let Some(added) = added {
            let inverse = model.nakajima_inverse_gram()?;
            let (distinct, equal) = pair_colors(&inverse);
            for left in 1..=added / 2 {
                let right = added - left;
                let colors = if left == right { &equal } else { &distinct };
                for (monomial, &coefficient) in &self.terms {
                    for &(a, b, inverse_entry) in colors {
                        model.spend(1 + monomial.len() + LIMBS)?;
                        let mut key = monomial.clone();
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
                        output.add(key, Integer::product(coefficient, inverse_entry))?;
                    }
                }
            }
        }
        output.finish()
    }
}

fn move_modes(
    model: &mut Model,
    output: &mut Accumulator,
    monomial: &[Factor],
    coefficient: i128,
    index: i128,
) -> Result<()> {
    model.spend(1 + monomial.len())?;
    for (position, factor) in monomial.iter().enumerate() {
        let length = factor
            .mode
            .length
            .checked_sub(index)
            .ok_or_else(grade_range)?;
        if length <= 0 {
            continue;
        }
        model.spend(1 + monomial.len() + LIMBS)?;
        let mut key = monomial.to_vec();
        remove_factor(&mut key, position);
        insert_factor(
            &mut key,
            Mode {
                length,
                color: factor.mode.color,
            },
        )?;
        output.add(
            key,
            Integer::product(coefficient, length_weight(factor, factor.exponent)?),
        )?;
    }
    Ok(())
}

fn contract_pairs(
    model: &mut Model,
    output: &mut Accumulator,
    gram: &Gram,
    monomial: &[Factor],
    coefficient: i128,
    index: i128,
) -> Result<()> {
    for (i, left) in monomial.iter().enumerate() {
        model.spend(1 + monomial.len().ilog2() as usize)?;
        let partner = index - left.mode.length;
        if partner < left.mode.length {
            continue;
        }
        let start = monomial.partition_point(|f| f.mode.length < partner).max(i);
        for (j, right) in monomial
            .iter()
            .enumerate()
            .skip(start)
            .take_while(|(_, f)| f.mode.length == partner)
        {
            model.spend(1)?;
            let pairing = gram.entry(left.mode.color, right.mode.color);
            if pairing == 0 || (i == j && left.exponent == 1) {
                continue;
            }
            model.spend(1 + monomial.len() + LIMBS)?;
            let (mut a, mut b) = (left.exponent, right.exponent);
            if i == j {
                // For identical keys take the number of unordered copies.
                b -= 1;
                if a.is_multiple_of(2) {
                    a /= 2;
                } else {
                    b /= 2;
                }
            }
            let weight = Integer::product(length_weight(left, a)?, length_weight(right, b)?)
                .multiply(pairing)
                .and_then(|term| term.multiply(coefficient))
                .map_err(|_| coordinate_range())?;
            let mut key = monomial.to_vec();
            remove_factor(&mut key, j);
            remove_factor(&mut key, i);
            output.add(key, weight)?;
        }
    }
    Ok(())
}

fn length_weight(factor: &Factor, count: u128) -> Result<i128> {
    // The input length grade already bounds the sum of all occupied lengths.
    let count = i128::try_from(count).map_err(|_| coordinate_range())?;
    factor
        .mode
        .length
        .checked_mul(count)
        .ok_or_else(coordinate_range)
}

type ColorPair = (u8, u8, i128);

fn pair_colors(inverse: &Gram) -> (Vec<ColorPair>, Vec<ColorPair>) {
    let mut distinct = Vec::new();
    let mut equal = Vec::new();
    for a in 0..COLORS as u8 {
        for b in 0..COLORS as u8 {
            let entry = inverse.entry(a, b);
            if entry != 0 {
                distinct.push((a, b, entry));
                if a < b {
                    equal.push((a, b, entry));
                } else if a == b {
                    equal.push((a, b, entry / 2));
                }
            }
        }
    }
    (distinct, equal)
}
