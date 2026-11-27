use std::collections::BTreeMap;
use std::rc::Rc;

use super::clifford::{self, NEGATIVE_SQUARES, VOLUME_MASK};
use super::rational_sparse::{Coefficients, Domain, Integer, WORK_LIMBS};
use crate::error::Result;
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

pub(crate) const DIMENSION: usize = 1 << 10;
const INDEX_MASK: u16 = DIMENSION as u16 - 1;
const LAST_GRADE: u16 = DIMENSION as u16;
const GRADE_BITS: u16 = 0b111 | (1 << 6) | LAST_GRADE;

pub(crate) enum Expression {
    Basis(Label, u8, u32),
    Scale(i128, u128, Box<ParsedExpr>),
    Act(Box<[ParsedExpr; 2]>),
}

pub(crate) struct Coefficient(pub(crate) u32, pub(crate) Box<ParsedExpr>);

pub(crate) enum Spinor {}

impl Domain for Spinor {
    const NAME: &'static str = "K3 spinor";
}

pub(crate) type Element = Coefficients<Spinor>;

impl Element {
    pub(crate) fn basis(index: u32) -> Self {
        Self::singleton(index, 1)
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        self.charge_storage(model)?;
        Ok(Value::K3Spinor(Rc::new(self)))
    }

    pub(crate) fn act(
        &self,
        model: &mut Model,
        actor: &clifford::Element,
        chirality: u8,
    ) -> Result<Self> {
        if self.terms().is_empty() || actor.terms().is_empty() {
            return Ok(Self::zero());
        }
        let actor_work: usize = actor
            .terms()
            .iter()
            .map(|&(mask, _)| 4 * WORK_LIMBS + 2 * mask.count_ones() as usize)
            .sum();
        model
            .spend(actor_work.saturating_mul(self.terms().len()))
            .map_err(|_| Self::work_limit())?;
        let states = self
            .terms()
            .iter()
            .map(|&(index, coefficient)| {
                Ok((
                    basis_state(chirality, index as u16),
                    coefficient.resize::<WORK_LIMBS>()?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut output = BTreeMap::new();
        for &(mask, coefficient) in actor.terms() {
            let coefficient: Integer = coefficient.resize()?;
            for &(state, value) in &states {
                let (image, phase) = monomial_action(mask, state);
                let coefficient = coefficient
                    .multiply_wide(value)?
                    .multiply(i128::from(phase))?;
                Self::accumulate(&mut output, u32::from(image & INDEX_MASK), coefficient)?;
            }
        }
        let denominator = self.product_denominator(actor)?;
        Self::normalized(
            model,
            output
                .iter()
                .map(|(&index, coefficient)| (index, coefficient)),
            denominator,
        )
    }
}

const fn grading(state: u16) -> u8 {
    ((state & GRADE_BITS).count_ones() & 1) as u8
}

// The highest state bit records the last block's grade. Omit it in the
// compact index, where the requested volume eigenvalue determines it.
const VOLUME_ORIENTATION: u8 = (monomial_action(VOLUME_MASK, 0).1 < 0) as u8;

const fn basis_state(chirality: u8, index: u16) -> u16 {
    if chirality ^ VOLUME_ORIENTATION ^ grading(index) == 0 {
        index
    } else {
        index | LAST_GRADE
    }
}

const fn monomial_action(mask: u32, mut state: u16) -> (u16, i16) {
    let mut phase = 1;
    let mut remaining = mask;
    while remaining != 0 {
        // The rightmost factor of an increasing, ordered product acts first.
        let slot = u32::BITS - 1 - remaining.leading_zeros();
        remaining ^= 1 << slot;
        let (image, coefficient) = gamma(slot, state);
        state = image;
        phase *= coefficient;
    }
    (state, phase)
}

const fn gamma(slot: u32, state: u16) -> (u16, i16) {
    let (image, mut phase, grade_bit) = if slot < 6 {
        let grade_bit = 1 << (slot / 2);
        let mut phase = 1;
        if state & grade_bit != 0 {
            phase = if NEGATIVE_SQUARES & (1 << slot) == 0 {
                2
            } else {
                -2
            };
        }
        (state ^ grade_bit, phase, grade_bit)
    } else {
        let (offset, unit) = if slot < 14 {
            (3, slot - 5)
        } else {
            (7, slot - 13)
        };
        let grade_bit = 1 << (offset + 3);
        if unit == 8 {
            // The last gamma sends the upper half to minus the lower half.
            let phase = if state & grade_bit == 0 { -1 } else { 1 };
            (state ^ grade_bit, phase, grade_bit)
        } else {
            let octonion = ((state >> offset) & 7) as usize;
            let phase = OCTONION_SIGN[unit as usize][octonion];
            (
                state ^ grade_bit ^ ((unit as u16) << offset),
                phase,
                grade_bit,
            )
        }
    };
    if (state & GRADE_BITS & (grade_bit - 1)).count_ones() & 1 != 0 {
        phase = -phase;
    }
    (image, phase)
}

// Octonion basis products have index i XOR j. Oriented triples specify
// the positive products of distinct imaginary units in cyclic order.
const OCTONION_SIGN: [[i16; 8]; 8] = {
    let triples = [
        [1, 2, 3],
        [1, 4, 5],
        [1, 7, 6],
        [2, 4, 6],
        [2, 5, 7],
        [3, 4, 7],
        [3, 6, 5],
    ];
    let mut signs = [[1; 8]; 8];
    let mut i = 1;
    while i < 8 {
        signs[i][i] = -1;
        i += 1;
    }
    let mut line = 0;
    while line < triples.len() {
        let [i, j, k] = triples[line];
        signs[i][j] = 1;
        signs[j][k] = 1;
        signs[k][i] = 1;
        signs[j][i] = -1;
        signs[k][j] = -1;
        signs[i][k] = -1;
        line += 1;
    }
    signs
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hyperelliptic_kummer_surface::unimodular::RANK;

    #[test]
    fn full_gamma_table_and_both_volume_eigenspaces() {
        assert_eq!(VOLUME_ORIENTATION, 0);
        let squares = [
            -2, 2, 2, -2, -2, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        ];
        let mut seen = [false; 2 * DIMENSION];
        for chirality in 0..2 {
            for index in 0..DIMENSION as u16 {
                let state = basis_state(chirality, index);
                assert_eq!(state & INDEX_MASK, index);
                assert!(!seen[state as usize]);
                seen[state as usize] = true;
                let eigenvalue = if chirality == 0 { 8 } else { -8 };
                assert_eq!(monomial_action(VOLUME_MASK, state), (state, eigenvalue));
                for i in 0..RANK as u32 {
                    let (after_i, coefficient_i) = gamma(i, state);
                    assert_eq!(basis_state(chirality ^ 1, after_i & INDEX_MASK), after_i);
                    let (back, again_i) = gamma(i, after_i);
                    assert_eq!(
                        (back, coefficient_i * again_i),
                        (state, squares[i as usize])
                    );
                    for j in (i + 1)..RANK as u32 {
                        let (after_j, coefficient_j) = gamma(j, state);
                        let (left, j_after_i) = gamma(j, after_i);
                        let (right, i_after_j) = gamma(i, after_j);
                        assert_eq!(left, right);
                        assert_eq!(j_after_i * coefficient_i, -i_after_j * coefficient_j);
                    }
                }
            }
        }
        assert!(seen.into_iter().all(|used| used));
    }
}
