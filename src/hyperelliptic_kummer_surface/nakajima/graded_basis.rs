use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;

use super::{Factor, Grade, Mode, Monomial, COLORS, POINT};

struct WeightedMode {
    mode: Mode,
    length: usize,
    degree: usize,
}

struct Enumeration {
    length: u128,
    complement: bool,
    modes: Vec<WeightedMode>,
    counts: Vec<usize>,
    length_stride: usize,
    mode_stride: usize,
    cells: usize,
    basis: Vec<Monomial>,
}

impl Enumeration {
    fn offset(&self, start: usize, length: usize, degree: usize) -> usize {
        start * self.mode_stride + length * self.length_stride + degree
    }

    fn count(&self, start: usize, length: usize, degree: usize) -> usize {
        self.counts[self.offset(start, length, degree)]
    }

    fn enumerate(
        &mut self,
        model: &mut Model,
        start: usize,
        length: usize,
        degree: usize,
        available: usize,
        key: &mut Monomial,
    ) -> Result<()> {
        if degree == 0 {
            let unused = self.length - (available - length) as u128;
            self.cells += 5 * (key.len() + usize::from(unused != 0));
            check_cells(self.cells)?;
            model.spend(1 + key.len())?;
            let mut monomial = Vec::with_capacity(key.len() + usize::from(unused != 0));
            if unused != 0 {
                monomial.push(Factor {
                    mode: Mode {
                        length: 1,
                        color: 0,
                    },
                    exponent: unused,
                });
            }
            monomial.extend_from_slice(key);
            if self.complement {
                complement_monomial(&mut monomial);
            }
            self.basis.push(monomial);
            return Ok(());
        }
        for slot in start..self.modes.len() {
            model.spend(1)?;
            if self.count(slot, length, degree) == 0 {
                break;
            }
            let WeightedMode {
                mode,
                length: weight,
                degree: order,
            } = self.modes[slot];
            for exponent in 1..=(length / weight).min(degree / order) {
                model.spend(1)?;
                let rest_length = length - exponent * weight;
                let rest_degree = degree - exponent * order;
                if self.count(slot + 1, rest_length, rest_degree) != 0 {
                    key.push(Factor {
                        mode,
                        exponent: exponent as u128,
                    });
                    self.enumerate(model, slot + 1, rest_length, rest_degree, available, key)?;
                    key.pop();
                }
            }
        }
        Ok(())
    }
}

pub(super) fn graded_basis(model: &mut Model, grade: &Grade) -> Result<Vec<Monomial>> {
    if !grade.supports_states() {
        return Ok(Vec::new());
    }
    let length = grade.length as u128;
    let degree = (grade.degree / 2) as u128;
    // Swapping the unit and point sends degree i to its complement 4n-i.
    let complement = degree > length;
    let degree = if complement {
        length - (degree - length)
    } else {
        degree
    };
    if degree == 0 {
        if length == 0 {
            return Ok(vec![Vec::new()]);
        }
        let color = if complement { POINT } else { 0 };
        return Ok(vec![vec![Factor {
            mode: Mode { length: 1, color },
            exponent: length,
        }]]);
    }

    // Each positive-half-degree oscillator uses at most twice its half-degree
    // in length. The remaining length is filled by the single q_1(1) mode.
    let available = usize::try_from(length.min(2 * degree)).map_err(|_| enumeration_limit())?;
    let degree = usize::try_from(degree).map_err(|_| enumeration_limit())?;
    let length_stride = degree.saturating_add(1);
    let mode_stride = available.saturating_add(1).saturating_mul(length_stride);
    let max_mode = available.min(length_stride);
    let mode_count = (max_mode - 1)
        .saturating_add((COLORS - 2).saturating_mul(available.min(degree)))
        .saturating_add(available.min(degree - 1));
    let cells = mode_count.saturating_add(1).saturating_mul(mode_stride);
    let scratch = cells.saturating_add(5usize.saturating_mul(mode_count));
    check_cells(scratch)?;
    model.spend(scratch)?;
    let mut modes = Vec::with_capacity(mode_count);
    for length in 1..=max_mode {
        for color in 0..COLORS as u8 {
            let mode = Mode {
                length: length as i128,
                color,
            };
            let order = length - 1 + (mode.color_degree() / 2) as usize;
            if (1..=degree).contains(&order) {
                modes.push(WeightedMode {
                    mode,
                    length,
                    degree: order,
                });
            }
        }
    }
    let mut enumeration = Enumeration {
        length,
        complement,
        modes,
        counts: vec![0; cells],
        length_stride,
        mode_stride,
        cells: scratch,
        basis: Vec::new(),
    };
    // The suffix table counts monomials with exact degree and at most the
    // available length. Either omit the current mode or remove one copy of it.
    for length in 0..=available {
        let offset = enumeration.offset(mode_count, length, 0);
        enumeration.counts[offset] = 1;
    }
    for slot in (0..mode_count).rev() {
        let WeightedMode {
            length: weight,
            degree: order,
            ..
        } = enumeration.modes[slot];
        for length in 0..=available {
            for degree in 0..=degree {
                let mut count = enumeration.count(slot + 1, length, degree);
                if length >= weight && degree >= order {
                    count += enumeration.count(slot, length - weight, degree - order);
                }
                let offset = enumeration.offset(slot, length, degree);
                enumeration.counts[offset] = count.min(MAX_PIVOT_CELLS + 1);
            }
        }
    }
    let dimension = enumeration.count(0, available, degree);
    enumeration.cells += 3 * dimension;
    check_cells(enumeration.cells)?;
    enumeration.basis = Vec::with_capacity(dimension);
    enumeration.enumerate(model, 0, available, degree, available, &mut Vec::new())?;
    model.spend(dimension * dimension.max(1).ilog2() as usize)?;
    enumeration.basis.sort_unstable();
    Ok(enumeration.basis)
}

fn check_cells(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(enumeration_limit());
    }
    Ok(())
}

fn enumeration_limit() -> Error {
    Error::plain("Nakajima calculation exceeds the graded enumeration storage limit")
}

pub(super) fn complement_monomial(monomial: &mut Monomial) {
    for factor in &mut *monomial {
        factor.mode.color = match factor.mode.color {
            0 => POINT,
            POINT => 0,
            middle => middle,
        };
    }
    monomial.sort_unstable();
}
