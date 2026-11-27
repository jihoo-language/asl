use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;

use super::super::linear::{self, Echelon, Matrix};
use super::super::Monomial;
use super::{indexing, storage, Chain, Form, Grade, Rational, Tensor};
use crate::hyperelliptic_kummer_surface::hodge_wick::rational::LIMBS;

type BasisElement = (u32, Monomial);

struct Axis {
    elements: Vec<BasisElement>,
    cells: usize,
}

impl Axis {
    fn new(model: &mut Model, form: &Form, grade: &Grade) -> Result<Self> {
        let mut axis = Self {
            elements: Vec::new(),
            cells: 0,
        };
        if form.supports_koszul(grade) {
            let variables = form.variables();
            let p = grade.homological as usize;
            let mut mask = (1u32 << p) - 1;
            while mask < (1u32 << variables) {
                model.spend(1)?;
                if let Some(coefficient) = form.coefficient_grade(grade, mask)? {
                    let slice = form.slice(model, &coefficient)?;
                    for index in slice.pivots() {
                        let monomial = &slice.generators[index];
                        axis.cells += 10 + 5 * monomial.len();
                        storage(axis.cells)?;
                        axis.elements.push((mask, monomial.clone()));
                    }
                }
                let Some(next) = next_mask(mask) else { break };
                mask = next;
            }
        }
        Ok(axis)
    }

    fn width(&self) -> usize {
        self.elements.len()
    }

    fn locate(&self, mask: u32, powers: &[super::Power]) -> Result<usize> {
        self.elements
            .binary_search_by(|(gen, monomial)| {
                gen.cmp(&mask).then_with(|| monomial.as_slice().cmp(powers))
            })
            .map_err(|_| indexing())
    }

    fn coordinates(&self, chain: &Chain) -> Result<Vec<Rational>> {
        let mut result = vec![Rational::zero(); self.width()];
        for (&mask, value) in &chain.terms {
            for (powers, &coefficient) in &value.terms {
                result[self.locate(mask, powers)?] = coefficient;
            }
        }
        Ok(result)
    }

    fn representative(&self, model: &mut Model, row: &[Rational]) -> Result<Chain> {
        let mut chain = Chain::default();
        for ((mask, powers), &coefficient) in self.elements.iter().zip(row) {
            if coefficient.is_zero() {
                continue;
            }
            let mut term = Tensor::default();
            term.insert(model, powers.clone(), coefficient)?;
            chain.insert(model, *mask, &term, false)?;
        }
        Ok(chain)
    }
}

// Gosper's increasing fixed-popcount enumeration.
fn next_mask(mask: u32) -> Option<u32> {
    if mask == 0 {
        return None;
    }
    let least = mask & mask.wrapping_neg();
    let higher = mask + least;
    let trailing = ((higher ^ mask) >> 2) / least;
    Some(higher | trailing)
}

fn differential(
    model: &mut Model,
    form: &Form,
    grade: &Grade,
    source: &Axis,
    target: &Axis,
    transpose: bool,
) -> Result<Matrix> {
    let (rows, columns) = if transpose {
        (source.width(), target.width())
    } else {
        (target.width(), source.width())
    };
    let mut matrix = linear::matrix(rows, columns)?;
    if rows == 0 || columns == 0 {
        return Ok(matrix);
    }
    for (column, (mask, powers)) in source.elements.iter().enumerate() {
        let chain = Chain::term(model, *mask, &Tensor::monomial(powers.clone()))?;
        let image = form.koszul_diff(model, grade, &chain)?;
        for (&mask, tensor) in &image.terms {
            for (powers, &coefficient) in &tensor.terms {
                let row = target.locate(mask, powers)?;
                let (i, j) = if transpose {
                    (column, row)
                } else {
                    (row, column)
                };
                matrix[i][j] = coefficient;
            }
        }
    }
    Ok(matrix)
}

pub(crate) struct Homology {
    axis: Axis,
    cycle_constraints: Echelon,
    boundaries: Echelon,
    basis: Matrix,
}

impl Homology {
    fn new(model: &mut Model, form: &Form, grade: &Grade) -> Result<Self> {
        let axis = Axis::new(model, form, grade)?;
        let columns = axis.width();
        if columns == 0 {
            return Ok(Self {
                axis,
                cycle_constraints: Echelon {
                    rows: Vec::new(),
                    pivots: Vec::new(),
                },
                boundaries: Echelon {
                    rows: Vec::new(),
                    pivots: Vec::new(),
                },
                basis: Vec::new(),
            });
        }
        let previous_grade = grade.differential()?;
        let previous = Axis::new(model, form, &previous_grade)?;
        let next_grade = Grade {
            homological: grade.homological + 1,
            ..grade.clone()
        };
        let next = Axis::new(model, form, &next_grade)?;
        let outgoing = differential(model, form, grade, &axis, &previous, false)?;
        let cycle_constraints = linear::rref(model, outgoing, columns)?;
        let incoming = differential(model, form, &next_grade, &next, &axis, true)?;
        let boundaries = linear::rref(model, incoming, columns)?;
        for row in &boundaries.rows {
            if !is_cycle(model, &cycle_constraints, row)? {
                return Err(Error::plain(
                    "Verbitsky Koszul differential does not square to zero",
                ));
            }
        }
        let mut cycles = cycle_constraints.kernel(columns)?;
        for row in &mut cycles {
            reduce_boundaries(model, &boundaries, row)?;
        }
        let basis = linear::rref(model, cycles, columns)?.rows;
        if cycle_constraints.rows.len() + boundaries.rows.len() + basis.len() != columns {
            return Err(Error::plain(
                "Verbitsky Koszul cycles and boundaries have inconsistent ranks",
            ));
        }
        Ok(Self {
            axis,
            cycle_constraints,
            boundaries,
            basis,
        })
    }

    fn cells(&self) -> usize {
        let width = self.axis.width() * (2 * LIMBS + 1) + 8;
        let rows =
            self.basis.len() + self.boundaries.rows.len() + self.cycle_constraints.rows.len();
        self.axis.cells + width * rows + 16
    }

    fn project(&self, model: &mut Model, chain: &Chain) -> Result<Chain> {
        let mut row = self.axis.coordinates(chain)?;
        if !is_cycle(model, &self.cycle_constraints, &row)? {
            return Err(Error::plain(
                "Verbitsky Tor projection requires a Koszul cycle",
            ));
        }
        reduce_boundaries(model, &self.boundaries, &mut row)?;
        self.axis.representative(model, &row)
    }
}

fn is_cycle(model: &mut Model, constraint: &Echelon, row: &[Rational]) -> Result<bool> {
    for test in &constraint.rows {
        let mut sum = Rational::zero();
        model.spend(test.len())?;
        for (&a, &b) in test.iter().zip(row) {
            let product = a.multiply(model, b)?;
            sum = sum.add(model, product)?;
        }
        if !sum.is_zero() {
            return Ok(false);
        }
    }
    Ok(true)
}

fn reduce_boundaries(model: &mut Model, boundaries: &Echelon, row: &mut [Rational]) -> Result<()> {
    for (&pivot, other) in boundaries.pivots.iter().zip(&boundaries.rows) {
        let correction = row[pivot].negate();
        linear::add_multiple(model, row, other, correction)?;
    }
    Ok(())
}

impl Form {
    fn homology(&self, model: &mut Model, grade: &Grade) -> Result<Rc<Homology>> {
        let key = (grade.length, grade.homological, grade.weight);
        if let Some(value) = self.koszul.borrow().get(&key).cloned() {
            return Ok(value);
        }
        let slice = Homology::new(model, self, grade)?;
        storage(slice.cells())?;
        model.allocate(slice.cells())?;
        let value = Rc::new(slice);
        self.koszul.borrow_mut().insert(key, value.clone());
        Ok(value)
    }

    pub(crate) fn tor_project(
        &self,
        model: &mut Model,
        grade: &Grade,
        cycle: &Chain,
    ) -> Result<Chain> {
        if cycle.terms.is_empty() {
            return Ok(Chain::default());
        }
        self.homology(model, grade)?.project(model, cycle)
    }

    pub(crate) fn tor_betti(&self, model: &mut Model, grade: &Grade) -> Result<usize> {
        if !self.supports_koszul(grade) {
            return Ok(0);
        }
        Ok(self.homology(model, grade)?.basis.len())
    }

    pub(crate) fn tor_basis(
        &self,
        model: &mut Model,
        grade: &Grade,
        index: usize,
    ) -> Result<Chain> {
        let homology = self.homology(model, grade)?;
        let row = homology
            .basis
            .get(index)
            .ok_or_else(|| Error::plain("Verbitsky Tor basis index is outside the Betti number"))?;
        homology.axis.representative(model, row)
    }
}
