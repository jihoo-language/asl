use std::collections::BTreeMap;

use crate::binary_field::BinaryField;
use crate::binary_linear::add_multiple;
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId};
use crate::rank_variety::{Column, Endpoint};

pub(crate) struct Report {
    pub(crate) dimension: usize,
    pub(crate) blocks: BTreeMap<usize, usize>,
}

pub(crate) fn compute(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    selected: &[usize],
    modulus: usize,
    order: usize,
    coefficients: &[Vec<usize>],
) -> Result<Report> {
    if order == 0 {
        return Err(Error::plain("rvjet requires a positive truncation order"));
    }
    let field = BinaryField::new(model, modulus)?;
    model.check_q_endpoints(source, target, selected)?;
    if selected.len() != coefficients.len() {
        return Err(Error::plain(
            "rvjet requires one jet per selected primitive",
        ));
    }
    let mut active = Vec::new();
    let mut arcs = Vec::new();
    let mut cells = 0;
    for (&primitive, polynomial) in selected.iter().zip(coefficients) {
        let mut arc = Vec::new();
        model.spend(1 + polynomial.len())?;
        for (power, &coefficient) in polynomial.iter().enumerate() {
            if coefficient >= field.order() {
                return Err(Error::plain(format!(
                    "field coefficient {coefficient} is out of range for F{}",
                    field.order()
                )));
            }
            if coefficient == 0 {
                continue;
            }
            if power >= order {
                return Err(Error::plain(
                    "rvjet has a nonzero coefficient beyond the truncation order",
                ));
            }
            cells += 2;
            check_storage(cells)?;
            arc.push((power, coefficient));
        }
        if !arc.is_empty() {
            cells += 2;
            check_storage(cells)?;
            active.push(primitive);
            arcs.push(arc);
        }
    }
    let from = Endpoint::build(model, source, &active)?;
    let other = if source == target {
        None
    } else {
        Some(Endpoint::build(model, target, &active)?)
    };
    let from = from.columns.as_slice();
    let to = other
        .as_ref()
        .map_or(from, |target| target.columns.as_slice());
    let hom = from
        .len()
        .checked_mul(to.len())
        .ok_or_else(dimension_overflow)?;
    let dimension = hom.checked_mul(order).ok_or_else(dimension_overflow)?;
    if dimension == 0 || from.iter().chain(to).all(Vec::is_empty) {
        return Ok(Report {
            dimension,
            blocks: if hom == 0 {
                BTreeMap::new()
            } else {
                BTreeMap::from([(order, hom)])
            },
        });
    }

    check_storage(dimension)?;
    // Store the weighted transposed source action for precomposition.
    cells += 2 * action_cells(from);
    if other.is_some() {
        cells += action_cells(to);
    }
    check_storage(cells)?;
    let mut transpose = vec![Vec::new(); from.len()];
    model.spend(from.len())?;
    for (domain, column) in from.iter().enumerate() {
        model.spend(column.len())?;
        for &(row, variable) in column {
            transpose[row].push((domain, variable));
        }
    }
    let action = HomAction {
        source: &transpose,
        target: to,
        arcs: &arcs,
        hom,
        order,
        dimension,
        cells,
    };
    action.homology(model, field)
}

fn action_cells(columns: &[Column]) -> usize {
    columns.len() + columns.iter().map(|column| 2 * column.len()).sum::<usize>()
}

struct HomAction<'a> {
    source: &'a [Column],
    target: &'a [Column],
    arcs: &'a [Column],
    hom: usize,
    order: usize,
    dimension: usize,
    cells: usize,
}

impl HomAction<'_> {
    fn column(&self, model: &mut Model, index: usize) -> Result<Column> {
        let power = index / self.hom;
        let source = (index % self.hom) / self.target.len();
        let target = index % self.target.len();
        let source_terms = self.source[source]
            .iter()
            .map(|&(row, variable)| (row * self.target.len() + target, variable));
        let target_terms = self.target[target]
            .iter()
            .map(|&(row, variable)| (source * self.target.len() + row, variable));

        let mut result = Vec::new();
        for (row, variable) in source_terms.chain(target_terms) {
            model.spend(1)?;
            for &(shift, coefficient) in &self.arcs[variable] {
                model.spend(1)?;
                if shift >= self.order - power {
                    break;
                }
                check_storage(self.cells + 2 * (result.len() + 1))?;
                result.push(((power + shift) * self.hom + row, coefficient));
            }
        }
        model.spend(result.len() * (1 + result.len().max(1).ilog2() as usize))?;
        // Endpoints change in disjoint ways under the two actions.
        // Distinct primitive degrees also land in distinct endpoint rows.
        result.sort_unstable();
        Ok(result)
    }

    fn homology(&self, model: &mut Model, field: BinaryField) -> Result<Report> {
        let mut image = Basis::default();
        let mut kernel = Vec::new();
        let mut kernel_cells = 0;
        for index in 0..self.dimension {
            model.spend(1)?;
            let column = self.column(model, index)?;
            if let Some(cycle) = image.insert(
                model,
                field,
                column,
                vec![(index, 1)],
                self.cells + kernel_cells,
            )? {
                kernel_cells += 1 + 2 * cycle.len();
                check_storage(self.cells + image.cells + kernel_cells)?;
                kernel.push(cycle);
            }
        }
        let dimension = self
            .dimension
            .checked_sub(2 * image.rank())
            .ok_or_else(incompatible_square)?;
        if dimension == 0 {
            return Ok(Report {
                dimension,
                blocks: BTreeMap::new(),
            });
        }
        image.drop_witnesses(model)?;
        let fixed = self.cells + image.cells;
        let mut span = Basis::default();
        while let Some(cycle) = kernel.pop() {
            kernel_cells -= 1 + 2 * cycle.len();
            let reduced = image.reduce(model, field, cycle)?;
            span.insert(model, field, reduced, Vec::new(), fixed + kernel_cells)?;
        }
        if span.rank() != dimension {
            return Err(incompatible_square());
        }

        // At step j, span represents u^j H in the quotient by im(D).
        let mut ranks = vec![dimension];
        for _ in 0..self.order {
            let mut next = Basis::default();
            for pivot in span.pivots.values() {
                model.spend(1 + pivot.column.len())?;
                let shifted = pivot
                    .column
                    .iter()
                    .take_while(|&&(row, _)| row < self.dimension - self.hom)
                    .map(|&(row, coefficient)| (row + self.hom, coefficient))
                    .collect();
                let reduced = image.reduce(model, field, shifted)?;
                next.insert(
                    model,
                    field,
                    reduced,
                    Vec::new(),
                    fixed + span.cells + ranks.len(),
                )?;
            }
            ranks.push(next.rank());
            if next.rank() == 0 {
                break;
            }
            span = next;
        }
        if ranks.last() != Some(&0) {
            return Err(incompatible_nilpotence());
        }
        ranks.push(0);
        let mut blocks = BTreeMap::new();
        for (length, triple) in ranks.windows(3).enumerate() {
            let [before, at, after] = *triple else {
                unreachable!()
            };
            let multiplicity = before
                .checked_sub(at)
                .zip(at.checked_sub(after))
                .and_then(|(entering, leaving)| entering.checked_sub(leaving))
                .ok_or_else(incompatible_nilpotence)?;
            if multiplicity != 0 {
                blocks.insert(length + 1, multiplicity);
            }
        }
        Ok(Report { dimension, blocks })
    }
}

#[derive(Default)]
struct Basis {
    pivots: BTreeMap<usize, Pivot>,
    cells: usize,
}

struct Pivot {
    column: Column,
    witness: Column,
}

impl Basis {
    fn rank(&self) -> usize {
        self.pivots.len()
    }

    fn insert(
        &mut self,
        model: &mut Model,
        field: BinaryField,
        mut column: Column,
        mut witness: Column,
        retained: usize,
    ) -> Result<Option<Column>> {
        while let Some(&(leading, coefficient)) = column.last() {
            model.spend(1)?;
            if let Some(pivot) = self.pivots.get(&leading) {
                column = add_multiple(model, field, &column, &pivot.column, coefficient)?;
                if !pivot.witness.is_empty() {
                    witness = add_multiple(model, field, &witness, &pivot.witness, coefficient)?;
                }
            } else {
                self.cells += 1 + 2 * (column.len() + witness.len());
                check_storage(retained + self.cells)?;
                if coefficient != 1 {
                    let scale = field.inverse(model, coefficient)?;
                    model.spend(field.degree() * (column.len() + witness.len()))?;
                    for (_, value) in column.iter_mut().chain(witness.iter_mut()) {
                        *value = field.multiply(scale, *value);
                    }
                }
                self.pivots.insert(leading, Pivot { column, witness });
                return Ok(None);
            }
        }
        Ok(Some(witness))
    }

    fn drop_witnesses(&mut self, model: &mut Model) -> Result<()> {
        model.spend(self.rank())?;
        for pivot in self.pivots.values_mut() {
            self.cells -= 2 * pivot.witness.len();
            pivot.witness = Vec::new();
        }
        Ok(())
    }

    fn reduce(&self, model: &mut Model, field: BinaryField, mut column: Column) -> Result<Column> {
        let mut remainder = Vec::new();
        while let Some(&(leading, coefficient)) = column.last() {
            model.spend(1)?;
            if let Some(pivot) = self.pivots.get(&leading) {
                column = add_multiple(model, field, &column, &pivot.column, coefficient)?;
            } else {
                remainder.push(column.pop().unwrap());
            }
        }
        remainder.reverse();
        Ok(remainder)
    }
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "Margolis jet calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

fn dimension_overflow() -> Error {
    Error::plain("Margolis jet dimension is out of range")
}

fn incompatible_square() -> Error {
    Error::plain("Margolis jet homology rank is incompatible with square zero")
}

fn incompatible_nilpotence() -> Error {
    Error::plain("Margolis jet homology ranks are incompatible with the truncation of u")
}
