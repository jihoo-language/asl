use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;

use super::super::{linear, monomials};
use super::{
    invalid_index, multiply_powers, storage, work_limit, Form, Grade, Monomial, Rational, Tensor,
    Vector,
};

pub(crate) struct Resolution {
    stages: Vec<Stage>,
}

#[derive(Default)]
struct Stage {
    shifts: Vec<usize>,
    columns: Vec<Vector>,
}

impl Resolution {
    pub(super) fn shifts(&self) -> Vec<&[usize]> {
        self.stages
            .iter()
            .map(|stage| stage.shifts.as_slice())
            .collect()
    }

    pub(super) fn shift(&self, p: usize, index: usize) -> Result<usize> {
        self.stages
            .get(p)
            .and_then(|stage| stage.shifts.get(index))
            .copied()
            .ok_or_else(invalid_index)
    }

    pub(super) fn column(&self, p: usize, index: usize) -> Result<&Vector> {
        self.shift(p, index)?;
        Ok(&self.stages[p].columns[index])
    }

    pub(super) fn entry(&self, p: usize, i: usize, j: usize) -> Result<Option<&Tensor>> {
        let prev = p.checked_sub(1).ok_or_else(invalid_index)?;
        self.shift(prev, i)?;
        self.shift(p, j)?;
        Ok(self.stages[p].columns[j].entries.get(&i))
    }

    pub(super) fn differential(
        &self,
        model: &mut Model,
        p: usize,
        input: &Vector,
    ) -> Result<Vector> {
        let stage = self.stages.get(p).ok_or_else(invalid_index)?;
        apply(model, &stage.columns, input)
    }

    fn cells(&self) -> usize {
        self.stages
            .iter()
            .map(|stage| {
                stage.shifts.len()
                    + 8
                    + stage
                        .columns
                        .iter()
                        .map(|column| column.cells + 8)
                        .sum::<usize>()
            })
            .sum()
    }
}

impl Form {
    pub(super) fn resolution(&self, model: &mut Model, grade: &Grade) -> Result<Rc<Resolution>> {
        self.validate_free(grade)?;
        if let Some(found) = self.resolutions.borrow().get(&grade.length).cloned() {
            return Ok(found);
        }
        let count = self.variable_count();
        let bound =
            (2 * grade.length as u128).saturating_add((count + self.field.degree()) as u128);
        let bound = usize::try_from(bound).map_err(|_| work_limit())?;
        storage((count + 1).saturating_mul(bound.saturating_add(1)))?;
        model
            .spend((count + 1) * (bound + 1))
            .map_err(|_| work_limit())?;

        let mut stages = vec![Stage {
            shifts: vec![0],
            columns: Vec::new(),
        }];
        let mut catalog = Catalog::default();
        for p in 1..=count {
            let mut stage = Stage::default();
            for weight in 0..=bound {
                let grade = Grade {
                    homological: p as i128,
                    weight: weight as i128,
                    ..grade.clone()
                };
                let koszul = grade.koszul_grade();
                let betti = self.tor_betti(model, &koszul)?;
                // A zero Tor slice certifies there are no minimal generators
                // of this weight. Avoid expanding their redundant multiples.
                if betti == 0 {
                    continue;
                }
                let columns = if p == 1 {
                    let mut columns = Vec::with_capacity(betti);
                    for index in 0..betti {
                        let cycle = self.tor_basis(model, &koszul, index)?;
                        let polynomial = self.lift_minimal_relation(model, &koszul, &cycle)?;
                        let mut column = Vector::default();
                        column.add_tensor(model, 0, &polynomial)?;
                        columns.push(column);
                    }
                    columns
                } else {
                    syzygies(
                        model,
                        self,
                        &mut catalog,
                        &stages[p - 1],
                        &stages[p - 2],
                        &stage,
                        weight,
                    )?
                };
                if columns.len() != betti {
                    return Err(Error::plain(format!(
                        "Verbitsky minimal syzygy count disagrees with Tor in degree ({p}, {weight})"
                    )));
                }
                for column in columns {
                    if column.entries.is_empty() {
                        return Err(Error::plain(
                            "Verbitsky minimal resolution has a zero generator",
                        ));
                    }
                    certify_minimality(self, &stages[p - 1], weight, &column)?;
                    if p > 1
                        && !apply(model, &stages[p - 1].columns, &column)?
                            .entries
                            .is_empty()
                    {
                        return Err(Error::plain(
                            "Verbitsky free differentials do not compose to zero",
                        ));
                    }
                    stage.shifts.push(weight);
                    stage.columns.push(column);
                }
                storage(stage.columns.iter().map(|column| column.cells + 8).sum())?;
            }
            stages.push(stage);
        }
        let resolution = Resolution { stages };
        storage(resolution.cells())?;
        model.allocate(resolution.cells())?;
        let resolution = Rc::new(resolution);
        self.resolutions
            .borrow_mut()
            .insert(grade.length, resolution.clone());
        Ok(resolution)
    }
}

fn apply(model: &mut Model, columns: &[Vector], input: &Vector) -> Result<Vector> {
    let mut result = Vector::default();
    for (&source, multiplier) in &input.entries {
        let column = columns.get(source).ok_or_else(invalid_index)?;
        for (&target, polynomial) in &column.entries {
            let product = multiplier.cup(model, polynomial)?;
            result.add_tensor(model, target, &product)?;
        }
    }
    Ok(result)
}

fn certify_minimality(form: &Form, prev: &Stage, weight: usize, column: &Vector) -> Result<()> {
    for (&index, polynomial) in &column.entries {
        let entry_degree = weight
            .checked_sub(prev.shifts[index])
            .filter(|&w| w > 0)
            .ok_or_else(|| {
                Error::plain("Verbitsky free resolution has a nonminimal constant entry")
            })?;
        for powers in polynomial.terms.keys() {
            if form.polynomial_weight(powers)? != entry_degree as u128 {
                return Err(Error::plain(
                    "Verbitsky free differential has inconsistent weight",
                ));
            }
        }
    }
    Ok(())
}

#[derive(Default)]
pub(super) struct Catalog {
    entries: BTreeMap<usize, Rc<Vec<Monomial>>>,
    cells: usize,
}

impl Catalog {
    pub(super) fn at(
        &mut self,
        model: &mut Model,
        form: &Form,
        weight: usize,
    ) -> Result<Rc<Vec<Monomial>>> {
        if let Some(monomials) = self.entries.get(&weight) {
            return Ok(Rc::clone(monomials));
        }
        let monomials = monomials(model, form, weight as u128)?;
        self.cells += monomials
            .iter()
            .map(|powers| powers.len() * 5 + 8)
            .sum::<usize>();
        storage(self.cells)?;
        let monomials = Rc::new(monomials);
        self.entries.insert(weight, Rc::clone(&monomials));
        Ok(monomials)
    }
}

pub(super) struct Axis(Vec<(usize, Monomial)>);

impl Axis {
    pub(super) fn new(
        model: &mut Model,
        form: &Form,
        catalog: &mut Catalog,
        shifts: &[usize],
        weight: usize,
    ) -> Result<Self> {
        let mut axis = Vec::new();
        let mut cells = 0;
        for (i, &shift) in shifts.iter().enumerate() {
            if let Some(remaining) = weight.checked_sub(shift) {
                let monomials = catalog.at(model, form, remaining)?;
                model.spend(monomials.len())?;
                for powers in monomials.iter() {
                    cells += 10 + 5 * powers.len();
                    storage(cells)?;
                    axis.push((i, powers.clone()));
                }
            }
        }
        Ok(Self(axis))
    }

    pub(super) fn entries(&self) -> &[(usize, Monomial)] {
        &self.0
    }

    pub(super) fn len(&self) -> usize {
        self.0.len()
    }

    pub(super) fn locate(&self, row: usize, monomial: &[super::Power]) -> Result<usize> {
        self.0
            .binary_search_by(|(i, powers)| {
                i.cmp(&row).then_with(|| powers.as_slice().cmp(monomial))
            })
            .map_err(|_| {
                Error::plain("Verbitsky free differential term lies outside its weighted module")
            })
    }

    pub(super) fn vector(&self, model: &mut Model, entries: &[Rational]) -> Result<Vector> {
        let mut result = Vector::default();
        for ((i, powers), &value) in self.0.iter().zip(entries) {
            result.insert(model, *i, powers.clone(), value)?;
        }
        Ok(result)
    }
}

fn visit_terms(
    model: &mut Model,
    column: &Vector,
    multiplier: &[super::Power],
    axis: &Axis,
    mut use_term: impl FnMut(usize, Rational),
) -> Result<()> {
    for (&row, polynomial) in &column.entries {
        model.spend(1 + polynomial.terms.len())?;
        for (powers, &value) in &polynomial.terms {
            let monomial = multiply_powers(powers, multiplier)?;
            use_term(axis.locate(row, &monomial)?, value);
        }
    }
    Ok(())
}

fn syzygies(
    model: &mut Model,
    form: &Form,
    catalog: &mut Catalog,
    source: &Stage,
    target: &Stage,
    known: &Stage,
    weight: usize,
) -> Result<Vec<Vector>> {
    let source_axis = Axis::new(model, form, catalog, &source.shifts, weight)?;
    if source_axis.len() == 0 {
        return Ok(Vec::new());
    }
    let target_axis = Axis::new(model, form, catalog, &target.shifts, weight)?;
    let mut outgoing = linear::matrix(target_axis.len(), source_axis.len())?;
    for (j, (i, powers)) in source_axis.0.iter().enumerate() {
        visit_terms(
            model,
            &source.columns[*i],
            powers,
            &target_axis,
            |row, value| outgoing[row][j] = value,
        )?;
    }
    let constraints = linear::rref(model, outgoing, source_axis.len())?;
    let mut kernel = constraints.kernel(source_axis.len())?;
    let lower_axis = Axis::new(model, form, catalog, &known.shifts, weight)?;
    let mut multiples = linear::matrix(lower_axis.len(), source_axis.len())?;
    for ((i, powers), row) in lower_axis.0.iter().zip(multiples.iter_mut()) {
        visit_terms(
            model,
            &known.columns[*i],
            powers,
            &source_axis,
            |index, value| row[index] = value,
        )?;
    }
    let lower = linear::rref(model, multiples, source_axis.len())?;
    for row in &mut kernel {
        for (&pivot, other) in lower.pivots.iter().zip(&lower.rows) {
            let coefficient = row[pivot].negate();
            linear::add_multiple(model, row, other, coefficient)?;
        }
    }
    let quotient = linear::rref(model, kernel, source_axis.len())?;
    if constraints.rows.len() + lower.rows.len() + quotient.rows.len() != source_axis.len() {
        return Err(Error::plain(
            "Verbitsky syzygy kernel and lower multiples have inconsistent ranks",
        ));
    }
    quotient
        .rows
        .iter()
        .map(|row| source_axis.vector(model, row))
        .collect()
}
