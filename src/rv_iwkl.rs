use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_iwahori::{self, Boundary};

const MAX_CELLS: usize = 300_000;

pub(crate) type Window = Rc<[i128]>;
type Polynomial = BTreeMap<i64, i128>;

pub(crate) fn polynomial(model: &mut Model, lower: &[i128], upper: &[i128]) -> Result<Vec<i128>> {
    if lower.len() != upper.len() {
        return Err(Error::plain("affine KL windows must have the same rank"));
    }
    rv_iwahori::validate_window(model, lower)?;
    rv_iwahori::validate_window(model, upper)?;
    if lower == upper {
        return Ok(vec![1]);
    }
    let component = rv_iwahori::component(lower)?;
    if component != rv_iwahori::component(upper)? {
        return Ok(vec![0]);
    }
    let lower = translate(model, lower, component)?;
    let upper = translate(model, upper, component)?;
    if !Boundary::Closure.contains(model, &lower, &upper)? {
        return Ok(vec![0]);
    }
    let upper_length = rv_iwahori::length(model, &upper)?;
    let difference = length_difference(model, &lower, upper_length)?;
    let basis = canonical_basis(model, &upper)?;
    let coefficient = basis
        .coefficients
        .get(&lower)
        .ok_or_else(|| Error::plain("comparable affine windows have no Hecke coefficient"))?;
    to_q(model, difference, coefficient)
}

fn length_difference(model: &mut Model, lower: &[i128], upper_length: u128) -> Result<i64> {
    let lower_length = rv_iwahori::length(model, lower)?;
    upper_length
        .checked_sub(lower_length)
        .and_then(|difference| i64::try_from(difference).ok())
        .ok_or_else(arithmetic_limit)
}

pub(crate) fn expansion(model: &mut Model, upper: &[i128]) -> Result<Vec<(Window, Vec<i128>)>> {
    rv_iwahori::validate_window(model, upper)?;
    let component = rv_iwahori::component(upper)?;
    let upper = translate(model, upper, component)?;
    let basis = canonical_basis(model, &upper)?;
    let upper_length = rv_iwahori::length(model, &upper)?;
    let mut output = Vec::new();
    let mut cells = 0;
    for (window, polynomial) in basis.coefficients {
        let difference = length_difference(model, &window, upper_length)?;
        let coefficient = to_q(model, difference, &polynomial)?;
        cells += window_cells(&window) + 2 * coefficient.len();
        check_storage(cells)?;
        let window = window
            .iter()
            .map(|&entry| entry.checked_add(component).ok_or_else(arithmetic_limit))
            .collect::<Result<Vec<_>>>()?;
        output.push((window.into(), coefficient));
    }
    Ok(output)
}

fn to_q(model: &mut Model, difference: i64, coefficient: &Polynomial) -> Result<Vec<i128>> {
    let mut output = Vec::new();
    for (&power, &value) in coefficient {
        let doubled_degree = difference.checked_sub(power).ok_or_else(arithmetic_limit)?;
        if doubled_degree < 0 || doubled_degree % 2 != 0 || value < 0 {
            return Err(Error::plain("invalid affine Kazhdan–Lusztig coefficient"));
        }
        let degree = usize::try_from(doubled_degree / 2).map_err(|_| arithmetic_limit())?;
        let cells = degree
            .checked_add(1)
            .and_then(|n| n.checked_mul(2))
            .ok_or_else(arithmetic_limit)?;
        check_storage(cells)?;
        if output.len() <= degree {
            model.spend(1 + degree - output.len())?;
            output.resize(degree + 1, 0);
        }
        output[degree] = value;
    }
    Ok(output)
}

pub(crate) fn translate(model: &mut Model, window: &[i128], component: i128) -> Result<Window> {
    model.spend(window.len())?;
    check_storage(2 * window.len())?;
    let result = window
        .iter()
        .map(|&entry| entry.checked_sub(component).ok_or_else(arithmetic_limit))
        .collect::<Result<Vec<_>>>()?;
    Ok(result.into())
}

fn arithmetic_limit() -> Error {
    Error::plain("affine Kazhdan–Lusztig arithmetic is out of range")
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain(
            "affine Kazhdan–Lusztig calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

pub(crate) fn descent(window: &[i128], simple: usize) -> Result<bool> {
    if let Some(&next) = window.get(simple + 1) {
        Ok(window[simple] > next)
    } else {
        let next = window[0]
            .checked_add(window.len() as i128)
            .ok_or_else(arithmetic_limit)?;
        Ok(window[simple] > next)
    }
}

pub(crate) fn right_simple(model: &mut Model, window: &[i128], simple: usize) -> Result<Window> {
    model.spend(window.len())?;
    let mut result = window.to_vec();
    if simple + 1 < window.len() {
        result.swap(simple, simple + 1);
    } else {
        let rank = window.len() as i128;
        result[0] = window[simple]
            .checked_sub(rank)
            .ok_or_else(arithmetic_limit)?;
        result[simple] = window[0].checked_add(rank).ok_or_else(arithmetic_limit)?;
    }
    Ok(result.into())
}

pub(crate) fn wall(model: &mut Model, window: &[i128]) -> Result<Option<usize>> {
    if window.len() > 1 {
        for simple in 0..window.len() {
            model.spend(1)?;
            if descent(window, simple)? {
                return Ok(Some(simple));
            }
        }
    }
    Ok(None)
}

fn window_cells(window: &[i128]) -> usize {
    1 + 2 * window.len()
}

#[derive(Default)]
struct Expansion {
    coefficients: BTreeMap<Window, Polynomial>,
    cells: usize,
}

impl Expansion {
    fn add(
        &mut self,
        model: &mut Model,
        window: &Window,
        polynomial: &Polynomial,
        shift: i64,
        multiple: i128,
        retained: usize,
    ) -> Result<()> {
        let width = 1 + self.coefficients.len().max(1).ilog2() as usize;
        model.spend(width * window_cells(window) + polynomial.len())?;
        if !self.coefficients.contains_key(window) {
            self.cells += window_cells(window);
            check_storage(retained + self.cells)?;
        }
        let coefficient = self.coefficients.entry(window.clone()).or_default();
        for (&degree, &value) in polynomial {
            let degree = degree.checked_add(shift).ok_or_else(arithmetic_limit)?;
            let value = value.checked_mul(multiple).ok_or_else(arithmetic_limit)?;
            if let Some(previous) = coefficient.get_mut(&degree) {
                *previous = previous.checked_add(value).ok_or_else(arithmetic_limit)?;
                if *previous == 0 {
                    coefficient.remove(&degree);
                    self.cells -= 3;
                }
            } else {
                self.cells += 3;
                check_storage(retained + self.cells)?;
                coefficient.insert(degree, value);
            }
        }
        if coefficient.is_empty() {
            self.coefficients.remove(window);
            self.cells -= window_cells(window);
        }
        Ok(())
    }
}

enum Job {
    Build(Window),
    Prepare(Window, Window, usize),
    Finish(Window, Window, usize, Vec<(Window, i128)>),
}

impl Job {
    fn cells(&self) -> usize {
        match self {
            Self::Build(window) => window_cells(window),
            Self::Prepare(window, parent, _) => window_cells(window) + window_cells(parent) + 1,
            Self::Finish(window, parent, _, corrections) => {
                window_cells(window)
                    + window_cells(parent)
                    + 1
                    + corrections
                        .iter()
                        .map(|(term, _)| window_cells(term) + 2)
                        .sum::<usize>()
            }
        }
    }
}

#[derive(Default)]
struct Jobs {
    pending: Vec<(Job, usize)>,
    cells: usize,
}

impl Jobs {
    fn push(&mut self, job: Job, retained: usize) -> Result<()> {
        let cells = job.cells();
        check_storage(retained + self.cells + cells)?;
        self.pending.push((job, cells));
        self.cells += cells;
        Ok(())
    }

    fn pop(&mut self) -> Option<Job> {
        let (job, cells) = self.pending.pop()?;
        self.cells -= cells;
        Some(job)
    }
}

#[derive(Default)]
struct Bases {
    elements: BTreeMap<Window, Expansion>,
    cells: usize,
}

impl Bases {
    fn insert(&mut self, window: Window, expansion: Expansion, retained: usize) -> Result<()> {
        self.cells += window_cells(&window) + expansion.cells;
        check_storage(self.cells + retained)?;
        self.elements.insert(window, expansion);
        Ok(())
    }

    fn product(
        &self,
        model: &mut Model,
        parent: &Window,
        simple: usize,
        corrections: &[(Window, i128)],
        retained: usize,
    ) -> Result<Expansion> {
        let mut result = Expansion::default();
        for (term, polynomial) in &self.elements[parent].coefficients {
            let reflected = right_simple(model, term, simple)?;
            result.add(model, &reflected, polynomial, 0, 1, retained)?;
            let shift = if descent(term, simple)? { -1 } else { 1 };
            result.add(model, term, polynomial, shift, 1, retained)?;
        }
        for (term, multiplicity) in corrections {
            for (lower, polynomial) in &self.elements[term].coefficients {
                result.add(
                    model,
                    lower,
                    polynomial,
                    0,
                    multiplicity.checked_neg().ok_or_else(arithmetic_limit)?,
                    retained,
                )?;
            }
        }
        Ok(result)
    }
}

fn canonical_basis(model: &mut Model, upper: &Window) -> Result<Expansion> {
    let mut bases = Bases::default();
    let mut jobs = Jobs::default();
    jobs.push(Job::Build(upper.clone()), 0)?;
    while let Some(job) = jobs.pop() {
        model.spend(1)?;
        let active_cells = job.cells();
        match job {
            Job::Build(window) => {
                if bases.elements.contains_key(&window) {
                    continue;
                }
                if let Some(simple) = wall(model, &window)? {
                    let parent = right_simple(model, &window, simple)?;
                    jobs.push(Job::Prepare(window, parent.clone(), simple), bases.cells)?;
                    if !bases.elements.contains_key(&parent) {
                        jobs.push(Job::Build(parent), bases.cells)?;
                    }
                } else {
                    let mut identity = Expansion::default();
                    identity.add(
                        model,
                        &window,
                        &BTreeMap::from([(0, 1)]),
                        0,
                        1,
                        bases.cells + jobs.cells,
                    )?;
                    bases.insert(window, identity, jobs.cells)?;
                }
            }
            Job::Prepare(window, parent, simple) => {
                let mut corrections = Vec::new();
                for (term, polynomial) in &bases.elements[&parent].coefficients {
                    model.spend(1)?;
                    if let Some(&coefficient) = polynomial.get(&1) {
                        if descent(term, simple)? {
                            corrections.push((term.clone(), coefficient));
                        }
                    }
                }
                let dependencies = corrections
                    .iter()
                    .filter(|(term, _)| !bases.elements.contains_key(term))
                    .map(|(term, _)| term.clone())
                    .collect::<Vec<_>>();
                jobs.push(
                    Job::Finish(window, parent, simple, corrections),
                    bases.cells,
                )?;
                for term in dependencies.into_iter().rev() {
                    jobs.push(Job::Build(term), bases.cells)?;
                }
            }
            Job::Finish(window, parent, simple, corrections) => {
                let expansion = bases.product(
                    model,
                    &parent,
                    simple,
                    &corrections,
                    bases.cells + jobs.cells + active_cells,
                )?;
                bases.insert(window, expansion, jobs.cells + active_cells)?;
            }
        }
    }
    Ok(bases.elements.remove(upper).unwrap())
}
