use std::cell::OnceCell;
use std::cmp::Ordering;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::may::Rows;
use crate::model::{Model, Value};
use crate::odd_ext;
use crate::prime_linear::{self, Basis, Coordinate};
use crate::sheaf::{self, Cochain, Cochains, Generator, Layer, Support, Term};
use crate::sheaf_hyper::Complex;

#[derive(Clone, Copy, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Grade {
    pub(crate) r: usize,
    pub(crate) p: usize,
    pub(crate) q: i128,
}

impl Grade {
    pub(crate) fn from_total(r: usize, p: usize, total: usize) -> Self {
        Self {
            r,
            p,
            q: total as i128 - p as i128,
        }
    }

    pub(crate) fn outgoing(self) -> Result<Self> {
        Ok(Self {
            r: self.r,
            p: self.p.checked_add(self.r).ok_or_else(grade_limit)?,
            q: self
                .q
                .checked_add(1 - self.r as i128)
                .ok_or_else(grade_limit)?,
        })
    }

    fn incoming(self) -> Option<Self> {
        Some(Self {
            r: self.r,
            p: self.p.checked_sub(self.r)?,
            q: self.q.checked_add(self.r as i128 - 1)?,
        })
    }

    pub(crate) fn next(self) -> Result<Self> {
        Ok(Self {
            r: self.r.checked_add(1).ok_or_else(page_limit)?,
            ..self
        })
    }

    pub(crate) fn previous(self) -> Result<Self> {
        Ok(Self {
            r: self
                .r
                .checked_sub(1)
                .ok_or_else(|| Error::plain("shEPrev requires a positive page index"))?,
            ..self
        })
    }

    fn cutoff(self) -> u128 {
        self.p as u128 + self.r as u128
    }

    pub(crate) fn total(self) -> Result<usize> {
        let total = self.q.checked_add(self.p as i128).ok_or_else(grade_limit)?;
        usize::try_from(total).map_err(|_| {
            Error::plain(
                "sheaf spectral sequence class has no representable nonnegative total degree",
            )
        })
    }

    pub(crate) fn has_leading(self, complex: &Complex) -> bool {
        usize::try_from(self.q).ok().is_some_and(|q| {
            complex
                .component_at(q)
                .is_some_and(|part| self.p < part.sheaf.height)
        })
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Differential,
    Representative,
    Next,
    Previous,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Differential => "shED",
            Self::Representative => "shERep",
            Self::Next => "shENext",
            Self::Previous => "shEPrev",
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Lift {
    Class,
    Correct,
}

impl Lift {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Class => "shEClass",
            Self::Correct => "shELiftTo",
        }
    }
}

pub(crate) struct Page {
    grade: Grade,
    layer: Layer,
    leading_start: usize,
    span: Basis,
    representatives: Vec<Box<[Coordinate]>>,
    transition: OnceCell<Basis>,
}

impl Page {
    pub(crate) fn build(model: &mut Model, complex: &Complex, grade: Grade) -> Result<Self> {
        let mut page = Self {
            grade,
            layer: Layer::default(),
            leading_start: 0,
            span: Basis::default(),
            representatives: Vec::new(),
            transition: OnceCell::new(),
        };
        if !grade.has_leading(complex) {
            return Ok(page);
        }
        let q = grade.q as usize;
        let total = grade.p.checked_add(q).ok_or_else(grade_limit)?;
        let high = if grade.r == 0 {
            grade.p
        } else {
            usize::try_from(grade.cutoff() - 1).unwrap_or(usize::MAX)
        };
        page.layer = complex.layer_window(model, total, grade.p, high, MAX_PIVOT_CELLS)?;
        page.leading_start = page.layer.vertical_start(q);
        if grade.r > 0 {
            page.quotient(model, complex, total)?;
        }
        model.allocate(5 + page.layer.cells())?;
        Ok(page)
    }

    pub(crate) fn rank(&self) -> usize {
        if self.grade.r == 0 {
            self.layer.size
        } else {
            self.representatives.len()
        }
    }

    fn leading_index(&self, model: &mut Model, term: &Generator) -> Result<usize> {
        Ok(self.layer.index(model, term)? - self.leading_start)
    }

    fn quotient(&mut self, model: &mut Model, complex: &Complex, total: usize) -> Result<()> {
        let prime = complex.prime();
        let mut span = self.boundaries(model, complex, total)?;
        let mut rows = Rows::new("sheaf hypercohomology");
        let mut cycles = Basis::default();
        let mut lift_cells = 0;
        for index in 0..self.layer.size {
            let generator = self.layer.generator(index);
            let column = below_column(model, complex, &generator, self.grade.cutoff(), &mut rows)?;
            // Higher filtration columns come first. Their kernel witnesses
            // have zero leading part but can correct a later lift.
            let lift = cycles.insert_with_witness(model, prime, column, vec![(index, 1)])?;
            sheaf::check_cells(lift_cells + span.cells() + cycles.cells() + self.layer.cells())?;
            let Some(lift) = lift else {
                continue;
            };
            if index < self.leading_start {
                continue;
            }
            let leading = lift
                .iter()
                .filter(|&&(i, _)| i >= self.leading_start)
                .map(|&(i, c)| (i - self.leading_start, c))
                .collect();
            let witness = vec![(self.representatives.len(), 1)];
            if span
                .insert_with_witness(model, prime, leading, witness)?
                .is_none()
            {
                lift_cells += 1 + 2 * lift.len();
                sheaf::check_cells(
                    lift_cells + span.cells() + cycles.cells() + self.layer.cells(),
                )?;
                self.representatives.push(lift.into_boxed_slice());
            }
        }
        model.allocate(lift_cells + span.cells())?;
        self.span = span;
        Ok(())
    }

    fn boundaries(&self, model: &mut Model, complex: &Complex, total: usize) -> Result<Basis> {
        let mut span = Basis::default();
        let Some(previous) = total.checked_sub(1) else {
            return Ok(span);
        };
        let p = self.grade.p;
        let low = p.saturating_sub(self.grade.r - 1);
        let layer = complex.layer_window(model, previous, low, p, MAX_PIVOT_CELLS)?;
        let prime = complex.prime();
        let mut lower_span = Basis::default();
        let mut rows = Rows::new("sheaf hypercohomology");
        for index in 0..layer.size {
            let mut lower = Vec::new();
            let mut leading = Vec::new();
            let generator = layer.generator(index);
            complex.cofaces(model, &generator, 1, &mut |model, term, coefficient| {
                match filtration(&term).cmp(&p) {
                    Ordering::Less => {
                        lower.push((rows.coordinate(model, cells(&term), term)?, coefficient));
                    }
                    Ordering::Equal => {
                        leading.push((self.leading_index(model, &term)?, coefficient));
                    }
                    Ordering::Greater => {}
                }
                sheaf::check_cells(2 * (lower.len() + leading.len()))
            })?;
            sort(model, &mut lower)?;
            sort(model, &mut leading)?;
            // A witness whose lower image cancels is a boundary that lies in F^p.
            if let Some(image) = lower_span.insert_with_witness(model, prime, lower, leading)? {
                span.insert(model, prime, image)?;
            }
            sheaf::check_cells(
                self.layer.cells() + layer.cells() + span.cells() + lower_span.cells(),
            )?;
        }
        Ok(span)
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        complex: &Complex,
        terms: &[Term],
    ) -> Result<Vec<Coordinate>> {
        let mut leading = Vec::new();
        for (term, coefficient) in terms {
            model.spend(1)?;
            match filtration(term).cmp(&self.grade.p) {
                Ordering::Less => return Err(filtration_error()),
                Ordering::Equal => leading.push((self.leading_index(model, term)?, *coefficient)),
                Ordering::Greater => {}
            }
        }
        if self.grade.r == 0 {
            return Ok(leading);
        }
        if !below_differential(model, complex, terms, self.grade.cutoff())?.is_empty() {
            return Err(Error::plain(
                "shEClass requires the sheaf page cycle condition",
            ));
        }
        let (rest, class) = self
            .span
            .reduce_with_witness(model, complex.prime(), leading)?;
        if !rest.is_empty() {
            return Err(Error::plain(
                "sheaf cycle is missing from the page quotient",
            ));
        }
        Ok(class)
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        complex: &Complex,
        class: &[Coordinate],
    ) -> Result<Vec<Term>> {
        let coordinates = if self.grade.r == 0 {
            class.to_vec()
        } else {
            let mut lift = Vec::new();
            for &(index, coefficient) in class {
                lift = prime_linear::add_scaled(
                    model,
                    complex.prime(),
                    &lift,
                    &self.representatives[index],
                    coefficient,
                )?;
            }
            lift
        };
        model.spend(coordinates.len())?;
        Ok(coordinates
            .into_iter()
            .map(|(i, coefficient)| (self.layer.generator(i), coefficient))
            .collect())
    }

    pub(crate) fn next(
        &self,
        model: &mut Model,
        complex: &Complex,
        class: &[Coordinate],
    ) -> Result<Value> {
        let transition = if let Some(transition) = self.transition.get() {
            transition
        } else {
            let data = self.build_transition(model, complex)?;
            self.transition.get_or_init(|| data)
        };
        let (rest, next) =
            transition.reduce_with_witness(model, complex.prime(), class.to_vec())?;
        if !rest.is_empty() {
            return Err(Error::plain(format!(
                "shENext requires a cycle for d_{}",
                self.grade.r
            )));
        }
        odd_ext::persist(model, next)
    }

    fn build_transition(&self, model: &mut Model, complex: &Complex) -> Result<Basis> {
        let mut cycles = Basis::default();
        let prime = complex.prime();
        if let Some(grade) = self.grade.incoming() {
            let incoming = complex.page(model, grade)?;
            for source in 0..incoming.rank() {
                let lift = incoming.representative(model, complex, &[(source, 1)])?;
                let image = complex.differential_support(model, &lift)?;
                let boundary = self.project(model, complex, &image)?;
                cycles.insert(model, prime, boundary)?;
            }
        }
        let next = complex.page(model, self.grade.next()?)?;
        for index in 0..next.rank() {
            let lift = next.representative(model, complex, &[(index, 1)])?;
            let earlier = self.project(model, complex, &lift)?;
            if cycles
                .insert_with_witness(model, prime, earlier, vec![(index, 1)])?
                .is_some()
            {
                return Err(Error::plain(
                    "sheaf page transition has dependent cycle lifts",
                ));
            }
        }
        model.allocate(cycles.cells())?;
        Ok(cycles)
    }
}

pub(crate) fn class(
    model: &mut Model,
    complex: &Complex,
    grade: Grade,
    raw: &[Term],
) -> Result<Value> {
    if raw.is_empty() {
        return Ok(odd_ext::zero());
    }
    let coordinates = complex.page(model, grade)?.project(model, complex, raw)?;
    odd_ext::persist(model, coordinates)
}

pub(crate) fn apply(
    model: &mut Model,
    complex: &Complex,
    grade: Grade,
    operation: Operation,
    coordinates: &[Coordinate],
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(match operation {
            Operation::Representative => Value::SheafC([].into()),
            _ => odd_ext::zero(),
        });
    }
    let source = complex.page(model, grade)?;
    match operation {
        Operation::Differential => {
            let lift = source.representative(model, complex, coordinates)?;
            let image = complex.differential_support(model, &lift)?;
            class(model, complex, grade.outgoing()?, &image)
        }
        Operation::Representative => {
            let lift = source.representative(model, complex, coordinates)?;
            sheaf::persist(model, lift)
        }
        Operation::Next => source.next(model, complex, coordinates),
        Operation::Previous => {
            let lift = source.representative(model, complex, coordinates)?;
            class(model, complex, grade.previous()?, &lift)
        }
    }
}

pub(crate) fn correct(
    model: &mut Model,
    complex: &Complex,
    grade: Grade,
    raw: &Support,
) -> Result<Value> {
    if raw.iter().any(|(term, _)| filtration(term) < grade.p) {
        return Err(filtration_error());
    }
    let defect = below_differential(model, complex, raw, grade.cutoff())?;
    if defect.is_empty() {
        return Ok(Value::SheafC(raw.clone()));
    }
    let prime = complex.prime();
    let correction = complex.layer_window(
        model,
        grade.total()?,
        grade.p.saturating_add(1),
        usize::try_from(grade.cutoff() - 1).unwrap_or(usize::MAX),
        MAX_PIVOT_CELLS,
    )?;
    let mut rows = Rows::new("sheaf hypercohomology correction");
    let mut span = Basis::default();
    for index in 0..correction.size {
        let generator = correction.generator(index);
        let column = below_column(model, complex, &generator, grade.cutoff(), &mut rows)?;
        span.insert_with_witness(model, prime, column, vec![(index, 1)])?;
        sheaf::check_cells(correction.cells() + span.cells())?;
    }
    let mut rhs = Vec::new();
    for (term, coefficient) in defect {
        rhs.push((
            rows.coordinate(model, cells(&term), term)?,
            prime.negate(coefficient),
        ));
    }
    sort(model, &mut rhs)?;
    let (rest, lift) = span.reduce_with_witness(model, prime, rhs)?;
    if !rest.is_empty() {
        return Err(Error::plain(format!(
            "shELiftTo cannot lift the prescribed leading part to page {}",
            grade.r
        )));
    }
    let lifted: Vec<_> = lift
        .into_iter()
        .map(|(index, coefficient)| (correction.generator(index), coefficient))
        .collect();
    sheaf::add(model, prime, raw, &lifted)
}

fn below_column(
    model: &mut Model,
    complex: &Complex,
    source: &Generator,
    cutoff: u128,
    rows: &mut Rows<Generator>,
) -> Result<Vec<Coordinate>> {
    let mut column = Vec::new();
    complex.cofaces(model, source, 1, &mut |model, term, coefficient| {
        if (filtration(&term) as u128) < cutoff {
            sheaf::check_cells(2 * (column.len() + 1))?;
            column.push((rows.coordinate(model, cells(&term), term)?, coefficient));
        }
        Ok(())
    })?;
    sort(model, &mut column)?;
    Ok(column)
}

fn below_differential(
    model: &mut Model,
    complex: &Complex,
    terms: &[Term],
    cutoff: u128,
) -> Result<Vec<Term>> {
    let mut output = Cochain::default();
    for (source, coefficient) in terms {
        if (filtration(source) as u128) < cutoff {
            complex.cofaces(model, source, *coefficient, &mut |model, term, scalar| {
                if (filtration(&term) as u128) < cutoff {
                    output.add(model, complex.prime(), term, scalar)?;
                }
                Ok(())
            })?;
        }
    }
    Ok(output.into_terms())
}

fn cells(generator: &Generator) -> usize {
    5 + generator.chain.len()
}

fn filtration(term: &Generator) -> usize {
    term.chain.len() - 1
}

fn sort(model: &mut Model, column: &mut [Coordinate]) -> Result<()> {
    model.spend(column.len() * (1 + column.len().max(1).ilog2() as usize))?;
    column.sort_unstable_by_key(|&(index, _)| index);
    Ok(())
}

fn grade_limit() -> Error {
    Error::plain("sheaf spectral sequence bidegree is out of range")
}

fn page_limit() -> Error {
    Error::plain("sheaf spectral sequence page is out of range")
}

fn filtration_error() -> Error {
    Error::plain("total cochain lies below the requested sheaf filtration")
}
