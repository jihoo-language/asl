use std::ops::RangeInclusive;
use std::rc::Rc;

use crate::cohomology;
use crate::error::{Error, Result};
use crate::may::{FilteredComplex, FilteredPage};
use crate::model::{Kind, Model, TypeId, Value};
use crate::motivic;
use crate::page_transition::Transition;
use crate::real_motivic::{self, Term};

const MAX_BASIS: usize = 32_768;

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct Grade {
    pub(crate) chain: motivic::Grade,
    pub(crate) r: usize,
    pub(crate) f: usize,
}

impl Grade {
    pub(crate) fn differential(self) -> Result<Self> {
        let chain = self.chain.differential()?;
        let f = self
            .f
            .checked_add(self.r)
            .ok_or_else(|| Error::plain("rhoD filtration is out of range"))?;
        Ok(Self { chain, f, ..self })
    }

    pub(crate) fn next(self) -> Result<Self> {
        let r = self
            .r
            .checked_add(1)
            .ok_or_else(|| Error::plain("rhoNext page is out of range"))?;
        Ok(Self { r, ..self })
    }

    pub(crate) fn previous(self) -> Result<Self> {
        let r = self
            .r
            .checked_sub(1)
            .ok_or_else(|| Error::plain("rhoPrev requires a positive page number"))?;
        Ok(Self { r, ..self })
    }

    pub(crate) fn product(self, other: Self) -> Result<Self> {
        if self.r != other.r {
            return Err(Error::plain(
                "rhoMul requires RhoBock classes on the same page",
            ));
        }
        let f = self
            .f
            .checked_add(other.f)
            .ok_or_else(|| Error::plain("rhoMul filtration is out of range"))?;
        Ok(Self {
            chain: self.chain.product(other.chain)?,
            f,
            ..self
        })
    }

    fn incoming(self) -> Option<Self> {
        let chain = self.chain.previous()?;
        let f = self.f.checked_sub(self.r)?;
        Some(Self { chain, f, ..self })
    }

    fn cutoff(self) -> u128 {
        self.f as u128 + self.r.max(1) as u128
    }
}

#[derive(Clone, Copy)]
struct Complex(motivic::Grade);

impl FilteredComplex for Complex {
    type Term = Term;
    const NAME: &'static str = "rho-Bockstein";
    const ELEMENT_NAME: &'static str = "real-motivic";

    fn basis(&mut self, model: &mut Model, s: usize, f: usize) -> Result<Vec<Term>> {
        self.basis_window(model, s, f..=usize::MAX)
    }

    fn basis_window(
        &mut self,
        model: &mut Model,
        s: usize,
        window: RangeInclusive<usize>,
    ) -> Result<Vec<Term>> {
        real_motivic::basis_window(model, motivic::Grade { s, ..self.0 }, window, MAX_BASIS)
    }

    fn differential(&mut self, model: &mut Model, term: &Term) -> Result<Vec<Term>> {
        self.differential_support(model, std::slice::from_ref(term))
    }

    fn differential_support(&mut self, model: &mut Model, terms: &[Term]) -> Result<Vec<Term>> {
        real_motivic::differential_support(model, terms)
    }

    fn differential_below(
        &mut self,
        model: &mut Model,
        terms: &[Term],
        cutoff: u128,
    ) -> Result<Vec<Term>> {
        real_motivic::differential_below(model, terms, cutoff)
    }

    fn filtration(&self, term: &Term) -> Result<usize> {
        Ok(term.rho)
    }

    fn cells(term: &Term) -> usize {
        term.cells()
    }
}

pub(crate) struct Page {
    grade: Grade,
    filtered: FilteredPage<Complex>,
}

impl Page {
    pub(crate) fn build(model: &mut Model, grade: Grade) -> Result<Self> {
        let mut complex = Complex(grade.chain);
        let filtered = FilteredPage::build(model, &mut complex, grade.r, grade.chain.s, grade.f)?;
        Ok(Self { grade, filtered })
    }

    pub(crate) fn rank(&self) -> usize {
        self.filtered.rank()
    }

    fn project(&self, model: &mut Model, terms: &[Term]) -> Result<Value> {
        let coords = self
            .filtered
            .project(model, &mut Complex(self.grade.chain), terms)?;
        model.allocate(coords.len())?;
        Ok(Value::Class(coords.into()))
    }

    fn representative(&self, model: &mut Model, coords: &[usize]) -> Result<Vec<Term>> {
        self.filtered.representative(model, coords)
    }

    fn project_differential(&self, model: &mut Model, lift: &[Term]) -> Result<Value> {
        let image = real_motivic::differential_below(model, lift, self.grade.cutoff())?;
        self.project(model, &image)
    }
}

pub(crate) fn class(model: &mut Model, output: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(cohomology::zero());
    }
    model.rho_page(output)?.project(model, terms)
}

pub(crate) fn representative(model: &mut Model, input: TypeId, coords: &[usize]) -> Result<Value> {
    if coords.is_empty() {
        return Ok(real_motivic::zero());
    }
    let terms = model.rho_page(input)?.representative(model, coords)?;
    real_motivic::persist(model, terms)
}

pub(crate) fn frobenius(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    coords: &[usize],
) -> Result<Value> {
    if coords.is_empty() {
        return Ok(cohomology::zero());
    }
    let Kind::RhoBock(grade) = *model.kind(output) else {
        unreachable!();
    };
    let lift = model.rho_page(input)?.representative(model, coords)?;
    let image = real_motivic::frobenius_below(model, &lift, grade.cutoff())?;
    class(model, output, &image)
}

pub(crate) fn differential(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    coords: &[usize],
) -> Result<Value> {
    if coords.is_empty() {
        return Ok(cohomology::zero());
    }
    let lift = model.rho_page(input)?.representative(model, coords)?;
    model.rho_page(output)?.project_differential(model, &lift)
}

pub(crate) fn previous(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    coords: &[usize],
) -> Result<Value> {
    if coords.is_empty() {
        return Ok(cohomology::zero());
    }
    let lift = model.rho_page(input)?.representative(model, coords)?;
    class(model, output, &lift)
}

pub(crate) fn product(
    model: &mut Model,
    left: (TypeId, &[usize]),
    right: (TypeId, &[usize]),
    output: TypeId,
) -> Result<Value> {
    if left.1.is_empty() || right.1.is_empty() {
        return Ok(cohomology::zero());
    }
    let a = model.rho_page(left.0)?.representative(model, left.1)?;
    let b = model.rho_page(right.0)?.representative(model, right.1)?;
    let terms = real_motivic::concat_support(model, &a, &b)?;
    class(model, output, &terms)
}

pub(crate) fn transition(model: &mut Model, ty: TypeId) -> Result<Transition> {
    let Kind::RhoBock(grade) = *model.kind(ty) else {
        unreachable!();
    };
    let source = model.rho_page(ty)?;
    let next_ty = model.rho_page_type(grade.next()?)?;
    let next = model.rho_page(next_ty)?;
    let incoming: Option<Rc<Page>> = grade
        .incoming()
        .map(|incoming| {
            let ty = model.rho_page_type(incoming)?;
            model.rho_page(ty)
        })
        .transpose()?;
    Transition::from_columns(
        model,
        grade.r,
        "rhoNext",
        incoming.as_ref().map_or(0, |page| page.rank()),
        next.rank(),
        |model, index| {
            let lift = incoming.as_ref().unwrap().representative(model, &[index])?;
            let image = real_motivic::differential_below(model, &lift, grade.cutoff())?;
            source
                .filtered
                .project(model, &mut Complex(grade.chain), &image)
        },
        |model, index| {
            let lift = next.representative(model, &[index])?;
            source
                .filtered
                .project(model, &mut Complex(grade.chain), &lift)
        },
    )
}
