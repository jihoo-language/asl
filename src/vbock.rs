use std::rc::Rc;

use crate::error::{Error, Result};
use crate::may::{FilteredComplex, FilteredPage};
use crate::model::{Kind, Model, TypeId, Value};
use crate::page_transition::Transition;
use crate::q_koszul::{self, Complex, Support, Term};
use crate::v_koszul;

const MAX_BASIS_TERMS: usize = 32_768;

#[derive(Clone, PartialEq, Eq, Hash)]
pub(crate) struct Grade {
    pub(crate) localized: v_koszul::Grade,
    pub(crate) r: usize,
    pub(crate) f: usize,
}

pub(crate) struct BocksteinPage {
    grade: v_koszul::Grade,
    length: i128,
    page: FilteredPage<LocalizedComplex>,
}

struct LocalizedComplex {
    grade: v_koszul::Grade,
    length: i128,
    differential: Complex,
}

impl LocalizedComplex {
    fn new(model: &mut Model, grade: &v_koszul::Grade, length: i128) -> Result<Self> {
        Ok(Self {
            grade: grade.clone(),
            length,
            differential: Complex::new(model, grade.source, grade.target, &grade.selected)?,
        })
    }
}

impl FilteredComplex for LocalizedComplex {
    type Term = Term;
    const NAME: &'static str = "v-Bockstein";

    fn basis(&mut self, model: &mut Model, s: usize, f: usize) -> Result<Vec<Term>> {
        // The filtered reducer only requests its current and preceding
        // lengths. Use 1 for the signed length and 0 for its predecessor.
        let s = self.length + s as i128 - 1;
        let basis = q_koszul::basis_localized(
            model,
            self.grade.source,
            self.grade.target,
            &self.grade.selected,
            s,
            self.grade.t,
            MAX_BASIS_TERMS,
        )?;
        let mut terms = Vec::new();
        for term in basis {
            model.spend(term.cells())?;
            if self.filtration(&term)? >= f {
                terms.push(term);
            }
        }
        Ok(terms)
    }

    fn differential(&mut self, model: &mut Model, term: &Term) -> Result<Vec<Term>> {
        self.differential_support(model, std::slice::from_ref(term))
    }

    fn differential_support(&mut self, model: &mut Model, terms: &[Term]) -> Result<Vec<Term>> {
        self.differential.differential(model, terms)
    }

    fn filtration(&self, term: &Term) -> Result<usize> {
        // Only the last variable is inverted. Its exponent can have either
        // sign and does not contribute to the lower-variable filtration.
        let lower = &term.powers[..term.powers.len() - 1];
        lower.iter().try_fold(0usize, |sum, &power| {
            let add = usize::try_from(power).map_err(|_| filtration_overflow())?;
            sum.checked_add(add).ok_or_else(filtration_overflow)
        })
    }

    fn cells(term: &Term) -> usize {
        term.cells()
    }
}

impl BocksteinPage {
    pub(crate) fn build(model: &mut Model, grade: &Grade) -> Result<Self> {
        Self::build_wide(model, grade, i128::from(grade.localized.s))
    }

    fn build_wide(model: &mut Model, grade: &Grade, length: i128) -> Result<Self> {
        let mut complex = LocalizedComplex::new(model, &grade.localized, length)?;
        let page = FilteredPage::build(model, &mut complex, grade.r, 1, grade.f)?;
        Ok(Self {
            grade: grade.localized.clone(),
            length,
            page,
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.page.rank()
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<Term>> {
        self.page.representative(model, class)
    }

    pub(crate) fn project(&self, model: &mut Model, terms: &[Term]) -> Result<Vec<usize>> {
        let mut complex = LocalizedComplex::new(model, &self.grade, self.length)?;
        self.page.project(model, &mut complex, terms)
    }
}

fn filtration_overflow() -> Error {
    Error::plain("v-Bockstein filtration is out of range")
}

pub(crate) fn differential_type(model: &mut Model, ty: TypeId) -> Result<TypeId> {
    let Kind::VBock(grade) = model.kind(ty) else {
        unreachable!()
    };
    let mut grade = grade.clone();
    grade.localized.s = grade
        .localized
        .s
        .checked_add(1)
        .ok_or_else(|| Error::plain("vbD cohomological degree is out of range"))?;
    grade.f = grade
        .f
        .checked_add(grade.r)
        .ok_or_else(filtration_overflow)?;
    model.vbock_type(grade)
}

pub(crate) fn composition_type(model: &mut Model, outer: TypeId, inner: TypeId) -> Result<TypeId> {
    let (Kind::VBock(outer), Kind::VBock(inner)) = (model.kind(outer), model.kind(inner)) else {
        return Err(Error::plain("ecomp requires two VBock classes"));
    };
    if outer.r != inner.r {
        return Err(Error::plain(
            "ecomp requires VBock classes on the same page",
        ));
    }
    let localized = v_koszul::composed(&outer.localized, &inner.localized)?;
    let f = outer
        .f
        .checked_add(inner.f)
        .ok_or_else(filtration_overflow)?;
    model.vbock_type(Grade {
        localized,
        r: outer.r,
        f,
    })
}

pub(crate) fn shift_type(model: &mut Model, ty: TypeId, k: i64) -> Result<TypeId> {
    let Kind::VBock(grade) = model.kind(ty) else {
        unreachable!()
    };
    let mut result = grade.clone();
    result.localized = v_koszul::shifted(&grade.localized, k)?;
    model.vbock_type(result)
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(zero());
    }
    let coordinates = model.vbock_page(ty)?.project(model, terms)?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

pub(crate) fn representative(model: &mut Model, ty: TypeId, class: &[usize]) -> Result<Support> {
    if class.is_empty() {
        return Ok([].into());
    }
    let terms = model.vbock_page(ty)?.representative(model, class)?;
    q_koszul::persist(model, terms)
}

pub(crate) fn differential(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    coordinates: &[usize],
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(zero());
    }
    let page = model.vbock_page(input)?;
    let lift = page.representative(model, coordinates)?;
    let mut complex = LocalizedComplex::new(model, &page.grade, page.length)?;
    let image = complex.differential_support(model, &lift)?;
    class(model, output, &image)
}

pub(crate) fn compose(
    model: &mut Model,
    (outer_ty, outer): (TypeId, &[usize]),
    (inner_ty, inner): (TypeId, &[usize]),
    output: TypeId,
) -> Result<Value> {
    if outer.is_empty() || inner.is_empty() {
        return Ok(zero());
    }
    let outer = model.vbock_page(outer_ty)?.representative(model, outer)?;
    let inner = model.vbock_page(inner_ty)?.representative(model, inner)?;
    let terms = q_koszul::compose_terms(model, &outer, &inner)?;
    class(model, output, &terms)
}

pub(crate) fn shift(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    k: i64,
    coordinates: &[usize],
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(zero());
    }
    let terms = model
        .vbock_page(input)?
        .representative(model, coordinates)?;
    let shifted = v_koszul::shift_terms(model, &terms, k)?;
    class(model, output, &shifted)
}

pub(crate) fn previous(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    coordinates: &[usize],
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(zero());
    }
    let lift = model
        .vbock_page(input)?
        .representative(model, coordinates)?;
    class(model, output, &lift)
}

pub(crate) fn transition(model: &mut Model, ty: TypeId) -> Result<Transition> {
    let Kind::VBock(grade) = model.kind(ty) else {
        unreachable!()
    };
    let grade = grade.clone();
    let source = model.vbock_page(ty)?;
    let incoming = incoming_page(model, &grade)?;
    let next_r = grade
        .r
        .checked_add(1)
        .ok_or_else(|| Error::plain("vbNext page is out of range"))?;
    let next = model.vbock_type(Grade {
        r: next_r,
        ..grade.clone()
    })?;
    let next = model.vbock_page(next)?;
    let mut complex = LocalizedComplex::new(model, &grade.localized, source.length)?;
    Transition::from_columns(
        model,
        grade.r,
        "vbNext",
        incoming.as_ref().map_or(0, |page| page.rank()),
        next.rank(),
        |model, index| {
            let terms = incoming.as_ref().unwrap().representative(model, &[index])?;
            let image = complex.differential_support(model, &terms)?;
            source.project(model, &image)
        },
        |model, index| {
            let terms = next.representative(model, &[index])?;
            source.project(model, &terms)
        },
    )
}

fn incoming_page(model: &mut Model, grade: &Grade) -> Result<Option<Rc<BocksteinPage>>> {
    let Some(f) = grade.f.checked_sub(grade.r) else {
        return Ok(None);
    };
    let mut incoming = grade.clone();
    incoming.f = f;
    if let Some(s) = incoming.localized.s.checked_sub(1) {
        incoming.localized.s = s;
        let ty = model.vbock_type(incoming)?;
        Ok(Some(model.vbock_page(ty)?))
    } else {
        // Internal elimination still sees an incoming degree below the
        // signed surface type range.
        let s = i128::from(incoming.localized.s) - 1;
        Ok(Some(Rc::new(BocksteinPage::build_wide(
            model, &incoming, s,
        )?)))
    }
}

fn zero() -> Value {
    Value::Class([].into())
}
