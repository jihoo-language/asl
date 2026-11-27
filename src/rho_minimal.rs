use std::ops::RangeInclusive;
use std::rc::Rc;

use crate::cobar::{self, Support as Sparse, Word};
use crate::complex::xor_sorted;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::motivic::{self, Grade};
use crate::real_motivic::{self, Term};

pub(crate) type Coordinate = (usize, usize);
pub(crate) type Support = Rc<[Coordinate]>;
type Cochain = Vec<Term>;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Inject(usize),
    Coefficient(usize),
    Differential,
    Lift,
    Project,
    Homotopy,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Inject(_) => "rhoMin",
            Self::Coefficient(_) => "rhoMinCoef",
            Self::Differential => "rhoMinD",
            Self::Lift => "rhoMinI",
            Self::Project => "rhoMinP",
            Self::Homotopy => "rhoMinH",
        }
    }
}

pub(crate) fn coefficient_grade(grade: Grade, f: usize) -> Result<Grade> {
    let t = grade.t.checked_add(f as i128).ok_or_else(grade_limit)?;
    let w = i64::try_from(i128::from(grade.w) + f as i128).map_err(|_| grade_limit())?;
    Ok(Grade { t, w, ..grade })
}

fn grade_limit() -> Error {
    Error::plain("rho-minimal coefficient grade exceeds the supported range")
}

pub(crate) fn zero() -> Value {
    Value::RhoMin([].into())
}

fn persist(model: &mut Model, support: Vec<Coordinate>) -> Result<Value> {
    check_cells(support.len())?;
    model.allocate(2 * support.len())?;
    Ok(Value::RhoMin(support.into()))
}

fn check_cells(terms: usize) -> Result<()> {
    cobar::check_storage(terms.saturating_mul(2))
}

pub(crate) fn inject(model: &mut Model, exponent: usize, coords: &[usize]) -> Result<Value> {
    model.spend(coords.len())?;
    persist(model, coords.iter().map(|&i| (exponent, i)).collect())
}

pub(crate) fn coefficient(model: &mut Model, support: &[Coordinate], f: usize) -> Result<Value> {
    model.spend(2 * (1 + support.len().max(1).ilog2() as usize))?;
    let start = support.partition_point(|&(e, _)| e < f);
    let end = support.partition_point(|&(e, _)| e <= f);
    model.spend(end - start)?;
    model.allocate(end - start)?;
    Ok(Value::Class(
        support[start..end].iter().map(|&(_, i)| i).collect(),
    ))
}

pub(crate) fn add(model: &mut Model, a: &[Coordinate], b: &[Coordinate]) -> Result<Value> {
    model.spend(a.len() + b.len())?;
    persist(model, xor_sorted(a, b))
}

pub(crate) fn shift(model: &mut Model, support: &Support, f: usize) -> Result<Value> {
    if f == 0 {
        return Ok(Value::RhoMin(support.clone()));
    }
    model.spend(support.len())?;
    let terms = support
        .iter()
        .map(|&(e, i)| {
            let e = e
                .checked_add(f)
                .ok_or_else(|| Error::plain("rho-minimal exponent is out of range"))?;
            Ok((e, i))
        })
        .collect::<Result<Vec<_>>>()?;
    persist(model, terms)
}

pub(crate) fn basis(model: &mut Model, grade: Grade, limit: usize) -> Result<Vec<Coordinate>> {
    let Some(window) = filtration_window(grade)? else {
        return Ok(Vec::new());
    };
    let mut terms = Vec::new();
    for f in window {
        model.spend(1)?;
        let coefficient = coefficient_grade(grade, f)?;
        let rank = if grade.s == 1 {
            usize::from(primitive_degree(coefficient)?.is_some())
        } else {
            let ty = model.mext_type(coefficient)?;
            model.adams(ty)?.rank()
        };
        if rank > limit - terms.len() {
            return Err(Error::plain(
                "rho-minimal basis exceeds the finite enumeration limit",
            ));
        }
        terms.extend((0..rank).map(|i| (f, i)));
    }
    model.allocate(2 * terms.len())?;
    Ok(terms)
}

fn filtration_window(grade: Grade) -> Result<Option<RangeInclusive<usize>>> {
    let w = i128::from(grade.w);
    if grade.t < w || (grade.s == 0 && grade.t > 0) {
        return Ok(None);
    }
    if grade.s == 0 {
        let f = usize::try_from(-grade.t).map_err(|_| grade_limit())?;
        return Ok(Some(f..=f));
    }
    let lower = (grade.s as i128)
        .checked_sub(grade.t)
        .ok_or_else(grade_limit)?
        .max(0);
    let upper = grade.t.checked_sub(2 * w).ok_or_else(grade_limit)?;
    if upper < lower {
        return Ok(None);
    }
    let lower = usize::try_from(lower).map_err(|_| grade_limit())?;
    let upper = usize::try_from(upper).map_err(|_| grade_limit())?;
    Ok(Some(lower..=upper))
}

pub(crate) fn lift(model: &mut Model, grade: Grade, support: &[Coordinate]) -> Result<Value> {
    let terms = match scalar_boundary(model, grade, support)? {
        Some(boundary) => boundary,
        None => lifted(model, grade, support)?,
    };
    real_motivic::persist(model, terms)
}

pub(crate) fn project(model: &mut Model, grade: Grade, source: &[Term]) -> Result<Value> {
    let corrected = correct_for_projection(model, grade, source)?;
    let terms = base_project(model, grade, &corrected)?;
    persist(model, terms)
}

pub(crate) fn homotopy(model: &mut Model, grade: Grade, source: &[Term]) -> Result<Value> {
    if grade.s <= 1 || source.is_empty() {
        return Ok(real_motivic::zero());
    }
    let corrected = correct_for_projection(model, grade, source)?;
    let terms = base_homotopy(model, grade, &corrected)?;
    real_motivic::persist(model, terms)
}

pub(crate) fn differential(
    model: &mut Model,
    grade: Grade,
    support: &[Coordinate],
) -> Result<Value> {
    if scalar_boundary(model, grade, support)?.is_some() {
        return Ok(zero());
    }
    let lifted = lifted(model, grade, support)?;
    let image = perturbation(model, &lifted)?;
    let terms = base_project(model, grade.differential()?, &image)?;
    persist(model, terms)
}

fn scalar_boundary(
    model: &mut Model,
    grade: Grade,
    support: &[Coordinate],
) -> Result<Option<Cochain>> {
    let [(exponent, 0)] = support else {
        return Ok(None);
    };
    if grade.s != 1 || grade.t > 0 {
        return Ok(None);
    }
    let rho = grade.t.checked_neg().and_then(|r| usize::try_from(r).ok());
    let tau = grade
        .t
        .checked_sub(i128::from(grade.w))
        .and_then(|t| usize::try_from(t).ok());
    let (Some(rho), Some(tau)) = (rho, tau) else {
        return Ok(None);
    };
    if !tau.is_power_of_two() || rho.checked_add(tau) != Some(*exponent) {
        return Ok(None);
    }

    // In scalar degree hQ is zero, so the scalar is already its own
    // perturbed lift and its minimal image is p(d scalar). Checking that
    // image identifies the canonical lift as d scalar by I Delta = d I.
    let scalar = Term {
        rho,
        tau,
        word: [].into(),
    };
    let boundary = real_motivic::differential_support(model, &[scalar])?;
    let projected = base_project(model, grade, &boundary)?;
    Ok((projected == support).then_some(boundary))
}

fn lifted(model: &mut Model, grade: Grade, support: &[Coordinate]) -> Result<Cochain> {
    let representative = base_lift(model, grade, support)?;
    if grade.s == 0 || representative.is_empty() {
        return Ok(representative);
    }
    let next = grade.differential()?;
    geometric_sum(model, representative, |model, terms| {
        let image = perturbation(model, terms)?;
        base_homotopy(model, next, &image)
    })
}

fn correct_for_projection(model: &mut Model, grade: Grade, source: &[Term]) -> Result<Cochain> {
    model.spend(source.iter().map(Term::cells).sum())?;
    let source = source.to_vec();
    if grade.s <= 1 {
        return Ok(source);
    }
    geometric_sum(model, source, |model, terms| {
        let primitive = base_homotopy(model, grade, terms)?;
        perturbation(model, &primitive)
    })
}

// Every successive summand raises the minimum rho exponent. The grade
// bounds all such exponents, and no dense filtration range is allocated.
fn geometric_sum(
    model: &mut Model,
    mut source: Cochain,
    mut correction: impl FnMut(&mut Model, &[Term]) -> Result<Cochain>,
) -> Result<Cochain> {
    let mut total = Vec::new();
    while let Some(first) = source.first() {
        model.spend(1)?;
        let next = correction(model, &source)?;
        if next.first().is_some_and(|term| term.rho <= first.rho) {
            return Err(Error::plain(
                "rho-minimal correction does not raise filtration",
            ));
        }
        total = if total.is_empty() {
            source
        } else {
            real_motivic::add(model, &total, &source)?
        };
        source = next;
    }
    Ok(total)
}

// The real differential preserves rho filtration. On each input level
// its terms at the same level are d0; retaining the higher ones gives Q.
fn perturbation(model: &mut Model, source: &[Term]) -> Result<Cochain> {
    let mut image = Sparse::new();
    for level in source.chunk_by(|a, b| a.rho == b.rho) {
        let f = level[0].rho;
        for term in real_motivic::differential_support(model, level)? {
            if term.rho > f {
                let cells = term.cells();
                image.toggle(model, term, cells)?;
            }
        }
    }
    Ok(image.into_vec())
}

fn base_lift(model: &mut Model, grade: Grade, support: &[Coordinate]) -> Result<Cochain> {
    let mut result = Vec::new();
    let mut cells = 0;
    for level in support.chunk_by(|a, b| a.0 == b.0) {
        let f = level[0].0;
        let grade = coefficient_grade(grade, f)?;
        model.spend(level.len())?;
        let words = if grade.s == 1 {
            let t = primitive_degree(grade)?.ok_or_else(grade_limit)?;
            let factor: cobar::Monomial = [t].into();
            vec![[factor].into()]
        } else {
            let ty = model.mext_type(grade)?;
            let coords = level.iter().map(|&(_, i)| i).collect::<Vec<_>>();
            model.adams(ty)?.representative(model, &coords)?
        };
        extend_level(model, &mut result, &mut cells, grade, f, words)?;
    }
    Ok(result)
}

fn base_project(model: &mut Model, grade: Grade, terms: &[Term]) -> Result<Vec<Coordinate>> {
    let mut result = Vec::new();
    for level in terms.chunk_by(|a, b| a.rho == b.rho) {
        let f = level[0].rho;
        let grade = coefficient_grade(grade, f)?;
        if grade.s == 1 {
            if let Some(t) = primitive_degree(grade)? {
                for term in level {
                    model.spend(term.cells())?;
                    if term.word[0].as_ref() == [t] {
                        check_cells(result.len() + 1)?;
                        result.push((f, 0));
                        break;
                    }
                }
            }
        } else {
            let words = words(model, level)?;
            let (_, coordinates) = complex_project(model, grade, words)?;
            check_cells(result.len() + coordinates.len())?;
            result.extend(coordinates.into_iter().map(|i| (f, i)));
        }
    }
    Ok(result)
}

fn primitive_degree(grade: Grade) -> Result<Option<usize>> {
    if grade.t <= 0 || i128::from(grade.w) > grade.t / 2 {
        return Ok(None);
    }
    let t = grade.classical_degree()?;
    // A degree-one complex class is a classical primitive satisfying the
    // weight bound. The classical Steenrod indecomposables Sq^(2^k) give
    // exactly xi_1^(2^k). Since this is a basis monomial, every other basis
    // monomial is in the chosen complement to the cycles.
    Ok(t.is_power_of_two().then_some(t))
}

fn base_homotopy(model: &mut Model, grade: Grade, terms: &[Term]) -> Result<Cochain> {
    // The complex-motivic differential from scalars is zero, so its
    // boundary section never has output in cohomological degree zero.
    if grade.s <= 1 {
        return Ok(Vec::new());
    }
    let mut result = Vec::new();
    let mut cells = 0;
    for level in terms.chunk_by(|a, b| a.rho == b.rho) {
        let f = level[0].rho;
        let grade = coefficient_grade(grade, f)?;
        let terms = words(model, level)?;
        let (closed, coords) = complex_project(model, grade, terms)?;
        if closed.is_empty() {
            continue;
        }
        let ty = model.mext_type(grade)?;
        let group = model.adams(ty)?;
        let representative = group.representative(model, &coords)?;
        let boundary = cobar::add_support(model, &closed, &representative)?;
        let primitive = group
            .primitive(model, &boundary)?
            .ok_or_else(missing_primitive)?;
        extend_level(
            model,
            &mut result,
            &mut cells,
            grade.previous().unwrap(),
            f,
            primitive,
        )?;
    }
    Ok(result)
}

fn words(model: &mut Model, level: &[Term]) -> Result<Vec<Word>> {
    let cells = level.iter().map(Term::cells).sum::<usize>();
    model.spend(cells * (1 + level.len().max(1).ilog2() as usize))?;
    let mut words: Vec<_> = level.iter().map(|term| term.word.clone()).collect();
    words.sort_unstable();
    Ok(words)
}

fn extend_level(
    model: &mut Model,
    result: &mut Cochain,
    cells: &mut usize,
    grade: Grade,
    rho: usize,
    words: Vec<Word>,
) -> Result<()> {
    let start = result.len();
    let mut new_cells = 0;
    let t = grade.classical_degree()?;
    for word in words {
        let weight = motivic::intrinsic_weight(&word, t);
        let tau =
            usize::try_from(weight as i128 - i128::from(grade.w)).map_err(|_| grade_limit())?;
        let term = Term { rho, tau, word };
        new_cells += term.cells();
        cobar::check_storage(*cells + new_cells)?;
        result.push(term);
    }
    model.spend(new_cells * (1 + (result.len() - start).max(1).ilog2() as usize))?;
    *cells += new_cells;
    result[start..].sort_unstable();
    Ok(())
}

fn complex_project(
    model: &mut Model,
    grade: Grade,
    source: Vec<Word>,
) -> Result<(Vec<Word>, Vec<usize>)> {
    let image = cobar::differential_support(model, &source)?;
    let closed = if image.is_empty() {
        source
    } else {
        let next_ty = model.mext_type(grade.differential()?)?;
        let primitive = model
            .adams(next_ty)?
            .primitive(model, &image)?
            .ok_or_else(missing_primitive)?;
        cobar::add_support(model, &source, &primitive)?
    };
    if closed.is_empty() {
        return Ok((closed, Vec::new()));
    }
    let ty = model.mext_type(grade)?;
    let coords = model.adams(ty)?.project(model, &closed)?;
    Ok((closed, coords))
}

fn missing_primitive() -> Error {
    Error::plain("complex-motivic contraction boundary has no primitive")
}
