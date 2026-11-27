use std::hash::Hash;

use crate::adams;
use crate::cobar::{self, Monomial, Word};
use crate::cohomology;
use crate::error::{Error, Result};
use crate::linear::{Basis, QuotientBasis};
use crate::may::Rows;
use crate::model::{Model, TypeId, Value};
use crate::motivic::{self, Grade};
use crate::real_motivic::{self, Term};
use crate::rho_localized;

const MAX_BASIS: usize = 32_768;

#[derive(Default)]
pub(crate) struct Group {
    terms: Vec<Term>,
    quotient: QuotientBasis,
}

impl Group {
    pub(crate) fn build(model: &mut Model, grade: Grade) -> Result<Self> {
        let terms = real_motivic::basis(model, grade, MAX_BASIS)?;
        if terms.is_empty() {
            return Ok(Self::default());
        }
        let mut boundaries = Basis::default();
        if let Some(previous) = grade.previous() {
            for term in real_motivic::basis(model, previous, MAX_BASIS)? {
                let image = real_motivic::differential_support(model, &[term])?;
                let column = indices(model, &terms, &image)?;
                boundaries.insert(model, column)?;
            }
        }
        let mut rows = Rows::new("real-motivic cobar");
        let quotient = QuotientBasis::build(model, boundaries, terms.len(), |model, index| {
            let image =
                real_motivic::differential_support(model, std::slice::from_ref(&terms[index]))?;
            column(model, &mut rows, image, Term::cells)
        })?;
        Ok(Self { terms, quotient })
    }

    pub(crate) fn rank(&self) -> usize {
        self.quotient.rank()
    }

    pub(crate) fn representative(&self, model: &mut Model, coords: &[usize]) -> Result<Vec<Term>> {
        Ok(self
            .quotient
            .representative(model, coords)?
            .into_iter()
            .map(|coordinate| self.terms[coordinate].clone())
            .collect())
    }

    pub(crate) fn project(&self, model: &mut Model, terms: &[Term]) -> Result<Vec<usize>> {
        let support = indices(model, &self.terms, terms)?;
        self.quotient
            .project(model, support)?
            .ok_or_else(not_closed)
    }
}

fn not_closed() -> Error {
    Error::plain("operation requires a closed real-motivic cobar cochain")
}

fn indices(model: &mut Model, basis: &[Term], terms: &[Term]) -> Result<Vec<usize>> {
    let steps = 1 + basis.len().max(1).ilog2() as usize;
    terms
        .iter()
        .map(|term| {
            model.spend(steps * term.cells())?;
            basis
                .binary_search(term)
                .map_err(|_| Error::plain("word is outside its real-motivic cobar tridegree"))
        })
        .collect()
}

fn column<T: Hash + Eq>(
    model: &mut Model,
    rows: &mut Rows<T>,
    support: impl IntoIterator<Item = T>,
    cells: fn(&T) -> usize,
) -> Result<Vec<usize>> {
    let mut column = support
        .into_iter()
        .map(|term| rows.coordinate(model, cells(&term), term))
        .collect::<Result<Vec<_>>>()?;
    model.spend(column.len() * (1 + column.len().max(1).ilog2() as usize))?;
    column.sort_unstable();
    Ok(column)
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(cohomology::zero());
    }
    let coords = model.real_ext(ty)?.project(model, terms)?;
    model.allocate(coords.len())?;
    Ok(Value::Class(coords.into()))
}

pub(crate) fn representative(model: &mut Model, input: TypeId, coords: &[usize]) -> Result<Value> {
    if coords.is_empty() {
        return Ok(real_motivic::zero());
    }
    let terms = model.real_ext(input)?.representative(model, coords)?;
    real_motivic::persist(model, terms)
}

pub(crate) fn yoneda(
    model: &mut Model,
    left: (TypeId, &[usize]),
    right: (TypeId, &[usize]),
    output: TypeId,
) -> Result<Value> {
    if left.1.is_empty() || right.1.is_empty() {
        return Ok(cohomology::zero());
    }
    let a = model.real_ext(left.0)?.representative(model, left.1)?;
    let b = model.real_ext(right.0)?.representative(model, right.1)?;
    let terms = real_motivic::concat_support(model, &a, &b)?;
    class(model, output, &terms)
}

pub(crate) fn shift(
    model: &mut Model,
    input: TypeId,
    coords: &[usize],
    output: TypeId,
    power: usize,
) -> Result<Value> {
    if coords.is_empty() {
        return Ok(cohomology::zero());
    }
    let lift = model.real_ext(input)?.representative(model, coords)?;
    let terms = real_motivic::shift_support(model, &lift, power)?;
    class(model, output, &terms)
}

pub(crate) fn specialize(
    model: &mut Model,
    input: TypeId,
    coords: &[usize],
    output: TypeId,
) -> Result<Value> {
    if coords.is_empty() {
        return Ok(cohomology::zero());
    }
    let lift = model.real_ext(input)?.representative(model, coords)?;
    let Value::Cobar(words) = real_motivic::specialize(model, &lift)? else {
        unreachable!();
    };
    adams::class(model, output, &words)
}

pub(crate) fn project_classical(
    model: &mut Model,
    input: TypeId,
    coords: &[usize],
    output: TypeId,
) -> Result<Value> {
    if coords.is_empty() {
        return Ok(cohomology::zero());
    }
    let lift = model.real_ext(input)?.representative(model, coords)?;
    let words = classical_words(model, &lift)?;
    adams::class(model, output, &words)
}

fn classical_words(model: &mut Model, terms: &[Term]) -> Result<Vec<Word>> {
    let mut projected = cobar::Support::new();
    for term in terms {
        model.spend(term.cells())?;
        // Each tau_i has zero image after killing the vertex and the U
        // variables. The remaining factors are the pure classical xi words.
        if term.tau != 0
            || term
                .word
                .iter()
                .flat_map(|factor| factor.iter())
                .any(|&power| power & 1 != 0)
        {
            continue;
        }
        let word: Word = term
            .word
            .iter()
            .map(|factor| factor.iter().map(|power| power / 2).collect::<Monomial>())
            .collect();
        let cells = cobar::word_cells(&word);
        projected.toggle(model, word, cells)?;
    }
    Ok(projected.into_vec())
}

pub(crate) fn primitive(model: &mut Model, grade: Grade, terms: &[Term]) -> Result<Value> {
    let witness = real_preimage(model, grade, terms, true)?
        .ok_or_else(|| Error::plain("extPrimitive requires a real-motivic cobar boundary"))?;
    real_motivic::persist(model, witness)
}

// Check the actual incoming differential. Stopping once the target
// reduces to zero avoids reducing unrelated later boundary columns.
fn preimage<T: Clone + Eq + Hash>(
    model: &mut Model,
    source: &[T],
    terms: &[T],
    cells: fn(&T) -> usize,
    mut differential: impl FnMut(&mut Model, &T) -> Result<Vec<T>>,
    witness: bool,
) -> Result<Option<Vec<T>>> {
    let mut rows = Rows::new("cobar boundary");
    let mut remainder = column(model, &mut rows, terms.iter().cloned(), cells)?;
    let mut coefficients = Vec::new();
    let mut span = Basis::default();
    for (index, term) in source.iter().enumerate() {
        let image = differential(model, term)?;
        let image = column(model, &mut rows, image, cells)?;
        let incoming = if witness { vec![index] } else { Vec::new() };
        if span.insert_with_witness(model, image, incoming)?.is_none()
            && remainder.last().is_some_and(|&j| span.contains_pivot(j))
        {
            (remainder, coefficients) = span.reduce_with_witness(model, remainder, coefficients)?;
            if remainder.is_empty() {
                return Ok(Some(
                    coefficients
                        .into_iter()
                        .map(|i| source[i].clone())
                        .collect(),
                ));
            }
        }
    }
    Ok(None)
}

pub(crate) fn real_preimage(
    model: &mut Model,
    grade: Grade,
    terms: &[Term],
    witness: bool,
) -> Result<Option<Vec<Term>>> {
    if terms.is_empty() {
        return Ok(Some(Vec::new()));
    }
    let Some(previous) = grade.previous() else {
        return Ok(None);
    };
    if previous.s == 0 {
        let source = real_motivic::basis(model, previous, 1)?;
        if let [scalar] = source.as_slice() {
            let image = real_motivic::differential_support(model, std::slice::from_ref(scalar))?;
            if real_motivic::compare(model, &image, terms)? == std::cmp::Ordering::Equal {
                return Ok(Some(if witness { source } else { Vec::new() }));
            }
        }
        return Ok(None);
    }

    if let Some(primitive) = left_tau_preimage(model, terms)? {
        return Ok(Some(if witness { primitive } else { Vec::new() }));
    }

    // A non-boundary after specialization is already a certificate
    // against a real primitive, without using the other rho levels.
    let complex = real_motivic::complex_support(model, terms)?;
    if !complex.is_empty() {
        let source = motivic::basis(model, previous, MAX_BASIS)?;
        if preimage(
            model,
            &source,
            &complex,
            cobar::word_cells,
            cobar::differential_word,
            false,
        )?
        .is_none()
        {
            return Ok(None);
        }
    }

    let source = real_motivic::basis(model, previous, MAX_BASIS)?;
    preimage(
        model,
        &source,
        terms,
        Term::cells,
        |model, term| real_motivic::differential_support(model, std::slice::from_ref(term)),
        witness,
    )
}

fn left_tau_preimage(model: &mut Model, terms: &[Term]) -> Result<Option<Vec<Term>>> {
    let mut candidate = Vec::with_capacity(terms.len());
    for term in terms {
        model.spend(term.cells())?;
        let (Some(rho), Some(tau), Some((factor, tail))) = (
            term.rho.checked_sub(1),
            term.tau.checked_add(1),
            term.word.split_first(),
        ) else {
            return Ok(None);
        };
        if factor.as_ref() != [1] {
            return Ok(None);
        }
        candidate.push(Term {
            rho,
            tau,
            word: tail.into(),
        });
    }
    // d(tau) = rho[tau_0]. Verify the differential of the proposed
    // left product before returning it as a primitive.
    let image = real_motivic::differential_support(model, &candidate)?;
    let exact = real_motivic::compare(model, &image, terms)? == std::cmp::Ordering::Equal;
    Ok(exact.then_some(candidate))
}

fn classical_boundary(model: &mut Model, grade: Grade, terms: &[Term]) -> Result<bool> {
    let words = classical_words(model, terms)?;
    if words.is_empty() {
        return Ok(true);
    }
    let Some(previous) = grade.s.checked_sub(1) else {
        return Ok(false);
    };
    let k = rho_localized::classical_degree(grade)?;
    let source = cobar::basis(model, previous, k, MAX_BASIS)?;
    Ok(preimage(
        model,
        &source,
        &words,
        cobar::word_cells,
        cobar::differential_word,
        false,
    )?
    .is_some())
}

fn shifted_preimage(
    model: &mut Model,
    grade: Grade,
    terms: &[Term],
    power: usize,
    witness: bool,
) -> Result<Option<Vec<Term>>> {
    if power == 0 {
        return real_preimage(model, grade, terms, witness);
    }
    let shifted = real_motivic::shifted(grade, power)?;
    let shifted_terms = real_motivic::shift_support(model, terms, power)?;
    real_preimage(model, shifted, &shifted_terms, witness)
}

fn verify_cycle(model: &mut Model, terms: &[Term], known_cycle: bool) -> Result<()> {
    if !known_cycle && !real_motivic::differential_support(model, terms)?.is_empty() {
        return Err(not_closed());
    }
    Ok(())
}

pub(crate) fn order(
    model: &mut Model,
    grade: Grade,
    terms: &[Term],
    known_cycle: bool,
) -> Result<Option<usize>> {
    if terms.is_empty() {
        return Ok(Some(0));
    }
    if grade.s == 0 {
        verify_cycle(model, terms, known_cycle)?;
        return Ok(None);
    }

    // For t<=0 the real cobar already equals the localized slice.
    if grade.t <= 0 && (known_cycle || grade.s > 1) {
        if known_cycle {
            return Ok(None);
        }
        verify_cycle(model, terms, false)?;
        return Ok(classical_boundary(model, grade, terms)?.then_some(0));
    }

    // In length one there is no possible scalar primitive before t<=0.
    if grade.s == 1 {
        let bound = usize::try_from(grade.t.max(0)).map_err(|_| exponent_overflow())?;
        if shifted_preimage(model, grade, terms, bound, false)?.is_some() {
            return Ok(Some(bound));
        }
        verify_cycle(model, terms, known_cycle)?;
        return Ok(None);
    }

    if !classical_boundary(model, grade, terms)? {
        verify_cycle(model, terms, known_cycle)?;
        return Ok(None);
    }
    let mut hi = usize::try_from(grade.t).map_err(|_| exponent_overflow())?;
    let mut lo = usize::from(known_cycle);
    let mut known_boundary = false;
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if shifted_preimage(model, grade, terms, mid, false)?.is_some() {
            hi = mid;
            known_boundary = true;
        } else {
            lo = mid + 1;
        }
    }
    if !known_boundary && shifted_preimage(model, grade, terms, hi, false)?.is_none() {
        verify_cycle(model, terms, known_cycle)?;
        return Err(Error::plain(
            "real-motivic cobar boundary disagrees with its stable classical image",
        ));
    }
    // An actual boundary rho^hi z certifies that z was closed, since
    // multiplication by rho is injective on the normalized cochains.
    Ok(Some(hi))
}

fn exponent_overflow() -> Error {
    Error::plain("rho exponent is outside the supported range")
}

pub(crate) fn null(
    model: &mut Model,
    grade: Grade,
    terms: &[Term],
    power: usize,
    known_cycle: bool,
) -> Result<Value> {
    let Some(witness) = shifted_preimage(model, grade, terms, power, true)? else {
        verify_cycle(model, terms, known_cycle)?;
        return Err(Error::plain("rhoNull power does not kill the class"));
    };
    if power != 0 && shifted_preimage(model, grade, terms, power - 1, false)?.is_some() {
        return Err(Error::plain("rhoNull power is not the minimal exponent"));
    }
    real_motivic::persist(model, witness)
}
