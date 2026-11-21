use std::cmp::Ordering;
use std::rc::Rc;

use crate::cobar::{self, Word};
use crate::comodule::Generator;
use crate::error::{Error, Result};
use crate::map_ext::{self, Term};
use crate::model::{Model, TypeId, Value};

// Terms are sorted by source generator, target generator, and cobar word.
pub(crate) type Support = Rc<[Term]>;

pub(crate) fn generator(
    model: &mut Model,
    base: TypeId,
    degree: usize,
    index: usize,
) -> Result<Generator> {
    let ty = model.cohomology_type(base, degree)?;
    let rank = model.cohomology(ty)?.rank();
    if index >= rank {
        return Err(Error::plain(format!(
            "mword generator index {index} is out of range for {} of rank {rank}",
            model.type_name(ty)
        )));
    }
    Ok(Generator { degree, index })
}

pub(crate) fn internal_degree(source: usize, target: usize, word: usize) -> Result<i64> {
    i64::try_from(target as i128 + word as i128 - source as i128)
        .map_err(|_| Error::plain("mword internal degree is out of range"))
}

pub(crate) fn elementary(
    model: &mut Model,
    source: Generator,
    target: Generator,
    words: &[Word],
) -> Result<Support> {
    let cells = words
        .iter()
        .map(|word| Term::GENERATOR_CELLS + cobar::word_cells(word))
        .sum();
    reserve(model, cells)?;
    Ok(words
        .iter()
        .cloned()
        .map(|word| Term {
            source,
            target,
            word,
        })
        .collect())
}

pub(crate) fn add(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Support> {
    let terms = map_ext::add_terms(model, left, right)?;
    persist(model, terms)
}

pub(crate) fn compare(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Ordering> {
    // Read the highest differing coordinate first, as in support enumeration.
    for (a, b) in left.iter().rev().zip(right.iter().rev()) {
        model.spend(a.cells())?;
        let comparison = a.cmp(b);
        if comparison != Ordering::Equal {
            return Ok(comparison);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

pub(crate) fn differential(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    terms: &[Term],
) -> Result<Support> {
    let terms = map_ext::differential_support(model, source, target, terms)?;
    persist(model, terms)
}

pub(crate) fn compose(model: &mut Model, outer: &[Term], inner: &[Term]) -> Result<Support> {
    if outer.is_empty() || inner.is_empty() {
        return Ok([].into());
    }
    let terms = map_ext::compose_terms(model, outer, inner)?;
    persist(model, terms)
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let coordinates = model.map_ext(ty)?.project(model, terms)?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

pub(crate) fn representative(model: &mut Model, ty: TypeId, class: &[usize]) -> Result<Support> {
    if class.is_empty() {
        return Ok([].into());
    }
    let terms = model.map_ext(ty)?.representative(model, class)?;
    persist(model, terms)
}

pub(crate) fn persist(model: &mut Model, terms: Vec<Term>) -> Result<Support> {
    reserve(model, terms.iter().map(Term::cells).sum())?;
    Ok(terms.into())
}

fn reserve(model: &mut Model, cells: usize) -> Result<()> {
    map_ext::check_storage(cells)?;
    model.spend(cells)?;
    model.allocate(cells)
}
