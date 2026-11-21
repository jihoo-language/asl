use std::hash::Hash;
use std::rc::Rc;

use crate::cohomology;
use crate::error::{Error, Result};
use crate::linear::{Basis, QuotientBasis};
use crate::may::{sort_support, Rows};
use crate::model::{Kind, Model, TypeId, Value};
use crate::rv_hh::{self, Differential, Grade, Word};

pub(crate) const MAX_BASIS: usize = 32_768;

pub(crate) struct Slice<T> {
    terms: Vec<T>,
    quotient: QuotientBasis,
    cells: fn(&T) -> usize,
}

impl<T: Clone + Hash + Ord> Slice<T> {
    pub(crate) fn build(
        model: &mut Model,
        terms: Vec<T>,
        incoming: impl FnOnce(&mut Model) -> Result<Vec<T>>,
        cells: fn(&T) -> usize,
        label: &'static str,
        mut differential: impl FnMut(&mut Model, &T) -> Result<Vec<T>>,
    ) -> Result<Self> {
        if terms.is_empty() {
            return Ok(Self {
                terms,
                quotient: QuotientBasis::default(),
                cells,
            });
        }
        let mut boundaries = Basis::default();
        for term in incoming(model)? {
            let image = differential(model, &term)?;
            let coordinates = term_indices(model, &terms, &image, cells)?;
            boundaries.insert(model, coordinates)?;
        }
        let mut rows = Rows::new(label);
        let quotient = QuotientBasis::build(model, boundaries, terms.len(), |model, i| {
            let mut image = differential(model, &terms[i])?
                .into_iter()
                .map(|term| rows.coordinate(model, cells(&term), term))
                .collect::<Result<Vec<_>>>()?;
            sort_support(model, &mut image)?;
            Ok(image)
        })?;
        model.allocate(terms.iter().map(cells).sum())?;
        Ok(Self {
            terms,
            quotient,
            cells,
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.quotient.rank()
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        terms: &[T],
        operation: &str,
    ) -> Result<Vec<usize>> {
        let coordinates = term_indices(model, &self.terms, terms, self.cells)?;
        self.quotient
            .project(model, coordinates)?
            .ok_or_else(|| Error::plain(format!("{operation} requires a cycle")))
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<T>> {
        let indices = self.quotient.representative(model, class)?;
        Ok(indices.into_iter().map(|i| self.terms[i].clone()).collect())
    }
}

fn term_indices<T: Ord>(
    model: &mut Model,
    basis: &[T],
    terms: &[T],
    cells: fn(&T) -> usize,
) -> Result<Vec<usize>> {
    let cost = 1 + basis.len().max(1).ilog2() as usize;
    terms
        .iter()
        .map(|term| {
            model.spend(cells(term) * cost)?;
            basis
                .binary_search(term)
                .map_err(|_| Error::plain("expected a term in this homology bidegree"))
        })
        .collect()
}

pub(crate) fn build(model: &mut Model, grade: &Grade) -> Result<Slice<Word>> {
    let terms = rv_hh::basis(model, grade, MAX_BASIS)?;
    let ideal = model.rv_ideal_named(&grade.ideal).unwrap();
    let mut differential = Differential::new(model, grade, &ideal)?;
    Slice::build(
        model,
        terms,
        |model| rv_hh::basis(model, &grade.connes()?, MAX_BASIS),
        |term| rv_hh::word_cells(term),
        "Hochschild homology",
        |model, term| differential.apply(model, std::slice::from_ref(term)),
    )
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Word]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(cohomology::zero());
    }
    let class = model
        .rv_hh_homology(ty)?
        .project(model, terms, "rvHHClass")?;
    model.allocate(class.len())?;
    Ok(Value::Class(class.into()))
}

pub(crate) fn representative(model: &mut Model, ty: TypeId, class: &[usize]) -> Result<Rc<[Word]>> {
    if class.is_empty() {
        return Ok([].into());
    }
    let terms = model.rv_hh_homology(ty)?.representative(model, class)?;
    rv_hh::persist(model, terms)
}

pub(crate) fn connes(
    model: &mut Model,
    input: TypeId,
    class: &[usize],
    output: TypeId,
) -> Result<Value> {
    if class.is_empty() {
        return Ok(cohomology::zero());
    }
    let Kind::RvHochschild(grade) = model.kind(input).clone() else {
        unreachable!()
    };
    let lift = model.rv_hh_homology(input)?.representative(model, class)?;
    let image = rv_hh::connes(model, &grade, &lift)?;
    self::class(model, output, &image)
}
