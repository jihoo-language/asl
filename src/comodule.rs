use std::collections::{hash_map::Entry, BTreeSet, HashMap};
use std::rc::Rc;

use crate::cobar::{self, Word};
use crate::cohomology;
use crate::error::{Error, Result};
use crate::linear::{Basis, QuotientBasis};
use crate::milnor;
use crate::model::{Kind, Model, TypeId, Value};

const MAX_BASIS_TERMS: usize = 32_768;
const MAX_COCHAIN_CELLS: usize = 300_000;

#[derive(Clone, Copy, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Generator {
    pub(crate) degree: usize,
    pub(crate) index: usize,
}

#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Term {
    pub(crate) generator: Generator,
    pub(crate) word: Word,
}

impl Term {
    fn cells(&self) -> usize {
        2 + cobar::word_cells(&self.word)
    }
}

#[derive(Default)]
pub(crate) struct Adams {
    terms: Vec<Term>,
    quotient: QuotientBasis,
}

impl Adams {
    pub(crate) fn build(model: &mut Model, base: TypeId, s: usize, t: usize) -> Result<Self> {
        let mut cochains = Cochains::new(base);
        let terms = cochains.basis(model, s, t)?;
        if terms.is_empty() {
            return Ok(Self::default());
        }

        let mut boundaries = Basis::default();
        if let Some(previous) = s.checked_sub(1) {
            for term in cochains.basis(model, previous, t)? {
                let image = cochains.differential(model, &term)?;
                let column = term_indices(model, &terms, &image)?;
                boundaries.insert(model, column)?;
            }
        }

        // Number only the rows that occur in the outgoing differential.
        let mut rows = HashMap::new();
        let mut row_cells = 0;
        let quotient = QuotientBasis::build(model, boundaries, terms.len(), |model, index| {
            let mut column = Vec::new();
            for term in cochains.differential(model, &terms[index])? {
                let cells = term.cells();
                model.spend(cells)?;
                let next = rows.len();
                let coordinate = match rows.entry(term) {
                    Entry::Occupied(entry) => *entry.get(),
                    Entry::Vacant(entry) => {
                        add_cells(&mut row_cells, cells)?;
                        *entry.insert(next)
                    }
                };
                column.push(coordinate);
            }
            model.spend(column.len() * (1 + column.len().max(1).ilog2() as usize))?;
            column.sort_unstable();
            Ok(column)
        })?;
        Ok(Self { terms, quotient })
    }

    pub(crate) fn rank(&self) -> usize {
        self.quotient.rank()
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<Term>> {
        Ok(self
            .quotient
            .representative(model, class)?
            .into_iter()
            .map(|index| self.terms[index].clone())
            .collect())
    }

    pub(crate) fn project(&self, model: &mut Model, terms: &[Term]) -> Result<Vec<usize>> {
        let column = term_indices(model, &self.terms, terms)?;
        self.quotient
            .project(model, column)?
            .ok_or_else(|| Error::plain("expected a closed coefficient cobar cochain"))
    }
}

pub(crate) struct Cochains {
    base: TypeId,
    coactions: HashMap<Generator, Vec<Term>>,
}

impl Cochains {
    pub(crate) fn new(base: TypeId) -> Self {
        Self {
            base,
            coactions: HashMap::new(),
        }
    }

    fn basis(&self, model: &mut Model, s: usize, t: usize) -> Result<Vec<Term>> {
        let Some(max_degree) = t.checked_sub(s) else {
            return Ok(Vec::new());
        };
        let dimension_bound = match model.height(self.base, false)? {
            0 => return Ok(Vec::new()),
            usize::MAX => usize::MAX,
            height => height - 1,
        };
        // Every factor has positive degree, and the length-zero word has degree zero.
        let min_degree = if s == 0 { t } else { 0 };
        let mut terms = Vec::new();
        let mut cells = 0;
        for degree in min_degree..=max_degree.min(dimension_bound) {
            let ty = model.cohomology_type(self.base, degree)?;
            let rank = model.cohomology(ty)?.rank();
            if rank == 0 {
                continue;
            }
            let available = (MAX_BASIS_TERMS - terms.len()) / rank;
            if available == 0 {
                return Err(Error::plain(
                    "coefficient cobar basis exceeds the enumeration limit",
                ));
            }
            let words = cobar::basis(model, s, t - degree, available)?;
            for index in 0..rank {
                let generator = Generator { degree, index };
                for word in &words {
                    let term = Term {
                        generator,
                        word: word.clone(),
                    };
                    add_cells(&mut cells, term.cells())?;
                    terms.push(term);
                }
            }
        }
        model.spend(cells)?;
        model.allocate(3 * terms.len())?;
        Ok(terms)
    }

    // The reduced coaction lowers the coefficient degree. Its word is
    // appended on the right, where the coefficient sits in the cobar complex.
    fn differential(&mut self, model: &mut Model, term: &Term) -> Result<Vec<Term>> {
        let mut image = Vec::new();
        let mut cells = 0;
        for component in self.coaction(model, term.generator)? {
            let size = term.cells() + cobar::word_cells(&component.word) - 1;
            add_cells(&mut cells, size)?;
            model.spend(size)?;
            image.push(Term {
                generator: component.generator,
                word: cobar::join(&term.word, &component.word),
            });
        }
        for word in cobar::differential_word(model, &term.word)? {
            let output = Term {
                generator: term.generator,
                word,
            };
            add_cells(&mut cells, output.cells())?;
            image.push(output);
        }
        Ok(image)
    }

    pub(crate) fn coaction(&mut self, model: &mut Model, generator: Generator) -> Result<&[Term]> {
        Ok(match self.coactions.entry(generator) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => {
                let mut terms = Vec::new();
                let mut cells = 0;
                // Positive-degree Steenrod operations vanish on H^0.
                for degree in 1..generator.degree {
                    let ty = model.cohomology_type(self.base, degree)?;
                    let rank = model.cohomology(ty)?.rank();
                    for index in 0..rank {
                        let words = milnor::coact(
                            model,
                            self.base,
                            generator.degree,
                            generator.index,
                            degree,
                            Rc::from([index]),
                        )?;
                        let generator = Generator { degree, index };
                        for word in words {
                            let term = Term { generator, word };
                            add_cells(&mut cells, term.cells())?;
                            terms.push(term);
                        }
                    }
                }
                model.allocate(1 + cells)?;
                entry.insert(terms)
            }
        })
    }
}

pub(crate) fn action(
    model: &mut Model,
    scalar_ty: TypeId,
    scalar: &[usize],
    coefficient_ty: TypeId,
    coefficient: &[usize],
    result_ty: TypeId,
) -> Result<Value> {
    if scalar.is_empty() || coefficient.is_empty() {
        return Ok(cohomology::zero());
    }
    let scalars = model.adams(scalar_ty)?.representative(model, scalar)?;
    let coefficients = model
        .xext(coefficient_ty)?
        .representative(model, coefficient)?;
    let mut terms = Vec::new();
    let mut cells = 0;
    // The common left word length determines a unique split in every term.
    for coefficient in coefficients {
        for scalar in &scalars {
            let word = cobar::concat_words(model, scalar, &coefficient.word)?;
            let term = Term {
                generator: coefficient.generator,
                word,
            };
            add_cells(&mut cells, term.cells())?;
            terms.push(term);
        }
    }
    model.spend(cells * (1 + terms.len().max(1).ilog2() as usize))?;
    terms.sort_unstable();
    class(model, result_ty, &terms)
}

pub(crate) fn map(
    model: &mut Model,
    function_ty: TypeId,
    table: &[Value],
    source_ty: TypeId,
    source: &[usize],
    result_ty: TypeId,
) -> Result<Value> {
    if source.is_empty() {
        return Ok(cohomology::zero());
    }
    let source = model.xext(source_ty)?.representative(model, source)?;
    let mut support = BTreeSet::new();
    let mut cells = 0;
    for terms in source.chunk_by(|a, b| a.generator.degree == b.generator.degree) {
        let degree = terms[0].generator.degree;
        let images = homology_map(model, function_ty, table, degree)?;
        for term in terms {
            for &index in &images[term.generator.index] {
                let generator = Generator { degree, index };
                let output = Term {
                    generator,
                    word: term.word.clone(),
                };
                let size = output.cells();
                model.spend(size * (1 + support.len().max(1).ilog2() as usize))?;
                if support.remove(&output) {
                    cells -= size;
                } else {
                    add_cells(&mut cells, size)?;
                    support.insert(output);
                }
            }
        }
    }
    let support: Vec<_> = support.into_iter().collect();
    class(model, result_ty, &support)
}

pub(crate) fn homology_map(
    model: &mut Model,
    function_ty: TypeId,
    table: &[Value],
    degree: usize,
) -> Result<Vec<Vec<usize>>> {
    let Kind::Arrow(domain, range) = *model.kind(function_ty) else {
        unreachable!()
    };
    let domain_ty = model.cohomology_type(domain, degree)?;
    let range_ty = model.cohomology_type(range, degree)?;
    let from_rank = model.cohomology(domain_ty)?.rank();
    let to_rank = model.cohomology(range_ty)?.rank();
    let mut images = vec![Vec::new(); from_rank];
    let mut cells = from_rank;
    // The rows of pullback are the columns of the induced homology map.
    for j in 0..to_rank {
        let Value::Class(pullback) =
            cohomology::pull(model, function_ty, table, range_ty, &[j], domain_ty)?
        else {
            unreachable!()
        };
        add_cells(&mut cells, pullback.len())?;
        for &i in pullback.iter() {
            images[i].push(j);
        }
    }
    Ok(images)
}

fn class(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(cohomology::zero());
    }
    let coordinates = model.xext(ty)?.project(model, terms)?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

fn term_indices(model: &mut Model, basis: &[Term], terms: &[Term]) -> Result<Vec<usize>> {
    let search_steps = 1 + basis.len().max(1).ilog2() as usize;
    terms
        .iter()
        .map(|term| {
            model.spend(search_steps * term.cells())?;
            basis
                .binary_search(term)
                .map_err(|_| Error::plain("term is outside its coefficient cobar bidegree"))
        })
        .collect()
}

fn add_cells(cells: &mut usize, extra: usize) -> Result<()> {
    *cells += extra;
    if *cells > MAX_COCHAIN_CELLS {
        return Err(Error::plain(
            "coefficient cobar calculation exceeds the storage limit",
        ));
    }
    Ok(())
}
