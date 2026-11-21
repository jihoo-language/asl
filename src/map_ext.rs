use std::collections::hash_map::Entry;
use std::collections::{BTreeSet, HashMap};
use std::ops::RangeInclusive;
use std::rc::Rc;

use crate::cobar::{self, Word};
use crate::comodule::{self, Cochains, Generator};
use crate::error::{Error, Result};
use crate::linear::{Basis, QuotientBasis};
use crate::milnor;
use crate::model::{Model, TypeId, Value};

const MAX_BASIS_TERMS: usize = 32_768;
const MAX_COCHAIN_CELLS: usize = 300_000;

// The term sends the source generator to word tensor target.
#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Term {
    pub(crate) source: Generator,
    pub(crate) target: Generator,
    pub(crate) word: Word,
}

impl Term {
    fn cells(&self) -> usize {
        4 + cobar::word_cells(&self.word)
    }
}

#[derive(Default)]
pub(crate) struct MapAdams {
    terms: Vec<Term>,
    quotient: QuotientBasis,
}

impl MapAdams {
    pub(crate) fn build(
        model: &mut Model,
        source: TypeId,
        target: TypeId,
        s: usize,
        t: i64,
    ) -> Result<Self> {
        let Some(source_max) = model.height(source, false)?.checked_sub(1) else {
            return Ok(Self::default());
        };
        let Some(target_max) = model.height(target, false)?.checked_sub(1) else {
            return Ok(Self::default());
        };
        let t = i128::from(t);
        if source_degree_bounds(s, t, source_max, target_max).is_none() {
            return Ok(Self::default());
        }
        let mut complex = MapCochains {
            source,
            target,
            source_max: refine_dimension_bound(model, source, source_max)?,
            target_max: refine_dimension_bound(model, target, target_max)?,
            t,
            source_coactions: HashMap::new(),
            source_coaction_cells: 0,
            target_cochains: Cochains::new(target),
        };
        let terms = complex.basis(model, s)?;
        if terms.is_empty() {
            return Ok(Self::default());
        }

        let mut boundaries = Basis::default();
        complex.for_each_boundary(model, s, &terms, |model, _, column| {
            boundaries.insert(model, column)?;
            Ok(())
        })?;

        let mut rows = HashMap::new();
        let mut row_cells = 0;
        let quotient = QuotientBasis::build(model, boundaries, terms.len(), |model, index| {
            let mut column = Vec::new();
            for term in complex.differential(model, &terms[index])? {
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
        self.project_indices(model, column)
    }

    fn project_indices(&self, model: &mut Model, column: Vec<usize>) -> Result<Vec<usize>> {
        self.quotient
            .project(model, column)?
            .ok_or_else(|| Error::plain("expected a closed map cobar cochain"))
    }
}

struct MapCochains {
    source: TypeId,
    target: TypeId,
    source_max: usize,
    target_max: usize,
    t: i128,
    source_coactions: HashMap<Generator, Vec<comodule::Term>>,
    source_coaction_cells: usize,
    target_cochains: Cochains,
}

impl MapCochains {
    fn basis(&self, model: &mut Model, s: usize) -> Result<Vec<Term>> {
        let length = s as i128;
        let Some((first_source, last_source)) =
            source_degree_bounds(s, self.t, self.source_max, self.target_max)
        else {
            return Ok(Vec::new());
        };
        let first_target = if s == 0 { first_source + self.t } else { 0 };
        let last_target = (last_source + self.t - length).min(self.target_max as i128);

        let source = graded_ranks(
            model,
            self.source,
            first_source as usize..=last_source as usize,
        )?;
        if source.is_empty() {
            return Ok(Vec::new());
        }
        let target = graded_ranks(
            model,
            self.target,
            first_target as usize..=last_target as usize,
        )?;
        let mut catalogue = HashMap::<usize, Vec<Word>>::new();
        let mut terms = Vec::new();
        let mut cells = 0;
        for (q, source_rank) in source {
            let total = q as i128 + self.t;
            let min = if s == 0 { total } else { 0 };
            let max = total - length;
            let first = target.partition_point(|&(degree, _)| (degree as i128) < min);
            for &(p, target_rank) in &target[first..] {
                model.spend(1)?;
                if p as i128 > max {
                    break;
                }
                let internal = usize::try_from(total - p as i128)
                    .map_err(|_| Error::plain("map cobar word degree is out of range"))?;
                let copies = source_rank.saturating_mul(target_rank);
                let available = (MAX_BASIS_TERMS - terms.len()) / copies;
                let words = match catalogue.entry(internal) {
                    Entry::Occupied(entry) => entry.into_mut(),
                    Entry::Vacant(entry) => {
                        if available == 0 {
                            return Err(basis_limit());
                        }
                        entry.insert(cobar::basis(model, s, internal, available)?)
                    }
                };
                if words.len() > available {
                    return Err(basis_limit());
                }
                for index in 0..source_rank {
                    let source = Generator { degree: q, index };
                    for index in 0..target_rank {
                        let target = Generator { degree: p, index };
                        for word in words.iter() {
                            let term = Term {
                                source,
                                target,
                                word: word.clone(),
                            };
                            add_cells(&mut cells, term.cells())?;
                            terms.push(term);
                        }
                    }
                }
            }
        }
        model.spend(cells * (1 + terms.len().max(1).ilog2() as usize))?;
        model.allocate(5 * terms.len())?;
        terms.sort_unstable();
        Ok(terms)
    }

    // Returning the ordered source terms lets a later preimage reduction
    // translate tracked column indices back into cochains.
    fn for_each_boundary(
        &mut self,
        model: &mut Model,
        s: usize,
        target: &[Term],
        mut use_column: impl FnMut(&mut Model, usize, Vec<usize>) -> Result<()>,
    ) -> Result<Vec<Term>> {
        let Some(previous) = s.checked_sub(1) else {
            return Ok(Vec::new());
        };
        let source = self.basis(model, previous)?;
        for (index, term) in source.iter().enumerate() {
            let image = self.differential(model, term)?;
            let column = term_indices(model, target, &image)?;
            use_column(model, index, column)?;
        }
        Ok(source)
    }

    fn differential(&mut self, model: &mut Model, term: &Term) -> Result<Vec<Term>> {
        let mut support = Support::default();
        for word in cobar::differential_word(model, &term.word)? {
            support.toggle(
                model,
                Term {
                    word,
                    ..term.clone()
                },
            )?;
        }
        for component in self.target_cochains.coaction(model, term.target)? {
            let word = cobar::concat_word_support(model, &term.word, &component.word)?;
            support.toggle(
                model,
                Term {
                    source: term.source,
                    target: component.generator,
                    word,
                },
            )?;
        }
        for component in self.transpose_coaction(model, term.source)? {
            let word = cobar::concat_word_support(model, &component.word, &term.word)?;
            support.toggle(
                model,
                Term {
                    source: component.generator,
                    target: term.target,
                    word,
                },
            )?;
        }
        Ok(support.into_vec())
    }

    fn transpose_coaction(
        &mut self,
        model: &mut Model,
        generator: Generator,
    ) -> Result<&[comodule::Term]> {
        // Positive operations annihilate degree zero, so no reduced coaction
        // can land in this component.
        if generator.degree == 0 || generator.degree >= self.source_max {
            return Ok(&[]);
        }
        Ok(match self.source_coactions.entry(generator) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => {
                let mut components = Vec::new();
                for (degree, rank) in
                    graded_ranks(model, self.source, generator.degree + 1..=self.source_max)?
                {
                    let class = Rc::<[usize]>::from([generator.index]);
                    for index in 0..rank {
                        let above = Generator { degree, index };
                        for word in milnor::coact(
                            model,
                            self.source,
                            degree,
                            index,
                            generator.degree,
                            class.clone(),
                        )? {
                            add_cells(
                                &mut self.source_coaction_cells,
                                3 + cobar::word_cells(&word),
                            )?;
                            components.push(comodule::Term {
                                generator: above,
                                word,
                            });
                        }
                    }
                }
                entry.insert(components)
            }
        })
    }
}

pub(crate) fn identity(model: &mut Model, ty: TypeId) -> Result<Value> {
    let adams = model.map_ext(ty)?;
    model.spend(adams.terms.len())?;
    let diagonal = adams
        .terms
        .iter()
        .enumerate()
        .filter_map(|(index, term)| (term.source == term.target).then_some(index));
    let coordinates = adams.project_indices(model, diagonal.collect())?;
    class_value(model, coordinates)
}

pub(crate) fn induced(
    model: &mut Model,
    function_ty: TypeId,
    table: &[Value],
    ty: TypeId,
) -> Result<Value> {
    let adams = model.map_ext(ty)?;
    let mut support = Vec::new();
    let mut coordinate = 0;
    for terms in adams
        .terms
        .chunk_by(|a, b| a.source.degree == b.source.degree)
    {
        let images = comodule::homology_map(model, function_ty, table, terms[0].source.degree)?;
        for term in terms {
            let image = &images[term.source.index];
            model.spend(1 + image.len().max(1).ilog2() as usize)?;
            if image.binary_search(&term.target.index).is_ok() {
                support.push(coordinate);
            }
            coordinate += 1;
        }
    }
    let coordinates = adams.project_indices(model, support)?;
    class_value(model, coordinates)
}

pub(crate) fn compose(
    model: &mut Model,
    (outer_ty, outer): (TypeId, &[usize]),
    (inner_ty, inner): (TypeId, &[usize]),
    ty: TypeId,
) -> Result<Value> {
    if outer.is_empty() || inner.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let outer = model.map_ext(outer_ty)?.representative(model, outer)?;
    let inner = model.map_ext(inner_ty)?.representative(model, inner)?;
    let terms = compose_terms(model, &outer, &inner)?;
    if terms.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let coordinates = model.map_ext(ty)?.project(model, &terms)?;
    class_value(model, coordinates)
}

pub(crate) fn compose_terms(
    model: &mut Model,
    outer: &[Term],
    inner: &[Term],
) -> Result<Vec<Term>> {
    let mut support = Support::default();
    let search = 1 + outer.len().max(1).ilog2() as usize;
    for f in inner {
        model.spend(search)?;
        let start = outer.partition_point(|g| g.source < f.target);
        for g in &outer[start..] {
            model.spend(1)?;
            if g.source != f.target {
                break;
            }
            let word = cobar::concat_word_support(model, &f.word, &g.word)?;
            support.toggle(
                model,
                Term {
                    source: f.source,
                    target: g.target,
                    word,
                },
            )?;
        }
    }
    Ok(support.into_vec())
}

fn source_degree_bounds(s: usize, t: i128, source: usize, target: usize) -> Option<(i128, i128)> {
    // The image of degree q has degree q+t and the word has degree at least s.
    // At length zero the word is the unit, so q+t is a target degree.
    let first = (s as i128 - t).max(0);
    let mut last = source as i128;
    if s == 0 {
        last = last.min(target as i128 - t);
    }
    (first <= last).then_some((first, last))
}

fn refine_dimension_bound(model: &mut Model, ty: TypeId, bound: usize) -> Result<usize> {
    if bound == usize::MAX - 1 {
        Ok(model.height(ty, true)? - 1)
    } else {
        Ok(bound)
    }
}

fn graded_ranks(
    model: &mut Model,
    base: TypeId,
    degrees: RangeInclusive<usize>,
) -> Result<Vec<(usize, usize)>> {
    let mut ranks = Vec::new();
    for degree in degrees {
        model.spend(1)?;
        let ty = model.cohomology_type(base, degree)?;
        let rank = model.cohomology(ty)?.rank();
        if rank > 0 {
            ranks.push((degree, rank));
        }
    }
    Ok(ranks)
}

fn term_indices(model: &mut Model, basis: &[Term], terms: &[Term]) -> Result<Vec<usize>> {
    let search = 1 + basis.len().max(1).ilog2() as usize;
    terms
        .iter()
        .map(|term| {
            model.spend(search * term.cells())?;
            basis
                .binary_search(term)
                .map_err(|_| Error::plain("term is outside its map cobar bidegree"))
        })
        .collect()
}

fn class_value(model: &mut Model, coordinates: Vec<usize>) -> Result<Value> {
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

#[derive(Default)]
struct Support {
    terms: BTreeSet<Term>,
    cells: usize,
}

impl Support {
    fn toggle(&mut self, model: &mut Model, term: Term) -> Result<()> {
        let size = term.cells();
        model.spend(size * (1 + self.terms.len().max(1).ilog2() as usize))?;
        if self.terms.remove(&term) {
            self.cells -= size;
        } else {
            add_cells(&mut self.cells, size)?;
            self.terms.insert(term);
        }
        Ok(())
    }

    fn into_vec(self) -> Vec<Term> {
        self.terms.into_iter().collect()
    }
}

fn add_cells(cells: &mut usize, size: usize) -> Result<()> {
    *cells += size;
    if *cells > MAX_COCHAIN_CELLS {
        return Err(Error::plain(
            "map cobar calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

fn basis_limit() -> Error {
    Error::plain("map cobar basis exceeds the enumeration limit")
}
