use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use crate::cohomology;
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::linear::{Basis, QuotientBasis};
use crate::may::{sort_support, Rows};
use crate::model::{Model, TypeId, Value};
use crate::rv_bar::{self, Factor, Grade, Path, PathId, Support, Word};

const MAX_BASIS: usize = 32_768;

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "dual bar calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

pub(crate) fn differential_grade(grade: Grade) -> Result<Grade> {
    let b = grade.b.checked_add(1).ok_or_else(degree_limit)?;
    Ok(Grade { b, ..grade })
}

pub(crate) fn cup_grade(left: Grade, right: Grade) -> Result<Grade> {
    if left.path != right.path || left.end != right.start {
        return Err(Error::plain(
            "rvPathCup requires intervals sharing a cut on the same bar path",
        ));
    }
    let b = left.b.checked_add(right.b).ok_or_else(degree_limit)?;
    Ok(Grade {
        end: right.end,
        b,
        ..left
    })
}

fn degree_limit() -> Error {
    Error::plain("dual bar cohomological degree is out of range")
}

pub(crate) fn pair(model: &mut Model, mut dual: &[Word], mut raw: &[Word]) -> Result<bool> {
    let mut result = false;
    while let (Some((a, tail_a)), Some((b, tail_b))) = (dual.split_first(), raw.split_first()) {
        model.spend(rv_bar::word_cells(a) + rv_bar::word_cells(b))?;
        match a.cmp(b) {
            Ordering::Less => dual = tail_a,
            Ordering::Greater => raw = tail_b,
            Ordering::Equal => {
                result = !result;
                dual = tail_a;
                raw = tail_b;
            }
        }
    }
    Ok(result)
}

pub(crate) fn cup_terms(model: &mut Model, left: &[Word], right: &[Word]) -> Result<Vec<Word>> {
    let count = left.len().saturating_mul(right.len());
    check_storage(count)?;
    let mut result = Vec::with_capacity(count);
    let mut cells = 0;
    for a in left {
        for b in right {
            let cost = rv_bar::word_cells(a) + rv_bar::word_cells(b);
            model.spend(cost)?;
            cells += cost;
            check_storage(cells)?;
            let mut word = a.to_vec();
            word.extend_from_slice(b);
            result.push(word.into());
        }
    }
    // A fixed cut determines the unique two factors. No left word is a
    // proper prefix of another at that cut, so nested sorted loops suffice.
    Ok(result)
}

struct Complex {
    path: PathId,
    ranks: Rc<Path>,
    expansions: BTreeMap<(usize, Factor), Support>,
    cells: usize,
}

impl Complex {
    fn new(model: &Model, path: PathId) -> Self {
        Self {
            path,
            ranks: model.rv_bar_path(path),
            expansions: BTreeMap::new(),
            cells: 0,
        }
    }

    fn expand(&mut self, model: &mut Model, start: usize, factor: Factor) -> Result<Support> {
        if factor.end - start < 2 {
            return Ok([].into());
        }
        let key = (start, factor);
        model.spend(1 + self.expansions.len().max(1).ilog2() as usize)?;
        if let Some(support) = self.expansions.get(&key) {
            return Ok(support.clone());
        }
        let b = factor.s.checked_add(2).ok_or_else(degree_limit)?;
        let grade = Grade {
            path: self.path,
            start,
            end: factor.end,
            b,
        };
        let mut expanded = Vec::new();
        for word in rv_bar::replacements(model, grade, MAX_BASIS)? {
            let image = self.ranks.multiply(model, start, &word, factor.s)?;
            if image.binary_search(&factor.generator).is_ok() {
                self.cells += rv_bar::word_cells(&word);
                check_storage(self.cells)?;
                expanded.push(word);
            }
        }
        self.cells += 5;
        check_storage(self.cells)?;
        let expanded: Support = expanded.into();
        self.expansions.insert(key, expanded.clone());
        Ok(expanded)
    }

    fn differential(
        &mut self,
        model: &mut Model,
        start: usize,
        dual: &[Word],
    ) -> Result<Vec<Word>> {
        let mut sum = BTreeSet::new();
        let mut cells = 0;
        for word in dual {
            let mut left = start;
            for (index, &factor) in word.iter().enumerate() {
                model.spend(1)?;
                // A transpose column expands just one of its supported
                // factors. All other letters and both endpoints are fixed.
                for piece in self.expand(model, left, factor)?.iter() {
                    let cost = 1 + 3 * (word.len() - 1 + piece.len());
                    model.spend(cost * (1 + sum.len().max(1).ilog2() as usize))?;
                    let mut source = Vec::with_capacity(word.len() - 1 + piece.len());
                    source.extend_from_slice(&word[..index]);
                    source.extend_from_slice(piece);
                    source.extend_from_slice(&word[index + 1..]);
                    let source: Word = source.into();
                    if sum.remove(&source) {
                        cells -= cost;
                    } else {
                        cells += cost;
                        check_storage(cells + self.cells)?;
                        sum.insert(source);
                    }
                }
                left = factor.end;
            }
        }
        Ok(sum.into_iter().collect())
    }
}

pub(crate) fn differential(model: &mut Model, grade: Grade, dual: &[Word]) -> Result<Support> {
    let terms = Complex::new(model, grade.path).differential(model, grade.start, dual)?;
    rv_bar::persist(model, terms)
}

#[derive(Default)]
pub(crate) struct Ext {
    terms: Vec<Word>,
    quotient: QuotientBasis,
}

impl Ext {
    pub(crate) fn build(model: &mut Model, grade: Grade) -> Result<Self> {
        let terms = rv_bar::basis(model, grade, grade.max_length(), MAX_BASIS)?;
        if terms.is_empty() {
            return Ok(Self::default());
        }
        let mut complex = Complex::new(model, grade.path);
        let mut boundaries = Basis::default();
        if let Ok(preceding) = grade.differential() {
            let incoming = rv_bar::basis(model, preceding, preceding.max_length(), MAX_BASIS)?;
            for term in &incoming {
                let image = complex.differential(model, grade.start, std::slice::from_ref(term))?;
                let coordinates = term_indices(model, &terms, &image)?;
                boundaries.insert(model, coordinates)?;
            }
        }
        let mut rows = Rows::new("path Ext");
        let quotient = QuotientBasis::build(model, boundaries, terms.len(), |model, i| {
            let image =
                complex.differential(model, grade.start, std::slice::from_ref(&terms[i]))?;
            let mut column = image
                .into_iter()
                .map(|word| rows.coordinate(model, rv_bar::word_cells(&word), word))
                .collect::<Result<Vec<_>>>()?;
            sort_support(model, &mut column)?;
            Ok(column)
        })?;
        model.allocate(terms.iter().map(|word| rv_bar::word_cells(word)).sum())?;
        Ok(Self { terms, quotient })
    }

    pub(crate) fn rank(&self) -> usize {
        self.quotient.rank()
    }

    fn project(&self, model: &mut Model, terms: &[Word]) -> Result<Vec<usize>> {
        let coordinates = term_indices(model, &self.terms, terms)?;
        self.quotient
            .project(model, coordinates)?
            .ok_or_else(|| Error::plain("rvPathClass requires a closed dual bar cochain"))
    }

    fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<Word>> {
        let coordinates = self.quotient.representative(model, class)?;
        Ok(coordinates
            .into_iter()
            .map(|i| self.terms[i].clone())
            .collect())
    }
}

fn term_indices(model: &mut Model, basis: &[Word], terms: &[Word]) -> Result<Vec<usize>> {
    let log = 1 + basis.len().max(1).ilog2() as usize;
    terms
        .iter()
        .map(|word| {
            model.spend(rv_bar::word_cells(word) * log)?;
            basis
                .binary_search(word)
                .map_err(|_| Error::plain("expected a word of this dual bar bidegree"))
        })
        .collect()
}

pub(crate) fn representative(model: &mut Model, ty: TypeId, class: &[usize]) -> Result<Support> {
    if class.is_empty() {
        return Ok([].into());
    }
    let words = model.rv_path_ext(ty)?.representative(model, class)?;
    rv_bar::persist(model, words)
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Word]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(cohomology::zero());
    }
    let class = model.rv_path_ext(ty)?.project(model, terms)?;
    model.allocate(class.len())?;
    Ok(Value::Class(class.into()))
}

pub(crate) fn cup(
    model: &mut Model,
    (left_ty, left): (TypeId, &[usize]),
    (right_ty, right): (TypeId, &[usize]),
    output: TypeId,
) -> Result<Value> {
    if left.is_empty() || right.is_empty() {
        return Ok(cohomology::zero());
    }
    let left = model.rv_path_ext(left_ty)?.representative(model, left)?;
    let right = model.rv_path_ext(right_ty)?.representative(model, right)?;
    let terms = cup_terms(model, &left, &right)?;
    class(model, output, &terms)
}
