use std::cell::OnceCell;
use std::cmp::Ordering;
use std::collections::{btree_map::Entry, BTreeMap};

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::odd_cobar::{self, word_cells, Prime, Term, Word};
use crate::prime_linear::{add_scaled, Basis, Coordinate, QuotientBasis};

const MAX_BASIS_WORDS: usize = 32_768;

pub(crate) struct OddExt {
    prime: Prime,
    s: usize,
    t: usize,
    words: Vec<Word>,
    quotient: QuotientBasis,
    preimages: OnceCell<BoundaryPreimages>,
}

struct BoundaryPreimages {
    words: Vec<Word>,
    images: Basis,
}

impl OddExt {
    pub(crate) fn build(model: &mut Model, prime: Prime, s: usize, t: usize) -> Result<Self> {
        let words = odd_cobar::basis(model, prime, s, t, MAX_BASIS_WORDS)?;
        let quotient = if words.is_empty() {
            QuotientBasis::default()
        } else {
            let mut boundaries = Basis::default();
            boundary_columns(model, prime, s, t, &words, |model, _, column| {
                boundaries.insert(model, prime, column)
            })?;

            let mut rows = WordRows::default();
            QuotientBasis::build(model, prime, boundaries, words.len(), |model, index| {
                let image = odd_cobar::differential_word(model, prime, &words[index])?;
                rows.encode(model, image)
            })?
        };
        Ok(Self {
            prime,
            s,
            t,
            words,
            quotient,
            preimages: OnceCell::new(),
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.quotient.rank()
    }

    pub(crate) fn project(&self, model: &mut Model, terms: &[Term]) -> Result<Vec<Coordinate>> {
        let column = word_indices(model, &self.words, terms)?;
        self.quotient
            .project(model, self.prime, column)?
            .ok_or_else(|| Error::plain("extClass requires a closed odd-prime cobar cochain"))
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        class: &[Coordinate],
    ) -> Result<Vec<Term>> {
        let support = self.quotient.representative(model, self.prime, class)?;
        Ok(decode(&self.words, support))
    }

    pub(crate) fn primitive(&self, model: &mut Model, terms: &[Term]) -> Result<Option<Vec<Term>>> {
        if terms.is_empty() {
            return Ok(Some(Vec::new()));
        }
        let column = word_indices(model, &self.words, terms)?;
        if !self
            .quotient
            .is_boundary(model, self.prime, column.clone())?
        {
            return Ok(None);
        }
        let preimages = if let Some(preimages) = self.preimages.get() {
            preimages
        } else {
            let mut images = Basis::default();
            let words = boundary_columns(
                model,
                self.prime,
                self.s,
                self.t,
                &self.words,
                |model, index, differential| {
                    images.insert_with_witness(
                        model,
                        self.prime,
                        differential,
                        vec![(index, 1)],
                    )?;
                    Ok(())
                },
            )?;
            model.allocate(images.cells())?;
            self.preimages
                .get_or_init(|| BoundaryPreimages { words, images })
        };
        let (remainder, support) = preimages
            .images
            .reduce_with_witness(model, self.prime, column)?;
        debug_assert!(remainder.is_empty());
        Ok(Some(decode(&preimages.words, support)))
    }
}

fn boundary_columns(
    model: &mut Model,
    prime: Prime,
    s: usize,
    t: usize,
    target: &[Word],
    mut process: impl FnMut(&mut Model, usize, Vec<Coordinate>) -> Result<()>,
) -> Result<Vec<Word>> {
    let Some(previous) = s.checked_sub(1) else {
        return Ok(Vec::new());
    };
    let source = odd_cobar::basis(model, prime, previous, t, MAX_BASIS_WORDS)?;
    for (index, word) in source.iter().enumerate() {
        let terms = odd_cobar::differential_word(model, prime, word)?;
        let column = word_indices(model, target, &terms)?;
        process(model, index, column)?;
    }
    Ok(source)
}

pub(crate) fn word_indices(
    model: &mut Model,
    basis: &[Word],
    terms: &[Term],
) -> Result<Vec<Coordinate>> {
    let search_steps = 1 + basis.len().max(1).ilog2() as usize;
    terms
        .iter()
        .map(|(word, coefficient)| {
            model.spend(search_steps * word_cells(word))?;
            let index = basis
                .binary_search(word)
                .map_err(|_| Error::plain("odd-prime word is outside its cobar bidegree"))?;
            Ok((index, *coefficient))
        })
        .collect()
}

pub(crate) fn decode(words: &[Word], coordinates: Vec<Coordinate>) -> Vec<Term> {
    coordinates
        .into_iter()
        .map(|(index, coefficient)| (words[index].clone(), coefficient))
        .collect()
}

// Only outgoing rows that occur in a differential need to be indexed.
#[derive(Default)]
pub(crate) struct WordRows {
    indices: BTreeMap<Word, usize>,
    cells: usize,
}

impl WordRows {
    pub(crate) fn encode(
        &mut self,
        model: &mut Model,
        terms: Vec<Term>,
    ) -> Result<Vec<Coordinate>> {
        let mut column = Vec::with_capacity(terms.len());
        for (word, coefficient) in terms {
            let size = 1 + word_cells(&word);
            let next = self.indices.len();
            model.spend(size * (1 + next.max(1).ilog2() as usize))?;
            let index = match self.indices.entry(word) {
                Entry::Occupied(entry) => *entry.get(),
                Entry::Vacant(entry) => {
                    self.cells += size;
                    if self.cells > MAX_PIVOT_CELLS {
                        return Err(Error::plain(
                            "odd-prime cobar reduction exceeds the row storage limit",
                        ));
                    }
                    *entry.insert(next)
                }
            };
            column.push((index, coefficient));
        }
        model.spend(column.len() * (1 + column.len().max(1).ilog2() as usize))?;
        column.sort_unstable();
        Ok(column)
    }
}

pub(crate) fn zero() -> Value {
    Value::OddClass([].into())
}

pub(crate) fn persist(model: &mut Model, coordinates: Vec<Coordinate>) -> Result<Value> {
    model.allocate(2 * coordinates.len())?;
    Ok(Value::OddClass(coordinates.into()))
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(zero());
    }
    let coordinates = model.odd_ext(ty)?.project(model, terms)?;
    persist(model, coordinates)
}

pub(crate) fn add(
    model: &mut Model,
    prime: Prime,
    left: &[Coordinate],
    right: &[Coordinate],
) -> Result<Value> {
    let coordinates = add_scaled(model, prime, left, right, 1)?;
    persist(model, coordinates)
}

pub(crate) fn scale(
    model: &mut Model,
    prime: Prime,
    scalar: usize,
    class: &[Coordinate],
) -> Result<Value> {
    let coordinates = add_scaled(model, prime, &[], class, prime.reduce(scalar))?;
    persist(model, coordinates)
}

pub(crate) fn elements(
    model: &mut Model,
    prime: Prime,
    rank: usize,
    limit: usize,
) -> Result<Vec<Value>> {
    let p = prime.get() as usize;
    let mut size = 1;
    for _ in 0..rank {
        if size > limit / p {
            return Err(Error::plain(format!(
                "type exceeds the enumeration limit of {limit} elements"
            )));
        }
        size *= p;
    }
    model.allocate(size)?;
    model.spend(size * (rank + 1))?;
    let mut values = Vec::with_capacity(size);
    for mut index in 0..size {
        let mut coordinates = Vec::new();
        for generator in 0..rank {
            let coefficient = (index % p) as u32;
            if coefficient != 0 {
                coordinates.push((generator, coefficient));
            }
            index /= p;
        }
        values.push(persist(model, coordinates)?);
    }
    Ok(values)
}

pub(crate) fn compare(
    model: &mut Model,
    left: &[Coordinate],
    right: &[Coordinate],
) -> Result<Ordering> {
    model.spend(1 + left.len().min(right.len()))?;
    Ok(left.iter().rev().cmp(right.iter().rev()))
}

pub(crate) fn yoneda(
    model: &mut Model,
    left_ty: TypeId,
    left: &[Coordinate],
    right_ty: TypeId,
    right: &[Coordinate],
    target: TypeId,
) -> Result<Value> {
    if left.is_empty() || right.is_empty() {
        return Ok(zero());
    }
    let left_group = model.odd_ext(left_ty)?;
    let left = left_group.representative(model, left)?;
    let right = model.odd_ext(right_ty)?.representative(model, right)?;
    let terms = odd_cobar::concat_support(model, left_group.prime, &left, &right)?;
    class(model, target, &terms)
}
