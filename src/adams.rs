use crate::cobar::{self, Word, WordRows};
use crate::error::{Error, Result};
use crate::linear::{Basis, QuotientBasis};
use crate::model::{Model, TypeId, Value};
use crate::profile::Profile;
use std::cell::OnceCell;

const MAX_BASIS_WORDS: usize = 32_768;

#[derive(Default)]
pub(crate) struct Adams {
    bidegree: (Profile, usize, usize),
    words: Vec<Word>,
    quotient: QuotientBasis,
    boundary_preimages: OnceCell<BoundaryPreimages>,
}

struct BoundaryPreimages {
    words: Vec<Word>,
    images: Basis,
}

impl Adams {
    pub(crate) fn build(model: &mut Model, profile: Profile, s: usize, t: usize) -> Result<Self> {
        let words = profile.basis(model, s, t, MAX_BASIS_WORDS)?;
        if words.is_empty() {
            return Ok(Self {
                bidegree: (profile, s, t),
                ..Self::default()
            });
        }

        let mut boundaries = Basis::default();
        for_each_boundary(model, profile, s, t, &words, |model, _, column| {
            boundaries.insert(model, column)?;
            Ok(())
        })?;

        let mut rows = WordRows::default();
        let quotient = QuotientBasis::build(model, boundaries, words.len(), |model, index| {
            let mut column = profile
                .differential_word(model, &words[index])?
                .into_iter()
                .map(|term| rows.coordinate(model, term))
                .collect::<Result<Vec<_>>>()?;
            model.spend(column.len() * (1 + column.len().max(1).ilog2() as usize))?;
            column.sort_unstable();
            Ok(column)
        })?;
        Ok(Self {
            bidegree: (profile, s, t),
            words,
            quotient,
            boundary_preimages: OnceCell::new(),
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.quotient.rank()
    }

    pub(crate) fn project(&self, model: &mut Model, words: &[Word]) -> Result<Vec<usize>> {
        let column = word_indices(model, &self.words, words)?;
        self.quotient
            .project(model, column)?
            .ok_or_else(|| Error::plain("extClass requires a closed cobar cochain"))
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<Word>> {
        Ok(self
            .quotient
            .representative(model, class)?
            .into_iter()
            .map(|index| self.words[index].clone())
            .collect())
    }

    pub(crate) fn primitive(&self, model: &mut Model, words: &[Word]) -> Result<Option<Vec<Word>>> {
        if words.is_empty() {
            return Ok(Some(Vec::new()));
        }
        let column = word_indices(model, &self.words, words)?;
        if !self.quotient.is_boundary(model, column.clone())? {
            return Ok(None);
        }
        let preimages = if let Some(preimages) = self.boundary_preimages.get() {
            preimages
        } else {
            let (profile, s, t) = self.bidegree;
            let mut images = Basis::default();
            let words = for_each_boundary(
                model,
                profile,
                s,
                t,
                &self.words,
                |model, index, differential| {
                    let _ = images.insert_with_witness(model, differential, vec![index])?;
                    Ok(())
                },
            )?;
            model.allocate(images.cells())?;
            self.boundary_preimages
                .get_or_init(|| BoundaryPreimages { words, images })
        };
        let (remainder, support) =
            preimages
                .images
                .reduce_with_witness(model, column, Vec::new())?;
        debug_assert!(remainder.is_empty());
        Ok(Some(
            support
                .into_iter()
                .map(|index| preimages.words[index].clone())
                .collect(),
        ))
    }
}

fn for_each_boundary(
    model: &mut Model,
    profile: Profile,
    s: usize,
    t: usize,
    target: &[Word],
    mut use_column: impl FnMut(&mut Model, usize, Vec<usize>) -> Result<()>,
) -> Result<Vec<Word>> {
    let Some(previous) = s.checked_sub(1) else {
        return Ok(Vec::new());
    };
    let source = profile.basis(model, previous, t, MAX_BASIS_WORDS)?;
    for (index, word) in source.iter().enumerate() {
        let differential = profile.differential_word(model, word)?;
        let column = word_indices(model, target, &differential)?;
        use_column(model, index, column)?;
    }
    Ok(source)
}

fn word_indices(model: &mut Model, basis: &[Word], words: &[Word]) -> Result<Vec<usize>> {
    let search_steps = 1 + basis.len().max(1).ilog2() as usize;
    words
        .iter()
        .map(|word| {
            model.spend(search_steps * cobar::word_cells(word))?;
            basis
                .binary_search(word)
                .map_err(|_| Error::plain("word is outside its cobar bidegree"))
        })
        .collect()
}

pub(crate) fn class(model: &mut Model, ty: TypeId, words: &[Word]) -> Result<Value> {
    if words.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let adams = model.adams(ty)?;
    let coordinates = adams.project(model, words)?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

pub(crate) fn yoneda(
    model: &mut Model,
    left_ty: TypeId,
    left: &[usize],
    right_ty: TypeId,
    right: &[usize],
    ty: TypeId,
) -> Result<Value> {
    if left.is_empty() || right.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let left = model.adams(left_ty)?.representative(model, left)?;
    let right = model.adams(right_ty)?.representative(model, right)?;
    let words = cobar::concat_support(model, &left, &right)?;
    class(model, ty, &words)
}
