use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::mixed::Complex;
use crate::mixed_page;
use crate::model::{Model, TypeId, Value};
use crate::odd_cobar::Prime;
use crate::odd_ext;
use crate::prime_linear::Coordinate;
use crate::rv_hh::{self, Differential, Grade, Word};
use crate::rv_hh_homology::{self, MAX_BASIS};

pub(crate) struct Origin {
    weight: Grade,
    words: Vec<Word>,
}

impl Origin {
    pub(crate) fn build(model: &mut Model, weight: Grade) -> Result<(Self, Complex)> {
        let words = rv_hh::weight_basis(model, &weight, MAX_BASIS)?;
        let ideal = model.rv_ideal_named(&weight.ideal).unwrap();
        let mut differential = Differential::new(model, &weight, &ideal)?;
        let mut retained = 0;
        let mut degrees = Vec::with_capacity(words.len());
        let mut down = Vec::with_capacity(words.len());
        let mut up = Vec::with_capacity(words.len());
        for word in &words {
            let singleton = std::slice::from_ref(word);
            let degree = word_degree(word)?;
            degrees.push(degree);
            let image = differential.apply(model, singleton)?;
            let connes = rv_hh::connes(model, &weight, singleton)?;
            let boundary = word_indices(model, &words, &image)?;
            let operator = word_indices(model, &words, &connes)?;
            retained += rv_hh::word_cells(word) + 2 * (boundary.len() + operator.len());
            if retained > MAX_PIVOT_CELLS {
                return Err(Error::plain(
                    "Hochschild mixed complex exceeds the materialization limit",
                ));
            }
            down.push(boundary);
            up.push(operator);
        }

        let complex = Complex::from_columns(model, Prime::validated(2), degrees, down, up)?;
        model.allocate(words.iter().map(|word| rv_hh::word_cells(word)).sum())?;
        Ok((Self { weight, words }, complex))
    }

    pub(crate) fn check_weight(&self, input: &Grade) -> Result<()> {
        if self.weight.ideal != input.ideal || self.weight.weight != input.weight {
            return Err(Error::plain(
                "rvHHToMixed requires the ideal and weight of the imported mixed complex",
            ));
        }
        Ok(())
    }

    pub(crate) fn original_grade(&self, b: i128) -> Result<Grade> {
        let b = i64::try_from(b).map_err(|_| degree_limit())?;
        Ok(Grade {
            b,
            ..self.weight.clone()
        })
    }

    pub(crate) fn chain_to_mixed(&self, model: &mut Model, words: &[Word]) -> Result<Value> {
        let support = word_indices(model, &self.words, words)?;
        odd_ext::persist(model, support)
    }

    pub(crate) fn chain_to_hh(&self, model: &mut Model, support: &[Coordinate]) -> Result<Value> {
        let terms = self.original_words(model, support)?;
        Ok(Value::RVHHChain(rv_hh::persist(model, terms)?))
    }

    pub(crate) fn class_to_mixed(
        &self,
        model: &mut Model,
        input: TypeId,
        output: TypeId,
        support: &[usize],
    ) -> Result<Value> {
        if support.is_empty() {
            return Ok(odd_ext::zero());
        }
        let lift = model
            .rv_hh_homology(input)?
            .representative(model, support)?;
        let terms = word_indices(model, &self.words, &lift)?;
        mixed_page::class(model, output, &terms)
    }

    pub(crate) fn class_to_hh(
        &self,
        model: &mut Model,
        input: TypeId,
        output: TypeId,
        support: &[Coordinate],
    ) -> Result<Value> {
        let Value::OddClass(terms) = mixed_page::representative(model, input, support)? else {
            unreachable!()
        };
        let terms = self.original_words(model, &terms)?;
        rv_hh_homology::class(model, output, &terms)
    }

    fn original_words(&self, model: &mut Model, support: &[Coordinate]) -> Result<Vec<Word>> {
        model.spend(
            support
                .iter()
                .map(|&(i, _)| rv_hh::word_cells(&self.words[i]))
                .sum(),
        )?;
        Ok(support
            .iter()
            .map(|&(i, coefficient)| {
                debug_assert_eq!(coefficient, 1);
                self.words[i].clone()
            })
            .collect())
    }
}

fn degree_limit() -> Error {
    Error::plain("Hochschild mixed complex degree is out of range")
}

fn word_degree(word: &Word) -> Result<i64> {
    let degree = word
        .iter()
        .try_fold(word.len() - 1, |b, factor| b.checked_add(factor.s))
        .ok_or_else(degree_limit)?;
    i64::try_from(degree).map_err(|_| degree_limit())
}

fn word_indices(model: &mut Model, basis: &[Word], words: &[Word]) -> Result<Vec<Coordinate>> {
    let comparisons = 1 + basis.len().max(1).ilog2() as usize;
    words
        .iter()
        .map(|word| {
            model.spend(rv_hh::word_cells(word) * comparisons)?;
            let index = basis.binary_search(word).map_err(|_| {
                Error::plain("Hochschild word lies outside the imported weight basis")
            })?;
            Ok((index, 1))
        })
        .collect()
}
