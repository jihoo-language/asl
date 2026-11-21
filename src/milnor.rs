use std::cell::OnceCell;
use std::collections::HashMap;
use std::rc::Rc;

use crate::cobar::{self, Monomial, Word as CobarWord};
use crate::cohomology;
use crate::error::{Error, Result};
use crate::linear::Basis;
use crate::model::{Model, TypeId, Value};
use crate::steenrod::{self, Word};

const MAX_BASIS_WORDS: usize = 512;
const MAX_SPLIT_CELLS: usize = 300_000;

pub(crate) struct MilnorBasis {
    cobar_words: Vec<CobarWord>,
    admissible_words: Vec<Word>,
    pairing: Basis,
    operations: Vec<OnceCell<Rc<[Word]>>>,
}

impl MilnorBasis {
    pub(crate) fn build(model: &mut Model, degree: usize) -> Result<Self> {
        let cobar_words = cobar::basis(model, usize::from(degree != 0), degree, MAX_BASIS_WORDS)?;
        let mut admissible_words = Vec::with_capacity(cobar_words.len());
        for word in &cobar_words {
            let Some(exponents) = word.first() else {
                admissible_words.push([].into());
                continue;
            };
            // The excesses of an admissible word recover the Milnor
            // exponents via r_k = i_k - 2*i_(k+1).
            model.spend(exponents.len())?;
            let mut indices = Vec::with_capacity(exponents.len());
            let mut degree = 0usize;
            for &exponent in exponents.iter().rev() {
                degree = degree
                    .checked_mul(2)
                    .and_then(|twice| twice.checked_add(exponent))
                    .ok_or_else(degree_overflow)?;
                indices.push(degree);
            }
            indices.reverse();
            admissible_words.push(indices.into());
        }
        let admissible_cells = admissible_words
            .iter()
            .map(|word: &Word| 1 + word.len())
            .sum::<usize>();
        model.spend(admissible_cells * (1 + cobar_words.len().ilog2() as usize))?;
        admissible_words.sort_unstable();

        let mut pairing = Basis::default();
        let mut evaluator = Pairing::default();
        for (column, word) in admissible_words.iter().enumerate() {
            let mut image = Vec::new();
            for (row, monomial) in cobar_words.iter().enumerate() {
                if evaluator.evaluate(model, word, monomial)? {
                    image.push(row);
                }
            }
            if pairing
                .insert_with_witness(model, image, vec![column])?
                .is_some()
            {
                return Err(Error::plain("Milnor pairing matrix is singular"));
            }
        }
        model.allocate(1 + admissible_cells + pairing.cells() + cobar_words.len())?;
        let operations = (0..cobar_words.len()).map(|_| OnceCell::new()).collect();
        Ok(Self {
            cobar_words,
            admissible_words,
            pairing,
            operations,
        })
    }

    fn operation(&self, model: &mut Model, index: usize) -> Result<Rc<[Word]>> {
        let cached = &self.operations[index];
        if let Some(operation) = cached.get() {
            return Ok(operation.clone());
        }
        let (remainder, coordinates) =
            self.pairing
                .reduce_with_witness(model, vec![index], Vec::new())?;
        if !remainder.is_empty() {
            return Err(Error::plain("Milnor pairing matrix is singular"));
        }
        model.spend(coordinates.len())?;
        model.allocate(coordinates.len())?;
        let operation = coordinates
            .into_iter()
            .map(|coordinate| self.admissible_words[coordinate].clone())
            .collect();
        Ok(cached.get_or_init(|| operation).clone())
    }
}

fn degree_overflow() -> Error {
    Error::plain("milnor degree is out of range")
}

pub(crate) fn literal(model: &mut Model, mut exponents: Vec<usize>) -> Result<(usize, Value)> {
    model.spend(1 + exponents.len())?;
    let Some(last) = exponents.iter().rposition(|&exponent| exponent != 0) else {
        return Ok((0, steenrod::square(model, 0)?));
    };
    exponents.truncate(last + 1);
    if let [degree] = exponents.as_slice() {
        return Ok((*degree, steenrod::square(model, *degree)?));
    }

    let word: CobarWord = [exponents.into()].into();
    let (_, degree) = cobar::word_metadata(&word).map_err(|_| degree_overflow())?;
    let basis = model.milnor_basis(degree)?;
    model.spend(cobar::word_cells(&word) * (1 + basis.cobar_words.len().ilog2() as usize))?;
    let index = basis
        .cobar_words
        .binary_search(&word)
        .map_err(|_| Error::plain("Milnor basis lacks a monomial in its degree"))?;
    Ok((degree, Value::Operation(basis.operation(model, index)?)))
}

pub(crate) fn pair(model: &mut Model, words: &[Word], cochains: &[CobarWord]) -> Result<bool> {
    if words.is_empty() || cochains.is_empty() {
        return Ok(false);
    }
    let mut evaluator = Pairing::default();
    let mut value = false;
    for word in words {
        for cochain in cochains {
            value ^= evaluator.evaluate(model, word, cochain)?;
        }
    }
    Ok(value)
}

#[derive(Default)]
struct Pairing {
    splits: HashMap<Monomial, Rc<[CobarWord]>>,
    split_cells: usize,
}

impl Pairing {
    fn evaluate(&mut self, model: &mut Model, squares: &[usize], word: &CobarWord) -> Result<bool> {
        model.spend(1)?;
        match word.as_ref() {
            [] => Ok(squares.is_empty()),
            [monomial] => self.monomial(model, squares, monomial),
            _ => Err(Error::plain(
                "Milnor pairing requires a single cobar factor or the unit",
            )),
        }
    }

    fn monomial(
        &mut self,
        model: &mut Model,
        squares: &[usize],
        monomial: &Monomial,
    ) -> Result<bool> {
        model.spend(1)?;
        match squares {
            [] => Ok(false),
            [index] => Ok(monomial.as_ref() == [*index]),
            [first, rest @ ..] => {
                let splits = self.square_left_terms(model, monomial)?;
                let search_steps = 1 + splits.len().checked_ilog2().unwrap_or(0) as usize;
                model.spend(search_steps)?;
                let start = splits.partition_point(|term| term[0][0] < *first);
                let mut value = false;
                for term in &splits[start..] {
                    if term[0][0] != *first {
                        break;
                    }
                    // The left coproduct factor pairs with the leftmost
                    // square. The remaining squares act on the right.
                    value ^= self.monomial(model, rest, &term[1])?;
                }
                Ok(value)
            }
        }
    }

    fn square_left_terms(
        &mut self,
        model: &mut Model,
        monomial: &Monomial,
    ) -> Result<Rc<[CobarWord]>> {
        model.spend(1 + monomial.len())?;
        if let Some(cached) = self.splits.get(monomial) {
            return Ok(cached.clone());
        }
        // Both inputs to the product pairing have positive degree, so
        // only the reduced terms of the existing cobar differential enter.
        let word: CobarWord = [monomial.clone()].into();
        let mut terms = cobar::differential_word(model, &word)?;
        model.spend(terms.len())?;
        terms.retain(|term| term[0].len() == 1);
        let cells = 1 + monomial.len() + terms.iter().map(cobar::word_cells).sum::<usize>();
        if self.split_cells + cells > MAX_SPLIT_CELLS {
            return Err(Error::plain(
                "Milnor pairing exceeds the coproduct storage limit",
            ));
        }
        self.split_cells += cells;
        let terms: Rc<[CobarWord]> = terms.into();
        self.splits.insert(monomial.clone(), terms.clone());
        Ok(terms)
    }
}

pub(crate) fn coact_bidegree(
    model: &mut Model,
    base: TypeId,
    n: usize,
    index: usize,
    q: usize,
) -> Result<(usize, usize)> {
    let difference = n.checked_sub(q).ok_or_else(|| {
        Error::plain(format!("coact target degree {n} is below input degree {q}"))
    })?;
    let target = model.cohomology_type(base, n)?;
    let rank = model.cohomology(target)?.rank();
    if index >= rank {
        return Err(Error::plain(format!(
            "coact generator index {index} is out of range for {} of rank {rank}",
            model.type_name(target)
        )));
    }
    Ok((usize::from(difference != 0), difference))
}

pub(crate) fn coact(
    model: &mut Model,
    base: TypeId,
    n: usize,
    index: usize,
    q: usize,
    class: Rc<[usize]>,
) -> Result<Vec<CobarWord>> {
    let (_, difference) = coact_bidegree(model, base, n, index, q)?;
    if difference == 0 {
        if has_coordinate(model, &class, index)? {
            model.allocate(1)?;
            return Ok(vec![[].into()]);
        }
        return Ok(Vec::new());
    }
    // Positive operations annihilate all degree-zero classes by instability.
    if q == 0 || class.is_empty() {
        return Ok(Vec::new());
    }
    let basis = model.milnor_basis(difference)?;
    let mut support = Vec::new();
    for (coordinate, monomial) in basis.cobar_words.iter().enumerate() {
        let operation = basis.operation(model, coordinate)?;
        let Value::Class(image) = cohomology::act(model, base, q, &operation, class.clone())?
        else {
            unreachable!()
        };
        if has_coordinate(model, &image, index)? {
            support.push(monomial.clone());
        }
    }
    model.allocate(support.len())?;
    Ok(support)
}

fn has_coordinate(model: &mut Model, class: &[usize], index: usize) -> Result<bool> {
    model.spend(1 + class.len().checked_ilog2().unwrap_or(0) as usize)?;
    Ok(class.binary_search(&index).is_ok())
}
