use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;

use super::{Chain, Field, Spec};

pub(super) const WORD_LIMIT: usize = 64;
const TERM_LIMIT: usize = 10_000;

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(super) struct Word {
    pub(super) start: usize,
    pub(super) letters: Vec<u8>,
}

pub(super) type Terms = BTreeMap<Word, u32>;

pub(super) fn display(terms: &Terms) -> String {
    let words: Vec<_> = terms
        .iter()
        .map(|(word, &c)| (word.start, &word.letters, c))
        .collect();
    format!("{words:?}")
}

pub(super) fn word_limit() -> Error {
    Error::plain("partition Ginzburg word length exceeds 64")
}

pub(super) fn accumulate(
    field: &Field,
    terms: &mut Terms,
    word: Word,
    coefficient: u32,
) -> Result<()> {
    if coefficient == 0 {
        return Ok(());
    }
    if word.letters.len() > WORD_LIMIT {
        return Err(word_limit());
    }
    let value = field.add(terms.get(&word).copied().unwrap_or(0), coefficient);
    if value == 0 {
        terms.remove(&word);
    } else {
        terms.insert(word, value);
        if terms.len() > TERM_LIMIT {
            return Err(Error::plain(
                "partition Ginzburg support exceeds 10,000 terms",
            ));
        }
    }
    Ok(())
}

impl Spec {
    pub(super) fn validate(&self, model: &mut Model, word: &Word) -> Result<i128> {
        let mut at = word.start;
        let mut degree = 0;
        model.spend(1 + word.letters.len())?;
        for &letter in &word.letters {
            let edge = self.symbols[usize::from(letter)];
            if edge.source != at {
                return Err(Error::plain(
                    "partition Ginzburg path has noncomposable symbols",
                ));
            }
            at = edge.target;
            degree += edge.degree;
        }
        Ok(degree)
    }

    fn end(&self, word: &Word) -> usize {
        word.letters
            .last()
            .map_or(word.start, |&a| self.symbols[usize::from(a)].target)
    }

    pub(crate) fn add(&self, model: &mut Model, a: &Chain, b: &Chain) -> Result<Chain> {
        let mut terms = a.terms.clone();
        model.spend(a.terms.len() + b.terms.len())?;
        for (w, &c) in &b.terms {
            accumulate(&self.hall.field, &mut terms, w.clone(), c)?;
        }
        Ok(self.chain(a.degree, terms))
    }

    pub(crate) fn scale(
        &self,
        model: &mut Model,
        chain: &Chain,
        coefficient: u128,
    ) -> Result<Chain> {
        let field = &self.hall.field;
        if coefficient >= u128::from(field.order()) {
            return Err(Error::plain(
                "partition Ginzburg coefficient is outside the finite field codes",
            ));
        }
        model.spend(chain.terms.len())?;
        let terms = if coefficient == 0 {
            Terms::new()
        } else {
            chain
                .terms
                .iter()
                .map(|(w, &c)| (w.clone(), field.mul(c, coefficient as u32)))
                .collect()
        };
        Ok(self.chain(chain.degree, terms))
    }

    pub(crate) fn negate(&self, model: &mut Model, chain: &Chain) -> Result<Chain> {
        self.scale(model, chain, u128::from(self.hall.field.sub(0, 1)))
    }

    pub(crate) fn product(&self, model: &mut Model, a: &Chain, b: &Chain) -> Result<Chain> {
        let field = &self.hall.field;
        let degree = a.degree.checked_add(b.degree).ok_or_else(degree_bound)?;
        let mut terms = Terms::new();
        for (x, &c) in &a.terms {
            let end = self.end(x);
            for (y, &d) in &b.terms {
                model.spend(1)?;
                if end == y.start {
                    model.spend(x.letters.len() + y.letters.len())?;
                    let word = Word {
                        start: x.start,
                        letters: x.letters.iter().chain(&y.letters).copied().collect(),
                    };
                    accumulate(field, &mut terms, word, field.mul(c, d))?;
                }
            }
        }
        Ok(self.chain(degree, terms))
    }

    pub(crate) fn differentiate(&self, model: &mut Model, chain: &Chain) -> Result<Chain> {
        let field = &self.hall.field;
        let degree = chain.degree.checked_add(1).ok_or_else(degree_bound)?;
        let mut terms = Terms::new();
        for (word, &c) in &chain.terms {
            let mut negative = false;
            for (i, &symbol) in word.letters.iter().enumerate() {
                model.spend(1)?;
                let symbol = usize::from(symbol);
                let coefficient = if negative { field.sub(0, c) } else { c };
                for (child, &value) in &self.differential[symbol] {
                    model.spend(word.letters.len() + child.letters.len())?;
                    let letters = word.letters[..i]
                        .iter()
                        .chain(&child.letters)
                        .chain(&word.letters[i + 1..])
                        .copied()
                        .collect();
                    accumulate(
                        field,
                        &mut terms,
                        Word {
                            start: word.start,
                            letters,
                        },
                        field.mul(coefficient, value),
                    )?;
                }
                negative ^= self.symbols[symbol].degree == -1;
            }
        }
        Ok(self.chain(degree, terms))
    }

    pub(crate) fn trace_report(&self, model: &mut Model, chain: &Chain) -> Result<String> {
        let field = &self.hall.field;
        let mut result = Terms::new();
        for (word, &coefficient) in &chain.terms {
            if self.end(word) != word.start {
                continue;
            }
            let odd = word
                .letters
                .iter()
                .filter(|&&a| self.symbols[usize::from(a)].degree == -1)
                .count();
            let mut minimum = word.clone();
            let mut c = coefficient;
            let mut cancelled = false;
            let mut prefix = 0;
            for (i, &a) in word.letters.iter().enumerate() {
                model.spend(1 + word.letters.len())?;
                let rotated = Word {
                    start: self.symbols[usize::from(a)].source,
                    letters: word.letters[i..]
                        .iter()
                        .chain(&word.letters[..i])
                        .copied()
                        .collect(),
                };
                let signed = if (prefix * (odd - prefix)).is_multiple_of(2) {
                    coefficient
                } else {
                    field.sub(0, coefficient)
                };
                match rotated.cmp(&minimum) {
                    std::cmp::Ordering::Less => {
                        minimum = rotated;
                        c = signed;
                        cancelled = false;
                    }
                    std::cmp::Ordering::Equal => cancelled |= c != signed,
                    _ => (),
                }
                prefix += usize::from(self.symbols[usize::from(a)].degree == -1);
            }
            if !cancelled {
                accumulate(field, &mut result, minimum, c)?;
            }
        }
        Ok(format!(
            "rvPartitionDGTrace[{}] = {}",
            self.name,
            display(&result)
        ))
    }
}

pub(crate) fn degree_bound() -> Error {
    Error::plain("partition Ginzburg degree overflows signed 128 bits")
}
