use std::borrow::Cow;
use std::collections::{hash_map::Entry, HashMap};

use crate::adams;
use crate::cobar::{self, Monomial, Support, Word};
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::profile::Profile;

pub(crate) fn cup_bidegree(
    i: usize,
    (s, t): (usize, usize),
    (u, v): (usize, usize),
) -> Result<(usize, usize)> {
    let length = (s as u128 + u as u128)
        .checked_sub(i as u128)
        .ok_or_else(|| Error::plain("cCup produces a negative cobar length"))?;
    let length =
        usize::try_from(length).map_err(|_| Error::plain("cCup cobar length is out of range"))?;
    let degree = t.checked_add(v).ok_or_else(degree_overflow)?;
    Ok((length, degree))
}

pub(crate) fn square_bidegree(k: usize, s: usize, t: usize) -> Result<(usize, usize)> {
    let length = s
        .checked_add(k)
        .ok_or_else(|| Error::plain("extSq cohomological degree is out of range"))?;
    let degree = t.checked_mul(2).ok_or_else(degree_overflow)?;
    Ok((length, degree))
}

fn degree_overflow() -> Error {
    Error::plain("cobar internal degree is out of range")
}

pub(crate) fn cup(
    model: &mut Model,
    profile: Profile,
    i: usize,
    left: &[Word],
    right: &[Word],
) -> Result<Value> {
    let terms = cup_support(model, profile, i, left, right)?;
    model.allocate(terms.iter().map(cobar::word_cells).sum())?;
    Ok(Value::Cobar(terms.into()))
}

pub(crate) fn cup_support(
    model: &mut Model,
    profile: Profile,
    i: usize,
    left: &[Word],
    right: &[Word],
) -> Result<Vec<Word>> {
    let (Some(a), Some(b)) = (left.first(), right.first()) else {
        return Ok(Vec::new());
    };
    if i == 0 {
        return cobar::concat_temporary(model, left, right);
    }
    let (length, _) = cup_bidegree(i, cobar::word_metadata(a)?, cobar::word_metadata(b)?)?;
    if i > a.len().min(b.len()) {
        return Ok(Vec::new());
    }
    let vertices = length
        .checked_add(1)
        .ok_or_else(|| Error::plain("cCup vertex count is out of range"))?;
    cobar::check_storage(length.saturating_mul(2).saturating_add(1))?;
    let size = length - i;
    let first_count = length - a.len();
    let mut subset: Vec<_> = (0..size).collect();
    let mut first = Vec::with_capacity(size);
    let mut second = Vec::with_capacity(size);
    let mut cofaces = Cofaces::new(profile);
    let mut result = Support::new();
    loop {
        model.spend(1 + size)?;
        first.clear();
        second.clear();
        for (j, &vertex) in subset.iter().enumerate() {
            // j is one-based in the coface parity convention.
            if (vertex + j + 1) & 1 == 0 {
                first.push(vertex);
            } else {
                second.push(vertex);
            }
        }
        if first.len() == first_count {
            let a = cofaces.lift(model, left, &first)?;
            let b = cofaces.lift(model, right, &second)?;
            for left in a.iter() {
                for right in b.iter() {
                    if let Some(word) = multiply_words(model, profile, left, right)? {
                        let cells = cobar::word_cells(&word);
                        result.toggle(model, word, cells)?;
                    }
                }
            }
        }
        if !next_subset(&mut subset, vertices) {
            return Ok(result.into_vec());
        }
    }
}

fn next_subset(subset: &mut [usize], vertices: usize) -> bool {
    let size = subset.len();
    let Some(i) = (0..size).rev().find(|&i| subset[i] < vertices - size + i) else {
        return false;
    };
    subset[i] += 1;
    for j in i + 1..size {
        subset[j] = subset[j - 1] + 1;
    }
    true
}

type Coproduct = Vec<(Monomial, Monomial)>;

struct Cofaces {
    profile: Profile,
    unit: Monomial,
    coproducts: HashMap<Monomial, Coproduct>,
    cells: usize,
}

impl Cofaces {
    fn new(profile: Profile) -> Self {
        Self {
            profile,
            unit: [].into(),
            coproducts: HashMap::new(),
            cells: 0,
        }
    }

    fn coproduct(&mut self, model: &mut Model, letter: &Monomial) -> Result<&Coproduct> {
        model.spend(1 + letter.len())?;
        Ok(match self.coproducts.entry(letter.clone()) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => {
                let terms = cobar::coproduct(model, self.profile, letter)?;
                self.cells += 1
                    + letter.len()
                    + terms
                        .iter()
                        .map(|(left, right)| 1 + left.len() + right.len())
                        .sum::<usize>();
                cobar::check_storage(self.cells)?;
                entry.insert(terms)
            }
        })
    }

    fn lift<'a>(
        &mut self,
        model: &mut Model,
        words: &'a [Word],
        indices: &[usize],
    ) -> Result<Cow<'a, [Word]>> {
        let mut words = Cow::Borrowed(words);
        for &i in indices {
            let mut next = Support::new();
            for word in words.iter() {
                let cells = cobar::word_cells(word);
                model.spend(cells)?;
                if i == 0 || i == word.len() + 1 {
                    cobar::check_storage(cells + 1)?;
                    let mut extended = Vec::with_capacity(word.len() + 1);
                    if i == 0 {
                        extended.push(self.unit.clone());
                    }
                    extended.extend_from_slice(word);
                    if i != 0 {
                        extended.push(self.unit.clone());
                    }
                    next.toggle(model, extended.into(), cells + 1)?;
                } else {
                    let letter = &word[i - 1];
                    for (left, right) in self.coproduct(model, letter)?.iter() {
                        let size = cells + 1 + left.len() + right.len() - letter.len();
                        cobar::check_storage(size)?;
                        let mut extended = Vec::with_capacity(word.len() + 1);
                        extended.extend_from_slice(&word[..i - 1]);
                        extended.extend([left.clone(), right.clone()]);
                        extended.extend_from_slice(&word[i..]);
                        next.toggle(model, extended.into(), size)?;
                    }
                }
            }
            words = Cow::Owned(next.into_vec());
        }
        Ok(words)
    }
}

fn multiply_words(
    model: &mut Model,
    profile: Profile,
    left: &Word,
    right: &Word,
) -> Result<Option<Word>> {
    debug_assert_eq!(left.len(), right.len());
    model.spend(1 + left.len())?;
    let mut term = Vec::with_capacity(left.len());
    let mut cells = 1;
    for (left, right) in left.iter().zip(right.iter()) {
        if left.is_empty() && right.is_empty() {
            // Projection kills any word with a unit factor.
            return Ok(None);
        }
        let Some(letter) = multiply_letters(model, profile, left, right)? else {
            return Ok(None);
        };
        cells += 1 + letter.len();
        cobar::check_storage(cells)?;
        term.push(letter);
    }
    Ok(Some(term.into()))
}

fn multiply_letters(
    model: &mut Model,
    profile: Profile,
    left: &Monomial,
    right: &Monomial,
) -> Result<Option<Monomial>> {
    let width = left.len().max(right.len());
    model.spend(1 + width)?;
    if left.is_empty() {
        return Ok(Some(right.clone()));
    }
    if right.is_empty() {
        return Ok(Some(left.clone()));
    }
    let mut product = Vec::with_capacity(width);
    for i in 0..width {
        let sum = left
            .get(i)
            .copied()
            .unwrap_or(0)
            .checked_add(right.get(i).copied().unwrap_or(0))
            .ok_or_else(degree_overflow)?;
        if !profile.allows_power(i, sum) {
            return Ok(None);
        }
        product.push(sum);
    }
    Ok(Some(product.into()))
}

pub(crate) fn square(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    k: usize,
    class: &[usize],
) -> Result<Value> {
    let (profile, s, _) = model.kind(source).ext_grade().unwrap();
    if class.is_empty() || k > s {
        return Ok(Value::Class([].into()));
    }
    let words = model.adams(source)?.representative(model, class)?;
    let square = if k == 0 {
        frobenius(model, profile, &words)?
    } else {
        cup_support(model, profile, s - k, &words, &words)?
    };
    adams::class(model, target, &square)
}

pub(crate) fn frobenius(model: &mut Model, profile: Profile, words: &[Word]) -> Result<Vec<Word>> {
    if let Some(word) = words.first() {
        let (s, t) = cobar::word_metadata(word)?;
        square_bidegree(0, s, t)?;
    }
    // Componentwise squaring is additive over F2. Doubling exponents
    // preserves the order of the surviving words, so no cross terms are needed.
    let mut squares = Vec::new();
    let mut cells = 0;
    for word in words {
        if let Some(square) = multiply_words(model, profile, word, word)? {
            cells += cobar::word_cells(&square);
            cobar::check_storage(cells)?;
            squares.push(square);
        }
    }
    Ok(squares)
}
