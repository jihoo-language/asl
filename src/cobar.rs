use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap};
use std::rc::Rc;

use crate::complex::xor_sorted;
use crate::error::{Error, Result, Span};
use crate::model::{Model, Value};

const MAX_WORD_CELLS: usize = 300_000;

// Exponents are in ascending generator order, with trailing zeros removed.
// Only positive-degree monomials occur in words; the empty word is the unit.
pub(crate) type Monomial = Rc<[usize]>;
pub(crate) type Word = Rc<[Monomial]>;
type CoproductTerm = (Monomial, Monomial);

// Number rows in encounter order without enumerating an entire codomain.
#[derive(Default)]
pub(crate) struct WordRows {
    indices: HashMap<Word, usize>,
    cells: usize,
}

impl WordRows {
    pub(crate) fn coordinate(&mut self, model: &mut Model, word: Word) -> Result<usize> {
        let size = 1 + word_cells(&word);
        model.spend(size)?;
        let next = self.indices.len();
        match self.indices.entry(word) {
            std::collections::hash_map::Entry::Occupied(entry) => Ok(*entry.get()),
            std::collections::hash_map::Entry::Vacant(entry) => {
                check_storage(self.cells + size)?;
                self.cells += size;
                Ok(*entry.insert(next))
            }
        }
    }
}

pub(crate) fn zero() -> Value {
    Value::Cobar([].into())
}

pub(crate) fn literal(
    model: &mut Model,
    factors: Vec<(Vec<usize>, Span)>,
) -> Result<((usize, usize), Value)> {
    let mut word = Vec::<Monomial>::with_capacity(factors.len());
    for (mut exponents, span) in factors {
        model.spend(1 + exponents.len())?;
        let width = exponents.iter().rposition(|&exponent| exponent != 0);
        let Some(last) = width else {
            return Err(Error::new(
                span,
                "cobar word factors must have positive degree",
            ));
        };
        exponents.truncate(last + 1);
        word.push(exponents.into());
    }
    let word: Word = word.into();
    let bidegree = word_metadata(&word)?;
    let cells = word_cells(&word);
    check_storage(cells)?;
    model.allocate(cells + 1)?;
    Ok((bidegree, Value::Cobar([word].into())))
}

fn generator_degree(index: usize) -> usize {
    usize::MAX >> (usize::BITS as usize - index)
}

fn degree_overflow() -> Error {
    Error::plain("cobar internal degree is out of range")
}

pub(crate) fn word_metadata(word: &Word) -> Result<(usize, usize)> {
    let mut degree = 0usize;
    for factor in word.iter() {
        if factor.len() > usize::BITS as usize {
            return Err(degree_overflow());
        }
        for (index, &exponent) in factor.iter().enumerate() {
            let contribution = exponent
                .checked_mul(generator_degree(index + 1))
                .ok_or_else(degree_overflow)?;
            degree = degree
                .checked_add(contribution)
                .ok_or_else(degree_overflow)?;
        }
    }
    Ok((word.len(), degree))
}

pub(crate) fn word_cells(word: &Word) -> usize {
    1 + word.iter().map(|factor| 1 + factor.len()).sum::<usize>()
}

pub(crate) fn add(model: &mut Model, left: &[Word], right: &[Word]) -> Result<Value> {
    model.spend(left.iter().chain(right).map(word_cells).sum())?;
    let words = xor_sorted(left, right);
    check_storage(words.iter().map(word_cells).sum())?;
    model.allocate(words.len())?;
    Ok(Value::Cobar(words.into()))
}

pub(crate) fn compare(model: &mut Model, left: &[Word], right: &[Word]) -> Result<Ordering> {
    // Compare from the most significant word to match enumeration by bitmask.
    for (left, right) in left.iter().rev().zip(right.iter().rev()) {
        model.spend(word_cells(left))?;
        let order = left.cmp(right);
        if order != Ordering::Equal {
            return Ok(order);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

fn join(left: &Word, right: &Word) -> Word {
    if left.is_empty() {
        return right.clone();
    }
    if right.is_empty() {
        return left.clone();
    }
    left.iter().chain(right.iter()).cloned().collect()
}

pub(crate) fn concat_words(model: &mut Model, left: &Word, right: &Word) -> Result<Word> {
    let cells = word_cells(left) + word_cells(right) - 1;
    check_storage(cells)?;
    model.spend(cells)?;
    let (_, left_degree) = word_metadata(left)?;
    let (_, right_degree) = word_metadata(right)?;
    left_degree
        .checked_add(right_degree)
        .ok_or_else(degree_overflow)?;
    model.allocate(cells)?;
    Ok(join(left, right))
}

pub(crate) fn concat(model: &mut Model, left: &[Word], right: &[Word]) -> Result<Value> {
    Ok(Value::Cobar(concat_support(model, left, right)?.into()))
}

pub(crate) fn concat_support(
    model: &mut Model,
    left: &[Word],
    right: &[Word],
) -> Result<Vec<Word>> {
    if right.is_empty() {
        return Ok(Vec::new());
    }
    let mut words = Vec::new();
    let mut cells = 0;
    // Each term has a unique split at the common left word length, so
    // the nested iteration gives distinct words in lexicographic order.
    for a in left {
        for b in right {
            cells += word_cells(a) + word_cells(b) - 1;
            check_storage(cells)?;
            words.push(concat_words(model, a, b)?);
        }
    }
    Ok(words)
}

pub(crate) fn differential_word(model: &mut Model, word: &Word) -> Result<Vec<Word>> {
    differential_support(model, std::slice::from_ref(word))
}

pub(crate) fn differential(model: &mut Model, words: &[Word]) -> Result<Value> {
    let support = match words {
        [word] => differential_word(model, word)?,
        _ => differential_support(model, words)?,
    };
    Ok(Value::Cobar(support.into()))
}

pub(crate) fn differential_support(model: &mut Model, words: &[Word]) -> Result<Vec<Word>> {
    let mut support = Support::new();
    for word in words {
        differentiate_into(model, word, &mut support)?;
    }
    support.finish(model)
}

fn differentiate_into(model: &mut Model, word: &Word, support: &mut Support<Word>) -> Result<()> {
    let cells = word_cells(word);
    model.spend(cells)?;
    for (position, factor) in word.iter().enumerate() {
        for (left, right) in reduced_coproduct(model, factor)? {
            let new_cells = cells + 1 + left.len() + right.len() - factor.len();
            check_storage(new_cells)?;
            model.spend(new_cells)?;
            let mut result = Vec::with_capacity(word.len() + 1);
            result.extend_from_slice(&word[..position]);
            result.extend([left, right]);
            result.extend_from_slice(&word[position + 1..]);
            support.toggle(model, result.into(), new_cells)?;
        }
    }
    Ok(())
}

fn reduced_coproduct(model: &mut Model, factor: &Monomial) -> Result<Vec<CoproductTerm>> {
    let unit: Monomial = [].into();
    let mut terms = Support::new();
    terms.toggle(model, (unit.clone(), unit), 1)?;
    for (index, &exponent) in factor.iter().enumerate() {
        let generator = index + 1;
        let mut remaining = exponent;
        // In characteristic two, each set bit contributes one Frobenius
        // power of the generator coproduct.
        while remaining != 0 {
            let power = 1usize << remaining.trailing_zeros();
            remaining ^= power;
            let mut next = Support::new();
            for (left, right) in terms.terms {
                for i in 0..=generator {
                    // Delta(xi_n) = sum_i xi_(n-i)^(2^i) tensor xi_i.
                    // The i = n term has a unit on the left.
                    let left_index = generator - i;
                    let left_power = if left_index == 0 {
                        0
                    } else {
                        power.checked_mul(1usize << i).ok_or_else(degree_overflow)?
                    };
                    model.spend(1 + left.len().max(left_index) + right.len().max(i))?;
                    let a = with_power(&left, left_index, left_power)?;
                    let b = with_power(&right, i, power)?;
                    let cells = 1 + a.len() + b.len();
                    next.toggle(model, (a, b), cells)?;
                }
            }
            terms = next;
        }
    }
    Ok(terms
        .terms
        .into_iter()
        .filter(|(left, right)| !left.is_empty() && !right.is_empty())
        .collect())
}

fn with_power(monomial: &Monomial, index: usize, power: usize) -> Result<Monomial> {
    if index == 0 {
        return Ok(monomial.clone());
    }
    let mut exponents = monomial.to_vec();
    exponents.resize(exponents.len().max(index), 0);
    exponents[index - 1] = exponents[index - 1]
        .checked_add(power)
        .ok_or_else(degree_overflow)?;
    Ok(exponents.into())
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_WORD_CELLS {
        return Err(Error::plain(
            "cobar calculation exceeds the word storage limit",
        ));
    }
    Ok(())
}

struct Support<T> {
    terms: BTreeSet<T>,
    cells: usize,
}

impl<T: Ord> Support<T> {
    fn new() -> Self {
        Self {
            terms: BTreeSet::new(),
            cells: 0,
        }
    }

    fn toggle(&mut self, model: &mut Model, term: T, cells: usize) -> Result<()> {
        model.spend(cells * (1 + (1 + self.terms.len()).ilog2() as usize))?;
        if self.terms.remove(&term) {
            self.cells -= cells;
        } else {
            check_storage(self.cells + cells)?;
            self.cells += cells;
            self.terms.insert(term);
        }
        Ok(())
    }

    fn finish(self, model: &mut Model) -> Result<Vec<T>> {
        model.allocate(self.cells)?;
        Ok(self.terms.into_iter().collect())
    }
}

pub(crate) fn basis(model: &mut Model, s: usize, t: usize, limit: usize) -> Result<Vec<Word>> {
    model.spend(1)?;
    if s == 0 {
        return if t == 0 {
            if limit == 0 {
                return Err(basis_limit(limit));
            }
            model.allocate(1)?;
            Ok(vec![[].into()])
        } else {
            Ok(Vec::new())
        };
    }
    let Some(excess) = t.checked_sub(s) else {
        return Ok(Vec::new());
    };
    // Each place for the excess gives a word. Putting it all in one factor
    // also gives at least floor((excess + 1) / 3) + 1 powers of xi_1 and xi_2.
    if (excess != 0 && s > limit) || (excess + 1) / 3 >= limit {
        return Err(basis_limit(limit));
    }
    check_storage(s.saturating_mul(2).saturating_add(1))?;

    let mut degrees = vec![1; s];
    degrees[s - 1] += excess;
    let mut catalogue = HashMap::new();
    let mut catalogue_cells = 0;
    let mut words = Vec::new();
    let mut cells = 0;
    loop {
        model.spend(s)?;
        for &degree in &degrees {
            if let std::collections::hash_map::Entry::Vacant(entry) = catalogue.entry(degree) {
                entry.insert(monomials(model, degree, limit, &mut catalogue_cells)?);
            }
        }
        let factors: Vec<_> = degrees
            .iter()
            .map(|degree| catalogue[degree].as_slice())
            .collect();
        let mut digits = vec![0; s];
        loop {
            if words.len() == limit {
                return Err(basis_limit(limit));
            }
            let word: Word = factors
                .iter()
                .zip(&digits)
                .map(|(choices, &digit)| choices[digit].clone())
                .collect();
            let size = word_cells(&word);
            cells += size;
            check_storage(cells + catalogue_cells)?;
            model.spend(size + s)?;
            words.push(word);
            if !next_choices(&factors, &mut digits) {
                break;
            }
        }
        if !next_degrees(&mut degrees) {
            break;
        }
    }
    model.spend(cells * (1 + words.len().ilog2() as usize))?;
    model.allocate(cells)?;
    words.sort_unstable();
    Ok(words)
}

fn basis_limit(limit: usize) -> Error {
    Error::plain(format!(
        "cobar basis exceeds the enumeration limit of {limit} words"
    ))
}

fn next_degrees(degrees: &mut [usize]) -> bool {
    let Some(index) = (1..degrees.len()).rev().find(|&index| degrees[index] > 1) else {
        return false;
    };
    let excess = degrees[index] - 1;
    degrees[index - 1] += 1;
    degrees[index] = 1;
    *degrees.last_mut().unwrap() = excess;
    true
}

fn next_choices(factors: &[&[Monomial]], digits: &mut [usize]) -> bool {
    for (choices, digit) in factors.iter().zip(digits).rev() {
        *digit += 1;
        if *digit < choices.len() {
            return true;
        }
        *digit = 0;
    }
    false
}

fn monomials(
    model: &mut Model,
    degree: usize,
    limit: usize,
    cells: &mut usize,
) -> Result<Vec<Monomial>> {
    let width = degree.checked_add(1).map_or(usize::BITS, usize::ilog2) as usize;
    let mut exponents = vec![0; width];
    exponents[0] = degree;
    let mut monomials = Vec::new();
    loop {
        if monomials.len() == limit {
            return Err(basis_limit(limit));
        }
        model.spend(width)?;
        let last = exponents
            .iter()
            .rposition(|&exponent| exponent != 0)
            .unwrap();
        *cells += last + 2;
        check_storage(*cells)?;
        monomials.push(Monomial::from(&exponents[..=last]));
        // The exponent of xi_1 pays for each increment and is restored when
        // a digit carries. It is the unique remainder for the other digits.
        let mut index = 1;
        while index < width {
            let weight = generator_degree(index + 1);
            if exponents[0] >= weight {
                exponents[0] -= weight;
                exponents[index] += 1;
                break;
            }
            exponents[0] += exponents[index] * weight;
            exponents[index] = 0;
            index += 1;
        }
        if index == width {
            return Ok(monomials);
        }
    }
}
