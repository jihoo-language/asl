use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::rc::Rc;

use crate::complex::xor_sorted;
use crate::error::{Error, Result};
use crate::model::{Model, Value};

const MAX_WORD_CELLS: usize = 300_000;

// A word has only positive entries. The empty word is Sq^0.
pub(crate) type Word = Rc<[usize]>;

pub(crate) fn zero() -> Value {
    Value::Operation([].into())
}

pub(crate) fn square(model: &mut Model, degree: usize) -> Result<Value> {
    model.allocate(2)?;
    let word: Word = if degree == 0 {
        [].into()
    } else {
        [degree].into()
    };
    Ok(Value::Operation([word].into()))
}

pub(crate) fn add(model: &mut Model, left: &[Word], right: &[Word]) -> Result<Value> {
    let cells = |words: &[Word]| words.iter().map(|word| 1 + word.len()).sum::<usize>();
    model.spend(cells(left) + cells(right))?;
    let words = xor_sorted(left, right);
    model.allocate(words.len())?;
    Ok(Value::Operation(words.into()))
}

pub(crate) fn compare(model: &mut Model, left: &[Word], right: &[Word]) -> Result<Ordering> {
    for (left, right) in left.iter().rev().zip(right.iter().rev()) {
        model.spend(1 + left.len().min(right.len()))?;
        let order = left.cmp(right);
        if order != Ordering::Equal {
            return Ok(order);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

pub(crate) fn multiply(model: &mut Model, left: &[Word], right: &[Word]) -> Result<Value> {
    let mut pending = BTreeSet::new();
    let mut cells = 0;
    for a in left {
        for b in right {
            let mut word = Vec::with_capacity(a.len() + b.len());
            word.extend_from_slice(a);
            word.extend_from_slice(b);
            toggle(model, &mut pending, &mut cells, word)?;
        }
    }

    let mut normal = Vec::new();
    while let Some(word) = pending.pop_first() {
        model.spend(1 + word.len())?;
        let violation = word.windows(2).position(|pair| pair[0] / 2 < pair[1]);
        if let Some(index) = violation {
            cells -= 1 + word.len();
            let (a, b) = (word[index], word[index + 1]);
            model.spend(1 + a / 2)?;
            for j in 0..=a / 2 {
                // Lucas's theorem modulo 2: each set bit of the lower
                // argument must also be set in the upper argument.
                let upper = b - j - 1;
                let lower = a - 2 * j;
                if lower & !upper == 0 {
                    let mut rewritten = Vec::with_capacity(word.len());
                    rewritten.extend_from_slice(&word[..index]);
                    rewritten.push(a + (b - j));
                    if j > 0 {
                        rewritten.push(j);
                    }
                    rewritten.extend_from_slice(&word[index + 2..]);
                    toggle(model, &mut pending, &mut cells, rewritten)?;
                }
            }
        } else {
            // Each rewrite strictly increases the word lexicographically, so
            // no later step can contribute to a word already in normal form.
            normal.push(Word::from(word));
        }
    }
    model.allocate(cells)?;
    Ok(Value::Operation(normal.into()))
}

fn toggle(
    model: &mut Model,
    pending: &mut BTreeSet<Vec<usize>>,
    cells: &mut usize,
    word: Vec<usize>,
) -> Result<()> {
    let size = 1 + word.len();
    model.spend(size * (1 + (1 + pending.len()).ilog2() as usize))?;
    if pending.remove(&word) {
        *cells -= size;
    } else {
        if *cells + size > MAX_WORD_CELLS {
            return Err(Error::plain(
                "Adem reduction exceeds the word storage limit",
            ));
        }
        *cells += size;
        pending.insert(word);
    }
    Ok(())
}

pub(crate) fn elements(model: &mut Model, degree: usize, limit: usize) -> Result<Vec<Value>> {
    let basis_limit = limit.ilog2() as usize;
    let overflow = || {
        Error::plain(format!(
            "type exceeds the enumeration limit of {limit} elements"
        ))
    };
    // The words (degree) and (degree-j, j) alone force this lower bound.
    // This also bounds the subsequent search before it starts.
    if degree > 0 && degree / 3 >= basis_limit {
        return Err(overflow());
    }

    let mut basis = Vec::new();
    if !admissible_words(
        model,
        degree,
        degree,
        &mut Vec::new(),
        &mut basis,
        basis_limit,
    )? {
        return Err(overflow());
    }
    let size = 1usize
        .checked_shl(basis.len() as u32)
        .ok_or_else(overflow)?;
    let word_cells: usize = basis.iter().map(|word| word.len() + 1).sum();
    model.allocate(word_cells + size + (size / 2) * basis.len())?;
    model.spend(size * (basis.len() + 1))?;
    Ok((0..size)
        .map(|mask| {
            let words = basis
                .iter()
                .enumerate()
                .filter(|(index, _)| (mask >> index) & 1 != 0)
                .map(|(_, word)| word.clone())
                .collect();
            Value::Operation(words)
        })
        .collect())
}

fn admissible_words(
    model: &mut Model,
    remaining: usize,
    bound: usize,
    word: &mut Vec<usize>,
    basis: &mut Vec<Word>,
    limit: usize,
) -> Result<bool> {
    model.spend(1)?;
    if remaining == 0 {
        if basis.len() == limit {
            return Ok(false);
        }
        basis.push(Word::from(word.as_slice()));
        return Ok(true);
    }
    for head in 1..=remaining.min(bound) {
        word.push(head);
        if !admissible_words(model, remaining - head, head / 2, word, basis, limit)? {
            return Ok(false);
        }
        word.pop();
    }
    Ok(true)
}
