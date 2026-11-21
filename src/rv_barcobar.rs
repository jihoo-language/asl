use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::rc::Rc;

use crate::complex::{xor_sorted, MAX_PIVOT_CELLS};
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_bar::{self, Grade};

pub(crate) type Word = Rc<[rv_bar::Word]>;
pub(crate) type Support = Rc<[Word]>;

pub(crate) fn suspend_down(grade: Grade, operation: &str) -> Result<Grade> {
    let b = grade.b.checked_sub(1).ok_or_else(|| {
        Error::plain(format!(
            "{operation} requires positive bar-cobar input degree"
        ))
    })?;
    Ok(Grade { b, ..grade })
}

pub(crate) fn homotopy_grade(grade: Grade) -> Result<Grade> {
    let b = grade.b.checked_add(1).ok_or_else(degree_limit)?;
    Ok(Grade { b, ..grade })
}

pub(crate) fn product_grade(left: Grade, right: Grade) -> Result<Grade> {
    if left.path != right.path || left.end != right.start {
        return Err(Error::plain(
            "rvBCMul requires intervals sharing an endpoint on the same bar path",
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
    Error::plain("bar-cobar degree is out of range")
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "bar-cobar calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

pub(crate) fn word_cells(word: &[rv_bar::Word]) -> usize {
    1 + word
        .iter()
        .map(|block| rv_bar::word_cells(block))
        .sum::<usize>()
}

fn support_cells(words: &[Word]) -> usize {
    1 + words.iter().map(|word| word_cells(word)).sum::<usize>()
}

fn persist(model: &mut Model, words: Vec<Word>) -> Result<Support> {
    model.allocate(support_cells(&words))?;
    Ok(words.into())
}

pub(crate) fn compare(model: &mut Model, a: &[Word], b: &[Word]) -> Result<Ordering> {
    model.spend(support_cells(a) + support_cells(b))?;
    Ok(a.iter().rev().cmp(b.iter().rev()))
}

pub(crate) fn add(model: &mut Model, a: &[Word], b: &[Word]) -> Result<Support> {
    model.spend(support_cells(a) + support_cells(b))?;
    persist(model, xor_sorted(a, b))
}

pub(crate) fn desuspend(model: &mut Model, chain: &[rv_bar::Word]) -> Result<Support> {
    let mut support = Vec::with_capacity(chain.len());
    for word in chain {
        model.spend(1 + rv_bar::word_cells(word))?;
        support.push([word.clone()].into());
    }
    persist(model, support)
}

pub(crate) fn product(model: &mut Model, left: &[Word], right: &[Word]) -> Result<Support> {
    let count = left.len().saturating_mul(right.len());
    check_storage(count)?;
    let mut support = Vec::with_capacity(count);
    let mut cells = 0;
    for a in left {
        for b in right {
            let cost = word_cells(a) + word_cells(b);
            model.spend(cost)?;
            cells += cost;
            check_storage(cells)?;
            let mut term = a.to_vec();
            term.extend_from_slice(b);
            support.push(term.into());
        }
    }
    // The prescribed join at the endpoint gives a unique decomposition.
    persist(model, support)
}

#[derive(Default)]
struct Sum {
    terms: BTreeSet<Word>,
    cells: usize,
}

impl Sum {
    fn toggle(&mut self, model: &mut Model, word: Vec<rv_bar::Word>) -> Result<()> {
        let cost = word_cells(&word);
        model.spend(cost * (1 + self.terms.len().max(1).ilog2() as usize))?;
        let word: Word = word.into();
        if self.terms.remove(&word) {
            self.cells -= cost;
        } else {
            self.cells += cost;
            check_storage(self.cells)?;
            self.terms.insert(word);
        }
        Ok(())
    }

    fn finish(self) -> Vec<Word> {
        self.terms.into_iter().collect()
    }
}

fn differential_terms(
    model: &mut Model,
    grade: Grade,
    words: &[Word],
    include_splits: bool,
) -> Result<Vec<Word>> {
    let mut result = Sum::default();
    for word in words {
        let mut left = grade.start;
        for (index, block) in word.iter().enumerate() {
            model.spend(1 + block.len())?;
            let end = block.last().unwrap().end;
            let internal = block
                .iter()
                .try_fold(0usize, |sum, factor| {
                    sum.checked_add(factor.s).and_then(|b| b.checked_add(1))
                })
                .ok_or_else(degree_limit)?;
            let interval = Grade {
                start: left,
                end,
                b: internal,
                ..grade
            };
            let cochains = rv_bar::differential(model, interval, std::slice::from_ref(block))?;
            for replacement in cochains {
                let mut term = word.to_vec();
                term[index] = replacement;
                result.toggle(model, term)?;
            }
            if include_splits {
                for gap in 1..block.len() {
                    let mut term = Vec::with_capacity(word.len() + 1);
                    term.extend_from_slice(&word[..index]);
                    term.push(block[..gap].into());
                    term.push(block[gap..].into());
                    term.extend_from_slice(&word[index + 1..]);
                    result.toggle(model, term)?;
                }
            }
            left = end;
        }
    }
    Ok(result.finish())
}

pub(crate) fn differential(model: &mut Model, grade: Grade, words: &[Word]) -> Result<Support> {
    let terms = differential_terms(model, grade, words, true)?;
    persist(model, terms)
}

fn merge_first_gap(word: &[rv_bar::Word]) -> Option<Vec<rv_bar::Word>> {
    if word.len() < 2 || word[0].len() != 1 {
        return None;
    }
    let mut block = Vec::with_capacity(word[1].len() + 1);
    block.push(word[0][0]);
    block.extend_from_slice(&word[1]);
    let mut merged = Vec::with_capacity(word.len() - 1);
    merged.push(block.into());
    merged.extend_from_slice(&word[2..]);
    Some(merged)
}

fn singleton(word: &[rv_bar::Word]) -> Option<usize> {
    let [block] = word else {
        return None;
    };
    let [factor] = block.as_ref() else {
        return None;
    };
    Some(factor.generator)
}

pub(crate) enum Contraction {
    Projection,
    Homotopy,
}

pub(crate) fn contract(
    model: &mut Model,
    grade: Grade,
    terms: &[Word],
    operation: Contraction,
) -> Result<crate::model::Value> {
    use crate::model::Value;
    let mut current = terms.to_vec();
    let mut classes = BTreeSet::new();
    let mut homotopy = Sum::default();
    while !current.is_empty() {
        let mut elementary = Sum::default();
        for word in &current {
            model.spend(word_cells(word))?;
            if let (Contraction::Projection, Some(generator)) = (&operation, singleton(word)) {
                model.spend(1 + classes.len().max(1).ilog2() as usize)?;
                if !classes.remove(&generator) {
                    classes.insert(generator);
                    check_storage(classes.len())?;
                }
            }
            if let Some(merged) = merge_first_gap(word) {
                if matches!(operation, Contraction::Homotopy) {
                    homotopy.toggle(model, merged.clone())?;
                }
                elementary.toggle(model, merged)?;
            }
        }
        // The bar part lowers flattened length, so the perturbation series
        // terminates on every word and uses only supported terms.
        current = differential_terms(model, grade, &elementary.finish(), false)?;
    }
    match operation {
        Contraction::Projection => {
            model.allocate(classes.len())?;
            Ok(Value::Class(classes.into_iter().collect()))
        }
        Contraction::Homotopy => Ok(Value::RVBarCobar(persist(model, homotopy.finish())?)),
    }
}

pub(crate) fn basis(model: &mut Model, grade: Grade, limit: usize) -> Result<Vec<Word>> {
    let span = grade.end - grade.start;
    if grade.b > span.saturating_mul(usize::BITS as usize + 1) {
        return Ok(Vec::new());
    }
    let mut plan = Vec::new();
    let mut size = 0;
    for outer in 1..=span {
        let b = grade.b.checked_add(outer).ok_or_else(degree_limit)?;
        let bar_grade = Grade { b, ..grade };
        for length in outer..=bar_grade.max_length() {
            model.spend(1)?;
            let count = rv_bar::dimension(model, bar_grade, Some(length))?;
            if count == 0 {
                continue;
            }
            if count > (limit - size) as u128 {
                return Err(basis_limit());
            }
            let count = count as usize;
            let partitions = binomial_bounded(length - 1, outer - 1, (limit - size) / count)?;
            size += count * partitions;
            plan.push((bar_grade, length, outer, count));
        }
    }
    let mut terms = Vec::with_capacity(size);
    let mut cells = 0;
    for (bar_grade, length, outer, count) in plan {
        for word in rv_bar::basis_at_length(model, bar_grade, length, count)? {
            partitions(model, &word, outer, &mut terms, &mut cells)?;
        }
    }
    model.spend(cells * (1 + terms.len().max(1).ilog2() as usize))?;
    terms.sort_unstable();
    Ok(terms)
}

fn partitions(
    model: &mut Model,
    word: &[rv_bar::Factor],
    outer: usize,
    output: &mut Vec<Word>,
    cells: &mut usize,
) -> Result<()> {
    let mut cuts: Vec<_> = (1..outer).collect();
    loop {
        model.spend(1 + 3 * word.len() + outer)?;
        *cells += 1 + 3 * word.len() + outer;
        check_storage(*cells)?;
        let mut start = 0;
        let mut blocks = Vec::with_capacity(outer);
        for end in cuts.iter().copied().chain(std::iter::once(word.len())) {
            blocks.push(word[start..end].into());
            start = end;
        }
        output.push(blocks.into());

        let Some(position) = (0..cuts.len())
            .rev()
            .find(|&i| cuts[i] < word.len() - outer + i + 1)
        else {
            return Ok(());
        };
        cuts[position] += 1;
        for i in position + 1..cuts.len() {
            cuts[i] = cuts[i - 1] + 1;
        }
    }
}

fn basis_limit() -> Error {
    Error::plain("bar-cobar basis exceeds the enumeration limit")
}

fn binomial_bounded(n: usize, k: usize, limit: usize) -> Result<usize> {
    let k = k.min(n - k);
    let mut count = 1u128;
    for j in 1..=k {
        count = count * (n - k + j) as u128 / j as u128;
        if count > limit as u128 {
            return Err(basis_limit());
        }
    }
    if count > limit as u128 {
        return Err(basis_limit());
    }
    Ok(count as usize)
}
