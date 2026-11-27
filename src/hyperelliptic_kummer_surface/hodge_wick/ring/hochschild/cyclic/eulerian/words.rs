use std::collections::btree_map::Entry;

use super::{storage, work_limit, BTreeMap, Model, Rational, Result, LIMBS};

pub(super) type Words = BTreeMap<Vec<usize>, Rational>;

pub(super) fn adams(model: &mut Model, word: &[usize], m: i128, transpose: bool) -> Result<Words> {
    // For a constant word the signed multinomial of each shuffle vanishes
    // when two blocks are odd. Summing the surviving block sizes gives
    // m^ceil(arity/2).
    if m == 0 || m == 1 || uniform(word) {
        let value = power(model, m, word.len().div_ceil(2))?;
        let mut result = Words::new();
        add(model, &mut result, word.to_vec(), value)?;
        return Ok(result);
    }
    let mut result = Words::new();
    visit(model, word.len(), m, |model, permutation, value| {
        let shuffled = if transpose {
            let mut inverse = vec![0; permutation.len()];
            for (output, &input) in permutation.iter().enumerate() {
                inverse[input] = word[output];
            }
            inverse
        } else {
            permutation.iter().map(|&i| word[i]).collect()
        };
        add(model, &mut result, shuffled, value)
    })?;
    Ok(result)
}

pub(super) fn project(
    model: &mut Model,
    word: &[usize],
    index: i128,
    transpose: bool,
) -> Result<Words> {
    let arity = word.len();
    if !supports(arity, index) {
        return Ok(Words::new());
    }
    if uniform(word) {
        return Ok(if index == arity.div_ceil(2) as i128 {
            Words::from([(word.to_vec(), Rational::integer(1))])
        } else {
            Words::new()
        });
    }
    let index = index as usize;
    let eigenvalue = power(model, 2, index)?;
    let mut output = Words::from([(word.to_vec(), Rational::integer(1))]);
    for other in (1..=arity).filter(|&j| j != index) {
        if output.is_empty() {
            break;
        }
        let other = power(model, 2, other)?;
        let difference = eigenvalue.add(model, other.negate())?.inverse()?;
        let mut filtered = Words::new();
        for (word, value) in output {
            let coefficient = value.multiply(model, difference)?;
            for (word, entry) in adams(model, &word, 2, transpose)? {
                let value = coefficient.multiply(model, entry)?;
                add(model, &mut filtered, word, value)?;
            }
            let correction = coefficient.multiply(model, other)?.negate();
            add(model, &mut filtered, word, correction)?;
        }
        output = filtered;
    }
    Ok(output)
}

fn uniform(word: &[usize]) -> bool {
    word.windows(2).all(|pair| pair[0] == pair[1])
}

pub(super) fn supports(arity: usize, index: i128) -> bool {
    if arity == 0 {
        index == 0
    } else {
        (1..=arity as i128).contains(&index)
    }
}

fn power(model: &mut Model, base: i128, mut exponent: usize) -> Result<Rational> {
    let mut value = Rational::integer(1);
    if base == 0 && exponent > 0 {
        return Ok(Rational::zero());
    }
    let mut base = Rational::integer(base);
    while exponent > 0 {
        if !exponent.is_multiple_of(2) {
            value = value.multiply(model, base)?;
        }
        exponent /= 2;
        if exponent > 0 {
            base = base.multiply(model, base)?;
        }
    }
    Ok(value)
}

fn add(model: &mut Model, words: &mut Words, word: Vec<usize>, value: Rational) -> Result<()> {
    if value.is_zero() {
        return Ok(());
    }
    let cells = (words.len() + 1).saturating_mul(word.len() + 2 * LIMBS + 8);
    match words.entry(word) {
        Entry::Vacant(entry) => {
            storage(cells)?;
            entry.insert(value);
        }
        Entry::Occupied(mut entry) => {
            let sum = entry.get().add(model, value)?;
            if sum.is_zero() {
                entry.remove();
            } else {
                entry.insert(sum);
            }
        }
    }
    Ok(())
}

fn visit(
    model: &mut Model,
    arity: usize,
    m: i128,
    mut emit: impl FnMut(&mut Model, &[usize], Rational) -> Result<()>,
) -> Result<()> {
    if m == 2 {
        // A signed binary shuffle of each pair of consecutive input blocks.
        // There are 2^arity schedules, including arity+1 copies of the identity.
        let size = u32::try_from(arity)
            .ok()
            .and_then(|exponent| 1usize.checked_shl(exponent))
            .ok_or_else(work_limit)?;
        model.spend(size.saturating_mul(arity + 1))?;
        let mut permutation = vec![0; arity];
        for subset in 0..size {
            let split = subset.count_ones() as usize;
            let (mut left, mut right) = (0, split);
            let mut negative = false;
            for (position, entry) in permutation.iter_mut().enumerate() {
                if subset & (1usize << position) == 0 {
                    *entry = right;
                    right += 1;
                } else {
                    *entry = left;
                    left += 1;
                    negative ^= !(right - split).is_multiple_of(2);
                }
            }
            emit(
                model,
                &permutation,
                Rational::integer(if negative { -1 } else { 1 }),
            )?;
        }
        return Ok(());
    }

    let count = (1..=arity).try_fold(1usize, |count, i| {
        count.checked_mul(i).ok_or_else(work_limit)
    })?;
    model.spend(count.saturating_mul(arity.saturating_mul(arity)))?;
    let coefficients = (0..arity)
        .map(|descent| binomial(model, m, arity, descent))
        .collect::<Result<Vec<_>>>()?;
    let mut inverse = vec![0; arity];
    let mut permutation = (0..arity).collect::<Vec<_>>();
    loop {
        for (i, &value) in permutation.iter().enumerate() {
            inverse[value] = i;
        }
        let descents = inverse.windows(2).filter(|pair| pair[0] > pair[1]).count();
        let value = coefficients[descents];
        if !value.is_zero() {
            let mut odd = false;
            for (i, &value) in permutation.iter().enumerate() {
                for &right in &permutation[i + 1..] {
                    odd ^= value > right;
                }
            }
            emit(
                model,
                &permutation,
                if odd { value.negate() } else { value },
            )?;
        }
        let Some(pivot) = permutation.windows(2).rposition(|pair| pair[0] < pair[1]) else {
            break;
        };
        let swap = permutation[pivot + 1..]
            .iter()
            .rposition(|&i| i > permutation[pivot])
            .unwrap()
            + pivot
            + 1;
        permutation.swap(pivot, swap);
        permutation[pivot + 1..].reverse();
    }
    Ok(())
}

fn binomial(model: &mut Model, m: i128, arity: usize, descents: usize) -> Result<Rational> {
    let first = m - descents as i128;
    if first <= 0 {
        return Ok(Rational::zero());
    }
    let first = Rational::integer(first);
    let mut result = Rational::integer(1);
    for j in 0..arity {
        let factor = first.add(model, Rational::integer(j as i128))?;
        let ratio = factor.multiply(model, Rational::integer(j as i128 + 1).inverse()?)?;
        result = result.multiply(model, ratio)?;
    }
    Ok(result)
}
