use std::collections::VecDeque;
use std::rc::Rc;

use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::{Model, Value};

// The twisted comparison can inspect all g! Levi permutations.
const MAX_LEVI_SEARCH_GENUS: usize = 8;

pub(crate) fn persist(model: &mut Model, sequence: Vec<usize>) -> Result<Value> {
    let mut previous = 0;
    model.spend(sequence.len())?;
    for (index, &entry) in sequence.iter().enumerate() {
        if !matches!(entry.checked_sub(previous), Some(0 | 1)) {
            let position = index + 1;
            return Err(Error::plain(format!(
                "EO final sequence must increase by zero or one at position {position}"
            )));
        }
        previous = entry;
    }
    model.allocate(1 + sequence.len())?;
    Ok(Value::EOType(Rc::from(sequence)))
}

pub(crate) fn dimension(model: &mut Model, sequence: &[usize]) -> Result<usize> {
    model.spend(sequence.len())?;
    sequence.iter().try_fold(0usize, |total, &entry| {
        total
            .checked_add(entry)
            .ok_or_else(|| Error::plain("EO stratum dimension is out of range"))
    })
}

pub(crate) fn bruhat(model: &mut Model, lower: &[usize], upper: &[usize]) -> Result<bool> {
    model.spend(lower.len())?;
    Ok(lower.iter().zip(upper).all(|(x, y)| x <= y))
}

fn weyl(model: &mut Model, sequence: &[usize]) -> Result<Vec<usize>> {
    let genus = sequence.len();
    let n = genus
        .checked_mul(2)
        .ok_or_else(|| Error::plain("EO Weyl permutation size is out of range"))?;
    model.spend(n)?;
    let mut permutation = vec![0; n];
    let mut fixed = 0;
    let mut jump = genus;
    let mut previous = 0;
    // Minimal left coset representative for W_J with J of type A_(g-1).
    // A jump places the corresponding image in the upper half of the range.
    for (i, &entry) in sequence.iter().enumerate() {
        let image = if entry == previous {
            let image = fixed;
            fixed += 1;
            image
        } else {
            let image = jump;
            jump += 1;
            image
        };
        permutation[i] = image;
        permutation[n - 1 - i] = n - 1 - image;
        previous = entry;
    }
    Ok(permutation)
}

// Return a zero-based Levi permutation u with u^{-1} w_low delta(u) <= w_high.
// For split GSp, delta is conjugation by the reversal in the Levi.
pub(crate) fn closure_witness(
    model: &mut Model,
    lower: &[usize],
    upper: &[usize],
) -> Result<Option<Vec<usize>>> {
    let g = lower.len();
    if bruhat(model, lower, upper)? {
        return Ok(Some((0..g).collect()));
    }
    if dimension(model, lower)? >= dimension(model, upper)? {
        return Ok(None);
    }
    if g > MAX_LEVI_SEARCH_GENUS {
        return Err(Error::plain(format!(
            "EO closure Levi search supports genus at most {MAX_LEVI_SEARCH_GENUS}"
        )));
    }
    let n = 2 * g;
    let lower = weyl(model, lower)?;
    let upper = weyl(model, upper)?;
    let prefixes = prefixes(&upper);
    let mut u: Vec<_> = (0..g).collect();
    let mut inverse = vec![0; n];
    let mut twist = vec![0; n];
    loop {
        model.spend(1 + g)?;
        for (i, &image) in u.iter().enumerate() {
            inverse[image] = i;
            inverse[n - 1 - image] = n - 1 - i;
            twist[i] = g - 1 - u[g - 1 - i];
            twist[n - 1 - i] = n - 1 - twist[i];
        }
        if prefix_bruhat(model, &prefixes, twist.iter().map(|&i| inverse[lower[i]]))? {
            return Ok(Some(u));
        }
        if !advance_permutation(&mut u) {
            return Ok(None);
        }
    }
}

fn prefixes(permutation: &[usize]) -> Vec<u32> {
    let mut seen = 0;
    permutation[..permutation.len() - 1]
        .iter()
        .map(|&value| {
            seen |= 1 << value;
            seen
        })
        .collect()
}

// The type-C Bruhat order is the induced order on balanced permutations
// in S_(2g). In S_(2g), compare the sorted values in every prefix.
fn prefix_bruhat(
    model: &mut Model,
    upper: &[u32],
    lower: impl Iterator<Item = usize>,
) -> Result<bool> {
    let mut seen = 0;
    for (&upper, value) in upper.iter().zip(lower) {
        seen |= 1u32 << value;
        model.spend(1 + seen.count_ones() as usize)?;
        let (mut a, mut b) = (seen, upper);
        while a != 0 {
            if a.trailing_zeros() > b.trailing_zeros() {
                return Ok(false);
            }
            a &= a - 1;
            b &= b - 1;
        }
    }
    Ok(true)
}

fn advance_permutation(values: &mut [usize]) -> bool {
    let Some(pivot) = values.windows(2).rposition(|pair| pair[0] < pair[1]) else {
        return false;
    };
    let successor = values
        .iter()
        .rposition(|&value| value > values[pivot])
        .unwrap();
    values.swap(pivot, successor);
    values[pivot + 1..].reverse();
    true
}

#[derive(Clone, Copy)]
pub(crate) enum Report {
    Dimension,
    Weyl,
    FirstSlope,
    SlopeCycles,
}

impl Report {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Dimension => "eoDim",
            Self::Weyl => "eoWeyl",
            Self::FirstSlope => "eoFirstSlope",
            Self::SlopeCycles => "eoSlopeCycles",
        }
    }

    pub(crate) fn evaluate(self, model: &mut Model, sequence: &[usize]) -> Result<String> {
        let value = match self {
            Self::Dimension => dimension(model, sequence)?.to_string(),
            Self::Weyl => {
                let mut permutation = weyl(model, sequence)?;
                for image in &mut permutation {
                    *image += 1;
                }
                format!("{permutation:?}")
            }
            Self::FirstSlope | Self::SlopeCycles => {
                let (high, count) = periodic_set(model, sequence)?;
                // Harashita: Psi on its ordered periodic set acts by rotation
                // through 'high' places, so this gcd also counts its cycles.
                let cycles = gcd(high, count);
                match self {
                    Self::FirstSlope if cycles == count => (high / cycles).to_string(),
                    Self::FirstSlope => format!("{}/{}", high / cycles, count / cycles),
                    Self::SlopeCycles => cycles.to_string(),
                    _ => unreachable!(),
                }
            }
        };
        Ok(format!("{} = {value}", self.name()))
    }
}

fn periodic_set(model: &mut Model, sequence: &[usize]) -> Result<(usize, usize)> {
    let genus = sequence.len();
    if genus == 0 {
        return Err(Error::plain(
            "EO first-slope invariants require positive genus",
        ));
    }
    let dimension = genus
        .checked_mul(2)
        .ok_or_else(|| Error::plain("EO first-slope dimension is out of range"))?;
    check_storage(dimension.saturating_mul(3))?;
    model.spend(5 * dimension)?;
    let mut successors = Vec::with_capacity(dimension);
    let mut indegrees = vec![0usize; dimension];
    let value = |i: usize| if i == 0 { 0 } else { sequence[i - 1] };
    for i in 1..=dimension {
        // Extend the elementary sequence to the full final flag by duality.
        let psi = if i <= genus {
            value(i)
        } else {
            let j = dimension - i;
            genus - j + value(j)
        };
        let image = if psi == 0 { genus + i } else { psi };
        let target = image - 1;
        indegrees[target] += 1;
        successors.push(target);
    }
    let mut queue = VecDeque::with_capacity(dimension);
    queue.extend(
        indegrees
            .iter()
            .enumerate()
            .filter_map(|(i, &n)| (n == 0).then_some(i)),
    );
    // The stable image of a self-map of a finite set is its periodic set.
    while let Some(i) = queue.pop_front() {
        let target = successors[i];
        indegrees[target] -= 1;
        if indegrees[target] == 0 {
            queue.push_back(target);
        }
    }
    let high = indegrees[genus..].iter().filter(|&&n| n != 0).count();
    let count = high + indegrees[..genus].iter().filter(|&&n| n != 0).count();
    Ok((high, count))
}

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}
