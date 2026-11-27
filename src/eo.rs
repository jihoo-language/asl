use std::collections::{BTreeSet, VecDeque};
use std::rc::Rc;

use crate::binary_linear::check_storage;
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};

// Preserve the lexicographically first Levi witness in small genus.
const DIRECT_LEVI_GENUS: usize = 8;

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
    let low_length = dimension(model, lower)?;
    let high_length = dimension(model, upper)?;
    if low_length >= high_length {
        return Ok(None);
    }
    let scratch = g.saturating_mul(32);
    check_interval_storage(scratch)?;
    let lower = weyl(model, lower)?;
    let upper = weyl(model, upper)?;
    if g > DIRECT_LEVI_GENUS {
        return interval_witness(model, &lower, upper, high_length - low_length, scratch);
    }
    direct_witness(model, &lower, &upper)
}

fn direct_witness(
    model: &mut Model,
    lower: &[usize],
    upper: &[usize],
) -> Result<Option<Vec<usize>>> {
    let n = lower.len();
    let g = n / 2;
    let prefixes = prefixes(upper);
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

fn interval_witness(
    model: &mut Model,
    lower: &[usize],
    upper: Vec<usize>,
    gap: usize,
    scratch: usize,
) -> Result<Option<Vec<usize>>> {
    // Minimality of lower in its left Levi coset gives
    // l(u^-1 lower) = l(u) + l(lower). Since delta preserves length,
    // l(u^-1 lower delta(u)) >= l(lower). The only possible targets
    // lie in the downward Bruhat interval truncated at that length.
    let source = signed_cycles(model, lower)?;
    let state_cells = upper.len() + 8;
    check_interval_storage(scratch + state_cells)?;
    let mut frontier = BTreeSet::from([upper]);
    for depth in 0..=gap {
        let mut next = BTreeSet::new();
        while let Some(mut target) = frontier.pop_first() {
            let cycles = signed_cycles(model, &target)?;
            if let Some(witness) = match_cycles(model, &source, &cycles)? {
                return Ok(Some(witness));
            }
            if depth < gap {
                let used = scratch + (1 + frontier.len()) * state_cells;
                lower_covers(model, &mut target, &mut next, used, state_cells)?;
            }
        }
        if next.is_empty() {
            return Ok(None);
        }
        frontier = next;
    }
    Ok(None)
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct SignedCycle {
    signs: Vec<bool>,
    points: Vec<usize>,
}

fn signed_cycles(model: &mut Model, permutation: &[usize]) -> Result<Vec<SignedCycle>> {
    let n = permutation.len();
    let g = n / 2;
    model.spend(g)?;
    let mut seen = vec![false; g];
    let mut cycles = Vec::new();
    for start in 0..g {
        if seen[start] {
            continue;
        }
        let mut points = Vec::new();
        let mut signs = Vec::new();
        let mut vertex = start;
        while !seen[vertex] {
            seen[vertex] = true;
            points.push(vertex);
            // This is the signed permutation w*r, for the longest
            // unsigned Levi element r. The second half represents
            // the oppositely signed members of the same g pairs.
            let image = permutation[g - 1 - vertex];
            let sign = image >= g;
            signs.push(sign);
            vertex = if sign { n - 1 - image } else { image };
        }
        let rotation = first_rotation(model, &signs)?;
        signs.rotate_left(rotation);
        points.rotate_left(rotation);
        cycles.push(SignedCycle { signs, points });
    }
    model.spend(g * (1 + cycles.len().max(1).ilog2() as usize))?;
    cycles.sort_unstable();
    Ok(cycles)
}

// Two candidate starts for the lexicographically least cyclic rotation.
// A first unequal symbol excludes all starts skipped by that candidate.
fn first_rotation(model: &mut Model, signs: &[bool]) -> Result<usize> {
    let n = signs.len();
    let (mut left, mut right, mut matched) = (0, 1, 0);
    while left < n && right < n && matched < n {
        model.spend(1)?;
        let a = signs[(left + matched) % n];
        let b = signs[(right + matched) % n];
        if a == b {
            matched += 1;
        } else {
            let (loser, winner) = if a {
                (&mut left, right)
            } else {
                (&mut right, left)
            };
            *loser += matched + 1;
            if *loser == winner {
                *loser += 1;
            }
            matched = 0;
        }
    }
    Ok(left.min(right))
}

fn match_cycles(
    model: &mut Model,
    source: &[SignedCycle],
    target: &[SignedCycle],
) -> Result<Option<Vec<usize>>> {
    // x = u^-1 lower (r*u*r) iff x*r = u^-1 (lower*r) u.
    // An unsigned conjugacy preserves the outgoing signs of every
    // directed cycle. Rotating or permuting matching cycles produces u.
    if source.len() != target.len() {
        return Ok(None);
    }
    let g: usize = source.iter().map(|cycle| cycle.points.len()).sum();
    model.spend(1 + g)?;
    if !source.iter().zip(target).all(|(a, b)| a.signs == b.signs) {
        return Ok(None);
    }
    let mut witness = vec![0; g];
    for (a, b) in source.iter().zip(target) {
        for (&x, &y) in a.points.iter().zip(&b.points) {
            witness[y] = x;
        }
    }
    Ok(Some(witness))
}

fn lower_covers(
    model: &mut Model,
    permutation: &mut [usize],
    next: &mut BTreeSet<Vec<usize>>,
    used: usize,
    state_cells: usize,
) -> Result<()> {
    let n = permutation.len();
    for i in 0..n / 2 {
        for j in i + 1..n - i {
            // The type-C reflections are mirrored transpositions of
            // positions, or a single swap of a complementary pair.
            let mirror_i = n - 1 - j;
            let mirror_j = n - 1 - i;
            let mut difference = swap_length_change(model, permutation, i, j)?;
            if i != mirror_i {
                difference += swap_length_change(model, permutation, mirror_i, mirror_j)?;
            }
            if difference == -2 {
                let search = 1 + next.len().max(1).ilog2() as usize;
                model.spend(n * search)?;
                if !next.contains(permutation) {
                    check_interval_storage(used + (next.len() + 1) * state_cells)?;
                    model.spend(n * (1 + search))?;
                    next.insert(permutation.to_vec());
                }
            }
            permutation.swap(i, j);
            if i != mirror_i {
                permutation.swap(mirror_i, mirror_j);
            }
        }
    }
    Ok(())
}

// In the balanced permutation, twice the type-C length is the
// inversion count plus the number of negative images in the first half.
// Compute their change for one swap without scanning the full matrix.
fn swap_length_change(
    model: &mut Model,
    permutation: &mut [usize],
    i: usize,
    j: usize,
) -> Result<isize> {
    model.spend(j - i)?;
    let (a, b) = (permutation[i], permutation[j]);
    let between = permutation[i + 1..j]
        .iter()
        .filter(|&&x| a.min(b) < x && x < a.max(b))
        .count();
    let gain = (1 + 2 * between) as isize;
    let inversions = if a < b { gain } else { -gain };
    let half = permutation.len() / 2;
    let signs = if i < half && j >= half {
        (b >= half) as isize - (a >= half) as isize
    } else {
        0
    };
    permutation.swap(i, j);
    Ok(inversions + signs)
}

fn check_interval_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "EO closure interval exceeds the storage limit",
        ));
    }
    Ok(())
}
