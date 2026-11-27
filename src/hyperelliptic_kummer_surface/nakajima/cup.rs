use crate::error::{Error, Result};
use crate::model::Model;

use super::boundary::{diagonal_of_unit, surface_product, ColorPair};
use super::{
    coordinate_range, insert_factor, Accumulator, Gram, Integer, Mode, State, LIMBS, POINT,
};

impl State {
    pub(crate) fn cup(&self, model: &mut Model, rhs: &Self, length: i128) -> Result<Self> {
        if self.terms.is_empty() || rhs.terms.is_empty() {
            return Ok(Self::zero());
        }
        if let Some(unit) = self.unit_coefficient(length) {
            return rhs.multiply_unit(model, length, unit);
        }
        if let Some(unit) = rhs.unit_coefficient(length) {
            return self.multiply_unit(model, length, unit);
        }
        let size = usize::try_from(length).map_err(|_| conjugator_limit())?;
        let factorial = (2..=size)
            .try_fold(1usize, usize::checked_mul)
            .ok_or_else(conjugator_limit)?;
        let work = factorial
            .saturating_mul(size.saturating_add(1))
            .saturating_mul(self.terms.len())
            .saturating_mul(rhs.terms.len());
        model.spend(work).map_err(|_| conjugator_limit())?;

        let gram = model.nakajima_gram()?;
        let inverse = model.nakajima_inverse_gram()?;
        let middle_diagonal: Vec<_> = diagonal_of_unit(model, &inverse)?
            .into_iter()
            .filter(|&(a, b, _)| a != 0 && a != POINT && b != 0 && b != POINT)
            .collect();
        let mut output = Accumulator::default();
        for (left, &a) in &self.terms {
            let left = Permutation::representative(left, size);
            for (right, &b) in &rhs.terms {
                let right = Permutation::representative(right, size);
                let product = Product {
                    gram: &gram,
                    middle_diagonal: &middle_diagonal,
                    left: &left,
                    right: &right,
                };
                let coefficient = Integer::product(a, b);
                let mut conjugator = (0..size).collect::<Vec<_>>();
                loop {
                    // Regroup the two full conjugator sums by their relative
                    // conjugator. The outer sum is the output's Fock basis vector.
                    product.contribute(model, &mut output, &conjugator, coefficient)?;
                    if !next_permutation(&mut conjugator) {
                        break;
                    }
                }
            }
        }
        output.finish()
    }

    fn unit_coefficient(&self, length: i128) -> Option<i128> {
        if self.terms.len() != 1 {
            return None;
        }
        let (key, &coefficient) = self.terms.first_key_value()?;
        match key.as_slice() {
            [] if length == 0 => Some(coefficient),
            [super::Factor {
                mode:
                    Mode {
                        length: 1,
                        color: 0,
                    },
                exponent,
            }] if *exponent == length as u128 => Some(coefficient),
            _ => None,
        }
    }

    fn multiply_unit(&self, model: &mut Model, length: i128, unit: i128) -> Result<Self> {
        // The full conjugator sum on the identity is n! times the unit.
        // Stop as soon as even the factorial itself exceeds the wide range.
        let scalar = (2..=length).try_fold(Integer::from_i128(unit), |value, k| {
            model.spend(LIMBS)?;
            scale(value, k)
        })?;
        let mut terms = std::collections::BTreeMap::new();
        for (key, &coefficient) in &self.terms {
            model.spend(1 + key.len() + LIMBS)?;
            let value = scale(scalar, coefficient)?
                .signed()
                .ok_or_else(coordinate_range)?;
            terms.insert(key.clone(), value);
        }
        Ok(Self { terms })
    }
}

struct Cycle {
    head: usize,
    color: u8,
}

struct Permutation {
    image: Vec<usize>,
    cycles: Vec<Cycle>,
}

impl Permutation {
    fn representative(monomial: &[super::Factor], size: usize) -> Self {
        let mut image = Vec::with_capacity(size);
        let mut cycles = Vec::new();
        for factor in monomial {
            for _ in 0..factor.exponent {
                let head = image.len();
                let end = head + factor.mode.length as usize;
                image.extend(head + 1..end);
                image.push(head);
                cycles.push(Cycle {
                    head,
                    color: factor.mode.color,
                });
            }
        }
        Self { image, cycles }
    }
}

struct Orbit {
    size: usize,
    left_cycles: usize,
    right_cycles: usize,
    output_cycles: Vec<usize>,
    color: u8,
    weight: i128,
}

impl Orbit {
    fn new(size: usize) -> Self {
        Self {
            size,
            left_cycles: 0,
            right_cycles: 0,
            output_cycles: Vec::new(),
            color: 0,
            weight: 1,
        }
    }

    fn multiply_color(&mut self, gram: &Gram, color: u8) -> bool {
        let (color, weight) = surface_product(gram, self.color, color);
        self.color = color;
        // A nonzero K3 color product uses at most one Gram entry.
        self.weight *= weight;
        weight != 0
    }
}

struct Product<'a> {
    gram: &'a Gram,
    middle_diagonal: &'a [ColorPair],
    left: &'a Permutation,
    right: &'a Permutation,
}

impl Product<'_> {
    fn contribute(
        &self,
        model: &mut Model,
        output: &mut Accumulator,
        conjugator: &[usize],
        mut weight: Integer,
    ) -> Result<()> {
        let size = conjugator.len();
        let left = &self.left.image;
        let mut right = vec![0; size];
        for (old, &new) in conjugator.iter().enumerate() {
            right[new] = conjugator[self.right.image[old]];
        }
        let (mut orbits, membership) = connected_orbits(left, &right);
        for cycle in &self.left.cycles {
            let orbit = &mut orbits[membership[cycle.head]];
            orbit.left_cycles += 1;
            if !orbit.multiply_color(self.gram, cycle.color) {
                return Ok(());
            }
        }
        for cycle in &self.right.cycles {
            let orbit = &mut orbits[membership[conjugator[cycle.head]]];
            orbit.right_cycles += 1;
            if !orbit.multiply_color(self.gram, cycle.color) {
                return Ok(());
            }
        }
        let lengths = product_cycles(left, &right, &membership, &mut orbits);
        let mut expansion = 1usize;
        for orbit in &mut orbits {
            let k = orbit.output_cycles.len();
            let defect = orbit.size + 2 - orbit.left_cycles - orbit.right_cycles - k;
            debug_assert!(defect.is_multiple_of(2));
            let excess = defect / 2;
            if excess > 0 {
                // The Euler element is 24[pt], and its square is zero.
                if excess > 1 || !orbit.multiply_color(self.gram, POINT) {
                    return Ok(());
                }
                orbit.weight *= 24;
            }
            weight = scale(weight, orbit.weight)?;
            let choices = match orbit.color {
                0 => k + k * (k - 1) / 2 * self.middle_diagonal.len(),
                POINT => 1,
                _ => k,
            };
            expansion = expansion.saturating_mul(choices);
        }
        let sign_defect = size + lengths.len() - self.left.cycles.len() - self.right.cycles.len();
        debug_assert!(sign_defect.is_multiple_of(2));
        if sign_defect % 4 == 2 {
            weight = weight.negate();
        }
        model
            .spend(expansion.saturating_mul(1 + lengths.len().pow(2) + LIMBS))
            .map_err(|_| {
                Error::plain("Nakajima cup product exceeds the diagonal expansion work limit")
            })?;
        let push = PushForward {
            orbits: &orbits,
            lengths: &lengths,
            middle_diagonal: self.middle_diagonal,
        };
        push.expand(output, 0, &mut vec![POINT; lengths.len()], weight)
    }
}

fn connected_orbits(left: &[usize], right: &[usize]) -> (Vec<Orbit>, Vec<usize>) {
    let size = left.len();
    let mut membership = vec![usize::MAX; size];
    let mut orbits = Vec::new();
    let mut queue = Vec::with_capacity(size);
    for head in 0..size {
        if membership[head] != usize::MAX {
            continue;
        }
        let index = orbits.len();
        queue.clear();
        queue.push(head);
        membership[head] = index;
        let mut next = 0;
        while next < queue.len() {
            let vertex = queue[next];
            for target in [left[vertex], right[vertex]] {
                if membership[target] == usize::MAX {
                    membership[target] = index;
                    queue.push(target);
                }
            }
            next += 1;
        }
        orbits.push(Orbit::new(queue.len()));
    }
    (orbits, membership)
}

fn product_cycles(
    left: &[usize],
    right: &[usize],
    membership: &[usize],
    orbits: &mut [Orbit],
) -> Vec<usize> {
    let mut seen = vec![false; left.len()];
    let mut lengths = Vec::new();
    for head in 0..left.len() {
        if seen[head] {
            continue;
        }
        let mut vertex = head;
        let mut length = 0;
        while !seen[vertex] {
            seen[vertex] = true;
            length += 1;
            vertex = left[right[vertex]];
        }
        orbits[membership[head]].output_cycles.push(lengths.len());
        lengths.push(length);
    }
    lengths
}

struct PushForward<'a> {
    orbits: &'a [Orbit],
    lengths: &'a [usize],
    middle_diagonal: &'a [ColorPair],
}

impl PushForward<'_> {
    fn expand(
        &self,
        output: &mut Accumulator,
        index: usize,
        colors: &mut [u8],
        weight: Integer,
    ) -> Result<()> {
        let Some(orbit) = self.orbits.get(index) else {
            let mut key = Vec::new();
            for (&length, &color) in self.lengths.iter().zip(colors.iter()) {
                insert_factor(
                    &mut key,
                    Mode {
                        length: length as i128,
                        color,
                    },
                )?;
            }
            return output.add(key, weight);
        };
        if orbit.color == POINT {
            return self.expand(output, index + 1, colors, weight);
        }
        for &position in &orbit.output_cycles {
            colors[position] = orbit.color;
            self.expand(output, index + 1, colors, weight)?;
            colors[position] = POINT;
        }
        if orbit.color == 0 {
            for (i, &left) in orbit.output_cycles.iter().enumerate() {
                for &right in &orbit.output_cycles[i + 1..] {
                    for &(a, b, entry) in self.middle_diagonal {
                        colors[left] = a;
                        colors[right] = b;
                        self.expand(output, index + 1, colors, scale(weight, entry)?)?;
                    }
                    colors[right] = POINT;
                }
                colors[left] = POINT;
            }
        }
        Ok(())
    }
}

fn next_permutation(values: &mut [usize]) -> bool {
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

fn scale(weight: Integer, scalar: i128) -> Result<Integer> {
    weight.multiply(scalar).map_err(|_| coordinate_range())
}

fn conjugator_limit() -> Error {
    Error::plain("Nakajima cup product exceeds the conjugator work limit")
}
