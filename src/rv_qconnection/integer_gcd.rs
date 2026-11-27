use std::cmp::Ordering;

use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;

use super::{arithmetic, checked_divide, gcd, storage, Division, Polynomial, Result};

// If the primitive gcd fits i128, scaling it by the gcd of the input
// leading coefficients gives magnitudes at most 2^254. Reconstruction
// needs a modulus above 2^255 and one last prime of up to 32 bits.
const LIMBS: usize = 9;

#[derive(Clone, Copy, Default)]
struct Integer([u32; LIMBS]);

impl Integer {
    fn from_small(mut n: u128) -> Self {
        let mut value = Self::default();
        for digit in &mut value.0[..4] {
            *digit = n as u32;
            n >>= 32;
        }
        value
    }

    fn small(self) -> Option<u128> {
        if self.0[4..].iter().any(|&digit| digit != 0) {
            return None;
        }
        Some(
            self.0[..4]
                .iter()
                .rev()
                .fold(0, |acc, &digit| (acc << 32) | u128::from(digit)),
        )
    }

    fn signed(self, negative: bool) -> Option<i128> {
        let value = self.small()?;
        if negative && value == 1u128 << 127 {
            return Some(i128::MIN);
        }
        let value = i128::try_from(value).ok()?;
        Some(if negative { -value } else { value })
    }

    fn product(left: i128, right: i128) -> Self {
        let left = Self::from_small(left.unsigned_abs());
        let right = Self::from_small(right.unsigned_abs());
        let mut result = Self::default();
        for (i, &a) in left.0[..4].iter().enumerate() {
            let mut carry = 0;
            for (j, &b) in right.0[..4].iter().enumerate() {
                let digit = &mut result.0[i + j];
                carry += u64::from(*digit) + u64::from(a) * u64::from(b);
                *digit = carry as u32;
                carry >>= 32;
            }
            result.0[i + 4] = carry as u32;
        }
        result
    }

    fn bits(self) -> usize {
        self.0
            .iter()
            .rposition(|&digit| digit != 0)
            .map_or(0, |i| 32 * i + (32 - self.0[i].leading_zeros()) as usize)
    }

    fn cmp(self, other: Self) -> Ordering {
        self.0.iter().rev().cmp(other.0.iter().rev())
    }

    fn add_multiple(&mut self, other: Self, multiplier: u32) -> Result<()> {
        let mut carry = 0;
        for (digit, other) in self.0.iter_mut().zip(other.0) {
            carry += u64::from(*digit) + u64::from(other) * u64::from(multiplier);
            *digit = carry as u32;
            carry >>= 32;
        }
        if carry != 0 {
            return Err(arithmetic());
        }
        Ok(())
    }

    fn difference(mut self, other: Self) -> Self {
        let mut borrow = 0i64;
        for (digit, other) in self.0.iter_mut().zip(other.0) {
            let remaining = i64::from(*digit) - i64::from(other) - borrow;
            *digit = remaining as u32;
            borrow = i64::from(remaining < 0);
        }
        debug_assert_eq!(borrow, 0);
        self
    }

    fn residue(self, p: u32) -> u32 {
        self.0.iter().rev().fold(0, |r, &digit| {
            ((u64::from(r) << 32 | u64::from(digit)) % u64::from(p)) as u32
        })
    }

    fn divide(self, model: &mut Model, divisor: u128) -> Result<(Self, u128)> {
        model.spend(LIMBS)?;
        if let Some(value) = self.small() {
            return Ok((Self::from_small(value / divisor), value % divisor));
        }
        if divisor == 1 {
            return Ok((self, 0));
        }
        // The divisor is at most 2^127. Doubling a remainder and
        // appending a bit therefore fits in u128.
        let width = self.bits();
        model.spend(width)?;
        let mut quotient = Self::default();
        let mut remainder = 0;
        for bit in (0..width).rev() {
            let (i, j) = (bit / 32, bit % 32);
            remainder = (remainder << 1) | u128::from((self.0[i] >> j) & 1);
            if remainder >= divisor {
                remainder -= divisor;
                quotient.0[i] |= 1 << j;
            }
        }
        Ok((quotient, remainder))
    }
}

// A coefficient of a product is a sum of at most the storage limit
// (300,000) of products of two i128 values, so 288 bits also suffice
// for exact long division as long as every quotient coefficient fits.
pub(super) fn wide_divide(model: &mut Model, input: &[i128], divisor: &[i128]) -> Result<Division> {
    let (&lead, lower) = divisor.split_last().unwrap();
    let degree = lower.len();
    let length = input.len().saturating_sub(degree);
    storage(length)?;
    let mut quotient = vec![0; length];
    for (power, &coefficient) in input.iter().enumerate().rev() {
        let first = power.checked_sub(degree);
        let mut magnitude = Integer::from_small(coefficient.unsigned_abs());
        let mut negative = coefficient < 0;
        let end = quotient.len().min(power + 1);
        for (j, &b) in quotient
            .iter()
            .enumerate()
            .take(end)
            .skip(first.map_or(0, |k| k + 1))
        {
            let a = divisor[power - j];
            model.spend(1)?;
            if a == 0 || b == 0 {
                continue;
            }
            model.spend(16 + 2 * LIMBS)?;
            let term = Integer::product(a, b);
            let subtract_negative = (a < 0) == (b < 0);
            if negative == subtract_negative {
                magnitude.add_multiple(term, 1)?;
            } else if magnitude.cmp(term).is_lt() {
                magnitude = term.difference(magnitude);
                negative = subtract_negative;
            } else {
                magnitude = magnitude.difference(term);
            }
        }
        if let Some(k) = first {
            let (value, remainder) = magnitude.divide(model, lead.unsigned_abs())?;
            if remainder != 0 {
                return Ok(Division::NotExact);
            }
            let Some(value) = value.signed(negative != (lead < 0)) else {
                return Ok(Division::OutOfRange);
            };
            quotient[k] = value;
        } else if magnitude.bits() != 0 {
            return Ok(Division::NotExact);
        }
    }
    Ok(Division::Exact(quotient))
}

struct Reconstruction {
    modulus: Integer,
    residues: Vec<Integer>,
}

impl Reconstruction {
    fn new(size: usize) -> Result<Self> {
        storage(size.saturating_mul(LIMBS + 2))?;
        Ok(Self {
            modulus: Integer::from_small(1),
            residues: vec![Integer::default(); size],
        })
    }

    fn update(&mut self, model: &mut Model, p: Prime, gcd: &[u32], scale: u32) -> Result<()> {
        let pending = self.modulus.residue(p.get());
        let inverse = invert(model, p, pending)?;
        for (current, &coefficient) in self.residues.iter_mut().zip(gcd) {
            model.spend(2 * LIMBS)?;
            let v = p.mul(coefficient, scale);
            let delta = p.add(v, p.negate(current.residue(p.get())));
            current.add_multiple(self.modulus, p.mul(delta, inverse))?;
        }
        let mut modulus = Integer::default();
        modulus.add_multiple(self.modulus, p.get())?;
        self.modulus = modulus;
        Ok(())
    }

    fn candidate(&self, model: &mut Model, lead: u128) -> Result<Option<Polynomial>> {
        // With the leading coefficient unambiguously centered, the
        // content is the gcd of that coefficient and all other entries.
        let leading = Integer::from_small(lead);
        let mut twice_lead = leading;
        twice_lead.add_multiple(leading, 1)?;
        if self.modulus.cmp(twice_lead).is_le() {
            return Ok(None);
        }
        let centered = |representative| {
            let opposite = self.modulus.difference(representative);
            if opposite.cmp(representative).is_lt() {
                (opposite, true)
            } else {
                (representative, false)
            }
        };
        let mut content = lead;
        for &representative in &self.residues {
            if content == 1 {
                break;
            }
            model.spend(2 * LIMBS)?;
            let (magnitude, _) = centered(representative);
            let (_, residue) = magnitude.divide(model, content)?;
            content = gcd(content, residue);
        }
        let mut answer = Vec::with_capacity(self.residues.len());
        for &representative in &self.residues {
            model.spend(2 * LIMBS)?;
            let (magnitude, negative) = centered(representative);
            let (quotient, _) = magnitude.divide(model, content)?;
            let Some(coefficient) = quotient.signed(negative) else {
                return Ok(None);
            };
            answer.push(coefficient);
        }
        Ok(Some(answer))
    }
}

fn primes() -> impl FnMut(&mut Model) -> Result<Prime> {
    let mut next = 65_537u32;
    move |model| loop {
        let p = next;
        next = next.checked_add(2).ok_or_else(arithmetic)?;
        let mut factor = 3;
        while factor <= p / factor {
            model.spend(1)?;
            if p.is_multiple_of(factor) {
                break;
            }
            factor += 2;
        }
        if factor > p / factor {
            return Ok(Prime::validated(p));
        }
    }
}

fn residue(coefficient: i128, p: Prime) -> u32 {
    coefficient.rem_euclid(i128::from(p.get())) as u32
}

fn image(model: &mut Model, polynomial: &[i128], p: Prime) -> Result<Vec<u32>> {
    model.spend(polynomial.len())?;
    let mut mapped: Vec<_> = polynomial.iter().map(|&c| residue(c, p)).collect();
    trim(&mut mapped);
    Ok(mapped)
}

fn trim(polynomial: &mut Vec<u32>) {
    while polynomial.last() == Some(&0) {
        polynomial.pop();
    }
}

fn remainder(model: &mut Model, a: &mut Vec<u32>, b: &[u32], p: Prime) -> Result<()> {
    let inverse = invert(model, p, *b.last().unwrap())?;
    while a.len() >= b.len() {
        let shift = a.len() - b.len();
        let scalar = p.negate(p.mul(*a.last().unwrap(), inverse));
        model.spend(b.len())?;
        for (entry, &coefficient) in a[shift..].iter_mut().zip(b) {
            *entry = p.add(*entry, p.mul(scalar, coefficient));
        }
        trim(a);
    }
    Ok(())
}

fn image_gcd(model: &mut Model, a: &[i128], b: &[i128], p: Prime) -> Result<Vec<u32>> {
    let mut a = image(model, a, p)?;
    let mut b = image(model, b, p)?;
    while !b.is_empty() {
        remainder(model, &mut a, &b, p)?;
        (a, b) = (b, a);
    }
    let inverse = invert(model, p, *a.last().unwrap())?;
    model.spend(a.len())?;
    for coefficient in &mut a {
        *coefficient = p.mul(*coefficient, inverse);
    }
    Ok(a)
}

fn height_bits(polynomial: &[i128]) -> usize {
    polynomial
        .iter()
        .map(|x| x.unsigned_abs())
        .max()
        .map_or(0, |x| (u128::BITS - x.leading_zeros()) as usize)
}

fn divides(model: &mut Model, input: &[i128], divisor: &[i128]) -> Result<bool> {
    match checked_divide(model, input.to_vec(), divisor)? {
        Division::Exact(_) => return Ok(true),
        Division::NotExact => return Ok(false),
        Division::OutOfRange => {}
    }

    // Every pseudo-division step scales the coefficient bound by at
    // most twice the divisor height. Certify every coefficient of the
    // integral pseudo-remainder by congruences whose product exceeds
    // that bound. This also handles a large shifted verification quotient.
    let steps = input.len() + 1 - divisor.len();
    let bound = height_bits(input) + steps * (1 + height_bits(divisor));
    let mut next = primes();
    let mut proven_bits = 0;
    while proven_bits < bound {
        let p = next(model)?;
        if residue(*divisor.last().unwrap(), p) == 0 {
            continue;
        }
        let mut remainder_mod_p = image(model, input, p)?;
        let divisor_mod_p = image(model, divisor, p)?;
        remainder(model, &mut remainder_mod_p, &divisor_mod_p, p)?;
        if !remainder_mod_p.is_empty() {
            return Ok(false);
        }
        proven_bits += p.get().ilog2() as usize;
    }
    Ok(true)
}

fn reconstruct(
    model: &mut Model,
    a: &[i128],
    b: &[i128],
    mut next: impl FnMut(&mut Model) -> Result<Prime>,
) -> Result<Polynomial> {
    storage(a.len().saturating_add(b.len()))?;
    let (a_lead, b_lead) = (*a.last().unwrap(), *b.last().unwrap());
    let scale = gcd(a_lead.unsigned_abs(), b_lead.unsigned_abs());
    let mut size = a.len().min(b.len());
    let mut lift = Reconstruction::new(0)?;
    loop {
        let p = next(model)?;
        // Degree preservation makes the modular gcd degree an upper
        // bound. Unlucky primes can only increase that degree.
        if residue(a_lead, p) == 0 || residue(b_lead, p) == 0 {
            continue;
        }
        let image = image_gcd(model, a, b, p)?;
        if image.len() > size {
            continue;
        }
        if image.len() == 1 {
            return Ok(vec![1]);
        }
        if image.len() < size || lift.residues.is_empty() {
            size = image.len();
            lift = Reconstruction::new(size)?;
        }
        if lift.modulus.bits() > 255 {
            // Further primes of this same degree cannot help recover
            // an i128 primitive gcd. A smaller degree can reset the lift.
            continue;
        }
        lift.update(model, p, &image, (scale % u128::from(p.get())) as u32)?;
        if let Some(candidate) = lift.candidate(model, scale)? {
            if divides(model, a, &candidate)? && divides(model, b, &candidate)? {
                // Exact divisibility plus the modular upper bound
                // certifies the candidate's maximal degree.
                return Ok(candidate);
            }
        }
    }
}

pub(super) fn gcd_fallback(model: &mut Model, a: &[i128], b: &[i128]) -> Result<Polynomial> {
    reconstruct(model, a, b, primes())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rv_qconnection::{exact_divide, polynomial_gcd};

    #[test]
    fn crt_waits_for_a_degree_preserving_good_prime() {
        let constant = (1i128 << 100) + 123;
        let expected = [constant, 2];
        let a = [constant, constant + 2, 2];
        let b = [16 * constant, constant + 32, 2];
        let mut small = [2, 3, 5, 7, 11].into_iter();
        let mut ordinary = primes();
        let mut calls = 0;
        let next = |model: &mut Model| {
            calls += 1;
            small
                .next()
                .map_or_else(|| ordinary(model), |p| Ok(Prime::validated(p)))
        };
        // Two drops the leading terms. Three and five give an extra
        // common factor. Seven starts the lift but has too few bits.
        let result =
            reconstruct(&mut Model::new(), &a, &b, next).unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(result, expected);
        assert!(calls > 5);

        let mut next = [7, 11].into_iter();
        let coprime = reconstruct(&mut Model::new(), &[7, 1], &[7, 2], |_| {
            Ok(Prime::validated(next.next().unwrap()))
        });
        assert_eq!(coprime.unwrap_or_else(|error| panic!("{error}")), [1]);
        assert_eq!(next.next(), None);
    }

    #[test]
    fn exact_divisibility_handles_cancelling_out_of_range_products() {
        let top = i128::MAX;
        let divisor = [1, 1, 1];
        let input = vec![top, 0, -top, -top, 0, top];
        let model = &mut Model::new();
        assert!(matches!(
            checked_divide(model, input.clone(), &divisor)
                .unwrap_or_else(|error| panic!("{error}")),
            Division::OutOfRange
        ));
        assert!(divides(model, &input, &divisor).unwrap_or_else(|error| panic!("{error}")));
        let mut damaged = input.clone();
        damaged[0] -= 1;
        assert!(!divides(model, &damaged, &divisor).unwrap_or_else(|error| panic!("{error}")));
        assert_eq!(
            exact_divide(model, input, &divisor).unwrap_or_else(|error| panic!("{error}")),
            [top, -top, -top, top]
        );
        assert!(exact_divide(model, vec![i128::MIN], &[-1]).is_err());
    }

    #[test]
    fn primitive_gcd_removes_signed_contents_including_zero_inputs() {
        let model = &mut Model::new();
        let a = vec![-3, -9, -6];
        let b = vec![5, 15, 10];
        assert_eq!(
            polynomial_gcd(model, a.clone(), b).unwrap_or_else(|error| panic!("{error}")),
            [1, 3, 2]
        );
        assert_eq!(
            polynomial_gcd(model, Vec::new(), a).unwrap_or_else(|error| panic!("{error}")),
            [1, 3, 2]
        );
        assert!(polynomial_gcd(model, Vec::new(), Vec::new())
            .unwrap_or_else(|error| panic!("{error}"))
            .is_empty());
    }
}
