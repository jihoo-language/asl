use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_ore::gcd;
use crate::syntax::ParsedExpr;
use crate::wide_integer::Natural;

use super::super::unimodular::{self, Vector, RANK};
use super::{grade_range, Accumulator, Factor, Integer, Mode, State};

pub(crate) enum Command {
    ChernDenominator(u128, u128, ParsedExpr),
    TopSegre(u128, ParsedExpr),
}

pub(crate) fn denominator(model: &mut Model, line: &Vector, n: u128, j: u128) -> Result<u128> {
    calculate(model, line, n, j, None).map(|coefficients| coefficients.denominator)
}

pub(crate) fn chern(
    model: &mut Model,
    line: &Vector,
    n: u128,
    j: u128,
    scale: i128,
) -> Result<State> {
    calculate(model, line, n, j, Some(scale))?.output.finish()
}

pub(crate) fn top_segre(model: &mut Model, line: &Vector, n: u128) -> Result<i128> {
    if n == 0 {
        return Ok(1);
    }
    let square = line.wide_pair(model, line)?;
    let (half_square, parity) = square.magnitude.divide_u32(2);
    debug_assert_eq!(parity, 0, "the marked K3 lattice is even");
    let length = unimodular::Integer {
        magnitude: Natural::from_small(n),
        negative: false,
    };
    let top = unimodular::Integer {
        magnitude: half_square,
        negative: square.negative,
    }
    .add(length.multiply(-2)?)?
    .add(unimodular::Integer::from_i128(2))?;
    let negative_top = top.negative && !top.is_zero();
    if !negative_top && top.magnitude.cmp(length.magnitude).is_lt() {
        return Ok(0);
    }
    if n >= 128 {
        return Err(segre_range());
    }
    // For top < 0 use binom(top,n) = (-1)^n binom(n-1-top,n).
    let mut factor = if negative_top {
        top.negate()
            .add(unimodular::Integer::from_i128(n as i128 - 1))?
            .magnitude
    } else {
        top.magnitude
    };
    let complement = factor.difference(length.magnitude).small();
    let steps = complement.map_or(n, |other| n.min(other));
    let mut value = Natural::from_small(1);
    for k in 1..=steps as u32 {
        model.spend(40)?;
        let product = value.checked_product(factor).ok_or_else(segre_range)?;
        let (quotient, remainder) = product.divide_u32(k);
        debug_assert_eq!(remainder, 0);
        value = quotient;
        // These binomial magnitudes are nondecreasing up to the requested k.
        if value.bits() as u128 + n > 128 {
            return Err(segre_range());
        }
        factor = factor.difference(Natural::from_small(1));
    }
    value
        .multiply_u128(1u128 << n)
        .and_then(|value| value.signed(negative_top && n % 2 != 0))
        .ok_or_else(segre_range)
}

#[derive(Clone, Copy)]
struct Part {
    length: u128,
    color: u8,
    count: u128,
}

impl Part {
    fn unit(length: u128, count: u128) -> Self {
        Self {
            length,
            color: 0,
            count,
        }
    }
}

struct Coefficients {
    scale: Option<i128>,
    denominator: u128,
    output: Accumulator,
}

fn calculate(
    model: &mut Model,
    line: &Vector,
    n: u128,
    j: u128,
    scale: Option<i128>,
) -> Result<Coefficients> {
    let mut coefficients = Coefficients {
        scale,
        denominator: 1,
        output: Accumulator::default(),
    };
    if scale == Some(0) || j > n {
        return Ok(coefficients);
    }
    let units = n - j;
    if units != 0 {
        // A purely unit-colored term has numerator ±1. Check a term with
        // one long part and all remaining parts of length one before enumeration.
        let witness = if j == 0 {
            vec![Part::unit(1, units)]
        } else {
            let mut parts = vec![Part::unit(j + 1, 1)];
            if units > 1 {
                parts.push(Part::unit(1, units - 1));
            }
            parts
        };
        coefficients.require_denominator(unit_denominator(&witness)?)?;
    }
    let coordinates = if j == 0 {
        [0; RANK]
    } else {
        line.basis_coordinates(model)?
    };
    let colors: Vec<_> = (1..=RANK as u8)
        .filter(|&color| coordinates[usize::from(color - 1)] != 0)
        .collect();
    let primes = if colors.is_empty() {
        Vec::new()
    } else {
        // Already one middle color gives all partitions built from lengths
        // 1,2,3,4, with up to floor(j/9) copies independently of 2,3,4.
        let lower = (j / 9 + 1).saturating_pow(3);
        reserve(model, lower)?;
        let n = usize::try_from(n).map_err(|_| enumeration_limit())?;
        primes_up_to(model, n)?
    };
    let mut generator = Generator {
        coordinates,
        colors,
        primes,
        coefficients,
    };
    let mut key = Vec::new();
    if units == 0 {
        if n == 0 || !generator.colors.is_empty() {
            generator.middle(
                model,
                n,
                n,
                generator.colors.len().saturating_sub(1),
                &mut key,
            )?;
        }
    } else if generator.colors.is_empty() {
        generator.units(model, n, units, n - units + 1, 0, &mut key)?;
    } else {
        for unit_total in units..=n {
            generator.units(
                model,
                unit_total,
                units,
                unit_total - units + 1,
                n - unit_total,
                &mut key,
            )?;
        }
    }
    Ok(generator.coefficients)
}

struct Generator {
    coordinates: [i128; RANK],
    colors: Vec<u8>,
    primes: Vec<u128>,
    coefficients: Coefficients,
}

impl Generator {
    fn units(
        &mut self,
        model: &mut Model,
        total: u128,
        count: u128,
        max: u128,
        middle: u128,
        key: &mut Vec<Part>,
    ) -> Result<()> {
        if count == 0 {
            debug_assert_eq!(total, 0);
            return self.middle(
                model,
                middle,
                middle,
                self.colors.len().saturating_sub(1),
                key,
            );
        }
        if count == 1 || total == count {
            key.push(Part::unit(total / count, count));
            self.units(model, 0, 0, 0, middle, key)?;
            key.pop();
            return Ok(());
        }
        let lower = total.div_ceil(count);
        let upper = max.min(total - count + 1);
        reserve(model, upper.saturating_sub(lower) + 1)?;
        for length in (lower..=upper).rev() {
            // All remaining parts lie between 1 and length-1.
            let least = total
                .saturating_sub(count.saturating_mul(length - 1))
                .max(1);
            let most = count.min((total - count) / (length - 1));
            for copies in least..=most {
                model.spend(1)?;
                key.push(Part::unit(length, copies));
                self.units(
                    model,
                    total - length * copies,
                    count - copies,
                    length - 1,
                    middle,
                    key,
                )?;
                key.pop();
            }
        }
        Ok(())
    }

    fn middle(
        &mut self,
        model: &mut Model,
        total: u128,
        max_length: u128,
        max_color: usize,
        key: &mut Vec<Part>,
    ) -> Result<()> {
        if total == 0 {
            return self.emit(model, key);
        }
        let last = self.colors.len() - 1;
        for length in (1..=total.min(max_length)).rev() {
            let upper = if length == max_length {
                max_color
            } else {
                last
            };
            for index in (0..=upper).rev() {
                let color = self.colors[index];
                let next = if index == 0 {
                    (length - 1, last)
                } else {
                    (length, index - 1)
                };
                let high = total / length;
                let low = if length == 1 && index == 0 { high } else { 1 };
                for count in low..=high {
                    model.spend(1)?;
                    key.push(Part {
                        length,
                        color,
                        count,
                    });
                    self.middle(model, total - length * count, next.0, next.1, key)?;
                    key.pop();
                }
            }
        }
        Ok(())
    }

    fn emit(&mut self, model: &mut Model, parts: &[Part]) -> Result<()> {
        model.spend(1 + parts.len())?;
        let coefficient = if parts.iter().all(|part| part.color == 0) {
            let denominator = unit_denominator(parts)?;
            self.coefficients.require_denominator(denominator)?;
            self.coefficients
                .scale
                .map(|scale| scale.unsigned_abs() / denominator)
        } else {
            self.colored_magnitude(model, parts)?
        };
        if let (Some(magnitude), Some(scale)) = (coefficient, self.coefficients.scale) {
            let mut negative = scale < 0;
            let mut key = Vec::with_capacity(parts.len());
            for &Part {
                length,
                color,
                count,
            } in parts
            {
                if count % 2 == 1 {
                    negative ^= length % 2 == 0;
                    negative ^= color != 0 && self.coordinates[usize::from(color - 1)] < 0;
                }
                let length = i128::try_from(length).map_err(|_| grade_range())?;
                key.push(Factor {
                    mode: Mode { length, color },
                    exponent: count,
                });
            }
            key.sort_unstable();
            let coefficient = Natural::<16>::from_small(magnitude)
                .signed(negative)
                .ok_or_else(chern_range)?;
            self.coefficients
                .output
                .add(key, Integer::from_i128(coefficient))?;
        }
        Ok(())
    }

    fn colored_magnitude(&mut self, model: &mut Model, parts: &[Part]) -> Result<Option<u128>> {
        let mut powers = [0; RANK];
        let mut max = 1;
        for &Part {
            length,
            color,
            count,
        } in parts
        {
            max = max.max(length).max(count);
            if color != 0 {
                powers[usize::from(color - 1)] += count;
            }
        }
        let mut bases = self.coordinates.map(i128::unsigned_abs);
        let mut denominator = 1u128;
        let mut surplus = Vec::new();
        // Cancel prime valuations before forming any powers, so the raw
        // numerator and centralizer may individually exceed the output range.
        for &prime in self.primes.iter().take_while(|&&prime| prime <= max) {
            model.spend(1 + parts.len() + self.colors.len())?;
            let needed: u128 = parts
                .iter()
                .map(
                    |&Part {
                         mut length, count, ..
                     }| {
                        count * remove_prime(&mut length, prime) + factorial_valuation(count, prime)
                    },
                )
                .sum();
            if needed == 0 {
                continue;
            }
            let mut supplied = 0;
            for &color in &self.colors {
                let index = usize::from(color - 1);
                if powers[index] != 0 {
                    supplied += powers[index] * remove_prime(&mut bases[index], prime);
                }
            }
            if needed > supplied {
                denominator = multiply_power(denominator, prime, needed - supplied)
                    .ok_or_else(denominator_range)?;
            } else if self.coefficients.scale.is_some() && supplied > needed {
                surplus.push((prime, supplied - needed));
            }
        }
        self.coefficients.require_denominator(denominator)?;
        let Some(scale) = self.coefficients.scale else {
            return Ok(None);
        };
        let mut value = scale.unsigned_abs() / denominator;
        for (base, power) in bases.into_iter().zip(powers).chain(surplus) {
            value = multiply_power(value, base, power).ok_or_else(chern_range)?;
        }
        Ok(Some(value))
    }
}

impl Coefficients {
    fn require_denominator(&mut self, denominator: u128) -> Result<()> {
        if let Some(scale) = self.scale {
            if !scale.unsigned_abs().is_multiple_of(denominator) {
                return Err(Error::plain(
                    "tautological Chern scale does not clear the denominator",
                ));
            }
        } else {
            self.denominator = (self.denominator / gcd(self.denominator, denominator))
                .checked_mul(denominator)
                .ok_or_else(denominator_range)?;
        }
        Ok(())
    }
}

fn unit_denominator(parts: &[Part]) -> Result<u128> {
    let mut value = 1;
    for &Part { length, count, .. } in parts {
        value = multiply_power(value, length, count).ok_or_else(denominator_range)?;
        value = (2..=count)
            .try_fold(value, u128::checked_mul)
            .ok_or_else(denominator_range)?;
    }
    Ok(value)
}

fn multiply_power(mut value: u128, mut base: u128, mut exponent: u128) -> Option<u128> {
    while exponent != 0 {
        if exponent % 2 == 1 {
            value = value.checked_mul(base)?;
        }
        exponent /= 2;
        if exponent != 0 {
            base = base.checked_mul(base)?;
        }
    }
    Some(value)
}

fn remove_prime(value: &mut u128, prime: u128) -> u128 {
    let mut count = 0;
    while value.is_multiple_of(prime) {
        *value /= prime;
        count += 1;
    }
    count
}

fn factorial_valuation(mut n: u128, prime: u128) -> u128 {
    let mut result = 0;
    while n >= prime {
        n /= prime;
        result += n;
    }
    result
}

fn primes_up_to(model: &mut Model, n: usize) -> Result<Vec<u128>> {
    model.spend(4 * n)?;
    let mut composite = vec![false; n + 1];
    let mut primes = Vec::new();
    for p in 2..=n {
        if !composite[p] {
            primes.push(p as u128);
            if p <= n / p {
                for multiple in (p * p..=n).step_by(p) {
                    composite[multiple] = true;
                }
            }
        }
    }
    Ok(primes)
}

fn reserve(model: &mut Model, lower: u128) -> Result<()> {
    model
        .spend(usize::try_from(lower).map_err(|_| enumeration_limit())?)
        .map_err(|_| enumeration_limit())
}

fn enumeration_limit() -> Error {
    Error::plain("tautological Chern partition enumeration exceeds the computation limit")
}

fn chern_range() -> Error {
    Error::plain("tautological Chern coefficient exceeds the signed 128-bit range")
}

fn denominator_range() -> Error {
    Error::plain("tautological Chern denominator exceeds the unsigned 128-bit range")
}

fn segre_range() -> Error {
    Error::plain("top tautological Segre number exceeds the signed 128-bit range")
}
