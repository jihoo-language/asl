use std::cmp::Ordering;
use std::rc::Rc;

use super::nikulin::{self, Action, ORDER};
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_lens::{self, Rational};
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};
use crate::wide_integer::Natural;

const SIZE: usize = ORDER as usize;
const RANK: usize = 16;
// A common-denominator sum can need one bit beyond two 128-bit factors.
type Integer = Wide<9>;
type Denominator = Natural<9>;

pub(crate) enum Expression {
    Basis(Box<ParsedExpr>),
    Apply(Operation, Box<ParsedExpr>),
    Frobenius(Label, usize, Box<ParsedExpr>),
    Shell(Label, usize),
}

pub(crate) struct Coefficient(pub(crate) Box<[ParsedExpr; 2]>);

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Scale(i128, u128),
    S,
    T,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Scale(..) => "rvKummerWeilScale",
            Self::S => "rvKummerWeilS",
            Self::T => "rvKummerWeilT",
        }
    }
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Vector(Rational);

impl Vector {
    pub(crate) fn zero() -> Self {
        Self(Rational::zero())
    }

    pub(crate) fn basis(model: &mut Model, mask: u8) -> Result<Self> {
        let mut coefficients = vec![0; usize::from(mask) + 1];
        coefficients[usize::from(mask)] = 1;
        Ok(Self(rv_lens::fraction(model, &coefficients, 1)?))
    }

    pub(crate) fn fraction(&self) -> &Rational {
        &self.0
    }

    pub(crate) fn coefficient(&self, model: &mut Model, index: u8) -> Result<Rational> {
        let value = self.numerator(usize::from(index));
        rv_lens::fraction(model, &[value], self.0.denominator)
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(4 + 2 * self.0.numerator.len())?;
        Ok(Value::KummerWeil(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        rv_lens::compare(model, &self.0, &other.0)
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let common = rv_lens::gcd(self.0.denominator, rhs.0.denominator);
        let left = Denominator::from_small(rhs.0.denominator / common);
        let right = Denominator::from_small(self.0.denominator / common);
        let denominator = left
            .multiply_u128(self.0.denominator)
            .ok_or_else(arithmetic)?;
        let length = self.0.numerator.len().max(rhs.0.numerator.len());
        model.spend(8 * length)?;
        let coefficients = (0..length)
            .map(|i| {
                unsigned_product(self.numerator(i), left)?
                    .add(unsigned_product(rhs.numerator(i), right)?)
            })
            .collect::<Result<Vec<_>>>()?;
        Self::normalized(model, &coefficients, denominator)
    }

    pub(crate) fn apply(&self, model: &mut Model, operation: Operation) -> Result<Self> {
        match operation {
            Operation::Scale(n, d) => self.scale(model, n, d),
            Operation::S => self.fourier(model),
            Operation::T => {
                model.spend(6 * self.0.numerator.len())?;
                let coefficients = self
                    .0
                    .numerator
                    .iter()
                    .enumerate()
                    .map(|(mask, &value)| {
                        if nikulin::quadratic(mask as u8) == 0 {
                            Ok(value)
                        } else {
                            value.checked_neg().ok_or_else(arithmetic)
                        }
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(Self(Rational {
                    numerator: coefficients.into(),
                    denominator: self.0.denominator,
                }))
            }
        }
    }

    pub(crate) fn frobenius(&self, model: &mut Model, action: &Action) -> Result<Self> {
        model.spend(6 * SIZE)?;
        let mut coefficients = [0; SIZE];
        for (i, &value) in self.0.numerator.iter().enumerate() {
            coefficients[usize::from(action.apply(i as u8))] = value;
        }
        Ok(Self(rv_lens::fraction(
            model,
            &coefficients,
            self.0.denominator,
        )?))
    }

    fn numerator(&self, index: usize) -> i128 {
        self.0.numerator.get(index).copied().unwrap_or(0)
    }

    fn scale(&self, model: &mut Model, n: i128, d: u128) -> Result<Self> {
        if d == 0 {
            return Err(Error::plain(
                "Kummer Weil scalar denominator must be positive",
            ));
        }
        if n == 0 || self.0.is_zero() {
            return Ok(Self::zero());
        }
        let denominator = Denominator::from_small(self.0.denominator)
            .multiply_u128(d)
            .ok_or_else(arithmetic)?;
        model.spend(4 * self.0.numerator.len())?;
        let coefficients: Vec<_> = self
            .0
            .numerator
            .iter()
            .map(|&v| Integer::product(v, n))
            .collect();
        Self::normalized(model, &coefficients, denominator)
    }

    fn fourier(&self, model: &mut Model) -> Result<Self> {
        if self.0.is_zero() {
            return Ok(Self::zero());
        }
        model.spend(12 * SIZE)?;
        let mut transformed: [_; SIZE] =
            std::array::from_fn(|i| Integer::from_i128(self.numerator(i)));
        // First compute the ordinary binary Walsh transform. Reindex the dual
        // characters by the discriminant pairing before dividing by sqrt(64).
        for bit in 0..6 {
            let half = 1 << bit;
            for block in transformed.chunks_exact_mut(2 * half) {
                let (left, right) = block.split_at_mut(half);
                for (a, b) in left.iter_mut().zip(right) {
                    (*a, *b) = (a.add(*b)?, a.add(b.negate())?);
                }
            }
        }
        let coefficients: [_; SIZE] = std::array::from_fn(|y| {
            let character = (0..6).fold(0, |mask, i| {
                mask | ((nikulin::pair(1 << i, y as u8) as usize) << i)
            });
            transformed[character]
        });
        let denominator = Denominator::from_small(self.0.denominator)
            .multiply_u128(8)
            .ok_or_else(arithmetic)?;
        Self::normalized(model, &coefficients, denominator)
    }

    fn normalized(
        model: &mut Model,
        coefficients: &[Integer],
        denominator: Denominator,
    ) -> Result<Self> {
        let mut factor = denominator;
        for entry in coefficients {
            if factor.small() == Some(1) {
                break;
            }
            factor = factor.gcd(model, entry.magnitude)?;
        }
        let divide = |model: &mut Model, value: Denominator| -> Result<Denominator> {
            value.divide_exact(model, factor)?.ok_or_else(arithmetic)
        };
        let denominator = divide(model, denominator)?.small().ok_or_else(arithmetic)?;
        let numerator = coefficients
            .iter()
            .map(|entry| {
                divide(model, entry.magnitude)?
                    .signed(entry.negative)
                    .ok_or_else(arithmetic)
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self(rv_lens::fraction(model, &numerator, denominator)?))
    }
}

fn unsigned_product(value: i128, factor: Denominator) -> Result<Integer> {
    Ok(Integer {
        magnitude: factor
            .multiply_u128(value.unsigned_abs())
            .ok_or_else(arithmetic)?,
        negative: value < 0,
    })
}

fn arithmetic() -> Error {
    Error::plain("Kummer Weil rational arithmetic exceeds the supported exact integer range")
}

pub(crate) fn theta_shell(model: &mut Model, norm: usize) -> Result<Vector> {
    let exponent = norm.checked_mul(2).ok_or_else(shell_limit)?;
    let width = exponent.checked_add(1).ok_or_else(shell_limit)?;
    let root = exponent.isqrt();
    let storage = width.saturating_mul(8 * (RANK + 1));
    if storage > crate::complex::MAX_PIVOT_CELLS {
        return Err(shell_limit());
    }
    let work = width.saturating_mul(2 * RANK * (root + 1) + RANK + 1) + SIZE * 32;
    model.spend(work).map_err(|_| shell_limit())?;
    let (even, odd): (Vec<_>, Vec<_>) = (0..=root).partition(|j| j.is_multiple_of(2));
    let even = powers(exponent, &even)?;
    let odd = powers(exponent, &odd)?;
    let mut by_weight = [0i128; RANK + 1];
    for (weight, count) in by_weight.iter_mut().enumerate() {
        for (&a, &b) in even[RANK - weight].iter().zip(odd[weight].iter().rev()) {
            let term = a.checked_mul(b).ok_or_else(arithmetic)?;
            *count = count.checked_add(term).ok_or_else(arithmetic)?;
        }
    }
    let mut coefficients = [0i128; SIZE];
    for (mask, coefficient) in coefficients.iter_mut().enumerate() {
        let table = nikulin::quadratic_table(mask as u8);
        let mut multiplicities = [0i128; RANK + 1];
        for affine in nikulin::affine_code() {
            multiplicities[(table ^ affine).count_ones() as usize] += 1;
        }
        for (&count, &multiplicity) in by_weight.iter().zip(&multiplicities) {
            let term = count.checked_mul(multiplicity).ok_or_else(arithmetic)?;
            *coefficient = coefficient.checked_add(term).ok_or_else(arithmetic)?;
        }
    }
    Ok(Vector(rv_lens::fraction(model, &coefficients, 1)?))
}

fn powers(exponent: usize, coordinates: &[usize]) -> Result<Vec<Vec<i128>>> {
    let mut identity = vec![0i128; exponent + 1];
    identity[0] = 1;
    let mut powers = vec![identity];
    for power in 0..RANK {
        let mut next = vec![0i128; exponent + 1];
        for &coordinate in coordinates {
            let square = coordinate * coordinate;
            let multiplicity = if coordinate == 0 { 1 } else { 2 };
            for (entry, &old) in next[square..].iter_mut().zip(&powers[power]) {
                let term = old.checked_mul(multiplicity).ok_or_else(arithmetic)?;
                *entry = entry.checked_add(term).ok_or_else(arithmetic)?;
            }
        }
        powers.push(next);
    }
    Ok(powers)
}

fn shell_limit() -> Error {
    Error::plain("Kummer theta shell exceeds the polynomial work or storage limit")
}
