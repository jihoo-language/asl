use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::marker::PhantomData;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::wide_integer::Natural;

const COEFFICIENT_LIMBS: usize = 8;
pub(super) const WORK_LIMBS: usize = 26;
const TERM_LIMIT: usize = crate::complex::MAX_PIVOT_CELLS / (WORK_LIMBS + 3);
pub(crate) type Scalar = Wide<COEFFICIENT_LIMBS>;
type Denominator = Natural<COEFFICIENT_LIMBS>;
// Holds a triple product and the sums allowed by the sparse work limit.
pub(super) type Integer = Wide<WORK_LIMBS>;
pub(super) type Magnitude = Natural<WORK_LIMBS>;

pub(crate) trait Domain {
    const NAME: &'static str;
}

pub(crate) struct Coefficients<D> {
    // Increasing indices, nonzero numerators, and joint gcd one with the denominator.
    terms: Vec<(u32, Scalar)>,
    denominator: Denominator,
    domain: PhantomData<D>,
}

impl<D: Domain> Coefficients<D> {
    pub(crate) fn zero() -> Self {
        Self {
            terms: Vec::new(),
            denominator: Denominator::from_small(1),
            domain: PhantomData,
        }
    }

    pub(super) fn singleton(index: u32, denominator: u128) -> Self {
        Self {
            terms: vec![(index, Scalar::from_i128(1))],
            denominator: Denominator::from_small(denominator),
            domain: PhantomData,
        }
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        model.spend(1 + self.terms.len().min(other.terms.len()) * COEFFICIENT_LIMBS)?;
        let order = self
            .terms
            .iter()
            .zip(&other.terms)
            .map(|(&(key, a), &(other_index, b))| {
                key.cmp(&other_index)
                    .then_with(|| a.negative.cmp(&b.negative))
                    .then_with(|| a.magnitude.cmp(b.magnitude))
            })
            .find(|order| !order.is_eq());
        Ok(order.unwrap_or_else(|| {
            self.terms
                .len()
                .cmp(&other.terms.len())
                .then_with(|| self.denominator.cmp(other.denominator))
        }))
    }

    pub(crate) fn terms(&self) -> &[(u32, Scalar)] {
        &self.terms
    }

    pub(crate) fn denominator_text(&self) -> String {
        decimal(self.denominator)
    }

    pub(super) fn charge_storage(&self, model: &mut Model) -> Result<()> {
        model.allocate(3 + COEFFICIENT_LIMBS + (COEFFICIENT_LIMBS + 2) * self.terms.len())
    }

    pub(super) fn product_denominator<E>(&self, other: &Coefficients<E>) -> Result<Magnitude> {
        self.denominator
            .resize::<WORK_LIMBS>()
            .ok_or_else(Self::arithmetic)?
            .checked_product(other.denominator.resize().ok_or_else(Self::arithmetic)?)
            .ok_or_else(Self::arithmetic)
    }

    pub(crate) fn coefficient(&self, model: &mut Model, key: u32) -> Result<String> {
        model.spend(1 + self.terms.len().max(1).ilog2() as usize)?;
        let Ok(index) = self.terms.binary_search_by_key(&key, |&(key, _)| key) else {
            return Ok("0".into());
        };
        let coefficient = self.terms[index].1;
        let factor = coefficient.magnitude.gcd(model, self.denominator)?;
        let denominator = self
            .denominator
            .divide_exact(model, factor)?
            .ok_or_else(Self::arithmetic)?;
        let magnitude = coefficient
            .magnitude
            .divide_exact(model, factor)?
            .ok_or_else(Self::arithmetic)?;
        let numerator = signed_text(Scalar {
            magnitude,
            negative: coefficient.negative,
        });
        if denominator.small() == Some(1) {
            Ok(numerator)
        } else {
            Ok(format!("{numerator}/{}", decimal(denominator)))
        }
    }

    pub(crate) fn add(&self, model: &mut Model, other: &Self) -> Result<Self> {
        let common = self.denominator.gcd(model, other.denominator)?;
        let left = other
            .denominator
            .divide_exact(model, common)?
            .ok_or_else(Self::arithmetic)?;
        let right = self
            .denominator
            .divide_exact(model, common)?
            .ok_or_else(Self::arithmetic)?;
        let left: Magnitude = left.resize().ok_or_else(Self::arithmetic)?;
        let right: Magnitude = right.resize().ok_or_else(Self::arithmetic)?;
        let denominator = left
            .checked_product(self.denominator.resize().ok_or_else(Self::arithmetic)?)
            .ok_or_else(Self::arithmetic)?;
        model.spend(WORK_LIMBS * (self.terms.len() + other.terms.len()))?;
        let mut terms = BTreeMap::new();
        for (source, scale) in [(self, left), (other, right)] {
            for &(key, coefficient) in &source.terms {
                let magnitude = scale
                    .checked_product(
                        coefficient
                            .magnitude
                            .resize()
                            .ok_or_else(Self::arithmetic)?,
                    )
                    .ok_or_else(Self::arithmetic)?;
                Self::accumulate(
                    &mut terms,
                    key,
                    Integer {
                        magnitude,
                        negative: coefficient.negative,
                    },
                )?;
            }
        }
        Self::normalized(
            model,
            terms.iter().map(|(&key, entry)| (key, entry)),
            denominator,
        )
    }

    pub(crate) fn signed_map(
        &self,
        model: &mut Model,
        negate: impl Fn(u32) -> bool,
    ) -> Result<Self> {
        model.spend(self.terms.len())?;
        let terms = self
            .terms
            .iter()
            .map(|&(key, value)| (key, if negate(key) { value.negate() } else { value }))
            .collect();
        Ok(Self {
            terms,
            denominator: self.denominator,
            domain: PhantomData,
        })
    }

    pub(crate) fn scale(
        &self,
        model: &mut Model,
        numerator: i128,
        denominator: u128,
    ) -> Result<Self> {
        if denominator == 0 {
            return Err(Error::plain(format!(
                "{} scalar denominator must be positive",
                D::NAME
            )));
        }
        if numerator == 0 || self.terms.is_empty() {
            return Ok(Self::zero());
        }
        let denominator = self
            .denominator
            .resize::<WORK_LIMBS>()
            .ok_or_else(Self::arithmetic)?
            .multiply_u128(denominator)
            .ok_or_else(Self::arithmetic)?;
        model.spend(WORK_LIMBS * self.terms.len())?;
        let terms = self
            .terms
            .iter()
            .map(|&(key, coefficient)| {
                Ok((
                    key,
                    coefficient.resize::<WORK_LIMBS>()?.multiply(numerator)?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        Self::normalized(
            model,
            terms.iter().map(|(key, entry)| (*key, entry)),
            denominator,
        )
    }

    pub(super) fn normalized<'a>(
        model: &mut Model,
        coefficients: impl Iterator<Item = (u32, &'a Integer)> + Clone,
        denominator: Magnitude,
    ) -> Result<Self> {
        let mut factor = denominator;
        for (_, entry) in coefficients.clone() {
            if factor.small() == Some(1) {
                break;
            }
            factor = factor.gcd(model, entry.magnitude)?;
        }
        let divide = |model: &mut Model, value: Magnitude| -> Result<Denominator> {
            value
                .divide_exact(model, factor)?
                .and_then(|d| d.resize())
                .ok_or_else(Self::arithmetic)
        };
        let denominator = divide(model, denominator)?;
        let mut terms = Vec::new();
        for (key, entry) in coefficients {
            if !entry.is_zero() {
                let magnitude = divide(model, entry.magnitude)?;
                terms.push((
                    key,
                    Scalar {
                        magnitude,
                        negative: entry.negative,
                    },
                ));
            }
        }
        Ok(Self {
            terms,
            denominator,
            domain: PhantomData,
        })
    }

    pub(super) fn accumulate(
        terms: &mut BTreeMap<u32, Integer>,
        key: u32,
        coefficient: Integer,
    ) -> Result<()> {
        let sum = if let Some(&existing) = terms.get(&key) {
            existing.add(coefficient)?
        } else {
            if terms.len() >= TERM_LIMIT {
                return Err(Self::work_limit());
            }
            coefficient
        };
        if sum.is_zero() {
            terms.remove(&key);
        } else {
            terms.insert(key, sum);
        }
        Ok(())
    }

    pub(super) fn arithmetic() -> Error {
        Error::plain(format!(
            "{} rational coefficient exceeds the supported exact integer range",
            D::NAME
        ))
    }

    pub(super) fn work_limit() -> Error {
        Error::plain(format!(
            "{} operation exceeds the sparse work or storage limit",
            D::NAME
        ))
    }
}

pub(super) fn decimal<const N: usize>(mut magnitude: Natural<N>) -> String {
    if let Some(value) = magnitude.small() {
        return value.to_string();
    }
    let mut digits = Vec::new();
    while !magnitude.is_zero() {
        let (quotient, remainder) = magnitude.divide_u32(1_000_000_000);
        digits.push(remainder);
        magnitude = quotient;
    }
    let mut output = digits.pop().expect("positive magnitude").to_string();
    for digit in digits.into_iter().rev() {
        write!(output, "{digit:09}").unwrap();
    }
    output
}

pub(crate) fn signed_text(value: Scalar) -> String {
    let magnitude = decimal(value.magnitude);
    if value.negative {
        format!("-{magnitude}")
    } else {
        magnitude
    }
}
