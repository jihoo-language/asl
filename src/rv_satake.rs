use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_sphhecke;

const MAX_CELLS: usize = 300_000;

type Exponents = Rc<[i128]>;

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Term {
    pub(crate) exponents: Exponents,
    pub(crate) power: i128,
    pub(crate) coefficient: i128,
}

pub(crate) type Support = Rc<[Term]>;

fn arithmetic() -> Error {
    Error::plain("Satake Laurent arithmetic is out of range")
}

fn storage() -> Error {
    Error::plain("Satake Laurent calculation exceeds the storage limit")
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(storage());
    }
    Ok(())
}

fn term_cells(rank: usize) -> Result<usize> {
    rank.checked_mul(2)
        .and_then(|result| result.checked_add(5))
        .ok_or_else(storage)
}

fn persist(model: &mut Model, terms: Vec<Term>) -> Result<Value> {
    let mut cells = 0;
    for term in &terms {
        cells += term_cells(term.exponents.len())?;
        check_storage(cells)?;
    }
    model.allocate(cells)?;
    Ok(Value::RVSatake(terms.into()))
}

pub(crate) fn zero() -> Value {
    Value::RVSatake([].into())
}

pub(crate) fn scalar(model: &mut Model, rank: usize, power: i128) -> Result<Value> {
    check_storage(term_cells(rank)?)?;
    model.spend(rank)?;
    persist(
        model,
        vec![Term {
            exponents: vec![0; rank].into(),
            power,
            coefficient: 1,
        }],
    )
}

pub(crate) fn compare(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Ordering> {
    for (a, b) in left.iter().zip(right) {
        model.spend(1 + a.exponents.len().min(b.exponents.len()))?;
        let order = a.cmp(b);
        if order != Ordering::Equal {
            return Ok(order);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

#[derive(Default)]
struct Sum {
    terms: BTreeMap<(Exponents, i128), i128>,
    cells: usize,
}

impl Sum {
    fn add(
        &mut self,
        model: &mut Model,
        exponents: &Exponents,
        power: i128,
        coefficient: i128,
        retained: usize,
    ) -> Result<()> {
        let cells = term_cells(exponents.len())?;
        model.spend(cells * (1 + self.terms.len().max(1).ilog2() as usize))?;
        let key = (exponents.clone(), power);
        if let Some(old) = self.terms.get_mut(&key) {
            *old = old.checked_add(coefficient).ok_or_else(arithmetic)?;
            if *old == 0 {
                self.terms.remove(&key);
                self.cells -= cells;
            }
        } else if coefficient != 0 {
            self.cells += cells;
            check_storage(self.cells + retained)?;
            self.terms.insert(key, coefficient);
        }
        Ok(())
    }

    fn finish(self) -> Vec<Term> {
        self.terms
            .into_iter()
            .map(|((exponents, power), coefficient)| Term {
                exponents,
                power,
                coefficient,
            })
            .collect()
    }
}

pub(crate) fn add(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Value> {
    let mut output = Sum::default();
    for term in left.iter().chain(right) {
        output.add(model, &term.exponents, term.power, term.coefficient, 0)?;
    }
    persist(model, output.finish())
}

pub(crate) fn neg(model: &mut Model, terms: &[Term]) -> Result<Value> {
    model.spend(terms.len())?;
    let output = terms
        .iter()
        .map(|term| {
            Ok(Term {
                coefficient: term.coefficient.checked_neg().ok_or_else(arithmetic)?,
                ..term.clone()
            })
        })
        .collect::<Result<Vec<_>>>()?;
    persist(model, output)
}

pub(crate) fn multiply(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Value> {
    let mut output = Sum::default();
    for a in left {
        for b in right {
            model.spend(1 + a.exponents.len())?;
            let exponents = a
                .exponents
                .iter()
                .zip(b.exponents.iter())
                .map(|(&x, &y)| x.checked_add(y).ok_or_else(arithmetic))
                .collect::<Result<Vec<_>>>()?;
            let power = a.power.checked_add(b.power).ok_or_else(arithmetic)?;
            let coefficient = a
                .coefficient
                .checked_mul(b.coefficient)
                .ok_or_else(arithmetic)?;
            output.add(model, &exponents.into(), power, coefficient, 0)?;
        }
    }
    persist(model, output.finish())
}

fn integer_sum(values: &[i128]) -> Result<i128> {
    values.iter().try_fold(0i128, |sum, &value| {
        sum.checked_add(value).ok_or_else(arithmetic)
    })
}

fn weight(partition: &[i128]) -> Result<i128> {
    let rank = partition.len();
    partition
        .windows(2)
        .enumerate()
        .try_fold(0i128, |sum, (index, pair)| {
            let multiplicity = (index + 1) as i128 * (rank - index - 1) as i128;
            let term = pair[0]
                .checked_sub(pair[1])
                .and_then(|gap| gap.checked_mul(multiplicity))
                .ok_or_else(arithmetic)?;
            sum.checked_add(term).ok_or_else(arithmetic)
        })
}

pub(crate) fn transform(model: &mut Model, input: &[rv_sphhecke::Term]) -> Result<Value> {
    let mut cache = HallLittlewood::default();
    let mut output = Sum::default();
    for source in input {
        let coweight = &source.coweight;
        let center = coweight.first().copied().unwrap_or(0);
        model.spend(coweight.len())?;
        let partition = coweight
            .iter()
            .rev()
            .map(|value| value.checked_sub(center).ok_or_else(arithmetic))
            .collect::<Result<Vec<_>>>()?;
        let normalization = weight(&partition)?;
        let image = cache.expand(model, &partition, output.cells)?;
        for term in image.iter() {
            model.spend(term.exponents.len())?;
            let exponents = term
                .exponents
                .iter()
                .map(|value| value.checked_add(center).ok_or_else(arithmetic))
                .collect::<Result<Vec<_>>>()?
                .into();
            let normalized = term
                .power
                .checked_add(normalization)
                .ok_or_else(arithmetic)?;
            for &(degree, scalar) in source.coefficient.iter() {
                let power = normalized
                    .checked_add(2 * degree as i128)
                    .ok_or_else(arithmetic)?;
                let coefficient = term
                    .coefficient
                    .checked_mul(scalar)
                    .ok_or_else(arithmetic)?;
                output.add(model, &exponents, power, coefficient, cache.cells)?;
            }
        }
    }
    persist(model, output.finish())
}

#[derive(Default)]
struct HallLittlewood {
    cache: BTreeMap<Exponents, Support>,
    cells: usize,
}

impl HallLittlewood {
    fn expand(&mut self, model: &mut Model, lambda: &[i128], retained: usize) -> Result<Support> {
        model.spend(lambda.len() * (1 + self.cache.len().max(1).ilog2() as usize))?;
        if let Some(polynomial) = self.cache.get(lambda) {
            return Ok(polynomial.clone());
        }
        let rank = lambda.len();
        let width = term_cells(rank)?;
        check_storage(width + retained + self.cells)?;
        let polynomial = if lambda.first() == lambda.last() {
            let term = Term {
                exponents: lambda.into(),
                power: 0,
                coefficient: 1,
            };
            self.cells += width;
            vec![term]
        } else {
            // Every permutation of the leading monomial occurs. This also
            // bounds the recursion rank before building interlacing partitions.
            check_storage(width.checked_mul(rank).ok_or_else(storage)?)?;
            let mut choices = 1usize;
            for pair in lambda.windows(2) {
                let width = pair[0]
                    .checked_sub(pair[1])
                    .and_then(|gap| gap.checked_add(1))
                    .and_then(|width| usize::try_from(width).ok())
                    .ok_or_else(storage)?;
                choices = choices.checked_mul(width).ok_or_else(storage)?;
                check_storage(choices)?;
            }
            // At v^-2=0 each child contributes its own leading monomial.
            check_storage(width.checked_mul(choices).ok_or_else(storage)?)?;
            model.spend(choices)?;
            let mut child = lambda[1..].to_vec();
            let total = integer_sum(lambda)?;
            let mut result = Sum::default();
            loop {
                let power = total
                    .checked_sub(integer_sum(&child)?)
                    .ok_or_else(arithmetic)?;
                let factor = branching_factor(model, lambda, &child)?;
                let polynomial = self.expand(model, &child, retained + result.cells)?;
                for term in polynomial.iter() {
                    model.spend(rank)?;
                    let exponents: Exponents = term
                        .exponents
                        .iter()
                        .copied()
                        .chain(std::iter::once(power))
                        .collect();
                    for &(power, coefficient) in &factor {
                        let power = power.checked_add(term.power).ok_or_else(arithmetic)?;
                        let coefficient = coefficient
                            .checked_mul(term.coefficient)
                            .ok_or_else(arithmetic)?;
                        result.add(model, &exponents, power, coefficient, retained + self.cells)?;
                    }
                }

                let Some(index) = (0..child.len()).rfind(|&index| child[index] < lambda[index])
                else {
                    break;
                };
                model.spend(child.len() - index)?;
                child[index] += 1;
                child[index + 1..].copy_from_slice(&lambda[index + 2..]);
            }
            self.cells += result.cells;
            result.finish()
        };
        self.cells += 1 + 2 * rank;
        check_storage(self.cells + retained)?;
        let result: Support = polynomial.into();
        self.cache.insert(lambda.into(), result.clone());
        Ok(result)
    }
}

fn add_scalar(map: &mut BTreeMap<i128, i128>, power: i128, coefficient: i128) -> Result<()> {
    let value = map
        .get(&power)
        .copied()
        .unwrap_or(0)
        .checked_add(coefficient)
        .ok_or_else(arithmetic)?;
    if value == 0 {
        map.remove(&power);
    } else {
        map.insert(power, value);
    }
    Ok(())
}

fn branching_factor(
    model: &mut Model,
    lambda: &[i128],
    child: &[i128],
) -> Result<Vec<(i128, i128)>> {
    let mut result = BTreeMap::from([(0i128, 1i128)]);
    for equal in child.chunk_by(|a, b| a == b) {
        let value = equal[0];
        if value == 0 {
            break;
        }
        model.spend(1 + 2 * lambda.len().ilog2() as usize)?;
        let first = lambda.partition_point(|entry| *entry > value);
        let end = lambda.partition_point(|entry| *entry >= value);
        if equal.len() > end - first {
            let shift = -2 * equal.len() as i128;
            let mut next = result.clone();
            for (power, coefficient) in result {
                model.spend(1)?;
                add_scalar(
                    &mut next,
                    power + shift,
                    coefficient.checked_neg().ok_or_else(arithmetic)?,
                )?;
            }
            result = next;
        }
    }
    Ok(result.into_iter().collect())
}

pub(crate) fn coefficient<'a>(
    model: &mut Model,
    terms: &'a [Term],
    exponents: &[i128],
) -> Result<&'a [Term]> {
    model.spend(2 * (1 + exponents.len()) * (1 + terms.len().max(1).ilog2() as usize))?;
    let begin = terms.partition_point(|term| term.exponents.as_ref() < exponents);
    let end = begin + terms[begin..].partition_point(|term| term.exponents.as_ref() == exponents);
    Ok(&terms[begin..end])
}

fn gcd(mut a: u128, mut b: u128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a as i128
}

struct Fraction {
    numerator: i128,
    denominator: i128,
}

impl Fraction {
    fn integer(value: i128) -> Self {
        Self {
            numerator: value,
            denominator: 1,
        }
    }

    fn add(self, other: Self) -> Result<Self> {
        let common = gcd(self.denominator as u128, other.denominator as u128);
        let left = self.denominator / common;
        let right = other.denominator / common;
        let numerator = self
            .numerator
            .checked_mul(right)
            .and_then(|a| {
                other
                    .numerator
                    .checked_mul(left)
                    .and_then(|b| a.checked_add(b))
            })
            .ok_or_else(arithmetic)?;
        let reduce = gcd(numerator.unsigned_abs(), common as u128);
        let denominator = left
            .checked_mul(other.denominator / reduce)
            .ok_or_else(arithmetic)?;
        Ok(Self {
            numerator: numerator / reduce,
            denominator,
        })
    }

    fn multiply(self, other: Self) -> Result<Self> {
        let first = gcd(self.numerator.unsigned_abs(), other.denominator as u128);
        let second = gcd(other.numerator.unsigned_abs(), self.denominator as u128);
        let numerator = (self.numerator / first)
            .checked_mul(other.numerator / second)
            .ok_or_else(arithmetic)?;
        let denominator = (self.denominator / second)
            .checked_mul(other.denominator / first)
            .ok_or_else(arithmetic)?;
        Ok(Self {
            numerator,
            denominator,
        })
    }
}

fn evaluate_power(model: &mut Model, mut base: i128, power: i128) -> Result<Fraction> {
    let mut exponent = power.unsigned_abs();
    let mut value = 1i128;
    while exponent != 0 {
        model.spend(1)?;
        if exponent & 1 != 0 {
            value = value.checked_mul(base).ok_or_else(arithmetic)?;
        }
        exponent >>= 1;
        if exponent != 0 {
            base = base.checked_mul(base).ok_or_else(arithmetic)?;
        }
    }
    if power >= 0 {
        Ok(Fraction::integer(value))
    } else if value == 0 {
        Err(Error::plain(
            "cannot specialize a negative Laurent power at zero",
        ))
    } else if value < 0 {
        Ok(Fraction {
            numerator: -1,
            denominator: value.checked_neg().ok_or_else(arithmetic)?,
        })
    } else {
        Ok(Fraction {
            numerator: 1,
            denominator: value,
        })
    }
}

pub(crate) fn evaluate(model: &mut Model, terms: &[Term], v: i128, x: &[i128]) -> Result<String> {
    let mut sum = Fraction::integer(0);
    for term in terms {
        model.spend(1 + x.len())?;
        let mut value =
            Fraction::integer(term.coefficient).multiply(evaluate_power(model, v, term.power)?)?;
        for (&base, &power) in x.iter().zip(term.exponents.iter()) {
            value = value.multiply(evaluate_power(model, base, power)?)?;
        }
        sum = sum.add(value)?;
    }
    if sum.denominator == 1 {
        Ok(sum.numerator.to_string())
    } else {
        Ok(format!("{}/{}", sum.numerator, sum.denominator))
    }
}
