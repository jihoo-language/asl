use std::collections::BTreeMap;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;
use crate::rv_as_characters::{coefficient_error, Characters, Query as CharacterQuery};
use crate::rv_cyclotomic as cyclo;

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Reduced,
    Genus,
    PRank,
    Character(CharacterQuery),
}

pub(crate) struct Cover {
    prime: Prime,
    genus: usize,
    constant: u32,
    head: Side,
    tail: Side,
    characters: Characters,
}

#[derive(Default)]
struct Side {
    degree: usize,
    terms: Vec<(usize, u32, u32)>, // exponent, coefficient, logarithmic weight
}

impl Side {
    fn push(&mut self, prime: Prime, exponent: u128, coefficient: u32) -> Result<()> {
        let exponent = usize::try_from(exponent).map_err(|_| {
            Error::plain("two-pole Artin–Schreier reduced pole order is out of range")
        })?;
        self.degree = self.degree.max(exponent);
        self.terms.push((
            exponent,
            coefficient,
            prime.mul(prime.reduce(exponent), coefficient),
        ));
        Ok(())
    }

    fn phase(
        &self,
        model: &mut Model,
        prime: Prime,
        series: &[u32],
        reciprocal: &[u32],
    ) -> Result<u32> {
        let mut total = 0;
        for &(exponent, coefficient, _) in &self.terms {
            let width = exponent.min(series.len() - 1);
            model.spend(1 + width)?;
            // The power sum is −[t^j] tL'(t)/L(t).
            let sum = series[1..=width]
                .iter()
                .zip(reciprocal[exponent - width..exponent].iter().rev())
                .enumerate()
                .fold(0, |sum, (i, (&a, &b))| {
                    prime.add(sum, prime.mul(prime.reduce(i + 1), prime.mul(a, b)))
                });
            total = prime.add(total, prime.negate(prime.mul(coefficient, sum)));
        }
        Ok(total)
    }

    fn linear_phase(
        &self,
        model: &mut Model,
        prime: Prime,
        exponent: usize,
        reciprocal: &[u32],
    ) -> Result<u32> {
        let mut total = 0;
        for &(j, _, weight) in self.terms.iter().rev() {
            model.spend(1)?;
            if j < exponent {
                break;
            }
            total = prime.add(total, prime.mul(weight, reciprocal[j - exponent]));
        }
        // The common negative sign is immaterial for the orthogonality test.
        Ok(total)
    }
}

impl Cover {
    pub(crate) fn build(model: &mut Model, p: usize, terms: Vec<(i128, i128)>) -> Result<Self> {
        if p < 3 || p.is_multiple_of(2) {
            return Err(Error::plain("rvastwo requires an odd prime"));
        }
        cyclo::check_order(p)?;
        let prime = Prime::new(model, p)?;
        let mut normalized = BTreeMap::new();
        let modulus = i128::from(prime.get());
        for (mut exponent, coefficient) in terms {
            model.spend(1 + normalized.len().max(1).ilog2() as usize)?;
            while exponent != 0 && exponent % modulus == 0 {
                exponent /= modulus;
                model.spend(1)?;
            }
            let entry = normalized.entry(exponent).or_insert(0);
            *entry = prime.add(*entry, coefficient.rem_euclid(modulus) as u32);
        }
        let mut constant = 0;
        let (mut head, mut tail) = (Side::default(), Side::default());
        for (exponent, coefficient) in normalized {
            model.spend(1)?;
            if coefficient == 0 {
                continue;
            }
            match exponent.cmp(&0) {
                std::cmp::Ordering::Less => {
                    tail.push(prime, exponent.unsigned_abs(), coefficient)?
                }
                std::cmp::Ordering::Equal => constant = coefficient,
                std::cmp::Ordering::Greater => {
                    head.push(prime, exponent.unsigned_abs(), coefficient)?
                }
            }
        }
        tail.terms.reverse();
        if head.degree == 0 || tail.degree == 0 {
            return Err(Error::plain(
                "rvastwo requires nonconstant terms at both zero and infinity after reduction",
            ));
        }
        let degree = head
            .degree
            .checked_add(tail.degree)
            .ok_or_else(degree_error)?;
        let genus = ((p - 1) / 2).checked_mul(degree).ok_or_else(degree_error)?;
        model.allocate(8 + 3 * (head.terms.len() + tail.terms.len()))?;
        Ok(Self {
            prime,
            genus,
            constant,
            head,
            tail,
            characters: Characters::new(prime, degree, genus),
        })
    }

    pub(crate) fn order(&self) -> usize {
        self.prime.get() as usize
    }

    pub(crate) fn coefficient(&self, model: &mut Model, j: usize, k: usize) -> Result<Value> {
        self.characters
            .coefficient(model, j, k, |model| self.factor(model))
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, query: Query) -> Result<String> {
        match query {
            Query::Reduced => {
                model.spend(1 + self.head.terms.len() + self.tail.terms.len())?;
                let mut reduced = self
                    .tail
                    .terms
                    .iter()
                    .rev()
                    .map(|&(e, c, _)| (-(e as i128), c))
                    .collect::<Vec<_>>();
                if self.constant != 0 {
                    reduced.push((0, self.constant));
                }
                reduced.extend(self.head.terms.iter().map(|&(e, c, _)| (e as i128, c)));
                Ok(format!("rvASTwoReduced[{name}] = {reduced:?}"))
            }
            Query::Genus => Ok(format!("rvASTwoGenus[{name}] = {}", self.genus)),
            Query::PRank => Ok(format!("rvASTwoPrank[{name}] = {}", self.order() - 1)),
            Query::Character(query) => {
                let operation = match query {
                    CharacterQuery::Character(_) => "rvASTwoCharacterL",
                    CharacterQuery::CharacterNewton(_) => "rvASTwoCharacterNewton",
                    CharacterQuery::LPolynomial => "rvASTwoLPolynomial",
                    CharacterQuery::Newton => "rvASTwoNewton",
                    CharacterQuery::Points(_) => "rvASTwoPoints",
                };
                self.characters
                    .report(model, name, operation, query, |model| self.factor(model))
            }
        }
    }

    fn factor(&self, model: &mut Model) -> Result<Vec<Vec<i128>>> {
        let p = self.order();
        let degree = self.head.degree + self.tail.degree;
        if degree.saturating_add(1).saturating_mul(2 * p + 12) > MAX_PIVOT_CELLS {
            return Err(Error::plain(
                "two-pole Artin–Schreier stationary phase exceeds the storage limit",
            ));
        }
        let mut result = Vec::with_capacity(degree + 1);
        result.push(vec![1]);
        let mut head_reciprocal = vec![0; self.head.degree];
        let mut tail_reciprocal = vec![0; self.tail.degree];
        for k in 1..=degree {
            // Fixed head and tail positions are disjoint here, including when
            // the two nominal blocks overlap in a low polynomial degree.
            let head = (self.head.degree / 2).min(k - 1);
            let tail = k.saturating_sub(self.tail.degree / 2).max(head + 1);
            let fixed = head + (k - tail);
            let free = tail - head - 1;
            let fixed_positions = (1..=head).chain(tail..k);
            let choices = p
                .checked_pow(fixed as u32)
                .and_then(|choices| choices.checked_mul(p - 1).map(|total| (choices, total)))
                .ok_or_else(|| {
                    Error::plain(
                        "two-pole Artin–Schreier stationary phase exceeds the computation limit",
                    )
                })?;
            model.spend(choices.1)?;
            let weight = (p as i128)
                .checked_pow(free as u32)
                .ok_or_else(coefficient_error)?;
            let mut counts = vec![0i128; p];
            let mut forward = vec![0; k + 1];
            let mut backward = vec![0; k + 1];
            forward[0] = 1;
            backward[0] = 1;
            for last in 1..self.prime.get() {
                let inv_last = invert(model, self.prime, last)?;
                forward[k] = last;
                backward[k] = inv_last;
                for mut code in 0..choices.0 {
                    model.spend(fixed)?;
                    for i in fixed_positions.clone() {
                        let a = (code % p) as u32;
                        code /= p;
                        forward[i] = a;
                        backward[k - i] = self.prime.mul(a, inv_last);
                    }
                    reciprocal(model, self.prime, &forward, &mut head_reciprocal)?;
                    reciprocal(model, self.prime, &backward, &mut tail_reciprocal)?;
                    // B(t)=t^k A(1/t)/a_k encodes the reciprocal roots.
                    // Each free a_i changes the phase by a sum of two
                    // logarithmic derivatives. Quadratic terms exceed
                    // the relevant pole order on each side. The B-side
                    // coefficient changes by a_i/a_k.
                    let mut stationary = true;
                    for i in head + 1..tail {
                        let a = self
                            .head
                            .linear_phase(model, self.prime, i, &head_reciprocal)?;
                        let b =
                            self.tail
                                .linear_phase(model, self.prime, k - i, &tail_reciprocal)?;
                        if self.prime.add(a, self.prime.mul(b, inv_last)) != 0 {
                            stationary = false;
                            break;
                        }
                    }
                    if stationary {
                        let a = self
                            .head
                            .phase(model, self.prime, &forward, &head_reciprocal)?;
                        let b = self
                            .tail
                            .phase(model, self.prime, &backward, &tail_reciprocal)?;
                        let constant = self.prime.mul(self.prime.reduce(k), self.constant);
                        let phase = self.prime.add(self.prime.add(a, b), constant) as usize;
                        counts[phase] = counts[phase]
                            .checked_add(weight)
                            .ok_or_else(coefficient_error)?;
                    }
                }
            }
            result.push(cyclo::reduced(model, p, counts)?);
        }
        Ok(result)
    }
}

fn reciprocal(model: &mut Model, prime: Prime, series: &[u32], result: &mut [u32]) -> Result<()> {
    result[0] = 1;
    for j in 1..result.len() {
        let width = j.min(series.len() - 1);
        model.spend(1 + width)?;
        let sum = series[1..=width]
            .iter()
            .zip(result[j - width..j].iter().rev())
            .fold(0, |sum, (&a, &b)| prime.add(sum, prime.mul(a, b)));
        result[j] = prime.negate(sum);
    }
    Ok(())
}

fn degree_error() -> Error {
    Error::plain("two-pole Artin–Schreier genus or pole-order sum is out of range")
}
