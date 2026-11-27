use crate::error::{Error, Result};
use crate::model::Model;

use super::rational::{arithmetic, Integer, Rational, LIMBS};
use super::storage;

pub(super) struct Etale {
    polynomial: Vec<i128>,
}

impl Etale {
    pub(super) fn new(model: &mut Model, polynomial: Vec<i128>) -> Result<Self> {
        // The shape and monicity are checked together with total dimension.
        let p: Vec<_> = polynomial.iter().copied().map(Rational::integer).collect();
        let q = p
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, &coefficient)| coefficient.multiply(model, Rational::integer(i as i128)))
            .collect::<Result<_>>()?;
        certify_totally_real(model, p, q)?;
        Ok(Self { polynomial })
    }

    pub(super) fn degree(&self) -> usize {
        self.polynomial.len() - 1
    }

    pub(super) fn extend_traces(
        &self,
        model: &mut Model,
        powers: &mut Vec<Integer>,
        cap: usize,
    ) -> Result<()> {
        storage(cap.saturating_add(1).saturating_mul(LIMBS + 1))?;
        let d = self.degree();
        if powers.is_empty() {
            powers.push(Integer::from_i128(d as i128));
        }
        while powers.len() <= cap {
            let k = powers.len();
            let mut value = Integer::from_i128(0);
            let last = d.min(k - 1);
            model.spend((last + 1) * LIMBS)?;
            for i in 1..=last {
                let term = powers[k - i]
                    .multiply(self.polynomial[d - i])
                    .map_err(|_| arithmetic())?;
                value = value.add(term.negate()).map_err(|_| arithmetic())?;
            }
            if k <= d {
                let term = Integer::product(k as i128, self.polynomial[d - k]);
                value = value.add(term.negate()).map_err(|_| arithmetic())?;
            }
            powers.push(value);
        }
        Ok(())
    }

    pub(super) fn element(
        &self,
        model: &mut Model,
        entries: &[(i128, i128)],
    ) -> Result<Vec<Rational>> {
        if entries.len() > self.degree() {
            return Err(Error::plain(
                "Hodge Wick element has too many power-basis coordinates",
            ));
        }
        entries
            .iter()
            .map(|&entry| Rational::input(model, entry))
            .collect()
    }

    pub(super) fn trace(&self, model: &mut Model, entries: &[(i128, i128)]) -> Result<Rational> {
        let element = self.element(model, entries)?;
        let mut powers = Vec::new();
        self.extend_traces(model, &mut powers, element.len().saturating_sub(1))?;
        let mut trace = Rational::zero();
        for (value, trace_power) in element.into_iter().zip(powers) {
            let contribution = value.multiply(model, Rational::wide(trace_power))?;
            trace = trace.add(model, contribution)?;
        }
        Ok(trace)
    }
}

fn certify_totally_real(
    model: &mut Model,
    mut p: Vec<Rational>,
    mut q: Vec<Rational>,
) -> Result<()> {
    let degree = p.len() - 1;
    let mut plus = Variations::default();
    let mut minus = Variations::default();
    append_signs(&p, &mut plus, &mut minus);
    while !q.is_empty() {
        append_signs(&q, &mut plus, &mut minus);
        let mut remainder = negative_remainder(model, p, &q)?;
        if remainder.is_empty() {
            if q.len() != 1 {
                return Err(Error::plain(
                    "Hodge Wick étale polynomial must be square-free",
                ));
            }
        } else {
            // A positive rescaling preserves both Sturm variation counts.
            let mut scale = remainder.last().unwrap().inverse()?;
            if scale.negative() {
                scale = scale.negate();
            }
            for coefficient in &mut remainder {
                *coefficient = coefficient.multiply(model, scale)?;
            }
        }
        p = q;
        q = remainder;
    }
    if minus.count.checked_sub(plus.count) != Some(degree) {
        return Err(Error::plain(
            "Hodge Wick étale polynomial must have all roots real",
        ));
    }
    Ok(())
}

fn negative_remainder(
    model: &mut Model,
    mut p: Vec<Rational>,
    q: &[Rational],
) -> Result<Vec<Rational>> {
    let inverse = q.last().unwrap().inverse()?;
    while p.len() >= q.len() {
        let shift = p.len() - q.len();
        let scalar = p.pop().unwrap().multiply(model, inverse)?;
        for (out, &source) in p[shift..].iter_mut().zip(q) {
            let subtraction = source.multiply(model, scalar)?.negate();
            *out = out.add(model, subtraction)?;
        }
        while p.last().is_some_and(|x| x.is_zero()) {
            p.pop();
        }
    }
    Ok(p.into_iter().map(Rational::negate).collect())
}

#[derive(Default)]
struct Variations {
    count: usize,
    last: Option<bool>,
}

impl Variations {
    fn sign(&mut self, negative: bool) {
        self.count += usize::from(self.last.is_some_and(|last| last != negative));
        self.last = Some(negative);
    }
}

fn append_signs(p: &[Rational], plus: &mut Variations, minus: &mut Variations) {
    let negative = p.last().unwrap().negative();
    plus.sign(negative);
    minus.sign(negative != p.len().is_multiple_of(2));
}
