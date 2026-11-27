use std::collections::BTreeMap;

use crate::error::Result;
use crate::model::Model;

use super::rational::{arithmetic, Integer, Rational, LIMBS};
use super::{storage, term_cells, work_limit, Form, Power, Tensor};

pub(super) fn integral(model: &mut Model, form: &Form, tensor: &Tensor) -> Result<Rational> {
    let mut moments = Moments::new(form);
    let mut sum = Rational::zero();
    for (powers, &coefficient) in &tensor.terms {
        let value = moments.monomial(model, powers)?;
        let product = coefficient.multiply(model, value)?;
        sum = sum.add(model, product)?;
    }
    Ok(sum)
}

pub(super) struct Moments<'a> {
    form: &'a Form,
    divisors: BTreeMap<Vec<Power>, Rational>,
    transcendental: BTreeMap<Vec<Power>, Integer>,
    traces: Vec<Integer>,
    factors: Vec<Option<Integer>>,
    cache_cells: usize,
}

impl<'a> Moments<'a> {
    pub(super) fn new(form: &'a Form) -> Self {
        Self {
            form,
            divisors: BTreeMap::new(),
            transcendental: BTreeMap::new(),
            traces: Vec::new(),
            factors: vec![Some(Integer::from_i128(form.multiplicity as i128))],
            cache_cells: 0,
        }
    }

    pub(super) fn monomial(&mut self, model: &mut Model, powers: &[Power]) -> Result<Rational> {
        let split =
            powers.partition_point(|power| usize::from(power.index) < self.form.gram.size());
        let (divisor, transcendental) = powers.split_at(split);
        let d = self.divisor(model, divisor)?;
        if d.is_zero() {
            return Ok(d);
        }
        let t = self.transcendental(model, transcendental)?;
        d.multiply(model, Rational::wide(t))
    }

    fn reserve(&mut self, powers: &[Power]) -> Result<()> {
        self.cache_cells = self.cache_cells.saturating_add(term_cells(powers));
        storage(self.cache_cells)
    }

    fn divisor(&mut self, model: &mut Model, powers: &[Power]) -> Result<Rational> {
        if let Some(&value) = self.divisors.get(powers) {
            return Ok(value);
        }
        if powers.is_empty() {
            return Ok(Rational::integer(1));
        }
        if !powers
            .iter()
            .filter(|power| !power.exponent.is_multiple_of(2))
            .count()
            .is_multiple_of(2)
        {
            return Ok(Rational::zero());
        }
        let axis = Counts::new(powers, 2 * LIMBS + 1)?;
        let mut moments = vec![Rational::zero(); axis.states];
        moments[0] = Rational::integer(1);
        let mut counts = vec![0; powers.len()];
        for state in 1..axis.states {
            model.spend(1 + powers.len())?;
            axis.decode(state, &mut counts);
            if !counts.iter().sum::<usize>().is_multiple_of(2) {
                continue;
            }
            let first = counts.iter().position(|&count| count != 0).unwrap();
            counts[first] -= 1;
            let rest = state - axis.steps[first];
            let mut sum = Rational::zero();
            for (j, &count) in counts.iter().enumerate().skip(first) {
                if count == 0 {
                    continue;
                }
                let pairing = *self
                    .form
                    .gram
                    .entry(powers[first].index.into(), powers[j].index.into());
                let term = pairing
                    .multiply(model, Rational::integer(count as i128))?
                    .multiply(model, moments[rest - axis.steps[j]])?;
                sum = sum.add(model, term)?;
            }
            moments[state] = sum;
        }
        let value = moments[axis.states - 1];
        self.reserve(powers)?;
        self.divisors.insert(powers.to_vec(), value);
        Ok(value)
    }

    fn transcendental(&mut self, model: &mut Model, powers: &[Power]) -> Result<Integer> {
        if let Some(&value) = self.transcendental.get(powers) {
            return Ok(value);
        }
        let rho = self.form.gram.size();
        let value = match powers {
            [] => Integer::from_i128(1),
            [power] if usize::from(power.index) == rho => self.radial(model, power.exponent)?,
            _ => self.mixed(model, powers)?,
        };
        self.reserve(powers)?;
        self.transcendental.insert(powers.to_vec(), value);
        Ok(value)
    }

    fn radial(&self, model: &mut Model, exponent: u128) -> Result<Integer> {
        let exponent = usize::try_from(exponent).map_err(|_| work_limit())?;
        model
            .spend(exponent.saturating_mul(LIMBS))
            .map_err(|_| work_limit())?;
        let rank = self.form.multiplicity * self.form.field.degree();
        let mut moment = Integer::from_i128(1);
        for j in 0..exponent {
            moment = moment
                .multiply(rank as i128 + 2 * j as i128)
                .map_err(|_| arithmetic())?;
        }
        Ok(moment)
    }

    fn mixed(&mut self, model: &mut Model, powers: &[Power]) -> Result<Integer> {
        let shape = Counts::new(powers, LIMBS + 1)?;
        let mut values = vec![Integer::from_i128(0); shape.states];
        values[0] = Integer::from_i128(1);
        let rho = self.form.gram.size();
        let exponents: Vec<_> = powers
            .iter()
            .map(|power| usize::from(power.index) - rho)
            .collect();
        let trace_cap = exponents
            .iter()
            .zip(&shape.maximum)
            .map(|(a, b)| a * b)
            .sum();
        self.form
            .field
            .extend_traces(model, &mut self.traces, trace_cap)?;

        let mut available = vec![0; powers.len()];
        let mut chosen = available.clone();
        for state in 1..shape.states {
            model.spend(1 + powers.len())?;
            shape.decode(state, &mut available);
            let first = available.iter().position(|&count| count != 0).unwrap();
            available[first] -= 1;
            let rest = state - shape.steps[first];
            chosen.fill(0);
            let mut result = Integer::from_i128(0);
            loop {
                model.spend(1 + LIMBS * powers.len())?;
                let block_id: usize = chosen
                    .iter()
                    .zip(&shape.steps)
                    .map(|(digit, step)| digit * step)
                    .sum();
                let tail = values[rest - block_id];
                if !tail.is_zero() {
                    let trace = exponents[first]
                        + chosen
                            .iter()
                            .zip(&exponents)
                            .map(|(count, power)| count * power)
                            .sum::<usize>();
                    if !self.traces[trace].is_zero() {
                        let size = 1 + chosen.iter().sum::<usize>();
                        let log = self
                            .cumulant_factor(size)?
                            .multiply_wide(self.traces[trace])
                            .map_err(|_| arithmetic())?;
                        let mut value = tail.multiply_wide(log).map_err(|_| arithmetic())?;
                        // Attach a block containing the first labeled factor.
                        // Counts recover exactly the number of indistinguishable
                        // choices for its other labels.
                        for (&count, &take) in available.iter().zip(&chosen) {
                            value = value
                                .multiply_wide(choose(model, count, take)?)
                                .map_err(|_| arithmetic())?;
                        }
                        result = result.add(value).map_err(|_| arithmetic())?;
                    }
                }
                if !advance(&mut chosen, &available) {
                    break;
                }
            }
            values[state] = result;
        }
        Ok(values[shape.states - 1])
    }

    fn cumulant_factor(&mut self, size: usize) -> Result<Integer> {
        storage(size.saturating_mul(LIMBS + 2))?;
        while self.factors.len() < size {
            let factor = 2 * self.factors.len() as i128;
            let term = self
                .factors
                .last()
                .unwrap()
                .and_then(|value| value.multiply(factor).ok());
            self.factors.push(term);
        }
        self.factors[size - 1].ok_or_else(arithmetic)
    }
}

fn choose(model: &mut Model, n: usize, k: usize) -> Result<Integer> {
    let mut value = Integer::from_i128(1);
    for i in 1..=k.min(n - k) {
        value = value
            .multiply((n + 1 - i) as i128)
            .map_err(|_| arithmetic())?;
        let (quotient, remainder) = value.magnitude.divide(model, i as u128)?;
        if remainder != 0 {
            return Err(arithmetic());
        }
        value.magnitude = quotient;
    }
    Ok(value)
}

fn advance(digits: &mut [usize], maximum: &[usize]) -> bool {
    for (digit, &maximum) in digits.iter_mut().zip(maximum) {
        if *digit < maximum {
            *digit += 1;
            return true;
        }
        *digit = 0;
    }
    false
}

struct Counts {
    maximum: Vec<usize>,
    steps: Vec<usize>,
    states: usize,
}

impl Counts {
    fn new(powers: &[Power], width: usize) -> Result<Self> {
        let mut result = Self {
            maximum: Vec::new(),
            steps: Vec::new(),
            states: 1,
        };
        for power in powers {
            let exponent = usize::try_from(power.exponent).map_err(|_| work_limit())?;
            result.maximum.push(exponent);
            result.steps.push(result.states);
            result.states = result.states.saturating_mul(exponent.saturating_add(1));
        }
        storage(result.states.saturating_mul(width))?;
        Ok(result)
    }

    fn decode(&self, state: usize, digits: &mut [usize]) {
        for (i, digit) in digits.iter_mut().enumerate() {
            *digit = state / self.steps[i] % (self.maximum[i] + 1);
        }
    }
}
