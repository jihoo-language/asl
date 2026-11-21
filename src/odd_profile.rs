use std::ops::Range;

use crate::error::Result;
use crate::model::Model;
use crate::odd_cobar::{self, Factor, Prime, Term, Word};
use crate::profile::Profile;

pub(crate) fn changes_basis(prime: Prime, n: usize, max_factor: usize) -> bool {
    // The lowest forbidden positive degree is that of xi_1^(p^n).
    odd_cobar::xi_degree(prime, 1)
        .ok()
        .zip(odd_cobar::power(prime, n).ok())
        .and_then(|(degree, exponent)| degree.checked_mul(exponent))
        .is_some_and(|cutoff| max_factor >= cutoff)
}

pub(crate) fn allows_power(prime: Prime, profile: Profile, index: usize, exponent: usize) -> bool {
    match profile {
        Profile::Full => true,
        Profile::An(n) => exponent <= exponent_bound(prime, n, index),
    }
}

fn exponent_bound(prime: Prime, n: usize, index: usize) -> usize {
    let Some(exponent) = n.checked_sub(index) else {
        return 0;
    };
    if exponent >= usize::BITS as usize {
        return usize::MAX;
    }
    (prime.get() as usize)
        .checked_pow(exponent as u32)
        .map_or(usize::MAX, |power| power - 1)
}

pub(crate) fn restrict(
    model: &mut Model,
    prime: Prime,
    n: usize,
    terms: &[Term],
) -> Result<Vec<Term>> {
    let mut projected = Vec::new();
    for term in terms {
        model.spend(odd_cobar::term_cells(term))?;
        if term.0.iter().all(|factor| {
            factor.taus.last().is_none_or(|&index| index <= n)
                && factor
                    .exponents
                    .iter()
                    .enumerate()
                    .all(|(index, &exponent)| allows_power(prime, Profile::An(n), index, exponent))
        }) {
            projected.push(term.clone());
        }
    }
    Ok(projected)
}

pub(crate) fn basis(
    model: &mut Model,
    prime: Prime,
    n: usize,
    s: usize,
    t: usize,
    limit: usize,
) -> Result<Vec<Word>> {
    model.spend(1)?;
    let maximum = t - s + 1;
    let (variables, top_degree) = variables(prime, n, maximum);
    let maximum = maximum.min(top_degree);
    if t > s.saturating_mul(maximum) {
        return Ok(Vec::new());
    }
    let tau_count = variables
        .iter()
        .filter(|var| matches!(var.kind, Kind::Tau(_)))
        .count();
    if let Ok(base) = odd_cobar::xi_degree(prime, 1) {
        // All polynomial weights are multiples of the base, while each
        // exterior generator contributes one to this degree residue.
        if t % base > s.saturating_mul(tau_count) {
            return Ok(Vec::new());
        }
    }
    odd_cobar::check_storage(s.saturating_mul(2).saturating_add(1))?;
    let minimum = t.saturating_sub((s - 1).saturating_mul(maximum)).max(1);
    let mut catalogue = Catalogue {
        variables: &variables,
        maximum,
        exponents: vec![0; variables.len()],
        taus: Vec::new(),
        factors: Vec::new(),
        cells: 0,
    };
    catalogue.enumerate(model, variables.len(), minimum, maximum)?;
    let factor_count = catalogue.factors.len();
    if factor_count == 0 {
        return Ok(Vec::new());
    }
    model.spend(catalogue.cells * (1 + factor_count.ilog2() as usize))?;
    catalogue
        .factors
        .sort_unstable_by_key(|factor| factor.degree);
    words(model, s, t, limit, catalogue)
}

#[derive(Clone, Copy)]
enum Kind {
    Xi(usize),
    Tau(usize),
}

#[derive(Clone, Copy)]
struct Variable {
    kind: Kind,
    weight: usize,
    bound: usize,
    lower_capacity: usize,
}

fn variables(prime: Prime, n: usize, maximum: usize) -> (Vec<Variable>, usize) {
    let mut variables = vec![Variable {
        kind: Kind::Tau(0),
        weight: 1,
        bound: 1,
        lower_capacity: 0,
    }];
    let mut capacity = 1usize;
    for i in 1..=n {
        let Ok(weight) = odd_cobar::xi_degree(prime, i) else {
            break;
        };
        if weight > maximum {
            break;
        }
        let bound = exponent_bound(prime, n, i - 1);
        variables.push(Variable {
            kind: Kind::Xi(i - 1),
            weight,
            bound,
            lower_capacity: capacity,
        });
        capacity = capacity.saturating_add(weight.saturating_mul(bound));
        if weight < maximum {
            let weight = weight + 1;
            variables.push(Variable {
                kind: Kind::Tau(i),
                weight,
                bound: 1,
                lower_capacity: capacity,
            });
            capacity = capacity.saturating_add(weight);
        }
    }
    (variables, capacity)
}

struct Catalogue<'a> {
    variables: &'a [Variable],
    maximum: usize,
    exponents: Vec<usize>,
    taus: Vec<usize>,
    factors: Vec<Factor>,
    cells: usize,
}

impl Catalogue<'_> {
    fn enumerate(
        &mut self,
        model: &mut Model,
        count: usize,
        minimum: usize,
        maximum: usize,
    ) -> Result<()> {
        model.spend(1)?;
        let Some(index) = count.checked_sub(1) else {
            debug_assert_eq!(minimum, 0);
            model.spend(self.exponents.len() + self.taus.len())?;
            let degree = self.maximum - maximum;
            let taus = self.taus.iter().rev().copied().collect();
            let factor = odd_cobar::factor(self.exponents.clone(), taus, degree);
            self.cells += factor.cells();
            odd_cobar::check_storage(self.cells)?;
            self.factors.push(factor);
            return Ok(());
        };
        let Variable {
            kind,
            weight,
            bound,
            lower_capacity,
        } = self.variables[index];
        let first = minimum.saturating_sub(lower_capacity).div_ceil(weight);
        let last = (maximum / weight).min(bound);
        for exponent in first..=last {
            let spent = exponent * weight;
            let (low, high) = (minimum.saturating_sub(spent), maximum - spent);
            match kind {
                Kind::Xi(index) => {
                    self.exponents[index] = exponent;
                    self.enumerate(model, count - 1, low, high)?;
                }
                Kind::Tau(index) => {
                    if exponent != 0 {
                        self.taus.push(index);
                    }
                    self.enumerate(model, count - 1, low, high)?;
                    if exponent != 0 {
                        self.taus.pop();
                    }
                }
            }
        }
        Ok(())
    }
}

struct Frame {
    remaining: usize,
    choices: Range<usize>,
}

impl Frame {
    fn new(model: &mut Model, factors: &[Factor], remaining: usize, slots: usize) -> Result<Self> {
        let later = slots - 1;
        let minimum = factors.first().unwrap().degree;
        let maximum = factors.last().unwrap().degree;
        let choices = if let Some(last) = later
            .checked_mul(minimum)
            .and_then(|rest| remaining.checked_sub(rest))
        {
            let first = remaining.saturating_sub(later.saturating_mul(maximum));
            model.spend(2 * (1 + factors.len().ilog2() as usize))?;
            let start = factors.partition_point(|factor| factor.degree < first);
            let end = factors.partition_point(|factor| factor.degree <= last);
            start..end
        } else {
            0..0
        };
        Ok(Self { remaining, choices })
    }
}

fn words(
    model: &mut Model,
    s: usize,
    t: usize,
    limit: usize,
    catalogue: Catalogue<'_>,
) -> Result<Vec<Word>> {
    let factors = catalogue.factors;
    let mut words = Vec::new();
    let mut cells = 0;
    let mut partial = Vec::new();
    let mut stack = vec![Frame::new(model, &factors, t, s)?];
    while let Some(frame) = stack.last_mut() {
        model.spend(1)?;
        if let Some(index) = frame.choices.next() {
            let factor = &factors[index];
            let remaining = frame.remaining - factor.degree;
            if partial.len() + 1 == s {
                if words.len() == limit {
                    return Err(odd_cobar::basis_limit(limit));
                }
                let word: Word = partial
                    .iter()
                    .chain(std::iter::once(factor))
                    .cloned()
                    .collect();
                let size = odd_cobar::word_cells(&word);
                cells += size;
                odd_cobar::check_storage(cells + catalogue.cells)?;
                model.spend(size)?;
                words.push(word);
            } else {
                partial.push(factor.clone());
                stack.push(Frame::new(model, &factors, remaining, s - partial.len())?);
            }
        } else {
            stack.pop();
            partial.pop();
        }
    }
    model.spend(cells * (1 + words.len().max(1).ilog2() as usize))?;
    words.sort_unstable();
    Ok(words)
}
