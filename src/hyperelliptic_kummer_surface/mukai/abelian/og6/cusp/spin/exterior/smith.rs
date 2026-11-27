use std::collections::BTreeMap;
use std::fmt::Write;

use super::super::{Algebra, EVEN_RANK};
use super::Form;
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::Limit;
use crate::hyperelliptic_kummer_surface::rational_sparse::{decimal, Magnitude, WORK_LIMBS};
use crate::model::Model;

const PRIMES: [u128; 11] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31];
type Factors = Vec<(Magnitude, u128)>;
type Diagonal = BTreeMap<u128, u128>;

pub(super) struct Spectrum {
    pub(super) dimension: u128,
    pub(super) base: Magnitude,
    pub(super) power: usize,
    factors: Factors,
}

impl Spectrum {
    pub(super) fn zero(dimension: u128) -> Self {
        Self {
            dimension,
            base: Magnitude::from_small(1),
            power: 0,
            factors: Vec::new(),
        }
    }

    pub(super) fn rank(&self) -> u128 {
        self.factors.iter().map(|(_, count)| count).sum()
    }

    pub(super) fn free(&self) -> u128 {
        self.dimension - self.rank()
    }

    pub(super) fn factor_text(&self, omit_unit: bool) -> String {
        let mut entries = String::new();
        for &(factor, count) in &self.factors {
            if omit_unit && factor.small() == Some(1) {
                continue;
            }
            if !entries.is_empty() {
                entries.push_str(", ");
            }
            write!(entries, "({}, {count})", decimal(factor)).unwrap();
        }
        entries
    }
}

impl Algebra {
    pub(crate) fn exterior_smith_report(
        &self,
        model: &mut Model,
        label: &str,
        limit: &Limit,
        degree: usize,
        order: usize,
    ) -> Result<String> {
        let spectrum = self.exterior_smith_spectrum(model, limit, degree, order)?;
        let dimension = spectrum.dimension;
        let rank = spectrum.rank();
        let free = spectrum.free();
        let power = spectrum.power;
        let entries = spectrum.factor_text(false);
        Ok(format!(
            "rvKSExteriorSmith[{label},{degree},{order}] = dimension={dimension}, rank={rank}, kernel={free}, cokernel_free={free}, common_base={}, common_power={power}, factors=[{entries}]",
            decimal(spectrum.base),
        ))
    }

    pub(super) fn exterior_smith_spectrum(
        &self,
        model: &mut Model,
        limit: &Limit,
        degree: usize,
        order: usize,
    ) -> Result<Spectrum> {
        if degree > EVEN_RANK {
            return Err(Error::plain(
                "Kuga–Satake cohomological degree must be between 0 and 64",
            ));
        }
        let degree_bound = degree.min(EVEN_RANK - degree);
        let dimension = exterior_dimension(degree);
        let one = Magnitude::from_small(1);
        let (base, power, factors) = if order == 0 {
            (one, 0, vec![(one, dimension)])
        } else if order > degree_bound {
            return Ok(Spectrum::zero(dimension));
        } else {
            let spin = self.action(model, limit)?.smith(model)?;
            if order > spin.len() {
                return Ok(Spectrum::zero(dimension));
            } else if spin.iter().all(|d| *d == spin[0]) {
                model.spend(EVEN_RANK * (EVEN_RANK + 1) / 2)?;
                let factors = uniform(model, &binomial(), degree_bound, order, spin.len())?;
                (spin[0], order, factors)
            } else if degree_bound <= 2 {
                (
                    spin[0],
                    order,
                    nonuniform(model, spin, degree_bound, order)?,
                )
            } else {
                return Err(Error::plain(
                    "nonuniform Kuga–Satake spin Smith factors support exterior Smith degrees 0, 1, 2, 62, 63, and 64, as well as identity or zero powers",
                ));
            }
        };
        Ok(Spectrum {
            dimension,
            base,
            power,
            factors,
        })
    }
}

pub(super) fn exterior_dimension(degree: usize) -> u128 {
    if degree > EVEN_RANK {
        return 0;
    }
    let degree = degree.min(EVEN_RANK - degree);
    (0..degree).fold(1, |value, i| {
        value * (EVEN_RANK - i) as u128 / (i + 1) as u128
    })
}

type Binomial = [[u128; EVEN_RANK + 1]; EVEN_RANK + 1];

fn binomial() -> Binomial {
    let mut table = [[0; EVEN_RANK + 1]; EVEN_RANK + 1];
    table[0][0] = 1;
    for n in 1..=EVEN_RANK {
        let (previous, rows) = table.split_at_mut(n);
        let row = &mut rows[0];
        row[0] = 1;
        for (k, entry) in row.iter_mut().enumerate().take(n + 1).skip(1) {
            *entry = previous[n - 1][k - 1] + previous[n - 1][k];
        }
    }
    table
}

fn uniform(
    model: &mut Model,
    binomial: &Binomial,
    degree: usize,
    order: usize,
    rank: usize,
) -> Result<Factors> {
    let diagonal = uniform_diagonal(model, binomial, degree, order, rank)?;
    canonical_valuations(model, &diagonal)
}

fn uniform_diagonal(
    model: &mut Model,
    binomial: &Binomial,
    degree: usize,
    order: usize,
    rank: usize,
) -> Result<Diagonal> {
    let inert = EVEN_RANK - 2 * rank;
    let mut diagonal = BTreeMap::<u128, u128>::new();
    for singles in order..=rank.min(degree) {
        let space = degree - singles;
        let placements: u128 = (0..=(rank - singles).min(space / 2))
            .map(|doubles| binomial[rank - singles][doubles] * binomial[inert][space - 2 * doubles])
            .sum();
        let placements = binomial[rank][singles] * placements;
        if placements == 0 {
            continue;
        }
        for source in order..=singles {
            // Transpose after taking complements when the target layer is larger.
            let smaller = (source - order).min(singles - source);
            model.spend(smaller + 1)?;
            let mut previous = 0;
            // Wilson's diagonal precedes sorting valuations across all layers.
            for (v, &subsets) in binomial[singles].iter().enumerate().take(smaller + 1) {
                let value = binomial[smaller + order - v][order];
                *diagonal.entry(value).or_default() += placements * (subsets - previous);
                previous = subsets;
            }
        }
    }
    Ok(diagonal)
}

pub(super) fn exterior_group_ranks_modulo(
    model: &mut Model,
    degree: usize,
    active_rank: usize,
    prime: u32,
) -> Result<Vec<u128>> {
    let mut ranks = vec![exterior_dimension(degree)];
    let degree = degree.min(EVEN_RANK - degree);
    let bound = degree.min(active_rank).min(prime as usize - 1);
    if bound == 0 {
        return Ok(ranks);
    }
    model.spend(EVEN_RANK * (EVEN_RANK + 1) / 2)?;
    let binomial = binomial();
    for order in 1..=bound {
        let diagonal = uniform_diagonal(model, &binomial, degree, order, active_rank)?;
        let rank = diagonal
            .into_iter()
            .filter(|&(value, _)| !value.is_multiple_of(u128::from(prime)))
            .map(|(_, count)| count)
            .sum();
        ranks.push(rank);
    }
    Ok(ranks)
}

fn canonical_valuations(model: &mut Model, diagonal: &BTreeMap<u128, u128>) -> Result<Factors> {
    let rank: u128 = diagonal.values().sum();
    let mut changes = BTreeMap::from([(0, 1u128), (rank, 1)]);
    for p in PRIMES {
        let mut valuations = BTreeMap::<u32, u128>::new();
        model.spend(diagonal.len())?;
        for (&value, &count) in diagonal {
            let mut value = value;
            let mut exponent = 0;
            while value.is_multiple_of(p) {
                value /= p;
                exponent += 1;
            }
            *valuations.entry(exponent).or_default() += count;
        }
        let mut index = 0;
        let mut previous = 0;
        for (exponent, count) in valuations {
            if exponent > previous {
                let multiplier = changes.entry(index).or_insert(1);
                *multiplier = multiplier
                    .checked_mul(p.pow(exponent - previous))
                    .ok_or_else(Form::arithmetic)?;
            }
            index += count;
            previous = exponent;
        }
    }
    let mut changes = changes.into_iter().peekable();
    let mut factor = 1u128;
    let mut result = Vec::new();
    while let Some((index, multiplier)) = changes.next() {
        if let Some(&(next, _)) = changes.peek() {
            factor = factor
                .checked_mul(multiplier)
                .ok_or_else(Form::arithmetic)?;
            result.push((Magnitude::from_small(factor), next - index));
        }
    }
    Ok(result)
}

fn nonuniform(
    model: &mut Model,
    spin: &[Magnitude],
    degree: usize,
    order: usize,
) -> Result<Factors> {
    let ratios = spin
        .iter()
        .map(|value| {
            value
                .divide_exact(model, spin[0])?
                .ok_or_else(|| Error::plain("spin Smith factors do not form a divisibility chain"))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut factors = Vec::new();
    if order == 1 {
        for (i, ratio) in ratios.into_iter().enumerate() {
            // A singly occupied pair with an inert vector gives one entry.
            // Two singly occupied pairs give two gcds, both equal to the earlier factor.
            let count = if degree == 1 {
                1
            } else {
                (EVEN_RANK - 2 - 2 * i) as u128
            };
            push_factor(&mut factors, ratio, count);
        }
    } else {
        let mut products = Vec::new();
        for (i, &d) in ratios.iter().enumerate() {
            for &e in &ratios[i + 1..] {
                products.push(Some(d.checked_product(e).ok_or_else(Form::arithmetic)?));
            }
        }
        divisibility_sort(model, &mut products)?;
        for factor in products.into_iter().flatten() {
            push_factor(&mut factors, factor, 1);
        }
    }
    Ok(factors)
}

fn push_factor(factors: &mut Factors, value: Magnitude, count: u128) {
    if count == 0 {
        return;
    }
    if let Some((prior, multiplicity)) = factors.last_mut() {
        if *prior == value {
            *multiplicity += count;
            return;
        }
    }
    factors.push((value, count));
}

// A bitonic sorting network sorts every prime valuation independently. Its
// comparison is gcd/lcm, and padding is an infinite valuation at every prime.
fn divisibility_sort(model: &mut Model, values: &mut Vec<Option<Magnitude>>) -> Result<()> {
    let length = values.len().next_power_of_two();
    values.resize(length, None);
    let mut width = 2;
    while width <= length {
        let mut stride = width / 2;
        while stride > 0 {
            model.spend(length * WORK_LIMBS)?;
            for i in 0..length {
                let j = i ^ stride;
                if j <= i || values[i] == values[j] {
                    continue;
                }
                let pair = match (values[i], values[j]) {
                    (Some(a), Some(b)) => {
                        let gcd = a.gcd(model, b)?;
                        let lcm = if gcd == a {
                            b
                        } else if gcd == b {
                            a
                        } else {
                            a.divide_exact(model, gcd)?
                                .expect("gcd divides its argument")
                                .checked_product(b)
                                .ok_or_else(Form::arithmetic)?
                        };
                        (Some(gcd), Some(lcm))
                    }
                    (a, b) => (a.or(b), None),
                };
                (values[i], values[j]) = if i & width == 0 {
                    pair
                } else {
                    (pair.1, pair.0)
                };
            }
            stride /= 2;
        }
        width *= 2;
    }
    Ok(())
}
