use crate::error::{Error, Result};
use crate::model::Model;
use crate::wide_integer::Natural;

use super::necklace::factorization;
use super::norm::{cycles, half, power};
use super::{
    add, degree_limit, determinant_from_ghosts, integrity, multiply, power_sums, require_degree,
    trace_at_product, zero, Integer, Spectrum, MAX_DEGREE,
};

const MAX_DIVISORS: usize = 65_536;

#[derive(Clone, Copy)]
pub(crate) enum CyclicPower {
    Isotypic(u128, u128),
    Lie(u128),
}

impl CyclicPower {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Isotypic(..) => "rvWittSpectrumCyclic",
            Self::Lie(_) => "rvWittSpectrumLie",
        }
    }

    pub(super) fn ghost_label(self, k: u128) -> String {
        match self {
            Self::Isotypic(n, d) => format!("rvWittSpectrumCyclicGhost[{n}, {d}, {k}]"),
            Self::Lie(n) => format!("rvWittSpectrumLieGhost[{n}, {k}]"),
        }
    }

    fn length(self) -> Result<u128> {
        let n = match self {
            Self::Isotypic(n, _) | Self::Lie(n) => n,
        };
        if n == 0 {
            return Err(Error::plain("cyclic tensor length must be positive"));
        }
        if matches!(self, Self::Isotypic(_, d) if d == 0 || !n.is_multiple_of(d)) {
            return Err(Error::plain(
                "cyclotomic order must be a positive divisor of the tensor length",
            ));
        }
        Ok(n)
    }

    fn line_supported(self, odd: bool) -> bool {
        match self {
            Self::Isotypic(n, d) => d == if odd && n.is_multiple_of(2) { 2 } else { 1 },
            Self::Lie(n) => n == 1 || odd && n == 2,
        }
    }
}

impl Spectrum {
    pub(crate) fn cyclic_power(&self, model: &mut Model, projection: CyclicPower) -> Result<Self> {
        let n = projection.length()?;
        if n == 1 {
            return Ok(self.clone());
        }
        if let Some(line) = self.cyclic_line(projection) {
            let Some((eigenvalue, odd)) = line else {
                return Ok(Self::zero());
            };
            let eigenvalue = power(model, eigenvalue, n)?;
            let spectrum = Self::effective(vec![Integer::from_i128(1), eigenvalue.negate()])?;
            return Ok(if odd && !n.is_multiple_of(2) {
                spectrum.negate()
            } else {
                spectrum
            });
        }
        check_primitive_lower_bound(n)?;
        let projector = Projector::new(model, projection, n)?;
        let ranks = projector.parity_traces(model, |_| {
            (dimension(&self.numerator), dimension(&self.denominator))
        })?;
        let even = require_rank(ranks[0])?;
        let odd = require_rank(ranks[1])?;
        let rank = even.max(odd);
        let limit = projector.last_divisor() as usize * rank;
        let pos = power_sums(model, &self.numerator, limit)?;
        let neg = power_sums(model, &self.denominator, limit)?;
        model.spend(even + odd)?;
        let mut even_ghosts = Vec::with_capacity(even);
        let mut odd_ghosts = Vec::with_capacity(odd);
        for k in 1..=rank {
            let [a, b] = projector.parity_traces(model, |h| {
                let i = h as usize * k;
                (pos[i], neg[i])
            })?;
            if k <= even {
                even_ghosts.push(a);
            }
            if k <= odd {
                odd_ghosts.push(b);
            }
        }
        let p = determinant_from_ghosts(model, &even_ghosts)?;
        let q = determinant_from_ghosts(model, &odd_ghosts)?;
        if p[even].is_zero() || q[odd].is_zero() {
            return Err(integrity());
        }
        Self::canonical(model, p, q)
    }

    pub(super) fn cyclic_ghost(
        &self,
        model: &mut Model,
        projection: CyclicPower,
        k: u128,
    ) -> Result<Integer> {
        let n = projection.length()?;
        if k == 0 {
            return Err(Error::plain("cyclic tensor ghost index must be positive"));
        }
        if n == 1 {
            return self.ghost(model, k);
        }
        if let Some(line) = self.cyclic_line(projection) {
            let Some((eigenvalue, odd)) = line else {
                return Ok(zero());
            };
            let initial = power(model, eigenvalue, n)?;
            let value = power(model, initial, k)?;
            return Ok(if odd && !n.is_multiple_of(2) {
                value.negate()
            } else {
                value
            });
        }
        // The h=1 term is always present. Check its arithmetic bound before
        // trying to factor a potentially large tensor length.
        let ghost = self.ghost(model, k)?;
        let initial = power(model, ghost, n)?;
        let projector = Projector::new(model, projection, n)?;
        projector.trace(model, |model, h| {
            if h == 1 {
                return Ok(initial);
            }
            let a = trace_at_product(model, &self.numerator, k, h)?;
            let b = trace_at_product(model, &self.denominator, k, h)?;
            power(model, add(a, b.negate())?, n / h)
        })
    }

    fn cyclic_line(&self, projection: CyclicPower) -> Option<Option<(Integer, bool)>> {
        let (coefficient, odd) = match (&*self.numerator, &*self.denominator) {
            ([1], [1]) => return Some(None),
            ([1, coefficient], [1]) => (coefficient, false),
            ([1], [1, coefficient]) => (coefficient, true),
            _ => return None,
        };
        Some(
            projection
                .line_supported(odd)
                .then(|| (Integer::from_i128(*coefficient).negate(), odd)),
        )
    }
}

struct Projector {
    length: u128,
    weights: Vec<(u128, Integer)>,
}

impl Projector {
    fn new(model: &mut Model, projection: CyclicPower, length: u128) -> Result<Self> {
        let factors = factorization(
            model,
            length,
            "cyclic tensor factorization exceeds the trial division limit",
        )?;
        let lie = matches!(projection, CyclicPower::Lie(_));
        let all = divisors(model, &factors, lie)?;
        let mut weights = Vec::with_capacity(all.len());
        let order_phi = match projection {
            CyclicPower::Isotypic(_, order) => totient(order, &factors),
            CyclicPower::Lie(_) => 1,
        };
        for h in all {
            model.spend(1 + 3 * factors.len())?;
            let (argument, degree) = match projection {
                CyclicPower::Isotypic(_, order) => {
                    let (g, _) = cycles(order, length / h)?;
                    (order / g, totient(h, &factors))
                }
                CyclicPower::Lie(_) => (h, 1),
            };
            if let Some(negative) = mobius_sign(argument, &factors) {
                // The primitive projector has φ(n) times the Lie character.
                // Its normalized coefficients are μ(h), with no φ factors.
                let quotient = if lie {
                    1
                } else {
                    order_phi / totient(argument, &factors)
                };
                let weight = multiply(model, natural(degree), natural(quotient))?;
                weights.push((h, if negative { weight.negate() } else { weight }));
            }
        }
        Ok(Self { length, weights })
    }

    fn last_divisor(&self) -> u128 {
        self.weights.iter().map(|&(h, _)| h).max().unwrap_or(1)
    }

    fn trace(
        &self,
        model: &mut Model,
        mut term: impl FnMut(&mut Model, u128) -> Result<Integer>,
    ) -> Result<Integer> {
        let mut sum = zero();
        for &(h, weight) in &self.weights {
            let value = term(model, h)?;
            sum = add(sum, multiply(model, weight, value)?)?;
        }
        let (magnitude, remainder) = sum.magnitude.divide(model, self.length)?;
        if remainder != 0 {
            return Err(Error::plain(
                "cyclic tensor projector trace division is inexact",
            ));
        }
        Ok(Integer {
            magnitude,
            negative: sum.negative,
        })
    }

    fn parity_traces(
        &self,
        model: &mut Model,
        parts: impl Fn(u128) -> (Integer, Integer),
    ) -> Result<[Integer; 2]> {
        let traces = |model: &mut Model, supertrace| {
            self.trace(model, |model, h| {
                let (a, b) = parts(h);
                let b = if supertrace || h.is_multiple_of(2) {
                    b.negate()
                } else {
                    b
                };
                power(model, add(a, b)?, self.length / h)
            })
        };
        let supertrace = traces(model, true)?;
        let trace = traces(model, false)?;
        Ok([
            half(add(trace, supertrace)?)?,
            half(add(trace, supertrace.negate())?)?,
        ])
    }
}

fn dimension(p: &[i128]) -> Integer {
    Integer::from_i128((p.len() - 1) as i128)
}

fn require_rank(rank: Integer) -> Result<usize> {
    if rank.negative && !rank.is_zero() {
        return Err(Error::plain("cyclic tensor projector rank is negative"));
    }
    let degree = rank
        .magnitude
        .small()
        .and_then(|value| usize::try_from(value).ok())
        .ok_or_else(degree_limit)?;
    require_degree(degree)?;
    Ok(degree)
}

fn check_primitive_lower_bound(n: u128) -> Result<()> {
    // With at least two input basis elements, the aperiodic binary words
    // give free cyclic orbits in the tensor basis. For n ≥ 4, fewer than
    // 2^(n−1) words have period ≤ n/2. Each free orbit contributes φ(d)
    // to the order-d projector, and one to the normalized Lie rank.
    if n >= 4 && (n > 128 || (1u128 << (n - 1)) / n > 2 * MAX_DEGREE as u128) {
        return Err(degree_limit());
    }
    Ok(())
}

fn divisors(model: &mut Model, factors: &[(u128, u32)], squarefree: bool) -> Result<Vec<u128>> {
    let mut values = vec![1];
    for &(prime, exponent) in factors {
        let width = values.len();
        let height = if squarefree { 1 } else { exponent };
        let count = width * (height as usize + 1);
        if count > MAX_DIVISORS {
            return Err(Error::plain(
                "cyclic tensor divisor count exceeds the storage limit",
            ));
        }
        model.spend(count - width)?;
        let mut factor = 1;
        for _ in 0..height {
            factor *= prime;
            for i in 0..width {
                values.push(values[i] * factor);
            }
        }
    }
    Ok(values)
}

fn totient(n: u128, factors: &[(u128, u32)]) -> u128 {
    factors
        .iter()
        .filter(|&&(prime, _)| n.is_multiple_of(prime))
        .fold(n, |value, &(prime, _)| value / prime * (prime - 1))
}

fn mobius_sign(mut n: u128, factors: &[(u128, u32)]) -> Option<bool> {
    let mut negative = false;
    for &(prime, _) in factors {
        if n.is_multiple_of(prime) {
            n /= prime;
            if n.is_multiple_of(prime) {
                return None;
            }
            negative = !negative;
        }
    }
    Some(negative)
}

fn natural(value: u128) -> Integer {
    Integer {
        magnitude: Natural::from_small(value),
        negative: false,
    }
}
