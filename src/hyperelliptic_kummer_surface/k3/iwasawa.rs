use std::rc::Rc;

use crate::elliptic_modular::decimal;
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::pencil::{sections::Section, Pencil};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_ore::gcd;
use crate::rv_qp_quad::Wide;
use crate::wide_integer::Natural;

use super::{abelian_spectrum, cross_quartic, node_orbits, AbelianSpectrum};

mod modular;

use modular::Newton;

type Valuation = Wide<5>;

pub(super) struct Spectrum<'a> {
    q: u64,
    abelian: AbelianSpectrum,
    node_orbits: &'a [usize],
}

pub(super) struct Growth {
    prime: u128,
    rank: usize,
    rank_one: usize,
    geometric_rank: usize,
    residual_degree: usize,
    breakpoint: Breakpoint,
}

impl<'a> Spectrum<'a> {
    pub(super) fn new(
        model: &mut Model,
        q: u64,
        weil: &[i128],
        orbits: &'a [usize],
    ) -> Result<Self> {
        let abelian = abelian_spectrum(model, q, cross_quartic(q, weil)?)?;
        Ok(Self {
            q,
            abelian,
            node_orbits: orbits,
        })
    }

    pub(super) fn growth(&self, model: &mut Model, base: u128, prime: usize) -> Result<Growth> {
        if base == 0 {
            return Err(Error::plain(
                "Kummer K3 Iwasawa base degree must be positive",
            ));
        }
        let prime = checked_prime(model, prime)?;
        if u128::from(self.q).is_multiple_of(prime) {
            return Err(Error::plain(
                "Kummer K3 Iwasawa prime must differ from the characteristic",
            ));
        }
        let rank = self.rank(base, 1);
        let rank_one = self.rank(1, 1);
        let geometric_rank = 18
            + self
                .abelian
                .cyclotomic
                .iter()
                .map(|&(_, d)| d)
                .sum::<usize>();
        let residual_degree = self.abelian.transcendental.len() - 1;
        if geometric_rank + residual_degree != 22 || rank > geometric_rank {
            return Err(inconsistent());
        }
        self.require_fixed_rank(base, prime, rank)?;
        // All nonfixed cyclotomic directions have a nontrivial prime-to-ell
        // part after this fixed-rank certificate. Their shifts have valuation 0.
        let breakpoint = Breakpoint::new(model, self.q, &self.abelian.transcendental, base, prime)?;
        Ok(Growth {
            prime,
            rank,
            rank_one,
            geometric_rank,
            residual_degree,
            breakpoint,
        })
    }

    fn rank(&self, degree: u128, multiplier: u128) -> usize {
        let algebraic = self
            .abelian
            .cyclotomic
            .iter()
            .filter(|&&(order, _)| multiplier.is_multiple_of(order / gcd(degree, order)))
            .map(|&(_, rank)| rank)
            .sum::<usize>();
        let exceptional = self
            .node_orbits
            .iter()
            .map(|&length| {
                let length = length as u128;
                let split = gcd(degree, length);
                (split * gcd(length / split, multiplier)) as usize
            })
            .sum::<usize>();
        2 + algebraic + exceptional
    }

    fn require_fixed_rank(&self, degree: u128, prime: u128, rank: usize) -> Result<()> {
        let maximum = self
            .abelian
            .cyclotomic
            .iter()
            .map(|&(order, _)| order)
            .chain(self.node_orbits.iter().map(|&length| length as u128))
            .max()
            .unwrap_or(1);
        let mut multiplier = 1;
        let mut level = 0;
        while multiplier < maximum {
            level += 1;
            multiplier *= prime;
            let next = self.rank(degree, multiplier);
            if next != rank {
                return Err(Error::plain(format!("Kummer K3 Iwasawa Néron–Severi rank jumps from {rank} to {next} at level k={level}")));
            }
        }
        // Every ell-part of a cyclotomic order or permutation cycle length
        // divides this multiplier, so later powers add no fixed directions.
        Ok(())
    }
}

impl Growth {
    pub(super) fn report(
        &self,
        operation: &str,
        group: &str,
        name: &str,
        degree: u128,
        level: u128,
    ) -> Result<String> {
        let start = self.breakpoint.index();
        let slopes = self
            .breakpoint
            .slopes
            .iter()
            .map(|s| s.report(self.breakpoint.offset))
            .collect::<Vec<_>>()
            .join(", ");
        let relative = self.breakpoint.relative(level)?;
        let from_base = self.breakpoint.from_base(level)?;
        let zeros = 22 - self.rank - self.breakpoint.lambda;
        Ok(format!(
            "{operation}[{name}, {degree}, {}, {level}] = NS rank={}, geometric NS rank={}, integral NS lattice fixed throughout tower, transcendental degree={}, zero valuations={zeros}, breakpoint k0={start}, positive Newton valuations=[{slopes}], determinant valuations through breakpoint={:?}, mu=0, lambda={}, relative v_ell({group}) from breakpoint={}, relative v_ell({group}) from base={}",
            self.prime, self.rank, self.geometric_rank, self.residual_degree,
            self.breakpoint.prefix, self.breakpoint.lambda, decimal(relative), decimal(from_base),
        ))
    }
}

impl Pencil {
    pub(crate) fn sha_iwasawa(
        &self,
        model: &mut Model,
        name: &str,
        degree: u128,
        prime: usize,
        level: u128,
        generators: Option<&[Rc<Section>]>,
    ) -> Result<String> {
        let q = self.curve.ring.size_of_residue_field()?;
        let orbits = node_orbits(model, &self.curve, q)?;
        let weil = self.curve.weil(model)?;
        let spectrum = Spectrum::new(model, q, weil, &orbits)?;
        let growth = spectrum.growth(model, degree, prime)?;
        let mut report = growth.report("rvK3PencilShaIwasawa", "Sha", name, degree, level)?;
        if let Some(generators) = generators {
            let base_order = self.regulator_sha_order(model, generators)?;
            if growth.rank_one != growth.rank {
                report.push_str(&format!(", certified Sha order over F_q={base_order}, NS rank over F_q={}, absolute tower valuation unresolved", growth.rank_one));
            } else {
                // Nested fixed lattices of equal rank coincide integrally.
                // The ell-part of the Artin–Tate q factor is trivial.
                let reference = modular::certify(
                    model,
                    q,
                    &spectrum.abelian.transcendental,
                    1,
                    growth.prime,
                    0,
                )?;
                let known =
                    modular::valuation(growth.prime, base_order).ok_or_else(inconsistent)?;
                let offset = known + growth.breakpoint.prefix[0] - reference.determinant;
                if offset < 0 {
                    return Err(inconsistent());
                }
                let absolute = growth
                    .breakpoint
                    .from_base(level)?
                    .add(Valuation::from_i128(i128::from(offset)))?;
                if absolute.negative && !absolute.is_zero() {
                    return Err(inconsistent());
                }
                report.push_str(&format!(
                    ", certified Sha order over F_q={base_order}, v_ell(Sha)={}",
                    decimal(absolute)
                ));
            }
        }
        Ok(report)
    }
}

struct Breakpoint {
    slopes: Vec<modular::Slope>,
    offset: usize,
    prefix: Vec<i64>,
    lambda: usize,
}

impl Breakpoint {
    fn new(
        model: &mut Model,
        q: u64,
        polynomial: &[i128],
        mut base: u128,
        prime: u128,
    ) -> Result<Self> {
        let mut removed = 0usize;
        while base.is_multiple_of(prime) {
            base /= prime;
            removed += 1;
        }
        let mut level = 0;
        let mut order = modular::certify(model, q, polynomial, base, prime, level)?;
        let lambda = order.positive.iter().map(|s| s.count).sum::<usize>();
        let mut prefix = vec![order.determinant];
        while !order
            .positive
            .iter()
            .all(|&slope| slope.past_threshold(prime))
        {
            // A positive slope is at least 1/d. Below the threshold it
            // multiplies by ell, and equality crosses it in one step.
            // Here d<=4, which bounds the number of initial levels by d.
            level += 1;
            if level >= polynomial.len() {
                return Err(inconsistent());
            }
            let next = modular::certify(model, q, polynomial, base, prime, level)?;
            if next.determinant < order.determinant
                || next.positive.iter().map(|s| s.count).sum::<usize>() != lambda
            {
                return Err(inconsistent());
            }
            prefix.push(next.determinant);
            order = next;
        }
        let Newton {
            determinant,
            positive: slopes,
        } = order;
        let offset = removed.saturating_sub(level);
        let prefix = if offset != 0 {
            vec![determinant + (lambda * offset) as i64]
        } else {
            prefix.split_off(removed)
        };
        Ok(Self {
            slopes,
            offset,
            prefix,
            lambda,
        })
    }

    fn index(&self) -> u128 {
        (self.prefix.len() - 1) as u128
    }

    fn relative(&self, target: u128) -> Result<Valuation> {
        if let Some(steps) = target.checked_sub(self.index()) {
            let magnitude = Natural::from_small(steps)
                .multiply_u128(self.lambda as u128)
                .ok_or_else(inconsistent)?;
            return Ok(Valuation {
                magnitude,
                negative: false,
            });
        }
        let value = self.prefix[target as usize] - self.prefix[self.index() as usize];
        Ok(Valuation::from_i128(i128::from(value)))
    }

    fn from_base(&self, target: u128) -> Result<Valuation> {
        let offset = self.prefix[self.index() as usize] - self.prefix[0];
        self.relative(target)?
            .add(Valuation::from_i128(i128::from(offset)))
    }
}

fn checked_prime(model: &mut Model, prime: usize) -> Result<u128> {
    if prime == 2 {
        return Ok(2);
    }
    if u32::try_from(prime).is_err() {
        return Err(Error::plain("Kummer K3 Iwasawa prime must fit in u32"));
    }
    if prime < 3 || prime.is_multiple_of(2) {
        return Err(Error::plain("Kummer K3 Iwasawa prime must be prime"));
    }
    Ok(u128::from(Prime::new(model, prime)?.get()))
}

fn inconsistent() -> Error {
    Error::plain("Kummer K3 Iwasawa data failed their exact slope or rank checks")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equality_at_the_strict_threshold_requires_an_exact_initial_power() {
        let mut model = Model::new();
        // Normalized roots have product 1 and sum 2/3. Their shifts start
        // at valuation 1 and squaring gives valuation 5/2, not 2.
        let polynomial = [1, -2, 9];
        let growth = Breakpoint::new(&mut model, 3, &polynomial, 1, 2)
            .unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(growth.index(), 1);
        assert_eq!(growth.prefix, [2, 5]);
        assert_eq!(growth.lambda, 2);
        assert_eq!(growth.slopes[0].report(growth.offset), "5/2: 2");
        assert_eq!(
            growth
                .relative(0)
                .unwrap_or_else(|error| panic!("{error}"))
                .signed(),
            Some(-3)
        );
        assert_eq!(
            growth
                .from_base(0)
                .unwrap_or_else(|error| panic!("{error}"))
                .signed(),
            Some(0)
        );
        assert_eq!(
            growth
                .from_base(1000)
                .unwrap_or_else(|error| panic!("{error}"))
                .signed(),
            Some(2001)
        );

        let later = Breakpoint::new(&mut model, 3, &polynomial, 1u128 << 126, 2)
            .unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(later.index(), 0);
        assert_eq!(later.prefix, [255]);
        assert_eq!(later.slopes[0].report(later.offset), "255/2: 2");
    }

    #[test]
    fn fractional_valuations_below_threshold_power_at_exact_finite_levels() {
        let mut model = Model::new();
        let polynomial = [1, 4, 9];
        let growth = Breakpoint::new(&mut model, 3, &polynomial, 1, 2)
            .unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(growth.index(), 2);
        assert_eq!(growth.prefix, [1, 2, 6]);
        assert_eq!(growth.slopes[0].report(growth.offset), "3: 2");
        assert_eq!(
            growth
                .relative(1)
                .unwrap_or_else(|error| panic!("{error}"))
                .signed(),
            Some(-4)
        );
    }

    #[test]
    fn exact_cyclotomic_and_permutation_orders_detect_future_rank_jumps() {
        let spectrum = Spectrum {
            q: 3,
            abelian: AbelianSpectrum {
                cyclotomic: vec![(4, 2)],
                transcendental: vec![1, -2, 9],
            },
            node_orbits: &[1; 16],
        };
        assert_eq!(spectrum.rank(1, 1), 18);
        assert_eq!(spectrum.rank(1, 2), 18);
        let result = spectrum.require_fixed_rank(1, 2, 18);
        assert!(result
            .err()
            .unwrap()
            .to_string()
            .contains("jumps from 18 to 20 at level k=2"));
        assert!(spectrum.require_fixed_rank(1, 3, 18).is_ok());
        let mixed = Spectrum {
            node_orbits: &[1, 1, 2, 2, 2, 2, 2, 2, 2],
            ..spectrum
        };
        assert_eq!(mixed.rank(1, 1), 11);
        let result = mixed.require_fixed_rank(1, 2, 11);
        assert!(result
            .err()
            .unwrap()
            .to_string()
            .contains("jumps from 11 to 18 at level k=1"));
    }
}
