use std::collections::BTreeSet;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::integer_smith::{self as smith, Matrix, Track};
use crate::model::Model;
use crate::rv_ore::gcd;

use super::super::super::super::sections::{covering::MordellImage, Section};
use super::super::super::heights::{signed_ratio, HeightData};
use super::{
    compatible_sha, divide_ratio, field, fraction_text, multiply_ratio, special_overflow, Factors,
    Pencil, Ratio, MAX_DIVISOR_SCAN,
};

mod sieve;

const MAX_SIEVED_PRIMES: usize = 8;

impl Pencil {
    pub(crate) fn regulator(
        &self,
        model: &mut Model,
        name: &str,
        generators: &[Rc<Section>],
    ) -> Result<String> {
        self.regulator_result(model, name, generators)
            .map(|(report, _)| report)
    }

    pub(crate) fn regulator_sha_order(
        &self,
        model: &mut Model,
        generators: &[Rc<Section>],
    ) -> Result<u128> {
        self.regulator_result(model, "", generators)?.1.ok_or_else(|| {
            Error::plain("K3 pencil Sha tower requires an exact base Sha order from the regulator certificate")
        })
    }

    fn regulator_result(
        &self,
        model: &mut Model,
        name: &str,
        generators: &[Rc<Section>],
    ) -> Result<(String, Option<u128>)> {
        let factors = Factors::new(model, self)?;
        let rank = factors.rank;
        if generators.len() != rank {
            return Err(Error::plain(format!(
                "K3 pencil regulator needs {rank} sections for the arithmetic Mordell-Weil rank",
            )));
        }
        let relative = self.relative(model)?;
        let heights = relative.heights(model, field(&self.curve))?;
        let d = heights.multiplier;
        if factors.height_denominator != d {
            return Err(inconsistent_lattice());
        }
        let gram = Gram::new(model, self, &heights, generators)?;
        let denominator = d.checked_pow(rank as u32).ok_or_else(special_overflow)?;
        let common_divisor = gcd(gram.determinant, denominator);
        let gamma = (
            gram.determinant / common_divisor,
            denominator / common_divisor,
        );
        let mod_two = MordellImage::dimension_of(model, self, generators)?;
        let deficiency = (rank + 2)
            .checked_sub(mod_two)
            .ok_or_else(inconsistent_lattice)?;
        if deficiency > rank {
            return Err(inconsistent_lattice());
        }
        let summary = format!(
            "{}, {}, Gram={}, subgroup determinant={}, D^r*det={}, mod-2 image dimension={mod_two}/{}",
            factors.summary(), heights.provenance(), gram.display(d), fraction_text(gamma),
            gram.determinant, rank + 2,
        );
        let maximum = gram.determinant.isqrt();
        if maximum > MAX_DIVISOR_SCAN {
            return Ok((format!(
                "rvK3PencilRegulator[{name}] = {summary}, index candidates exceed divisor scan limit",
            ), None));
        }
        let candidates = index_candidates(model, &factors, gram.determinant, gamma, deficiency)?;
        let primes = odd_primes(model, &candidates)?;
        let mut excluded = BTreeSet::new();
        let mut evidence = Vec::new();
        for &prime in primes.iter().take(MAX_SIEVED_PRIMES) {
            let result = sieve::test(model, self, relative, generators, prime)?;
            if result.excluded {
                excluded.insert(prime);
            }
            evidence.push(result.report);
        }
        if primes.len() > MAX_SIEVED_PRIMES {
            evidence.push(format!(
                "prime limit leaves untested {:?}",
                &primes[MAX_SIEVED_PRIMES..]
            ));
        }
        let remaining: Vec<_> = candidates
            .iter()
            .copied()
            .filter(|(index, _)| excluded.iter().all(|&p| !index.is_multiple_of(p)))
            .collect();
        let indices: Vec<_> = candidates.iter().map(|&(index, _)| index).collect();
        let exact_sha = if let [(_, sha)] = remaining.as_slice() {
            Some(*sha)
        } else {
            None
        };
        let report = match remaining.as_slice() {
            [] => return Err(inconsistent_lattice()),
            [(index, sha)] => {
                let regulator = fraction_text(divide_ratio(gamma, index * index)?);
                let structure = if *sha == 1 {
                    ", Sha=0".to_owned()
                } else if Some(*sha) == 1u128.checked_shl(factors.sha_dimension as u32) {
                    format!(", Sha=(Z/2)^{}", factors.sha_dimension)
                } else {
                    String::new()
                };
                format!("exact index={index}, regulator={regulator}, Sha order={sha}{structure}")
            }
            _ => {
                let mut choices = Vec::new();
                for (index, sha) in remaining {
                    let regulator = fraction_text(divide_ratio(gamma, index * index)?);
                    choices.push(format!("({index}, {regulator}, {sha})"));
                }
                format!(
                    "unresolved (index, regulator, Sha order)=[{}]",
                    choices.join(", ")
                )
            }
        };
        Ok((format!(
            "rvK3PencilRegulator[{name}] = {summary}, arithmetic index candidates={indices:?}, good-fiber sieve=[{}], {report}",
            evidence.join("; "),
        ), exact_sha))
    }
}

struct Gram {
    entries: Vec<Vec<i128>>,
    determinant: u128,
}

impl Gram {
    fn new(
        model: &mut Model,
        pencil: &Pencil,
        heights: &HeightData,
        sections: &[Rc<Section>],
    ) -> Result<Self> {
        let rank = sections.len();
        let d = heights.multiplier as i128;
        let self_heights = sections
            .iter()
            .map(|section| heights.numerator(model, pencil, section))
            .collect::<Result<Vec<_>>>()?;
        let mut entries = vec![vec![0; rank]; rank];
        let mut matrix = Matrix::new(model, rank, rank)?;
        for j in 0..rank {
            let diagonal = self_heights[j] / d;
            entries[j][j] = diagonal;
            matrix.write(model, j, j, diagonal)?;
            for i in 0..j {
                let numerator = heights.twice_pair(
                    model,
                    pencil,
                    &sections[i],
                    &sections[j],
                    self_heights[i],
                    self_heights[j],
                )?;
                let entry = numerator / (2 * d);
                entries[i][j] = entry;
                entries[j][i] = entry;
                matrix.write(model, i, j, entry)?;
                matrix.write(model, j, i, entry)?;
            }
        }
        let invariants = smith::reduce(model, matrix, Track::None, Track::None)?.factors;
        // The height form on sections modulo torsion is positive definite;
        // for arbitrary lists its Gram is positive semidefinite.
        if invariants.len() != rank {
            return Err(Error::plain(
                "K3 pencil regulator sections are dependent modulo torsion",
            ));
        }
        let determinant = invariants
            .into_iter()
            .try_fold(1u128, |product, invariant| {
                product
                    .checked_mul(invariant.unsigned_abs())
                    .ok_or_else(special_overflow)
            })?;
        Ok(Self {
            entries,
            determinant,
        })
    }

    fn display(&self, denominator: u128) -> String {
        let rows = self
            .entries
            .iter()
            .map(|row| {
                let entries = row
                    .iter()
                    .map(|&n| signed_ratio(n, denominator))
                    .collect::<Vec<_>>();
                format!("[{}]", entries.join(", "))
            })
            .collect::<Vec<_>>();
        format!("[{}]", rows.join(", "))
    }
}

fn index_candidates(
    model: &mut Model,
    factors: &Factors,
    determinant: u128,
    gamma: Ratio,
    two_deficiency: usize,
) -> Result<Vec<(u128, u128)>> {
    let maximum = determinant.isqrt();
    model.spend(maximum as usize)?;
    let product_over_gamma =
        multiply_ratio(divide_ratio(factors.sha_regulator, gamma.0)?, gamma.1)?;
    let mut candidates = Vec::new();
    for index in 1..=maximum {
        let square = index * index;
        let valuation = index.trailing_zeros() as usize;
        let valid_two = if two_deficiency == 0 {
            valuation == 0
        } else {
            valuation >= two_deficiency
        };
        if !valid_two || !determinant.is_multiple_of(square) {
            continue;
        }
        let (sha, denominator) = multiply_ratio(product_over_gamma, square)?;
        if denominator == 1 && compatible_sha(sha, factors.sha_dimension) {
            candidates.push((index, sha));
        }
    }
    if candidates.is_empty() {
        return Err(inconsistent_lattice());
    }
    Ok(candidates)
}

fn odd_primes(model: &mut Model, candidates: &[(u128, u128)]) -> Result<Vec<u128>> {
    let mut primes = BTreeSet::new();
    for &(mut index, _) in candidates {
        index >>= index.trailing_zeros();
        let mut divisor = 3;
        while divisor * divisor <= index {
            model.spend(1)?;
            if index.is_multiple_of(divisor) {
                primes.insert(divisor);
                while index.is_multiple_of(divisor) {
                    index /= divisor;
                }
            }
            divisor += 2;
        }
        if index > 1 {
            primes.insert(index);
        }
    }
    Ok(primes.into_iter().collect())
}

fn inconsistent_lattice() -> Error {
    Error::plain(
        "K3 pencil Mordell-Weil lattice is inconsistent with the height, Selmer and BSD data",
    )
}
