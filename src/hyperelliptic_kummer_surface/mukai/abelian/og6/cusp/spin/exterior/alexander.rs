use std::rc::Rc;

use super::smith::{exterior_dimension, exterior_group_ranks_modulo};
use super::{character, jordan_from_character, residue, top_pair, Algebra, Form, Terms, EVEN_RANK};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::Limit;
use crate::hyperelliptic_kummer_surface::rational_sparse::{decimal, Integer, WORK_LIMBS};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::syntax::{Label, ParsedExpr};

pub(crate) const MAX_LIFT: usize = EVEN_RANK / 2 + 1;
pub(crate) const FIBER_DEGREE: usize = EVEN_RANK;

pub(crate) enum Command {
    Data(Label, i128, usize, Option<usize>),
    Jets(Label, i128, usize, u128, Option<usize>),
    Witness(Label, i128, Option<usize>, Box<[Vec<ParsedExpr>; 2]>),
}

pub(crate) struct Witness {
    pub(crate) degree: usize,
    pub(crate) lifts: [Vec<Rc<Form>>; 2],
}

#[derive(Clone, Copy)]
pub(crate) struct LocalSystem {
    pub(crate) cover: i128,
    pub(crate) prime: Option<Prime>,
}

impl LocalSystem {
    pub(crate) fn new(
        model: &mut Model,
        cover: i128,
        characteristic: Option<usize>,
    ) -> Result<Self> {
        let prime = match characteristic {
            None => None,
            Some(2) => Some(Prime::validated(2)),
            Some(p) if p < 3 || p.is_multiple_of(2) => {
                return Err(Error::plain(
                    "Kuga–Satake Alexander coefficient characteristic must be prime",
                ));
            }
            Some(p) => Some(Prime::new(model, p)?),
        };
        Ok(Self { cover, prime })
    }

    fn field_name(self) -> String {
        self.prime
            .map_or_else(|| "Q".into(), |p| p.get().to_string())
    }
}

struct Profile {
    active_rank: usize,
    prime: Option<Prime>,
}

impl Profile {
    fn exponents(&self, model: &mut Model, degree: usize) -> Result<Vec<(usize, u128)>> {
        if degree > EVEN_RANK {
            return Ok(Vec::new());
        }
        if self.active_rank == 0 || matches!(degree, 0 | EVEN_RANK) {
            return Ok(vec![(1, exterior_dimension(degree))]);
        }
        let rank = self.active_rank;
        if let Some(p) = self.prime {
            let mut powers = exterior_group_ranks_modulo(model, degree, rank, p.get())?;
            // Nilpotent powers beyond the bound are zero, including the pth.
            powers.extend([0, 0]);
            Ok(powers
                .windows(3)
                .enumerate()
                .filter_map(|(i, triple)| {
                    let count = triple[0] + triple[2] - 2 * triple[1];
                    (count != 0).then_some((i + 1, count))
                })
                .collect())
        } else {
            let weights = character(model, rank, degree)?;
            Ok(jordan_from_character(&weights, rank))
        }
    }

    fn jet_dimension(&self, model: &mut Model, degree: usize, length: u128) -> Result<u128> {
        Ok(self
            .exponents(model, degree)?
            .into_iter()
            .map(|(size, count)| (size as u128).min(length) * count)
            .sum())
    }
}

impl Algebra {
    fn alexander_profile(
        &self,
        model: &mut Model,
        limit: &Limit,
        system: LocalSystem,
    ) -> Result<Profile> {
        let active_rank = match system.prime {
            _ if system.cover == 0 => 0,
            Some(p)
                if system
                    .cover
                    .unsigned_abs()
                    .is_multiple_of(u128::from(p.get())) =>
            {
                0
            }
            Some(p) => {
                let factors = self.action(model, limit)?.smith(model)?;
                // In the integral spin-pair basis the pairs with p|d_i become
                // inert, and the remaining d_i can be rescaled to units.
                factors.iter().filter(|d| d.residue(p.get()) != 0).count()
            }
            None => self.action(model, limit)?.smith(model)?.len(),
        };
        Ok(Profile {
            active_rank,
            prime: system.prime,
        })
    }

    pub(crate) fn alexander_report(
        &self,
        model: &mut Model,
        label: &str,
        limit: &Limit,
        system: LocalSystem,
        degree: usize,
    ) -> Result<String> {
        if degree > EVEN_RANK {
            return Err(Error::plain(
                "Kuga–Satake Alexander fiber degree must be between 0 and 64",
            ));
        }
        let profile = self.alexander_profile(model, limit, system)?;
        let exponents = profile.exponents(model, degree)?;
        let dimension = exterior_dimension(degree);
        let field = system.field_name();
        let cover = system.cover;
        Ok(format!("rvKSAlexander[{label},{cover},{degree},{field}] = dimension={dimension}, exponents={exponents:?}, differential_ranks={exponents:?}"))
    }

    pub(crate) fn alexander_jets_report(
        &self,
        model: &mut Model,
        label: &str,
        limit: &Limit,
        system: LocalSystem,
        degree: usize,
        length: u128,
    ) -> Result<String> {
        if degree > EVEN_RANK + 1 {
            return Err(Error::plain(
                "Kuga–Satake Alexander jet degree must be between 0 and 65",
            ));
        }
        if length == 0 {
            return Err(Error::plain(
                "Kuga–Satake Alexander jet length must be positive",
            ));
        }
        let profile = self.alexander_profile(model, limit, system)?;
        let fiber_kernel = profile.jet_dimension(model, degree, length)?;
        let fiber_cokernel = if let Some(lower) = degree.checked_sub(1) {
            profile.jet_dimension(model, lower, length)?
        } else {
            0
        };
        let dimension = fiber_kernel + fiber_cokernel;
        let field = system.field_name();
        let cover = system.cover;
        Ok(format!("rvKSAlexanderJets[{label},{cover},{degree},{length},{field}] = dimension={dimension}, fiber_kernel={fiber_kernel}, fiber_cokernel={fiber_cokernel}"))
    }

    pub(crate) fn alexander_witness_report(
        &self,
        model: &mut Model,
        label: &str,
        limit: &Limit,
        system: LocalSystem,
        witness: Witness,
    ) -> Result<String> {
        let [left, dual] = &witness.lifts;
        let order = left.len();
        let forward = Integer::from_i128(system.cover);
        let obstruction = self
            .check_alexander_lift(model, limit, system.prime, forward, witness.degree, left)
            .map_err(|error| Error::plain(format!("forward Alexander lift: {}", error.message)))?;
        self.check_alexander_lift(
            model,
            limit,
            system.prime,
            forward.negate(),
            EVEN_RANK - witness.degree,
            dual,
        )
        .map_err(|error| Error::plain(format!("adjoint Alexander lift: {}", error.message)))?;
        let terms = obstruction.into_iter().map(Ok);
        let pairing = top_pair(model, terms, &dual[0])?;
        let (value, certified) = if let Some(p) = system.prime {
            let value = residue(pairing, p);
            (value.to_string(), value != 0)
        } else {
            let sign = if pairing.negative { "-" } else { "" };
            (
                format!("{sign}{}", decimal(pairing.magnitude)),
                !pairing.is_zero(),
            )
        };
        let cover = system.cover;
        let field = system.field_name();
        Ok(format!("rvKSAlexanderWitness[{label},{cover},{field}] = order={order}, obstruction_pairing={value}, certified_nonzero={certified}"))
    }

    fn check_alexander_lift(
        &self,
        model: &mut Model,
        limit: &Limit,
        prime: Option<Prime>,
        exponent: Integer,
        degree: usize,
        coefficients: &[Rc<Form>],
    ) -> Result<Terms> {
        let action = if matches!(degree, 0 | EVEN_RANK)
            || exponent.is_zero()
            || prime.is_some_and(|p| residue(exponent, p) == 0)
        {
            None
        } else {
            Some(self.action(model, limit)?)
        };
        let mut previous = Terms::new();
        for (i, coefficient) in coefficients.iter().enumerate() {
            let image = if let Some(action) = action {
                action.field_exponential_terms(
                    model,
                    exponent,
                    coefficient,
                    degree > EVEN_RANK / 2,
                    prime,
                )?
            } else {
                coefficient.widen()?
            };
            let mut residual = image.clone();
            model.spend(
                (previous.len() + residual.len() + coefficient.terms().len()) * WORK_LIMBS,
            )?;
            for (mask, value) in previous {
                Form::accumulate(&mut residual, mask, value)?;
            }
            for &(mask, value) in coefficient.terms() {
                Form::accumulate(&mut residual, mask, value.resize()?.negate())?;
            }
            if residual
                .values()
                .any(|&value| prime.is_none_or(|p| residue(value, p) != 0))
            {
                return Err(Error::plain(format!(
                    "coefficient {i} of the Alexander lift equation is nonzero"
                )));
            }
            previous = image;
        }
        Ok(previous)
    }
}
