use crate::ec_poly::Field;
use crate::elliptic_modular::decimal;
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hyperelliptic_kummer_surface::mordell::Mordell;
use crate::model::Model;
use crate::rv_ore::gcd;
use crate::wide_integer::Natural;

use super::super::field;
use super::{fraction_text, special_overflow, unsigned, BadPlace, Integer, Kodaira, Pencil, LIMBS};

mod regulator;

const MAX_DIVISOR_SCAN: u128 = 100_000;

type Magnitude = Natural<LIMBS>;
type Ratio = (u128, u128);

struct Factors {
    rank: usize,
    sha_dimension: usize,
    sha_regulator: Ratio,
    height_denominator: u128,
    additive_place: String,
    multiplicative_place: String,
}

impl Pencil {
    pub(crate) fn bsd_factors(&self, model: &mut Model, name: &str) -> Result<String> {
        let factors = Factors::new(model, self)?;
        let bound = factors.sha_divisibility_bound(model)?;
        let common = factors.summary();
        let separation = if factors.rank == 0 {
            let order = bound.small().ok_or_else(special_overflow)?;
            if !compatible_sha(order, factors.sha_dimension) {
                return Err(inconsistent_factors());
            }
            format!("Sha order={order}, regulator=1")
        } else if let Some(limit) = bound.small().filter(|n| n.isqrt() <= MAX_DIVISOR_SCAN) {
            let orders = sha_divisors(model, limit, factors.sha_dimension)?;
            if let [unique] = orders.as_slice() {
                let regulator = fraction_text(divide_ratio(factors.sha_regulator, *unique)?);
                format!("Sha order={unique}, regulator={regulator}")
            } else {
                format!("Sha order candidates={orders:?}")
            }
        } else {
            "Sha order candidates exceed divisor scan limit".into()
        };
        Ok(format!(
            "rvK3PencilBSDFactors[{name}] = {common}, Sha order divides={}, {separation}",
            display(bound),
        ))
    }
}

impl Mordell {
    pub(crate) fn bsd_factors(
        &self,
        model: &mut Model,
        name: &str,
        trope: u8,
        multiple: i128,
    ) -> Result<String> {
        if multiple == 0 {
            return Err(Error::plain(
                "K3 pencil BSD section multiple must be nonzero",
            ));
        }
        let factors = Factors::new(model, &self.pencil)?;
        if factors.rank != 1 {
            return Err(Error::plain(
                "K3 pencil BSD section refinement requires arithmetic Mordell-Weil rank one",
            ));
        }
        self.certified_trope_height(model, trope)?;
        let candidates = factors.section_candidates(model)?;
        let step = Magnitude::from_small(multiple.unsigned_abs());
        let height = Magnitude::product(multiple, multiple);
        let mut outcomes = Vec::with_capacity(candidates.len());
        for (index, sha) in candidates {
            // The height-one trope difference has free index index. Its
            // multiple has index |multiple|*index, with no need to factor n.
            let free_index = step.multiply_u128(index).ok_or_else(special_overflow)?;
            let regulator = fraction_text(divide_ratio(factors.sha_regulator, sha)?);
            outcomes.push((display(free_index), regulator, sha));
        }
        let separation = if let [(index, regulator, sha)] = outcomes.as_slice() {
            format!("generator index={index}, regulator={regulator}, Sha order={sha}")
        } else {
            let rendered = outcomes
                .iter()
                .map(|(index, regulator, sha)| format!("({index}, {regulator}, {sha})"))
                .collect::<Vec<_>>()
                .join(", ");
            format!("candidates (generator index, regulator, Sha order)=[{rendered}]")
        };
        Ok(format!(
            "rvK3MordellBSDFactors[{name}, {trope}, {multiple}] = {}, zero trope={}, section height={}, {separation}",
            factors.summary(), self.zero_mask(), display(height),
        ))
    }
}

impl Factors {
    fn new(model: &mut Model, pencil: &Pencil) -> Result<Self> {
        if pencil.torsion_roots(model)?.len() != 3 {
            return Err(Error::plain(
                "K3 pencil BSD factors require all three affine two-torsion points",
            ));
        }
        let relative = pencil.relative(model)?;
        let data = relative.bsd_data(model, pencil)?;
        let rank = data.l_function.rank;
        let field = field(&pencil.curve);
        let places = &data.l_function.places;
        let (additive, multiplicative) = torsion_places(places, field.ring.prime())?;
        let height_denominator = height_denominator(places)?;
        let additive_place = describe_place(model, field, additive)?;
        let multiplicative_place = describe_place(model, field, multiplicative)?;
        let sha_dimension = pencil.selmer_space(model)?.sha_dimension(rank)?;
        let sha_regulator = multiply_ratio(data.quotient, 16)?;
        Ok(Self {
            rank,
            sha_dimension,
            sha_regulator,
            height_denominator,
            additive_place,
            multiplicative_place,
        })
    }

    fn summary(&self) -> String {
        format!(
            "rank={}, torsion=(Z/2)^2, prime-to-characteristic certificate={}, characteristic certificate={}, height denominator={}, dim Sha[2]={}, Sha*Reg={}",
            self.rank, self.additive_place, self.multiplicative_place, self.height_denominator,
            self.sha_dimension, fraction_text(self.sha_regulator),
        )
    }

    fn sha_divisibility_bound(&self, model: &mut Model) -> Result<Magnitude> {
        // D times the arithmetic height Gram is integral, so D^r Reg
        // is a positive integer. The exact BSD product gives the bound.
        let (numerator, denominator) = self.sha_regulator;
        let mut bound = Magnitude::from_small(numerator);
        model.spend(self.rank * LIMBS * LIMBS)?;
        for _ in 0..self.rank {
            bound = bound
                .multiply_u128(self.height_denominator)
                .ok_or_else(special_overflow)?;
        }
        let (bound, remainder) = bound.divide(model, denominator)?;
        if remainder != 0 || bound.is_zero() {
            return Err(inconsistent_factors());
        }
        if self.sha_dimension != 0 {
            let order = 1u128
                .checked_shl(self.sha_dimension as u32)
                .ok_or_else(special_overflow)?;
            if bound.divide(model, order)?.1 != 0 {
                return Err(inconsistent_factors());
            }
        }
        Ok(bound)
    }

    fn section_candidates(&self, model: &mut Model) -> Result<Vec<(u128, u128)>> {
        let denominator = self.height_denominator;
        let maximum = denominator.isqrt();
        model.spend(usize::try_from(maximum).map_err(|_| special_overflow())?)?;
        let mut candidates = Vec::new();
        for index in 1..=maximum {
            let square = index * index;
            if !denominator.is_multiple_of(square) {
                continue;
            }
            let (sha, denominator) = multiply_ratio(self.sha_regulator, square)?;
            if denominator == 1 && compatible_sha(sha, self.sha_dimension) {
                candidates.push((index, sha));
            }
        }
        if candidates.is_empty() {
            return Err(inconsistent_factors());
        }
        Ok(candidates)
    }
}

fn torsion_places(places: &[BadPlace], prime: u64) -> Result<(&BadPlace, &BadPlace)> {
    // In I0* the smooth identity component is additive. The reduction
    // kernel is pro-p for prime-to-p torsion, so the torsion bounded by
    // D4's component group of order four is exactly the rational E[2].
    let additive = places
        .iter()
        .find(|p| matches!(p.kind, Kodaira::Star(0)))
        .ok_or_else(|| {
            Error::plain("K3 pencil torsion certificate requires a rational I0* place")
        })?;
    if additive.action.tamagawa != 4 || additive.action.cycles != [1, 1, 1, 1] {
        return Err(Error::plain(
            "K3 pencil rational I0* components disagree with split two-torsion",
        ));
    }
    // Over a split unramified completion a multiplicative I_n fiber has
    // toric identity and the Tate identity subgroup has no p-torsion;
    // p not dividing n also excludes p-torsion in the components.
    let multiplicative = places.iter().find(|p| matches!(p.kind, Kodaira::Multiplicative(n) if !(n as u64).is_multiple_of(prime)))
        .ok_or_else(|| Error::plain("K3 pencil torsion certificate requires a rational I_n place with n prime to the characteristic"))?;
    Ok((additive, multiplicative))
}

fn describe_place(model: &mut Model, field: Field<'_>, place: &BadPlace) -> Result<String> {
    let kodaira = place.kind.symbol();
    if let Some(factor) = &place.uniformizer {
        let parameter = field.subtract(model, &field.zero(), &factor[0])?;
        Ok(format!("{kodaira} at t={parameter:?}"))
    } else {
        Ok(format!("{kodaira} at infinity"))
    }
}

fn height_denominator(places: &[BadPlace]) -> Result<u128> {
    places.iter().try_fold(1, |lcm, place| {
        let denominator = match place.kind {
            Kodaira::Good => 1,
            Kodaira::Multiplicative(n) => n as u128,
            Kodaira::Star(0) => 2, // The inverse D4 Cartan matrix has half-integral entries.
            _ => {
                return Err(Error::plain(
                    "K3 pencil BSD factors lack a height denominator for this fiber",
                ))
            }
        };
        (lcm / gcd(lcm, denominator))
            .checked_mul(denominator)
            .ok_or_else(special_overflow)
    })
}

fn sha_divisors(model: &mut Model, bound: u128, dimension: usize) -> Result<Vec<u128>> {
    let maximum = bound.isqrt();
    model.spend(maximum as usize)?;
    let mut results = Vec::new();
    for divisor in 1..=maximum {
        if bound.is_multiple_of(divisor) {
            let partner = bound / divisor;
            if compatible_sha(divisor, dimension) {
                results.push(divisor);
            }
            if partner != divisor && compatible_sha(partner, dimension) {
                results.push(partner);
            }
        }
    }
    if results.is_empty() {
        return Err(inconsistent_factors());
    }
    results.sort_unstable();
    Ok(results)
}

fn compatible_sha(order: u128, dimension: usize) -> bool {
    // The canonical principal polarization comes from the rational zero
    // divisor. Cassels–Tate is alternating and perfect on this finite Sha.
    let square_root = order.isqrt();
    order != 0
        && square_root * square_root == order
        && if dimension == 0 {
            order % 2 == 1
        } else {
            order.trailing_zeros() as usize >= dimension
        }
}

fn multiply_ratio((numerator, denominator): Ratio, multiple: u128) -> Result<Ratio> {
    let common = gcd(multiple, denominator);
    let numerator = numerator
        .checked_mul(multiple / common)
        .ok_or_else(special_overflow)?;
    Ok((numerator, denominator / common))
}

fn divide_ratio((numerator, denominator): Ratio, divisor: u128) -> Result<Ratio> {
    let common = gcd(numerator, divisor);
    let denominator = denominator
        .checked_mul(divisor / common)
        .ok_or_else(special_overflow)?;
    Ok((numerator / common, denominator))
}

fn display(magnitude: Magnitude) -> String {
    // Use the existing exact decimal formatter for a wide positive integer.
    decimal(Integer {
        magnitude,
        ..unsigned(0)
    })
}

fn inconsistent_factors() -> Error {
    Error::plain("K3 pencil BSD factors are inconsistent with the height and Selmer certificates")
}

#[cfg(test)]
mod tests {
    use super::super::Components;
    use super::*;

    #[test]
    fn local_torsion_certificate_checks_both_places() {
        let place = |kind, tamagawa, cycles| BadPlace {
            uniformizer: None,
            kind,
            action: Components {
                cycles,
                tamagawa,
                split: None,
            },
        };
        let additive = || place(Kodaira::Star(0), 4, vec![1; 4]);
        let multiplicative = |n| place(Kodaira::Multiplicative(n), 1, vec![]);
        assert!(torsion_places(&[additive(), multiplicative(7)], 7).is_err());
        assert!(torsion_places(&[multiplicative(2)], 7).is_err());
        assert!(torsion_places(&[additive(), multiplicative(2)], 7).is_ok());
        let inconsistent = place(Kodaira::Star(0), 2, vec![1, 1, 2]);
        assert!(torsion_places(&[inconsistent, multiplicative(2)], 7).is_err());
    }

    #[test]
    fn index_candidates_filter_height_integrality_and_sha_two_rank() {
        let mut model = Model::new();
        let mut factors = Factors {
            rank: 1,
            sha_dimension: 2,
            sha_regulator: (1, 1),
            height_denominator: 16,
            additive_place: String::new(),
            multiplicative_place: String::new(),
        };
        assert_eq!(
            factors
                .section_candidates(&mut model)
                .unwrap_or_else(|error| panic!("{error}")),
            [(2, 4), (4, 16)]
        );
        factors.sha_dimension = 0;
        factors.sha_regulator = (9, 4);
        assert_eq!(
            factors
                .section_candidates(&mut model)
                .unwrap_or_else(|error| panic!("{error}")),
            [(2, 9)]
        );
        factors.rank = 0;
        factors.sha_regulator = (9, 1);
        let rank_zero = factors
            .sha_divisibility_bound(&mut model)
            .unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(rank_zero.small(), Some(9));
    }
}
