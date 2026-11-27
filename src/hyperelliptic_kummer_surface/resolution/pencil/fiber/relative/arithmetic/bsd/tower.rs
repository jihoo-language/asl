use std::rc::Rc;

use crate::elliptic_modular::decimal;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_ore::gcd;

use super::super::super::sections::Section;
use super::super::{polynomial_overflow, root_multiplicity};
use super::{
    central_value, checked_bsd_quotient, fraction, fraction_text, k3, special_overflow,
    trivial_lattice, BadPlace, Components, Integer, Kodaira, LFunction, Pencil, LIMBS,
};

const MAX_SQUAREFREE_TRIAL: u128 = 100_000;

struct TowerTerm {
    q: i128,
    polynomial: Vec<Integer>,
    rank: usize,
    ns_rank: u128,
    trivial_rank: usize,
    artin_tate: u128,
    root_number: i128,
    leading: String,
    tamagawa: u128,
    orbit_product: u128,
    trivial_det: u128,
    bsd: (u128, u128),
}

impl Pencil {
    pub(crate) fn l_power(&self, model: &mut Model, name: &str, degree: u128) -> Result<String> {
        positive(degree)?;
        let relative = self.relative(model)?;
        let base = relative.l_data(model, self)?;
        let term = TowerTerm::new(model, self, &base, degree)?;
        let coefficients = term
            .polynomial
            .iter()
            .map(|&c| decimal(c))
            .collect::<Vec<_>>()
            .join(", ");
        let sign = if term.ns_rank.is_multiple_of(2) {
            "-"
        } else {
            ""
        };
        Ok(format!(
            "rvK3PencilLPower[{name}, {degree}] = Q={}, coefficients=[{coefficients}], NS rank={}, trivial rank={}, Mordell-Weil rank={}, root number={}, L*={}, Tamagawa product={}, component orbit product={}, arithmetic trivial |det|={}, signed NS discriminant times Brauer order={sign}{}, BSD quotient={}",
            term.q, term.ns_rank, term.trivial_rank, term.rank, term.root_number, term.leading,
            term.tamagawa, term.orbit_product, term.trivial_det, term.artin_tate, fraction_text(term.bsd),
        ))
    }

    pub(crate) fn artin_tate_ratio(
        &self,
        model: &mut Model,
        name: &str,
        lower: u128,
        upper: u128,
    ) -> Result<String> {
        positive(lower)?;
        positive(upper)?;
        if !upper.is_multiple_of(lower) {
            return Err(Error::plain(
                "K3 pencil constant-field degrees must satisfy m dividing n",
            ));
        }
        let base = self.relative(model)?.l_data(model, self)?;
        let lower_term = TowerTerm::new(model, self, &base, lower)?;
        let upper_term = TowerTerm::new(model, self, &base, upper)?;
        let result = comparison(&lower_term, &upper_term)?.map_or_else(
            || "NS rank jump, integral NS discriminants require separate data".into(),
            |ratio| {
                format!(
                    "integral NS lattice fixed, Brauer and Sha order ratio={}",
                    fraction_text(ratio)
                )
            },
        );
        Ok(format!(
            "rvK3PencilArtinTateRatio[{name}, {lower}, {upper}] = NS ranks=[{}, {}], Artin-Tate magnitudes=[{}, {}], {result}",
            lower_term.ns_rank, upper_term.ns_rank, lower_term.artin_tate, upper_term.artin_tate,
        ))
    }

    pub(crate) fn sha_power(
        &self,
        model: &mut Model,
        name: &str,
        degree: u128,
        generators: &[Rc<Section>],
    ) -> Result<String> {
        positive(degree)?;
        let base = self.relative(model)?.l_data(model, self)?;
        let lower = TowerTerm::new(model, self, &base, 1)?;
        let upper = TowerTerm::new(model, self, &base, degree)?;
        let ratio = comparison(&lower, &upper)?;
        let base_order = self.regulator_sha_order(model, generators)?;
        let result = if let Some(ratio) = ratio {
            let order = exact_order(base_order, ratio)?;
            let structure = elementary_divisors(model, order)?
                .map_or_else(String::new, |structure| format!(", Sha={structure}"));
            format!("integral NS lattice fixed, Sha order ratio={}, extension Sha order={order}{structure}", fraction_text(ratio))
        } else {
            "NS rank jump, extension Sha order unresolved".into()
        };
        Ok(format!(
            "rvK3PencilShaPower[{name}, {degree}] = base Sha order={base_order}, NS ranks=[{}, {}], Artin-Tate magnitudes=[{}, {}], {result}",
            lower.ns_rank, upper.ns_rank, lower.artin_tate, upper.artin_tate,
        ))
    }
}

impl TowerTerm {
    fn new(model: &mut Model, pencil: &Pencil, base: &LFunction, degree: u128) -> Result<Self> {
        let weil = pencil.curve.weil(model)?;
        let artin_tate = k3::artin_tate::data(model, base.q, weil, &base.node_orbits, degree)?;
        let q = artin_tate.field_order;
        let degree_usize = usize::try_from(degree).map_err(|_| polynomial_overflow())?;
        let polynomial = power_roots(model, &base.polynomial, degree_usize)?;
        let rank = root_multiplicity(model, &polynomial, q)?;
        let places = base
            .places
            .iter()
            .map(|place| power_place(place, degree))
            .collect::<Result<Vec<_>>>()?;
        let trivial_rank = 2 + places.iter().map(|p| p.action.cycles.len()).sum::<usize>();
        if (rank + trivial_rank) as u128 != artin_tate.rank {
            return Err(Error::plain(
                "K3 pencil tower failed the powered Shioda–Tate rank check",
            ));
        }
        let (tamagawa, orbit_product, trivial_det) = trivial_lattice(model, pencil, &places)?;
        let central = central_value(model, &polynomial, q, rank)?;
        let numerator = central.numerators[rank];
        let leading = fraction(model, numerator, central.denominator)?;
        let bsd = checked_bsd_quotient(
            model,
            q,
            numerator,
            central.denominator,
            (tamagawa, orbit_product, trivial_det),
            artin_tate.magnitude,
        )?;
        Ok(Self {
            q,
            polynomial,
            rank,
            ns_rank: artin_tate.rank,
            trivial_rank,
            artin_tate: artin_tate.magnitude,
            root_number: central.root_number,
            leading,
            tamagawa,
            orbit_product,
            trivial_det,
            bsd,
        })
    }
}

fn power_place(place: &BadPlace, degree: u128) -> Result<BadPlace> {
    let mut cycles = Vec::new();
    for &length in &place.action.cycles {
        if length == 0 {
            return Err(inconsistent());
        }
        let divisions = gcd(length as u128, degree) as usize;
        cycles.extend(std::iter::repeat_n(length / divisions, divisions));
    }
    cycles.sort_unstable();
    let fixed = cycles.iter().take_while(|&&length| length == 1).count();
    let (tamagawa, split) = match place.kind {
        Kodaira::Good => (1, None),
        Kodaira::Multiplicative(_) => (
            1 + fixed,
            Some(place.action.split.ok_or_else(inconsistent)? || degree.is_multiple_of(2)),
        ),
        // The D4 center is always fixed. Counting it instead of the
        // identity leaf gives exactly 1 + the remaining fixed leaves.
        Kodaira::Star(0) => (fixed, None),
        _ => return Err(inconsistent()),
    };
    Ok(BadPlace {
        uniformizer: place.uniformizer.clone(),
        kind: place.kind,
        action: Components {
            cycles,
            tamagawa,
            split,
        },
    })
}

fn power_roots(model: &mut Model, polynomial: &[Integer], exponent: usize) -> Result<Vec<Integer>> {
    let n = polynomial.len() - 1;
    let maximum = n.checked_mul(exponent).ok_or_else(polynomial_overflow)?;
    model.spend(maximum.saturating_mul(n).saturating_mul(LIMBS * LIMBS))?;
    let zero = Integer::from_i128(0);
    let mut traces = vec![zero; maximum + 1];
    // Newton's formulas transition to the companion recurrence once the
    // power index exceeds the degree. The base q^exponent bound keeps
    // the trace table small even for enormous extension fields.
    for k in 1..=maximum {
        let mut sum = zero;
        for j in 1..=n.min(k - 1) {
            sum = add_product(sum, polynomial[j], traces[k - j])?;
        }
        if k <= n {
            let term = polynomial[k]
                .multiply(k as i128)
                .map_err(|_| polynomial_overflow())?;
            sum = sum.add(term).map_err(|_| polynomial_overflow())?;
        }
        traces[k] = sum.negate();
    }
    model.spend(n * n * LIMBS * LIMBS)?;
    let mut powered = vec![Integer::from_i128(1)];
    for k in 1..=n {
        let mut sum = zero;
        for j in 1..=k {
            sum = add_product(sum, powered[k - j], traces[exponent * j])?;
        }
        let (magnitude, remainder) = sum.magnitude.divide_u32(k as u32);
        if remainder != 0 {
            return Err(Error::plain(
                "K3 pencil powered Newton sum did not divide exactly",
            ));
        }
        powered.push(Integer {
            magnitude,
            negative: !sum.negative && !magnitude.is_zero(),
        });
    }
    Ok(powered)
}

fn add_product(sum: Integer, a: Integer, b: Integer) -> Result<Integer> {
    sum.add(a.multiply_wide(b).map_err(|_| polynomial_overflow())?)
        .map_err(|_| polynomial_overflow())
}

fn comparison(lower: &TowerTerm, upper: &TowerTerm) -> Result<Option<(u128, u128)>> {
    if lower.ns_rank > upper.ns_rank {
        return Err(inconsistent());
    }
    if lower.ns_rank != upper.ns_rank {
        return Ok(None);
    }
    let divisor = gcd(lower.artin_tate, upper.artin_tate);
    Ok(Some((
        upper.artin_tate / divisor,
        lower.artin_tate / divisor,
    )))
}

fn exact_order(order: u128, (numerator, denominator): (u128, u128)) -> Result<u128> {
    if !order.is_multiple_of(denominator) {
        return Err(Error::plain(
            "K3 pencil tower ratio and certified base Sha do not give an integer order",
        ));
    }
    let order = (order / denominator)
        .checked_mul(numerator)
        .ok_or_else(special_overflow)?;
    if order == 0 || order.isqrt().pow(2) != order {
        return Err(Error::plain(
            "K3 pencil tower did not give a positive square Sha order",
        ));
    }
    Ok(order)
}

fn elementary_divisors(model: &mut Model, order: u128) -> Result<Option<String>> {
    if order == 1 {
        return Ok(Some("0".into()));
    }
    let half = order.isqrt();
    let maximum = half.isqrt();
    if maximum > MAX_SQUAREFREE_TRIAL {
        return Ok(None);
    }
    model.spend(maximum as usize)?;
    let squarefree = (2..=maximum).all(|factor| !half.is_multiple_of(factor * factor));
    Ok(squarefree.then(|| format!("(Z/{half})^2")))
}

fn positive(degree: u128) -> Result<()> {
    if degree == 0 {
        return Err(Error::plain(
            "K3 pencil constant-field degree must be positive",
        ));
    }
    Ok(())
}

fn inconsistent() -> Error {
    Error::plain("K3 pencil tower data have inconsistent ranks or component permutations")
}

#[cfg(test)]
mod tests {
    use super::super::equal;
    use super::*;

    #[test]
    fn powering_roots_composes_exactly_with_wide_coefficients() {
        let mut model = Model::new();
        let polynomial = [1, -5, -414, -9522, -60835, 6436343].map(Integer::from_i128);
        let block = |model: &mut Model, p: &[Integer], power| {
            power_roots(model, p, power).unwrap_or_else(|error| panic!("{error}"))
        };
        let two = block(&mut model, &polynomial, 2);
        let three = block(&mut model, &polynomial, 3);
        let six = block(&mut model, &polynomial, 6);
        let two_then_three = block(&mut model, &two, 3);
        let three_then_two = block(&mut model, &three, 2);
        for ((a, b), c) in six.into_iter().zip(two_then_three).zip(three_then_two) {
            assert!(equal(a, b) && equal(a, c));
        }
    }

    #[test]
    fn actual_group_orders_accept_rational_ratios_and_only_forced_divisors() {
        let mut model = Model::new();
        assert!(matches!(exact_order(4, (1, 4)), Ok(1)));
        assert!(exact_order(4, (1, 3)).is_err());
        assert!(exact_order(4, (2, 1)).is_err());
        assert_eq!(
            elementary_divisors(&mut model, 100).unwrap_or_else(|error| panic!("{error}")),
            Some("(Z/10)^2".into())
        );
        assert!(elementary_divisors(&mut model, 16)
            .unwrap_or_else(|error| panic!("{error}"))
            .is_none());
    }

    #[test]
    fn local_cycles_and_component_points_power_for_triality_and_nonsplit_fibers() {
        let check = |kind, cycles, split, exponent, expected, tamagawa| {
            let place = BadPlace {
                uniformizer: None,
                kind,
                action: Components {
                    cycles,
                    tamagawa: 1,
                    split,
                },
            };
            let power = power_place(&place, exponent).unwrap_or_else(|error| panic!("{error}"));
            assert_eq!(power.action.cycles, expected);
            assert_eq!(power.action.tamagawa, tamagawa);
        };
        check(Kodaira::Star(0), vec![1, 3], None, 2, vec![1, 3], 1);
        check(Kodaira::Star(0), vec![1, 3], None, 3, vec![1; 4], 4);
        check(
            Kodaira::Multiplicative(4),
            vec![1, 2],
            Some(false),
            2,
            vec![1; 3],
            4,
        );
    }
}
