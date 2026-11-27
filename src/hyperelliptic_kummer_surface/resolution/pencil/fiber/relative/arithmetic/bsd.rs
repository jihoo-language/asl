use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::polarized::{Divisor, Generator};
use crate::integer_smith::{self as smith, Matrix, Track};
use crate::model::Model;
use crate::wide_integer::Natural;

use super::{k3, Components, Integer, Kodaira, LFunction, Pencil, RelativeJacobian, LIMBS};

impl RelativeJacobian {
    pub(crate) fn bsd(&self, model: &mut Model, pencil: &Pencil, name: &str) -> Result<String> {
        let LFunction {
            q,
            polynomial,
            rank,
            ns_rank,
            node_orbits,
            places,
            ..
        } = self.l_data(model, pencil)?;
        let powers = powers_of_q(q, polynomial.len() - 1)?;
        let root_number = reciprocity(model, &polynomial, &powers)?;
        if root_number != if rank.is_multiple_of(2) { 1 } else { -1 } {
            return Err(Error::plain(
                "K3 pencil root number disagrees with central vanishing parity",
            ));
        }
        let numerators = taylor_numerators(model, &polynomial, &powers)?;
        let numerator = numerators[rank];
        if numerators[..rank].iter().any(|c| !c.is_zero())
            || numerator.negative
            || numerator.is_zero()
        {
            return Err(Error::plain(
                "K3 pencil L polynomial failed its central order and positivity checks",
            ));
        }
        let denominator = powers[polynomial.len() - 1];

        let (tamagawa, orbit_product, trivial_det) = trivial_lattice(model, pencil, &places)?;
        let weil = pencil.curve.weil(model)?;
        let artin_tate = k3::artin_tate::data(model, q, weil, &node_orbits, 1)?;
        let special_with_q = numerator
            .multiply(i128::from(q))
            .map_err(|_| special_overflow())?;
        let predicted = special_with_q
            .multiply_wide(unsigned(orbit_product))
            .map_err(|_| special_overflow())?;
        let actual = denominator
            .multiply_wide(unsigned(artin_tate.magnitude))
            .map_err(|_| special_overflow())?;
        if artin_tate.rank != ns_rank as u128 || !equal(predicted, actual) {
            return Err(Error::plain(
                "K3 pencil BSD and Artin–Tate products disagree",
            ));
        }

        let leading = fraction(model, numerator, denominator.magnitude)?;
        let quotient_denominator = denominator
            .magnitude
            .multiply_u128(tamagawa)
            .ok_or_else(special_overflow)?;
        let quotient = fraction(model, special_with_q, quotient_denominator)?;
        let artin_quotient = fraction(
            model,
            unsigned(artin_tate.magnitude),
            Natural::from_small(trivial_det),
        )?;
        if quotient != artin_quotient {
            return Err(Error::plain(
                "K3 pencil BSD quotient failed the trivial-lattice comparison",
            ));
        }

        let sign = if ns_rank.is_multiple_of(2) { "-" } else { "" };
        let product = artin_tate.magnitude;
        let central = numerators
            .into_iter()
            .map(|c| fraction(model, c, denominator.magnitude))
            .collect::<Result<Vec<_>>>()?
            .join(", ");
        Ok(format!(
            "rvK3PencilBSD[{name}] = q={q}, rank={rank}, root number={root_number}, L*={leading}, Tamagawa product={tamagawa}, component orbit product={orbit_product}, arithmetic trivial |det|={trivial_det}, signed NS discriminant times Brauer order={sign}{product}, BSD quotient={quotient}, Taylor coefficients=[{central}]"
        ))
    }
}

fn trivial_lattice(
    model: &mut Model,
    pencil: &Pencil,
    places: &[(Kodaira, Components)],
) -> Result<(u128, u128, u128)> {
    let fiber = pencil.fiber_class(model)?;
    let section = Divisor::generator(
        model,
        &pencil.basis,
        Generator::Trope(pencil.section_masks.trailing_zeros() as usize),
    )?;
    let fiber_square = fiber.pair(model, &fiber)?;
    let section_square = section.pair(model, &section)?;
    let pairing = fiber.pair(model, &section)?;
    if (fiber_square, section_square, pairing) != (0, -2, 1) {
        return Err(inconsistent_lattice());
    }
    // F,O generate a unimodular hyperbolic plane, orthogonal to the
    // components that do not meet O.
    let mut determinant = (fiber_square * section_square - pairing * pairing).unsigned_abs();
    let mut tamagawa = 1u128;
    let mut orbit_product = 1u128;
    for (kind, action) in places {
        let component_det = component_determinant(model, *kind, action)?;
        let local_orbits = action.cycles.iter().try_fold(1u128, |product, &length| {
            product
                .checked_mul(length as u128)
                .ok_or_else(special_overflow)
        })?;
        let local_tamagawa = action.tamagawa as u128;
        if local_orbits.checked_mul(local_tamagawa) != Some(component_det) {
            return Err(inconsistent_lattice());
        }
        determinant = determinant
            .checked_mul(component_det)
            .ok_or_else(special_overflow)?;
        tamagawa = tamagawa
            .checked_mul(local_tamagawa)
            .ok_or_else(special_overflow)?;
        orbit_product = orbit_product
            .checked_mul(local_orbits)
            .ok_or_else(special_overflow)?;
    }
    Ok((tamagawa, orbit_product, determinant))
}

fn component_determinant(model: &mut Model, kind: Kodaira, action: &Components) -> Result<u128> {
    let orbits = match kind {
        Kodaira::Good => Vec::new(),
        Kodaira::Multiplicative(n) => {
            let split = action.split.ok_or_else(inconsistent_lattice)?;
            let mut orbits = Vec::new();
            for component in 1..n {
                let image = if split { component } else { n - component };
                if image == component {
                    orbits.push(vec![component]);
                } else if image > component {
                    orbits.push(vec![component, image]);
                }
            }
            orbits
        }
        Kodaira::Star(0) => {
            // The D4 center is fixed. The remaining three vertices are
            // leaves, so the cubic's cycle lengths determine the folded
            // Gram up to permuting those leaves.
            let mut first = 0;
            let orbits = action
                .cycles
                .iter()
                .map(|&length| {
                    let last = first + length;
                    let orbit = (first..last).collect::<Vec<_>>();
                    first = last;
                    orbit
                })
                .collect::<Vec<_>>();
            if first != 4 || orbits.first().map(Vec::as_slice) != Some(&[0]) {
                return Err(inconsistent_lattice());
            }
            orbits
        }
        _ => return Err(inconsistent_lattice()),
    };
    let mut actual_lengths: Vec<_> = orbits.iter().map(Vec::len).collect();
    let mut expected_lengths = action.cycles.clone();
    actual_lengths.sort_unstable();
    expected_lengths.sort_unstable();
    if actual_lengths != expected_lengths {
        return Err(inconsistent_lattice());
    }

    let rank = orbits.len();
    let mut gram = Matrix::new(model, rank, rank)?;
    for (i, left) in orbits.iter().enumerate() {
        for (j, right) in orbits[..=i].iter().enumerate() {
            model.spend(left.len() * right.len())?;
            let mut pairing = 0;
            for &a in left {
                for &b in right {
                    pairing += if a == b {
                        -2
                    } else {
                        i128::from(match kind {
                            Kodaira::Multiplicative(_) => a.abs_diff(b) == 1,
                            Kodaira::Star(0) => a == 0 || b == 0,
                            _ => false,
                        })
                    };
                }
            }
            gram.write(model, i, j, pairing)?;
            if i != j {
                gram.write(model, j, i, pairing)?;
            }
        }
    }
    let result = smith::reduce(model, gram, Track::None, Track::None)?;
    if result.factors.len() != rank {
        return Err(inconsistent_lattice());
    }
    result
        .factors
        .into_iter()
        .try_fold(1u128, |product, factor| {
            product
                .checked_mul(factor.unsigned_abs())
                .ok_or_else(special_overflow)
        })
}

fn powers_of_q(q: u64, degree: usize) -> Result<Vec<Integer>> {
    let mut powers = Vec::with_capacity(degree + 1);
    let mut power = Integer::from_i128(1);
    powers.push(power);
    for _ in 0..degree {
        power = power
            .multiply(i128::from(q))
            .map_err(|_| special_overflow())?;
        powers.push(power);
    }
    Ok(powers)
}

fn reciprocity(model: &mut Model, polynomial: &[Integer], powers: &[Integer]) -> Result<i128> {
    let degree = polynomial.len() - 1;
    let leading = polynomial[degree];
    let sign = if leading.negative { -1 } else { 1 };
    // L(T)=sign*(qT)^degree*L(1/(q²T)).
    for i in 0..=degree / 2 {
        model.spend(LIMBS * LIMBS)?;
        let partner = polynomial[i]
            .multiply_wide(powers[degree - 2 * i])
            .and_then(|c| c.multiply(sign))
            .map_err(|_| special_overflow())?;
        if !equal(partner, polynomial[degree - i]) {
            return Err(Error::plain(
                "K3 pencil L polynomial failed functional equation reciprocity",
            ));
        }
    }
    Ok(sign)
}

fn taylor_numerators(
    model: &mut Model,
    polynomial: &[Integer],
    powers: &[Integer],
) -> Result<Vec<Integer>> {
    let degree = polynomial.len() - 1;
    let mut output = vec![Integer::from_i128(0); polynomial.len()];
    for (i, coefficient) in polynomial.iter().enumerate() {
        model.spend(LIMBS * LIMBS)?;
        let normalized = coefficient
            .multiply_wide(powers[degree - i])
            .map_err(|_| special_overflow())?;
        let mut binomial = 1i128;
        for (j, out) in output[..=i].iter_mut().enumerate() {
            model.spend(4 * LIMBS)?;
            let term = normalized
                .multiply(binomial)
                .map_err(|_| special_overflow())?;
            let term = if j.is_multiple_of(2) {
                term
            } else {
                term.negate()
            };
            *out = out.add(term).map_err(|_| special_overflow())?;
            binomial = binomial * (i - j) as i128 / (j + 1) as i128;
        }
    }
    Ok(output)
}

fn equal(a: Integer, b: Integer) -> bool {
    a.magnitude == b.magnitude && (a.negative == b.negative || a.is_zero())
}

fn unsigned(value: u128) -> Integer {
    Integer {
        magnitude: Natural::from_small(value),
        negative: false,
    }
}

fn fraction(model: &mut Model, numerator: Integer, denominator: Natural<LIMBS>) -> Result<String> {
    if numerator.is_zero() {
        return Ok("0".into());
    }
    let gcd = numerator.magnitude.gcd(model, denominator)?;
    let reduced = |model: &mut Model, value: Natural<LIMBS>| -> Result<u128> {
        value
            .divide_exact(model, gcd)?
            .and_then(Natural::small)
            .ok_or_else(special_overflow)
    };
    let magnitude = reduced(model, numerator.magnitude)?;
    let denominator = reduced(model, denominator)?;
    let sign = if numerator.negative { "-" } else { "" };
    if denominator == 1 {
        Ok(format!("{sign}{magnitude}"))
    } else {
        Ok(format!("{sign}{magnitude}/{denominator}"))
    }
}

fn inconsistent_lattice() -> Error {
    Error::plain("K3 pencil has inconsistent arithmetic trivial lattice data")
}

fn special_overflow() -> Error {
    Error::plain("K3 pencil BSD result exceeds the exact arithmetic range")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn determinant_of_triality_folded_d4_is_three() {
        let mut model = Model::new();
        let action = Components {
            cycles: vec![1, 3],
            tamagawa: 1,
            split: None,
        };
        let determinant = component_determinant(&mut model, Kodaira::Star(0), &action)
            .unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(determinant, 3);
    }
}
