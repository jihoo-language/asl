use crate::error::{Error, Result};
use crate::legendre::{inverse, split_prime};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_unram::{add_mod, mul_mod, reduce, require_scalar, sub_mod, Ring};

use super::numerator_factors;

pub(super) fn unit_root(
    model: &mut Model,
    prime: Prime,
    psi: i128,
    precision: usize,
) -> Result<String> {
    let p = u64::from(prime.get());
    if precision == 0 {
        return Err(Error::plain(
            "Dwork quartic unit-root precision must be positive",
        ));
    }
    let modulus = u32::try_from(precision)
        .ok()
        .and_then(|s| p.checked_pow(s))
        .ok_or_else(|| Error::plain("Dwork quartic unit-root modulus p^s must fit in u64"))?;

    let hasse_witt = coefficient(model, p, p, p, psi)?;
    if hasse_witt == 0 {
        return Ok("Hasse–Witt=0, nonordinary, no Frobenius unit root".into());
    }
    let degree = modulus;
    let previous_degree = degree / p;
    let current = if precision == 1 {
        hasse_witt
    } else {
        coefficient(model, p, modulus, degree, psi)?
    };
    let previous = coefficient(model, p, modulus, previous_degree, psi)?;
    let unit = mul_mod(current, inverse(model, previous, modulus)?, modulus);
    if precision > 1 {
        let earlier = coefficient(model, p, previous_degree, previous_degree / p, psi)?;
        let prior = mul_mod(
            previous % previous_degree,
            inverse(model, earlier, previous_degree)?,
            previous_degree,
        );
        if unit % previous_degree != prior {
            return Err(ratio_error());
        }
    }
    if unit % p != hasse_witt {
        return Err(Error::plain(
            "Dwork quartic unit root disagrees with Hasse–Witt",
        ));
    }
    Ok(format!(
        "mod {modulus}, Hasse–Witt={hasse_witt}, b_{degree}={current}, b_{previous_degree}={previous}, unit={unit}"
    ))
}

pub(crate) fn closed_unit(model: &mut Model, ring: &Ring, input: &[u64]) -> Result<String> {
    let p = ring.prime();
    if p == 2 {
        return Err(super::odd_characteristic());
    }
    let field = ring.residue_ring(model)?;
    let residue = field.residue(input);
    if field.power(model, &residue, 4)? == field.scalar(1) {
        return Err(super::singular());
    }
    let hasse_witt = coefficient_in_ring(model, &field, p, &residue)?;
    let modulus = ring.modulus();
    let q = ring.size_of_residue_field()?;
    let psi = ring.teichmuller(model, input)?;
    let head = format!(
        "rvDworkK3ClosedUnit mod {modulus} = q={q}, Teichmüller psi={psi:?}, Hasse–Witt={hasse_witt:?}"
    );
    if hasse_witt == field.scalar(0) {
        return Ok(format!("{head}, nonordinary, no Frobenius unit root"));
    }
    let current = if ring.precision() == 1 {
        hasse_witt
    } else {
        coefficient_in_ring(model, ring, modulus, &psi)?
    };
    let previous_degree = modulus / p;
    let previous = coefficient_in_ring(model, ring, previous_degree, &psi)?;
    let ratio = twisted_quotient(model, ring, &current, &previous)?;
    if ring.precision() > 1 {
        let earlier = coefficient_in_ring(model, ring, previous_degree / p, &psi)?;
        let prior = twisted_quotient(model, ring, &previous, &earlier)?;
        if ratio
            .iter()
            .zip(prior)
            .any(|(a, b)| a % previous_degree != b % previous_degree)
        {
            return Err(ratio_error());
        }
    }
    let norm = ring.conjugate_sum_or_product(model, &ratio, true)?;
    let unit = require_scalar(&norm)?;
    Ok(format!("{head}, Cartier factor={ratio:?}, unit={unit}"))
}

fn twisted_quotient(model: &mut Model, ring: &Ring, a: &[u64], b: &[u64]) -> Result<Vec<u64>> {
    let shifted = ring.frobenius(model, b, 1)?;
    let inverse = ring.inverse(model, &shifted)?;
    ring.multiply(model, a, &inverse)
}

fn coefficient_in_ring(
    model: &mut Model,
    ring: &Ring,
    degree: u64,
    psi: &[u64],
) -> Result<Vec<u64>> {
    let n = degree - 1;
    let modulus = ring.modulus();
    let z = ring.scale(model, psi, sub_mod(0, 4 % modulus, modulus))?;
    let fourth = ring.power(model, &z, 4)?;
    let mut sum = ring.scalar(0);
    multinomial_terms(model, ring.prime(), modulus, n, |model, term| {
        sum = ring.multiply(model, &sum, &fourth)?;
        if term != 0 {
            sum = ring.add(model, &sum, &ring.scalar(term))?;
        }
        Ok(())
    })?;
    let tail = ring.power(model, &z, n % 4)?;
    ring.multiply(model, &sum, &tail)
}

fn coefficient(model: &mut Model, p: u64, modulus: u64, degree: u64, psi: i128) -> Result<u64> {
    let n = degree - 1;
    let z = mul_mod(
        sub_mod(0, 4 % modulus, modulus),
        reduce(psi, modulus),
        modulus,
    );
    let square = mul_mod(z, z, modulus);
    let fourth = mul_mod(square, square, modulus);
    let tails = [1, z, square, mul_mod(square, z, modulus)];
    let mut sum = 0;
    multinomial_terms(model, p, modulus, n, |_, term| {
        sum = add_mod(mul_mod(sum, fourth, modulus), term, modulus);
        Ok(())
    })?;
    Ok(mul_mod(sum, tails[(n % 4) as usize], modulus))
}

// In descending powers of z^4, emit the coefficients of b_(n+1)/z^(n mod 4).
fn multinomial_terms(
    model: &mut Model,
    p: u64,
    modulus: u64,
    n: u64,
    mut emit: impl FnMut(&mut Model, u64) -> Result<()>,
) -> Result<()> {
    let last = n / 4;
    let work = usize::try_from(last + 1)
        .map_err(|_| Error::plain("Dwork quartic period truncation is out of range"))?;
    model.spend(work)?;
    let mut powers = vec![1];
    while *powers.last().unwrap() < modulus / p {
        powers.push(p * powers.last().unwrap());
    }
    let mut unit = 1;
    let mut valuation = 0usize;
    for j in 0..=last {
        let term = powers
            .get(valuation)
            .map_or(0, |&factor| mul_mod(unit, factor, modulus));
        emit(model, term)?;
        if j != last {
            // u_(j+1)/u_j is the product of four consecutive
            // factors divided by (j+1)^4. Retain the unit even when
            // this term vanishes at the requested precision.
            for factor in numerator_factors(n, j) {
                let (value, order) = split_prime(factor, p);
                valuation += order;
                unit = mul_mod(unit, value % modulus, modulus);
            }
            let (divisor, order) = split_prime(j + 1, p);
            valuation = valuation.checked_sub(4 * order).ok_or_else(|| {
                Error::plain("Dwork quartic multinomial valuation became negative")
            })?;
            let divisor_inverse = inverse(model, divisor % modulus, modulus)?;
            let divisor_inverse_squared = mul_mod(divisor_inverse, divisor_inverse, modulus);
            let divisor_inverse_fourth =
                mul_mod(divisor_inverse_squared, divisor_inverse_squared, modulus);
            unit = mul_mod(unit, divisor_inverse_fourth, modulus);
        }
    }
    Ok(())
}

fn ratio_error() -> Error {
    Error::plain("Dwork quartic successive unit-root residues disagree")
}
