use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;

// `remaining` is monic and squarefree, possibly the constant 1.
pub(crate) fn distinct_degrees(
    model: &mut Model,
    field: Field<'_>,
    q: Option<u64>,
    mut remaining: Polynomial,
    mut collect: impl FnMut(&mut Model, Polynomial, usize) -> Result<()>,
) -> Result<()> {
    let x = field.polynomial(&[0, 1]);
    let mut x_power = x.clone();
    let mut degree = 1;
    while degree <= (remaining.len() - 1) / 2 {
        x_power = frobenius(model, field, q, x_power, &remaining)?;
        let difference = poly::sum(model, field, &x_power, &x, true)?;
        model.spend(remaining.len() * field.cells())?;
        let block = poly::monic_gcd(model, field, remaining.clone(), difference)?;
        if block.len() > 1 {
            remaining = poly::exact_divide(model, field, &remaining, &block)?;
            x_power = poly::remainder(model, field, x_power, &remaining)?;
            collect(model, block, degree)?;
        }
        degree += 1;
    }
    // All remaining factors have degree at least `degree`, so a nonconstant
    // remainder here must be irreducible.
    if remaining.len() > 1 {
        let degree = remaining.len() - 1;
        collect(model, remaining, degree)?;
    }
    Ok(())
}

pub(crate) fn factor_squarefree(
    model: &mut Model,
    field: Field<'_>,
    monic: Polynomial,
) -> Result<Vec<Polynomial>> {
    // The total degree of the factors held in the stack and output is bounded
    // by the input degree. Each factor has one extra leading coefficient.
    field.check_polynomial(monic.len().saturating_mul(2))?;
    let q = field.ring.size_of_residue_field().ok();
    let mut factors = Vec::new();
    let mut candidates = Candidates(0);
    distinct_degrees(model, field, q, monic, |model, block, degree| {
        let mut pending = vec![block];
        while let Some(block) = pending.pop() {
            if block.len() == degree + 1 {
                factors.push(block);
                continue;
            }
            // Cantor–Zassenhaus separates by quadratic character or, in
            // characteristic two, by absolute trace. Charge every attempt.
            loop {
                model.spend(block.len() * field.cells())?;
                let candidate = candidates.polynomial(field, block.len() - 1);
                let splitter = if field.ring.prime() == 2 {
                    binary_trace(model, field, &block, degree, candidate)?
                } else {
                    let sign = quadratic_sign(model, field, q, &block, degree, &candidate)?;
                    poly::sum(model, field, &sign, &[field.one()], true)?
                };
                model.spend(block.len() * field.cells())?;
                let factor = poly::monic_gcd(model, field, block.clone(), splitter)?;
                if factor.len() > 1 && factor.len() < block.len() {
                    pending.push(poly::exact_divide(model, field, &block, &factor)?);
                    pending.push(factor);
                    break;
                }
            }
        }
        Ok(())
    })?;
    Ok(factors)
}

fn binary_trace(
    model: &mut Model,
    field: Field<'_>,
    block: &[Scalar],
    degree: usize,
    mut power: Polynomial,
) -> Result<Polynomial> {
    let mut trace = power.clone();
    for _ in 1..field.ring.degree() * degree {
        power = poly::multiply_mod(model, field, &power, &power, block)?;
        trace = poly::sum(model, field, &trace, &power, false)?;
    }
    Ok(trace)
}

struct Candidates(u64);

impl Candidates {
    fn polynomial(&mut self, field: Field<'_>, size: usize) -> Polynomial {
        let values = (0..size)
            .map(|_| {
                (0..field.ring.degree())
                    .map(|_| self.next() % field.ring.prime())
                    .collect()
            })
            .collect();
        poly::trim(values)
    }

    // SplitMix64 provides a reproducible stream even when |K| exceeds u64.
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut value = self.0;
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^ (value >> 31)
    }
}

fn frobenius(
    model: &mut Model,
    field: Field<'_>,
    q: Option<u64>,
    mut value: Polynomial,
    modulus: &[Scalar],
) -> Result<Polynomial> {
    if let Some(q) = q {
        return poly::power_mod(model, field, &value, q, modulus);
    }
    for _ in 0..field.ring.degree() {
        value = poly::power_mod(model, field, &value, field.ring.prime(), modulus)?;
    }
    Ok(value)
}

pub(crate) fn quadratic_sign(
    model: &mut Model,
    field: Field<'_>,
    q: Option<u64>,
    block: &[Scalar],
    degree: usize,
    unit: &[Scalar],
) -> Result<Polynomial> {
    if let Some(size) = q.and_then(|q| q.checked_pow(degree as u32)) {
        return poly::power_mod(model, field, unit, (size - 1) / 2, block);
    }
    // The base-p digits of (p^(e*d)−1)/2 all equal (p−1)/2.
    let p = field.ring.prime();
    let digit_power = poly::power_mod(model, field, unit, (p - 1) / 2, block)?;
    model.spend(digit_power.len() * field.cells())?;
    let mut sign = digit_power.clone();
    for _ in 1..field.ring.degree() * degree {
        let raised = poly::power_mod(model, field, &sign, p, block)?;
        sign = poly::multiply_mod(model, field, &raised, &digit_power, block)?;
    }
    Ok(sign)
}

type Factor = (usize, usize, Polynomial);

pub(crate) fn factor_with_multiplicity(
    model: &mut Model,
    field: Field<'_>,
    p: &[Scalar],
) -> Result<Vec<Factor>> {
    let mut remaining = poly::monic(model, field, p.to_vec())?;
    if remaining.is_empty() {
        return Err(Error::plain(
            "cannot factor the zero polynomial over the residue field",
        ));
    }
    let mut places = Vec::new();
    let mut scale = 1usize;
    while remaining.len() > 1 {
        let derivative = poly::derivative(model, field, &remaining)?;
        if derivative.is_empty() {
            // Over a finite (perfect) field, an inseparable polynomial is
            // a pth power. Recover its pth root even in extension fields.
            let characteristic = field.ring.prime() as usize;
            remaining = remaining
                .iter()
                .step_by(characteristic)
                .map(|c| field.ring.frobenius(model, c, field.ring.degree() - 1))
                .collect::<Result<_>>()?;
            scale *= characteristic;
            continue;
        }
        let common = poly::monic_gcd(model, field, remaining.clone(), derivative)?;
        // f/gcd(f,f') contains precisely the irreducibles whose orders
        // are nonzero mod p, each once. Strip their entire multiplicity.
        let squarefree = poly::exact_divide(model, field, &remaining, &common)?;
        for factor in factor_squarefree(model, field, squarefree)? {
            let multiplicity = remove_factor(model, field, &mut remaining, &factor)?;
            if multiplicity == 0 {
                return Err(Error::plain("finite-field factor failed exact division"));
            }
            places.push((factor.len() - 1, multiplicity * scale, factor));
        }
    }
    places.sort_unstable_by(|(d, _, p), (e, _, q)| d.cmp(e).then_with(|| p.cmp(q)));
    Ok(places)
}

pub(crate) fn remove_factor(
    model: &mut Model,
    field: Field<'_>,
    p: &mut Polynomial,
    factor: &[Scalar],
) -> Result<usize> {
    if factor.len() < 2 {
        return Err(Error::plain(
            "finite-field factor must have positive degree",
        ));
    }
    let mut multiplicity = 0;
    while p.len() >= factor.len() {
        let (quotient, remainder) = poly::divide(model, field, p.clone(), factor)?;
        if !remainder.is_empty() {
            break;
        }
        *p = quotient;
        multiplicity += 1;
    }
    Ok(multiplicity)
}
