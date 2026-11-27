use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::Result;
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
            // Cantor–Zassenhaus separates the factors on which a candidate
            // has quadratic character +1. Charge each attempt to the work
            // budget, including attempts that leave the block intact.
            loop {
                model.spend(block.len() * field.cells())?;
                let candidate = candidates.polynomial(field, block.len() - 1);
                let sign = quadratic_sign(model, field, q, &block, degree, &candidate)?;
                let positive = poly::sum(model, field, &sign, &[field.one()], true)?;
                model.spend(block.len() * field.cells())?;
                let factor = poly::monic_gcd(model, field, block.clone(), positive)?;
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
