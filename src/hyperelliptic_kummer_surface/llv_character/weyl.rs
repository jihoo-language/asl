use crate::error::{Error, Result};
use crate::model::Model;

use super::{arithmetic, Magnitude, Signed, LARGEST_RHO, RANK, WEYL_LIMBS};

pub(super) fn dimension(model: &mut Model, shape: &[u128]) -> Result<Magnitude> {
    let mut weight = [0; RANK];
    weight[..shape.len()].copy_from_slice(shape);
    let mut ratio = WeylProduct {
        numerator: Magnitude::from_small(1),
        denominator: Magnitude::from_small(1),
    };
    for (i, &part) in weight.iter().enumerate() {
        let rho = (2 * (RANK - 1 - i) + 1) as u128;
        ratio.multiply(model, natural_sum([part, part, rho])?, rho)?;
        for (j, &other) in weight.iter().enumerate().skip(i + 1) {
            let difference = (j - i) as u128;
            let sum = (2 * RANK - 1 - i - j) as u128;
            ratio.multiply(model, natural_sum([part - other, difference])?, difference)?;
            ratio.multiply(model, natural_sum([part, other, sum])?, sum)?;
        }
    }
    if ratio.denominator.small() != Some(1) {
        return Err(Error::plain(
            "type B12 Weyl dimension did not reduce to an integer",
        ));
    }
    Ok(ratio.numerator)
}

struct WeylProduct {
    numerator: Magnitude,
    denominator: Magnitude,
}

impl WeylProduct {
    fn multiply(
        &mut self,
        model: &mut Model,
        mut factor: Magnitude,
        denominator: u128,
    ) -> Result<()> {
        let mut divisor = Magnitude::from_small(denominator);
        if factor == divisor {
            return Ok(());
        }
        model.spend(WEYL_LIMBS)?;
        cancel(model, &mut factor, &mut divisor)?;
        cancel(model, &mut self.numerator, &mut divisor)?;
        cancel(model, &mut factor, &mut self.denominator)?;
        self.numerator = self
            .numerator
            .checked_product(factor)
            .ok_or_else(arithmetic)?;
        self.denominator = self
            .denominator
            .checked_product(divisor)
            .ok_or_else(arithmetic)?;
        Ok(())
    }
}

fn cancel(model: &mut Model, numerator: &mut Magnitude, denominator: &mut Magnitude) -> Result<()> {
    if denominator.small() == Some(1) {
        return Ok(());
    }
    let common = numerator.gcd(model, *denominator)?;
    if common.small() == Some(1) {
        return Ok(());
    }
    *numerator = numerator
        .divide_exact(model, common)?
        .ok_or_else(arithmetic)?;
    *denominator = denominator
        .divide_exact(model, common)?
        .ok_or_else(arithmetic)?;
    Ok(())
}

fn natural_sum<const N: usize>(terms: [u128; N]) -> Result<Magnitude> {
    let mut result = Magnitude::default();
    for term in terms {
        result
            .add_multiple(Magnitude::from_small(term), 1)
            .ok_or_else(arithmetic)?;
    }
    Ok(result)
}

pub(super) fn specialized_numerator(model: &mut Model, shape: &[u128]) -> Result<Vec<Signed>> {
    let max = shape.first().copied().unwrap_or(0) as usize + LARGEST_RHO;
    let mut numerator = vec![Signed::from_i128(0); max + 1];
    let highest: [usize; RANK] =
        std::array::from_fn(|i| shape.get(i).copied().unwrap_or(0) as usize + LARGEST_RHO - i);
    let doubled: [u128; RANK] = std::array::from_fn(|i| (2 * highest[i] + 1) as u128);
    // With x₂,…,x₁₂ tending to 1, their determinant rows become odd
    // powers of the shifted weights. The remaining row is the interval
    // character x⁻ᵏ+⋯+xᵏ. The cofactor is an odd Vandermonde product.
    for (excluded, &degree) in highest.iter().enumerate() {
        let cofactor = odd_vandermonde(model, &doubled, &[excluded])?;
        let term = Signed {
            magnitude: cofactor,
            negative: !excluded.is_multiple_of(2),
        };
        for coefficient in numerator.iter_mut().take(degree + 1) {
            *coefficient = coefficient.add(term).map_err(|_| arithmetic())?;
        }
    }
    Ok(numerator)
}

pub(super) fn odd_vandermonde(
    model: &mut Model,
    doubled: &[u128; RANK],
    excluded: &[usize],
) -> Result<Magnitude> {
    let mut product = Magnitude::from_small(1);
    for (i, &left) in doubled.iter().enumerate() {
        if excluded.contains(&i) {
            continue;
        }
        model.spend(WEYL_LIMBS * (RANK - i))?;
        product = product.multiply_u128(left).ok_or_else(arithmetic)?;
        for (j, &right) in doubled.iter().enumerate().skip(i + 1) {
            if !excluded.contains(&j) {
                product = product
                    .multiply_u128(left * left - right * right)
                    .ok_or_else(arithmetic)?;
            }
        }
    }
    Ok(product)
}

pub(super) fn specialize(
    model: &mut Model,
    shape: &[u128],
    denominator: &[Signed],
) -> Result<Vec<Magnitude>> {
    let numerator = specialized_numerator(model, shape)?;
    let highest = numerator.len() - 1 - LARGEST_RHO;
    let mut weights = vec![Magnitude::default(); highest + 1];
    let divisor = denominator[LARGEST_RHO].magnitude;
    for charge in (0..=highest).rev() {
        let degree = charge + LARGEST_RHO;
        let mut coefficient = numerator[degree];
        for (next, &weight) in weights.iter().enumerate().skip(charge + 1) {
            let distance = degree.abs_diff(next);
            if weight.is_zero() || distance > LARGEST_RHO {
                continue;
            }
            model.spend(WEYL_LIMBS)?;
            let correction = denominator[distance]
                .multiply_wide(Signed {
                    magnitude: weight,
                    negative: true,
                })
                .map_err(|_| arithmetic())?;
            coefficient = coefficient.add(correction).map_err(|_| arithmetic())?;
        }
        if coefficient.negative && !coefficient.is_zero() {
            return Err(Error::plain(
                "LLV Weyl specialization produced a negative weight multiplicity",
            ));
        }
        weights[charge] = coefficient
            .magnitude
            .divide_exact(model, divisor)?
            .ok_or_else(|| Error::plain("LLV Weyl specialization did not divide exactly"))?;
    }
    Ok(weights)
}
