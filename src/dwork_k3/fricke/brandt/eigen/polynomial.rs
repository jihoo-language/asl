use crate::ec_poly::Field;
use crate::error::Result;
use crate::finite_field_factor;
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_unram::Ring;

use super::linear::{self as lin, add, mul, negate, Matrix};
use super::{certificate, FactorKind, Fraction};

pub(super) fn multiply(left: &[i128], right: &[i128]) -> Result<Vec<i128>> {
    let mut product = vec![0; left.len() + right.len() - 1];
    for (i, &a) in left.iter().enumerate() {
        for (j, &b) in right.iter().enumerate() {
            product[i + j] = add(product[i + j], mul(a, b)?)?;
        }
    }
    Ok(product)
}

pub(super) fn evaluate(poly: &[i128], value: i128) -> Result<i128> {
    poly.iter()
        .rev()
        .try_fold(0, |prior, &c| add(mul(value, prior)?, c))
}

// Return the quotient and the fixed-length residue modulo a monic polynomial.
pub(super) fn divide(mut poly: Vec<i128>, monic: &[i128]) -> Result<(Vec<i128>, Vec<i128>)> {
    let degree = monic.len() - 1;
    let mut quotient = vec![0; poly.len().saturating_sub(degree)];
    for (j, coefficient) in quotient.iter_mut().enumerate().rev() {
        *coefficient = poly[j + degree];
        for (entry, &c) in poly[j..j + degree].iter_mut().zip(monic) {
            *entry = add(*entry, negate(mul(c, *coefficient)?)?)?;
        }
    }
    poly.resize(degree, 0);
    Ok((quotient, poly))
}

pub(super) fn projector(
    model: &mut Model,
    total: &[i128],
    factor: &[i128],
    powers: &[Matrix],
) -> Result<(Matrix, i128)> {
    let (cofactor, residual) = divide(total.to_vec(), factor)?;
    if residual.iter().any(|&c| c != 0) {
        return Err(certificate("exact factor division"));
    }
    let degree = factor.len() - 1;
    let mut columns = Vec::with_capacity(degree);
    for j in 0..degree {
        let shifted = vec![0; j]
            .into_iter()
            .chain(cofactor.iter().copied())
            .collect();
        columns.push(divide(shifted, factor)?.1);
    }
    let mut unit = vec![0; degree];
    unit[0] = 1;
    let inverse = lin::solve(model, &columns, &unit)?
        .ok_or_else(|| crate::error::Error::plain("Brandt factors are not pairwise coprime"))?;
    let (inverse, denominator) = lin::clear(&inverse)?;
    let crt = multiply(&inverse, &cofactor)?;
    let numerator = lin::combine(model, &crt, powers)?;
    let fraction = Fraction::new(lin::flat(&numerator), denominator);
    let n = numerator.len();
    let matrix = fraction
        .numerator
        .chunks(n)
        .map(|row| row.to_vec())
        .collect();
    Ok((matrix, fraction.denominator))
}

pub(super) fn squarefree(model: &mut Model, total: &[i128]) -> Result<bool> {
    let degree = total.len() - 1;
    let mut derivative = Vec::with_capacity(degree);
    for (j, &value) in total.iter().enumerate().skip(1) {
        derivative.push(mul(j as i128, value)?);
    }
    let mut columns = Vec::with_capacity(degree);
    for j in 0..degree {
        let shifted = vec![0; j]
            .into_iter()
            .chain(derivative.iter().copied())
            .collect();
        columns.push(divide(shifted, total)?.1);
    }
    let mut unit = vec![0; degree];
    unit[0] = 1;
    Ok(lin::solve(model, &columns, &unit)?.is_some())
}

pub(super) fn factor_kind(
    model: &mut Model,
    polynomial: &[i128],
    eigen_bound: i128,
) -> Result<FactorKind> {
    let degree = polynomial.len() - 1;
    if degree == 1 {
        return Ok(FactorKind::Irreducible);
    }
    for root in -eigen_bound..=eigen_bound {
        if evaluate(polynomial, root)? == 0 {
            return Ok(FactorKind::Reducible);
        }
    }
    if degree <= 3 {
        return Ok(FactorKind::Irreducible);
    }
    // Irreducibility of one monic prime reduction proves rational irreducibility.
    for p in [3, 5, 7, 11, 13, 17, 19, 23] {
        let ring = Ring::field_of_degree(model, Prime::validated(p), 1)?;
        let field = Field { ring: &ring };
        let reduction = field.polynomial(polynomial);
        let factors = finite_field_factor::factor_with_multiplicity(model, field, &reduction)?;
        if matches!(factors.as_slice(), [(d, 1, _)] if *d == degree) {
            return Ok(FactorKind::Irreducible);
        }
    }
    Ok(FactorKind::Uncertified)
}
