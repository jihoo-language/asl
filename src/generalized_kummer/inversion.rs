use crate::error::{Error, Result};
use crate::hw_algebra::{self, Matrix};
use crate::hw_closed::Curve;
use crate::hyperelliptic::Residues;
use crate::hyperelliptic_jacobian::Divisor;
use crate::hyperelliptic_kummer_surface::hilbert::{
    checked_product, intermediate_overflow, inversion_lefschetz, unsigned, Integer,
};
use crate::jacobian_level::Basis;
use crate::model::Model;

use super::cohomology::torsion_correction;
use super::{affine_fixed, require_level_basis, require_translation_order};

pub(crate) fn euler(model: &mut Model, n: u128, h: u128) -> Result<u128> {
    require_translation_order(n, h)?;
    // At Q=1 the ungraded inverted Hilbert character has generating series
    // ∏_j ((1+z^j)/(1−z^j))^8. The inverted surface Euler term equals 16.
    let hilbert = inversion_lefschetz(model, n, None)?;
    let n = n as usize;
    let mut numerator = hilbert[n][0];
    if n.is_multiple_of(2) {
        let weight = 16 * i128::from(((n as u128) / h).is_multiple_of(2)) - 1;
        let term = checked_product(model, Integer::from_i128(weight), hilbert[n / 2][0])?;
        numerator = numerator.add(term).map_err(|_| intermediate_overflow())?;
    }
    exact_unsigned_quotient(model, numerator, Integer::from_i128(16))
}

pub(crate) fn points(
    model: &mut Model,
    basis: &Basis,
    curve: &Curve,
    r: u128,
    translation: &Divisor,
) -> Result<(u128, u128)> {
    let coordinate = coordinates(model, basis, curve, translation)?;
    let n = u128::from(basis.prime.get());
    let hilbert = inversion_lefschetz(model, n, Some((curve, r)))?;
    let invariant = quotient(model, hilbert[n as usize][0], hilbert[1][0])?;
    let fixed = fixed_torsion(model, basis, curve, r, &coordinate)?;
    let correction = torsion_correction(model, curve, r, n, fixed)?;
    let total = invariant
        .add(correction)
        .map_err(|_| intermediate_overflow())?;
    let count = unsigned(total).ok_or_else(|| {
        Error::plain("generalized Kummer inversion fixed count exceeds unsigned 128-bit range")
    })?;
    Ok((count, fixed))
}

pub(super) fn coordinates(
    model: &mut Model,
    basis: &Basis,
    base: &Curve,
    translation: &Divisor,
) -> Result<Vec<u32>> {
    require_level_basis(basis)?;
    basis.require_base_extension(model, base)?;
    basis.translation_coordinates(model, translation)
}

pub(super) fn fixed_torsion(
    model: &mut Model,
    basis: &Basis,
    base: &Curve,
    mut r: u128,
    coordinate: &[u32],
) -> Result<u128> {
    let prime = basis.prime;
    let residues = Residues::new(prime, 1)?;
    let frobenius = basis.frobenius(model, base, 1)?.matrix;
    // Homogeneous coordinates encode x ↦ t−F(x). Taking powers also
    // evolves the constant term when Frobenius moves the translating point.
    let mut entries = Vec::with_capacity(25);
    for (row, &component) in coordinate.iter().enumerate() {
        entries
            .extend((0..4).map(|col| u64::from(prime.negate(*frobenius.entry(row, col) as u32))));
        entries.push(u64::from(component));
    }
    entries.extend([0, 0, 0, 0, 1]);
    let mut factor = Matrix::from_entries(5, entries)?;
    let mut matrix = hw_algebra::identity(model, &residues, 5)?;
    while r != 0 {
        if !r.is_multiple_of(2) {
            matrix = hw_algebra::multiply(model, &residues, &matrix, &factor)?;
        }
        r >>= 1;
        if r != 0 {
            factor = hw_algebra::multiply(model, &residues, &factor, &factor)?;
        }
    }
    let translation: [u32; 4] = std::array::from_fn(|row| *matrix.entry(row, 4) as u32);
    affine_fixed(
        model,
        prime,
        |row, col| *matrix.entry(row, col) as u32,
        Some(&translation),
    )
}

fn exact_unsigned_quotient(model: &mut Model, value: Integer, divisor: Integer) -> Result<u128> {
    unsigned(quotient(model, value, divisor)?).ok_or_else(|| {
        Error::plain(
            "generalized Kummer inversion Euler characteristic exceeds unsigned 128-bit range",
        )
    })
}

fn quotient(model: &mut Model, value: Integer, divisor: Integer) -> Result<Integer> {
    value
        .divide_exact(model, divisor)?
        .ok_or_else(|| Error::plain("generalized Kummer inversion Lefschetz division is inexact"))
}
