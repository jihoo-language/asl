use std::fmt::Write;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_qp_quad::Wide;
use crate::rv_unram::Ring;
use crate::syntax::ParsedExpr;

const LIMBS: usize = 40;
pub(crate) type Integer = Wide<LIMBS>;

pub(crate) enum Statement {
    Degree(usize),
    Coefficient(usize, usize, usize),
}

pub(crate) enum Expression {
    Evaluate(usize, Box<[ParsedExpr; 2]>),
    WittEvaluate(usize, Box<[ParsedExpr; 2]>),
    SpecializeCoefficient(usize, usize, Box<ParsedExpr>),
}

pub(crate) struct Polynomial {
    // The row indexed by j contains the coefficients of X^i Y^j.
    rows: Vec<Vec<Integer>>,
}

pub(crate) fn prime(model: &mut Model, level: usize) -> Result<Prime> {
    if level == 2 {
        return Ok(Prime::validated(2));
    }
    if level < 2 || level.is_multiple_of(2) {
        return Err(Error::plain(
            "classical modular polynomial level must be prime",
        ));
    }
    Prime::new(model, level)
}

pub(crate) fn degree(level: Prime) -> usize {
    level.get() as usize + 1
}

impl Polynomial {
    pub(crate) fn new(model: &mut Model, level: Prime) -> Result<Self> {
        Self::from_hauptmodul(model, level, 0, |_, cutoff| shifted_j(cutoff))
    }

    pub(crate) fn from_hauptmodul(
        model: &mut Model,
        level: Prime,
        audit: usize,
        shifted: impl FnOnce(&mut Model, usize) -> Result<Vec<Integer>>,
    ) -> Result<Self> {
        Self::reconstruct(
            model,
            level.get() as usize,
            degree(level),
            None,
            audit,
            shifted,
        )
    }

    pub(crate) fn from_cyclic_square(
        model: &mut Model,
        prime: Prime,
        audit: usize,
        shifted: impl FnOnce(&mut Model, usize) -> Result<Vec<Integer>>,
    ) -> Result<Self> {
        let ell = prime.get() as usize;
        Self::reconstruct(model, ell * ell, ell * (ell + 1), Some(ell), audit, shifted)
    }

    fn reconstruct(
        model: &mut Model,
        step: usize,
        r: usize,
        middle: Option<usize>,
        audit: usize,
        shifted: impl FnOnce(&mut Model, usize) -> Result<Vec<Integer>>,
    ) -> Result<Self> {
        let m = step.saturating_mul(r);
        // The trace of the k-th power has cusp pole at most ell*k.
        // The extra widths admit the positive trace coefficients and
        // direct substitution on the separate F(q^ell) branch.
        let low_extra = if audit == 0 {
            0
        } else {
            (step * audit).max(m + audit)
        };
        let width = (m + audit).max(r + low_extra);
        let cubic = (m as u128)
            .saturating_mul(m.saturating_add(1) as u128)
            .saturating_mul(m.saturating_add(2) as u128)
            / 6;
        let extra = (m as u128)
            .saturating_mul(audit as u128)
            .saturating_mul(width as u128)
            .saturating_add(
                (r as u128)
                    .saturating_mul(low_extra as u128)
                    .saturating_mul(width as u128),
            );
        let work = cubic
            .saturating_add(extra)
            .saturating_add(
                8u128
                    .saturating_mul(width as u128)
                    .saturating_mul(width as u128),
            )
            // Primitive prime-square levels have a separate bounded
            // width envelope and a separate wide arithmetic budget.
            .saturating_mul(if middle.is_some() {
                24
            } else {
                (LIMBS * 2) as u128
            });
        let work = usize::try_from(work).map_err(|_| computation_limit())?;
        model.spend(work).map_err(|_| computation_limit())?;
        model.allocate((r + 1) * (r + 1) * (LIMBS + 1))?;

        let series = shifted(model, width)?;
        if series.len() != width + 1 || series[0].signed() != Some(1) {
            return Err(Error::plain(
                "modular Hauptmodul must have simple pole with leading coefficient one",
            ));
        }
        let precision = (1..=m).map(|d| d + if d <= r { low_extra } else { audit });
        let powers = shifted_powers(&series, precision)?;
        let mut traces = vec![Vec::new()];
        for k in 1..=r {
            traces.push(trace_polynomial(step, k, middle, audit, &powers)?);
        }

        let mut elementary = vec![vec![Integer::from_i128(1)]];
        for k in 1..=r {
            let mut next = vec![Integer::from_i128(0); step * k + 1];
            for i in 1..=k {
                add_convolution(&mut next, &elementary[k - i], &traces[i], i % 2 == 0)?;
            }
            if next.iter().skip(r + 1).any(|c| !c.is_zero()) {
                return Err(Error::plain(
                    "classical modular polynomial has an inconsistent bidegree",
                ));
            }
            next.resize(r + 1, Integer::from_i128(0));
            for coefficient in &mut next {
                *coefficient = divide_small(*coefficient, k)?;
            }
            elementary.push(next);
        }

        let mut rows = Vec::with_capacity(r + 1);
        for (k, mut row) in elementary.into_iter().enumerate().rev() {
            if k % 2 == 1 {
                for coefficient in &mut row {
                    *coefficient = coefficient.negate();
                }
            }
            row.resize(r + 1, Integer::from_i128(0));
            rows.push(row);
        }
        for (i, row) in rows.iter().enumerate() {
            for (j, other) in rows[..i].iter().enumerate() {
                let a = row[j];
                let b = other[i];
                if a.magnitude != b.magnitude || (!a.is_zero() && a.negative != b.negative) {
                    return Err(Error::plain(
                        "classical modular polynomial failed the symmetry check",
                    ));
                }
            }
        }
        if audit > 0 {
            check_neighbor_branch(model, step, audit, &rows, &powers)?;
        }
        Ok(Self { rows })
    }

    pub(crate) fn coefficient_text(&self, i: usize, j: usize) -> String {
        self.rows
            .get(j)
            .and_then(|row| row.get(i))
            .copied()
            .map_or_else(|| "0".into(), decimal)
    }

    pub(crate) fn wide_coefficient<const N: usize>(&self, i: usize, j: usize) -> Result<Wide<N>> {
        self.rows
            .get(j)
            .and_then(|row| row.get(i))
            .copied()
            .unwrap_or_else(|| Integer::from_i128(0))
            .resize()
    }

    pub(crate) fn field_coefficient(
        &self,
        model: &mut Model,
        field: &Ring,
        i: usize,
        j: usize,
    ) -> Result<Vec<u64>> {
        let coefficient = self.wide_coefficient::<LIMBS>(i, j)?;
        Ok(field.scalar(residue(model, coefficient, field.modulus())?))
    }

    pub(crate) fn specialization(
        &self,
        model: &mut Model,
        field: &Ring,
        x: &[u64],
        j: usize,
    ) -> Result<Vec<u64>> {
        let Some(row) = self.rows.get(j) else {
            return Ok(field.scalar(0));
        };
        let mut result = field.scalar(0);
        let p = field.modulus();
        for coefficient in row.iter().rev() {
            result = field.multiply(model, &result, x)?;
            let value = residue(model, *coefficient, p)?;
            result[0] = crate::rv_unram::add_mod(result[0], value, p);
        }
        Ok(result)
    }

    pub(crate) fn evaluate(
        &self,
        model: &mut Model,
        field: &Ring,
        x: &[u64],
        y: &[u64],
    ) -> Result<Vec<u64>> {
        let mut result = field.scalar(0);
        for j in (0..self.rows.len()).rev() {
            result = field.multiply(model, &result, y)?;
            let row = self.specialization(model, field, x, j)?;
            result = field.add(model, &result, &row)?;
        }
        Ok(result)
    }
}

fn residue(model: &mut Model, coefficient: Integer, modulus: u64) -> Result<u64> {
    if coefficient.is_zero() {
        return Ok(0);
    }
    let (_, remainder) = coefficient.magnitude.divide(model, u128::from(modulus))?;
    let remainder = remainder as u64;
    Ok(if coefficient.negative && remainder != 0 {
        modulus - remainder
    } else {
        remainder
    })
}

pub(crate) fn special_char_evaluate(
    model: &mut Model,
    field: &Ring,
    x: &[u64],
    y: &[u64],
) -> Result<Vec<u64>> {
    let x_ell = field.frobenius(model, x, 1)?;
    let y_ell = field.frobenius(model, y, 1)?;
    let minus_y = field.negate(model, y)?;
    let minus_y_ell = field.negate(model, &y_ell)?;
    let left = field.add(model, &x_ell, &minus_y)?;
    let right = field.add(model, x, &minus_y_ell)?;
    field.multiply(model, &left, &right)
}

pub(crate) fn special_char_coefficient(
    model: &mut Model,
    field: &Ring,
    x: &[u64],
    j: usize,
) -> Result<Vec<u64>> {
    let ell = field.prime() as usize;
    match j {
        0 => {
            let x_ell = field.frobenius(model, x, 1)?;
            field.multiply(model, x, &x_ell)
        }
        1 => field.negate(model, x),
        j if j == ell => {
            let x_ell = field.frobenius(model, x, 1)?;
            field.negate(model, &x_ell)
        }
        j if j == ell + 1 => Ok(field.scalar(1)),
        _ => Ok(field.scalar(0)),
    }
}

fn trace_polynomial(
    ell: usize,
    k: usize,
    middle: Option<usize>,
    audit: usize,
    powers: &[Vec<Integer>],
) -> Result<Vec<Integer>> {
    let pole = ell * k;
    let mut expansion = vec![Integer::from_i128(0); pole + audit + 1];
    for (i, &coefficient) in powers[k].iter().enumerate() {
        let exponent = i as i64 - k as i64;
        let stretched = exponent * ell as i64;
        if stretched <= audit as i64 {
            let n = (pole as i64 + stretched) as usize;
            expansion[n] = expansion[n].add(coefficient)?;
        }
        if exponent % ell as i64 == 0 {
            let contracted = exponent / ell as i64;
            if contracted <= audit as i64 {
                let n = (pole as i64 + contracted) as usize;
                expansion[n] = expansion[n].add(coefficient.multiply(ell as i128)?)?;
            }
        }
        if let Some(p) = middle {
            if exponent <= audit as i64 {
                // The middle cyclic cosets are q times the nontrivial
                // p-th roots of unity, so the sum is c_p(exponent).
                let ramanujan = if exponent % p as i64 == 0 {
                    p as i128 - 1
                } else {
                    -1
                };
                let n = (pole as i64 + exponent) as usize;
                expansion[n] = expansion[n].add(coefficient.multiply(ramanujan)?)?;
            }
        }
    }
    // 1,F,...,F^pole form a monic triangular basis for principal parts.
    let mut trace = vec![Integer::from_i128(0); pole + 1];
    for d in (0..=pole).rev() {
        let coefficient = expansion[pole - d];
        trace[d] = coefficient;
        if !coefficient.is_zero() {
            for (output, &b) in expansion[pole - d..].iter_mut().zip(&powers[d]) {
                *output = output.add(coefficient.multiply_wide(b)?.negate())?;
            }
        }
    }
    if expansion.iter().any(|c| !c.is_zero()) {
        return Err(Error::plain(
            "modular Hauptmodul trace failed the extra q-coefficient check",
        ));
    }
    Ok(trace)
}

// The fully reconstructed polynomial is checked independently on the
// neighbor F(q^ell), with all terms of exponent <= audit retained.
fn check_neighbor_branch(
    model: &mut Model,
    ell: usize,
    audit: usize,
    rows: &[Vec<Integer>],
    powers: &[Vec<Integer>],
) -> Result<()> {
    let pole = (rows.len() - 1) * (ell + 1);
    let mut residual = vec![Integer::from_i128(0); pole + audit + 1];
    for (j, row) in rows.iter().enumerate() {
        for (i, &coefficient) in row.iter().enumerate() {
            if coefficient.is_zero() {
                continue;
            }
            let offset = pole - i - ell * j;
            for (b, &y) in powers[j].iter().enumerate() {
                let start = offset + ell * b;
                if start >= residual.len() {
                    break;
                }
                if y.is_zero() {
                    continue;
                }
                let cy = coefficient.multiply_wide(y)?;
                let count = (residual.len() - start).min(powers[i].len());
                model.spend(count)?;
                for (out, &x) in residual[start..start + count].iter_mut().zip(&powers[i]) {
                    *out = out.add(cy.multiply_wide(x)?)?;
                }
            }
        }
    }
    if residual.iter().any(|c| !c.is_zero()) {
        return Err(Error::plain(
            "modular polynomial failed direct q-branch substitution",
        ));
    }
    Ok(())
}

fn shifted_j(m: usize) -> Result<Vec<Integer>> {
    let mut sigma = vec![0i128; m + 1];
    let mut e4 = vec![Integer::from_i128(0); m + 1];
    e4[0] = Integer::from_i128(1);
    for divisor in 1..=m {
        let d = divisor as i128;
        let cube = Integer::from_i128(240 * d * d * d);
        for index in (divisor..=m).step_by(divisor) {
            e4[index] = e4[index].add(cube)?;
            sigma[index] += d;
        }
    }
    let mut square = vec![Integer::from_i128(0); m + 1];
    add_convolution(&mut square, &e4, &e4, false)?;
    let mut cube = vec![Integer::from_i128(0); m + 1];
    add_convolution(&mut cube, &square, &e4, false)?;

    // qj = E4^3 P, where P = product (1-q^n)^(-24).
    // Its logarithmic derivative qP'/P has coefficients 24 sigma_1(n).
    let mut reciprocal = vec![Integer::from_i128(1)];
    for n in 1..=m {
        let mut coefficient = Integer::from_i128(0);
        for i in 1..=n {
            coefficient = coefficient.add(reciprocal[n - i].multiply(24 * sigma[i])?)?;
        }
        reciprocal.push(divide_small(coefficient, n)?);
    }
    let mut j = vec![Integer::from_i128(0); m + 1];
    add_convolution(&mut j, &cube, &reciprocal, false)?;
    Ok(j)
}

fn shifted_powers(
    j: &[Integer],
    precision: impl Iterator<Item = usize>,
) -> Result<Vec<Vec<Integer>>> {
    let mut logarithm = vec![Integer::from_i128(0); j.len()];
    for n in 1..j.len() {
        let mut coefficient = j[n].multiply(n as i128)?;
        for k in 1..n {
            coefficient = coefficient.add(logarithm[k].multiply_wide(j[n - k])?.negate())?;
        }
        logarithm[n] = coefficient;
    }
    let mut powers = vec![vec![Integer::from_i128(1)]];
    for (index, bound) in precision.enumerate() {
        let d = index + 1;
        let mut power = vec![Integer::from_i128(1)];
        for n in 1..=bound {
            let mut coefficient = Integer::from_i128(0);
            for k in 1..=n {
                coefficient = coefficient.add(logarithm[k].multiply_wide(power[n - k])?)?;
            }
            power.push(divide_small(coefficient.multiply(d as i128)?, n)?);
        }
        powers.push(power);
    }
    Ok(powers)
}

pub(crate) fn divide_small(value: Integer, divisor: usize) -> Result<Integer> {
    let (magnitude, remainder) = value.magnitude.divide_u32(divisor as u32);
    if remainder != 0 {
        return Err(Error::plain(
            "classical modular polynomial recurrence failed exact division",
        ));
    }
    Ok(Integer {
        magnitude,
        negative: value.negative,
    })
}

fn add_convolution(
    output: &mut [Integer],
    a: &[Integer],
    b: &[Integer],
    negative: bool,
) -> Result<()> {
    for (i, &a) in a.iter().take(output.len()).enumerate() {
        if a.is_zero() {
            continue;
        }
        for (slot, &b) in output[i..].iter_mut().zip(b) {
            if !b.is_zero() {
                let term = a.multiply_wide(b)?;
                *slot = slot.add(if negative { term.negate() } else { term })?;
            }
        }
    }
    Ok(())
}

pub(crate) fn decimal<const N: usize>(value: Wide<N>) -> String {
    if value.is_zero() {
        return "0".into();
    }
    let mut magnitude = value.magnitude;
    let mut digits = Vec::new();
    while !magnitude.is_zero() {
        let (quotient, remainder) = magnitude.divide_u32(1_000_000_000);
        digits.push(remainder);
        magnitude = quotient;
    }
    let sign = if value.negative { "-" } else { "" };
    let mut result = format!(
        "{sign}{}",
        digits.pop().expect("nonzero modular coefficient")
    );
    for digit in digits.into_iter().rev() {
        write!(result, "{digit:09}").expect("writing to a string");
    }
    result
}

fn computation_limit() -> Error {
    Error::plain("exact classical modular polynomial exceeds the computation limit")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_exact_coefficients_satisfy_kronecker_congruence() -> std::result::Result<(), String> {
        verify_kronecker().map_err(|error| error.to_string())
    }

    fn verify_kronecker() -> Result<()> {
        let mut model = Model::new();
        for ell in [2, 3, 5, 7] {
            model.reset_work();
            let polynomial = Polynomial::new(&mut model, Prime::validated(ell))?;
            let e = ell as usize;
            for (j, row) in polynomial.rows.iter().enumerate() {
                for (i, value) in row.iter().enumerate() {
                    let residue = i64::from(value.magnitude.residue(ell));
                    let residue = if value.negative { -residue } else { residue };
                    let expected = if (i, j) == (1, 1) || (i, j) == (e, e) {
                        i64::from(ell - 1)
                    } else {
                        i64::from((i, j) == (e + 1, 0) || (i, j) == (0, e + 1))
                    };
                    assert_eq!(
                        residue.rem_euclid(i64::from(ell)),
                        expected,
                        "level {ell}, i={i}, j={j}"
                    );
                }
            }
        }
        Ok(())
    }
}
