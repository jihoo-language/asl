use std::fmt::Write;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_qp_quad::Wide;
use crate::rv_unram::Ring;
use crate::syntax::ParsedExpr;

const LIMBS: usize = 40;
type Integer = Wide<LIMBS>;

pub(crate) enum Statement {
    Degree(usize),
    Coefficient(usize, usize, usize),
}

pub(crate) enum Expression {
    Evaluate(usize, Box<[ParsedExpr; 2]>),
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
        let ell = level.get() as usize;
        let r = degree(level);
        let m = ell.saturating_mul(r);
        // Construct qj-powers through m, then solve r triangular systems
        // and take r Newton steps. Charge for the cubic work before allocation.
        let work = (m as u128)
            .saturating_mul(m.saturating_add(1) as u128)
            .saturating_mul(m.saturating_add(2) as u128)
            / 6;
        let work = work
            .saturating_add(8u128.saturating_mul(m as u128).saturating_mul(m as u128))
            .saturating_mul((LIMBS * 2) as u128);
        let work = usize::try_from(work).map_err(|_| computation_limit())?;
        model.spend(work).map_err(|_| computation_limit())?;
        model.allocate((r + 1) * (r + 1) * (LIMBS + 1))?;

        let j = shifted_j(m)?;
        let powers = shifted_powers(&j)?;
        let mut traces = vec![Vec::new()];
        for k in 1..=r {
            let mut series = vec![Integer::from_i128(0); ell * k + 1];
            for d in (0..series.len()).step_by(ell) {
                series[d] = powers[k][k - d / ell];
            }
            for d in 0..=k / ell {
                series[d] = series[d].add(powers[k][k - ell * d].multiply(ell as i128)?)?;
            }
            // The distinct leading poles of 1,j,j^2,... give a unit
            // triangular basis for the nonpositive Laurent coefficients.
            let mut trace = vec![Integer::from_i128(0); series.len()];
            for d in (0..series.len()).rev() {
                let coefficient = series[d];
                trace[d] = coefficient;
                if !coefficient.is_zero() {
                    for k in 1..=d {
                        series[d - k] =
                            series[d - k].add(coefficient.multiply_wide(powers[d][k])?.negate())?;
                    }
                }
            }
            traces.push(trace);
        }

        let mut elementary = vec![vec![Integer::from_i128(1)]];
        for k in 1..=r {
            let mut next = vec![Integer::from_i128(0); ell * k + 1];
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
        let p = field.prime();
        for coefficient in row.iter().rev() {
            result = field.multiply(model, &result, x)?;
            if !coefficient.is_zero() {
                let (_, remainder) = coefficient.magnitude.divide(model, u128::from(p))?;
                let remainder = remainder as u64;
                let residue = if coefficient.negative && remainder != 0 {
                    p - remainder
                } else {
                    remainder
                };
                result[0] = crate::rv_unram::add_mod(result[0], residue, p);
            }
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

fn shifted_powers(j: &[Integer]) -> Result<Vec<Vec<Integer>>> {
    let mut logarithm = vec![Integer::from_i128(0); j.len()];
    for n in 1..j.len() {
        let mut coefficient = j[n].multiply(n as i128)?;
        for k in 1..n {
            coefficient = coefficient.add(logarithm[k].multiply_wide(j[n - k])?.negate())?;
        }
        logarithm[n] = coefficient;
    }
    let mut powers = vec![vec![Integer::from_i128(1)]];
    for d in 1..j.len() {
        let mut power = vec![Integer::from_i128(1)];
        for n in 1..=d {
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

fn divide_small(value: Integer, divisor: usize) -> Result<Integer> {
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
