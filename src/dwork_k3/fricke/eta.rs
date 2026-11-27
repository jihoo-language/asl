use crate::error::Result;
use crate::lee_linear::Rational;
use crate::model::Model;

use super::super::mirror::series::{
    compose, degree_ok, display_row, inverse_mirror, product, quotient_by_unit, theta, Series,
};
use super::{arithmetic, certificate};

const ONE: Rational = Rational::integer(1);
const ZERO: Rational = Rational::ZERO;

pub(super) fn report(model: &mut Model, degree: usize) -> Result<String> {
    degree_ok(degree)?;
    let (f, inverse) = inverse_mirror(model, degree)?;
    let t = eta_product(model, degree)?;
    let mut plus = t.clone();
    let mut minus = t.clone();
    if degree > 0 {
        plus[1] = plus[1]
            .add(Rational::integer(64))
            .map_err(|_| arithmetic())?;
        minus[1] = minus[1]
            .add(Rational::integer(-64))
            .map_err(|_| arithmetic())?;
    }
    let denominator = product(model, &plus, &plus)?;
    let z_over_q = quotient_by_unit(model, &t, &denominator)?;
    let z = std::iter::once(ZERO)
        .chain(z_over_q.iter().copied())
        .collect::<Vec<_>>();
    if z != inverse {
        return Err(certificate());
    }
    let period = compose(model, &f, &z[..=degree])?;
    let mut log_derivative = quotient_by_unit(model, &theta(&t)?, &t)?;
    for entry in &mut log_derivative {
        entry.numerator = entry.numerator.checked_neg().ok_or_else(arithmetic)?;
    }
    log_derivative[0] = ONE;

    let eisenstein = eisenstein_difference(degree);
    let square_root = quotient_by_unit(model, &minus, &plus)?;
    let mut z_derivative = quotient_by_unit(model, &theta(&z_over_q)?, &z_over_q)?;
    z_derivative[0] = ONE;
    if period != log_derivative
        || period != eisenstein
        || product(model, &period, &square_root)? != z_derivative
    {
        return Err(certificate());
    }
    Ok(format!(
        "rvDworkFrickeSeries[{degree}] = q*t={}, z(q)={}, f(z(q))={}, certified inverse mirror map, Eisenstein difference and Wronskian normalization",
        display_row(&t), display_row(&z), display_row(&period),
    ))
}

fn eta_product(model: &mut Model, degree: usize) -> Result<Series> {
    let mut product_so_far = vec![ZERO; degree + 1];
    product_so_far[0] = ONE;
    for n in 1..=degree {
        let mut factor = vec![ZERO; degree + 1];
        factor[0] = ONE;
        let mut coefficient = 1i128;
        for j in 1..=degree / n {
            // (-1)^j binom(j+23,j), with the exact division checked.
            let term = coefficient * -(j as i128 + 23);
            if term % j as i128 != 0 {
                return Err(certificate());
            }
            coefficient = term / j as i128;
            factor[n * j] = Rational::integer(coefficient);
        }
        product_so_far = product(model, &product_so_far, &factor)?;
    }
    Ok(product_so_far)
}

pub(super) fn eisenstein_difference(degree: usize) -> Series {
    let mut sums = vec![0i128; degree + 1];
    for d in 1..=degree {
        for n in (d..=degree).step_by(d) {
            sums[n] += d as i128;
        }
    }
    let coefficient = |n: usize| {
        if n.is_multiple_of(2) {
            sums[n] - 2 * sums[n / 2]
        } else {
            sums[n]
        }
    };
    std::iter::once(ONE)
        .chain((1..=degree).map(|n| Rational::integer(24 * coefficient(n))))
        .collect()
}
