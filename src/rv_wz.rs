use std::collections::{btree_map::Entry, BTreeMap};

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_ore::{self, Powers};

// Keep the same matrix/row storage order of magnitude for each sparse
// polynomial. The verification holds a bounded number of these at once.
const MAX_TERMS: usize = MAX_PIVOT_CELLS / 8;

pub(crate) struct Certificate {
    pub(crate) u: usize,
    pub(crate) v: usize,
    pub(crate) z: i128,
    pub(crate) start: usize,
    pub(crate) numerator: Vec<(usize, usize, i128)>,
    pub(crate) denominator: Vec<Factor>,
}

pub(crate) struct Factor {
    pub(crate) n: i128,
    pub(crate) k: i128,
    pub(crate) constant: i128,
    pub(crate) multiplicity: usize,
}

#[derive(Default)]
struct Polynomial(BTreeMap<(usize, usize), i128>);

impl Polynomial {
    fn one() -> Self {
        Self(BTreeMap::from([((0, 0), 1)]))
    }

    fn affine(n: i128, k: i128, constant: i128) -> Self {
        Self(
            [((1, 0), n), ((0, 1), k), ((0, 0), constant)]
                .into_iter()
                .filter(|&(_, c)| c != 0)
                .collect(),
        )
    }

    fn coefficient(
        &mut self,
        model: &mut Model,
        powers: (usize, usize),
        value: i128,
        subtract: bool,
    ) -> Result<()> {
        let size = self.0.len();
        model.spend(1 + size.max(1).ilog2() as usize)?;
        if value == 0 {
            return Ok(());
        }
        match self.0.entry(powers) {
            Entry::Vacant(entry) => {
                if size == MAX_TERMS {
                    return Err(Error::plain("WZ polynomial exceeds the term storage limit"));
                }
                let value = if subtract {
                    value.checked_neg().ok_or_else(arithmetic)?
                } else {
                    value
                };
                entry.insert(value);
            }
            Entry::Occupied(mut entry) => {
                let old = *entry.get();
                let sum = if subtract {
                    old.checked_sub(value)
                } else {
                    old.checked_add(value)
                }
                .ok_or_else(arithmetic)?;
                if sum == 0 {
                    entry.remove();
                } else {
                    *entry.get_mut() = sum;
                }
            }
        }
        Ok(())
    }

    fn accumulate(&mut self, model: &mut Model, source: &Self, subtract: bool) -> Result<()> {
        for (&powers, &value) in &source.0 {
            self.coefficient(model, powers, value, subtract)?;
        }
        Ok(())
    }

    fn scale(mut self, model: &mut Model, value: i128) -> Result<Self> {
        model.spend(self.0.len())?;
        if value == 0 {
            return Ok(Self::default());
        }
        for c in self.0.values_mut() {
            *c = c.checked_mul(value).ok_or_else(arithmetic)?;
        }
        Ok(self)
    }

    fn product(&self, model: &mut Model, right: &Self) -> Result<Self> {
        let mut product = Self::default();
        for (&(i, j), &a) in &self.0 {
            for (&(h, l), &b) in &right.0 {
                let n = i.checked_add(h).ok_or_else(degree_error)?;
                let k = j.checked_add(l).ok_or_else(degree_error)?;
                let value = a.checked_mul(b).ok_or_else(arithmetic)?;
                product.coefficient(model, (n, k), value, false)?;
            }
        }
        Ok(product)
    }

    fn power(mut self, model: &mut Model, mut exponent: usize) -> Result<Self> {
        let mut result = Self::one();
        while exponent != 0 {
            if !exponent.is_multiple_of(2) {
                result = result.product(model, &self)?;
            }
            exponent /= 2;
            if exponent != 0 {
                self = self.product(model, &self)?;
            }
        }
        Ok(result)
    }

    fn shift_k(&self, model: &mut Model) -> Result<Self> {
        let mut result = Self::default();
        for (&(n, k), &value) in &self.0 {
            model.spend(k)?;
            let mut choose = 1;
            for i in 0..=k {
                let contribution = value.checked_mul(choose).ok_or_else(arithmetic)?;
                result.coefficient(model, (n, i), contribution, false)?;
                if i != k {
                    choose = binomial_step(choose, (k - i) as i128, (i + 1) as i128)?;
                }
            }
        }
        Ok(result)
    }
}

pub(crate) fn verify(
    model: &mut Model,
    operator: &[(Powers, i128)],
    certificate: &Certificate,
) -> Result<()> {
    if certificate.u == 0 || certificate.z == 0 {
        return Err(Error::plain("WZ summand requires u >= 1 and nonzero z"));
    }
    let max_shift = operator.last().map_or(0, |(powers, _)| powers.q);
    if certificate.start == 0 || certificate.start < max_shift {
        return Err(Error::plain(format!(
            "WZ start must be positive and at least the largest q-shift {max_shift}"
        )));
    }
    let vertical_degree = certificate
        .u
        .checked_add(certificate.v)
        .ok_or_else(degree_error)?;
    model.spend(max_shift)?;
    let (denominator, shifted_denominator) = denominators(model, certificate)?;
    let mut numerator = Polynomial::default();
    for &(n, k, value) in &certificate.numerator {
        numerator.coefficient(model, (n, k), value, false)?;
    }
    if numerator.0.keys().any(|&(_, k)| k == 0) {
        return Err(Error::plain("WZ certificate must have R(n,0) = 0"));
    }

    let (horizontal_num, horizontal_den) =
        horizontal(model, operator, max_shift, certificate.u, certificate.v)?;
    let rising = Polynomial::affine(1, 1, 1).power(model, certificate.v)?;
    let vertical_num = Polynomial::affine(1, -1, 0)
        .power(model, certificate.u)?
        .product(model, &rising)?
        .scale(model, certificate.z)?;
    let vertical_den = Polynomial::affine(0, 1, 1).power(model, vertical_degree)?;

    // Common denominator: horizontal_den * vertical_den *
    // denominator(n,k) * denominator(n,k+1).
    let scaled_left = vertical_den.product(model, &horizontal_num)?;
    let common_right = shifted_denominator.product(model, &vertical_den)?;
    let mut difference = scaled_left
        .product(model, &denominator)?
        .product(model, &shifted_denominator)?;
    let forward = numerator
        .shift_k(model)?
        .product(model, &vertical_num)?
        .product(model, &horizontal_den)?
        .product(model, &denominator)?;
    difference.accumulate(model, &forward, true)?;
    let back = numerator
        .product(model, &horizontal_den)?
        .product(model, &common_right)?;
    difference.accumulate(model, &back, false)?;
    if let Some((&(n, k), &value)) = difference.0.first_key_value() {
        return Err(Error::plain(format!(
            "WZ cleared bivariate identity fails at n^{n} k^{k} (coefficient {value})"
        )));
    }

    verify_initial(model, operator, certificate)
}

fn denominators(model: &mut Model, certificate: &Certificate) -> Result<(Polynomial, Polynomial)> {
    let n0 = certificate.start as i128;
    let mut denominator = Polynomial::one();
    let mut shifted = Polynomial::one();
    for factor in &certificate.denominator {
        if factor.multiplicity == 0 {
            return Err(Error::plain("WZ denominator multiplicity must be positive"));
        }
        let sum = factor.n.checked_add(factor.k).ok_or_else(arithmetic)?;
        let next_constant = factor
            .constant
            .checked_add(factor.k)
            .ok_or_else(arithmetic)?;
        // Across 0<=k<=n+1, an affine factor attains its minimum at
        // an endpoint. Both endpoint functions are increasing in n.
        if factor.n < 0
            || sum < 0
            || endpoint(factor.n, n0, factor.constant)? <= 0
            || endpoint(sum, n0, next_constant)? <= 0
        {
            return Err(Error::plain(
                "WZ denominator factor must be positive on n>=n0, 0<=k<=n+1",
            ));
        }
        let a = Polynomial::affine(factor.n, factor.k, factor.constant)
            .power(model, factor.multiplicity)?;
        let b = Polynomial::affine(factor.n, factor.k, next_constant)
            .power(model, factor.multiplicity)?;
        denominator = denominator.product(model, &a)?;
        shifted = shifted.product(model, &b)?;
    }
    Ok((denominator, shifted))
}

fn endpoint(slope: i128, n0: i128, constant: i128) -> Result<i128> {
    slope
        .checked_mul(n0)
        .and_then(|value| value.checked_add(constant))
        .ok_or_else(arithmetic)
}

fn horizontal(
    model: &mut Model,
    operator: &[(Powers, i128)],
    max_shift: usize,
    u: usize,
    v: usize,
) -> Result<(Polynomial, Polynomial)> {
    let mut parts = operator.chunk_by(|a, b| a.0.q == b.0.q).peekable();
    let mut numerator = Polynomial::default();
    let mut denominator = Polynomial::one();
    let mut product = Polynomial::one();
    for a in 0..=max_shift {
        if a != 0 {
            let step = 1 - a as i128;
            // Cancel common powers of n-a+1 from the horizontal ratio.
            let excess_v = Polynomial::affine(1, 0, step).power(model, v.saturating_sub(u))?;
            let top = Polynomial::affine(1, -1, step)
                .power(model, u)?
                .product(model, &excess_v)?;
            let rising = Polynomial::affine(1, 1, step).power(model, v)?;
            let bottom = Polynomial::affine(1, 0, step)
                .power(model, u.saturating_sub(v))?
                .product(model, &rising)?;
            product = product.product(model, &top)?;
            denominator = denominator.product(model, &bottom)?;
            numerator = numerator.product(model, &bottom)?;
        }
        if let Some(part) = parts.next_if(|part| part[0].0.q == a) {
            let mut shifted = Polynomial::default();
            for &(powers, c) in part {
                let term = Polynomial::affine(1, 0, -(a as i128))
                    .power(model, powers.theta)?
                    .scale(model, c)?;
                shifted.accumulate(model, &term, false)?;
            }
            let term = shifted.product(model, &product)?;
            numerator.accumulate(model, &term, false)?;
        }
    }
    Ok((numerator, denominator))
}

fn verify_initial(
    model: &mut Model,
    operator: &[(Powers, i128)],
    certificate: &Certificate,
) -> Result<()> {
    let start = certificate.start;
    // The total counts all original hypergeometric summands before n0.
    model.spend(start.saturating_mul(start.saturating_add(1)) / 2)?;
    if start.saturating_mul(2) > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "WZ initial coefficients exceed the storage limit",
        ));
    }
    let mut sums = Vec::with_capacity(start);
    for n in 0..start {
        sums.push(sum(model, certificate, n)?);
        let mut residue = 0i128;
        for &(powers, c) in operator.iter().take_while(|(powers, _)| powers.q <= n) {
            model.spend(1)?;
            let index = n - powers.q;
            if sums[index] != 0 {
                let factor = scalar_power(model, index as i128, powers.theta)?;
                let contribution = c
                    .checked_mul(factor)
                    .and_then(|c| c.checked_mul(sums[index]))
                    .ok_or_else(arithmetic)?;
                residue = residue.checked_add(contribution).ok_or_else(arithmetic)?;
            }
        }
        if residue != 0 {
            return Err(Error::plain(format!(
                "WZ initial coefficient of L(sum) at q^{n} is {residue}"
            )));
        }
    }
    Ok(())
}

fn sum(model: &mut Model, certificate: &Certificate, degree: usize) -> Result<i128> {
    let (mut ordinary, mut rising, mut signed_power) = (1i128, 1i128, 1i128);
    let mut total = 0i128;
    for k in 0..=degree {
        let left = scalar_power(model, ordinary, certificate.u)?;
        let right = scalar_power(model, rising, certificate.v)?;
        let term = left
            .checked_mul(right)
            .and_then(|x| x.checked_mul(signed_power))
            .ok_or_else(arithmetic)?;
        total = total.checked_add(term).ok_or_else(arithmetic)?;
        if k != degree {
            ordinary = binomial_step(ordinary, (degree - k) as i128, (k + 1) as i128)?;
            if certificate.v != 0 {
                rising = binomial_step(rising, degree as i128 + k as i128 + 1, (k + 1) as i128)?;
            }
            signed_power = signed_power
                .checked_mul(certificate.z)
                .ok_or_else(arithmetic)?;
        }
    }
    Ok(total)
}

fn scalar_power(model: &mut Model, mut value: i128, mut exponent: usize) -> Result<i128> {
    let mut result = 1i128;
    while exponent != 0 {
        model.spend(1)?;
        if !exponent.is_multiple_of(2) {
            result = result.checked_mul(value).ok_or_else(arithmetic)?;
        }
        exponent /= 2;
        if exponent != 0 {
            value = value.checked_mul(value).ok_or_else(arithmetic)?;
        }
    }
    Ok(result)
}

fn binomial_step(coefficient: i128, numerator: i128, denominator: i128) -> Result<i128> {
    let common = rv_ore::gcd(numerator as u128, denominator as u128) as i128;
    let remaining = denominator / common;
    if coefficient % remaining != 0 {
        return Err(Error::plain("WZ binomial update is not integral"));
    }
    (coefficient / remaining)
        .checked_mul(numerator / common)
        .ok_or_else(arithmetic)
}

fn arithmetic() -> Error {
    Error::plain("WZ exact integer arithmetic is out of range")
}

fn degree_error() -> Error {
    Error::plain("WZ polynomial degree is out of range")
}
