use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_ore::gcd;

use super::Field;

pub(super) struct Distribution {
    pub(super) framed: Vec<u128>,
    total: u128,
    group: u128,
}

impl Distribution {
    pub(super) fn new(framed: Vec<u128>, total: u128, group: u128) -> Result<Self> {
        let actual = framed
            .iter()
            .try_fold(0u128, |sum, &n| sum.checked_add(n).ok_or_else(bound))?;
        if actual != total {
            return Err(Error::plain(
                "partition quiver trace histogram did not exhaust the ambient matrix space",
            ));
        }
        Ok(Self {
            framed,
            total,
            group,
        })
    }

    pub(super) fn render(&self, name: &str) -> String {
        let stack = self
            .framed
            .iter()
            .map(|&n| {
                let d = gcd(n, self.group);
                format!("{}/{}", n / d, self.group / d)
            })
            .collect::<Vec<_>>()
            .join(", ");
        format!("rvPartitionQuiverDistribution[{name}] = q={}, framed={:?}, total={}, groupOrder={}, stack=[{stack}]", self.framed.len(), self.framed, self.total, self.group)
    }
}

pub(super) fn bound() -> Error {
    Error::plain("partition quiver counts or group orders exceed unsigned 128-bit arithmetic")
}

pub(super) fn power(q: u128, exponent: usize) -> Result<u128> {
    q.checked_pow(u32::try_from(exponent).map_err(|_| bound())?)
        .ok_or_else(bound)
}

pub(super) fn mass(field: &Field, dims: &[usize], entries: usize) -> Result<(u128, u128)> {
    let q = u128::from(field.order());
    if q > 256 {
        return Err(Error::plain(
            "partition quiver trace distributions require a field of order at most 256",
        ));
    }
    let total = power(q, entries)?;
    let group = dims.iter().try_fold(1u128, |product, &d| {
        product.checked_mul(gl(q, d)?).ok_or_else(bound)
    })?;
    Ok((total, group))
}

pub(super) fn gl(q: u128, n: usize) -> Result<u128> {
    let qn = power(q, n)?;
    (0..n).try_fold(1u128, |product, i| {
        product.checked_mul(qn - power(q, i)?).ok_or_else(bound)
    })
}

pub(super) fn unit(field: &Field) -> Vec<u128> {
    let mut result = vec![0; field.order() as usize];
    result[0] = 1;
    result
}

pub(super) fn convolve(
    model: &mut Model,
    field: &Field,
    a: &[u128],
    b: &[u128],
    scale: u128,
) -> Result<Vec<u128>> {
    let mut out = vec![0; field.order() as usize];
    accumulate_convolution(model, field, &mut out, a, b, scale)?;
    Ok(out)
}

pub(super) fn accumulate_convolution(
    model: &mut Model,
    field: &Field,
    out: &mut [u128],
    a: &[u128],
    b: &[u128],
    scale: u128,
) -> Result<()> {
    model.spend(a.len())?;
    for (i, &x) in a.iter().enumerate().filter(|(_, x)| **x != 0) {
        model.spend(b.len())?;
        let scaled = x.checked_mul(scale).ok_or_else(bound)?;
        for (j, &y) in b.iter().enumerate().filter(|(_, y)| **y != 0) {
            let destination = &mut out[field.add(i as u32, j as u32) as usize];
            let contribution = scaled.checked_mul(y).ok_or_else(bound)?;
            *destination = destination.checked_add(contribution).ok_or_else(bound)?;
        }
    }
    Ok(())
}

fn rank_count(q: u128, r: usize, s: usize, k: usize) -> Result<u128> {
    let (qr, qs) = (power(q, r)?, power(q, s)?);
    let numerator = (0..k).try_fold(1u128, |product, i| {
        let qi = power(q, i)?;
        product
            .checked_mul(qr - qi)
            .and_then(|t| t.checked_mul(qs - qi))
            .ok_or_else(bound)
    })?;
    exact_divide(numerator, gl(q, k)?)
}

pub(super) fn exact_divide(numerator: u128, denominator: u128) -> Result<u128> {
    if !numerator.is_multiple_of(denominator) {
        return Err(Error::plain(
            "partition quiver integral orbit count failed exact division",
        ));
    }
    Ok(numerator / denominator)
}

pub(super) fn two_cycle(
    model: &mut Model,
    field: &Field,
    [r, s]: [usize; 2],
    half: u32,
) -> Result<Vec<u128>> {
    let q = u128::from(field.order());
    let mut total = vec![0; q as usize];
    if r == 0 || s == 0 {
        total[0] = 1;
        return Ok(total);
    }
    let mut squares = total.clone();
    for x in 0..q as u32 {
        let value = field.mul(half, field.mul(x, x));
        squares[value as usize] += 1;
    }
    // In Tr(C²)/2, each off-diagonal pair contributes xy. There are
    // 2q−1 zeros and q−1 pairs for each nonzero target, also over extensions.
    let mut hyperbolic = vec![q - 1; q as usize];
    hyperbolic[0] = 2 * q - 1;
    let mut active = unit(field);
    for k in 0..=r.min(s) {
        if k != 0 {
            active = convolve(model, field, &active, &squares, 1)?;
            for _ in 1..k {
                active = convolve(model, field, &active, &hyperbolic, 1)?;
            }
        }
        let free = power(q, r * s - k * k)?;
        let count = rank_count(q, r, s, k)?
            .checked_mul(free)
            .ok_or_else(bound)?;
        for (total, &multiplicity) in total.iter_mut().zip(&active) {
            *total = total
                .checked_add(count.checked_mul(multiplicity).ok_or_else(bound)?)
                .ok_or_else(bound)?;
        }
    }
    Ok(total)
}
