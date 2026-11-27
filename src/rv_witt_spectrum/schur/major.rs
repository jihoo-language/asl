use crate::error::{Error, Result};
use crate::model::Model;

use super::super::{add, multiply, one, zero, Integer};
use super::Shape;

pub(crate) enum Branch {
    Major(Shape),
    Cyclic(Shape, u128),
    Lie(Shape),
}

impl Branch {
    pub(crate) fn report(self, model: &mut Model) -> Result<String> {
        match self {
            Self::Major(shape) => {
                let counts = residues(model, &shape)?
                    .into_iter()
                    .map(magnitude)
                    .collect::<Result<Vec<_>>>()?;
                Ok(format!(
                    "rvWittSchurMajorIndex[{:?}] = {counts:?}",
                    shape.rows()
                ))
            }
            Self::Cyclic(shape, d) => {
                let n = shape.weight;
                if d == 0 || n == 0 || !n.is_multiple_of(d) {
                    return Err(Error::plain("cyclotomic order must be a positive divisor of a nonempty Schur shape's weight"));
                }
                let single = multiplicity(model, &shape, (n / d) % n)?;
                let rational = if single == 0 {
                    0
                } else {
                    let order = usize::try_from(d).map_err(|_| overflow())?;
                    let phi = model.cyclotomic(order)?.len() - 1;
                    single.checked_mul(phi as u128).ok_or_else(overflow)?
                };
                Ok(format!("rvWittCyclicSchurMultiplicity[{:?}, {d}] = single={single}, rational={rational}", shape.rows()))
            }
            Self::Lie(shape) => {
                if shape.weight == 0 {
                    return Err(Error::plain(
                        "Lie Schur multiplicity requires a nonempty shape",
                    ));
                }
                let single = multiplicity(model, &shape, 1 % shape.weight)?;
                Ok(format!(
                    "rvWittLieSchurMultiplicity[{:?}] = {single}",
                    shape.rows()
                ))
            }
        }
    }
}

fn multiplicity(model: &mut Model, shape: &Shape, index: u128) -> Result<u128> {
    if let Some(residue) = shape.elementary_major() {
        return Ok(u128::from(residue == index));
    }
    magnitude(residues(model, shape)?[index as usize])
}

fn residues(model: &mut Model, shape: &Shape) -> Result<Vec<Integer>> {
    let n = shape.bounded_weight()?.max(1);
    if let Some(index) = shape.elementary_major() {
        model.spend(n)?;
        let mut counts = vec![zero(); n];
        counts[index as usize] = one();
        return Ok(counts);
    }
    let mut heights = vec![0; shape.rows[0] as usize];
    for &row in shape.rows() {
        for height in &mut heights[..row as usize] {
            *height += 1;
        }
    }
    let mut hooks = Vec::with_capacity(n);
    let mut shift = 0;
    for (i, &row) in shape.rows().iter().enumerate() {
        let row = row as usize;
        shift += i * row;
        hooks.extend(
            heights[..row]
                .iter()
                .enumerate()
                .map(|(j, &height)| row - j + height - i - 1),
        );
    }
    let mut values = vec![zero(); n];
    values[shift % n] = one();
    // Cancel q-integers in Z[q] by cyclotomic valuation before reducing
    // modulo q^n-1. In particular no division occurs at a root of unity.
    for order in 2..=n {
        model.spend(n)?;
        let count = hooks
            .iter()
            .filter(|&&hook| hook.is_multiple_of(order))
            .count();
        let exponent = (n / order)
            .checked_sub(count)
            .ok_or_else(|| Error::plain("Schur q-hook cyclotomic division is inexact"))?;
        if exponent != 0 {
            let factor = model.cyclotomic(order)?;
            let terms: Vec<_> = factor
                .iter()
                .enumerate()
                .filter(|(_, &c)| c != 0)
                .map(|(degree, &c)| (degree, Integer::from_i128(c)))
                .collect();
            for _ in 0..exponent {
                model.spend(n)?;
                let mut product = vec![zero(); n];
                for (i, &value) in values.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
                    for &(j, coefficient) in &terms {
                        let target = &mut product[(i + j) % n];
                        *target = add(*target, multiply(model, value, coefficient)?)?;
                    }
                }
                values = product;
            }
        }
    }
    if values.iter().any(|v| v.negative && !v.is_zero()) {
        return Err(Error::plain(
            "Schur major-index multiplicity became negative",
        ));
    }
    Ok(values)
}

fn magnitude(value: Integer) -> Result<u128> {
    value.magnitude.small().ok_or_else(overflow)
}

fn overflow() -> Error {
    Error::plain("Schur major-index multiplicity exceeds unsigned 128-bit range")
}
