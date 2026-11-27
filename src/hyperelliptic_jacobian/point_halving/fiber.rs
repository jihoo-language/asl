use std::cmp::Ordering;

use super::{field, from_selection, lift, require_affine_point, Divisor};
use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::finite_field_factor::factor_squarefree;
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::model::Model;

type Pair = [Polynomial; 2];

pub(crate) struct Fiber {
    free_pairs: Vec<Pair>,
    pivot: Option<Pair>,
}

impl Fiber {
    pub(crate) fn new(
        model: &mut Model,
        curve: &Curve,
        a: &[u64],
        b: &[u64],
        operation: &str,
    ) -> Result<Option<Self>> {
        require_affine_point(model, curve, a, b, operation)?;
        let field = field(curve);
        let reflection = [a.to_vec(), field.signed(-1)];
        let branch = lift(model, field, curve.polynomial(), &reflection)?;
        let branch = poly::scale(model, field, &branch, &reflection[1])?;
        let ramified = poly::zero(b);
        let branch = if ramified {
            poly::exact_divide(model, field, &branch, &field.polynomial(&[0, 0, 1]))?
        } else {
            branch
        };
        let degree = branch.len() - 1;
        let mut pairs = Vec::new();
        // The factors are permuted by z ↦ −z. A self-paired factor prevents a
        // rational choice of one square root at each of its branch points.
        for factor in factor_squarefree(model, field, branch)? {
            let opposite = opposite(model, field, &factor)?;
            match factor.cmp(&opposite) {
                Ordering::Equal => return Ok(None),
                Ordering::Less => pairs.push([factor, opposite]),
                Ordering::Greater => (),
            }
        }
        if 2 * pairs
            .iter()
            .map(|[factor, _]| factor.len() - 1)
            .sum::<usize>()
            != degree
        {
            return Err(Error::plain(
                "incomplete factor pairs in rational point halving",
            ));
        }
        let widest = pairs
            .iter()
            .map(|[factor, _]| factor.len())
            .max()
            .unwrap_or(0);
        model.spend(
            pairs
                .len()
                .saturating_mul(1 + pairs.len().max(1).ilog2() as usize)
                .saturating_mul(widest * field.cells()),
        )?;
        pairs.sort_unstable();
        let pivot = if ramified {
            None
        } else {
            // One odd-degree factor is necessary because the degrees of the
            // selected factors sum to 2g+1. Flipping it changes R(0)'s sign.
            let index = pairs
                .iter()
                .position(|[factor, _]| factor.len().is_multiple_of(2))
                .ok_or_else(|| {
                    Error::plain("rational point halving has no odd-degree sign pivot")
                })?;
            Some(pairs.remove(index))
        };
        Ok(Some(Self {
            free_pairs: pairs,
            pivot,
        }))
    }

    pub(crate) fn dimension(&self) -> usize {
        self.free_pairs.len()
    }

    pub(crate) fn at(
        self,
        model: &mut Model,
        curve: &Curve,
        a: &[u64],
        b: &[u64],
        index: &str,
    ) -> Result<Divisor> {
        let mut choices = index_bits(model, index, self.dimension())?.into_iter();
        let field = field(curve);
        let mut selected = if self.pivot.is_some() {
            vec![field.one()]
        } else {
            // The ramified branch point contributes z once.
            field.polynomial(&[0, 1])
        };
        for [first, second] in self.free_pairs {
            let factor = if choices.next().unwrap_or(false) {
                second
            } else {
                first
            };
            selected = poly::multiply(model, field, &selected, &factor)?;
        }
        if let Some([first, second]) = self.pivot {
            let constant = field.multiply(model, &selected[0], &first[0])?;
            let factor = if constant == b { first } else { second };
            selected = poly::multiply(model, field, &selected, &factor)?;
        }
        if selected.len() != curve.polynomial().len() || selected[0] != b {
            return Err(Error::plain(
                "rvJacHalfAt obtained inconsistent branch factors",
            ));
        }
        from_selection(model, curve, a, b, selected, "rvJacHalfAt")
    }
}

fn opposite(model: &mut Model, field: Field<'_>, factor: &[Scalar]) -> Result<Polynomial> {
    model.spend(factor.len() * field.cells())?;
    let degree = factor.len() - 1;
    factor
        .iter()
        .enumerate()
        .map(|(exponent, coefficient)| {
            if (degree - exponent).is_multiple_of(2) {
                Ok(coefficient.clone())
            } else {
                field.ring.negate(model, coefficient)
            }
        })
        .collect()
}

fn index_bits(model: &mut Model, index: &str, dimension: usize) -> Result<Vec<bool>> {
    let out_of_range = || Error::plain("rvJacHalfAt index is outside the rational half fiber");
    let index = index.trim_start_matches('0');
    // Every three bits fit in one decimal digit. Reject oversized inputs
    // before converting and check the exact bound while extracting bits.
    if index.len() > dimension.div_ceil(3) {
        return Err(out_of_range());
    }
    let mut digits: Vec<_> = index.bytes().map(|digit| digit - b'0').collect();
    let mut digits = digits.as_mut_slice();
    let mut bits = Vec::new();
    while !digits.is_empty() {
        if bits.len() == dimension {
            return Err(out_of_range());
        }
        model.spend(digits.len())?;
        let mut carry = 0;
        for digit in digits.iter_mut() {
            let value = 10 * carry + *digit;
            *digit = value / 2;
            carry = value % 2;
        }
        bits.push(carry != 0);
        if digits[0] == 0 {
            digits = &mut digits[1..];
        }
    }
    Ok(bits)
}
