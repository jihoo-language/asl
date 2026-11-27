use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hyperelliptic_kummer_surface::square_root;
use crate::model::Model;

use super::ResidueSquares;

#[derive(Clone, Eq, PartialEq)]
pub(super) struct Laurent {
    pub(super) order: i64,
    pub(super) coefficients: Polynomial,
}

impl Laurent {
    pub(super) fn new(mut order: i64, coefficients: Polynomial) -> Self {
        let mut coefficients = poly::trim(coefficients);
        if let Some(first) = coefficients.iter().position(|c| !poly::zero(c)) {
            coefficients.drain(..first);
            order += first as i64;
        } else {
            order = 0;
        }
        Self {
            order,
            coefficients,
        }
    }

    pub(super) fn is_zero(&self) -> bool {
        self.coefficients.is_empty()
    }

    pub(super) fn shift(mut self, amount: i64) -> Self {
        if !self.is_zero() {
            self.order += amount;
        }
        self
    }

    pub(super) fn symbol(&self, squares: &ResidueSquares) -> u8 {
        squares.symbol(&self.coefficients[0], self.order.rem_euclid(2) as usize)
    }

    pub(super) fn product(
        &self,
        model: &mut Model,
        field: Field<'_>,
        right: &Self,
    ) -> Result<Self> {
        let coefficients = poly::multiply(model, field, &self.coefficients, &right.coefficients)?;
        Ok(Self::new(self.order + right.order, coefficients))
    }

    pub(super) fn sum(
        &self,
        model: &mut Model,
        field: Field<'_>,
        right: &Self,
        subtract: bool,
    ) -> Result<Self> {
        let start = if self.is_zero() {
            right.order
        } else if right.is_zero() {
            self.order
        } else {
            self.order.min(right.order)
        };
        let align = |model: &mut Model, term: &Laurent| -> Result<Polynomial> {
            if term.is_zero() {
                return Ok(Vec::new());
            }
            let shift = (term.order - start) as usize;
            field.check_polynomial(shift + term.coefficients.len())?;
            model.spend(shift + term.coefficients.len())?;
            let mut padded = vec![field.zero(); shift];
            padded.extend_from_slice(&term.coefficients);
            Ok(padded)
        };
        let a = align(model, self)?;
        let b = align(model, right)?;
        Ok(Self::new(start, poly::sum(model, field, &a, &b, subtract)?))
    }

    fn coefficient(&self, order: i64) -> Option<&Scalar> {
        usize::try_from(order - self.order)
            .ok()
            .and_then(|i| self.coefficients.get(i))
    }
}

pub(super) struct Root<'a> {
    radicand: &'a Laurent,
    divisor: &'a Laurent,
    order: i64,
    inverse_divisor: Scalar,
    inverse_derivative: Scalar,
    quotient: Polynomial,
    coefficients: Polynomial,
}

impl<'a> Root<'a> {
    pub(super) fn new(
        model: &mut Model,
        field: Field<'_>,
        radicand: &'a Laurent,
        divisor: &'a Laurent,
        weight: i64,
    ) -> Result<Self> {
        // The input x-chart was scaled by s^weight at infinity.
        let value_order = radicand.order - divisor.order - weight;
        if value_order.rem_euclid(2) != 0 {
            return Err(Error::plain(
                "Cassels local covering square root has odd valuation",
            ));
        }
        let inverse_divisor = field.inverse_unit(model, &divisor.coefficients[0])?;
        let angular = field.multiply(model, &radicand.coefficients[0], &inverse_divisor)?;
        let root = square_root(model, field, &angular)?.ok_or_else(|| {
            Error::plain("Cassels local covering square root has nonsquare angular component")
        })?;
        let derivative = field.add(model, &root, &root)?;
        let inverse_derivative = field.inverse_unit(model, &derivative)?;
        Ok(Self {
            radicand,
            divisor,
            order: value_order / 2,
            inverse_divisor,
            inverse_derivative,
            quotient: vec![angular],
            coefficients: vec![root],
        })
    }

    pub(super) fn extend(
        &mut self,
        model: &mut Model,
        field: Field<'_>,
        precision: usize,
    ) -> Result<()> {
        field.check_polynomial(precision)?;
        while self.coefficients.len() < precision {
            let n = self.coefficients.len();
            let mut value = self
                .radicand
                .coefficients
                .get(n)
                .cloned()
                .unwrap_or_else(|| field.zero());
            // First divide unit series, then use the simple-root Hensel
            // coefficient equation with derivative 2*root(0).
            for (j, c) in self
                .divisor
                .coefficients
                .iter()
                .enumerate()
                .take(n + 1)
                .skip(1)
            {
                let term = field.multiply(model, c, &self.quotient[n - j])?;
                value = field.subtract(model, &value, &term)?;
            }
            value = field.multiply(model, &value, &self.inverse_divisor)?;
            self.quotient.push(value.clone());
            for j in 1..n {
                let term =
                    field.multiply(model, &self.coefficients[j], &self.coefficients[n - j])?;
                value = field.subtract(model, &value, &term)?;
            }
            self.coefficients
                .push(field.multiply(model, &value, &self.inverse_derivative)?);
        }
        Ok(())
    }
}

pub(super) fn certified_lead(
    model: &mut Model,
    field: Field<'_>,
    coefficients: &[Laurent; 3],
    roots: [&Root<'_>; 2],
    precision: usize,
) -> Result<Option<(i64, Scalar)>> {
    let leading_order = coefficients
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(i, c)| c.order + roots.get(i).map_or(0, |root| root.order))
        .min()
        .ok_or_else(|| Error::plain("Cassels tangent linear form vanished"))?;
    // Every coefficient in this interval is exact. Each unknown root
    // coefficient starts at least precision higher than its term's order.
    for k in 0..precision {
        let order = leading_order + k as i64;
        let mut value = coefficients[2]
            .coefficient(order)
            .cloned()
            .unwrap_or_else(|| field.zero());
        for (c, root) in coefficients.iter().zip(roots) {
            if c.is_zero() {
                continue;
            }
            let Ok(n) = usize::try_from(order - c.order - root.order) else {
                continue;
            };
            for (j, a) in c.coefficients.iter().enumerate().take(n + 1) {
                let term = field.multiply(model, a, &root.coefficients[n - j])?;
                value = field.add(model, &value, &term)?;
            }
        }
        if !field.is_zero(&value) {
            return Ok(Some((order, value)));
        }
    }
    Ok(None)
}
