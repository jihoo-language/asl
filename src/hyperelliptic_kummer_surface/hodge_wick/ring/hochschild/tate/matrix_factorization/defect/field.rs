use crate::hw_algebra::CoefficientRing;
use crate::rv_cyclotomic as cyclo;
use crate::rv_lens::{self, Rational as Scalar};

use super::{Error, Model, Result};

pub(super) struct Field {
    pub(super) d: usize,
    degree: usize,
}

impl Field {
    pub(super) fn new(model: &mut Model, d: usize) -> Result<Self> {
        let degree = model.cyclotomic(d)?.len() - 1;
        Ok(Self { d, degree })
    }
}

impl CoefficientRing for Field {
    type Element = Scalar;

    fn cells(&self) -> usize {
        4 + 2 * self.degree
    }

    fn zero(&self) -> Scalar {
        Scalar::zero()
    }

    fn one(&self) -> Scalar {
        Scalar::integer(1)
    }

    fn is_zero(&self, value: &Scalar) -> bool {
        value.is_zero()
    }

    fn is_unit(&self, value: &Scalar) -> bool {
        !value.is_zero()
    }

    fn add(&self, model: &mut Model, a: &Scalar, b: &Scalar) -> Result<Scalar> {
        rv_lens::sum(model, a, b)
    }

    fn subtract(&self, model: &mut Model, a: &Scalar, b: &Scalar) -> Result<Scalar> {
        let negative = rv_lens::negative(model, b)?;
        rv_lens::sum(model, a, &negative)
    }

    fn multiply(&self, model: &mut Model, a: &Scalar, b: &Scalar) -> Result<Scalar> {
        rv_lens::multiply(model, self.d, a, b)
    }

    fn inverse_unit(&self, model: &mut Model, value: &Scalar) -> Result<Scalar> {
        let bound = || Error::plain("defect field inversion exceeds the signed coefficient range");
        let denominator: i128 = value.denominator.try_into().map_err(|_| bound())?;
        let mut terms = value.numerator.iter().enumerate().filter(|(_, &c)| c != 0);
        let Some((exponent, &coefficient)) = terms.next() else {
            return Err(Error::plain(
                "division by zero in the defect coefficient field",
            ));
        };
        if terms.next().is_none() {
            let signed = if coefficient < 0 {
                -denominator
            } else {
                denominator
            };
            let inverse = rv_lens::fraction(model, &[signed], coefficient.unsigned_abs())?;
            return rv_lens::shifted(model, self.d, &inverse, -(exponent as i128));
        }

        // The product of all other Galois conjugates times the element is its rational norm.
        let mut complement = vec![1];
        model.spend(self.d)?;
        for k in (2..self.d).filter(|&k| rv_lens::gcd(k as u128, self.d as u128) == 1) {
            let conjugate = cyclo::substitute_root(model, self.d, &value.numerator, k)?;
            complement = cyclo::product_coefficients(model, self.d, &complement, &conjugate)?;
        }
        let norm = cyclo::product_coefficients(model, self.d, &complement, &value.numerator)?;
        let [n] = *norm else {
            return Err(Error::plain("defect field norm is not a nonzero integer"));
        };
        let inverse = rv_lens::fraction(model, &complement, n.unsigned_abs())?;
        let multiplier = Scalar::integer(if n < 0 { -denominator } else { denominator });
        let inverse = self.multiply(model, &inverse, &multiplier)?;
        if self.multiply(model, &inverse, value)? != self.one() {
            return Err(Error::plain("defect field inverse failed verification"));
        }
        Ok(inverse)
    }
}
