use crate::error::{Error, Result};
use crate::model::Model;

#[derive(Clone, Copy)]
pub(crate) struct BinaryField {
    modulus: usize,
    order: usize,
    degree: usize,
}

impl BinaryField {
    pub(crate) fn new(model: &mut Model, modulus: usize) -> Result<Self> {
        if modulus < 2 {
            return Err(Error::plain(
                "field modulus must have positive polynomial degree",
            ));
        }
        let degree = modulus.ilog2() as usize;
        let field = Self {
            modulus,
            order: 1 << degree,
            degree,
        };
        let x = polynomial_remainder(2, modulus);
        let mut power = x;
        for exponent in 1..=degree {
            model.spend(1 + 2 * degree)?;
            power = field.multiply(power, power);
            // Any proper factor of degree k <= m/2 divides X^(2^k)-X
            // when it is irreducible of degree k.
            if exponent <= degree / 2 && polynomial_gcd(modulus, power ^ x) != 1 {
                return Err(Error::plain("field modulus is reducible over F2"));
            }
        }
        if power != x {
            return Err(Error::plain("field modulus is reducible over F2"));
        }
        Ok(field)
    }

    pub(crate) fn order(self) -> usize {
        self.order
    }

    pub(crate) fn degree(self) -> usize {
        self.degree
    }

    pub(crate) fn multiply(self, mut a: usize, mut b: usize) -> usize {
        let mut output = 0;
        while b != 0 {
            if b & 1 != 0 {
                output ^= a;
            }
            b >>= 1;
            a <<= 1;
            if a & self.order != 0 {
                a ^= self.modulus;
            }
        }
        output
    }

    pub(crate) fn inverse(self, model: &mut Model, mut value: usize) -> Result<usize> {
        debug_assert!(value != 0 && value < self.order);
        let mut exponent = self.order - 2;
        let mut output = 1;
        while exponent != 0 {
            model.spend(2 * self.degree)?;
            if exponent & 1 != 0 {
                output = self.multiply(output, value);
            }
            value = self.multiply(value, value);
            exponent >>= 1;
        }
        Ok(output)
    }
}

fn polynomial_remainder(mut value: usize, divisor: usize) -> usize {
    while value != 0 && value.ilog2() >= divisor.ilog2() {
        value ^= divisor << (value.ilog2() - divisor.ilog2());
    }
    value
}

fn polynomial_gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        (a, b) = (b, polynomial_remainder(a, b));
    }
    a
}
