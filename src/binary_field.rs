use crate::error::{Error, Result};
use crate::model::Model;

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
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
        if !field.is_irreducible(model)? {
            return Err(Error::plain("field modulus is reducible over F2"));
        }
        Ok(field)
    }

    pub(crate) fn of_degree(model: &mut Model, degree: usize) -> Result<Self> {
        if degree == 0 || degree >= usize::BITS as usize {
            return Err(Error::plain("binary extension degree is out of range"));
        }
        let order = 1 << degree;
        // Scan monic polynomials with nonzero constant term in encoding order.
        for lower in (1..order).step_by(2) {
            let field = Self {
                modulus: order | lower,
                order,
                degree,
            };
            if field.is_irreducible(model)? {
                return Ok(field);
            }
        }
        Err(Error::plain(
            "no irreducible modulus found for binary extension",
        ))
    }

    fn is_irreducible(self, model: &mut Model) -> Result<bool> {
        let x = polynomial_remainder(2, self.modulus);
        let mut power = x;
        for exponent in 1..=self.degree {
            model.spend(1 + 2 * self.degree)?;
            power = self.multiply(power, power);
            // Any proper factor of degree k <= m/2 divides X^(2^k)-X
            // when it is irreducible of degree k.
            if exponent <= self.degree / 2 && polynomial_gcd(self.modulus, power ^ x) != 1 {
                return Ok(false);
            }
        }
        Ok(power == x)
    }

    pub(crate) fn order(self) -> usize {
        self.order
    }

    pub(crate) fn modulus(self) -> usize {
        self.modulus
    }

    pub(crate) fn check_element(self, value: usize) -> Result<()> {
        if value >= self.order {
            return Err(Error::plain(format!(
                "field coefficient {value} is out of range for F{}",
                self.order
            )));
        }
        Ok(())
    }

    pub(crate) fn check_point(
        self,
        coefficients: &[usize],
        expected: usize,
        name: &str,
    ) -> Result<()> {
        if coefficients.len() != expected {
            return Err(Error::plain(format!(
                "{name} requires one field coefficient per selected primitive"
            )));
        }
        for &value in coefficients {
            self.check_element(value)?;
        }
        Ok(())
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

    pub(crate) fn inverse(self, model: &mut Model, value: usize) -> Result<usize> {
        debug_assert!(value != 0 && value < self.order);
        self.power(model, value, (self.order - 2) as u128)
    }

    pub(crate) fn power(
        self,
        model: &mut Model,
        mut value: usize,
        mut exponent: u128,
    ) -> Result<usize> {
        debug_assert!(value < self.order);
        if exponent == 0 || value == 1 {
            return Ok(1);
        }
        if value == 0 {
            return Ok(0);
        }
        exponent %= (self.order - 1) as u128;
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
