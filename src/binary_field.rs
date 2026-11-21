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

    pub(crate) fn frobenius(
        self,
        model: &mut Model,
        mut value: usize,
        iterations: usize,
    ) -> Result<usize> {
        debug_assert!(value < self.order);
        if value > 1 {
            let iterations = iterations % self.degree;
            model.spend(iterations * self.degree)?;
            for _ in 0..iterations {
                value = self.multiply(value, value);
            }
        }
        Ok(value)
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

pub(crate) struct Embedding {
    small: BinaryField,
    large: BinaryField,
    powers: Vec<usize>,
    pivots: Vec<(usize, usize)>,
}

impl Embedding {
    pub(crate) fn new(
        model: &mut Model,
        small: BinaryField,
        large: BinaryField,
        root: usize,
    ) -> Result<Self> {
        if !large.degree.is_multiple_of(small.degree) {
            return Err(Error::plain(
                "field embedding requires the smaller degree to divide the larger degree",
            ));
        }
        large.check_element(root)?;
        model.spend((1 + small.degree) * large.degree)?;
        let mut value = 0;
        for degree in (0..=small.degree).rev() {
            value = large.multiply(value, root) ^ ((small.modulus >> degree) & 1);
        }
        if value != 0 {
            return Err(Error::plain(
                "field embedding root does not satisfy the smaller modulus",
            ));
        }

        let mut powers = Vec::with_capacity(small.degree);
        let mut pivots = vec![(0, 0); large.degree];
        let mut power = 1usize;
        for degree in 0..small.degree {
            powers.push(power);
            let mut image = power;
            let mut preimage = 1 << degree;
            while image != 0 {
                model.spend(1)?;
                let row = &mut pivots[image.ilog2() as usize];
                if row.0 == 0 {
                    *row = (image, preimage);
                    break;
                }
                image ^= row.0;
                preimage ^= row.1;
            }
            if image == 0 {
                return Err(Error::plain("field embedding basis is linearly dependent"));
            }
            model.spend(large.degree)?;
            power = large.multiply(power, root);
        }
        model.allocate(powers.len() + 2 * pivots.len())?;
        Ok(Self {
            small,
            large,
            powers,
            pivots,
        })
    }

    pub(crate) fn small(&self) -> BinaryField {
        self.small
    }

    pub(crate) fn large(&self) -> BinaryField {
        self.large
    }

    pub(crate) fn image(&self, model: &mut Model, mut value: usize) -> Result<usize> {
        debug_assert!(value < self.small.order);
        model.spend(value.count_ones() as usize)?;
        let mut image = 0;
        while value != 0 {
            image ^= self.powers[value.trailing_zeros() as usize];
            value &= value - 1;
        }
        Ok(image)
    }

    pub(crate) fn preimage(&self, model: &mut Model, mut value: usize) -> Result<Option<usize>> {
        debug_assert!(value < self.large.order);
        let mut preimage = 0;
        while value != 0 {
            model.spend(1)?;
            let (image, witness) = self.pivots[value.ilog2() as usize];
            if image == 0 {
                return Ok(None);
            }
            value ^= image;
            preimage ^= witness;
        }
        Ok(Some(preimage))
    }

    pub(crate) fn trace(&self, model: &mut Model, mut value: usize) -> Result<usize> {
        let mut trace = value;
        for _ in 1..self.large.degree / self.small.degree {
            model.spend(1)?;
            value = self.large.frobenius(model, value, self.small.degree)?;
            trace ^= value;
        }
        self.preimage(model, trace)?.ok_or_else(|| {
            Error::plain("relative field trace did not land in the specified subfield")
        })
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
