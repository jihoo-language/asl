use std::cmp::Ordering;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;

const MAX_DEGREE: usize = 16;

pub(crate) type Coordinates = Rc<[u64]>;
type Element = Vec<u64>;

pub(crate) struct Ring {
    prime: Prime,
    precision: usize,
    modulus: u64,
    lower: Box<[u64]>,
    frobenius_images: Vec<Element>,
}

impl Ring {
    pub(crate) fn new(
        model: &mut Model,
        prime: usize,
        precision: usize,
        polynomial: Vec<i128>,
    ) -> Result<Self> {
        if prime < 2 || (prime != 2 && prime.is_multiple_of(2)) {
            return Err(Error::plain(
                "unramified coefficient characteristic must be prime",
            ));
        }
        let p = if prime == 2 {
            Prime::validated(2)
        } else {
            Prime::new(model, prime)?
        };
        let Some((&1, lower)) = polynomial.split_last() else {
            return Err(Error::plain("unramified modulus polynomial must be monic"));
        };
        if lower.is_empty() || lower.len() > MAX_DEGREE {
            return Err(Error::plain(format!(
                "unramified modulus polynomial degree must be in 1..={MAX_DEGREE}"
            )));
        }
        if precision == 0 {
            return Err(Error::plain("unramified precision must be positive"));
        }
        let modulus = u32::try_from(precision)
            .ok()
            .and_then(|e| u64::from(p.get()).checked_pow(e))
            .ok_or_else(|| Error::plain("unramified coefficient modulus p^s exceeds u64"))?;
        model.spend(lower.len())?;
        let lower: Box<[_]> = lower.iter().map(|&c| reduce(c, modulus)).collect();
        let ring = Self {
            prime: p,
            precision,
            modulus,
            lower,
            frobenius_images: Vec::new(),
        };
        if !ring.is_irreducible(model)? {
            return Err(Error::plain(
                "unramified modulus polynomial must be irreducible modulo p",
            ));
        }
        ring.finish(model)
    }

    pub(crate) fn field_of_degree(model: &mut Model, prime: Prime, degree: usize) -> Result<Self> {
        if !(1..=MAX_DEGREE).contains(&degree) {
            return Err(Error::plain(format!(
                "residue field degree must be in 1..={MAX_DEGREE}"
            )));
        }
        let p = u64::from(prime.get());
        let count = p
            .checked_pow(degree as u32)
            .ok_or_else(|| Error::plain("residue field size exceeds u64"))?;
        let mut ring = Self {
            prime,
            precision: 1,
            modulus: p,
            lower: vec![0; degree].into(),
            frobenius_images: Vec::new(),
        };
        if degree == 1 {
            return ring.finish(model);
        }
        for mut index in 1..count {
            model.spend(1)?;
            if index.is_multiple_of(p) {
                continue;
            }
            model.spend(degree)?;
            for coefficient in &mut ring.lower {
                *coefficient = index % p;
                index /= p;
            }
            if ring.is_irreducible(model)? {
                return ring.finish(model);
            }
        }
        Err(Error::plain(
            "could not find a monic irreducible residue polynomial",
        ))
    }

    fn finish(mut self, model: &mut Model) -> Result<Self> {
        self.frobenius_images = self.lift_frobenius(model)?;
        model.allocate(8 + self.degree() * (4 + self.degree()))?;
        Ok(self)
    }

    pub(crate) fn prime(&self) -> u64 {
        u64::from(self.prime.get())
    }

    pub(crate) fn precision(&self) -> usize {
        self.precision
    }

    pub(crate) fn modulus(&self) -> u64 {
        self.modulus
    }

    pub(crate) fn degree(&self) -> usize {
        self.lower.len()
    }

    pub(crate) fn count(&self) -> usize {
        usize::try_from(self.modulus)
            .unwrap_or(usize::MAX)
            .saturating_pow(self.degree() as u32)
    }

    pub(crate) fn size_of_residue_field(&self) -> Result<u64> {
        self.prime()
            .checked_pow(self.degree() as u32)
            .ok_or_else(|| Error::plain("unramified residue field size exceeds u64"))
    }

    pub(crate) fn scalar(&self, value: u64) -> Element {
        self.algebra(self.modulus).constant(value)
    }

    pub(crate) fn literal(&self, model: &mut Model, coefficients: Vec<i128>) -> Result<Element> {
        if coefficients.len() > self.degree() {
            return Err(Error::plain(format!(
                "unramified element has at most {} coordinates",
                self.degree()
            )));
        }
        model.spend(self.degree())?;
        let mut value: Element = coefficients
            .into_iter()
            .map(|c| reduce(c, self.modulus))
            .collect();
        value.resize(self.degree(), 0);
        Ok(value)
    }

    fn algebra(&self, modulus: u64) -> Algebra<'_> {
        Algebra {
            lower: &self.lower,
            modulus,
        }
    }

    fn is_irreducible(&self, model: &mut Model) -> Result<bool> {
        let p = self.prime();
        let alg = self.algebra(p);
        let x = alg.generator();
        let mut power = x.clone();
        for _ in 1..=self.degree() / 2 {
            power = alg.pow(model, &power, p)?;
            let defect = alg.subtract(model, &power, &x)?;
            let mut polynomial: Vec<_> = self.lower.iter().map(|x| (x % p) as u32).collect();
            polynomial.push(1);
            if !coprime(
                model,
                self.prime,
                polynomial,
                defect.into_iter().map(|x| x as u32).collect(),
            )? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn lift_frobenius(&self, model: &mut Model) -> Result<Vec<Element>> {
        let p = self.prime();
        let residue = self.algebra(p);
        let mut root = residue.pow(model, &residue.generator(), p)?;
        let mut polynomial = self.lower.to_vec();
        polynomial.push(1);
        if self.precision > 1 {
            let derivative: Vec<_> = polynomial
                .iter()
                .enumerate()
                .skip(1)
                .map(|(i, &c)| mul_mod(i as u64, c, p))
                .collect();
            let derivative = residue.evaluate(model, &derivative, &root)?;
            let inverse = self.residue_inverse(model, &derivative)?;
            let mut previous = p;
            for _ in 1..self.precision {
                let next = previous * p;
                let defect = self.algebra(next).evaluate(model, &polynomial, &root)?;
                let digit: Vec<_> = defect
                    .into_iter()
                    .map(|c| neg_mod(c / previous, p))
                    .collect();
                let correction = residue.multiply(model, &digit, &inverse)?;
                model.spend(self.degree())?;
                for (value, c) in root.iter_mut().zip(correction) {
                    *value += c * previous;
                }
                previous = next;
            }
        }
        let algebra = self.algebra(self.modulus);
        if algebra
            .evaluate(model, &polynomial, &root)?
            .iter()
            .any(|&x| x != 0)
        {
            return Err(Error::plain(
                "unramified Frobenius root failed its polynomial identity",
            ));
        }
        let generator = algebra.generator();
        let mut images = Vec::with_capacity(self.degree());
        let mut value = generator.clone();
        for _ in 0..self.degree() {
            let next = algebra.evaluate(model, &value, &root)?;
            images.push(value);
            value = next;
        }
        if value != generator {
            return Err(Error::plain(
                "unramified Frobenius failed its residue-degree identity",
            ));
        }
        Ok(images)
    }

    // In the residue field, multiply all conjugates except the input and
    // invert the resulting scalar norm in the prime field.
    fn residue_inverse(&self, model: &mut Model, x: &[u64]) -> Result<Element> {
        let p = self.prime();
        let algebra = self.algebra(p);
        let x: Element = x.iter().map(|&c| c % p).collect();
        if x.iter().all(|&c| c == 0) {
            return Err(Error::plain("unramified inverse requires a unit modulo p"));
        }
        let mut cofactor = algebra.constant(1);
        let mut sigma = x.clone();
        for _ in 1..self.degree() {
            sigma = algebra.pow(model, &sigma, p)?;
            cofactor = algebra.multiply(model, &cofactor, &sigma)?;
        }
        let norm = algebra.multiply(model, &x, &cofactor)?;
        let norm = require_scalar(&norm)?;
        let inverse = u64::from(invert(model, self.prime, norm as u32)?);
        algebra.scale(model, &cofactor, inverse)
    }

    pub(crate) fn inverse(&self, model: &mut Model, x: &[u64]) -> Result<Element> {
        let mut inverse = self.residue_inverse(model, x)?;
        let algebra = self.algebra(self.modulus);
        let two = algebra.constant(2);
        let mut precision = 1;
        while precision < self.precision {
            let product = algebra.multiply(model, x, &inverse)?;
            let correction = algebra.subtract(model, &two, &product)?;
            inverse = algebra.multiply(model, &inverse, &correction)?;
            precision *= 2;
        }
        Ok(inverse)
    }

    pub(crate) fn add(&self, model: &mut Model, x: &[u64], y: &[u64]) -> Result<Element> {
        self.algebra(self.modulus).add(model, x, y)
    }

    pub(crate) fn negate(&self, model: &mut Model, x: &[u64]) -> Result<Element> {
        self.algebra(self.modulus).scale(model, x, self.modulus - 1)
    }

    pub(crate) fn scale(&self, model: &mut Model, x: &[u64], c: u64) -> Result<Element> {
        self.algebra(self.modulus).scale(model, x, c)
    }

    pub(crate) fn multiply(&self, model: &mut Model, x: &[u64], y: &[u64]) -> Result<Element> {
        self.algebra(self.modulus).multiply(model, x, y)
    }

    pub(crate) fn power(&self, model: &mut Model, x: &[u64], k: u64) -> Result<Element> {
        self.algebra(self.modulus).pow(model, x, k)
    }

    pub(crate) fn frobenius(&self, model: &mut Model, x: &[u64], k: usize) -> Result<Element> {
        let k = k % self.degree();
        if k == 0 {
            model.spend(self.degree())?;
            return Ok(x.to_vec());
        }
        self.algebra(self.modulus)
            .evaluate(model, x, &self.frobenius_images[k])
    }

    pub(crate) fn teichmuller(&self, model: &mut Model, x: &[u64]) -> Result<Element> {
        let p = self.prime();
        let mut value: Element = x.iter().map(|&c| c % p).collect();
        let algebra = self.algebra(self.modulus);
        let degree = self.degree();
        // A p-power increases the congruence with the corresponding
        // Teichmüller lift by one digit. Use complete residue orbits so
        // the power returns to the specified residue class.
        for _ in 0..(self.precision - 1).div_ceil(degree) * degree {
            value = algebra.pow(model, &value, p)?;
        }
        Ok(value)
    }

    pub(crate) fn conjugate_sum_or_product(
        &self,
        model: &mut Model,
        x: &[u64],
        product: bool,
    ) -> Result<Element> {
        let mut result = self.scalar(u64::from(product));
        for i in 0..self.degree() {
            let conjugate = self.frobenius(model, x, i)?;
            result = if product {
                self.multiply(model, &result, &conjugate)?
            } else {
                self.add(model, &result, &conjugate)?
            };
        }
        require_scalar(&result)?;
        Ok(result)
    }

    pub(crate) fn valuation(&self, model: &mut Model, x: &[u64]) -> Result<usize> {
        let p = self.prime();
        let mut valuation = self.precision;
        for mut c in x.iter().copied().filter(|&c| c != 0) {
            model.spend(1)?;
            let mut order = 0;
            while c.is_multiple_of(p) && order < valuation {
                model.spend(1)?;
                c /= p;
                order += 1;
            }
            valuation = valuation.min(order);
        }
        Ok(valuation)
    }

    pub(crate) fn residue(&self, x: &[u64]) -> Element {
        x.iter().map(|&c| c % self.prime()).collect()
    }

    pub(crate) fn residue_ring(&self, model: &mut Model) -> Result<Self> {
        let p = self.prime();
        let h = self.degree();
        model.spend(h * (h + 1))?;
        model.allocate(8 + h * (4 + h))?;
        Ok(Self {
            prime: self.prime,
            precision: 1,
            modulus: p,
            lower: self.residue(&self.lower).into(),
            frobenius_images: self
                .frobenius_images
                .iter()
                .map(|image| self.residue(image))
                .collect(),
        })
    }

    pub(crate) fn reduces_to(&self, target: &Self) -> bool {
        self.prime() == target.prime()
            && self.precision >= target.precision
            && self.degree() == target.degree()
            && self
                .lower
                .iter()
                .zip(&target.lower)
                .all(|(&a, &b)| a % target.modulus == b)
    }

    pub(crate) fn evaluate_polynomial(
        &self,
        model: &mut Model,
        coefficients: &[Element],
        point: &[u64],
    ) -> Result<Element> {
        let Some((leading, rest)) = coefficients.split_last() else {
            return Ok(self.scalar(0));
        };
        model.spend(self.degree())?;
        let mut output = leading.clone();
        for term in rest.iter().rev() {
            let product = self.multiply(model, &output, point)?;
            output = self.add(model, &product, term)?;
        }
        Ok(output)
    }

    pub(crate) fn embed_residue_coefficients(
        &self,
        model: &mut Model,
        target: &Self,
        coefficients: &[Element],
    ) -> Result<Vec<Element>> {
        if self.prime() != target.prime()
            || target.precision != 1
            || !target.degree().is_multiple_of(self.degree())
        {
            return Err(Error::plain(
                "residue embedding requires an extension of the coefficient field",
            ));
        }
        if self.reduces_to(target) {
            model.spend(coefficients.len().saturating_mul(self.degree()))?;
            return Ok(coefficients
                .iter()
                .map(|value| self.residue(value))
                .collect());
        }
        if self.degree() == 1 {
            model.spend(coefficients.len().saturating_mul(target.degree()))?;
            return Ok(coefficients
                .iter()
                .map(|value| target.scalar(value[0]))
                .collect());
        }

        let mut defining = self.residue(&self.lower);
        defining.push(1);
        let root = target.root_in_residue(model, &defining)?;
        let algebra = target.algebra(target.prime());
        coefficients
            .iter()
            .map(|value| algebra.evaluate(model, &self.residue(value), &root))
            .collect()
    }

    fn root_in_residue(&self, model: &mut Model, polynomial: &[u64]) -> Result<Element> {
        let p = self.prime();
        let algebra = self.algebra(p);
        let mut root = algebra.constant(0);
        for mut index in 0..self.size_of_residue_field()? {
            model.spend(self.degree())?;
            for c in &mut root {
                *c = index % p;
                index /= p;
            }
            if algebra
                .evaluate(model, polynomial, &root)?
                .iter()
                .all(|&c| c == 0)
            {
                return Ok(root);
            }
        }
        Err(Error::plain(
            "defining polynomial has no root in the proposed residue extension",
        ))
    }

    pub(crate) fn legendre_trace(&self, model: &mut Model, lambda: &[u64]) -> Result<i64> {
        let alg = self.algebra(self.prime());
        let lambda = self.residue(lambda);
        let one = alg.constant(1);
        self.residue_quadratic_trace(model, |model, x| {
            if x == lambda || x == one || x.iter().all(|&c| c == 0) {
                return Ok(alg.constant(0));
            }
            let x_one = alg.subtract(model, x, &one)?;
            let x_lambda = alg.subtract(model, x, &lambda)?;
            let rhs = alg.multiply(model, x, &x_one)?;
            alg.multiply(model, &rhs, &x_lambda)
        })
    }

    pub(crate) fn residue_polynomial_trace(
        &self,
        model: &mut Model,
        polynomial: &[u64],
    ) -> Result<i64> {
        let alg = self.algebra(self.prime());
        self.residue_quadratic_trace(model, |model, x| alg.evaluate(model, polynomial, x))
    }

    pub(crate) fn residue_quadratic_trace(
        &self,
        model: &mut Model,
        mut evaluate: impl FnMut(&mut Model, &[u64]) -> Result<Element>,
    ) -> Result<i64> {
        let q = self.size_of_residue_field()?;
        let count = usize::try_from(q)
            .map_err(|_| Error::plain("quadratic character point field is too large"))?;
        model.spend(count)?;
        let p = self.prime();
        if p == 2 {
            return Err(Error::plain("quadratic character requires an odd prime"));
        }
        let alg = self.algebra(p);
        let exponent = (q - 1) / 2;
        let mut x = alg.constant(0);
        let mut trace = 0;
        for mut index in 0..q {
            model.spend(self.degree())?;
            for c in &mut x {
                *c = index % p;
                index /= p;
            }
            let rhs = evaluate(model, &x)?;
            if rhs.iter().all(|c| c.is_multiple_of(p)) {
                continue;
            }
            let character = alg.pow(model, &rhs, exponent)?;
            trace += match require_scalar(&character)? {
                1 => -1,
                c if c == p - 1 => 1,
                _ => {
                    return Err(Error::plain(
                        "nonzero residue field value must have quadratic character +1 or -1",
                    ));
                }
            };
        }
        Ok(trace)
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Negate,
    Inverse,
    Power(u64),
    Frobenius(usize),
    Teichmuller,
    Trace,
    Norm,
}

impl Unary {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Negate => "neg",
            Self::Inverse => "rvUnramInv",
            Self::Power(_) => "rvUnramPow",
            Self::Frobenius(_) => "rvUnramFrob",
            Self::Teichmuller => "rvUnramTeich",
            Self::Trace => "rvUnramTrace",
            Self::Norm => "rvUnramNorm",
        }
    }

    pub(crate) fn apply(self, model: &mut Model, ring: &Ring, value: &[u64]) -> Result<Element> {
        match self {
            Self::Negate => ring.negate(model, value),
            Self::Inverse => ring.inverse(model, value),
            Self::Power(k) => ring.power(model, value, k),
            Self::Frobenius(k) => ring.frobenius(model, value, k),
            Self::Teichmuller => ring.teichmuller(model, value),
            Self::Trace => ring.conjugate_sum_or_product(model, value, false),
            Self::Norm => ring.conjugate_sum_or_product(model, value, true),
        }
    }
}

pub(crate) fn persist(model: &mut Model, values: Element) -> Result<Value> {
    model.allocate(3 + values.len())?;
    Ok(Value::RVUnram(values.into()))
}

pub(crate) fn compare(model: &mut Model, left: &[u64], right: &[u64]) -> Result<Ordering> {
    model.spend(left.len().min(right.len()))?;
    Ok(left.cmp(right))
}

pub(crate) fn enumerate(model: &mut Model, ring: &Ring, count: usize) -> Result<Vec<Value>> {
    (0..count)
        .map(|mut i| {
            model.spend(ring.degree())?;
            let mut row = ring.scalar(0);
            for c in row.iter_mut().rev() {
                *c = (i as u64) % ring.modulus;
                i /= ring.modulus as usize;
            }
            persist(model, row)
        })
        .collect()
}

pub(crate) fn require_scalar(x: &[u64]) -> Result<u64> {
    if x[1..].iter().any(|&coefficient| coefficient != 0) {
        return Err(Error::plain("expected descent to the unramified base ring"));
    }
    Ok(x[0])
}

struct Algebra<'a> {
    lower: &'a [u64],
    modulus: u64,
}

impl Algebra<'_> {
    fn constant(&self, c: u64) -> Element {
        let mut values = vec![0; self.lower.len()];
        values[0] = c % self.modulus;
        values
    }

    fn generator(&self) -> Element {
        if self.lower.len() == 1 {
            return self.constant(neg_mod(self.lower[0] % self.modulus, self.modulus));
        }
        let mut value = self.constant(0);
        value[1] = 1;
        value
    }

    fn add(&self, model: &mut Model, x: &[u64], y: &[u64]) -> Result<Element> {
        model.spend(self.lower.len())?;
        Ok(x.iter()
            .zip(y)
            .map(|(&a, &b)| add_mod(a, b, self.modulus))
            .collect())
    }

    fn subtract(&self, model: &mut Model, x: &[u64], y: &[u64]) -> Result<Element> {
        model.spend(self.lower.len())?;
        Ok(x.iter()
            .zip(y)
            .map(|(&a, &b)| sub_mod(a, b, self.modulus))
            .collect())
    }

    fn scale(&self, model: &mut Model, x: &[u64], c: u64) -> Result<Element> {
        model.spend(self.lower.len())?;
        Ok(x.iter().map(|&a| mul_mod(a, c, self.modulus)).collect())
    }

    fn multiply(&self, model: &mut Model, x: &[u64], y: &[u64]) -> Result<Element> {
        let n = self.lower.len();
        let modulus = self.modulus;
        model.spend(2 * n * (n + 1))?;
        let mut product = vec![0; 2 * n - 1];
        for (i, &a) in x.iter().enumerate().filter(|&(_, &a)| a != 0) {
            for (j, &b) in y.iter().enumerate().filter(|&(_, &b)| b != 0) {
                product[i + j] = add_mod(product[i + j], mul_mod(a, b, modulus), modulus);
            }
        }
        for i in (n..product.len()).rev() {
            let coefficient = product[i];
            if coefficient != 0 {
                for (j, &h) in self.lower.iter().enumerate() {
                    product[i - n + j] = sub_mod(
                        product[i - n + j],
                        mul_mod(coefficient, h, modulus),
                        modulus,
                    );
                }
            }
        }
        product.truncate(n);
        Ok(product)
    }

    fn pow(&self, model: &mut Model, a: &[u64], mut k: u64) -> Result<Element> {
        let mut result = self.constant(1);
        let mut power = a.to_vec();
        while k != 0 {
            if !k.is_multiple_of(2) {
                result = self.multiply(model, &result, &power)?;
            }
            k >>= 1;
            if k != 0 {
                power = self.multiply(model, &power, &power)?;
            }
        }
        Ok(result)
    }

    fn evaluate(&self, model: &mut Model, polynomial: &[u64], x: &[u64]) -> Result<Element> {
        let Some((leading, lower)) = polynomial.split_last() else {
            return Ok(self.constant(0));
        };
        let mut output = self.constant(*leading);
        for &coefficient in lower.iter().rev() {
            output = self.multiply(model, &output, x)?;
            output[0] = add_mod(output[0], coefficient % self.modulus, self.modulus);
        }
        Ok(output)
    }
}

pub(crate) fn coprime(
    model: &mut Model,
    prime: Prime,
    mut a: Vec<u32>,
    mut b: Vec<u32>,
) -> Result<bool> {
    trim(&mut b);
    while let Some(&leading) = b.last() {
        let inverse = invert(model, prime, leading)?;
        while a.len() >= b.len() {
            model.spend(b.len())?;
            let shift = a.len() - b.len();
            let scale = prime.negate(prime.mul(*a.last().expect("nonzero remainder"), inverse));
            for (i, &c) in b.iter().enumerate() {
                a[i + shift] = prime.add(a[i + shift], prime.mul(scale, c));
            }
            trim(&mut a);
        }
        (a, b) = (b, a);
    }
    Ok(a.len() == 1)
}

fn trim(polynomial: &mut Vec<u32>) {
    while polynomial.last() == Some(&0) {
        polynomial.pop();
    }
}

pub(crate) fn reduce(x: i128, modulus: u64) -> u64 {
    x.rem_euclid(i128::from(modulus)) as u64
}

pub(crate) fn add_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((u128::from(a) + u128::from(b)) % u128::from(modulus)) as u64
}

pub(crate) fn mul_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((u128::from(a) * u128::from(b)) % u128::from(modulus)) as u64
}

fn neg_mod(x: u64, modulus: u64) -> u64 {
    if x == 0 {
        0
    } else {
        modulus - x
    }
}

pub(crate) fn sub_mod(a: u64, b: u64, modulus: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        modulus - (b - a)
    }
}
