use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_ideal::Natural;
use crate::syntax::{Label, ParsedExpr};

use super::cubic_algebra;
use super::cubic_hilbert::{self as cubic, Class};
use super::kummer::divide_by_lambda;
use super::{Algebra as BaseAlgebra, Element as BaseElement, Extension as Base};

pub(crate) struct Declaration(
    pub(crate) Label,
    pub(crate) Label,
    pub(crate) usize,
    pub(crate) usize,
);

pub(crate) enum Expression {
    Literal(Label, Vec<Vec<i128>>),
    Generator(Label, Generator),
    Convert(Label, Conversion, Box<ParsedExpr>),
    Binary(Binary, Box<[ParsedExpr; 2]>),
    Unary(Unary, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Generator {
    Pi,
    Theta,
}

#[derive(Clone, Copy)]
pub(crate) enum Conversion {
    Symbol,
    Reduce,
}

#[derive(Clone, Copy)]
pub(crate) enum Binary {
    Multiply,
    Commutator,
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Negate,
    Inverse,
    Norm(usize),
    IsUnit,
    ConjugatePi,
}

pub(crate) enum Report {
    Data(Label),
    Valuation(Box<ParsedExpr>),
    Commutator(Box<[ParsedExpr; 2]>),
}

type Coefficient = [BaseElement; 3];

struct Unramified<'a> {
    base: BaseAlgebra<'a>,
    relation: Coefficient,
    images: [[Coefficient; 2]; 2],
}

impl<'a> Unramified<'a> {
    fn new(model: &mut Model, field: &'a Base, precision: usize) -> Result<Self> {
        let base = field.algebra(precision)?;
        let one = base.scalar(1);
        let lambda = base.pi();
        let linear = base.add(model, &base.scalar(2), &lambda)?;
        let quadratic = base.add(model, &base.scalar(3), &lambda)?;
        let zeta = base.add(model, &one, &lambda)?;
        let zeta2 = base.multiply(model, &zeta, &zeta)?;
        let theta1 = [one.clone(), zeta.clone(), base.zero()];
        let theta2 = [base.add(model, &one, &zeta)?, zeta2, base.zero()];
        let images = [[theta1, zero(&base)], [theta2, zero(&base)]];
        let mut field = Self {
            relation: [one, base.negate(model, &linear)?, quadratic],
            base,
            images,
        };
        for k in 0..2 {
            field.images[k][1] = field.multiply(model, &field.images[k][0], &field.images[k][0])?;
        }
        Ok(field)
    }

    fn add(&self, model: &mut Model, a: &Coefficient, b: &Coefficient) -> Result<Coefficient> {
        Ok([
            self.base.add(model, &a[0], &b[0])?,
            self.base.add(model, &a[1], &b[1])?,
            self.base.add(model, &a[2], &b[2])?,
        ])
    }

    fn negate(&self, model: &mut Model, x: &Coefficient) -> Result<Coefficient> {
        Ok([
            self.base.negate(model, &x[0])?,
            self.base.negate(model, &x[1])?,
            self.base.negate(model, &x[2])?,
        ])
    }

    fn scale(&self, model: &mut Model, x: &Coefficient, a: &BaseElement) -> Result<Coefficient> {
        Ok([
            self.base.multiply(model, &x[0], a)?,
            self.base.multiply(model, &x[1], a)?,
            self.base.multiply(model, &x[2], a)?,
        ])
    }

    fn multiply(&self, model: &mut Model, a: &Coefficient, b: &Coefficient) -> Result<Coefficient> {
        let mut terms = [(); 5].map(|()| self.base.zero());
        for (i, a) in a.iter().enumerate().filter(|(_, a)| !base_zero(a)) {
            for (j, b) in b.iter().enumerate().filter(|(_, b)| !base_zero(b)) {
                let term = self.base.multiply(model, a, b)?;
                terms[i + j] = self.base.add(model, &terms[i + j], &term)?;
            }
        }
        for k in (3..5).rev() {
            if !base_zero(&terms[k]) {
                for (j, c) in self.relation.iter().enumerate() {
                    let term = self.base.multiply(model, &terms[k], c)?;
                    terms[k - 3 + j] = self.base.add(model, &terms[k - 3 + j], &term)?;
                }
            }
        }
        let [a, b, c, _, _] = terms;
        Ok([a, b, c])
    }

    fn frobenius(&self, model: &mut Model, power: usize, x: &Coefficient) -> Result<Coefficient> {
        let power = power % 3;
        if power == 0 {
            return Ok(x.clone());
        }
        let [theta, theta_sq] = &self.images[power - 1];
        let linear = self.scale(model, theta, &x[1])?;
        let quadratic = self.scale(model, theta_sq, &x[2])?;
        let mut result = self.add(model, &linear, &quadratic)?;
        result[0] = self.base.add(model, &result[0], &x[0])?;
        Ok(result)
    }

    fn scalar(&self, n: u64) -> Coefficient {
        [self.base.scalar(n), self.base.zero(), self.base.zero()]
    }

    fn theta(&self) -> Coefficient {
        [self.base.zero(), self.base.scalar(1), self.base.zero()]
    }

    fn check(&self, model: &mut Model) -> Result<()> {
        // X^3-X-1 has no F_3 root, and the indicated integral generator
        // and its Frobenius are compatible with the original Kummer chart.
        for digit in 0i32..3 {
            if (digit.pow(3) - digit - 1).rem_euclid(3) == 0 {
                return Err(inconsistent());
            }
        }
        let theta = self.theta();
        let sigma = self.frobenius(model, 1, &theta)?;
        let sigma2 = self.frobenius(model, 1, &sigma)?;
        if sigma2 != self.frobenius(model, 2, &theta)?
            || self.frobenius(model, 1, &sigma2)? != theta
        {
            return Err(inconsistent());
        }
        for image in [&sigma, &sigma2] {
            let squared = self.multiply(model, image, image)?;
            let cubed = self.multiply(model, &squared, image)?;
            let linear = self.scale(model, image, &self.relation[1])?;
            let quadratic = self.scale(model, &squared, &self.relation[2])?;
            let tail = self.add(model, &linear, &quadratic)?;
            if cubed != self.add(model, &self.scalar(1), &tail)? {
                return Err(inconsistent());
            }
        }
        let lambda = self.base.pi();
        let zeta = self.base.add(model, &self.base.scalar(1), &lambda)?;
        let tail = self.scale(model, &theta, &lambda)?;
        let y = self.add(model, &self.scalar(1), &tail)?;
        let y2 = self.multiply(model, &y, &y)?;
        let y3 = self.multiply(model, &y, &y2)?;
        let lambda3 = self.base.power(model, lambda, 3)?;
        let unit3 = self.base.add(model, &self.base.scalar(1), &lambda3)?;
        if y3 != [unit3, self.base.zero(), self.base.zero()]
            || self.frobenius(model, 1, &y)? != self.scale(model, &y, &zeta)?
        {
            return Err(inconsistent());
        }
        Ok(())
    }
}

pub(crate) struct Extension {
    pub(crate) base_name: Box<str>,
    pub(crate) precision: usize,
    pub(crate) base: Rc<Base>,
    invariant: usize,
    frobenius: usize,
}

impl Extension {
    pub(crate) fn new(
        model: &mut Model,
        base_name: Box<str>,
        base: Rc<Base>,
        invariant: usize,
        precision: usize,
    ) -> Result<Self> {
        cubic::require_field(&base)?;
        if !matches!(invariant, 1 | 2) {
            return Err(Error::plain(
                "cubic division invariant numerator must be 1 or 2",
            ));
        }
        let coefficient_precision = precision.div_ceil(3);
        base.algebra(coefficient_precision)?;
        let a = Class([0, 0, 0, 1]);
        let b = Class([invariant as u8, 0, 0, 0]);
        if usize::from(model.cubic_certificate()?.exponent(a, b)) != invariant {
            return Err(inconsistent());
        }
        Unramified::new(model, &base, base.degree() * base.base.precision())?.check(model)?;
        model.allocate(24)?;
        Ok(Self {
            base_name,
            precision,
            base,
            invariant,
            frobenius: invariant, // Both nonzero classes are self-inverse in F_3.
        })
    }

    pub(crate) fn ring(&self, model: &mut Model) -> Result<Order<'_>> {
        Ok(Order {
            extension: self,
            coefficients: Unramified::new(model, &self.base, self.precision.div_ceil(3))?,
        })
    }

    pub(crate) fn same_presentation(&self, source: &Self) -> Result<()> {
        if self.base_name != source.base_name
            || self.invariant != source.invariant
            || self.precision > source.precision
        {
            return Err(Error::plain(
                "cubic division reduction requires the same named base and invariant with sufficient source precision",
            ));
        }
        Ok(())
    }

    pub(crate) fn symbol_compatible(&self, source: &cubic_algebra::Extension) -> Result<()> {
        let required = self
            .precision
            .div_ceil(3)
            .max(cubic::CLASS_PRECISION + self.invariant);
        if source.base_name != self.base_name || source.precision < required {
            return Err(Error::plain(
                "symbol map requires the same named base and sufficient source precision",
            ));
        }
        if usize::from(source.exponent()) != self.invariant {
            return Err(Error::plain(
                "symbol map requires the canonical division invariant",
            ));
        }
        Ok(())
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str) -> Result<String> {
        let mut preceding = Natural::unit();
        for _ in 1..self.precision {
            let mut next = Natural::default();
            next.add_multiple(model, &preceding, 27)?;
            preceding = next;
        }
        let mut size = Natural::default();
        size.add_multiple(model, &preceding, 27)?;
        let mut units = Natural::default();
        units.add_multiple(model, &preceding, 26)?;
        let bridge = if self.invariant == 1 {
            "v"
        } else {
            "v^2/lambda"
        };
        let precision = self.precision;
        Ok(format!(
            "rvCubicDivData[{name}] = base={}, invariant={}/3, frobenius_power={}, ramification_index=3, residue_degree=3, residue_size=27, precision={precision}, coefficient_precisions={:?}, order_size={}, units={}, Pi={bridge}",
            self.base_name, self.invariant, self.frobenius, [0, 1, 2].map(|j| self.digits(j)),
            size.decimal(), units.decimal()
        ))
    }

    pub(crate) fn norm_limit(&self) -> usize {
        self.precision
            .min(self.base.degree() * self.base.base.precision())
    }

    fn digits(&self, j: usize) -> usize {
        self.precision.saturating_sub(j).div_ceil(3)
    }
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Element([Coefficient; 3]);

impl Element {
    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(45)?;
        Ok(Value::CubicDiv(Rc::new(self)))
    }

    pub(crate) fn rows(&self) -> [[u64; 2]; 9] {
        std::array::from_fn(|k| coordinates(&self.0[k / 3][k % 3]))
    }
}

pub(crate) struct Order<'a> {
    extension: &'a Extension,
    coefficients: Unramified<'a>,
}

impl Order<'_> {
    pub(crate) fn zero(&self) -> Element {
        Element([(); 3].map(|()| zero(&self.coefficients.base)))
    }

    pub(super) fn unit(&self, n: u64) -> Element {
        let mut result = self.zero();
        result.0[0] = self.coefficients.scalar(n);
        result
    }

    pub(crate) fn generator(&self, generator: Generator) -> Element {
        let mut result = self.zero();
        match generator {
            Generator::Pi if self.extension.precision > 1 => {
                result.0[1] = self.coefficients.scalar(1);
            }
            Generator::Pi => (),
            Generator::Theta => result.0[0] = self.coefficients.theta(),
        }
        result
    }

    pub(crate) fn literal(&self, model: &mut Model, rows: Vec<Vec<i128>>) -> Result<Element> {
        if rows.len() > 9 {
            return Err(Error::plain(
                "cubic division literals accept up to nine basis coefficients",
            ));
        }
        let mut result = self.zero();
        for (slot, row) in result.0.iter_mut().flatten().zip(rows) {
            *slot = self
                .coefficients
                .base
                .literal(model, row.into_iter().map(|x| vec![x]).collect())?;
        }
        self.reduce(result)
    }

    pub(crate) fn reduce(&self, x: Element) -> Result<Element> {
        self.truncate(x, self.extension.precision)
    }

    pub(super) fn truncate(&self, mut x: Element, precision: usize) -> Result<Element> {
        for (j, w) in x.0.iter_mut().enumerate() {
            let digits = precision.saturating_sub(j).div_ceil(3);
            if digits == 0 {
                *w = zero(&self.coefficients.base);
            } else {
                let base = self.extension.base.algebra(digits)?;
                for c in w {
                    *c = base.reduce(c.clone());
                }
            }
        }
        Ok(x)
    }

    pub(crate) fn add(&self, model: &mut Model, a: &Element, b: &Element) -> Result<Element> {
        let mut result = self.zero();
        for (out, (a, b)) in result.0.iter_mut().zip(a.0.iter().zip(&b.0)) {
            *out = self.coefficients.add(model, a, b)?;
        }
        self.reduce(result)
    }

    pub(crate) fn negate(&self, model: &mut Model, x: &Element) -> Result<Element> {
        let mut result = self.zero();
        for (out, w) in result.0.iter_mut().zip(&x.0) {
            *out = self.coefficients.negate(model, w)?;
        }
        self.reduce(result)
    }

    pub(crate) fn multiply(&self, model: &mut Model, a: &Element, b: &Element) -> Result<Element> {
        model.spend(9)?;
        let mut result = self.zero();
        let w = &self.coefficients;
        let twist = self.extension.frobenius;
        let lambda = w.base.pi();
        for (i, a) in a.0.iter().enumerate().filter(|(_, a)| !coefficient_zero(a)) {
            for (j, b) in b.0.iter().enumerate().filter(|(_, b)| !coefficient_zero(b)) {
                let sigma_b = w.frobenius(model, twist * i, b)?;
                let mut product = w.multiply(model, a, &sigma_b)?;
                if i + j >= 3 {
                    product = w.scale(model, &product, &lambda)?;
                }
                let out = &mut result.0[(i + j) % 3];
                *out = w.add(model, out, &product)?;
            }
        }
        self.reduce(result)
    }

    fn power(&self, model: &mut Model, x: &Element, mut n: usize) -> Result<Element> {
        let mut base = x.clone();
        let mut power = self.unit(1);
        while n != 0 {
            if !n.is_multiple_of(2) {
                power = self.multiply(model, &power, &base)?;
            }
            n /= 2;
            if n != 0 {
                base = self.multiply(model, &base, &base)?;
            }
        }
        Ok(power)
    }

    pub(crate) fn is_unit(&self, x: &Element) -> bool {
        residue(&x.0[0]) != [0, 0, 0]
    }

    pub(crate) fn inverse(&self, model: &mut Model, x: &Element) -> Result<Element> {
        if !self.is_unit(x) {
            return Err(Error::plain(
                "cubic division inverse in the maximal order requires a unit",
            ));
        }
        let one = self.unit(1);
        let two = self.unit(2);
        let mut seed = self.zero();
        let residue = self.extension.base.algebra(1)?;
        seed.0[0] = x.0[0].clone().map(|c| residue.reduce(c));
        let mut inverse = self.power(model, &seed, 25)?;
        let mut precision = 1;
        while precision < self.extension.precision {
            let product = self.multiply(model, x, &inverse)?;
            let negative = self.negate(model, &product)?;
            let correction = self.add(model, &two, &negative)?;
            inverse = self.multiply(model, &inverse, &correction)?;
            precision *= 2;
        }
        if self.multiply(model, x, &inverse)? != one || self.multiply(model, &inverse, x)? != one {
            return Err(inconsistent());
        }
        Ok(inverse)
    }

    pub(crate) fn conjugate_pi(&self, model: &mut Model, x: &Element) -> Result<Element> {
        self.conjugate_pi_power(model, 1, x)
    }

    pub(super) fn conjugate_pi_power(
        &self,
        model: &mut Model,
        power: usize,
        x: &Element,
    ) -> Result<Element> {
        let mut result = self.zero();
        for (out, coefficient) in result.0.iter_mut().zip(&x.0) {
            *out = self.coefficients.frobenius(
                model,
                self.extension.frobenius * power,
                coefficient,
            )?;
        }
        self.reduce(result)
    }

    pub(super) fn divide_right_pi(&self, model: &mut Model, x: &Element) -> Result<Element> {
        if self.is_unit(x) {
            return Err(Error::plain(
                "right division by Pi requires an element of the maximal ideal",
            ));
        }
        if self.extension.precision <= 1 {
            return Ok(self.zero());
        }
        let digits = self.extension.digits(0);
        let mut wrapped = zero(&self.coefficients.base);
        if digits > 1 {
            for (out, c) in wrapped.iter_mut().zip(&x.0[0]) {
                *out = divide_by_lambda(model, &self.extension.base, digits, c.clone(), 1)?;
            }
        }
        let shifted = Element([x.0[1].clone(), x.0[2].clone(), wrapped]);
        self.truncate(shifted, self.extension.precision - 1)
    }

    pub(super) fn graded_digit(
        &self,
        model: &mut Model,
        grade: usize,
        digit: usize,
    ) -> Result<Element> {
        let mut out = self.zero();
        if grade >= self.extension.precision || digit == 0 {
            return Ok(out);
        }
        let w = &self.coefficients;
        let d = [digit / 9, (digit / 3) % 3, digit % 3].map(|a| w.base.scalar(a as u64));
        let lambda = w.base.power(model, w.base.pi(), (grade / 3) as u64)?;
        out.0[grade % 3] = w.scale(model, &d, &lambda)?;
        self.reduce(out)
    }

    pub(crate) fn commutator(
        &self,
        model: &mut Model,
        a: &Element,
        b: &Element,
    ) -> Result<Element> {
        let x = self.inverse(model, a)?;
        let y = self.inverse(model, b)?;
        let ab = self.multiply(model, a, b)?;
        let abx = self.multiply(model, &ab, &x)?;
        self.multiply(model, &abx, &y)
    }

    fn valuation(&self, model: &mut Model, x: &Element) -> Result<usize> {
        let mut order = self.extension.precision;
        for (j, row) in x.0.iter().enumerate() {
            for c in row {
                order = order.min(3 * self.coefficients.base.valuation(model, c)? + j);
            }
        }
        Ok(order)
    }

    fn leading(&self, model: &mut Model, x: &Element, val: usize) -> Result<Coefficient> {
        if val >= self.extension.precision {
            return Err(inconsistent());
        }
        let mut leading = zero(&self.coefficients.base);
        let digits = self.extension.digits(val % 3);
        let residue = self.extension.base.algebra(1)?;
        for (out, c) in leading.iter_mut().zip(&x.0[val % 3]) {
            let divided =
                divide_by_lambda(model, &self.extension.base, digits, c.clone(), val / 3)?;
            *out = residue.reduce(divided);
        }
        Ok(leading)
    }

    fn norm_precision(&self, valuation: usize) -> usize {
        (valuation + (self.extension.precision - valuation).div_ceil(3))
            .min(self.extension.norm_limit())
    }

    pub(crate) fn norm(
        &self,
        model: &mut Model,
        x: &Element,
        digits: usize,
    ) -> Result<BaseElement> {
        let val = self.valuation(model, x)?;
        let available = self.norm_precision(val);
        if digits == 0 || digits > available {
            return Err(Error::plain(format!(
                "reduced norm requires positive lambda precision at most {available} for this cubic division residue"
            )));
        }
        self.determinant(model, x, digits, val)
    }

    fn determinant(
        &self,
        model: &mut Model,
        x: &Element,
        digits: usize,
        val: usize,
    ) -> Result<BaseElement> {
        let w = Unramified::new(model, &self.extension.base, digits)?;
        let lambda = w.base.pi();
        let mut matrix = [[(); 3]; 3].map(|row| row.map(|()| zero(&w.base)));
        // Left multiplication on the right W basis 1, Pi, Pi^2.
        for (j, row) in x.0.iter().enumerate() {
            let row = row.clone().map(|c| w.base.reduce(c));
            for (target, matrix_row) in matrix.iter_mut().enumerate() {
                let column = (target + 3 - j) % 3;
                let power = self.extension.frobenius * (3 - target);
                let mut coefficient = w.frobenius(model, power, &row)?;
                if j + column >= 3 {
                    coefficient = w.scale(model, &coefficient, &lambda)?;
                }
                matrix_row[column] = coefficient;
            }
        }
        let mut norm = zero(&w.base);
        for [a, b, c] in [
            [0, 1, 2],
            [1, 2, 0],
            [2, 0, 1],
            [0, 2, 1],
            [1, 0, 2],
            [2, 1, 0],
        ] {
            let pair = w.multiply(model, &matrix[0][a], &matrix[1][b])?;
            let term = w.multiply(model, &pair, &matrix[2][c])?;
            let term = if (a + 1) % 3 == b {
                term
            } else {
                w.negate(model, &term)?
            };
            norm = w.add(model, &norm, &term)?;
        }
        let [norm, linear, square] = norm;
        if !base_zero(&linear)
            || !base_zero(&square)
            || w.base.valuation(model, &norm)? != val.min(digits)
        {
            return Err(inconsistent());
        }
        Ok(norm)
    }

    pub(crate) fn valuation_report(&self, model: &mut Model, x: &Element) -> Result<String> {
        let val = self.valuation(model, x)?;
        let digits = self.norm_precision(val);
        let norm = coordinates(&self.determinant(model, x, digits, val)?);
        let leading = if val == self.extension.precision {
            format!("at_least={val}, leading=unresolved")
        } else {
            format!(
                "valuation={val}, leading={:?}",
                residue(&self.leading(model, x, val)?)
            )
        };
        Ok(format!(
            "rvCubicDivVal = {leading}, norm_precision={digits}, reduced_norm={norm:?}"
        ))
    }

    pub(crate) fn commutator_report(
        &self,
        model: &mut Model,
        a: &Element,
        b: &Element,
    ) -> Result<String> {
        let minus_one = self.negate(model, &self.unit(1))?;
        let x = self.add(model, a, &minus_one)?;
        let y = self.add(model, b, &minus_one)?;
        let i = self.valuation(model, &x)?;
        let j = self.valuation(model, &y)?;
        if i == 0 || j == 0 {
            return Err(Error::plain(
                "graded commutator check requires principal units",
            ));
        }
        let k = i + j;
        let commutator = self.commutator(model, a, b)?;
        let delta = self.add(model, &commutator, &minus_one)?;
        let val = self.valuation(model, &delta)?;
        if k >= self.extension.precision {
            if val != self.extension.precision {
                return Err(inconsistent());
            }
            return Ok(format!(
                "rvCubicDivCommutatorCheck = layers=[{i}, {j}], at_least={val}, grade=unresolved"
            ));
        }
        let w = Unramified::new(model, &self.extension.base, 1)?;
        let a = self.leading(model, &x, i)?;
        let b = self.leading(model, &y, j)?;
        let twist = self.extension.frobenius;
        let sigma_b = w.frobenius(model, twist * i, &b)?;
        let sigma_a = w.frobenius(model, twist * j, &a)?;
        let ab = w.multiply(model, &a, &sigma_b)?;
        let ba = w.multiply(model, &b, &sigma_a)?;
        let negative = w.negate(model, &ba)?;
        let expected = residue(&w.add(model, &ab, &negative)?);
        let actual = if val == k {
            residue(&self.leading(model, &delta, k)?)
        } else {
            [0; 3]
        };
        if val < k || actual != expected {
            return Err(inconsistent());
        }
        let filtration = if val == self.extension.precision {
            format!("at_least={val}")
        } else {
            format!("valuation={val}")
        };
        Ok(format!("rvCubicDivCommutatorCheck = layers=[{i}, {j}], grade={k}, leading={actual:?}, {filtration}"))
    }

    pub(crate) fn embed_symbol(
        &self,
        model: &mut Model,
        source: &cubic_algebra::Extension,
        x: &cubic_algebra::Element,
    ) -> Result<Element> {
        self.extension.symbol_compatible(source)?;
        let source_ring = self.extension.base.algebra(source.precision)?;
        let lambda = source_ring.pi();
        let lambda3 = source_ring.power(model, lambda.clone(), 3)?;
        let required_a = source_ring.add(model, &source_ring.scalar(1), &lambda3)?;
        let required_b = source_ring.power(model, lambda, self.extension.invariant as u64)?;
        let [a, b] = source.parameters();
        if *a != required_a || *b != required_b {
            return Err(Error::plain(
                "symbol map requires parameters U_3 and lambda^t at the source chart precision",
            ));
        }
        let base = &self.coefficients.base;
        let lambda = base.pi();
        let one = self.unit(1);
        let mut y = one.clone();
        y.0[0][1] = lambda;
        let v = self.power(
            model,
            &self.generator(Generator::Pi),
            self.extension.invariant,
        )?;
        let ys = [one.clone(), y.clone(), self.multiply(model, &y, &y)?];
        let vs = [one, v.clone(), self.multiply(model, &v, &v)?];
        let mut result = self.zero();
        for (i, y) in ys.iter().enumerate() {
            for (j, v) in vs.iter().enumerate() {
                let coefficient = base.reduce(x.coefficient(3 * i + j));
                let mut basis = self.multiply(model, y, v)?;
                for row in &mut basis.0 {
                    *row = self.coefficients.scale(model, row, &coefficient)?;
                }
                result = self.add(model, &result, &basis)?;
            }
        }
        Ok(result)
    }
}

fn zero(base: &BaseAlgebra<'_>) -> Coefficient {
    [(); 3].map(|()| base.zero())
}

fn coordinates(x: &BaseElement) -> [u64; 2] {
    [x.0[0][0], x.0[1][0]]
}

fn base_zero(x: &BaseElement) -> bool {
    coordinates(x) == [0, 0]
}

fn coefficient_zero(x: &Coefficient) -> bool {
    x.iter().all(base_zero)
}

fn residue(x: &Coefficient) -> [u64; 3] {
    x.each_ref().map(|c| c.0[0][0] % 3)
}

fn inconsistent() -> Error {
    Error::plain("inconsistent cubic division valuation, norm or ramification identity")
}
