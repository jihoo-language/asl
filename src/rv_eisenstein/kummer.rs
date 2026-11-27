use std::rc::Rc;

use crate::ec_poly::{self as poly, Field};
use crate::error::{Error, Result};
use crate::finite_field_factor::factor_with_multiplicity;
use crate::hw_algebra::{self, Matrix};
use crate::model::{Model, Value};
use crate::rv_unram::sub_mod;
use crate::syntax::{Label, ParsedExpr};

use super::{binomial, Algebra as BaseAlgebra, Element as BaseElement, Extension as BaseExtension};

const MAX_P: usize = 7;

pub(crate) struct Declaration(
    pub(crate) Label,
    pub(crate) Label,
    pub(crate) usize,
    pub(crate) Vec<Vec<i128>>,
);
pub(crate) struct Spec(pub(crate) Label, pub(crate) usize);

pub(crate) enum Generator {
    Y,
    Theta,
    Uniformizer,
}

pub(crate) enum Expression {
    Literal(Spec, Vec<Vec<Vec<i128>>>),
    Generator(Spec, Generator),
    Base(Spec, Box<ParsedExpr>),
    Multiply(Box<[ParsedExpr; 2]>),
    Evaluate(Box<[ParsedExpr; 2]>),
    Sigma(i128, Box<ParsedExpr>),
}

pub(crate) enum Unary {
    Negate,
    Sigma(usize),
}

pub(crate) enum Report {
    Valuation,
    Trace,
    Norm,
}

impl Report {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            Self::Valuation => "rvKummerLocalVal",
            Self::Trace => "rvKummerLocalTrace",
            Self::Norm => "rvKummerLocalNorm",
        }
    }
}

pub(crate) enum FieldReport {
    Data(Spec),
    Polynomial(Spec),
    Boundary(Label),
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub(crate) struct Grade {
    pub(crate) name: Box<str>,
    pub(crate) precision: usize,
}

impl Grade {
    pub(crate) fn resolve(
        model: &Model,
        Spec(name, precision): Spec,
    ) -> Result<(Self, Rc<Extension>)> {
        let field = model.local_kummer(&name)?;
        field.work_precision(precision)?;
        Ok((
            Self {
                name: name.text.into(),
                precision,
            },
            field,
        ))
    }
}

pub(crate) struct Extension {
    pub(crate) base_name: Box<str>,
    base: Rc<BaseExtension>,
    p: usize,
    b: usize,
    unit: BaseElement,
    denominators: Vec<usize>,
    orders: Vec<usize>,
    clearing: usize,
    e: usize,
}

impl Extension {
    pub(crate) fn new(
        model: &mut Model,
        base_name: Box<str>,
        base: Rc<BaseExtension>,
        b: usize,
        unit: Vec<Vec<i128>>,
    ) -> Result<Self> {
        let p = base.prime() as usize;
        if p == 2 || p > MAX_P || base.cyclotomic_order != Some(p as u64) {
            return Err(Error::plain(format!("relative Kummer base must be a named zeta_p extension for an odd prime up to {MAX_P}")));
        }
        if b == 0 || b > p {
            return Err(Error::plain("relative Kummer unit level must lie in 1..=p"));
        }
        let full = base.algebra(base.degree() * base.base.precision())?;
        let unit = full.literal(model, unit)?;
        if full.valuation(model, &unit)? != 0 {
            return Err(Error::plain("relative Kummer parameter u must be a unit"));
        }
        let denominators: Vec<_> = (0..p).map(|i| b * i / p).collect();
        let orders = (0..p).map(|i| b * i % p).collect();
        let clearing = denominators[p - 1];
        let e = if b == p { 1 } else { p };
        model.allocate(8 + p * 2 + base.degree() * base.base.degree())?;
        Ok(Self {
            base_name,
            base,
            p,
            b,
            unit,
            denominators,
            orders,
            clearing,
            e,
        })
    }

    pub(crate) fn prime(&self) -> usize {
        self.p
    }

    pub(crate) fn absolute_e(&self) -> usize {
        self.e * self.base.degree()
    }

    pub(crate) fn base_digits(&self, precision: usize) -> usize {
        precision.div_ceil(self.e)
    }

    fn work_precision(&self, precision: usize) -> Result<usize> {
        let cap = self.base.degree() * self.base.base.precision();
        self.base_digits(precision)
            .checked_add(2 * self.clearing)
            .filter(|&working| precision != 0 && working <= cap)
            .ok_or_else(precision_error)
    }

    pub(crate) fn algebra(&self, model: &mut Model, precision: usize) -> Result<Algebra<'_>> {
        let working = self.work_precision(precision)?;
        let base = self.base.algebra(working)?;
        let mut lambda = vec![base.scalar(1)];
        let generator = base.pi();
        for _ in 1..=self.b.max(2 * self.clearing) {
            lambda.push(base.multiply(model, lambda.last().expect("initial power"), &generator)?);
        }
        // The coefficient a in y^p=a is represented by w^p + sum binom(p,i)w^i + 1.
        let mut relation =
            vec![base.multiply(model, &base.reduce(self.unit.clone()), &lambda[self.b])?];
        for i in 1..self.p {
            relation.push(base.negate(model, &base.scalar(binomial(self.p, i) as u64))?);
        }
        Ok(Algebra {
            field: self,
            precision,
            base,
            lambda,
            relation,
        })
    }

    pub(crate) fn sigma_exponent(&self, c: i128) -> usize {
        c.rem_euclid(self.p as i128) as usize
    }

    fn uniformizer_index(&self) -> Result<usize> {
        if self.b == self.p {
            return Err(Error::plain(
                "relative Kummer uniformizer requires a totally ramified unit level b<p",
            ));
        }
        Ok(self
            .orders
            .iter()
            .position(|&r| r == 1)
            .expect("orders are all residues"))
    }

    pub(crate) fn boundary_report(&self, model: &mut Model, name: &str) -> Result<String> {
        if self.b != self.p {
            return Err(Error::plain("Artin–Schreier boundary report requires b=p"));
        }
        let ring = self.base.base.residue_ring(model)?;
        let residue = ring.residue(&self.unit.0[0]);
        let trace = ring.conjugate_sum_or_product(model, &residue, false)?[0];
        let mut equation = vec![ring.scalar(0); self.p + 1];
        equation[0] = ring.negate(model, &residue)?;
        equation[1] = ring.scalar(self.p as u64 - 1);
        equation[self.p] = ring.scalar(1);
        let factors = factor_with_multiplicity(model, Field { ring: &ring }, &equation)?;
        let (kind, roots) = if trace == 0 {
            if factors.len() != self.p
                || factors
                    .iter()
                    .any(|(degree, mult, _)| (*degree, *mult) != (1, 1))
            {
                return Err(inconsistent());
            }
            let roots = factors
                .iter()
                .map(|(_, _, q)| ring.negate(model, &q[0]))
                .collect::<Result<Vec<_>>>()?;
            ("kind=split".to_owned(), roots)
        } else {
            if factors.len() != 1 || (factors[0].0, factors[0].1) != (self.p, 1) {
                return Err(inconsistent());
            }
            (
                format!("kind=unramified, residue_degree={}, different=0", self.p),
                Vec::new(),
            )
        };
        let verification = self.boundary_lifts(model, &equation, &roots)?;
        Ok(format!("rvKummerLocalBoundary[{name}] = trace={trace}, {kind}, residue_polynomial={equation:?}, residual_roots={roots:?}{verification}"))
    }

    fn boundary_lifts(
        &self,
        model: &mut Model,
        residue: &[Vec<u64>],
        roots: &[Vec<u64>],
    ) -> Result<String> {
        let available = self.base.degree() * self.base.base.precision();
        let precision = available - (self.p - 1);
        if precision == 0 {
            return Ok(", theta_precision=0".to_owned());
        }
        let full = self.base.algebra(available)?;
        let base = self.base.algebra(precision)?;
        let mut equation = vec![base.negate(model, &base.reduce(self.unit.clone()))?];
        for i in 1..self.p {
            let coefficient = full.scalar(binomial(self.p, i) as u64);
            let divided = divide_by_lambda(model, &self.base, available, coefficient, self.p - i)?;
            equation.push(base.reduce(divided));
        }
        equation.push(base.scalar(1));
        for (coefficient, expected) in equation.iter().zip(residue) {
            if self.base.base.residue(&coefficient.0[0]) != *expected {
                return Err(inconsistent());
            }
        }
        let derivative: Vec<_> = equation
            .iter()
            .enumerate()
            .skip(1)
            .map(|(j, coefficient)| base.multiply(model, coefficient, &base.scalar(j as u64)))
            .collect::<Result<_>>()?;
        let mut lifts = Vec::new();
        for digits in roots {
            let mut root = base.zero();
            root.0[0].clone_from(digits);
            let mut known = 1;
            while known < precision {
                let error = evaluate_base(model, &base, &equation, &root)?;
                let slope = evaluate_base(model, &base, &derivative, &root)?;
                let inverse = inverse_base_unit(model, &base, &slope)?;
                let correction = base.multiply(model, &error, &inverse)?;
                let negative = base.negate(model, &correction)?;
                root = base.add(model, &root, &negative)?;
                known *= 2;
            }
            if !zero(&evaluate_base(model, &base, &equation, &root)?) {
                return Err(inconsistent());
            }
            lifts.push(root);
        }
        let equation: Vec<_> = equation.iter().map(BaseElement::rows).collect();
        let lifts: Vec<_> = lifts.iter().map(BaseElement::rows).collect();
        Ok(format!(
            ", theta_precision={precision}, theta_polynomial={equation:?}, split_lifts={lifts:?}"
        ))
    }

    pub(crate) fn data_report(
        &self,
        model: &mut Model,
        name: &str,
        precision: usize,
    ) -> Result<String> {
        let i = self.uniformizer_index()?;
        let algebra = self.algebra(model, precision)?;
        let q = self.denominators[i];
        let pi = algebra.basis(i);
        let w = algebra.w();
        if algebra.valuation(model, &w)? != self.b.min(precision)
            || algebra.valuation(model, &pi)? != 1.min(precision)
        {
            return Err(inconsistent());
        }
        let power = algebra.power(model, &w, i)?;
        let base_lambda = algebra.embed(algebra.lambda[q].clone());
        if algebra.multiply(model, &pi, &base_lambda)? != power {
            return Err(inconsistent());
        }
        let expected = 1 + self.p - self.b;
        let mut valuations = Vec::new();
        let negative = algebra.negate(model, &pi)?;
        for c in 1..self.p {
            let image = algebra.sigma(model, &pi, c)?;
            let delta = algebra.add(model, &image, &negative)?;
            let measured = algebra.valuation(model, &delta)?;
            if measured != expected.min(precision) {
                return Err(inconsistent());
            }
            valuations.push(measured);
        }
        let info = format!("rvKummerLocalData[{name},{precision}] = b={}, w_order={}, uniformizer=(w^{i})/lambda^{q}", self.b, self.b);
        if precision <= expected {
            return Ok(format!(
                "{info}, distances_at_least={precision}, unresolved"
            ));
        }
        let lower = expected - 1;
        let relative = (self.p - 1) * expected;
        let absolute = self.p * (self.p - 2) + relative;
        Ok(format!("{info}, distances={valuations:?}, lower_break={lower}, upper_break={lower}, relative_different={relative}, absolute_different={absolute}"))
    }

    pub(crate) fn polynomial_report(
        &self,
        model: &mut Model,
        name: &str,
        precision: usize,
    ) -> Result<String> {
        let index = self.uniformizer_index()?;
        let algebra = self.algebra(model, precision)?;
        let degree = self.absolute_e();
        let digits = precision / degree;
        if digits == 0 {
            return Ok(format!(
                "rvKummerLocalPolynomial[{name},{precision}] = precision=0, unresolved"
            ));
        }
        let modulus = (self.p as u64).pow(digits as u32);
        let ring = &self.base.base;
        let pi = algebra.basis(index);
        let basis_width = self.base.degree();
        let mut entries = vec![ring.scalar(0); degree * degree];
        let lambda = algebra.embed(algebra.base.pi());
        for i in 0..self.p {
            let mut value = algebra.multiply(model, &pi, &algebra.basis(i))?;
            for j in 0..basis_width {
                let column = i * basis_width + j;
                for (ii, coordinate) in value.0.iter().enumerate() {
                    for (jj, row) in coordinate.0.iter().enumerate() {
                        entries[(ii * basis_width + jj) * degree + column].clone_from(row);
                    }
                }
                if j + 1 < basis_width {
                    value = algebra.multiply(model, &value, &lambda)?;
                }
            }
        }
        let matrix = Matrix::from_entries(degree, entries)?;
        let mut polynomial = hw_algebra::characteristic(model, &Field { ring }, &matrix)?;
        polynomial.reverse();
        for coordinate in &mut polynomial {
            for c in coordinate {
                *c %= modulus;
            }
        }
        let relative = (self.p - 1) * (1 + self.p - self.b);
        let expected = self.p * (self.p - 2) + relative;
        let certified = digits * degree;
        let value = algebra.evaluate_unramified(model, &polynomial, &pi)?;
        if algebra.valuation(model, &value)? < certified {
            return Err(inconsistent());
        }
        let derivative = if certified > expected {
            let derivative: Vec<_> = polynomial
                .iter()
                .enumerate()
                .skip(1)
                .map(|(i, coefficient)| ring.scale(model, coefficient, i as u64))
                .collect::<Result<_>>()?;
            let value = algebra.evaluate_unramified(model, &derivative, &pi)?;
            let valuation = algebra.valuation(model, &value)?;
            if valuation != expected {
                return Err(inconsistent());
            }
            format!(", derivative_valuation={valuation}")
        } else {
            String::new()
        };
        let trace: Vec<_> = polynomial[degree - 1]
            .iter()
            .map(|&c| sub_mod(0, c, modulus))
            .collect();
        let norm = &polynomial[0]; // p(p-1) is even.
        let text = if ring.degree() == 1 {
            format!(
                "{:?}",
                polynomial.iter().map(|row| row[0]).collect::<Vec<_>>()
            )
        } else {
            format!("{polynomial:?}")
        };
        Ok(format!("rvKummerLocalPolynomial[{name},{precision}] = degree={degree}, precision={digits}, modulus={modulus}, coefficients={text}, base_trace={trace:?}, base_norm={norm:?}{derivative}"))
    }
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Element(Vec<BaseElement>);

impl Element {
    pub(crate) fn rows(&self) -> Vec<&[Vec<u64>]> {
        self.0.iter().map(BaseElement::rows).collect()
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + self.0.len() * (1 + self.0[0].0.len() * (1 + self.0[0].0[0].len())))?;
        Ok(Value::LocalKummer(Rc::new(self)))
    }
}

pub(crate) struct Algebra<'a> {
    field: &'a Extension,
    precision: usize,
    base: BaseAlgebra<'a>,
    lambda: Vec<BaseElement>,
    relation: Vec<BaseElement>,
}

impl Algebra<'_> {
    pub(crate) fn zero(&self) -> Element {
        Element(vec![self.base.zero(); self.field.p])
    }

    fn basis(&self, i: usize) -> Element {
        let mut value = self.zero();
        value.0[i] = self.base.scalar(1);
        self.reduce(value)
    }

    pub(crate) fn embed(&self, base: BaseElement) -> Element {
        let mut value = self.zero();
        value.0[0] = base;
        self.reduce(value)
    }

    pub(crate) fn literal(&self, model: &mut Model, rows: Vec<Vec<Vec<i128>>>) -> Result<Element> {
        if rows.len() > self.field.p {
            return Err(Error::plain("too many relative Kummer basis rows"));
        }
        let mut value = self.zero();
        for (target, row) in value.0.iter_mut().zip(rows) {
            *target = self.base.literal(model, row)?;
        }
        Ok(self.reduce(value))
    }

    pub(crate) fn generator(&self, model: &mut Model, generator: Generator) -> Result<Element> {
        match generator {
            Generator::Y => self.add(model, &self.basis(0), &self.w()),
            Generator::Uniformizer => Ok(self.basis(self.field.uniformizer_index()?)),
            Generator::Theta if self.field.b == self.field.p => Ok(self.basis(1)),
            Generator::Theta => Err(Error::plain("relative Kummer theta generator requires b=p")),
        }
    }

    fn w(&self) -> Element {
        let mut value = self.zero();
        value.0[1] = self.lambda[self.field.denominators[1]].clone();
        self.reduce(value)
    }

    pub(crate) fn add(
        &self,
        model: &mut Model,
        left: &Element,
        right: &Element,
    ) -> Result<Element> {
        let values = left
            .0
            .iter()
            .zip(&right.0)
            .map(|(a, b)| self.base.add(model, a, b))
            .collect::<Result<_>>()?;
        Ok(self.reduce(Element(values)))
    }

    pub(crate) fn negate(&self, model: &mut Model, value: &Element) -> Result<Element> {
        let values = value
            .0
            .iter()
            .map(|a| self.base.negate(model, a))
            .collect::<Result<_>>()?;
        Ok(self.reduce(Element(values)))
    }

    pub(crate) fn multiply(
        &self,
        model: &mut Model,
        left: &Element,
        right: &Element,
    ) -> Result<Element> {
        let left = self.clear_denominators(model, left)?;
        let right = self.clear_denominators(model, right)?;
        let product = self.ambient_product(model, &left, &right)?;
        self.restore_lattice(model, product, 2 * self.field.clearing)
    }

    fn ambient_product(
        &self,
        model: &mut Model,
        left: &[BaseElement],
        right: &[BaseElement],
    ) -> Result<Vec<BaseElement>> {
        let p = self.field.p;
        let mut product = vec![self.base.zero(); 2 * p - 1];
        for (i, a) in left.iter().enumerate().filter(|(_, a)| !zero(a)) {
            for (j, b) in right.iter().enumerate().filter(|(_, b)| !zero(b)) {
                let term = self.base.multiply(model, a, b)?;
                product[i + j] = self.base.add(model, &product[i + j], &term)?;
            }
        }
        for k in (p..product.len()).rev() {
            let leading = product.pop().expect("ambient product degree");
            if zero(&leading) {
                continue;
            }
            for (i, c) in self.relation.iter().enumerate().filter(|(_, c)| !zero(c)) {
                let term = self.base.multiply(model, &leading, c)?;
                product[k - p + i] = self.base.add(model, &product[k - p + i], &term)?;
            }
        }
        Ok(product)
    }

    fn clear_denominators(&self, model: &mut Model, value: &Element) -> Result<Vec<BaseElement>> {
        value
            .0
            .iter()
            .zip(&self.field.denominators)
            .map(|(a, q)| {
                self.base
                    .multiply(model, a, &self.lambda[self.field.clearing - q])
            })
            .collect()
    }

    fn restore_lattice(
        &self,
        model: &mut Model,
        value: Vec<BaseElement>,
        clearing: usize,
    ) -> Result<Element> {
        let output = value
            .into_iter()
            .zip(&self.field.denominators)
            .map(|(a, q)| {
                divide_by_lambda(
                    model,
                    &self.field.base,
                    self.base.precision,
                    a,
                    clearing - q,
                )
            })
            .collect::<Result<_>>()?;
        Ok(self.reduce(Element(output)))
    }

    pub(crate) fn sigma(&self, model: &mut Model, value: &Element, c: usize) -> Result<Element> {
        if c == 0 {
            return Ok(Element(value.0.clone()));
        }
        let zeta = self
            .base
            .add(model, &self.base.scalar(1), &self.base.pi())?;
        let power = self.base.power(model, zeta, c as u64)?;
        let mut point = vec![self.base.zero(); self.field.p];
        let negative_one = self.base.negate(model, &self.base.scalar(1))?;
        point[0] = self.base.add(model, &power, &negative_one)?;
        point[1] = power;
        let clear = self.clear_denominators(model, value)?;
        let mut output = vec![self.base.zero(); self.field.p];
        for coefficient in clear.iter().rev() {
            output = self.ambient_product(model, &output, &point)?;
            output[0] = self.base.add(model, &output[0], coefficient)?;
        }
        self.restore_lattice(model, output, self.field.clearing)
    }

    fn power(&self, model: &mut Model, value: &Element, power: usize) -> Result<Element> {
        let mut result = self.basis(0);
        for _ in 0..power {
            result = self.multiply(model, &result, value)?;
        }
        Ok(result)
    }

    pub(crate) fn evaluate_rational(
        &self,
        model: &mut Model,
        values: &[u64],
        point: &Element,
    ) -> Result<Element> {
        let values: Vec<_> = values
            .iter()
            .map(|&a| self.field.base.base.scalar(a))
            .collect();
        self.evaluate_unramified(model, &values, point)
    }

    fn evaluate_unramified(
        &self,
        model: &mut Model,
        values: &[Vec<u64>],
        point: &Element,
    ) -> Result<Element> {
        let mut result = self.zero();
        for coefficient in values.iter().rev() {
            let mut constant = self.base.zero();
            constant.0[0].clone_from(coefficient);
            let tail = self.multiply(model, &result, point)?;
            result = self.add(model, &tail, &self.embed(constant))?;
        }
        Ok(result)
    }

    fn valuation(&self, model: &mut Model, value: &Element) -> Result<usize> {
        let mut order = self.precision;
        for (a, &r) in value.0.iter().zip(&self.field.orders) {
            order = order.min(self.field.e * self.base.valuation(model, a)? + r);
        }
        Ok(order)
    }

    pub(crate) fn report(
        &self,
        model: &mut Model,
        value: &Element,
        report: Report,
    ) -> Result<String> {
        let name = report.name();
        if matches!(report, Report::Valuation) {
            let order = self.valuation(model, value)?;
            return Ok(if order == self.precision {
                format!("{name} = at_least={order}")
            } else {
                format!("{name} = {order}")
            });
        }
        let result = self.conjugate_invariant(model, value, matches!(report, Report::Norm))?;
        Ok(format!(
            "{name} = base={}, lambda_precision={}, coordinates={:?}",
            self.field.base_name,
            self.field.base_digits(self.precision),
            result.0
        ))
    }

    pub(super) fn conjugate_invariant(
        &self,
        model: &mut Model,
        value: &Element,
        product: bool,
    ) -> Result<BaseElement> {
        let mut result = if product { self.basis(0) } else { self.zero() };
        for c in 0..self.field.p {
            let image = self.sigma(model, value, c)?;
            result = if product {
                self.multiply(model, &result, &image)?
            } else {
                self.add(model, &result, &image)?
            };
        }
        if result.0[1..].iter().any(|value| !zero(value)) {
            return Err(inconsistent());
        }
        Ok(result.0.remove(0))
    }

    fn reduce(&self, mut value: Element) -> Element {
        for (coefficient, &order) in value.0.iter_mut().zip(&self.field.orders) {
            let digits = self.precision.saturating_sub(order).div_ceil(self.field.e);
            let e = self.field.base.degree();
            let p = self.field.p as u64;
            for (j, row) in coefficient.0.iter_mut().enumerate() {
                let modulus = p.pow(digits.saturating_sub(j).div_ceil(e) as u32);
                for c in row {
                    *c %= modulus;
                }
            }
        }
        value
    }
}

// The cyclotomic equation gives p/lambda = -lambda^(e-1) - sum a_i lambda^(i-1).
// A division is integral precisely when the constant row is divisible by p.
pub(super) fn divide_by_lambda(
    model: &mut Model,
    field: &BaseExtension,
    precision: usize,
    mut value: BaseElement,
    count: usize,
) -> Result<BaseElement> {
    let ring = &field.base;
    let p = ring.prime();
    let e = field.degree();
    for lost in 0..count {
        if value.0[0].iter().any(|x| !x.is_multiple_of(p)) {
            return Err(inconsistent());
        }
        let quotient: Vec<_> = value.0[0].iter().map(|x| x / p).collect();
        let mut rows = Vec::with_capacity(e);
        for i in 1..e {
            let c = sub_mod(0, field.coefficients[i], ring.modulus());
            let term = ring.scale(model, &quotient, c)?;
            rows.push(ring.add(model, &value.0[i], &term)?);
        }
        rows.push(ring.negate(model, &quotient)?);
        value = field
            .algebra(precision - lost - 1)?
            .reduce(BaseElement(rows));
    }
    Ok(value)
}

fn evaluate_base(
    model: &mut Model,
    base: &BaseAlgebra<'_>,
    coefficients: &[BaseElement],
    point: &BaseElement,
) -> Result<BaseElement> {
    let mut result = base.zero();
    for coefficient in coefficients.iter().rev() {
        let tail = base.multiply(model, &result, point)?;
        result = base.add(model, &tail, coefficient)?;
    }
    Ok(result)
}

pub(super) fn inverse_base_unit(
    model: &mut Model,
    base: &BaseAlgebra<'_>,
    unit: &BaseElement,
) -> Result<BaseElement> {
    let mut inverse = base.zero();
    inverse.0[0] = base.extension.base.inverse(model, &unit.0[0])?;
    let mut approximation = base.reduce(inverse);
    let two = base.scalar(2);
    let mut digits = 1;
    while digits < base.precision {
        let product = base.multiply(model, unit, &approximation)?;
        let negative = base.negate(model, &product)?;
        let correction = base.add(model, &two, &negative)?;
        approximation = base.multiply(model, &approximation, &correction)?;
        digits *= 2;
    }
    Ok(approximation)
}

fn zero(a: &BaseElement) -> bool {
    a.0.iter().all(|row| poly::zero(row))
}

fn precision_error() -> Error {
    Error::plain("relative Kummer precision exceeds the base digits after clearing denominators")
}

fn inconsistent() -> Error {
    Error::plain("inconsistent relative Kummer lattice or ramification identity")
}
