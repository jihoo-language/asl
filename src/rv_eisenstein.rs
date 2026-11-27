use std::rc::Rc;

use crate::ec_poly::{self as poly, Field as Coefficients};
use crate::error::{Error, Result};
use crate::hw_algebra::{self, Matrix};
use crate::model::{Model, Value};
use crate::rv_ore::gcd;
use crate::rv_unram::{mul_mod, reduce, sub_mod, Ring};
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod kummer;

const MAX_DEGREE: usize = 32;

pub(crate) struct Spec(pub(crate) Label, pub(crate) usize);

pub(crate) enum Declaration {
    Polynomial(Label, Label, Vec<i128>),
    Cyclotomic(Label, Label, usize),
}

pub(crate) enum Expression {
    Literal(Spec, Vec<Vec<i128>>),
    Uniformizer(Spec),
    Multiply(Box<[ParsedExpr; 2]>),
    Evaluate(Box<[ParsedExpr; 2]>),
    Automorphism(i128, Box<ParsedExpr>),
}

pub(crate) enum Unary {
    Negate,
    Automorphism(u64),
}

pub(crate) enum ExtensionReport {
    Polynomial(Label),
    Polygon(Label),
    Herbrand(Spec),
}

pub(crate) enum ElementReport {
    Valuation,
    Trace,
    Norm,
}

impl ElementReport {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            Self::Valuation => "rvEisensteinVal",
            Self::Trace => "rvEisensteinTrace",
            Self::Norm => "rvEisensteinNorm",
        }
    }
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) name: Box<str>,
    pub(crate) precision: usize,
}

impl Grade {
    pub(crate) fn resolve(
        model: &Model,
        Spec(name, precision): Spec,
    ) -> Result<(Self, Rc<Extension>)> {
        let extension = model.eisenstein(&name)?;
        extension.algebra(precision)?;
        Ok((
            Self {
                name: name.text.into(),
                precision,
            },
            extension,
        ))
    }
}

pub(crate) struct Extension {
    base: Rc<Ring>,
    exact: Vec<i128>,
    coefficients: Vec<u64>,
    cyclotomic_order: Option<u64>,
}

impl Extension {
    pub(crate) fn new(
        model: &mut Model,
        base: Rc<Ring>,
        polynomial: Vec<i128>,
        level: Option<usize>,
    ) -> Result<Self> {
        let (exact, cyclotomic_order) = if let Some(level) = level {
            let (polynomial, order) = cyclotomic_polynomial(model, base.prime(), level)?;
            (polynomial, Some(order))
        } else {
            (polynomial, None)
        };
        let Some((&1, lower)) = exact.split_last() else {
            return Err(Error::plain("Eisenstein polynomial must be monic"));
        };
        if !(1..=MAX_DEGREE).contains(&lower.len()) {
            return Err(degree_limit());
        }
        let p = i128::from(base.prime());
        if lower.iter().any(|&a| a % p != 0) || lower[0] % (p * p) == 0 {
            return Err(Error::plain("Eisenstein coefficients require p to divide every lower coefficient and v_p(a_0)=1"));
        }
        let coefficients = exact.iter().map(|&a| reduce(a, base.modulus())).collect();
        model.allocate(4 + exact.len() * 4)?;
        Ok(Self {
            base,
            exact,
            coefficients,
            cyclotomic_order,
        })
    }

    pub(crate) fn prime(&self) -> u64 {
        self.base.prime()
    }

    pub(crate) fn degree(&self) -> usize {
        self.exact.len() - 1
    }

    pub(crate) fn algebra(&self, precision: usize) -> Result<Algebra<'_>> {
        let e = self.degree();
        if precision == 0 || precision > e * self.base.precision() {
            return Err(Error::plain("Eisenstein uniformizer precision must be positive and supported by the named unramified ring"));
        }
        let moduli = (0..e)
            .map(|j| {
                self.prime()
                    .pow(precision.saturating_sub(j).div_ceil(e) as u32)
            })
            .collect();
        Ok(Algebra {
            extension: self,
            precision,
            moduli,
        })
    }

    pub(crate) fn automorphism_exponent(&self, exponent: i128) -> Result<u64> {
        let order = self.cyclotomic_order.ok_or_else(|| {
            Error::plain("cyclotomic automorphisms require a named p-power cyclotomic extension")
        })?;
        let residue = exponent.rem_euclid(i128::from(order)) as u64;
        if residue.is_multiple_of(self.prime()) {
            return Err(Error::plain(
                "cyclotomic automorphism exponent must be a unit modulo p^n",
            ));
        }
        Ok(residue)
    }

    pub(crate) fn polynomial_report(&self, name: &str) -> String {
        format!(
            "rvEisensteinPolynomial[{name}] = p={}, degree={}, coefficients={:?}",
            self.prime(),
            self.degree(),
            self.exact
        )
    }

    pub(crate) fn polygon_report(&self, model: &mut Model, name: &str) -> Result<String> {
        let polygon = self.ramification_polygon(model)?;
        Ok(format!("rvRamificationPolygon[{name}] = valuations={:?}, vertices={:?}, embedding_distances={:?}, different={}", polygon.valuations, polygon.vertices, polygon.distances, polygon.different))
    }

    fn ramification_polygon(&self, model: &mut Model) -> Result<Polygon> {
        let p = self.prime();
        let e = self.degree();
        model.spend(e * e)?;
        let valuations: Vec<_> = (1..=e)
            .map(|j| {
                (j..=e)
                    .filter(|&i| self.exact[i] != 0)
                    .map(|i| {
                        e * (order(self.exact[i].unsigned_abs(), p) + order(binomial(i, j), p)) + i
                            - e
                    })
                    .min()
                    .expect("monic term in every ramification coefficient")
            })
            .collect();
        let different = (1..=e)
            .filter(|&i| self.exact[i] != 0)
            .map(|i| e * (order(i as u128, p) + order(self.exact[i].unsigned_abs(), p)) + i - 1)
            .min()
            .expect("nonzero monic derivative");
        if different != valuations[0] + e - 1 {
            return Err(inconsistent());
        }
        let check_precision = different + 1;
        if check_precision <= e * self.base.precision() {
            let algebra = self.algebra(check_precision)?;
            let derivative: Vec<_> = self
                .coefficients
                .iter()
                .enumerate()
                .skip(1)
                .map(|(i, &c)| mul_mod(i as u64, c, self.base.modulus()))
                .collect();
            let value = algebra.evaluate_rational(model, &derivative, &algebra.pi())?;
            if algebra.valuation(model, &value)? != different {
                return Err(inconsistent());
            }
        }
        let mut vertices = Vec::<(usize, usize)>::new();
        for (j, &v) in (1..=e).zip(&valuations) {
            while let [.., (a, b), (c, d)] = vertices.as_slice() {
                let left = (*d as isize - *b as isize) * (j - c) as isize;
                let right = (v as isize - *d as isize) * (c - a) as isize;
                if left < right {
                    break;
                }
                vertices.pop();
            }
            vertices.push((j, v));
        }
        let distances = vertices
            .windows(2)
            .rev()
            .map(|edge| {
                let [(a, b), (c, d)] = [edge[0], edge[1]];
                let width = c - a;
                let (n, den) = ratio(width + b - d, width);
                (n, den, width)
            })
            .collect();
        Ok(Polygon {
            valuations,
            vertices,
            distances,
            different,
        })
    }

    pub(crate) fn herbrand_report(
        &self,
        model: &mut Model,
        name: &str,
        precision: usize,
    ) -> Result<String> {
        let order_of_cyclotomic_root = self.cyclotomic_order.ok_or_else(|| {
            Error::plain("Hasse–Herbrand filtration requires a named p-power cyclotomic extension")
        })?;
        let algebra = self.algebra(precision)?;
        let p = self.prime();
        let pi = algebra.pi();
        let negative_pi = algebra.negate(model, &pi)?;
        let mut distances = Vec::new();
        let mut unresolved = Vec::new();
        for a in (2..order_of_cyclotomic_root).filter(|a| !a.is_multiple_of(p)) {
            let image = algebra.image_of_pi(model, a)?;
            let difference = algebra.add(model, &image, &negative_pi)?;
            let observed = algebra.valuation(model, &difference)?;
            let expected = p.pow(order(u128::from(a - 1), p) as u32) as usize;
            if observed != expected.min(precision) {
                return Err(inconsistent());
            }
            if observed == precision {
                unresolved.push(a);
            } else {
                distances.push((a, observed));
            }
        }
        if !unresolved.is_empty() {
            return Ok(format!("rvHerbrand[{name},{precision}] = distances={distances:?}, valuation_at_least={precision}, unresolved={unresolved:?}"));
        }
        let last = distances.iter().map(|(_, v)| *v).max().unwrap_or(0);
        let polygon = self.ramification_polygon(model)?;
        let groups: Vec<_> = (0..=last)
            .map(|i| {
                let mut group = vec![1];
                group.extend(distances.iter().filter(|(_, v)| *v > i).map(|(a, _)| *a));
                group
            })
            .collect();
        let sizes: Vec<_> = groups.iter().map(Vec::len).collect();
        let different = sizes.iter().map(|&size| size - 1).sum::<usize>();
        let mut polygon_distances = vec![0; last + 1];
        for (n, d, count) in polygon.distances {
            if d != 1 || n > last {
                return Err(inconsistent());
            }
            polygon_distances[n] += count;
        }
        for &(_, v) in &distances {
            polygon_distances[v] = polygon_distances[v]
                .checked_sub(1)
                .ok_or_else(inconsistent)?;
        }
        if different != polygon.different || polygon_distances.iter().any(|&n| n != 0) {
            return Err(inconsistent());
        }

        let mut upper_numerator = 0;
        let mut lower_breaks = Vec::new();
        let mut upper_breaks = Vec::new();
        for i in 0..last {
            if i != 0 {
                upper_numerator += sizes[i];
            }
            if sizes[i] > sizes[i + 1] {
                lower_breaks.push(i);
                upper_breaks.push(ratio(upper_numerator, self.degree()));
            }
        }
        Ok(format!("rvHerbrand[{name},{precision}] = distances={distances:?}, lower_groups={groups:?}, sizes={sizes:?}, lower_breaks={lower_breaks:?}, upper_breaks={upper_breaks:?}, different={different}"))
    }
}

struct Polygon {
    valuations: Vec<usize>,
    vertices: Vec<(usize, usize)>,
    distances: Vec<(usize, usize, usize)>,
    different: usize,
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Element(Vec<Vec<u64>>);

impl Element {
    pub(crate) fn rows(&self) -> &[Vec<u64>] {
        &self.0
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + self.0.len() * (1 + self.0[0].len()))?;
        Ok(Value::Eisenstein(Rc::new(self)))
    }
}

pub(crate) struct Algebra<'a> {
    extension: &'a Extension,
    precision: usize,
    moduli: Vec<u64>,
}

impl Algebra<'_> {
    pub(crate) fn zero(&self) -> Element {
        Element(vec![self.extension.base.scalar(0); self.extension.degree()])
    }

    pub(crate) fn pi(&self) -> Element {
        let mut result = self.zero();
        let extension = self.extension;
        if extension.degree() == 1 {
            result.0[0] = extension.base.scalar(sub_mod(
                0,
                extension.coefficients[0],
                extension.base.modulus(),
            ));
        } else {
            result.0[1] = extension.base.scalar(1);
        }
        self.reduce(result)
    }

    fn scalar(&self, value: u64) -> Element {
        let mut result = self.zero();
        result.0[0] = self.extension.base.scalar(value);
        self.reduce(result)
    }

    pub(crate) fn literal(
        &self,
        model: &mut Model,
        coefficients: Vec<Vec<i128>>,
    ) -> Result<Element> {
        if coefficients.len() > self.extension.degree() {
            return Err(Error::plain(
                "Eisenstein literal has more rows than the extension degree",
            ));
        }
        let mut output = self.zero();
        for (row, coefficients) in output.0.iter_mut().zip(coefficients) {
            *row = self.extension.base.literal(model, coefficients)?;
        }
        Ok(self.reduce(output))
    }

    pub(crate) fn add(&self, model: &mut Model, a: &Element, b: &Element) -> Result<Element> {
        let ring = &self.extension.base;
        let rows =
            a.0.iter()
                .zip(&b.0)
                .map(|(x, y)| ring.add(model, x, y))
                .collect::<Result<_>>()?;
        Ok(self.reduce(Element(rows)))
    }

    pub(crate) fn negate(&self, model: &mut Model, value: &Element) -> Result<Element> {
        let rows = value
            .0
            .iter()
            .map(|row| self.extension.base.negate(model, row))
            .collect::<Result<_>>()?;
        Ok(self.reduce(Element(rows)))
    }

    pub(crate) fn multiply(&self, model: &mut Model, a: &Element, b: &Element) -> Result<Element> {
        let e = self.extension.degree();
        let ring = &self.extension.base;
        model.spend(e * e * ring.degree())?;
        let length = (2 * e - 1).min(self.precision).max(e);
        let mut coefficients = vec![ring.scalar(0); length];
        for (i, x) in a.0.iter().enumerate().filter(|(_, x)| !poly::zero(x)) {
            for (j, y) in
                b.0.iter()
                    .enumerate()
                    .take(length - i)
                    .filter(|(_, y)| !poly::zero(y))
            {
                let product = ring.multiply(model, x, y)?;
                coefficients[i + j] = ring.add(model, &coefficients[i + j], &product)?;
            }
        }
        // Reduction by an Eisenstein relation cannot lower a term's
        // uniformizer valuation, so truncation above is valid.
        for k in (e..length).rev() {
            let leading = coefficients.pop().expect("current polynomial degree");
            if poly::zero(&leading) {
                continue;
            }
            for (i, &c) in self.extension.coefficients[..e]
                .iter()
                .enumerate()
                .filter(|(_, c)| **c != 0)
            {
                let product = ring.scale(model, &leading, ring.modulus() - c)?;
                coefficients[k - e + i] = ring.add(model, &coefficients[k - e + i], &product)?;
            }
        }
        Ok(self.reduce(Element(coefficients)))
    }

    pub(crate) fn evaluate_rational(
        &self,
        model: &mut Model,
        coefficients: &[u64],
        point: &Element,
    ) -> Result<Element> {
        let mut output = self.zero();
        for &coefficient in coefficients.iter().rev() {
            let tail = self.multiply(model, &output, point)?;
            output = self.add(model, &tail, &self.scalar(coefficient))?;
        }
        Ok(output)
    }

    fn image_of_pi(&self, model: &mut Model, exponent: u64) -> Result<Element> {
        let one = self.scalar(1);
        let generator = self.add(model, &one, &self.pi())?;
        let power = self.power(model, generator, exponent)?;
        let negative_one = self.negate(model, &one)?;
        self.add(model, &power, &negative_one)
    }

    fn power(&self, model: &mut Model, mut base: Element, mut exponent: u64) -> Result<Element> {
        let mut value = self.scalar(1);
        while exponent != 0 {
            if !exponent.is_multiple_of(2) {
                value = self.multiply(model, &value, &base)?;
            }
            exponent /= 2;
            if exponent != 0 {
                base = self.multiply(model, &base, &base)?;
            }
        }
        Ok(value)
    }

    pub(crate) fn automorphism(
        &self,
        model: &mut Model,
        value: &Element,
        exponent: u64,
    ) -> Result<Element> {
        let point = self.image_of_pi(model, exponent)?;
        let mut output = self.zero();
        for row in value.0.iter().rev() {
            let mut constant = self.zero();
            constant.0[0].clone_from(row);
            let tail = self.multiply(model, &output, &point)?;
            output = self.add(model, &tail, &constant)?;
        }
        Ok(output)
    }

    fn valuation(&self, model: &mut Model, value: &Element) -> Result<usize> {
        let e = self.extension.degree();
        let mut valuation = self.precision;
        for (j, c) in value.0.iter().enumerate() {
            valuation = valuation.min(e * self.extension.base.valuation(model, c)? + j);
        }
        Ok(valuation)
    }

    pub(crate) fn report(
        &self,
        model: &mut Model,
        value: &Element,
        report: ElementReport,
    ) -> Result<String> {
        let name = report.name();
        if matches!(report, ElementReport::Valuation) {
            let valuation = self.valuation(model, value)?;
            return Ok(if valuation == self.precision {
                format!("{name} = at_least={valuation}")
            } else {
                format!("{name} = {valuation}")
            });
        }
        let e = self.extension.degree();
        let precision = self.precision / e;
        if precision == 0 {
            return Ok(format!("{name} = precision=0, unresolved"));
        }
        let ring = &self.extension.base;
        let modulus = ring.prime().pow(precision as u32);
        let pi = self.pi();
        let mut column = value.0.clone();
        let mut trace = ring.scalar(0);
        let mut matrix = matches!(report, ElementReport::Norm).then(|| vec![ring.scalar(0); e * e]);
        for j in 0..e {
            trace = ring.add(model, &trace, &column[j])?;
            if let Some(matrix) = &mut matrix {
                for (i, c) in column.iter().enumerate() {
                    matrix[i * e + j].clone_from(c);
                }
            }
            if j + 1 < e {
                column = self.multiply(model, &Element(column), &pi)?.0;
            }
        }
        let coordinates = if let Some(matrix) = matrix {
            let matrix = Matrix::from_entries(e, matrix)?;
            let characteristic =
                hw_algebra::characteristic(model, &Coefficients { ring }, &matrix)?;
            let determinant = characteristic
                .last()
                .expect("characteristic polynomial with constant slot");
            if e.is_multiple_of(2) {
                determinant.clone()
            } else {
                ring.negate(model, determinant)?
            }
        } else {
            trace
        };
        let coordinates: Vec<_> = coordinates.into_iter().map(|x| x % modulus).collect();
        Ok(format!(
            "{name} = precision={precision}, modulus={modulus}, coordinates={coordinates:?}"
        ))
    }

    fn reduce(&self, mut value: Element) -> Element {
        for (row, &modulus) in value.0.iter_mut().zip(&self.moduli) {
            for c in row {
                *c %= modulus;
            }
        }
        value
    }
}

fn cyclotomic_polynomial(model: &mut Model, p: u64, level: usize) -> Result<(Vec<i128>, u64)> {
    if level == 0 {
        return Err(Error::plain("p-power cyclotomic level must be positive"));
    }
    let exponent = u32::try_from(level - 1).map_err(|_| degree_limit())?;
    let step = p.checked_pow(exponent).ok_or_else(degree_limit)?;
    let degree = (p - 1)
        .checked_mul(step)
        .filter(|&e| e <= MAX_DEGREE as u64)
        .ok_or_else(degree_limit)?;
    let mut result = vec![0; degree as usize + 1];
    for j in 0..p {
        let power = (j * step) as usize;
        model.spend(power + 1)?;
        for (k, coefficient) in result.iter_mut().enumerate().take(power + 1) {
            *coefficient += binomial(power, k) as i128;
        }
    }
    Ok((result, p * step))
}

fn binomial(n: usize, k: usize) -> u128 {
    let mut coefficient = 1;
    for j in 0..k.min(n - k) {
        coefficient = coefficient * (n - j) as u128 / (j + 1) as u128;
    }
    coefficient
}

fn order(mut value: u128, p: u64) -> usize {
    debug_assert_ne!(value, 0);
    let p = u128::from(p);
    let mut exponent = 0;
    while value.is_multiple_of(p) {
        value /= p;
        exponent += 1;
    }
    exponent
}

fn ratio(n: usize, d: usize) -> (usize, usize) {
    let gcd = gcd(n as u128, d as u128) as usize;
    (n / gcd, d / gcd)
}

fn degree_limit() -> Error {
    Error::plain(format!(
        "Eisenstein extension degree must be in 1..={MAX_DEGREE}"
    ))
}

fn inconsistent() -> Error {
    Error::plain("inconsistent Eisenstein ramification identity")
}
