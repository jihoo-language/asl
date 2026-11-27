use std::rc::Rc;

use crate::ec_poly::{self as poly, Field};
use crate::error::{Error, Result};
use crate::finite_field_factor::factor_with_multiplicity;
use crate::model::{Model, Value};
use crate::rv_unram::Ring;
use crate::syntax::{Label, ParsedExpr};

use super::{newton_polygon, order, Preparation};

const MAX_E: usize = 32;

pub(crate) struct Spec(pub(crate) Label, pub(crate) usize, pub(crate) usize);

pub(crate) enum Expression {
    Literal(Spec, Vec<Vec<i128>>),
    Teichmuller(Spec, Box<ParsedExpr>),
    Multiply(Box<[ParsedExpr; 2]>),
    Evaluate(Box<[ParsedExpr; 2]>),
    Frobenius(usize, Box<ParsedExpr>),
    Inertia(Box<ParsedExpr>),
}

pub(crate) enum Unary {
    Negate,
    Teichmuller,
    Frobenius(usize),
    Inertia,
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) name: Box<str>,
    pub(crate) e: usize,
    pub(crate) precision: usize,
}

impl Grade {
    pub(crate) fn resolve(
        model: &Model,
        Spec(name, e, precision): Spec,
    ) -> Result<(Self, Rc<Ring>)> {
        let ring = model.rv_unram(&name)?;
        Algebra::new(&ring, e, precision)?;
        Ok((
            Self {
                name: name.text.into(),
                e,
                precision,
            },
            ring,
        ))
    }

    pub(crate) fn algebra<'a>(&self, ring: &'a Ring) -> Result<Algebra<'a>> {
        Algebra::new(ring, self.e, self.precision)
    }
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Element(Vec<Vec<u64>>);

impl Element {
    pub(crate) fn rows(&self) -> &[Vec<u64>] {
        &self.0
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + self.0.len() * (1 + self.0[0].len()))?;
        Ok(Value::Tame(Rc::new(self)))
    }
}

pub(crate) struct Algebra<'a> {
    ring: &'a Ring,
    e: usize,
    precision: usize,
    moduli: Vec<u64>,
}

impl<'a> Algebra<'a> {
    fn new(ring: &'a Ring, e: usize, precision: usize) -> Result<Self> {
        if e == 0 || e > MAX_E {
            return Err(Error::plain(format!(
                "tame ramification index must be in 1..={MAX_E}"
            )));
        }
        if e.is_multiple_of(ring.prime() as usize) {
            return Err(Error::plain(
                "ramification index divisible by the residue characteristic is wild",
            ));
        }
        if precision == 0 || precision > e * ring.precision() {
            return Err(Error::plain("tame uniformizer precision must be positive and supported by the named unramified ring"));
        }
        let moduli = (0..e)
            .map(|j| {
                ring.prime()
                    .pow(precision.saturating_sub(j).div_ceil(e) as u32)
            })
            .collect();
        Ok(Self {
            ring,
            e,
            precision,
            moduli,
        })
    }

    pub(crate) fn zero(&self) -> Element {
        Element(vec![self.ring.scalar(0); self.e])
    }

    pub(crate) fn literal(
        &self,
        model: &mut Model,
        coefficients: Vec<Vec<i128>>,
    ) -> Result<Element> {
        if coefficients.len() > self.e {
            return Err(Error::plain(
                "tame literal has more rows than the ramification index",
            ));
        }
        let mut value = self.zero();
        for (row, coefficients) in value.0.iter_mut().zip(coefficients) {
            *row = self.ring.literal(model, coefficients)?;
        }
        Ok(self.reduce(value))
    }

    pub(crate) fn teichmuller(&self, model: &mut Model, value: &[u64]) -> Result<Element> {
        Ok(self.constant(self.ring.teichmuller(model, value)?))
    }

    fn constant(&self, row: Vec<u64>) -> Element {
        let mut value = self.zero();
        value.0[0] = row;
        self.reduce(value)
    }

    pub(crate) fn add(&self, model: &mut Model, a: &Element, b: &Element) -> Result<Element> {
        let rows =
            a.0.iter()
                .zip(&b.0)
                .map(|(x, y)| self.ring.add(model, x, y))
                .collect::<Result<_>>()?;
        Ok(self.reduce(Element(rows)))
    }

    pub(crate) fn negate(&self, model: &mut Model, a: &Element) -> Result<Element> {
        self.scale(model, a, self.ring.modulus() - 1)
    }

    fn scale(&self, model: &mut Model, a: &Element, scalar: u64) -> Result<Element> {
        let rows =
            a.0.iter()
                .map(|x| self.ring.scale(model, x, scalar))
                .collect::<Result<_>>()?;
        Ok(self.reduce(Element(rows)))
    }

    pub(crate) fn multiply(&self, model: &mut Model, a: &Element, b: &Element) -> Result<Element> {
        model.spend(self.e * self.e * self.ring.degree())?;
        let mut answer = self.zero();
        for (i, x) in a.0.iter().enumerate().filter(|(_, x)| !poly::zero(x)) {
            for (j, y) in b.0.iter().enumerate().filter(|(_, y)| !poly::zero(y)) {
                let degree = i + j;
                let row = degree % self.e;
                if degree >= self.precision {
                    continue;
                }
                let mut product = self.ring.multiply(model, x, y)?;
                if degree >= self.e {
                    product = self.ring.scale(model, &product, self.ring.prime())?;
                }
                answer.0[row] = self.ring.add(model, &answer.0[row], &product)?;
            }
        }
        Ok(self.reduce(answer))
    }

    fn monomial(&self, model: &mut Model, coefficient: u64, exponent: usize) -> Result<Element> {
        let mut result = self.zero();
        if exponent < self.precision {
            let power = self.ring.prime().pow((exponent / self.e) as u32);
            result.0[exponent % self.e] =
                self.ring
                    .scale(model, &self.ring.scalar(coefficient), power)?;
        }
        Ok(self.reduce(result))
    }

    fn shift(&self, model: &mut Model, value: &Element, exponent: usize) -> Result<Element> {
        let mut result = self.zero();
        for (j, row) in value.0.iter().enumerate() {
            let degree = j + exponent;
            if degree < self.precision {
                let power = self.ring.prime().pow((degree / self.e) as u32);
                result.0[degree % self.e] = self.ring.scale(model, row, power)?;
            }
        }
        Ok(self.reduce(result))
    }

    pub(crate) fn evaluate_rational(
        &self,
        model: &mut Model,
        polynomial: &[u64],
        point: &Element,
    ) -> Result<Element> {
        let mut output = self.zero();
        for &coefficient in polynomial.iter().rev() {
            let next = self.multiply(model, &output, point)?;
            output = self.add(model, &next, &self.constant(self.ring.scalar(coefficient)))?;
        }
        Ok(output)
    }

    fn evaluate(
        &self,
        model: &mut Model,
        polynomial: &[Element],
        point: &Element,
    ) -> Result<Element> {
        let mut output = self.zero();
        for coefficient in polynomial.iter().rev() {
            let next = self.multiply(model, &output, point)?;
            output = self.add(model, &next, coefficient)?;
        }
        Ok(output)
    }

    fn inverse(&self, model: &mut Model, unit: &Element) -> Result<Element> {
        let mut inverse = self.constant(self.ring.inverse(model, &unit.0[0])?);
        let mut known = 1;
        let two = self.constant(self.ring.scalar(2));
        while known < self.precision {
            let product = self.multiply(model, unit, &inverse)?;
            let negative = self.negate(model, &product)?;
            let correction = self.add(model, &two, &negative)?;
            inverse = self.multiply(model, &inverse, &correction)?;
            known *= 2;
        }
        Ok(inverse)
    }

    pub(crate) fn frobenius(
        &self,
        model: &mut Model,
        value: &Element,
        exponent: usize,
    ) -> Result<Element> {
        let rows = value
            .0
            .iter()
            .map(|row| self.ring.frobenius(model, row, exponent))
            .collect::<Result<_>>()?;
        Ok(self.reduce(Element(rows)))
    }

    pub(crate) fn inertia(&self, model: &mut Model, value: &Element) -> Result<Element> {
        let root = primitive_unity(model, self.ring, self.e)?.ok_or_else(|| {
            Error::plain("named residue field does not contain the required tame roots of unity")
        })?;
        self.apply_inertia(model, value, &root)
    }

    fn apply_inertia(&self, model: &mut Model, value: &Element, root: &[u64]) -> Result<Element> {
        let mut coefficient = self.ring.scalar(1);
        let mut rows = Vec::with_capacity(self.e);
        for row in &value.0 {
            rows.push(self.ring.multiply(model, row, &coefficient)?);
            coefficient = self.ring.multiply(model, &coefficient, root)?;
        }
        Ok(self.reduce(Element(rows)))
    }

    pub(crate) fn valuation(&self, model: &mut Model, value: &Element) -> Result<usize> {
        let mut best = self.precision;
        for (j, row) in value.0.iter().enumerate() {
            best = best.min(self.e * self.ring.valuation(model, row)? + j);
        }
        Ok(best)
    }

    fn reduce(&self, mut value: Element) -> Element {
        for (row, &modulus) in value.0.iter_mut().zip(&self.moduli) {
            for entry in row {
                *entry %= modulus;
            }
        }
        value
    }
}

fn primitive_unity(model: &mut Model, ring: &Ring, e: usize) -> Result<Option<Vec<u64>>> {
    if e == 1 {
        return Ok(Some(ring.scalar(1)));
    }
    let mut q_mod_e = 1;
    for _ in 0..ring.degree() {
        q_mod_e = (q_mod_e * (ring.prime() % e as u64)) % e as u64;
    }
    if q_mod_e != 1 {
        return Ok(None);
    }
    let finite = ring.residue_ring(model)?;
    // Prefer the named power-basis generator when it already has order e.
    if ring.degree() > 1 {
        let generator = finite.literal(model, vec![0, 1])?;
        if exact_order(model, &finite, &generator, e)? {
            return lifted_unity(model, ring, &generator, e).map(Some);
        }
    }
    let mut equation = vec![finite.scalar(0); e + 1];
    equation[0] = finite.scalar(finite.prime() - 1);
    equation[e] = finite.scalar(1);
    for (degree, _, factor) in factor_with_multiplicity(model, Field { ring: &finite }, &equation)?
    {
        if degree == 1 {
            let root = finite.negate(model, &factor[0])?;
            if exact_order(model, &finite, &root, e)? {
                return lifted_unity(model, ring, &root, e).map(Some);
            }
        }
    }
    Err(Error::plain(
        "tame primitive root-of-unity factorization is inconsistent",
    ))
}

fn lifted_unity(model: &mut Model, ring: &Ring, residue: &[u64], order: usize) -> Result<Vec<u64>> {
    let lift = ring.teichmuller(model, residue)?;
    if ring.power(model, &lift, order as u64)? != ring.scalar(1)
        || ring.frobenius(model, &lift, 1)? != ring.power(model, &lift, ring.prime())?
    {
        return Err(inconsistent());
    }
    Ok(lift)
}

fn exact_order(model: &mut Model, ring: &Ring, value: &[u64], expected: usize) -> Result<bool> {
    let one = ring.scalar(1);
    let mut power = one.clone();
    for exponent in 1..=expected {
        power = ring.multiply(model, &power, value)?;
        if power == one {
            return Ok(exponent == expected);
        }
    }
    Ok(false)
}

pub(super) fn report(
    model: &mut Model,
    prepared: &Preparation,
    name: &str,
    ring: &Ring,
    index: usize,
    header: String,
) -> Result<String> {
    model.spend(prepared.polynomial.len().pow(2))?;
    let (segments, unknown) = newton_polygon(prepared.ring, &prepared.polynomial);
    let Some(&(h, e, multiplicity)) = segments.get(index) else {
        if let Some(floor) = unknown.filter(|_| index == segments.len()) {
            return Ok(format!(
                "{header}, segment={index}, valuation_floor={floor:?}, unresolved"
            ));
        }
        return Err(Error::plain(format!(
            "Newton segment index {index} is out of range"
        )));
    };
    let header =
        format!("{header}, segment={index}, slope=({h}, {e}, {multiplicity}), ring={name}");
    if e.is_multiple_of(ring.prime() as usize) {
        return Ok(format!("{header}, wild, unresolved"));
    }
    if e > MAX_E {
        return Ok(format!(
            "{header}, ramification_index_exceeds={MAX_E}, unresolved"
        ));
    }

    let weighted_orders: Vec<_> = prepared
        .polynomial
        .iter()
        .enumerate()
        .map(|(i, &a)| e * order(a, prepared.ring) + h * i)
        .collect();
    let m = *weighted_orders
        .iter()
        .min()
        .expect("monic polynomial has a nonempty coefficient vector");
    let verification = e * prepared.ring.digits.min(ring.precision());
    if verification <= m {
        return Ok(format!("{header}, normalization={m}, available_uniformizer_precision={verification}, insufficient_precision"));
    }
    let precision = verification - m;
    let algebra = Algebra::new(ring, e, precision)?;
    let alpha = Algebra::new(ring, e, precision + h)?;
    let verifier = Algebra::new(ring, e, verification)?;
    let p = ring.prime();
    let equation: Vec<_> = prepared
        .polynomial
        .iter()
        .zip(weighted_orders)
        .map(|(&a, weight)| {
            let unit = a / p.pow(order(a, prepared.ring) as u32);
            algebra.monomial(model, unit, weight - m)
        })
        .collect::<Result<_>>()?;
    let derivative: Vec<_> = equation
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, a)| algebra.scale(model, a, i as u64))
        .collect::<Result<_>>()?;
    let reduction: Vec<_> = equation.iter().map(|x| ring.residue(&x.0[0])).collect();
    let first = reduction
        .iter()
        .position(|x| !poly::zero(x))
        .ok_or_else(inconsistent)?;
    let residual = poly::trim(reduction[first..].to_vec());
    if residual.len() != multiplicity + 1 {
        return Err(inconsistent());
    }

    let finite = ring.residue_ring(model)?;
    let factors = factor_with_multiplicity(model, Field { ring: &finite }, &residual)?;
    let mut y_roots = Vec::new();
    let mut roots = Vec::new();
    let mut repeated = Vec::new();
    let mut outside = Vec::new();
    for (degree, count, factor) in factors {
        if degree != 1 {
            outside.push((factor, count));
            continue;
        }
        if count != 1 {
            repeated.push((factor, count));
            continue;
        }
        let residue = finite.negate(model, &factor[0])?;
        let y = lift(model, &algebra, &equation, &derivative, &residue)?;
        let lifted = alpha.shift(model, &y, h)?;
        let point = verifier.shift(model, &y, h)?;
        let value = verifier.evaluate_rational(model, &prepared.polynomial, &point)?;
        if value != verifier.zero() || alpha.valuation(model, &lifted)? != h {
            return Err(inconsistent());
        }
        y_roots.push(y);
        roots.push(lifted);
    }
    model.spend(roots.len().pow(2) * e * ring.degree())?;
    let f = permutation(&roots, |point| alpha.frobenius(model, point, 1))?;
    let inertia = match primitive_unity(model, ring, e)? {
        None => "inertia=unavailable".to_owned(),
        Some(xi) => {
            model.spend(roots.len().pow(2) * e * ring.degree())?;
            let t = permutation(&roots, |point| alpha.apply_inertia(model, point, &xi))?;
            for j in 0..roots.len() {
                let mut power = f[j];
                for _ in 0..p % e as u64 {
                    power = t[power];
                }
                if f[t[j]] != power {
                    return Err(inconsistent());
                }
            }
            format!("zeta={xi:?}, inertia={t:?}, relation=verified")
        }
    };
    let y_rows: Vec<_> = y_roots.iter().map(Element::rows).collect();
    let root_rows: Vec<_> = roots.iter().map(Element::rows).collect();
    Ok(format!("{header}, normalization={m}, y_precision={precision}, y_moduli={:?}, root_precision={}, root_moduli={:?}, verified_mod_pi={verification}, removed_zero_power={first}, residual={residual:?}, y={y_rows:?}, roots={root_rows:?}, repeated={repeated:?}, outside_field={outside:?}, frobenius={f:?}, {inertia}", algebra.moduli, alpha.precision, alpha.moduli))
}

fn lift(
    model: &mut Model,
    algebra: &Algebra<'_>,
    equation: &[Element],
    derivative: &[Element],
    residue: &[u64],
) -> Result<Element> {
    let mut y = algebra.teichmuller(model, residue)?;
    let mut known = 1;
    while known < algebra.precision {
        let value = algebra.evaluate(model, equation, &y)?;
        let slope = algebra.evaluate(model, derivative, &y)?;
        let inverse = algebra.inverse(model, &slope)?;
        let correction = algebra.multiply(model, &value, &inverse)?;
        let negative = algebra.negate(model, &correction)?;
        y = algebra.add(model, &y, &negative)?;
        known *= 2;
    }
    if algebra.evaluate(model, equation, &y)? != algebra.zero() {
        return Err(inconsistent());
    }
    Ok(y)
}

fn permutation(
    roots: &[Element],
    mut map: impl FnMut(&Element) -> Result<Element>,
) -> Result<Vec<usize>> {
    let mut images = Vec::with_capacity(roots.len());
    let mut used = vec![false; roots.len()];
    for root in roots {
        let transformed = map(root)?;
        let target = roots
            .iter()
            .position(|root| root == &transformed)
            .ok_or_else(inconsistent)?;
        if used[target] {
            return Err(inconsistent());
        }
        used[target] = true;
        images.push(target);
    }
    Ok(images)
}

fn inconsistent() -> Error {
    Error::plain("inconsistent tame Newton lift or Galois permutation")
}
