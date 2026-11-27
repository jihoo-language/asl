use std::rc::Rc;

use super::{horizontal, Certificate};
use crate::ec_morphism::Morphism;
use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::elliptic::{Curve, Isogeny, Point};
use crate::elliptic_cm_recognition::Recognition;
use crate::error::{Error, Result};
use crate::finite_field_factor;
use crate::hw_algebra::CoefficientRing;
use crate::hyperelliptic_kummer_surface::cm_ring::{Order, GROUP_LIMIT};
use crate::legendre_cm::{Element, Family};
use crate::model::{Model, MAX_SHOWN_BYTES};
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Statement {
    Define(Label, horizontal::Address),
    Data(Label),
    Step(Label, usize),
    Verify(Label),
}

pub(crate) enum Expression {
    Apply(Label, Box<ParsedExpr>),
    Element(Label),
}

pub(crate) struct Cycle {
    pub(crate) certificate_name: Label,
    certificate: Rc<Certificate>,
    pub(crate) curve_name: Rc<str>,
    vertices: Vec<Rc<Curve>>,
    steps: Vec<Isogeny>,
    j_path: Vec<Scalar>,
    chart: Morphism,
    scale: Scalar,
    chart_count: usize,
    differential: Scalar,
    degree: u128,
    conjugate: bool,
    pub(crate) element: Element,
    trace: i128,
}

impl Cycle {
    pub(crate) fn new(
        model: &mut Model,
        name: &Label,
        address: horizontal::Address,
    ) -> Result<Self> {
        let action = horizontal::Action::new(model, &address)?;
        let certificate = Rc::clone(&action.certificate);
        let order = Order::new(certificate.discriminant)?;
        let length = order.class_order(model, action.ideal.class)?;
        if length > GROUP_LIMIT {
            return Err(Error::plain(
                "CM isogeny cycle exceeds the finite-value limit",
            ));
        }
        let degree = u128::from(action.ell)
            .checked_pow(u32::try_from(length).map_err(|_| norm_limit())?)
            .ok_or_else(norm_limit)?;
        let four_degree = degree.checked_mul(4).ok_or_else(norm_limit)?;
        let bound = (four_degree / certificate.discriminant.unsigned_abs()).isqrt();
        model
            .spend(
                usize::try_from(bound)
                    .map_err(|_| norm_limit())?
                    .checked_mul(32)
                    .and_then(|work| work.checked_add(32))
                    .ok_or_else(norm_limit)?,
            )
            .map_err(|_| norm_limit())?;

        let source_name: Rc<str> = certificate.curve_name().into();
        let mut vertices = vec![Rc::clone(&action.curve)];
        let mut steps = Vec::with_capacity(length);
        let mut j_path = vec![action.curve.j_invariant(model)?];
        let mut current: Rc<str> = Rc::clone(&source_name);
        for i in 0..length {
            let domain = &vertices[i];
            let kernel = Rc::new(action.kernel_on(model, domain)?);
            let target: Rc<str> = format!("{}::{}", name.text, i + 1).into();
            let (map, quotient) = Isogeny::new(model, domain, current, Rc::clone(&target), kernel)?;
            action.check_quotient(model, domain, &quotient)?;
            j_path.push(quotient.j_invariant(model)?);
            vertices.push(Rc::new(quotient));
            steps.push(map);
            current = target;
        }
        if j_path.first() != j_path.last() {
            return Err(Error::plain(
                "principal CM ideal cycle did not return to its initial j-invariant",
            ));
        }
        let start = &vertices[0];
        let end = &vertices[length];
        let (scale, chart_count) = chart_scale(model, end, start)?;
        let chart = Morphism::chart(model, current, Rc::clone(&source_name), end, start, &scale)?;
        // Normalized Vélu maps pull back the invariant differential
        // unchanged. A short chart (u^2*x,u^3*y) has multiplier 1/u.
        let differential = start.field().inverse_unit(model, &scale)?;
        let element = generator(model, &action, order, &differential, four_degree, bound)?;
        let norm = certificate.norm_value(model, element)?;
        if norm.negative || norm.magnitude.small() != Some(degree) {
            return Err(Error::plain("CM cycle generator has the wrong norm"));
        }
        let trace = certificate
            .trace_value(model, element)?
            .signed()
            .ok_or_else(norm_limit)?;
        model.allocate(
            16 + (vertices.len() + steps.len() + j_path.len()) * (2 + start.field().cells()),
        )?;
        let cycle = Self {
            certificate_name: address.0,
            certificate,
            curve_name: source_name,
            vertices,
            steps,
            j_path,
            chart,
            scale,
            chart_count,
            differential,
            degree,
            conjugate: address.2,
            element,
            trace,
        };
        Ok(cycle)
    }

    pub(crate) fn apply(&self, model: &mut Model, mut point: Point) -> Result<Point> {
        for (i, map) in self.steps.iter().enumerate() {
            point = map.apply(model, &self.vertices[i], &self.vertices[i + 1], &point)?;
        }
        self.chart.apply(
            model,
            self.vertices.last().expect("cycle endpoint"),
            &self.vertices[0],
            &point,
        )
    }

    pub(crate) fn report(&self, name: &str) -> Result<String> {
        let orientation = usize::from(self.conjugate);
        let Element {
            constant,
            generator,
        } = self.element;
        let output = format!(
            "rvECCMCycleData[{name}] = curve={}, ell={}, orientation={orientation}, length={}, degree={}, trace={}, element=[{constant}, {generator}], chart_scale={:?}, chart_count={}, differential={:?}, j_path={:?}",
            self.curve_name,
            self.steps[0].kernel.ell,
            self.steps.len(),
            self.degree,
            self.trace,
            self.scale,
            self.chart_count,
            self.differential,
            self.j_path,
        );
        display(output)
    }

    pub(crate) fn step(&self, name: &str, index: usize) -> Result<String> {
        let map = self
            .steps
            .get(index)
            .ok_or_else(|| Error::plain("CM cycle step index is out of range"))?;
        let domain = &self.vertices[index];
        let range = &self.vertices[index + 1];
        display(format!(
            "rvECCMCycleStep[{name},{index}] = source_a={:?}, source_b={:?}, target_a={:?}, target_b={:?}, kernel={:?}",
            domain.a, domain.b, range.a, range.b, map.kernel.polynomial,
        ))
    }

    pub(crate) fn verify_expanded(&self, model: &mut Model, name: &str) -> Result<String> {
        let domain = &self.vertices[0];
        let degree = usize::try_from(self.degree).map_err(|_| norm_limit())?;
        domain
            .field()
            .check_polynomial(degree.checked_add(1).ok_or_else(norm_limit)?)?;
        let mut map = Morphism::from_velu(model, domain, &self.steps[0])?;
        for (i, step) in self.steps.iter().enumerate().skip(1) {
            let outer = Morphism::from_velu(model, &self.vertices[i], step)?;
            map = outer.compose(model, domain, &map)?;
        }
        map = self.chart.compose(model, domain, &map)?;
        if map.degree() != degree {
            return Err(Error::plain("expanded CM cycle has the wrong degree"));
        }
        let (q, trace) = domain.trace(model)?;
        let family = Family::from_frobenius(
            domain.ring.prime() as usize,
            i128::from(trace),
            i128::from(q),
        )?;
        let witness = Recognition::check(model, family, &map)?;
        if !witness.matches_order_element(self.certificate.frobenius, self.element)? {
            return Err(Error::plain(
                "expanded CM cycle has different endomorphism coordinates",
            ));
        }
        Ok(format!(
            "rvECCMCycleVerify[{name}] = degree={degree}, rational_map_agrees=true"
        ))
    }
}

fn display(output: String) -> Result<String> {
    if output.len() > MAX_SHOWN_BYTES {
        return Err(Error::plain(
            "CM isogeny cycle report exceeds the display limit",
        ));
    }
    Ok(output)
}

fn generator(
    model: &mut Model,
    action: &horizontal::Action,
    order: Order,
    differential: &[u64],
    four_degree: u128,
    bound: u128,
) -> Result<Element> {
    let field = action.curve.field();
    let certificate = &action.certificate;
    let pi = certificate.frobenius;
    let inverse = field.inverse_unit(model, &field.signed(pi.generator))?;
    let tangent_omega = field.multiply(model, &field.signed(-pi.constant), &inverse)?;
    let d = certificate.discriminant.unsigned_abs();
    let parity = certificate.generator_trace;
    let ell = i128::from(action.ell);
    let omega = action.ideal.generator_residue(order);
    let bound = i128::try_from(bound).map_err(|_| norm_limit())?;
    let mut answer = None;
    // An element a + b*omega of trace t and norm N satisfies
    // t^2 + |D| b^2 = 4N. The bound on |b| keeps this product in range.
    for b in -bound..=bound {
        let square = four_degree - d * b.unsigned_abs().pow(2);
        let root = square.isqrt();
        if root * root != square {
            continue;
        }
        let root = root as i128;
        for trace in [root, -root]
            .into_iter()
            .take(if root == 0 { 1 } else { 2 })
        {
            if (trace - parity * b) % 2 != 0 {
                continue;
            }
            let a = (trace - parity * b) / 2;
            let eval =
                |value: i128| (a.rem_euclid(ell) + b.rem_euclid(ell) * value).rem_euclid(ell);
            if eval(omega) != 0 || eval(parity - omega) == 0 {
                continue;
            }
            let linear = field.multiply(model, &field.signed(b), &tangent_omega)?;
            if field.add(model, &field.signed(a), &linear)? == differential
                && answer.replace(Element::new(a, b)).is_some()
            {
                return Err(Error::plain(
                    "CM cycle differential did not determine a unique generator",
                ));
            }
        }
    }
    answer.ok_or_else(|| {
        Error::plain("CM cycle has no principal generator matching the chart differential")
    })
}

fn chart_scale(model: &mut Model, domain: &Curve, range: &Curve) -> Result<(Scalar, usize)> {
    let field = domain.field();
    let for_a = chart_polynomial(model, field, &domain.a, &range.a, 4)?;
    let for_b = chart_polynomial(model, field, &domain.b, &range.b, 6)?;
    let condition = match (for_a, for_b) {
        (Some(a), Some(b)) => poly::monic_gcd(model, field, a, b)?,
        (Some(a), None) | (None, Some(a)) => a,
        (None, None) => return Err(no_chart()),
    };
    if condition.len() <= 1 {
        return Err(no_chart());
    }
    let q = field.ring.size_of_residue_field()?;
    let x = field.polynomial(&[0, 1]);
    let power = poly::power_mod(model, field, &x, q, &condition)?;
    let difference = poly::sum(model, field, &power, &x, true)?;
    let rational = poly::monic_gcd(model, field, condition, difference)?;
    if rational.len() <= 1 {
        return Err(no_chart());
    }
    let factors = finite_field_factor::factor_squarefree(model, field, rational)?;
    let mut roots = factors
        .into_iter()
        .map(|f| {
            if f.len() != 2 {
                return Err(no_chart());
            }
            field.ring.negate(model, &f[0])
        })
        .collect::<Result<Vec<_>>>()?;
    roots.sort_unstable();
    Ok((roots.first().ok_or_else(no_chart)?.clone(), roots.len()))
}

fn chart_polynomial(
    model: &mut Model,
    field: Field<'_>,
    domain: &Scalar,
    range: &Scalar,
    exponent: usize,
) -> Result<Option<Polynomial>> {
    match (poly::zero(domain), poly::zero(range)) {
        (true, true) => Ok(None),
        (true, false) | (false, true) => Err(no_chart()),
        (false, false) => {
            let inverse = field.inverse_unit(model, domain)?;
            let ratio = field.multiply(model, &inverse, range)?;
            let mut polynomial = vec![field.zero(); exponent + 1];
            polynomial[0] = field.ring.negate(model, &ratio)?;
            polynomial[exponent] = field.one();
            Ok(Some(polynomial))
        }
    }
}

fn no_chart() -> Error {
    Error::plain("principal CM cycle has no returning chart over the named residue field")
}

fn norm_limit() -> Error {
    Error::plain(
        "CM principal generator norm or enumeration exceeds the supported work or integer range",
    )
}
