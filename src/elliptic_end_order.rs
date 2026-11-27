use crate::cm_class_poly::Polynomial;
use crate::elliptic::{self, Curve};
use crate::elliptic_modular::decimal;
use crate::error::{Error, Result};
use crate::legendre_cm::Element;
use crate::model::{Model, Value};
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod canonical;

type Integer = Wide<13>;

pub(crate) enum Statement {
    Define(Label, Label, Label),
    Check(Label, Label),
    Data(Label),
    Invariant(bool, Box<ParsedExpr>),
}

pub(crate) enum Expression {
    Literal(Label, i128, i128),
    Frobenius(Label),
    JInvariant(Label),
    Multiply(Box<[ParsedExpr; 2]>),
    Conjugate(Box<ParsedExpr>),
}

pub(crate) struct Evidence {
    discriminant: i128,
    p: u64,
    q: u64,
    trace: i128,
    delta: i128,
    j: Vec<u64>,
    residual: Vec<u64>,
    ordinary: bool,
    split: bool,
    index: Option<i128>,
}

impl Evidence {
    pub(crate) fn read(model: &mut Model, curve: &Curve, polynomial: &Polynomial) -> Result<Self> {
        let p = curve.ring.prime();
        let (q, trace) = curve.trace(model)?;
        let trace = i128::from(trace);
        let delta = trace * trace - 4 * i128::from(q);
        let discriminant = polynomial.discriminant;
        let j = curve.j_invariant(model)?;
        let residual = polynomial.evaluate(model, &curve.ring, &j)?;
        let ordinary = trace.rem_euclid(i128::from(p)) != 0;
        let residue = discriminant.rem_euclid(i128::from(p)) as u64;
        // A nonzero square is the splitting criterion and also rules out
        // p dividing the discriminant, including p dividing the conductor.
        let split =
            residue != 0 && elliptic::integer_power_mod(model, residue, (p - 1) / 2, p)? == 1;
        let index = if delta < 0 && delta % discriminant == 0 {
            let ratio = (delta / discriminant) as u128;
            let root = ratio.isqrt();
            (root > 0 && root * root == ratio).then_some(root as i128)
        } else {
            None
        };
        Ok(Self {
            discriminant,
            p,
            q,
            trace,
            delta,
            j,
            residual,
            ordinary,
            split,
            index,
        })
    }

    fn is_root(&self) -> bool {
        self.residual.iter().all(|&c| c == 0)
    }

    pub(crate) fn certified(&self) -> bool {
        self.ordinary && self.split && self.index.is_some() && self.is_root()
    }

    pub(crate) fn report(&self, curve: &str, polynomial: &str) -> String {
        let index = self.index.map_or_else(|| "none".into(), |n| n.to_string());
        format!(
            "rvECEndCheck[{curve},{polynomial}] = p={}, q={}, trace={}, frobenius_discriminant={}, candidate_discriminant={}, j={:?}, residual={:?}, ordinary={}, split={}, index={index}, certified={}",
            self.p, self.q, self.trace, self.delta, self.discriminant, self.j, self.residual, self.ordinary,
            self.split, self.certified(),
        )
    }

    pub(crate) fn certify(self, model: &mut Model, curve: String) -> Result<Certificate> {
        if !self.ordinary {
            return Err(Error::plain(
                "elliptic endomorphism-order certification requires an ordinary curve",
            ));
        }
        if !self.is_root() {
            return Err(Error::plain(
                "elliptic j-invariant is not a root of the candidate CM class polynomial",
            ));
        }
        let index = self.index.ok_or_else(|| {
            Error::plain("Frobenius discriminant divided by the candidate CM discriminant must be a positive square")
        })?;
        if !self.split {
            return Err(Error::plain(
                "elliptic characteristic must split and be prime to the candidate CM discriminant",
            ));
        }
        let parity = self.discriminant.rem_euclid(2);
        let generator_norm = -((self.discriminant - parity) / 4);
        let numerator = self.trace - parity * index;
        if numerator % 2 != 0 {
            return Err(Error::plain(
                "certified Frobenius coordinates must be integral",
            ));
        }
        let frobenius = Element::new(numerator / 2, index);
        let certificate = Certificate {
            curve,
            discriminant: self.discriminant,
            q: self.q,
            trace: self.trace,
            frobenius_discriminant: self.delta,
            j: self.j,
            generator_trace: parity,
            generator_norm,
            frobenius,
            index,
        };
        let a = certificate.trace;
        let q = i128::from(certificate.q);
        let pi2 = certificate.multiply(model, frobenius, frobenius)?;
        let minus_a_pi = certificate.multiply(model, frobenius, Element::new(-a, 0))?;
        let relation = pi2.add(model, minus_a_pi)?.add(model, Element::new(q, 0))?;
        if relation != Element::ZERO
            || certificate.trace_value(model, frobenius)?.signed() != Some(a)
            || certificate.norm_value(model, frobenius)?.signed() != Some(q)
        {
            return Err(Error::plain(
                "certified Frobenius failed the quadratic relation",
            ));
        }
        model.allocate(16 + certificate.curve.len() + certificate.j.len())?;
        Ok(certificate)
    }
}

pub(crate) struct Certificate {
    curve: String,
    discriminant: i128,
    q: u64,
    trace: i128,
    frobenius_discriminant: i128,
    j: Vec<u64>,
    generator_trace: i128,
    generator_norm: i128,
    pub(crate) frobenius: Element,
    index: i128,
}

impl Certificate {
    pub(crate) fn discriminant(&self) -> i128 {
        self.discriminant
    }

    pub(crate) fn report(&self, name: &str) -> String {
        let Element {
            constant,
            generator,
        } = self.frobenius;
        format!(
            "rvECEndData[{name}] = curve={}, q={}, trace={}, frobenius_discriminant={}, discriminant={}, j={:?}, index={}, frobenius=[{constant}, {generator}]",
            self.curve, self.q, self.trace, self.frobenius_discriminant, self.discriminant, self.j, self.index,
        )
    }

    pub(crate) fn multiply(&self, model: &mut Model, x: Element, y: Element) -> Result<Element> {
        x.product(model, y, self.generator_trace, self.generator_norm)
    }

    pub(crate) fn conjugate(&self, model: &mut Model, x: Element) -> Result<Element> {
        model.spend(4)?;
        let constant = Integer::from_i128(x.constant)
            .add(Integer::product(self.generator_trace, x.generator))?
            .signed()
            .ok_or_else(arithmetic)?;
        Ok(Element::new(
            constant,
            x.generator.checked_neg().ok_or_else(arithmetic)?,
        ))
    }

    fn trace_value(&self, model: &mut Model, x: Element) -> Result<Integer> {
        model.spend(40)?;
        Integer::product(2, x.constant).add(Integer::product(self.generator_trace, x.generator))
    }

    fn norm_value(&self, model: &mut Model, x: Element) -> Result<Integer> {
        model.spend(80)?;
        Integer::product(x.constant, x.constant)
            .add(Integer::product(x.constant, x.generator).multiply(self.generator_trace)?)?
            .add(Integer::product(x.generator, x.generator).multiply(self.generator_norm)?)
    }

    pub(crate) fn invariant(&self, model: &mut Model, norm: bool, x: Element) -> Result<String> {
        let value = if norm {
            self.norm_value(model, x)?
        } else {
            self.trace_value(model, x)?
        };
        Ok(decimal(value))
    }
}

pub(crate) fn persist(model: &mut Model, value: Element) -> Result<Value> {
    model.allocate(4)?;
    Ok(Value::ECEndOrder(value))
}

fn arithmetic() -> Error {
    Error::plain("CM order arithmetic exceeds signed i128 range")
}
