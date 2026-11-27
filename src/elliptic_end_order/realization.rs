use std::rc::Rc;

use super::{cycle, Certificate, Integer};
use crate::ec_morphism::Morphism;
use crate::elliptic::{Curve, Point};
use crate::elliptic_cm_recognition::{frobenius_over, sum};
use crate::error::{Error, Result};
use crate::legendre_cm::Element;
use crate::model::{Model, MAX_SHOWN_BYTES};
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Statement {
    Define(Label, Label, Vec<Label>),
    BaseChange(Label, Label, Label),
    Data(Label),
}

pub(crate) enum Expression {
    Apply(Label, Box<[ParsedExpr; 2]>),
    Expand(Label, Box<ParsedExpr>),
    Generator(Label),
}

struct CycleTerm {
    name: String,
    cycle: Rc<cycle::Cycle>,
    coefficient: i128,
}

struct Witness {
    certificate_name: Label,
    certificate: Rc<Certificate>,
    base_curve: Rc<Curve>,
    index: i128,
    pi_coefficient: i128,
    cycles: Vec<CycleTerm>,
}

pub(crate) struct Realization {
    witness: Rc<Witness>,
    pub(crate) curve_name: Rc<str>,
    curve: Rc<Curve>,
    maps: Vec<Rc<cycle::PointMap>>,
}

impl Realization {
    pub(crate) fn new(model: &mut Model, name: Label, cycles: Vec<Label>) -> Result<Self> {
        let certificate = model.ec_end_order(&name.text)?;
        let curve_name: Rc<str> = certificate.curve_name().into();
        let base_curve = model.ec_curve(&curve_name)?;
        let mut index = certificate.frobenius.generator;
        let mut pi_coefficient = 1i128;
        let mut terms: Vec<CycleTerm> = Vec::new();
        for label in cycles {
            model.spend(1)?;
            let cycle = model.ec_cm_cycle(&label)?;
            if cycle.certificate_name.text != name.text {
                return Err(Error::new(
                    label.span,
                    "CM action cycle must use the same named endomorphism-order certificate",
                ));
            }
            let b = cycle.element.generator;
            if b % index == 0 {
                continue;
            }
            let (gcd, old, new) = bezout(model, index, b)?;
            pi_coefficient = pi_coefficient.checked_mul(old).ok_or_else(range)?;
            for term in &mut terms {
                term.coefficient = term.coefficient.checked_mul(old).ok_or_else(range)?;
            }
            terms.retain(|term| term.coefficient != 0);
            terms.push(CycleTerm {
                name: label.text,
                cycle,
                coefficient: new,
            });
            index = gcd;
        }
        let mut check = Integer::product(pi_coefficient, certificate.frobenius.generator);
        for term in &terms {
            check = check.add(Integer::product(
                term.coefficient,
                term.cycle.element.generator,
            ))?;
        }
        if check.signed() != Some(index) {
            return Err(Error::plain(
                "CM action failed its integral Bézout identity",
            ));
        }

        model.allocate(
            16 + 8 * terms.len() + terms.iter().map(|term| term.name.len()).sum::<usize>(),
        )?;
        let maps = terms
            .iter()
            .map(|term| Rc::clone(&term.cycle.map))
            .collect();
        let curve = Rc::clone(&base_curve);
        let witness = Rc::new(Witness {
            certificate_name: name,
            certificate,
            base_curve,
            index,
            pi_coefficient,
            cycles: terms,
        });
        Ok(Self {
            witness,
            curve_name,
            curve,
            maps,
        })
    }

    pub(crate) fn base_change(&self, model: &mut Model, target: &Label) -> Result<Self> {
        let curve_name: Rc<str> = target.text.as_str().into();
        let curve = model.ec_curve(&curve_name)?;
        let base = &self.witness.base_curve;
        let embedding = base.ring.residue_embedding(model, &curve.ring)?;
        let coefficients = embedding.coefficients(model, &[base.a.clone(), base.b.clone()])?;
        if coefficients != [curve.a.clone(), curve.b.clone()] {
            return Err(Error::new(
                target.span,
                "target elliptic equation is not the declared CM base extension",
            ));
        }
        model.allocate(8 + self.witness.cycles.len())?;
        let maps = self
            .witness
            .cycles
            .iter()
            .map(|term| {
                term.cycle
                    .base_change(model, &embedding, Rc::clone(&curve_name), Rc::clone(&curve))
            })
            .collect::<Result<_>>()?;
        Ok(Self {
            witness: Rc::clone(&self.witness),
            curve_name,
            curve,
            maps,
        })
    }

    pub(crate) fn certificate_name(&self) -> &Label {
        &self.witness.certificate_name
    }

    pub(crate) fn generator(&self) -> Element {
        Element::new(0, self.witness.index)
    }

    pub(crate) fn report(&self, name: &Label) -> Result<String> {
        let witness = &self.witness;
        let conductor_ratio = witness.certificate.frobenius.generator / witness.index;
        let discriminant =
            witness.certificate.frobenius_discriminant / (conductor_ratio * conductor_ratio);
        let cycles = witness
            .cycles
            .iter()
            .map(|term| (term.name.as_str(), term.coefficient))
            .collect::<Vec<_>>();
        let output = format!(
            "rvECCMActionData[{}] = certificate={}, base={}, curve={}, index={}, discriminant={discriminant}, full={}, generator=[0, {}], frobenius_weight={}, cycle_weights={cycles:?}",
            name.text, witness.certificate_name.text, witness.certificate.curve_name(),
            self.curve_name, witness.index, witness.index == 1, witness.index,
            witness.pi_coefficient,
        );
        if output.len() > MAX_SHOWN_BYTES {
            return Err(Error::new(
                name.span,
                "CM action report exceeds the display limit",
            ));
        }
        Ok(output)
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        element: Element,
        point: &Point,
    ) -> Result<Point> {
        let witness = &self.witness;
        let multiplier = self.multiplier(element)?;
        let linear = self.curve.multiple(model, point, element.constant)?;
        if multiplier == 0 || *point == Point::Infinity {
            return Ok(linear);
        }
        let mut generator = Point::Infinity;
        if witness.pi_coefficient != 0 {
            let p_power = witness.base_curve.ring.degree();
            let ring = &self.curve.ring;
            let Point::Affine(x, y) = point else {
                unreachable!()
            };
            let pi_point = Point::Affine(
                ring.frobenius(model, x, p_power)?,
                ring.frobenius(model, y, p_power)?,
            );
            generator = self.shifted(
                model,
                point,
                &pi_point,
                witness.certificate.frobenius.constant,
                witness.pi_coefficient,
            )?;
        }
        for (term, map) in witness.cycles.iter().zip(&self.maps) {
            let image = map.apply(model, point.clone())?;
            let summand = self.shifted(
                model,
                point,
                &image,
                term.cycle.element.constant,
                term.coefficient,
            )?;
            generator = self.curve.add(model, &generator, &summand)?;
        }
        let nonlinear = self.curve.multiple(model, &generator, multiplier)?;
        self.curve.add(model, &linear, &nonlinear)
    }

    pub(crate) fn expand(&self, model: &mut Model, element: Element) -> Result<Morphism> {
        let witness = &self.witness;
        let multiplier = self.multiplier(element)?;
        let degree = witness
            .certificate
            .norm_value(model, element)?
            .magnitude
            .small()
            .and_then(|n| usize::try_from(n).ok())
            .ok_or_else(expansion_limit)?;
        self.curve
            .field()
            .check_polynomial(degree.checked_add(1).ok_or_else(expansion_limit)?)?;
        let linear = self.scalar_map(model, element.constant)?;
        if multiplier == 0 {
            return Ok(linear);
        }

        let mut generator =
            Morphism::zero(Rc::clone(&self.curve_name), Rc::clone(&self.curve_name));
        if witness.pi_coefficient != 0 {
            let pi = frobenius_over(
                model,
                Rc::clone(&self.curve_name),
                &self.curve,
                witness.certificate.q,
            )?;
            generator = self.shifted_map(
                model,
                &pi,
                witness.certificate.frobenius.constant,
                witness.pi_coefficient,
            )?;
        }
        for (term, view) in witness.cycles.iter().zip(&self.maps) {
            let map = term.cycle.expanded(model, view)?;
            let summand =
                self.shifted_map(model, &map, term.cycle.element.constant, term.coefficient)?;
            generator = sum(model, &self.curve, &self.curve, &generator, &summand)?;
        }
        let nonlinear =
            self.scalar_map(model, multiplier)?
                .compose(model, &self.curve, &generator)?;
        let map = sum(model, &self.curve, &self.curve, &linear, &nonlinear)?;
        if map.degree() != degree {
            return Err(Error::plain(
                "expanded CM action has degree different from the CM norm",
            ));
        }
        Ok(map)
    }

    fn scalar_map(&self, model: &mut Model, scalar: i128) -> Result<Morphism> {
        Morphism::multiple(model, Rc::clone(&self.curve_name), &self.curve, scalar)
    }

    fn shifted_map(
        &self,
        model: &mut Model,
        map: &Morphism,
        constant: i128,
        coefficient: i128,
    ) -> Result<Morphism> {
        let shift = self.scalar_map(model, constant.checked_neg().ok_or_else(range)?)?;
        let difference = sum(model, &self.curve, &self.curve, map, &shift)?;
        self.scalar_map(model, coefficient)?
            .compose(model, &self.curve, &difference)
    }

    fn multiplier(&self, element: Element) -> Result<i128> {
        let index = self.witness.index;
        if element.generator % index != 0 {
            return Err(Error::plain(
                "CM element lies outside the realized suborder",
            ));
        }
        Ok(element.generator / index)
    }

    fn shifted(
        &self,
        model: &mut Model,
        point: &Point,
        image: &Point,
        constant: i128,
        coefficient: i128,
    ) -> Result<Point> {
        let shift = self
            .curve
            .multiple(model, point, constant.checked_neg().ok_or_else(range)?)?;
        let difference = self.curve.add(model, image, &shift)?;
        self.curve.multiple(model, &difference, coefficient)
    }
}

fn bezout(model: &mut Model, a: i128, b: i128) -> Result<(i128, i128, i128)> {
    let mut next = b.checked_abs().ok_or_else(range)?;
    let mut order = a;
    let mut coefficients = [1i128, 0];
    let mut next_coefficients = [0i128, if b < 0 { -1 } else { 1 }];
    while next != 0 {
        model.spend(1)?;
        let quotient = order / next;
        (order, next) = (next, order % next);
        for (left, right) in coefficients.iter_mut().zip(&mut next_coefficients) {
            let product = quotient.checked_mul(*right).ok_or_else(range)?;
            (*left, *right) = (*right, left.checked_sub(product).ok_or_else(range)?);
        }
    }
    let [left, right] = coefficients;
    Ok((order, left, right))
}

fn range() -> Error {
    Error::plain("CM action Bézout coefficients exceed the supported integer range")
}

fn expansion_limit() -> Error {
    Error::plain("expanded CM action degree exceeds the supported integer range")
}
