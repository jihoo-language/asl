use std::rc::Rc;

use super::cm_cycles::{pair_gram, Cycle};
use crate::ec_morphism::Morphism;
use crate::elliptic;
use crate::elliptic_cm_recognition as maps;
use crate::elliptic_end_order::{conductor, Certificate};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::cm_ring::{Class, Order};
use crate::model::Model;
use crate::rv_ore::gcd;
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod composition;
pub(crate) mod marked;

pub(crate) enum Statement {
    Define(Label, Option<[Label; 2]>, Box<[ParsedExpr; 2]>, bool),
    Composition(Box<composition::Statement>),
    Mark(Label, Label, Label, Label, [Label; 2]),
    Data(Label),
    MarkedData(Label),
    Norm(Label, [i128; 2]),
    CyclePair(Box<[ParsedExpr; 2]>),
    WedgePair(Box<[ParsedExpr; 2]>),
}

pub(crate) enum Expression {
    Class(Label, Label),
    Composition(Box<composition::Expression>),
    Map(Label, [i128; 2]),
    Cycle(Label, Cycle),
    Wedge(Label, [i128; 6]),
    Lift(Label, Box<ParsedExpr>),
    Project(bool, Box<ParsedExpr>),
    Frobenius(usize, Box<ParsedExpr>),
    Change(Label, Box<ParsedExpr>),
}

pub(crate) struct Hom {
    pub(crate) maps: [Rc<Morphism>; 2],
    pub(crate) endpoints: [Rc<elliptic::Curve>; 2],
    pub(crate) q: u64,
    pub(crate) trace: i128,
    j: [crate::ec_poly::Scalar; 2],
    form: [i128; 3],
    reduced: Class,
    fundamental: i128,
    conductor: u64,
    span_discriminant: i128,
    primitive_discriminant: i128,
    content: i128,
    maximal: bool,
    index: u128,
    plus: usize,
    minus: usize,
    pub(crate) gram: [[i128; 4]; 4],
    generators: [[i128; 2]; 2],
}

impl Hom {
    pub(crate) fn new(
        model: &mut Model,
        pair: [Rc<Morphism>; 2],
        witnesses: Option<[Rc<Certificate>; 2]>,
        nonprincipal_only: bool,
    ) -> Result<Self> {
        let [phi, psi] = &pair;
        if phi.source != psi.source || phi.target != psi.target {
            return Err(Error::plain(
                "CM Hom maps require the same named source and target",
            ));
        }
        if phi.degree() == 0 || psi.degree() == 0 {
            return Err(Error::plain("CM Hom maps must have positive degree"));
        }
        let endpoints = [model.ec_curve(&phi.source)?, model.ec_curve(&phi.target)?];
        let [source, target] = endpoints.each_ref().map(Rc::as_ref);
        phi.certify(model, source, target)?;
        psi.certify(model, source, target)?;
        let j = [source.j_invariant(model)?, target.j_invariant(model)?];
        if nonprincipal_only && j[0] == j[1] {
            return Err(Error::plain(
                "CM Hom correspondences require distinct elliptic j-invariants",
            ));
        }
        let (q, trace) = source.trace(model)?;
        if trace.rem_euclid(source.ring.prime() as i64) == 0 || (q, trace) != target.trace(model)? {
            return Err(Error::plain(
                "CM Hom correspondences require an ordinary isogenous pair over the common field",
            ));
        }
        let trace = i128::from(trace);
        let frobenius_discriminant = trace * trace - 4 * i128::from(q);
        let conductor = conductor::frobenius_conductor(model, frobenius_discriminant)?;
        let fundamental = frobenius_discriminant / i128::from(conductor).pow(2);
        let maximal = if let Some([e, f]) = witnesses {
            if e.curve_name() != phi.source.as_ref()
                || f.curve_name() != phi.target.as_ref()
                || e.discriminant() != fundamental
                || f.discriminant() != fundamental
            {
                return Err(Error::plain("CM Hom endpoint certificates must prove the common maximal order on the actual map endpoints"));
            }
            true
        } else {
            conductor == 1
        };
        let negative = maps::negate(model, source, target, psi)?;
        let plus = maps::sum(model, source, target, phi, psi)?.degree();
        let minus = maps::sum(model, source, target, phi, &negative)?.degree();
        let a = phi.degree() as i128;
        let c = psi.degree() as i128;
        let b = plus as i128 - a - c;
        if a + c - minus as i128 != b {
            return Err(inconsistent());
        }
        let span_discriminant = b
            .checked_mul(b)
            .and_then(|square| {
                a.checked_mul(c)?
                    .checked_mul(4)
                    .and_then(|x| square.checked_sub(x))
            })
            .ok_or_else(overflow)?;
        if span_discriminant >= 0 {
            return Err(Error::plain(
                "CM Hom witnesses must be independent over the integers",
            ));
        }
        let discriminants = (span_discriminant.unsigned_abs(), fundamental.unsigned_abs());
        let ratio = discriminants.0 / discriminants.1;
        let index = ratio.isqrt();
        if index == 0 || !discriminants.0.is_multiple_of(discriminants.1) || ratio != index * index
        {
            return Err(inconsistent());
        }
        let k = i128::try_from(index).map_err(|_| overflow())?;
        let parity = fundamental.rem_euclid(2);
        let second = b - k * parity;
        if second % 2 != 0 {
            return Err(inconsistent());
        }
        // Composition on the left by phi-dagger sends (phi,psi) to
        // (deg phi,(b+sqrt(discriminant))/2), up to conjugate orientation.
        let generators = [[a, 0], [second / 2, k]];
        let form = [a, b, c];
        let content = gcd(gcd(a.unsigned_abs(), b.unsigned_abs()), c.unsigned_abs()) as i128;
        let primitive_discriminant = span_discriminant / content / content;
        let reduced =
            Order::new(primitive_discriminant)?.class(model, form.map(|x| x / content))?;
        if maximal && index == 1 && j[0] != j[1] && reduced.0[0] * content <= 1 {
            return Err(inconsistent());
        }
        let gram = [
            [0, 1, a, c],
            [1, 0, 1, 1],
            [a, 1, 0, minus as i128],
            [c, 1, minus as i128, 0],
        ];
        model.allocate(96)?;
        Ok(Self {
            maps: pair,
            endpoints,
            q,
            trace,
            j,
            form,
            reduced,
            fundamental,
            conductor,
            span_discriminant,
            primitive_discriminant,
            content,
            maximal,
            index,
            plus,
            minus,
            gram,
            generators,
        })
    }

    pub(crate) fn degree(&self, model: &mut Model, [a, b]: [i128; 2]) -> Result<i128> {
        model.spend(3)?;
        let [l, m, r] = self.form;
        Wide::<16>::product(a, a)
            .multiply(l)?
            .add(Wide::product(a, b).multiply(m)?)?
            .add(Wide::product(b, b).multiply(r)?)?
            .signed()
            .ok_or_else(overflow)
    }

    pub(crate) fn pair(&self, model: &mut Model, a: Cycle, b: Cycle) -> Result<i128> {
        pair_gram(model, &self.gram, a, b)
    }

    pub(crate) fn class(&self, model: &mut Model, order: Order) -> Result<Class> {
        if !self.maximal || self.index != 1 || order.discriminant() != self.fundamental {
            return Err(Error::plain("CM Hom ideal class requires a full certified maximal-order span and matching order"));
        }
        let [a, b, c] = self.form;
        // The class-group convention for a form with middle coefficient b
        // represents the ideal <a,(-b+sqrt(D))/2>.
        order.class(model, [a, -b, c])
    }

    pub(crate) fn combination(&self, model: &mut Model, coeffs: [i128; 2]) -> Result<Morphism> {
        let endpoints = self.endpoints.each_ref().map(Rc::as_ref);
        let a = composition::scale_map(model, endpoints, &self.maps[0], coeffs[0])?;
        let b = composition::scale_map(model, endpoints, &self.maps[1], coeffs[1])?;
        let [domain, range] = endpoints;
        let map = maps::sum(model, domain, range, &a, &b)?;
        if map.degree() as i128 != self.degree(model, coeffs)? {
            return Err(inconsistent());
        }
        Ok(map)
    }

    pub(crate) fn report(&self, name: &str) -> String {
        let min = self.reduced.0[0] * self.content;
        let full_minimum = if self.maximal && self.index == 1 {
            min.to_string()
        } else {
            "unproved".into()
        };
        let status = if self.maximal {
            let nonprincipal = if self.index == 1 {
                if min == 1 {
                    "false"
                } else {
                    "true"
                }
            } else if self.j[0] != self.j[1] {
                "true"
            } else if self.maps[0].source == self.maps[0].target {
                "false"
            } else {
                "unproved"
            };
            format!(
                "span_index={}, full_lattice={}, nonprincipal={nonprincipal}",
                self.index,
                self.index == 1
            )
        } else {
            "span_index=unproved, full_lattice=unproved, nonprincipal=unproved".into()
        };
        format!("rvCMHomData[{name}] = j={:?}, trace={}, frobenius_conductor={}, fundamental_discriminant={}, maximal_order_certified={}, degrees=[{}, {}], degree_sum={}, degree_difference={}, norm_form={:?}, primitive_reduced={:?}, content={}, primitive_discriminant={}, span_minimum={min}, full_minimum={full_minimum}, graph_gram={:?}, graph_discriminant={}, image_basis_1_omega={:?}, {status}",
            self.j, self.trace, self.conductor, self.fundamental, self.maximal,
            self.form[0], self.form[2], self.plus, self.minus, self.form, self.reduced.0, self.content,
            self.primitive_discriminant, self.gram, self.span_discriminant, self.generators)
    }
}

pub(super) fn inconsistent() -> Error {
    Error::plain("CM Hom norm, gluing, or marked intersection certificates disagree")
}

fn overflow() -> Error {
    Error::plain("CM Hom exact degree arithmetic exceeds signed 128-bit range")
}
