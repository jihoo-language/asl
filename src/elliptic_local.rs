use std::rc::Rc;

use crate::ec_poly::{self as poly, Polynomial};
use crate::elliptic::{self, Curve, DivisionPolynomials, Isogeny, Kernel};
use crate::elliptic_end_order::conductor;
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::legendre_cm;
use crate::model::Model;
use crate::rv_qp_quad::{self, Element as QuadElement};
use crate::syntax::Label;

pub(crate) enum Command {
    Certificate(Label, Test),
    Order(Label, usize),
    Field(Label, Label),
    Edge(Label),
    Neighborhood(Label, usize),
}

pub(crate) enum Expression {
    Coordinates(Label, Test),
    Generator(Label, Label, usize),
}

#[derive(Clone, Copy)]
pub(crate) struct Test {
    pub(crate) prime: usize,
    pub(crate) level: usize,
    pub(crate) candidate: Option<i128>,
}

pub(crate) struct Certificate {
    modulus: Polynomial,
    frobenius_x: Polynomial,
    frobenius_y: Polynomial,
    residual_x: Polynomial,
    residual_y: Polynomial,
    order: u64,
    candidate: u64,
    representative: i128,
}

impl Certificate {
    fn build(
        model: &mut Model,
        curve: &Curve,
        division: &mut DivisionPolynomials<'_>,
        q: u64,
        ell: u64,
        order: u64,
        candidate: u64,
    ) -> Result<Self> {
        let field = curve.field();
        if candidate.is_multiple_of(ell) {
            return Err(Error::plain(
                "scalar candidate must be prime to the torsion order",
            ));
        }
        let order_i = i128::from(order);
        let mut representative = i128::from(candidate);
        if candidate > order / 2 {
            representative -= order_i;
        }
        if representative % i128::from(field.ring.prime()) == 0 {
            representative += if representative < 0 {
                order_i
            } else {
                -order_i
            };
        }
        let square = u128::from(order).pow(2);
        let degree = if ell == 2 {
            (square + 2) / 2
        } else {
            (square - 1) / 2
        };
        let size = usize::try_from(degree + 1).map_err(|_| {
            Error::plain("geometric torsion division polynomial degree is out of range")
        })?;
        field.check_polynomial(size)?;
        let psi = division.get(model, order)?;
        let branch_degree = if ell == 2 { 3 } else { 0 };
        if psi.len() + branch_degree != size {
            return Err(Error::plain(
                "geometric torsion division polynomial has unexpected degree",
            ));
        }
        let unit = field.inverse_unit(model, psi.last().unwrap())?;
        let nonbranch = poly::scale(model, field, &psi, &unit)?;
        // For even order psi_m/(2y) omits the three nonzero points
        // of order two. They belong to the x test but not to y^q/y.
        let full = if ell == 2 {
            Some(poly::multiply(
                model,
                field,
                curve.polynomial(),
                &nonbranch,
            )?)
        } else {
            None
        };
        let modulus = full.as_deref().unwrap_or(&nonbranch);
        let n = u64::try_from(representative.unsigned_abs())
            .expect("representative fits torsion order");
        let (numerator, denominator) = division.x_image(model, n)?;
        // Differentiate the affine rational x-map before reducing any of its
        // polynomials. The invariant differential gives the signed y-map.
        let numerator_prime = poly::derivative(model, field, &numerator)?;
        let denominator_prime = poly::derivative(model, field, &denominator)?;
        let numerator = poly::remainder(model, field, numerator, modulus)?;
        let denominator = poly::remainder(model, field, denominator, modulus)?;

        let variable = field.polynomial(&[0, 1]);
        let frobenius_x = poly::power_mod(model, field, &variable, q, modulus)?;
        let product = poly::multiply_mod(model, field, &frobenius_x, &denominator, modulus)?;
        let residual_x = poly::sum(model, field, &product, &numerator, true)?;

        let (frobenius_y, residual_y) = if nonbranch.len() == 1 {
            (Vec::new(), Vec::new())
        } else {
            let modulus = &nonbranch;
            let numerator_prime = poly::remainder(model, field, numerator_prime, modulus)?;
            let denominator_prime = poly::remainder(model, field, denominator_prime, modulus)?;
            let frobenius_y =
                poly::power_mod(model, field, curve.polynomial(), (q - 1) / 2, modulus)?;
            let square = poly::multiply_mod(model, field, &denominator, &denominator, modulus)?;
            let scaled = poly::scale(model, field, &square, &field.signed(representative))?;
            let target = poly::multiply_mod(model, field, &frobenius_y, &scaled, modulus)?;
            let left = poly::multiply_mod(model, field, &numerator_prime, &denominator, modulus)?;
            let right = poly::multiply_mod(model, field, &numerator, &denominator_prime, modulus)?;
            let derivative = poly::sum(model, field, &left, &right, true)?;
            let residual = poly::sum(model, field, &target, &derivative, true)?;
            (frobenius_y, residual)
        };
        Ok(Self {
            modulus: full.unwrap_or(nonbranch),
            frobenius_x,
            frobenius_y,
            residual_x,
            residual_y,
            order,
            candidate,
            representative,
        })
    }

    pub(crate) fn checks(&self) -> (bool, bool) {
        (self.residual_x.is_empty(), self.residual_y.is_empty())
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, test: Test) -> Result<String> {
        let (x, y) = self.checks();
        let mut output = format!(
            "rvEllFrobCert[{name}, {}, {}] = m={}, a={}, representative={}, deg(psi)={}, x={x}, y={y}",
            test.prime, test.level,self.order,self.candidate,self.representative,self.modulus.len()-1,
        );
        for (label, terms) in [
            ("psi", &self.modulus),
            ("frob_x", &self.frobenius_x),
            ("frob_y", &self.frobenius_y),
            ("residual_x", &self.residual_x),
            ("residual_y", &self.residual_y),
        ] {
            append_polynomial(model, &mut output, label, terms)?;
        }
        if output.len() > crate::model::MAX_SHOWN_BYTES {
            return Err(Error::plain("value exceeds the display size limit"));
        }
        Ok(output)
    }
}

fn append_polynomial(
    model: &mut Model,
    output: &mut String,
    label: &str,
    polynomial: &Polynomial,
) -> Result<()> {
    output.push_str(&format!(", {label}=["));
    let mut first = true;
    for (i, coefficient) in polynomial.iter().enumerate() {
        model.spend(1 + coefficient.len())?;
        if poly::zero(coefficient) {
            continue;
        }
        if !first {
            output.push_str(", ");
        }
        first = false;
        let term = format!("({i}, {coefficient:?})");
        if output.len() + term.len() > crate::model::MAX_SHOWN_BYTES {
            return Err(Error::plain("value exceeds the display size limit"));
        }
        output.push_str(&term);
    }
    output.push(']');
    Ok(())
}

fn scalar_from_trace(trace: i128, order: u64) -> u64 {
    let order = i128::from(order);
    let trace = trace.rem_euclid(order);
    let half = if trace & 1 == 0 {
        trace / 2
    } else {
        (trace + order) / 2
    };
    half as u64
}

fn scalar_candidates(trace: i128, ell: u64, order: u64) -> Vec<u64> {
    if ell != 2 {
        return vec![scalar_from_trace(trace, order)];
    }
    if trace % 2 != 0 {
        return Vec::new();
    }
    let half = order / 2;
    let reduced = (trace / 2).rem_euclid(i128::from(half)) as u64;
    [reduced, reduced + half]
        .into_iter()
        .filter(|value| value % 2 != 0)
        .collect()
}

fn divided_characteristic(q: u64, trace: i128, scalar: u64, modulus: u64) -> Option<(i128, i128)> {
    let m = i128::from(modulus);
    let a = i128::from(scalar);
    let trace_quotient = trace - 2 * a;
    let norm_quotient = i128::from(q)
        .checked_sub(trace.checked_mul(a)?)?
        .checked_add(a.checked_mul(a)?)?;
    let square = m.checked_mul(m)?;
    (trace_quotient % m == 0 && norm_quotient % square == 0)
        .then(|| (trace_quotient / m, norm_quotient / square))
}

pub(crate) fn certificate(model: &mut Model, curve: &Curve, test: Test) -> Result<Certificate> {
    let ell = elliptic::validate_ell(model, curve, test.prime)?;
    let order = u32::try_from(test.level)
        .ok()
        .and_then(|e| ell.checked_pow(e))
        .ok_or_else(|| Error::plain("geometric torsion order ell^k exceeds u64"))?;
    if order == 1 {
        return Ok(Certificate {
            modulus: vec![curve.field().one()],
            frobenius_x: Vec::new(),
            frobenius_y: Vec::new(),
            residual_x: Vec::new(),
            residual_y: Vec::new(),
            order,
            candidate: 0,
            representative: 1,
        });
    }
    let q = curve.ring.size_of_residue_field()?;
    let candidate = if let Some(value) = test.candidate {
        value.rem_euclid(i128::from(order)) as u64
    } else {
        let trace = i128::from(curve.trace(model)?.1);
        let candidates = scalar_candidates(trace, ell, order);
        candidates
            .iter()
            .copied()
            .find(|&a| divided_characteristic(q, trace, a, order).is_some())
            .or_else(|| candidates.first().copied())
            .unwrap_or(1)
    };
    let mut division = DivisionPolynomials::new(model, curve, None)?;
    Certificate::build(model, curve, &mut division, q, ell, order, candidate)
}

#[derive(Clone, Copy, Eq, PartialEq)]
struct Frobenius {
    prime: usize,
    q: u64,
    trace: i128,
    discriminant: i128,
}

impl Frobenius {
    fn ordinary(model: &mut Model, curve: &Curve) -> Result<Self> {
        let prime = curve.ring.prime() as usize;
        let (q, trace) = curve.trace(model)?;
        let trace = i128::from(trace);
        if trace % prime as i128 == 0 {
            return Err(Error::plain(
                "local geometric endomorphism order requires an ordinary curve",
            ));
        }
        let discriminant = trace * trace - 4 * i128::from(q);
        if discriminant >= 0 {
            return Err(Error::plain(
                "ordinary elliptic Frobenius discriminant must be negative",
            ));
        }
        Ok(Self {
            prime,
            q,
            trace,
            discriminant,
        })
    }

    fn family(self) -> Result<legendre_cm::Family> {
        legendre_cm::Family::from_frobenius(self.prime, self.trace, i128::from(self.q))
    }
}

pub(crate) fn cm_family(model: &mut Model, curve: &Curve) -> Result<legendre_cm::Family> {
    Frobenius::ordinary(model, curve)?.family()
}

pub(crate) struct LocalOrder {
    frobenius: Frobenius,
    level: usize,
    bound: usize,
    modulus: u64,
    scalar: u64,
    trace: i128,
    norm: i128,
    discriminant: i128,
}

impl LocalOrder {
    pub(crate) fn compute(model: &mut Model, curve: &Curve, prime: usize) -> Result<Self> {
        let ell = elliptic::validate_ell(model, curve, prime)?;
        let frobenius = Frobenius::ordinary(model, curve)?;
        let bound = if ell == 2 {
            conductor::frobenius_conductor(model, frobenius.discriminant)?.trailing_zeros() as usize
        } else {
            // Odd prime valuations in the fundamental discriminant
            // are at most one.
            let mut discriminant = frobenius.discriminant.unsigned_abs();
            let square = u128::from(ell).pow(2);
            let mut bound = 0;
            while discriminant.is_multiple_of(square) {
                model.spend(1)?;
                bound += 1;
                discriminant /= square;
            }
            bound
        };
        let mut observed = (0, 1, 0);
        if bound > 0 {
            let mut order = 1u64;
            let mut division = DivisionPolynomials::new(model, curve, None)?;
            for level in 1..=bound {
                order = order
                    .checked_mul(ell)
                    .ok_or_else(|| Error::plain("geometric torsion order ell^k exceeds u64"))?;
                let candidates = scalar_candidates(frobenius.trace, ell, order);
                let mut next = None;
                for candidate in candidates {
                    if candidate % observed.1 != observed.2
                        || divided_characteristic(frobenius.q, frobenius.trace, candidate, order)
                            .is_none()
                    {
                        continue;
                    }
                    let proof = Certificate::build(
                        model,
                        curve,
                        &mut division,
                        frobenius.q,
                        ell,
                        order,
                        candidate,
                    )?;
                    if proof.checks() == (true, true) {
                        next = Some(candidate);
                        break;
                    }
                }
                let Some(candidate) = next else { break };
                observed = (level, order, candidate);
            }
        }
        let (level, modulus, scalar) = observed;
        let (trace, norm) = divided_characteristic(frobenius.q, frobenius.trace, scalar, modulus)
            .ok_or_else(|| {
            Error::plain("Frobenius division produced nonintegral local order data")
        })?;
        let modulus_square = i128::from(modulus).pow(2);
        if frobenius.discriminant % modulus_square != 0 {
            return Err(Error::plain(
                "Frobenius division produced nonintegral local order data",
            ));
        }
        Ok(Self {
            frobenius,
            level,
            bound,
            modulus,
            scalar,
            trace,
            norm,
            discriminant: frobenius.discriminant / modulus_square,
        })
    }

    pub(crate) fn report(&self, name: &str, prime: usize) -> String {
        format!(
            "rvEllLocalOrder[{name}, {prime}] = r={}, bound={}, m={}, a={}, trace={}, norm={}, discriminant={}",
            self.level,self.bound,self.modulus,self.scalar,self.trace,self.norm,self.discriminant,
        )
    }

    pub(crate) fn generator(
        &self,
        model: &mut Model,
        family: legendre_cm::Family,
    ) -> Result<crate::model::Value> {
        if family != self.frobenius.family()? {
            return Err(Error::plain(
                "local endomorphism generator requires the curve's Frobenius field and embedding",
            ));
        }
        let constant = self.frobenius.trace - 2 * i128::from(self.scalar);
        let denominator = 2 * i128::from(self.modulus);
        rv_qp_quad::persist(model, QuadElement::new(constant, 1, denominator)?)
    }
}

pub(crate) fn edge(
    model: &mut Model,
    name: &str,
    source: &Curve,
    target: &Curve,
    prime: usize,
) -> Result<String> {
    let src = LocalOrder::compute(model, source, prime)?;
    let dst = LocalOrder::compute(model, target, prime)?;
    let orientation = orientation(&src, &dst)?.name();
    Ok(format!(
        "rvEllEdge[{name}] = {orientation}, r_source={}, r_target={}",
        src.level, dst.level
    ))
}

#[derive(Clone, Copy)]
enum Orientation {
    Ascending,
    Horizontal,
    Descending,
}

impl Orientation {
    fn name(self) -> &'static str {
        match self {
            Self::Ascending => "ascending",
            Self::Horizontal => "horizontal",
            Self::Descending => "descending",
        }
    }
}

fn orientation(src: &LocalOrder, dst: &LocalOrder) -> Result<Orientation> {
    if src.frobenius != dst.frobenius {
        return Err(Error::plain(
            "isogeny endpoints have different Frobenius characteristic polynomials",
        ));
    }
    match dst.level as i128 - src.level as i128 {
        -1 => Ok(Orientation::Descending),
        0 => Ok(Orientation::Horizontal),
        1 => Ok(Orientation::Ascending),
        _ => Err(Error::plain(
            "a prime-degree edge changes the local conductor by at most one",
        )),
    }
}

pub(crate) fn neighborhood(
    model: &mut Model,
    name: &str,
    source: &Curve,
    prime: usize,
) -> Result<String> {
    let src = LocalOrder::compute(model, source, prime)?;
    let kernels = elliptic::enumerate_kernels(model, source, prime)?;
    let source_name: Rc<str> = name.into();
    let quotient_name: Rc<str> = format!("{name}::neighbor").into();
    let mut counts = [0; 3];
    let mut edges = String::new();
    for polynomial in kernels {
        let kernel = Rc::new(Kernel::from_certified(
            model,
            source,
            prime as u64,
            polynomial,
        )?);
        let (_, target) = Isogeny::new(
            model,
            source,
            Rc::clone(&source_name),
            Rc::clone(&quotient_name),
            Rc::clone(&kernel),
        )?;
        let dst = LocalOrder::compute(model, &target, prime)?;
        let direction = orientation(&src, &dst)?;
        counts[direction as usize] += 1;
        if !edges.is_empty() {
            edges.push_str(", ");
        }
        let j = target.j_invariant(model)?;
        edges.push_str(&format!(
            "(kernel={:?}, a={:?}, b={:?}, j={j:?}, {}, r={})",
            kernel.polynomial,
            target.a,
            target.b,
            direction.name(),
            dst.level
        ));
        if edges.len() > crate::model::MAX_SHOWN_BYTES {
            return Err(Error::plain(
                "isogeny neighborhood exceeds the display size limit",
            ));
        }
    }
    let [up, horizontal, down] = counts;
    let output = format!(
        "rvEllNeighborhood[{name}, {prime}] = r_source={}, ascending={up}, horizontal={horizontal}, descending={down}, edges=[{edges}]",
        src.level
    );
    if output.len() > crate::model::MAX_SHOWN_BYTES {
        return Err(Error::plain(
            "isogeny neighborhood exceeds the display size limit",
        ));
    }
    Ok(output)
}
