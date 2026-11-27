use std::rc::Rc;

use crate::ec_morphism::{self, AffineMap, Morphism, Rational};
use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::elliptic::Curve;
use crate::error::{Error, Result};
use crate::finite_field_factor::factor_with_multiplicity;
use crate::hw_algebra::CoefficientRing;
use crate::legendre::closed_point_hensel;
use crate::model::{Model, Value};
use crate::rv_unram::{mul_mod, Ring};
use crate::syntax::Label;

use super::super::{period, quartic_prime};
use super::certificate;

#[derive(Clone, Copy)]
pub(crate) enum Invariant {
    Coordinate,
    Partner,
    J(usize),
    Twist,
}

impl Invariant {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Coordinate => "rvDworkFrickeZ",
            Self::Partner => "rvDworkFrickeT",
            Self::J(_) => "rvDworkFrickeJ",
            Self::Twist => "rvDworkFrickeTwist",
        }
    }

    pub(crate) fn apply(self, model: &mut Model, ring: &Ring, t: &[u64]) -> Result<Scalar> {
        let fiber = Fiber::new(model, ring, t)?;
        match self {
            Self::Coordinate => Ok(fiber.z),
            Self::Partner => Ok(fiber.partner),
            Self::Twist => Ok(fiber.twist),
            Self::J(1) => Ok(fiber.j1),
            Self::J(2) => Ok(fiber.j2),
            Self::J(_) => Err(Error::plain("Dwork Fricke j index must be 1 or 2")),
        }
    }
}

fn check_ring(ring: &Ring) -> Result<()> {
    if ring.precision() != 1 || ring.prime() == 2 {
        return Err(Error::plain(
            "Dwork Fricke specialization requires an odd residue field at precision one",
        ));
    }
    Ok(())
}

struct Fiber<'a> {
    field: Field<'a>,
    z: Scalar,
    j1: Scalar,
    j2: Scalar,
    partner: Scalar,
    twist: Scalar,
    u: Scalar,
}

impl<'a> Fiber<'a> {
    fn new(model: &mut Model, ring: &'a Ring, t: &[u64]) -> Result<Self> {
        check_ring(ring)?;
        let field = Field { ring };
        let t_plus_64 = ring.add(model, t, &field.signed(64))?;
        if poly::zero(t) || poly::zero(&t_plus_64) {
            return Err(Error::plain(
                "Dwork Fricke fiber requires t different from 0 and -64",
            ));
        }
        let inverse = ring.inverse(model, t)?;
        let partner = field.multiply(model, &field.signed(4096), &inverse)?;
        let u = parameter_u(model, ring, t)?;
        let twist = scale_fraction(model, ring, t, -1, 32)?;
        let j1 = j_value(model, ring, t, &inverse, 256, 2)?;
        let j2 = j_value(model, ring, t, &inverse, 16, 1)?;
        let numerator = ring.power(model, &t_plus_64, 2)?;
        let v = field.multiply(model, &numerator, &inverse)?;
        let z = ring.inverse(model, &v)?;
        Ok(Self {
            field,
            z,
            j1,
            j2,
            partner,
            twist,
            u,
        })
    }

    fn long_models(&self, model: &mut Model) -> Result<[[Scalar; 2]; 2]> {
        let f = self.field;
        let a = f.multiply(model, &f.signed(-2), &self.u)?;
        let u_minus_4 = f.subtract(model, &self.u, &f.signed(4))?;
        let b = f.multiply(model, &self.u, &u_minus_4)?;
        Ok([[self.u.clone(), self.u.clone()], [a, b]])
    }

    fn short_models(&self, model: &mut Model) -> Result<[[Scalar; 2]; 2]> {
        require_short(self.field.ring)?;
        let [left, right] = self.long_models(model)?;
        Ok([
            short(model, self.field, &left)?,
            short(model, self.field, &right)?,
        ])
    }

    fn traces(&self, model: &mut Model) -> Result<(u64, i64, i64, i64)> {
        let ring = self.field.ring;
        let q = ring.size_of_residue_field()?;
        let [left, right] = self.long_models(model)?;
        let partner_u = parameter_u(model, ring, &self.partner)?;
        // Verify the actual quadratic twist of the long model at W_2(t).
        let twist2 = ring.power(model, &self.twist, 2)?;
        if right[0] != ring.multiply(model, &partner_u, &self.twist)?
            || right[1] != ring.multiply(model, &partner_u, &twist2)?
        {
            return Err(certificate());
        }
        let first = trace(model, self.field, &left)?;
        let second = trace(model, self.field, &right)?;
        let third = trace(model, self.field, &[partner_u.clone(), partner_u])?;
        let character = ring.power(model, &self.twist, (q - 1) / 2)?;
        let sign = if character == self.field.one() {
            1
        } else if character == self.field.signed(-1) {
            -1
        } else {
            return Err(certificate());
        };
        if first != second || second != sign * third {
            return Err(certificate());
        }
        Ok((q, first, third, sign))
    }
}

pub(crate) fn fiber_report(model: &mut Model, ring: &Ring, t: &[u64]) -> Result<String> {
    let fiber = Fiber::new(model, ring, t)?;
    let [[a, b], [aa, bb]] = fiber.long_models(model)?;
    let (q, trace, partner_trace, character) = fiber.traces(model)?;
    Ok(format!(
        "rvDworkFrickeFiber = q={q}, t={t:?}, z={:?}, W2(t)={:?}, j1={:?}, j2={:?}, long E=({a:?}, {b:?}), long E' = ({aa:?}, {bb:?}), trace E=trace E'={trace}, trace E(W2(t))={partner_trace}, twist={:?}, twist character={character}",
        fiber.z, fiber.partner, fiber.j1, fiber.j2, fiber.twist,
    ))
}

fn parameter_u(model: &mut Model, ring: &Ring, t: &[u64]) -> Result<Scalar> {
    let shifted = scale_fraction(model, ring, t, 1, 16)?;
    let four = Field { ring }.signed(4);
    ring.add(model, &shifted, &four)
}

fn scale_fraction(
    model: &mut Model,
    ring: &Ring,
    t: &[u64],
    numerator: i128,
    denominator: i128,
) -> Result<Scalar> {
    let f = Field { ring };
    let inverse = ring.inverse(model, &f.signed(denominator))?;
    let quot = ring.multiply(model, &inverse, &f.signed(numerator))?;
    ring.multiply(model, t, &quot)
}

fn j_value(
    model: &mut Model,
    ring: &Ring,
    t: &[u64],
    inverse: &[u64],
    shift: i128,
    power: u64,
) -> Result<Scalar> {
    let shifted = ring.add(model, t, &Field { ring }.signed(shift))?;
    let cube = ring.power(model, &shifted, 3)?;
    let inverse = ring.power(model, inverse, power)?;
    ring.multiply(model, &cube, &inverse)
}

fn trace(model: &mut Model, field: Field<'_>, coefficients: &[Scalar; 2]) -> Result<i64> {
    let [a, b] = coefficients;
    let cubic = [field.zero(), b.clone(), a.clone(), field.one()];
    field
        .ring
        .residue_quadratic_trace(model, |model, x| poly::evaluate(model, field, &cubic, x))
}

fn require_short(ring: &Ring) -> Result<()> {
    if ring.prime() < 5 {
        return Err(Error::plain(
            "typed short Fricke elliptic curves require p at least 5",
        ));
    }
    Ok(())
}

fn short(model: &mut Model, f: Field<'_>, [a, b]: &[Scalar; 2]) -> Result<[Scalar; 2]> {
    let s = scale_fraction(model, f.ring, a, 1, 3)?;
    let s2 = f.multiply(model, &s, &s)?;
    let a2_over_3 = f.multiply(model, a, &s)?;
    let short_a = f.subtract(model, b, &a2_over_3)?;
    let b_over = f.subtract(model, &s2, b)?;
    let reduced = f.add(model, &s2, &b_over)?;
    let short_b = f.multiply(model, &s, &reduced)?;
    Ok([short_a, short_b])
}

pub(crate) fn define(
    model: &mut Model,
    names: [Label; 2],
    ring_name: Rc<str>,
    ring: Rc<Ring>,
    t: &[u64],
) -> Result<()> {
    let fiber = Fiber::new(model, &ring, t)?;
    let [[a, b], [aa, bb]] = fiber.short_models(model)?;
    let e = Curve::new(model, ring_name.clone(), ring.clone(), a, b)?;
    let ep = Curve::new(model, ring_name, ring.clone(), aa, bb)?;
    if e.j_invariant(model)? != fiber.j2 || ep.j_invariant(model)? != fiber.j1 {
        return Err(certificate());
    }
    model.define_ec_pair(names, [e, ep])
}

pub(crate) fn map(
    model: &mut Model,
    source: Rc<str>,
    target: Rc<str>,
    t: &[u64],
    reverse: bool,
) -> Result<Value> {
    let original = model.ec_curve(&source)?;
    let quotient = model.ec_curve(&target)?;
    let f = original.field();
    let fiber = Fiber::new(model, f.ring, t)?;
    let [[a, b], [aa, bb]] = fiber.short_models(model)?;
    if original.ring_name != quotient.ring_name
        || original.a != a
        || original.b != b
        || quotient.a != aa
        || quotient.b != bb
    {
        return Err(Error::plain(
            "elliptic curves do not match the ordered Fricke models for t",
        ));
    }
    let s = scale_fraction(model, f.ring, &fiber.u, 1, 3)?;
    let (shift, kernel_b, x_scale, y_scale) = if reverse {
        let [_, [_, b]] = fiber.long_models(model)?;
        (f.multiply(model, &s, &f.signed(2))?, b, 4, 8)
    } else {
        (f.ring.negate(model, &s)?, fiber.u, 1, 1)
    };
    // With w=X+shift the quotient has X' = (X+b/w)/sx and
    // Y' = Y*(1-b/w²)/sy. The reverse quotient is rescaled by (4,8).
    let field = f.ring;
    let shift2 = f.multiply(model, &shift, &shift)?;
    let double_shift = f.multiply(model, &shift, &f.signed(2))?;
    let y_constant = f.subtract(model, &shift2, &kernel_b)?;
    let rx = Rational::new(
        model,
        f,
        vec![kernel_b, shift.clone(), f.one()],
        vec![shift, f.one()],
    )?;
    let ry = Rational::new(
        model,
        f,
        vec![y_constant, double_shift.clone(), f.one()],
        vec![shift2, double_shift, f.one()],
    )?;
    let x_scale = field.inverse(model, &f.signed(x_scale))?;
    let y_scale = field.inverse(model, &f.signed(y_scale))?;
    let coordinates = AffineMap {
        x: rx.scale(model, f, &x_scale)?,
        y_multiplier: ry.scale(model, f, &y_scale)?,
    };
    let morphism = if reverse {
        Morphism::from_affine(model, target, source, &quotient, &original, coordinates)?
    } else {
        Morphism::from_affine(model, source, target, &original, &quotient, coordinates)?
    };
    if morphism.degree() != 2 {
        return Err(certificate());
    }
    ec_morphism::persist(model, morphism)
}

enum Parameters {
    Split(Vec<(Scalar, usize)>),
    Nonsplit(Polynomial),
}

fn parameters(model: &mut Model, ring: &Ring, z: &[u64]) -> Result<Parameters> {
    check_ring(ring)?;
    if poly::zero(z) {
        return Err(Error::plain("Dwork Fricke base requires z nonzero"));
    }
    let f = Field { ring };
    let inverse = ring.inverse(model, z)?;
    let coefficient = f.subtract(model, &f.signed(128), &inverse)?;
    let quadratic = vec![f.signed(4096), coefficient, f.one()];
    let factors = factor_with_multiplicity(model, f, &quadratic)?;
    if factors.len() == 1 && factors[0].0 == 2 {
        return Ok(Parameters::Nonsplit(quadratic));
    }
    let mut roots = factors
        .into_iter()
        .map(|(degree, multiplicity, factor)| {
            if degree != 1 {
                return Err(certificate());
            }
            Ok((ring.negate(model, &factor[0])?, multiplicity))
        })
        .collect::<Result<Vec<_>>>()?;
    roots.sort_unstable();
    if roots
        .iter()
        .map(|(_, multiplicity)| multiplicity)
        .sum::<usize>()
        != 2
    {
        return Err(certificate());
    }
    Ok(Parameters::Split(roots))
}

pub(crate) fn base_report(model: &mut Model, ring: &Ring, z: &[u64]) -> Result<String> {
    match parameters(model, ring, z)? {
        Parameters::Split(parameters) => {
            let rows = parameters
                .into_iter()
                .map(|(t, multiplicity)| {
                    let fiber = Fiber::new(model, ring, &t)?;
                    if fiber.z != z {
                        return Err(certificate());
                    }
                    Ok(format!(
                        "t={t:?} (multiplicity {multiplicity}), j1={:?}, j2={:?}",
                        fiber.j1, fiber.j2
                    ))
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(format!(
                "rvDworkFrickeBase = z={z:?}, split, {}",
                rows.join(" | ")
            ))
        }
        Parameters::Nonsplit(quadratic) => {
            let f = Field { ring };
            let q = ring.size_of_residue_field()?;
            let t = vec![f.zero(), f.one()];
            let inverse = poly::inverse_mod(model, f, &t, &quadratic)?;
            let partner = poly::scale(model, f, &inverse, &f.signed(4096))?;
            let conjugate = poly::power_mod(model, f, &t, q, &quadratic)?;
            let j1 = j_mod(model, f, &quadratic, &t, &inverse, 256, 2)?;
            let j2 = j_mod(model, f, &quadratic, &t, &inverse, 16, 1)?;
            if conjugate != partner
                || poly::power_mod(model, f, &j1, q, &quadratic)? != j2
                || poly::power_mod(model, f, &j2, q, &quadratic)? != j1
            {
                return Err(certificate());
            }
            let detail = if j1 == j2 {
                "coincident j-values"
            } else {
                "distinct j-values exchanged"
            };
            Ok(format!(
                "rvDworkFrickeBase = z={z:?}, nonsplit quadratic={quadratic:?}, j1={j1:?}, j2={j2:?}, Frobenius t->W2(t), {detail}"
            ))
        }
    }
}

fn j_mod(
    model: &mut Model,
    f: Field<'_>,
    quadratic: &[Scalar],
    t: &[Scalar],
    inverse: &[Scalar],
    shift: i128,
    power: u64,
) -> Result<Polynomial> {
    let shifted = poly::sum(model, f, t, &[f.signed(shift)], false)?;
    let numerator = poly::power_mod(model, f, &shifted, 3, quadratic)?;
    let denominator = poly::power_mod(model, f, inverse, power, quadratic)?;
    poly::multiply_mod(model, f, &numerator, &denominator, quadratic)
}

pub(crate) fn compare(
    model: &mut Model,
    prime: usize,
    psi: i128,
    precision: usize,
) -> Result<String> {
    let p = quartic_prime(model, prime, psi)?;
    if psi.rem_euclid(i128::from(p.get())) == 0 {
        return Err(Error::plain(
            "Dwork Fricke comparison requires psi nonzero modulo p",
        ));
    }
    let ring = Ring::new(model, prime, 1, vec![0, 1])?;
    let f = Field { ring: &ring };
    let cross = f.multiply(model, &f.signed(4), &f.signed(psi))?;
    let power = ring.power(model, &cross, 4)?;
    let z = ring.inverse(model, &power)?;
    let unit = period::prime_unit_data(model, p, psi, precision)?;
    let known = unit.unit();
    let head = format!(
        "rvDworkFrickeUnit[{prime}, {psi}, {precision}] = z={}, mod {}, Hasse–Witt={}, Dwork unit={}",
        z[0], unit.modulus, unit.hasse_witt, known.map_or_else(|| "none".into(), |i| i.to_string()),
    );
    let Parameters::Split(roots) = parameters(model, &ring, &z)? else {
        return Ok(format!(
            "{head}, nonsplit Fricke parameters over the prime field"
        ));
    };
    let mut results = Vec::new();
    for (t, _) in roots {
        let fiber = Fiber::new(model, &ring, &t)?;
        let (q, trace, _, _) = fiber.traces(model)?;
        let elliptic = closed_point_hensel(model, ring.prime(), precision, trace, q)?;
        let square = elliptic.map(|a| mul_mod(a, a, unit.modulus));
        let result = match (known, square) {
            (None, None) => "both nonordinary",
            (Some(a), Some(b)) if a == b => "equal",
            (Some(a), Some(b)) if a == unit.modulus - b => "opposite residues",
            _ => "different",
        };
        let display = |value: Option<u64>| value.map_or_else(|| "none".into(), |i| i.to_string());
        results.push(format!(
            "t={}, trace={trace}, elliptic unit={}, square={}, comparison={result}",
            t[0],
            display(elliptic),
            display(square)
        ));
    }
    Ok(format!("{head}, {}", results.join(" | ")))
}
