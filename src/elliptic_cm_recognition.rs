use std::borrow::Cow;
use std::cmp::Ordering;
use std::rc::Rc;

use crate::ec_morphism::{AffineMap, Morphism, Rational};
use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::elliptic::Curve;
use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing};
use crate::legendre_cm::Family;
use crate::model::{Model, Value};
use crate::rv_qp_quad::{self, Element as QuadElement};
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Expression {
    Frobenius(Label, bool),
    Negate(Box<ParsedExpr>),
    Sum(Box<[ParsedExpr; 2]>),
    Extend(Label, Label, Box<ParsedExpr>),
    Recognize(Label, Box<ParsedExpr>),
}

pub(crate) enum Command {
    Certify(Label, Box<ParsedExpr>),
}

#[derive(Clone, Copy, Eq, PartialEq)]
struct Coordinates {
    constant: i128,
    frobenius: i128,
    denominator: i128,
}

impl Coordinates {
    // The trace pairing in (1, pi) is [[2,t],[t,t²-2q]].
    fn recover(family: Family, q: u64, trace: i128, trace_product: i128) -> Result<Self> {
        let t = family.trace;
        let second_diagonal = subtract(multiply(t, t)?, 2 * i128::from(q))?;
        let a = subtract(
            multiply(second_diagonal, trace)?,
            multiply(t, trace_product)?,
        )?;
        let b = subtract(multiply(2, trace_product)?, multiply(t, trace)?)?;
        let denominator = family.discriminant.unsigned_abs();
        let common = gcd(gcd(denominator, a.unsigned_abs()), b.unsigned_abs());
        let divisor = -i128::try_from(common).map_err(|_| overflow())?;
        Ok(Self {
            constant: a.checked_div(divisor).ok_or_else(overflow)?,
            frobenius: b.checked_div(divisor).ok_or_else(overflow)?,
            denominator: i128::try_from(denominator / common).map_err(|_| overflow())?,
        })
    }

    fn field_element(self, model: &mut Model, family: Family) -> Result<Value> {
        let a = add(
            multiply(2, self.constant)?,
            multiply(family.trace, self.frobenius)?,
        )?;
        let d = multiply(2, self.denominator)?;
        rv_qp_quad::persist(model, QuadElement::new(a, self.frobenius, d)?)
    }
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn add(a: i128, b: i128) -> Result<i128> {
    a.checked_add(b).ok_or_else(overflow)
}

fn subtract(a: i128, b: i128) -> Result<i128> {
    a.checked_sub(b).ok_or_else(overflow)
}

fn multiply(a: i128, b: i128) -> Result<i128> {
    a.checked_mul(b).ok_or_else(overflow)
}

fn overflow() -> Error {
    Error::plain("elliptic CM trace arithmetic exceeds the exact integer range")
}

pub(crate) fn frobenius(model: &mut Model, name: Rc<str>, curve: &Curve) -> Result<Morphism> {
    let q = curve.ring.size_of_residue_field()?;
    frobenius_over(model, name, curve, q)
}

pub(crate) fn frobenius_over(
    model: &mut Model,
    name: Rc<str>,
    curve: &Curve,
    q: u64,
) -> Result<Morphism> {
    let x_bound = usize::try_from(q)
        .ok()
        .and_then(|q| q.checked_add(1))
        .ok_or_else(|| Error::plain("Frobenius polynomial degree exceeds the exact range"))?;
    let field = curve.field();
    field.check_polynomial(x_bound)?;
    model.spend(x_bound)?;
    let mut x = vec![field.zero(); x_bound];
    x[x_bound - 1] = field.one();
    let exponent = (q - 1) / 2;
    let y_bound = usize::try_from(exponent)
        .ok()
        .and_then(|power| power.checked_mul(3))
        .and_then(|degree| degree.checked_add(1))
        .ok_or_else(|| Error::plain("Frobenius polynomial degree exceeds the exact range"))?;
    let y = hw_algebra::pow_polynomial(model, &field, curve.polynomial(), exponent, y_bound)?;
    Morphism::from_affine(
        model,
        name.clone(),
        name,
        curve,
        curve,
        AffineMap {
            x: Rational::polynomial(field, x),
            y_multiplier: Rational::polynomial(field, y),
        },
    )
}

pub(crate) fn negate(
    model: &mut Model,
    source: &Curve,
    target: &Curve,
    map: &Morphism,
) -> Result<Morphism> {
    let Some((x, y)) = coordinates(map) else {
        return Ok(map.clone());
    };
    let field = source.field();
    let y = y.scale(model, field, &field.signed(-1))?;
    Morphism::from_affine(
        model,
        map.source.clone(),
        map.target.clone(),
        source,
        target,
        AffineMap {
            x: x.clone(),
            y_multiplier: y,
        },
    )
}

fn coordinates(map: &Morphism) -> Option<(&Rational, &Rational)> {
    map.coordinates
        .as_ref()
        .map(|coordinates| (&coordinates.x, &coordinates.y_multiplier))
}

pub(crate) fn sum(
    model: &mut Model,
    source: &Curve,
    target: &Curve,
    a: &Morphism,
    b: &Morphism,
) -> Result<Morphism> {
    if a.source != b.source || a.target != b.target {
        return Err(Error::plain(
            "elliptic map addition requires the same named source and target",
        ));
    }
    let (Some(ca), Some(cb)) = (coordinates(a), coordinates(b)) else {
        return Ok(if a.coordinates.is_none() { b } else { a }.clone());
    };
    let law = GroupLaw {
        source,
        target,
        a: ca,
        b: cb,
    };
    let Some(slope) = law.slope(model)? else {
        return Ok(Morphism::zero(a.source.clone(), a.target.clone()));
    };
    let x = law.x_sum(model, &slope)?;
    let y = law.y_sum(model, &slope, &x)?;
    Morphism::from_affine(
        model,
        a.source.clone(),
        a.target.clone(),
        source,
        target,
        AffineMap { x, y_multiplier: y },
    )
}

struct GroupLaw<'a> {
    source: &'a Curve,
    target: &'a Curve,
    a: (&'a Rational, &'a Rational),
    b: (&'a Rational, &'a Rational),
}

impl GroupLaw<'_> {
    fn slope(&self, model: &mut Model) -> Result<Option<Rational>> {
        let field = self.source.field();
        let (x, y) = self.a;
        let (u, v) = self.b;
        if x != u {
            let dx = u.subtract(model, field, x)?;
            let dy = v.subtract(model, field, y)?;
            return Ok(Some(dy.divide(model, field, &dx)?));
        }
        if y != v {
            // Certified points with the same x and unequal y are opposites.
            return Ok(None);
        }
        let numerator = self.tangent_numerator(model)?;
        let denominator = self.source_polynomial().multiply(model, field, y)?;
        let denominator = denominator.scale(model, field, &field.signed(2))?;
        Ok(Some(numerator.divide(model, field, &denominator)?))
    }

    fn x_sum(&self, model: &mut Model, slope: &Rational) -> Result<Rational> {
        let field = self.source.field();
        let slope_squared = slope.multiply(model, field, slope)?;
        let sum_x = self.a.0.add(model, field, self.b.0)?;
        self.source_polynomial()
            .multiply(model, field, &slope_squared)?
            .subtract(model, field, &sum_x)
    }

    fn source_polynomial(&self) -> Rational {
        Rational::polynomial(self.source.field(), self.source.polynomial().to_vec())
    }

    fn tangent_numerator(&self, model: &mut Model) -> Result<Rational> {
        let field = self.source.field();
        let (x, _) = self.a;
        let a = Rational::polynomial(field, poly::trim(vec![self.target.a.clone()]));
        let three = Rational::polynomial(field, field.polynomial(&[3]));
        x.multiply(model, field, x)?
            .multiply(model, field, &three)?
            .add(model, field, &a)
    }

    fn y_sum(&self, model: &mut Model, slope: &Rational, sum_x: &Rational) -> Result<Rational> {
        let field = self.source.field();
        let difference = self.a.0.subtract(model, field, sum_x)?;
        slope
            .multiply(model, field, &difference)?
            .subtract(model, field, self.a.1)
    }
}

pub(crate) fn lang(model: &mut Model, name: Rc<str>, curve: &Curve) -> Result<Morphism> {
    let identity = Morphism::identity(model, name.clone(), curve)?;
    let pi = frobenius(model, name, curve)?;
    let negative = negate(model, curve, curve, &pi)?;
    sum(model, curve, curve, &identity, &negative)
}

pub(crate) fn extend(
    model: &mut Model,
    map: &Morphism,
    source_name: Rc<str>,
    target_name: Rc<str>,
) -> Result<Morphism> {
    let old_source = model.ec_curve(&map.source)?;
    let old_target = model.ec_curve(&map.target)?;
    let new_source = model.ec_curve(&source_name)?;
    let new_target = model.ec_curve(&target_name)?;
    extend_over(
        model,
        map,
        [source_name, target_name],
        [&old_source, &old_target],
        [&new_source, &new_target],
    )
}

pub(crate) fn extend_over(
    model: &mut Model,
    map: &Morphism,
    [source_name, target_name]: [Rc<str>; 2],
    [old_source, old_target]: [&Curve; 2],
    [new_source, new_target]: [&Curve; 2],
) -> Result<Morphism> {
    if old_source.ring_name != old_target.ring_name || new_source.ring_name != new_target.ring_name
    {
        return Err(Error::plain(
            "elliptic base extension requires common fields at both endpoints",
        ));
    }
    let source = &old_source.ring;
    let target = &new_source.ring;
    let ab = [
        old_source.a.clone(),
        old_source.b.clone(),
        old_target.a.clone(),
        old_target.b.clone(),
    ];
    let embedded = source.embed_residue_coefficients(model, target, &ab)?;
    if embedded
        != [
            new_source.a.clone(),
            new_source.b.clone(),
            new_target.a.clone(),
            new_target.b.clone(),
        ]
    {
        return Err(Error::plain(
            "elliptic base extension does not give the declared endpoint curves",
        ));
    }
    let Some((x, y)) = coordinates(map) else {
        return Ok(Morphism::zero(source_name, target_name));
    };
    let field = new_source.field();
    let mut transport = |coordinate: &Rational| -> Result<Rational> {
        field.check_polynomial(coordinate.numerator.len().max(coordinate.denominator.len()))?;
        let numerator = source.embed_residue_coefficients(model, target, &coordinate.numerator)?;
        let denominator =
            source.embed_residue_coefficients(model, target, &coordinate.denominator)?;
        // An embedding preserves coprimality and the leading unit.
        Ok(Rational {
            numerator,
            denominator,
        })
    };
    let coordinates = AffineMap {
        x: transport(x)?,
        y_multiplier: transport(y)?,
    };
    Morphism::from_affine(
        model,
        source_name,
        target_name,
        new_source,
        new_target,
        coordinates,
    )
}

pub(crate) struct Recognition {
    family: Family,
    coordinates: Coordinates,
    norm: i128,
    trace: i128,
    trace_product: i128,
}

impl Recognition {
    pub(crate) fn check(model: &mut Model, family: Family, map: &Morphism) -> Result<Self> {
        if map.source != map.target {
            return Err(Error::plain(
                "elliptic CM recognition requires an endomorphism of a named curve",
            ));
        }
        let curve = model.ec_curve(&map.source)?;
        let (q, t) = curve.trace(model)?;
        let t = i128::from(t);
        let actual = Family::from_frobenius(curve.ring.prime() as usize, t, i128::from(q))?;
        if actual != family {
            return Err(Error::plain(
                "CM field does not match the endomorphism's ordinary Frobenius data",
            ));
        }
        let identity = Morphism::identity(model, map.source.clone(), &curve)?;
        let pi = frobenius(model, map.source.clone(), &curve)?;
        let negative = negate(model, &curve, &curve, map)?;
        let norm = map.degree() as i128;
        let lang_degree = degree_of_sum(model, &curve, &identity, &negative)?;
        let trace = subtract(add(1, norm)?, lang_degree)?;
        // Degree polarization gives tr(alpha*pi_bar) from the actual sum map.
        // pi_bar = t-pi, so this also gives tr(alpha*pi) without forming a
        // composite of degree norm*q just to recover its trace.
        let sum_degree = degree_of_sum(model, &curve, map, &pi)?;
        let n_plus_q = add(norm, i128::from(q))?;
        let trace_product = subtract(add(multiply(t, trace)?, n_plus_q)?, sum_degree)?;
        let coordinates = Coordinates::recover(family, q, trace, trace_product)?;
        let name = &map.source;
        let denominator = Morphism::multiple(model, name.clone(), &curve, coordinates.denominator)?;
        let lhs = denominator.compose(model, &curve, map)?;
        let rational_part = Morphism::multiple(model, name.clone(), &curve, coordinates.constant)?;
        let pi_factor = Morphism::multiple(model, name.clone(), &curve, coordinates.frobenius)?;
        let pi_part = pi_factor.compose(model, &curve, &pi)?;
        if !sum_equals(model, &curve, &rational_part, &pi_part, &lhs)? {
            return Err(Error::plain(
                "CM coordinates failed the cleared rational-map identity",
            ));
        }
        Ok(Self {
            family,
            coordinates,
            norm,
            trace,
            trace_product,
        })
    }

    pub(crate) fn into_value(self, model: &mut Model) -> Result<Value> {
        self.coordinates.field_element(model, self.family)
    }

    pub(crate) fn matches_order_element(
        &self,
        pi: crate::legendre_cm::Element,
        element: crate::legendre_cm::Element,
    ) -> Result<bool> {
        let c = self.coordinates;
        Ok(add(c.constant, multiply(c.frobenius, pi.constant)?)?
            == multiply(c.denominator, element.constant)?
            && multiply(c.frobenius, pi.generator)? == multiply(c.denominator, element.generator)?)
    }

    pub(crate) fn report(&self, name: &str) -> String {
        format!(
            "rvECCMCertify[{name}] = degree={}, trace={}, trace_product={}, [D]alpha=[A]+[B]pi with (A,B,D)=({},{},{})",
            self.norm,
            self.trace,
            self.trace_product,
            self.coordinates.constant,
            self.coordinates.frobenius,
            self.coordinates.denominator,
        )
    }
}

fn degree_of_sum(model: &mut Model, curve: &Curve, a: &Morphism, b: &Morphism) -> Result<i128> {
    let (Some(ca), Some(cb)) = (coordinates(a), coordinates(b)) else {
        return Ok(a.degree().max(b.degree()) as i128);
    };
    let law = GroupLaw {
        source: curve,
        target: curve,
        a: ca,
        b: cb,
    };
    let Some(slope) = law.slope(model)? else {
        return Ok(0);
    };
    let x = law.x_sum(model, &slope)?;
    Ok((x.numerator.len().max(x.denominator.len()) - 1) as i128)
}

fn sum_equals(
    model: &mut Model,
    curve: &Curve,
    a: &Morphism,
    b: &Morphism,
    expected: &Morphism,
) -> Result<bool> {
    let (Some(ca), Some(cb)) = (coordinates(a), coordinates(b)) else {
        let nonzero = if a.coordinates.is_none() { b } else { a };
        return Ok(nonzero.compare(model, expected)? == Ordering::Equal);
    };
    let law = GroupLaw {
        source: curve,
        target: curve,
        a: ca,
        b: cb,
    };
    let Some(slope) = law.slope(model)? else {
        return Ok(expected.coordinates.is_none());
    };
    let Some((x, y)) = coordinates(expected) else {
        return Ok(false);
    };
    // Comparing both normalized chart coordinates verifies the cleared
    // identity without constructing or recertifying an additional map.
    let result_x = law.x_sum(model, &slope)?;
    model.spend(
        x.numerator.len().min(result_x.numerator.len())
            + x.denominator.len().min(result_x.denominator.len()),
    )?;
    if &result_x != x {
        return Ok(false);
    }
    let result_y = law.y_sum(model, &slope, x)?;
    model.spend(
        y.numerator.len().min(result_y.numerator.len())
            + y.denominator.len().min(result_y.denominator.len()),
    )?;
    Ok(&result_y == y)
}

impl Rational {
    fn polynomial(field: Field<'_>, numerator: Polynomial) -> Self {
        Self {
            numerator: poly::trim(numerator),
            denominator: vec![field.one()],
        }
    }

    fn add(&self, model: &mut Model, field: Field<'_>, other: &Self) -> Result<Self> {
        self.combine(model, field, other, false)
    }

    fn subtract(&self, model: &mut Model, field: Field<'_>, other: &Self) -> Result<Self> {
        self.combine(model, field, other, true)
    }

    fn combine(
        &self,
        model: &mut Model,
        field: Field<'_>,
        other: &Self,
        negative: bool,
    ) -> Result<Self> {
        if self.denominator == other.denominator {
            let numerator = poly::sum(model, field, &self.numerator, &other.numerator, negative)?;
            return Self::new(model, field, numerator, self.denominator.clone());
        }
        let (left, right) = (&self.denominator, &other.denominator);
        let (left_cofactor, right_cofactor) = cancel(model, field, left, right)?;
        let a = poly::multiply(model, field, &self.numerator, &right_cofactor)?;
        let b = poly::multiply(model, field, &other.numerator, &left_cofactor)?;
        let numerator = poly::sum(model, field, &a, &b, negative)?;
        let denominator = poly::multiply(model, field, left, &right_cofactor)?;
        Self::new(model, field, numerator, denominator)
    }

    fn divide(&self, model: &mut Model, field: Field<'_>, other: &Self) -> Result<Self> {
        if other.numerator.is_empty() {
            return Err(Error::plain("rational map has a zero denominator"));
        }
        if self.numerator.is_empty() {
            return Ok(Self::scalar(field, field.zero()));
        }
        // Cancel the pairs of like factors before multiplying across: the
        // individual rational coordinates already have coprime fractions.
        let (a, c) = cancel(model, field, &self.numerator, &other.numerator)?;
        let (b, d) = cancel(model, field, &self.denominator, &other.denominator)?;
        let numerator = poly::multiply(model, field, &a, &d)?;
        let denominator = poly::multiply(model, field, &b, &c)?;
        Self::new(model, field, numerator, denominator)
    }
}

type Cofactors<'a> = (Cow<'a, [Scalar]>, Cow<'a, [Scalar]>);

fn cancel<'a>(
    model: &mut Model,
    field: Field<'_>,
    left: &'a [Scalar],
    right: &'a [Scalar],
) -> Result<Cofactors<'a>> {
    if left.len().min(right.len()) > 1 {
        model.spend(left.len() + right.len())?;
        let common = poly::gcd(model, field, left.to_vec(), right.to_vec())?;
        if common.len() > 1 {
            return Ok((
                Cow::Owned(poly::exact_divide(model, field, left, &common)?),
                Cow::Owned(poly::exact_divide(model, field, right, &common)?),
            ));
        }
    }
    Ok((Cow::Borrowed(left), Cow::Borrowed(right)))
}
