use std::cmp::Ordering;
use std::rc::Rc;

use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::finite_field_factor::{factor_squarefree, quadratic_sign};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic_jacobian::{relative_frobenius_power, require_extension, Divisor};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod arithmetic;
pub(crate) mod k3;
pub(crate) mod resolution;

pub(crate) enum Expression {
    Coordinates(Label, Box<[ParsedExpr; 4]>),
    Exceptional(Box<ParsedExpr>, Box<[ParsedExpr; 4]>),
    DifferentialAdd(Box<[ParsedExpr; 3]>),
    Richelot(Label, Box<ParsedExpr>),
    RichelotDescend(Label, Box<ParsedExpr>),
    Base(BaseOperation, Label, Box<ParsedExpr>),
    Apply(Operation, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum BaseOperation {
    Extend,
    Descend,
    Frobenius(usize),
}

impl BaseOperation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Extend => "rvKummerExtend",
            Self::Descend => "rvKummerDescend",
            Self::Frobenius(_) => "rvKummerFrob",
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Project,
    Coordinate(usize),
    Liftable,
    Lift,
    Scale(i128),
    Smooth,
    Blowdown,
    IsExceptional,
    Direction(usize),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Project => "rvJacKummerPoint",
            Self::Coordinate(_) => "rvKummerCoord",
            Self::Liftable => "rvKummerLiftable",
            Self::Lift => "rvKummerLift",
            Self::Scale(_) => "rvKummerScale",
            Self::Smooth => "rvKummerSmooth",
            Self::Blowdown => "rvKummerBlowdown",
            Self::IsExceptional => "rvKummerIsExceptional",
            Self::Direction(_) => "rvKummerDirection",
        }
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Point([Scalar; 4]);

pub(crate) enum Recovery {
    Identity,
    Linear {
        x: Scalar,
        y_squared: Scalar,
    },
    Quadratic {
        u: Polynomial,
        a: Scalar,
        b: Scalar,
        c: Scalar,
    },
}

pub(crate) fn require_curve(curve: &Curve) -> Result<()> {
    if curve.ring.precision() != 1
        || curve.polynomial().len() != 6
        || curve.polynomial().last() != Some(&curve.ring.scalar(1))
    {
        return Err(Error::plain(
            "RVKummer requires a smooth monic quintic over a named residue field of precision one",
        ));
    }
    Ok(())
}

fn field(curve: &Curve) -> Field<'_> {
    Field { ring: &curve.ring }
}

impl Point {
    pub(crate) fn new(
        model: &mut Model,
        curve: &Curve,
        mut coordinates: [Scalar; 4],
    ) -> Result<Self> {
        require_curve(curve)?;
        let field = field(curve);
        let scale = coordinates
            .iter()
            .find(|c| !poly::zero(c))
            .ok_or_else(|| Error::plain("rvKummer requires a nonzero projective quadruple"))?;
        if scale != &field.one() {
            let inverse = field.inverse_unit(model, scale)?;
            for coefficient in &mut coordinates {
                *coefficient = field.multiply(model, coefficient, &inverse)?;
            }
        }
        let [k0, k1, k2, k3] = &coordinates;
        let valid = if poly::zero(k0) {
            field.multiply(model, k1, k3)? == field.multiply(model, k2, k2)?
        } else {
            let data = QuadraticData::new(model, curve, k1, k2)?;
            data.matrix(model, field, k2, k3)?
                .on_quartic(model, field)?
        };
        if !valid {
            return Err(Error::plain(
                "rvKummer coordinates do not lie on the Kummer quartic",
            ));
        }
        Ok(Self(coordinates))
    }

    pub(crate) fn identity(curve: &Curve) -> Self {
        let field = field(curve);
        Self([field.zero(), field.zero(), field.zero(), field.one()])
    }

    pub(crate) fn project(model: &mut Model, curve: &Curve, divisor: &Divisor) -> Result<Self> {
        require_curve(curve)?;
        let field = field(curve);
        let coordinates = match divisor.u.as_slice() {
            [_] => return Ok(Self::identity(curve)),
            [constant, _] => {
                let x = field.ring.negate(model, constant)?;
                let square = field.multiply(model, &x, &x)?;
                [field.zero(), field.one(), x, square]
            }
            [p, linear, _] => {
                let s = field.ring.negate(model, linear)?;
                let q = divisor.norm_quotient(model, curve)?;
                let pq2 = field.multiply(model, p, &q[2])?;
                let sum = field.add(model, &q[0], &pq2)?;
                let t = field.ring.negate(model, &sum)?;
                [field.one(), s, p.clone(), t]
            }
            _ => {
                return Err(Error::plain(
                    "Kummer projection requires a reduced genus-two divisor",
                ))
            }
        };
        Ok(Self(coordinates))
    }

    pub(crate) fn recover(&self, model: &mut Model, curve: &Curve) -> Result<Recovery> {
        require_curve(curve)?;
        let [k0, k1, k2, k3] = &self.0;
        let field = field(curve);
        if poly::zero(k0) {
            if poly::zero(k1) {
                return Ok(Recovery::Identity);
            }
            let y_squared = poly::evaluate(model, field, curve.polynomial(), k2)?;
            return Ok(Recovery::Linear {
                x: k2.clone(),
                y_squared,
            });
        }
        let data = QuadraticData::new(model, curve, k1, k2)?;
        let SquareMatrix { a, b, c } = data.matrix(model, field, k2, k3)?;
        Ok(Recovery::Quadratic { u: data.u, a, b, c })
    }

    pub(crate) fn is_node(&self, model: &mut Model, curve: &Curve) -> Result<bool> {
        Ok(match self.recover(model, curve)? {
            Recovery::Identity => true,
            Recovery::Linear { y_squared, .. } => poly::zero(&y_squared),
            Recovery::Quadratic { a, b, c, .. } => {
                poly::zero(&a) && poly::zero(&b) && poly::zero(&c)
            }
        })
    }

    pub(crate) fn liftable(&self, model: &mut Model, curve: &Curve) -> Result<bool> {
        let diagonal = match self.recover(model, curve)? {
            Recovery::Identity => return Ok(true),
            Recovery::Linear { y_squared, .. } => y_squared,
            Recovery::Quadratic { a, c, .. } => {
                if poly::zero(&a) {
                    c
                } else {
                    a
                }
            }
        };
        if poly::zero(&diagonal) {
            return Ok(true);
        }
        let field = field(curve);
        let q = field.ring.size_of_residue_field().ok();
        let sign = quadratic_sign(model, field, q, &field.polynomial(&[0, 1]), 1, &[diagonal])?;
        Ok(sign == [field.one()])
    }

    pub(crate) fn lift(&self, model: &mut Model, curve: &Curve) -> Result<Divisor> {
        let field = field(curve);
        match self.recover(model, curve)? {
            Recovery::Identity => Ok(Divisor::zero(curve)),
            Recovery::Linear { x, y_squared } => {
                let y = lift_root(model, field, &y_squared)?;
                Divisor::from_point(model, curve, &x, y)
            }
            Recovery::Quadratic { u, a, b, c } => {
                let v = if poly::zero(&a) {
                    vec![field.zero(), lift_root(model, field, &c)?]
                } else {
                    let constant = lift_root(model, field, &a)?;
                    let inverse = field.inverse_unit(model, &constant)?;
                    vec![constant, field.multiply(model, &b, &inverse)?]
                };
                Divisor::new(model, curve, u, v)
            }
        }
    }

    pub(crate) fn coordinates(&self) -> &[Scalar; 4] {
        &self.0
    }

    pub(crate) fn extend(&self, model: &mut Model, base: &Curve, target: &Curve) -> Result<Self> {
        require_extension(model, base, target)?;
        let coordinates = base
            .ring
            .embed_residue_coefficients(model, &target.ring, &self.0)?;
        // The embedding preserves the quartic and the first nonzero coordinate 1.
        Ok(Self(
            coordinates.try_into().expect("four embedded coordinates"),
        ))
    }

    pub(crate) fn descend(&self, model: &mut Model, source: &Curve, base: &Curve) -> Result<Self> {
        require_extension(model, base, source)?;
        let preimage = base.ring.residue_preimage(model, &source.ring)?;
        let coordinates = preimage
            .coefficients(model, &self.0)?
            .ok_or_else(|| Error::plain("Kummer coordinate lies outside the requested subfield"))?;
        Ok(Self(
            coordinates.try_into().expect("four descended coordinates"),
        ))
    }

    pub(crate) fn frobenius(
        &self,
        model: &mut Model,
        base: &Curve,
        target: &Curve,
        power: usize,
    ) -> Result<Self> {
        let p_power = relative_frobenius_power(model, base, target, power)?;
        let [k0, k1, k2, k3] = &self.0;
        Ok(Self([
            target.ring.frobenius(model, k0, p_power)?,
            target.ring.frobenius(model, k1, p_power)?,
            target.ring.frobenius(model, k2, p_power)?,
            target.ring.frobenius(model, k3, p_power)?,
        ]))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        model.spend(4 * (1 + self.0[0].len()))?;
        Ok(self.cmp(other))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(4 * (3 + self.0[0].len()))?;
        Ok(Value::Kummer(Rc::new(self)))
    }
}

struct QuadraticData {
    u: Polynomial,
    f0: Scalar,
    f1: Scalar,
    correction: Scalar,
    one_half: Scalar,
    s: Scalar,
}

struct SquareMatrix {
    a: Scalar,
    b: Scalar,
    c: Scalar,
}

impl QuadraticData {
    fn new(model: &mut Model, curve: &Curve, s: &Scalar, p: &Scalar) -> Result<Self> {
        let field = field(curve);
        let f = curve.polynomial();
        let u = vec![p.clone(), field.ring.negate(model, s)?, field.one()];
        let mut remainder = poly::remainder(model, field, f.to_vec(), &u)?;
        remainder.resize_with(2, || field.zero());
        let f1 = remainder.pop().unwrap();
        let f0 = remainder.pop().unwrap();
        // For q=(f−V²)/U, q2=f4+s and q1=f3+s q2−p. Thus
        // V's linear coefficient squared is t + f2 + s q1.
        let q2 = field.add(model, &f[4], s)?;
        let sq2 = field.multiply(model, s, &q2)?;
        let q1 = field.add(model, &f[3], &sq2)?;
        let q1 = field.subtract(model, &q1, p)?;
        let sq1 = field.multiply(model, s, &q1)?;
        let correction = field.add(model, &f[2], &sq1)?;
        let one_half = field.ring.scalar(field.ring.prime().div_ceil(2));
        Ok(Self {
            u,
            f0,
            f1,
            correction,
            one_half,
            s: s.clone(),
        })
    }

    fn matrix(
        &self,
        model: &mut Model,
        field: Field<'_>,
        p: &Scalar,
        t: &Scalar,
    ) -> Result<SquareMatrix> {
        let c = field.add(model, t, &self.correction)?;
        let sc = field.multiply(model, &self.s, &c)?;
        let twice_b = field.subtract(model, &self.f1, &sc)?;
        let b = field.multiply(model, &twice_b, &self.one_half)?;
        let pc = field.multiply(model, p, &c)?;
        let a = field.add(model, &self.f0, &pc)?;
        Ok(SquareMatrix { a, b, c })
    }
}

impl SquareMatrix {
    fn on_quartic(&self, model: &mut Model, field: Field<'_>) -> Result<bool> {
        Ok(field.multiply(model, &self.a, &self.c)? == field.multiply(model, &self.b, &self.b)?)
    }
}

fn lift_root(model: &mut Model, field: Field<'_>, value: &Scalar) -> Result<Scalar> {
    square_root(model, field, value)?.ok_or_else(|| {
        Error::plain("rvKummerLift has nonsquare recovery data over the named field")
    })
}

fn square_root(model: &mut Model, field: Field<'_>, value: &Scalar) -> Result<Option<Scalar>> {
    if poly::zero(value) {
        return Ok(Some(field.zero()));
    }
    let constant = field.ring.negate(model, value)?;
    let factors = factor_squarefree(model, field, vec![constant, field.zero(), field.one()])?;
    let [left, right] = factors.as_slice() else {
        return Ok(None);
    };
    let left = field.ring.negate(model, &left[0])?;
    let right = field.ring.negate(model, &right[0])?;
    Ok(Some(left.min(right)))
}

pub(crate) fn elements(model: &mut Model, curve: &Curve, limit: usize) -> Result<Vec<Value>> {
    let field = field(curve);
    let too_many = || Error::plain("Kummer surface exceeds the finite element limit");
    let q = usize::try_from(curve.ring.size_of_residue_field().map_err(|_| too_many())?)
        .map_err(|_| too_many())?;
    if q >= limit {
        return Err(too_many());
    }
    field.check_polynomial(q)?;
    let volume = q.checked_pow(3).ok_or_else(too_many)?;
    model.spend(volume + q * field.cells())?;
    let coordinates: Vec<_> = (0..q).map(|i| field.ring.coordinate_row(i)).collect();
    let mut points = Vec::new();
    let mut insert = |model: &mut Model, point: Point| -> Result<()> {
        if points.len() == limit {
            return Err(too_many());
        }
        points.push(point.persist(model)?);
        Ok(())
    };
    // These loops visit the normalized projective coordinates in lexicographic order.
    insert(model, Point::identity(curve))?;
    for x in &coordinates {
        let square = field.multiply(model, x, x)?;
        insert(model, Point([field.zero(), field.one(), x.clone(), square]))?;
    }
    for s in &coordinates {
        for p in &coordinates {
            let data = QuadraticData::new(model, curve, s, p)?;
            for t in &coordinates {
                if data.matrix(model, field, p, t)?.on_quartic(model, field)? {
                    insert(model, Point([field.one(), s.clone(), p.clone(), t.clone()]))?;
                }
            }
        }
    }
    Ok(points)
}
