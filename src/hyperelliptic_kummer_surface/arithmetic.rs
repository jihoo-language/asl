use std::rc::Rc;

use super::{field, square_root, Point, Recovery};
use crate::ec_poly::{self as poly, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic_jacobian::Divisor;
use crate::model::Model;
use crate::richelot::correspondence::Correspondence;

pub(crate) fn scale(model: &mut Model, curve: &Curve, n: i128, point: &Point) -> Result<Point> {
    if n == 0 {
        return Ok(Point::identity(curve));
    }
    let lift = GeometricLift::recover(model, curve, point)?;
    let twist = Twist::new(model, curve, lift.preferred_radicand(curve))?;
    let divisor = twist.divisor(model, lift)?;
    let multiple = divisor.scale(model, &twist.curve, n)?;
    project_scaled(model, &twist.curve, curve, &multiple, &twist.k)
}

pub(crate) fn differential_add(
    model: &mut Model,
    curve: &Curve,
    left: &Point,
    right: &Point,
    difference: &Point,
) -> Result<Point> {
    let left = GeometricLift::recover(model, curve, left)?;
    let right = GeometricLift::recover(model, curve, right)?;
    let k = left
        .radicand
        .as_ref()
        .or(right.radicand.as_ref())
        .cloned()
        .unwrap_or_else(|| field(curve).one());
    let twist = Twist::new(model, curve, k)?;
    let left = twist.divisor(model, left)?;
    let right = twist.divisor(model, right)?;
    let sum = left.add(model, &twist.curve, &right)?;
    let sum = project_scaled(model, &twist.curve, curve, &sum, &twist.k)?;
    let opposite = right.negate(model, &twist.curve)?;
    let candidate = left.add(model, &twist.curve, &opposite)?;
    let candidate = project_scaled(model, &twist.curve, curve, &candidate, &twist.k)?;
    if sum.compare(model, difference)?.is_eq() {
        Ok(candidate)
    } else if candidate.compare(model, difference)?.is_eq() {
        Ok(sum)
    } else {
        Err(Error::plain(
            "rvKummerDiffAdd difference is incompatible with the two input points",
        ))
    }
}

pub(crate) fn richelot(
    model: &mut Model,
    correspondence: &Correspondence,
    point: &Point,
) -> Result<Point> {
    let source = &correspondence.source;
    let lift = GeometricLift::recover(model, source, point)?;
    let k = lift.preferred_radicand(source);
    if k == field(source).one() {
        // Both the original curve and its selected target already have the
        // required coordinates when the geometric ordinate is over K.
        let divisor = Divisor::new(model, source, lift.u, lift.w)?;
        let image = correspondence.push(model, &divisor)?;
        return Point::project(model, &correspondence.target, &image);
    }
    let twist = Twist::new(model, source, k)?;
    let divisor = twist.divisor(model, lift)?;
    let (normalized, target_scale) =
        correspondence.twisted(model, twist.curve.clone(), &twist.k, &twist.inverse_powers)?;
    let image = normalized.push(model, &divisor)?;
    project_scaled(
        model,
        &normalized.target,
        &correspondence.target,
        &image,
        &target_scale,
    )
}

struct GeometricLift {
    u: Polynomial,
    // A nonnodal ordinate is sqrt(radicand) times this polynomial.
    w: Polynomial,
    radicand: Option<Scalar>,
}

impl GeometricLift {
    fn preferred_radicand(&self, curve: &Curve) -> Scalar {
        self.radicand.clone().unwrap_or_else(|| field(curve).one())
    }

    fn recover(model: &mut Model, curve: &Curve, point: &Point) -> Result<Self> {
        let field = field(curve);
        let (u, w, radicand) = match point.recover(model, curve)? {
            Recovery::Identity => (vec![field.one()], Vec::new(), None),
            Recovery::Linear { x, y_squared } => {
                let u = vec![field.ring.negate(model, &x)?, field.one()];
                if poly::zero(&y_squared) {
                    (u, Vec::new(), None)
                } else {
                    (u, vec![field.one()], Some(y_squared))
                }
            }
            Recovery::Quadratic { u, a, b, c } => {
                if !poly::zero(&a) {
                    let inverse = field.inverse_unit(model, &a)?;
                    let linear = field.multiply(model, &b, &inverse)?;
                    (u, vec![field.one(), linear], Some(a))
                } else if !poly::zero(&c) {
                    (u, vec![field.zero(), field.one()], Some(c))
                } else {
                    (u, Vec::new(), None)
                }
            }
        };
        Ok(Self { u, w, radicand })
    }
}

struct Twist {
    curve: Rc<Curve>,
    k: Scalar,
    inverse_powers: [Scalar; 6],
}

impl Twist {
    fn new(model: &mut Model, curve: &Curve, k: Scalar) -> Result<Self> {
        let field = field(curve);
        let inverse = field.inverse_unit(model, &k)?;
        let mut inverse_powers = std::array::from_fn(|_| field.one());
        for index in 1..inverse_powers.len() {
            inverse_powers[index] = field.multiply(model, &inverse_powers[index - 1], &inverse)?;
        }
        // The twist has coordinates X=x/k and Y=y/(sqrt(k) k²).
        let f = poly::scale_homogeneous(model, field, curve.polynomial(), &inverse_powers, 5)?;
        let normalized = Curve::new(model, curve.grade.ring.clone(), curve.ring.clone(), f)?;
        Ok(Self {
            curve: Rc::new(normalized),
            k,
            inverse_powers,
        })
    }

    fn divisor(&self, model: &mut Model, lift: GeometricLift) -> Result<Divisor> {
        let field = field(&self.curve);
        let u = poly::scale_homogeneous(
            model,
            field,
            &lift.u,
            &self.inverse_powers,
            lift.u.len() - 1,
        )?;
        let mut v = poly::scale_homogeneous(model, field, &lift.w, &self.inverse_powers, 2)?;
        if let Some(radicand) = lift.radicand.filter(|radicand| radicand != &self.k) {
            let ratio = field.multiply(model, &radicand, &self.inverse_powers[1])?;
            // Opposite Frobenius signs make the sum and difference nonrational
            // on the quartic unless one of the inputs is a node.
            let ordinate_scale = square_root(model, field, &ratio)?.ok_or_else(|| {
                Error::plain(
                    "rvKummerDiffAdd nonnodal inputs have different quadratic twist classes",
                )
            })?;
            v = poly::scale(model, field, &v, &ordinate_scale)?;
        }
        Divisor::new(model, &self.curve, u, v)
    }
}

fn project_scaled(
    model: &mut Model,
    source: &Curve,
    target: &Curve,
    divisor: &Divisor,
    scale: &Scalar,
) -> Result<Point> {
    let point = Point::project(model, source, divisor)?;
    let field = field(target);
    let mut coordinates = point.coordinates().clone();
    // The two terms of −q₀−p q₂ have weight three under x = kX.
    let square = field.multiply(model, scale, scale)?;
    let cube = field.multiply(model, &square, scale)?;
    for (coordinate, weight) in coordinates[1..].iter_mut().zip([scale, &square, &cube]) {
        *coordinate = field.multiply(model, coordinate, weight)?;
    }
    Point::new(model, target, coordinates)
}
