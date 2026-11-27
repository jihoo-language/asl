use super::{field, square_root, Point, Recovery};
use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic_jacobian::Divisor;
use crate::model::Model;

pub(crate) fn scale(model: &mut Model, curve: &Curve, n: i128, point: &Point) -> Result<Point> {
    if n == 0 {
        return Ok(Point::identity(curve));
    }
    let lift = GeometricLift::recover(model, curve, point)?;
    let k = lift.radicand.clone().unwrap_or_else(|| field(curve).one());
    let twist = Twist::new(model, curve, k)?;
    let divisor = twist.divisor(model, lift)?;
    let multiple = divisor.scale(model, &twist.curve, n)?;
    twist.project(model, curve, &multiple)
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
    let sum = twist.project(model, curve, &sum)?;
    let opposite = right.negate(model, &twist.curve)?;
    let candidate = left.add(model, &twist.curve, &opposite)?;
    let candidate = twist.project(model, curve, &candidate)?;
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

struct GeometricLift {
    u: Polynomial,
    // A nonnodal ordinate is sqrt(radicand) times this polynomial.
    w: Polynomial,
    radicand: Option<Scalar>,
}

impl GeometricLift {
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
    curve: Curve,
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
        let f = reweight(model, field, curve.polynomial(), &inverse_powers, 5)?;
        let normalized = Curve::new(model, curve.grade.ring.clone(), curve.ring.clone(), f)?;
        Ok(Self {
            curve: normalized,
            k,
            inverse_powers,
        })
    }

    fn divisor(&self, model: &mut Model, lift: GeometricLift) -> Result<Divisor> {
        let field = field(&self.curve);
        let u = reweight(
            model,
            field,
            &lift.u,
            &self.inverse_powers,
            lift.u.len() - 1,
        )?;
        let mut v = reweight(model, field, &lift.w, &self.inverse_powers, 2)?;
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

    fn project(&self, model: &mut Model, curve: &Curve, divisor: &Divisor) -> Result<Point> {
        let point = Point::project(model, &self.curve, divisor)?;
        let field = field(curve);
        let mut coordinates = point.coordinates().clone();
        // The two terms of −q₀−p q₂ have weight three under x = kX.
        let square = field.multiply(model, &self.k, &self.k)?;
        let cube = field.multiply(model, &square, &self.k)?;
        for (coordinate, weight) in coordinates[1..].iter_mut().zip([&self.k, &square, &cube]) {
            *coordinate = field.multiply(model, coordinate, weight)?;
        }
        Point::new(model, curve, coordinates)
    }
}

fn reweight(
    model: &mut Model,
    field: Field<'_>,
    polynomial: &[Scalar],
    inverse_powers: &[Scalar; 6],
    degree: usize,
) -> Result<Polynomial> {
    polynomial
        .iter()
        .enumerate()
        .map(|(index, coefficient)| {
            field.multiply(model, coefficient, &inverse_powers[degree - index])
        })
        .collect()
}
