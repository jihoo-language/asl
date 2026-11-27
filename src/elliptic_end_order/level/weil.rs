use crate::ec_poly::{zero, Field, Scalar};
use crate::elliptic::{Curve, Point};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;

pub(super) fn pairing(
    model: &mut Model,
    curve: &Curve,
    order: u64,
    p: &Point,
    q: &Point,
) -> Result<Scalar> {
    let field = curve.field();
    if p == q || *p == Point::Infinity || *q == Point::Infinity {
        return Ok(field.one());
    }
    let first = miller_unit(model, curve, order, p, q)?;
    let second = miller_unit(model, curve, order, q, p)?;
    let inverse = field.inverse_unit(model, &second)?;
    let value = field.multiply(model, &first, &inverse)?;
    // The two degree-zero divisors have one common point, infinity.
    // Tame reciprocity gives (-1)^N for that shared support.
    if order.is_multiple_of(2) {
        Ok(value)
    } else {
        curve.ring.negate(model, &value)
    }
}

fn miller_unit(
    model: &mut Model,
    curve: &Curve,
    order: u64,
    point: &Point,
    at: &Point,
) -> Result<Scalar> {
    let field = curve.field();
    let mut sum = point.clone();
    let mut result = field.one();
    for bit in (0..order.ilog2()).rev() {
        let (double, unit) = addition_unit(model, curve, &sum, &sum, at)?;
        let square = field.multiply(model, &result, &result)?;
        result = field.multiply(model, &square, &unit)?;
        sum = double;
        if (order >> bit) & 1 != 0 {
            let (next, unit) = addition_unit(model, curve, &sum, point, at)?;
            result = field.multiply(model, &result, &unit)?;
            sum = next;
        }
    }
    if sum != Point::Infinity {
        return Err(Error::plain(
            "elliptic Weil pairing requires points killed by the level",
        ));
    }
    Ok(result)
}

// All line factors have leading unit 1 at infinity for the parameter
// t=x/y. At finite points use t=x-x_Q off the branch points, and t=y
// when y_Q=0. Multiplying leading units also handles shared supports.
fn addition_unit(
    model: &mut Model,
    curve: &Curve,
    a: &Point,
    b: &Point,
    at: &Point,
) -> Result<(Point, Scalar)> {
    let field = curve.field();
    let (sum, slope) = curve.add_with_slope(model, a, b)?;
    let (Point::Affine(xa, ya), Point::Affine(..), Point::Affine(x, y)) = (a, b, at) else {
        return Ok((sum, field.one()));
    };
    let Some(slope) = slope else {
        return Ok((sum, vertical_unit(model, curve, xa, x, y)?));
    };
    let delta = field.subtract(model, x, xa)?;
    let linear = field.multiply(model, &slope, &delta)?;
    let line = field.add(model, ya, &linear)?;
    let value = field.subtract(model, y, &line)?;
    let numerator = if !zero(&value) {
        value
    } else if zero(y) {
        field.one()
    } else {
        line_unit(model, curve, &slope, x, y)?
    };
    let Point::Affine(root, _) = &sum else {
        unreachable!("a nonvertical line has finite sum")
    };
    let vertical = vertical_unit(model, curve, root, x, y)?;
    let inverse = field.inverse_unit(model, &vertical)?;
    let unit = field.multiply(model, &numerator, &inverse)?;
    Ok((sum, unit))
}

fn vertical_unit(
    model: &mut Model,
    curve: &Curve,
    root: &Scalar,
    x: &Scalar,
    y: &Scalar,
) -> Result<Scalar> {
    let field = curve.field();
    if root != x {
        return field.subtract(model, x, root);
    }
    if !zero(y) {
        return Ok(field.one());
    }
    let derivative = derivative(model, curve, x)?;
    field.inverse_unit(model, &derivative)
}

fn derivative(model: &mut Model, curve: &Curve, x: &Scalar) -> Result<Scalar> {
    let field = curve.field();
    let square = field.multiply(model, x, x)?;
    let triple = field.multiply(model, &square, &field.signed(3))?;
    field.add(model, &triple, &curve.a)
}

fn line_unit(
    model: &mut Model,
    curve: &Curve,
    slope: &Scalar,
    x: &Scalar,
    y: &Scalar,
) -> Result<Scalar> {
    let field: Field<'_> = curve.field();
    // y-s(x)=(f(x)-s(x)^2)/(y+s(x)); the numerator is a monic cubic.
    let double_y = field.add(model, y, y)?;
    let linear_s = field.multiply(model, &double_y, slope)?;
    let cubic_derivative = derivative(model, curve, x)?;
    let linear = field.subtract(model, &cubic_derivative, &linear_s)?;
    let first_nonzero = if !zero(&linear) {
        linear
    } else {
        let triple_x = field.multiply(model, x, &field.signed(3))?;
        let square_s = field.multiply(model, slope, slope)?;
        let quadratic = field.subtract(model, &triple_x, &square_s)?;
        if zero(&quadratic) {
            field.one()
        } else {
            quadratic
        }
    };
    let inverse = field.inverse_unit(model, &double_y)?;
    field.multiply(model, &first_nonzero, &inverse)
}
