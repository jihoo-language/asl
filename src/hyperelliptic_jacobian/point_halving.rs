use super::{field, Divisor};
use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::model::Model;

pub(crate) fn witness(
    model: &mut Model,
    curve: &Curve,
    coordinates: Vec<Vec<i128>>,
) -> Result<Polynomial> {
    let field = field(curve);
    field.check_polynomial(coordinates.len())?;
    let root = coordinates
        .into_iter()
        .map(|coefficient| field.ring.literal(model, coefficient))
        .collect::<Result<_>>()?;
    let root = poly::trim(root);
    model.allocate(root.len() * field.cells())?;
    Ok(root)
}

pub(crate) fn half(
    model: &mut Model,
    curve: &Curve,
    x: &[u64],
    y: &[u64],
    root: &[Scalar],
) -> Result<Divisor> {
    let field = field(curve);
    let f = curve.polynomial();
    if field.ring.multiply(model, y, y)? != poly::evaluate(model, field, f, x)? {
        return Err(Error::plain(
            "rvJacHalfPoint requires an affine point on the curve",
        ));
    }
    let reflection = [x.to_vec(), field.signed(-1)];
    let expected_square = poly::remainder(model, field, reflection.to_vec(), f)?;
    model.spend(root.len() * field.cells())?;
    let root = poly::remainder(model, field, root.to_vec(), f)?;
    if poly::multiply_mod(model, field, &root, &root, f)? != expected_square {
        return Err(Error::plain(
            "rvJacHalfPoint witness must square to x−X modulo f",
        ));
    }

    // The common roots are the chosen values of the square root at the branch
    // points. At a ramified input, z = 0 is simple in z − r(x − z²).
    let branch = reflect(model, field, f, &reflection)?;
    let branch = at_square(model, field, branch)?;
    let lifted_root = reflect(model, field, &root, &reflection)?;
    let lifted_root = at_square(model, field, lifted_root)?;
    let z = [field.zero(), field.one()];
    let relation = poly::sum(model, field, &z, &lifted_root, true)?;
    let selected = poly::gcd(model, field, branch, relation)?;
    if selected.len() != f.len() {
        return Err(Error::plain(
            "rvJacHalfPoint witness gives the wrong degree",
        ));
    }
    let selected = monic(model, field, selected)?;
    if selected[0].as_slice() != y {
        return Err(Error::plain("rvJacHalfPoint witness norm must equal −y"));
    }

    // R(z) = E(z²) + z O(z²) implies f(X) − E(x − X)² = (X − x) O(x − X)².
    // The divisor of y − E(x − X) identifies the half with the conjugate
    // of the points cut out by O(x − X) and y = E(x − X).
    model.spend(selected.len())?;
    let mut even = Vec::with_capacity(selected.len() / 2);
    let mut odd = Vec::with_capacity(selected.len() / 2);
    for (index, coefficient) in selected.into_iter().enumerate() {
        if index.is_multiple_of(2) {
            even.push(coefficient);
        } else {
            odd.push(coefficient);
        }
    }
    let u = reflect(model, field, &odd, &reflection)?;
    let u = monic(model, field, u)?;
    let even = reflect(model, field, &even, &reflection)?;
    let v = poly::scale(model, field, &even, &field.signed(-1))?;
    let half = Divisor::new(model, curve, u, v)?;
    let point = Divisor::from_point(model, curve, x, y.to_vec())?;
    let double = half.add(model, curve, &half)?;
    if double != point {
        return Err(Error::plain(
            "rvJacHalfPoint failed Cantor doubling verification",
        ));
    }
    Ok(half)
}

fn reflect(
    model: &mut Model,
    field: Field<'_>,
    polynomial: &[Scalar],
    reflection: &[Scalar],
) -> Result<Polynomial> {
    let mut result = Vec::new();
    for coefficient in polynomial.iter().rev() {
        let product = poly::multiply(model, field, &result, reflection)?;
        result = poly::sum(
            model,
            field,
            &product,
            std::slice::from_ref(coefficient),
            false,
        )?;
    }
    Ok(result)
}

fn at_square(model: &mut Model, field: Field<'_>, polynomial: Polynomial) -> Result<Polynomial> {
    let length = polynomial.len().saturating_mul(2).saturating_sub(1);
    field.check_polynomial(length)?;
    model.spend(length * field.cells())?;
    let mut result = Vec::with_capacity(length);
    let mut coefficients = polynomial.into_iter();
    if let Some(constant) = coefficients.next() {
        result.push(constant);
    }
    for coefficient in coefficients {
        result.push(field.zero());
        result.push(coefficient);
    }
    Ok(result)
}

fn monic(model: &mut Model, field: Field<'_>, polynomial: Polynomial) -> Result<Polynomial> {
    let leading = polynomial
        .last()
        .ok_or_else(|| Error::plain("rvJacHalfPoint obtained a zero polynomial"))?;
    if leading == &field.one() {
        return Ok(polynomial);
    }
    let inverse = field.inverse_unit(model, leading)?;
    poly::scale(model, field, &polynomial, &inverse)
}
