use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::syntax::Label;
use crate::wide_integer::Natural;

use super::super::{k3, nikulin::NodePermutation};
use super::{spectrum, Word};

pub(super) fn fixed(model: &mut Model, name: &Label, base: &Label, power: u128) -> Result<String> {
    let (word, curve, frobenius) = data(model, name, base)?;
    let q = curve.ring.size_of_residue_field()?;
    let exponent = u32::try_from(power).map_err(|_| count_overflow())?;
    // The marked isometry preserves the positive cone. Its eigenvalues off
    // the unit circle, if present, are a positive reciprocal pair. Thus its
    // trace is at least -17, and the primitive trace has magnitude at most
    // 5Q by the Weil bounds. If Q does not fit, the count cannot fit either.
    let field_order = u128::from(q)
        .checked_pow(exponent)
        .ok_or_else(count_overflow)?;
    let field_order = Wide {
        magnitude: Natural::from_small(field_order),
        negative: false,
    };
    let matrix = word.marked_matrix(model, Some(&frobenius))?;
    let marked_trace = spectrum::power_trace(model, matrix, exponent)?;
    let degree = exponent as usize;
    let weil = curve.weil(model)?;
    let traces = k3::h1_power_sums(model, weil, 2 * degree, count_overflow)?;
    let exterior_trace = k3::exterior_power_sum(model, &traces, degree, count_overflow)?;
    let primitive_trace = exterior_trace
        .add(field_order.negate())
        .map_err(|_| count_overflow())?;
    let primitive_trace = if word.reverses_primitive() && power % 2 != 0 {
        primitive_trace.negate()
    } else {
        primitive_trace
    };

    // Combine signed terms before testing the unsigned result. Near the
    // output bound the intermediate Q^2 and Q tr(A^r) may require wide limbs.
    let count = field_order
        .add(marked_trace)
        .and_then(|sum| field_order.multiply_wide(sum))
        .and_then(|sum| sum.add(primitive_trace))
        .and_then(|sum| sum.add(Wide::from_i128(1)))
        .map_err(|_| count_overflow())?;
    if count.negative && !count.is_zero() {
        return Err(Error::plain(
            "K3 Frobenius fixed-point count became negative",
        ));
    }
    let count = count.magnitude.small().ok_or_else(count_overflow)?;
    Ok(format!(
        "rvK3FrobFixed[{}, {}, {power}] = {count}",
        name.text, base.text
    ))
}

pub(super) fn zeta(model: &mut Model, name: &Label, base: &Label) -> Result<String> {
    let (word, curve, frobenius) = data(model, name, base)?;
    let q = curve.ring.size_of_residue_field()?;
    let q2 = u128::from(q) * u128::from(q);
    let cross = k3::cross_quartic(q, curve.weil(model)?)?;
    let matrix = word.marked_matrix(model, Some(&frobenius))?;
    // On divisor cohomology the pullback is q times the inverse of this
    // push action. The preserved intersection form gives the same spectrum.
    let marked = spectrum::characteristic_coefficients(model, matrix)?;
    let (primitive_sign, cross_argument) = if word.reverses_primitive() {
        ('+', "-T")
    } else {
        ('-', "T")
    };
    Ok(format!(
        "rvK3FrobZeta[{}, {}] = 1/((1-T)(1-{q2}*T)D({q}*T)(1{primitive_sign}{q}*T)h({cross_argument})), D={marked:?}, h={cross:?}",
        name.text, base.text,
    ))
}

fn data(
    model: &mut Model,
    name: &Label,
    base: &Label,
) -> Result<(Rc<Word>, Rc<Curve>, NodePermutation)> {
    let word = model.k3_automorphism(name)?;
    let curve = model
        .kummer_curve(&base.text)
        .map_err(|error| error.at(base.span))?;
    let frobenius =
        NodePermutation::new(model, &word.basis, &curve, 1).map_err(|error| error.at(base.span))?;
    Ok((word, curve, frobenius))
}

fn count_overflow() -> Error {
    Error::plain("K3 Frobenius fixed-point count exceeds unsigned 128-bit range")
}
