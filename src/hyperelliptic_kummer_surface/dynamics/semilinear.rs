use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::syntax::Label;
use crate::wide_integer::Natural;

use super::super::{hilbert, k3, nikulin::NodePermutation};
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
    let field_order = Wide::<{ spectrum::LIMBS }> {
        magnitude: Natural::from_small(field_order),
        negative: false,
    };
    let fields = [Natural::from_small(1), field_order.magnitude];
    let middle = middle_traces(
        model,
        &word,
        &curve,
        &frobenius,
        power,
        &fields,
        count_overflow,
    )?;

    // Combine signed terms before testing the unsigned result. Near the
    // output bound the intermediate Q^2 and Q tr(A^r) may require wide limbs.
    let count = field_order
        .multiply_wide(field_order)
        .and_then(|sum| sum.add(middle[1]))
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

pub(super) fn hilbert(
    model: &mut Model,
    name: &Label,
    base: &Label,
    n: u128,
    r: u128,
    degree: Option<u128>,
) -> Result<String> {
    let (word, curve) = word_and_base(model, name, base)?;
    word.basis
        .require_base_extension(model, &curve)
        .map_err(|error| error.at(base.span))?;
    let coefficient = |model: &mut Model, length, degree| {
        hilbert::semilinear_coefficient(model, &curve, length, degree, r, |model, fields| {
            let frobenius = node_permutation(model, &word, &curve, base)?;
            middle_traces(
                model,
                &word,
                &curve,
                &frobenius,
                r,
                fields,
                hilbert::intermediate_overflow,
            )
        })
    };
    if let Some(i) = degree {
        let trace = hilbert::graded_trace(model, &curve, n, r, i, |model, length, degree| {
            coefficient(model, length, Some(degree))
        })?
        .signed()
        .ok_or_else(|| {
            Error::plain("K3 Frobenius Hilbert-scheme trace exceeds signed 128-bit range")
        })?;
        return Ok(format!(
            "rvK3FrobHilbTrace[{}, {}, {n}, {r}, {i}] = {trace}",
            name.text, base.text,
        ));
    }
    let value = if n == 0 {
        Wide::from_i128(1)
    } else {
        coefficient(model, n, None)?
    };
    if value.negative && !value.is_zero() {
        return Err(Error::plain(
            "K3 Frobenius Hilbert-scheme fixed-point count became negative",
        ));
    }
    let count = value.magnitude.small().ok_or_else(|| {
        Error::plain("K3 Frobenius Hilbert-scheme fixed-point count exceeds unsigned 128-bit range")
    })?;
    Ok(format!(
        "rvK3FrobHilbFixed[{}, {}, {n}, {r}] = {count}",
        name.text, base.text,
    ))
}

fn middle_traces<const N: usize>(
    model: &mut Model,
    word: &Word,
    curve: &Curve,
    frobenius: &NodePermutation,
    r: u128,
    fields: &[Natural<N>],
    overflow: fn() -> Error,
) -> Result<Vec<Wide<N>>> {
    let mut middle = vec![Wide::from_i128(0)];
    if fields.len() <= 1 {
        return Ok(middle);
    }
    let step = usize::try_from(r).map_err(|_| overflow())?;
    let max = step.checked_mul(fields.len() - 1).ok_or_else(overflow)?;
    let h1_limit = max.checked_mul(2).ok_or_else(overflow)?;
    // Compute the lattice characteristic polynomial of GP at most once.
    // The sequence supplies every power of the same product, even when
    // Frobenius moves a marked generator in the word.
    let matrix = word.marked_matrix(model, Some(frobenius))?;
    let marked = spectrum::power_sums(model, matrix, max)?;
    let weil = curve.weil(model)?;
    let h1 = k3::h1_power_sums::<{ spectrum::LIMBS }>(model, weil, h1_limit, overflow)?;
    let reversal = word.reverses_primitive();
    for (l, field) in fields.iter().enumerate().skip(1) {
        model.spend(spectrum::LIMBS)?;
        let degree = step * l;
        let q = Wide {
            magnitude: field.resize().ok_or_else(overflow)?,
            negative: false,
        };
        let exterior = k3::exterior_power_sum(model, &h1, degree, overflow)?;
        let primitive = exterior.add(q.negate()).map_err(|_| overflow())?;
        let primitive = if reversal && !degree.is_multiple_of(2) {
            primitive.negate()
        } else {
            primitive
        };
        let value = marked[degree]
            .multiply_wide(q)
            .and_then(|trace| trace.add(primitive))
            .map_err(|_| overflow())?;
        middle.push(Wide {
            magnitude: value.magnitude.resize().ok_or_else(overflow)?,
            negative: value.negative,
        });
    }
    Ok(middle)
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
    let (word, curve) = word_and_base(model, name, base)?;
    let frobenius = node_permutation(model, &word, &curve, base)?;
    Ok((word, curve, frobenius))
}

fn word_and_base(model: &Model, name: &Label, base: &Label) -> Result<(Rc<Word>, Rc<Curve>)> {
    let word = model.k3_automorphism(name)?;
    let curve = model
        .kummer_curve(&base.text)
        .map_err(|error| error.at(base.span))?;
    Ok((word, curve))
}

fn node_permutation(
    model: &mut Model,
    word: &Word,
    curve: &Curve,
    base: &Label,
) -> Result<NodePermutation> {
    NodePermutation::new(model, &word.basis, curve, 1).map_err(|error| error.at(base.span))
}

fn count_overflow() -> Error {
    Error::plain("K3 Frobenius fixed-point count exceeds unsigned 128-bit range")
}
