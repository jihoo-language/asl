use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::model::Model;
use crate::rv_witt_spectrum::{self as witt, HilbertBridge, Spectrum};
use crate::syntax::Label;

use super::super::hilbert::{self, Integer, Magnitude, SurfaceKind};
use super::{semilinear, Word};

struct HilbertPowers {
    length: u128,
    fields: Vec<Magnitude>,
    middle: Vec<Integer>,
}

impl HilbertPowers {
    fn new(
        model: &mut Model,
        word: &Word,
        curve: &Curve,
        base: &Label,
        (length, r, steps): (u128, u128, u128),
    ) -> Result<Self> {
        let q = curve.ring.size_of_residue_field()?;
        let through = length
            .checked_mul(steps)
            .ok_or_else(hilbert::intermediate_overflow)?;
        let last = r
            .checked_mul(through)
            .ok_or_else(hilbert::intermediate_overflow)?;
        // Bound the series allocation before materializing any of its terms.
        hilbert::power(model, q, last)?;
        let limit = usize::try_from(through).map_err(|_| hilbert::intermediate_overflow())?;
        let fields = hilbert::field_powers(model, q, r, limit)?;
        let frobenius = semilinear::node_permutation(model, word, curve, base)?;
        let middle = semilinear::middle_traces(
            model,
            word,
            curve,
            &frobenius,
            r,
            &fields,
            hilbert::intermediate_overflow,
        )?;
        Ok(Self {
            length,
            fields,
            middle,
        })
    }

    fn fixed(&self, model: &mut Model, power: usize) -> Result<Integer> {
        let value = hilbert::k3_lefschetz_coefficient(
            model,
            self.length,
            &self.fields,
            &self.middle,
            power,
        )?;
        if value.negative && !value.is_zero() {
            return Err(Error::plain(
                "K3 Frobenius Hilbert-scheme fixed-point count became negative",
            ));
        }
        Ok(value)
    }
}

pub(super) fn report(
    model: &mut Model,
    name: &Label,
    base: &Label,
    length: u128,
    period: u128,
) -> Result<String> {
    let (word, curve) = semilinear::word_and_base(model, name, base)?;
    word.basis
        .require_base_extension(model, &curve)
        .map_err(|error| error.at(base.span))?;
    if period == 0 {
        return Err(Error::plain("primitive period must be positive"));
    }
    let count = if length == 0 {
        u128::from(period == 1)
    } else {
        let powers = HilbertPowers::new(model, &word, &curve, base, (length, 1, period))?;
        let value = witt::primitive_coordinate(model, period, |model, k| {
            let power = usize::try_from(k).map_err(|_| hilbert::intermediate_overflow())?;
            powers.fixed(model, power)
        })?;
        if value.negative && !value.is_zero() {
            return Err(Error::plain(
                "K3 Frobenius Hilbert-scheme primitive orbit count became negative",
            ));
        }
        value.magnitude.small().ok_or_else(|| {
            Error::plain(
                "K3 Frobenius Hilbert-scheme primitive orbit count exceeds unsigned 128-bit range",
            )
        })?
    };
    Ok(format!(
        "rvK3FrobHilbOrbits[{}, {}, {length}, {period}] = {count}",
        name.text, base.text
    ))
}

pub(crate) fn hilbert_spectrum(model: &mut Model, query: &HilbertBridge) -> Result<Spectrum> {
    let (word, curve) = semilinear::word_and_base(model, &query.word, &query.base)?;
    word.basis
        .require_base_extension(model, &curve)
        .map_err(|error| error.at(query.base.span))?;
    if query.r == 0 {
        return Err(Error::plain(
            "K3 Frobenius iteration count must be positive",
        ));
    }
    if query.length == 0 {
        return Spectrum::line(1);
    }
    // The unit and the top-degree class in mode one already give n + 1
    // independent classes in length n.
    let length = usize::try_from(query.length).unwrap_or(usize::MAX);
    witt::require_degree(length.saturating_add(1))?;
    let rows = hilbert::coefficients(model, query.length, None, None, SurfaceKind::K3)?;
    let rank = hilbert::unsigned(rows[length][0])
        .and_then(|rank| usize::try_from(rank).ok())
        .unwrap_or(usize::MAX);
    witt::require_degree(rank)?;
    let powers = HilbertPowers::new(
        model,
        &word,
        &curve,
        &query.base,
        (query.length, query.r, rank as u128),
    )?;
    let mut ghosts = Vec::with_capacity(rank);
    for j in 1..=rank {
        ghosts.push(powers.fixed(model, j)?);
    }
    Spectrum::effective(witt::determinant_from_ghosts(model, &ghosts)?)
}
