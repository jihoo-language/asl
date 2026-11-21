use crate::cobar::{self, Word};
use crate::cobar_squares::{cup_support, frobenius};
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::profile::Profile;

pub(crate) fn bidegree(k: usize, s: usize, t: usize) -> Result<(usize, usize)> {
    let length = s
        .checked_add(k)
        .ok_or_else(|| Error::plain("cNakamura cobar length is out of range"))?;
    let degree = t
        .checked_mul(2)
        .ok_or_else(|| Error::plain("cNakamura internal degree is out of range"))?;
    Ok((length, degree))
}

pub(crate) fn apply(
    model: &mut Model,
    profile: Profile,
    k: usize,
    s: usize,
    words: &[Word],
) -> Result<Value> {
    if words.is_empty() || k.saturating_sub(s) > 1 {
        return Ok(cobar::zero());
    }

    let result = if k == 0 {
        frobenius(model, profile, words)?
    } else if let Some(i) = s.checked_sub(k) {
        let square = cup_support(model, profile, i, words, words)?;
        let differential = profile.differential_support(model, words)?;
        let correction = cup_support(model, profile, i + 1, words, &differential)?;
        cobar::add_support(model, &square, &correction)?
    } else {
        // At k=s+1 only the correction remains, with cup index zero.
        let differential = profile.differential_support(model, words)?;
        cobar::concat_temporary(model, words, &differential)?
    };
    model.allocate(result.iter().map(cobar::word_cells).sum())?;
    Ok(Value::Cobar(result.into()))
}
