use crate::cobar::{self, Monomial, Word};
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::FactorLiteral;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) s: usize,
    pub(crate) t: i128,
    pub(crate) w: i64,
}

impl Grade {
    pub(crate) fn differential(self) -> Result<Self> {
        let s = self
            .s
            .checked_add(1)
            .ok_or_else(|| Error::plain("motivic cobar length is out of range"))?;
        Ok(Self { s, ..self })
    }

    pub(crate) fn previous(self) -> Option<Self> {
        Some(Self {
            s: self.s.checked_sub(1)?,
            ..self
        })
    }

    pub(crate) fn product(self, right: Self) -> Result<Self> {
        let s = self
            .s
            .checked_add(right.s)
            .ok_or_else(|| Error::plain("motivic cobar length is out of range"))?;
        let t = self
            .t
            .checked_add(right.t)
            .ok_or_else(|| Error::plain("motivic internal degree is out of range"))?;
        let w = self.w.checked_add(right.w).ok_or_else(weight_overflow)?;
        Ok(Self { s, t, w })
    }

    pub(crate) fn shift(self, power: usize) -> Result<Self> {
        let w = i64::try_from(i128::from(self.w) - power as i128).map_err(|_| weight_overflow())?;
        Ok(Self { w, ..self })
    }

    pub(crate) fn classical_degree(self) -> Result<usize> {
        usize::try_from(self.t).map_err(|_| {
            Error::plain(
                "realize target is outside the supported classical topological degree range",
            )
        })
    }

    pub(crate) fn tau_power(self, word: &Word) -> u128 {
        let intrinsic = (self.t - odd_exponents(word) as i128) / 2;
        (intrinsic - i128::from(self.w)) as u128
    }
}

fn weight_overflow() -> Error {
    Error::plain("motivic weight is out of range")
}

pub(crate) fn intrinsic_weight(word: &Word, t: usize) -> usize {
    (t - odd_exponents(word)) / 2
}

pub(crate) fn odd_exponents(word: &Word) -> usize {
    word.iter().map(|factor| cobar::odd_exponents(factor)).sum()
}

pub(crate) fn literal(
    model: &mut Model,
    coefficient: usize,
    factors: Vec<FactorLiteral>,
) -> Result<(Grade, Value)> {
    let mut carry = coefficient as i128;
    let factors = factors
        .into_iter()
        .map(|(xi, tau, span)| {
            let width = xi.len().max(tau.len());
            model.spend(1 + width)?;
            let mut exponents = Vec::with_capacity(width);
            for i in 0..width {
                let x = xi.get(i).copied().unwrap_or_default();
                let a = tau.get(i).copied().unwrap_or_default();
                carry += (a / 2) as i128;
                let encoded = x
                    .checked_mul(2)
                    .and_then(|value| value.checked_add(a))
                    .ok_or_else(|| Error::new(span, "motivic monomial exponent is out of range"))?;
                exponents.push(encoded);
            }
            Ok((exponents, span))
        })
        .collect::<Result<Vec<_>>>()?;
    let ((s, t), value) = cobar::literal(model, factors)?;
    let Value::Cobar(words) = &value else {
        unreachable!();
    };
    let w = i64::try_from(intrinsic_weight(&words[0], t) as i128 - carry)
        .map_err(|_| weight_overflow())?;
    Ok((Grade { s, t: t as i128, w }, value))
}

pub(crate) fn basis(model: &mut Model, grade: Grade, limit: usize) -> Result<Vec<Word>> {
    let Grade { s, t, w } = grade;
    if t < 0 {
        return Ok(Vec::new());
    }
    let t = usize::try_from(t).map_err(|_| {
        Error::plain("motivic topological basis degree exceeds the supported range")
    })?;
    if w <= 0 {
        return cobar::basis(model, s, t, limit);
    }
    let minimum = 2 * i128::from(w);
    if minimum > t as i128 {
        return Ok(Vec::new());
    }
    let max_odd = (t as i128 - minimum) as usize;
    let Some(excess) = t.checked_sub(s) else {
        return Ok(Vec::new());
    };
    if excess < s.saturating_sub(max_odd) {
        return Ok(Vec::new());
    }
    if max_odd != 0 {
        return cobar::basis_with_odd_bound(model, s, t, max_odd, limit);
    }
    // At maximal weight all exponents are even. Halving them identifies
    // the whole word basis with the classical basis in degree t/2.
    let mut words = cobar::basis(model, s, t / 2, limit)?;
    for word in &mut words {
        model.spend(cobar::word_cells(word))?;
        *word = word
            .iter()
            .map(|factor| factor.iter().map(|power| 2 * power).collect::<Monomial>())
            .collect();
    }
    Ok(words)
}
