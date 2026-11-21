use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::q_koszul::{self, Powers, Term};

pub(crate) type Grade = q_koszul::Grade<i64>;

pub(crate) fn top(selected: &[usize]) -> Result<usize> {
    selected
        .last()
        .copied()
        .ok_or_else(|| Error::plain("localization requires a nonempty Milnor primitive selection"))
}

pub(crate) fn localization_type(model: &mut Model, ty: TypeId) -> Result<TypeId> {
    let (grade, class) = match model.kind(ty) {
        Kind::QKoszul(grade) => (grade, false),
        Kind::QExt(grade) => (grade, true),
        _ => {
            return Err(Error::plain(
                "vLocal requires a QKoszul cochain or QExt class",
            ))
        }
    };
    let s = i64::try_from(grade.s).map_err(|_| degree_overflow("vLocal"))?;
    let grade = Grade {
        source: grade.source,
        target: grade.target,
        selected: grade.selected.clone(),
        s,
        t: grade.t,
    };
    if class {
        model.v_ext_type(grade)
    } else {
        model.v_koszul_type(grade)
    }
}

pub(crate) fn localize(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    value: Value,
) -> Result<Value> {
    match value {
        Value::QKoszul(_) => Ok(value),
        Value::Class(class) if class.is_empty() => Ok(Value::Class([].into())),
        Value::Class(class) => {
            let terms = model.q_ext(source)?.representative(model, &class)?;
            q_koszul::class(model, target, &terms)
        }
        _ => unreachable!(),
    }
}

pub(crate) fn shift_type(model: &mut Model, ty: TypeId, k: i64) -> Result<TypeId> {
    let (grade, class) = match model.kind(ty) {
        Kind::VKoszul(grade) => (grade, false),
        Kind::VExt(grade) => (grade, true),
        _ => {
            return Err(Error::plain(
                "vShift requires a VKoszul cochain or VExt class",
            ))
        }
    };
    let grade = shifted(grade, k)?;
    if class {
        model.v_ext_type(grade)
    } else {
        model.v_koszul_type(grade)
    }
}

pub(crate) fn shifted(grade: &Grade, k: i64) -> Result<Grade> {
    let weight = q_koszul::weight(top(&grade.selected)?) as i128;
    let s = grade
        .s
        .checked_add(k)
        .ok_or_else(|| degree_overflow("vShift"))?;
    let t = i128::from(k)
        .checked_mul(weight)
        .and_then(|delta| i128::from(grade.t).checked_add(delta))
        .and_then(|sum| i64::try_from(sum).ok())
        .ok_or_else(|| degree_overflow("vShift"))?;
    Ok(Grade {
        s,
        t,
        ..grade.clone()
    })
}

pub(crate) fn composed(outer: &Grade, inner: &Grade) -> Result<Grade> {
    q_koszul::matching_endpoints(outer, inner)?;
    let s = outer
        .s
        .checked_add(inner.s)
        .ok_or_else(|| degree_overflow("ecomp"))?;
    let t = outer
        .t
        .checked_add(inner.t)
        .ok_or_else(|| degree_overflow("ecomp"))?;
    Ok(Grade {
        source: inner.source,
        target: outer.target,
        selected: outer.selected.clone(),
        s,
        t,
    })
}

pub(crate) fn shift(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    k: i64,
    value: Value,
) -> Result<Value> {
    if k == 0 {
        return Ok(value);
    }
    multiply(model, source, target, value, |model, terms| {
        shift_terms(model, terms, k)
    })
}

pub(crate) fn multiply(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    value: Value,
    transform: impl FnOnce(&mut Model, &[Term]) -> Result<Vec<Term>>,
) -> Result<Value> {
    match value {
        Value::QKoszul(terms) => {
            let shifted = transform(model, &terms)?;
            Ok(Value::QKoszul(q_koszul::persist(model, shifted)?))
        }
        Value::Class(class) if class.is_empty() => Ok(Value::Class([].into())),
        Value::Class(class) => {
            let terms = model.q_ext(source)?.representative(model, &class)?;
            let shifted = transform(model, &terms)?;
            q_koszul::class(model, target, &shifted)
        }
        _ => unreachable!(),
    }
}

pub(crate) fn shift_terms(model: &mut Model, terms: &[Term], k: i64) -> Result<Vec<Term>> {
    model.spend(terms.iter().map(Term::cells).sum())?;
    let mut result = Vec::with_capacity(terms.len());
    for term in terms {
        let mut powers = term.powers.to_vec();
        let top = powers.last_mut().unwrap();
        *top = top
            .checked_add(i128::from(k))
            .ok_or_else(|| Error::plain("localized Koszul exponent is out of range"))?;
        // Translation by one monomial preserves the canonical term order.
        result.push(Term {
            powers: powers.into(),
            ..term.clone()
        });
    }
    Ok(result)
}

fn degree_overflow(operation: &str) -> Error {
    Error::plain(format!("{operation} signed bidegree is out of range"))
}

pub(crate) fn laurent_basis(
    model: &mut Model,
    weights: &[usize],
    s: i128,
    degree: i128,
    limit: usize,
) -> Result<Vec<Powers>> {
    let (&top, weights) = weights.split_last().unwrap();
    let Some(gap) = degree_gap(s, top, degree)? else {
        return Ok(Vec::new());
    };
    let lower: Vec<_> = weights.iter().map(|&weight| top - weight).collect();
    let mut divisor = 0;
    let divisors = lower
        .iter()
        .map(|&weight| {
            divisor = q_koszul::gcd(divisor, weight);
            divisor
        })
        .collect();
    let mut basis = LaurentBasis {
        s,
        powers: vec![0; weights.len() + 1],
        words: Vec::new(),
        divisors,
        limit,
    };
    basis.enumerate(model, &lower, gap, 0)?;
    Ok(basis.words)
}

fn degree_gap(s: i128, weight: usize, degree: i128) -> Result<Option<u128>> {
    let top = s
        .unsigned_abs()
        .checked_mul(weight as u128)
        .ok_or_else(|| degree_overflow("localized Koszul basis"))?;
    Ok(match (s < 0, degree < 0) {
        (true, true) => degree.unsigned_abs().checked_sub(top),
        (true, false) => None,
        (false, true) => Some(
            top.checked_add(degree.unsigned_abs())
                .ok_or_else(|| degree_overflow("localized Koszul basis"))?,
        ),
        (false, false) => top.checked_sub(degree as u128),
    })
}

struct LaurentBasis {
    s: i128,
    powers: Vec<i128>,
    words: Vec<Powers>,
    divisors: Vec<usize>,
    limit: usize,
}

impl LaurentBasis {
    fn enumerate(
        &mut self,
        model: &mut Model,
        weights: &[usize],
        remaining: u128,
        used: i128,
    ) -> Result<()> {
        model.spend(1)?;
        let Some((&weight, rest)) = weights.split_last() else {
            if remaining == 0 {
                return self.push(model, used);
            }
            return Ok(());
        };
        let common = self.divisors[rest.len()];
        if !remaining.is_multiple_of(common as u128) {
            return Ok(());
        }
        let weight_wide = weight as u128;
        if rest.is_empty() {
            let exponent = i128::try_from(remaining / weight_wide)
                .map_err(|_| degree_overflow("localized Koszul basis"))?;
            self.powers[0] = exponent;
            return self.push(model, add_exponents(used, exponent)?);
        }

        // Skip exponents that cannot leave a multiple of the remaining gcd.
        let step = (self.divisors[rest.len() - 1] / common) as u128;
        let mut exponent = if step == 1 {
            0
        } else {
            let inverse = inverse_mod(model, (weight / common) as u128 % step, step)?;
            ((remaining / common as u128) % step) * inverse % step
        };
        let max = remaining / weight_wide;
        while exponent <= max {
            let signed =
                i128::try_from(exponent).map_err(|_| degree_overflow("localized Koszul basis"))?;
            self.powers[rest.len()] = signed;
            self.enumerate(
                model,
                rest,
                remaining - weight_wide * exponent,
                add_exponents(used, signed)?,
            )?;
            exponent += step;
        }
        Ok(())
    }

    fn push(&mut self, model: &mut Model, used: i128) -> Result<()> {
        if self.words.len() == self.limit {
            return Err(Error::plain(
                "localized Koszul basis exceeds the enumeration limit",
            ));
        }
        let top = self
            .s
            .checked_sub(used)
            .ok_or_else(|| degree_overflow("localized Koszul basis"))?;
        *self.powers.last_mut().unwrap() = top;
        let cells = 1 + q_koszul::power_cells(self.powers.len());
        q_koszul::check_storage(cells * (1 + self.words.len()))?;
        model.spend(cells)?;
        self.words.push(self.powers.as_slice().into());
        Ok(())
    }
}

fn add_exponents(a: i128, b: i128) -> Result<i128> {
    a.checked_add(b)
        .ok_or_else(|| degree_overflow("localized Koszul basis"))
}

fn inverse_mod(model: &mut Model, value: u128, modulus: u128) -> Result<u128> {
    let modulus = modulus as i128;
    let (mut r, mut next_r) = (modulus, value as i128);
    let (mut x, mut next_x) = (0i128, 1i128);
    while next_r != 0 {
        model.spend(1)?;
        let quotient = r / next_r;
        (r, next_r) = (next_r, r - quotient * next_r);
        (x, next_x) = (next_x, x - quotient * next_x);
    }
    debug_assert_eq!(r, 1);
    Ok(x.rem_euclid(modulus) as u128)
}
