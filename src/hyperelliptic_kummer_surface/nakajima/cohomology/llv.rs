use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;

use super::super::{grade_range, Factor, Grade, Mode, POINT};
use super::lefschetz::{divisor_pair, Operation};
use super::{
    add_scaled, cancel, divide, factorial, normalize, positive, pure, Accumulator, Class,
    Denominator, Integer, State,
};

pub(crate) fn beauville_report(
    model: &mut Model,
    a: &Class,
    b: &Class,
    length: i128,
) -> Result<String> {
    let ratio = beauville(model, a, b, length)?;
    let numerator = ratio.value.signed().ok_or_else(arithmetic)?;
    let denominator = ratio.divisor.small().ok_or_else(arithmetic)?;
    if denominator == 1 {
        Ok(numerator.to_string())
    } else {
        Ok(format!("{numerator}/{denominator}"))
    }
}

fn beauville(model: &mut Model, a: &Class, b: &Class, length: i128) -> Result<Ratio> {
    if is_zero(a) || is_zero(b) || length < 1 {
        return Ok(Ratio::zero());
    }
    let raw = divisor_pair(model, a, b, length)?;
    if raw.is_zero() {
        return Ok(Ratio::zero());
    }
    // For n ≥ 2 the compressed raw pairing omits the common factor
    // (n−1)(n−2)!². For n=1 it already is the surface Gram pairing.
    let n = length as u128;
    let fact = factorial(model, n.saturating_sub(2))?;
    let scale = Integer {
        magnitude: fact,
        negative: false,
    };
    let common = positive((n - 1).max(1));
    Ratio::new(
        model,
        [raw, common, scale, scale],
        [small(a.denominator), small(b.denominator)],
    )
}

pub(crate) fn lower(model: &mut Model, b: &Class, x: &Class, grade: &Grade) -> Result<Class> {
    if !grade.supports_states() || grade.degree == 0 || is_zero(b) || is_zero(x) {
        return Ok(Class::zero());
    }
    let n = grade.length as u128;
    if grade.degree == 2 {
        // F_b(x) = (2n)q(b,x)·1. Cancel n! symbolically in raw coordinates.
        let pair = divisor_pair(model, b, x, grade.length)?;
        if pair.is_zero() {
            return Ok(Class::zero());
        }
        let fact = Integer {
            magnitude: factorial(model, n.saturating_sub(2))?,
            negative: false,
        };
        let coefficient = Ratio::new(
            model,
            [pair, positive(2), fact],
            [small(b.denominator), small(x.denominator)],
        )?;
        return sum(model, [(&pure(n, 0), coefficient)]);
    }
    if !divisor_pair(model, b, b, grade.length)?.is_zero() {
        let (lowered, coefficient) = component(model, b, x, grade)?;
        return sum(model, [(&lowered.numerator, coefficient)]);
    }
    let (anchor, translated) = nonisotropic_pair(model, b, grade.length)?;
    let (left, a) = component(model, &translated, x, grade)?;
    let (right, t) = component(model, &anchor, x, grade)?;
    sum(
        model,
        [(&left.numerator, a), (&right.numerator, t.negate())],
    )
}

fn component(model: &mut Model, b: &Class, x: &Class, grade: &Grade) -> Result<(Class, Ratio)> {
    let lowered = Operation::Lower.apply(model, b, x, grade)?;
    if is_zero(&lowered) {
        return Ok((lowered, Ratio::zero()));
    }
    let norm = beauville(model, b, b, grade.length)?;
    let coefficient = norm.adjust(model, positive(1), lowered.denominator, 1)?;
    Ok((lowered, coefficient))
}

fn nonisotropic_pair(model: &mut Model, b: &Class, length: i128) -> Result<(Class, Class)> {
    let gram = model.nakajima_gram()?;
    let color = (1..POINT)
        .find(|&color| gram.entry(color, color) != 0)
        .ok_or_else(|| Error::plain("no nonisotropic surface anchor for Hilbert LLV lowering"))?;
    let mut monomial = Vec::new();
    if length > 1 {
        monomial.push(Factor {
            mode: Mode {
                length: 1,
                color: 0,
            },
            exponent: length as u128 - 1,
        });
    }
    monomial.push(Factor {
        mode: Mode { length: 1, color },
        exponent: 1,
    });
    // Use the same denominator and add in the direction of zero at this axis.
    let sign = if b
        .numerator
        .terms
        .get(&monomial)
        .is_some_and(|&value| value > 0)
    {
        -1
    } else {
        1
    };
    let axis = Class {
        numerator: State {
            terms: BTreeMap::from([(monomial, sign)]),
        },
        denominator: b.denominator,
    };
    let translated = b.add(model, &axis)?;
    if !divisor_pair(model, &translated, &translated, length)?.is_zero() {
        return Ok((axis, translated));
    }
    // q(b)=0 and q(axis)≠0, so q(b+c·axis) has at most one nonzero root.
    let twice = axis.scale(model, 2, 1)?;
    let translated = b.add(model, &twice)?;
    if divisor_pair(model, &translated, &translated, length)?.is_zero() {
        return Err(Error::plain(
            "failed to certify a nonisotropic Hilbert LLV decomposition",
        ));
    }
    Ok((twice, translated))
}

pub(crate) fn rotate(
    model: &mut Model,
    a: &Class,
    b: &Class,
    x: &Class,
    grade: &Grade,
) -> Result<Class> {
    if !grade.supports_states()
        || grade.degree == 0
        || grade.length.checked_mul(4) == Some(grade.degree)
        || is_zero(a)
        || is_zero(b)
        || is_zero(x)
    {
        return Ok(Class::zero());
    }
    let n = grade.length;
    if grade.degree == 2 {
        let bx = beauville(model, b, x, n)?.adjust(model, positive(1), a.denominator, 1)?;
        let ax = beauville(model, a, x, n)?.adjust(model, positive(1), b.denominator, 1)?;
        return sum(model, [(&a.numerator, bx), (&b.numerator, ax.negate())]);
    }
    let lowered = lower(model, b, x, grade)?;
    let first = a.cup(model, &lowered, n)?;
    let raised = a.cup(model, x, n)?;
    let up = Grade {
        degree: grade.degree.checked_add(2).ok_or_else(grade_range)?,
        ..grade.clone()
    };
    let second = lower(model, b, &raised, &up)?;
    let grading = Integer::from_i128(grade.degree).add(Integer::product(-2, n))?;
    let correction = if grading.is_zero() {
        Ratio::zero()
    } else {
        beauville(model, a, b, n)?.adjust(model, grading, x.denominator, 2)?
    };
    let left = Ratio::new(model, [positive(1)], [small(2), small(first.denominator)])?;
    let right = Ratio::new(
        model,
        [positive(1).negate()],
        [small(2), small(second.denominator)],
    )?;
    sum(
        model,
        [
            (&first.numerator, left),
            (&second.numerator, right),
            (&x.numerator, correction.negate()),
        ],
    )
}

#[derive(Clone, Copy)]
struct Ratio {
    value: Integer,
    divisor: Denominator,
}

impl Ratio {
    fn zero() -> Self {
        Self {
            value: Integer::from_i128(0),
            divisor: small(1),
        }
    }

    fn new<const A: usize, const B: usize>(
        model: &mut Model,
        mut numerator: [Integer; A],
        mut denominator: [Denominator; B],
    ) -> Result<Self> {
        if numerator.iter().any(|value| value.is_zero()) {
            return Ok(Self::zero());
        }
        for value in &mut numerator {
            for divisor in &mut denominator {
                cancel(model, &mut value.magnitude, divisor)?;
            }
        }
        let value = numerator.into_iter().try_fold(positive(1), |left, right| {
            left.multiply_wide(right).map_err(|_| arithmetic())
        })?;
        let divisor = denominator.into_iter().try_fold(small(1), |left, right| {
            left.checked_product(right).ok_or_else(arithmetic)
        })?;
        Ok(Self { value, divisor })
    }

    fn adjust(self, model: &mut Model, value: Integer, divisor: u128, extra: u128) -> Result<Self> {
        Self::new(
            model,
            [self.value, value],
            [self.divisor, small(divisor), small(extra)],
        )
    }

    fn negate(mut self) -> Self {
        self.value = self.value.negate();
        self
    }
}

fn sum<const N: usize>(model: &mut Model, terms: [(&State, Ratio); N]) -> Result<Class> {
    let mut denominator = small(1);
    let used =
        |(state, factor): &(&State, Ratio)| !state.terms.is_empty() && !factor.value.is_zero();
    for (_, factor) in terms.iter().filter(|entry| used(entry)) {
        let overlap = denominator.gcd(model, factor.divisor)?;
        denominator = divide(model, denominator, overlap)?
            .checked_product(factor.divisor)
            .ok_or_else(arithmetic)?;
    }
    let mut output = Accumulator::default();
    for (state, factor) in terms.into_iter().filter(used) {
        let quotient = divide(model, denominator, factor.divisor)?;
        let scale = factor
            .value
            .multiply_wide(Integer {
                magnitude: quotient,
                negative: false,
            })
            .map_err(|_| arithmetic())?;
        add_scaled(model, &mut output, state, scale)?;
    }
    normalize(model, output, denominator)
}

fn small(value: u128) -> Denominator {
    Denominator::from_small(value)
}

fn is_zero(class: &Class) -> bool {
    class.numerator.terms.is_empty()
}

fn arithmetic() -> Error {
    Error::plain("Hilbert Beauville or LLV coefficient exceeds the exact arithmetic range")
}
