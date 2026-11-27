use std::cmp::Ordering;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_cyclotomic::{self, Polynomial};
use crate::rv_wzw::{self, Grade, Support};

pub(crate) const MAX_WORD: usize = 256;

#[derive(Clone, Copy)]
pub(crate) enum Step {
    S,
    T(i128),
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Rational {
    pub(crate) numerator: Polynomial,
    pub(crate) denominator: u128,
}

fn arithmetic() -> Error {
    Error::plain("rational cyclotomic arithmetic is out of range")
}

pub(crate) fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn content(terms: &[i128]) -> u128 {
    terms.iter().map(|x| x.unsigned_abs()).fold(0, gcd)
}

fn divide_signed(x: i128, factor: u128) -> i128 {
    debug_assert!(factor > 0);
    if let Ok(factor) = i128::try_from(factor) {
        x / factor
    } else if x == i128::MIN {
        -1
    } else {
        0
    }
}

fn multiply_signed(x: i128, factor: u128) -> Result<i128> {
    if x == 0 {
        return Ok(0);
    }
    if let Ok(factor) = i128::try_from(factor) {
        x.checked_mul(factor).ok_or_else(arithmetic)
    } else if x == -1 && factor == 1u128 << 127 {
        Ok(i128::MIN)
    } else {
        Err(arithmetic())
    }
}

impl Rational {
    pub(crate) fn zero() -> Self {
        Self {
            numerator: [].into(),
            denominator: 1,
        }
    }

    pub(crate) fn integer(value: i128) -> Self {
        if value == 0 {
            Self::zero()
        } else {
            Self {
                numerator: [value].into(),
                denominator: 1,
            }
        }
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.numerator.is_empty()
    }
}

pub(crate) fn rational_value(model: &mut Model, ratio: Rational) -> Result<Value> {
    if !ratio.is_zero() {
        model.allocate(3 + 2 * ratio.numerator.len())?;
    }
    Ok(Value::RVCycloQ(ratio))
}

fn normalized(mut numerator: Vec<i128>, denominator: u128) -> Rational {
    rv_cyclotomic::trim(&mut numerator);
    if numerator.is_empty() {
        Rational::zero()
    } else {
        Rational {
            numerator: numerator.into(),
            denominator,
        }
    }
}

fn persist(model: &mut Model, numerator: Vec<i128>, denominator: u128) -> Result<Value> {
    rational_value(model, normalized(numerator, denominator))
}

pub(crate) fn zero() -> Value {
    Value::RVCycloQ(Rational::zero())
}

pub(crate) fn fraction(
    model: &mut Model,
    numerator: &[i128],
    denominator: u128,
) -> Result<Rational> {
    if denominator == 0 {
        return Err(Error::plain(
            "cyclotomic rational denominator must be positive",
        ));
    }
    model.spend(numerator.len())?;
    let factor = gcd(content(numerator), denominator);
    let numerator = numerator
        .iter()
        .map(|&x| divide_signed(x, factor))
        .collect();
    Ok(normalized(numerator, denominator / factor))
}

pub(crate) fn ratio(model: &mut Model, numerator: &Polynomial, denominator: u128) -> Result<Value> {
    let ratio = fraction(model, numerator, denominator)?;
    rational_value(model, ratio)
}

pub(crate) fn compare(model: &mut Model, a: &Rational, b: &Rational) -> Result<Ordering> {
    if !Rc::ptr_eq(&a.numerator, &b.numerator) {
        model.spend(1 + a.numerator.len().min(b.numerator.len()))?;
    }
    Ok(a.cmp(b))
}

pub(crate) fn negative(model: &mut Model, ratio: &Rational) -> Result<Rational> {
    model.spend(ratio.numerator.len())?;
    let numerator = ratio
        .numerator
        .iter()
        .map(|x| x.checked_neg().ok_or_else(arithmetic))
        .collect::<Result<Vec<_>>>()?;
    Ok(normalized(numerator, ratio.denominator))
}

pub(crate) fn sum(model: &mut Model, a: &Rational, b: &Rational) -> Result<Rational> {
    let length = a.numerator.len().max(b.numerator.len());
    model.spend(length)?;
    let common = gcd(a.denominator, b.denominator);
    let left = a.denominator / common;
    let right = b.denominator / common;
    let mut numerator = vec![0; length];
    for (i, out) in numerator.iter_mut().enumerate() {
        let x = multiply_signed(a.numerator.get(i).copied().unwrap_or(0), right)?;
        let y = multiply_signed(b.numerator.get(i).copied().unwrap_or(0), left)?;
        *out = x.checked_add(y).ok_or_else(arithmetic)?;
    }
    // For reduced a,b, the common denominator can only cancel with
    // gcd(a.denominator,b.denominator). Cancel it before multiplying.
    let factor = gcd(content(&numerator), common);
    let denominator = left
        .checked_mul(b.denominator / factor)
        .ok_or_else(arithmetic)?;
    for entry in &mut numerator {
        *entry = divide_signed(*entry, factor);
    }
    Ok(normalized(numerator, denominator))
}

pub(crate) fn negate(model: &mut Model, a: &Rational) -> Result<Value> {
    let value = negative(model, a)?;
    rational_value(model, value)
}

pub(crate) fn add(model: &mut Model, a: &Rational, b: &Rational) -> Result<Value> {
    let value = sum(model, a, b)?;
    rational_value(model, value)
}

pub(crate) fn multiply(
    model: &mut Model,
    order: usize,
    a: &Rational,
    b: &Rational,
) -> Result<Rational> {
    if a.is_zero() || b.is_zero() {
        return Ok(Rational::zero());
    }
    model.spend(a.numerator.len() + b.numerator.len())?;
    let left_cross = gcd(content(&a.numerator), b.denominator);
    let right_cross = gcd(content(&b.numerator), a.denominator);
    let left: Vec<_> = a
        .numerator
        .iter()
        .map(|&x| divide_signed(x, left_cross))
        .collect();
    let right: Vec<_> = b
        .numerator
        .iter()
        .map(|&x| divide_signed(x, right_cross))
        .collect();
    let mut numerator = rv_cyclotomic::product_coefficients(model, order, &left, &right)?;
    let da = a.denominator / right_cross;
    let db = b.denominator / left_cross;
    let c = content(&numerator);
    let ga = gcd(c, da);
    let gb = gcd(c / ga, db);
    let denominator = (da / ga).checked_mul(db / gb).ok_or_else(arithmetic)?;
    for x in &mut numerator {
        *x = divide_signed(divide_signed(*x, ga), gb);
    }
    Ok(normalized(numerator, denominator))
}

pub(crate) fn shifted(
    model: &mut Model,
    order: usize,
    a: &Rational,
    power: i128,
) -> Result<Rational> {
    let numerator = rv_cyclotomic::shifted(model, order, &a.numerator, power)?;
    Ok(normalized(numerator, a.denominator))
}

pub(crate) fn galois(model: &mut Model, order: usize, a: &Rational, power: i128) -> Result<Value> {
    let power = power.rem_euclid(order as i128) as usize;
    if gcd(power as u128, order as u128) != 1 {
        return Err(Error::plain(
            "cyclotomic Galois power must be coprime to the root order",
        ));
    }
    let numerator = rv_cyclotomic::substitute_root(model, order, &a.numerator, power)?;
    persist(model, numerator, a.denominator)
}

#[derive(Default)]
struct Factors {
    numerator: Vec<u128>,
    denominator: Vec<u128>,
}

fn cancel(model: &mut Model, mut value: u128, opposite: &mut Vec<u128>) -> Result<u128> {
    model.spend(opposite.len())?;
    for factor in opposite.iter_mut() {
        if value == 1 {
            break;
        }
        let common = gcd(value, *factor);
        value /= common;
        *factor /= common;
    }
    opposite.retain(|&factor| factor > 1);
    Ok(value)
}

impl Factors {
    fn denominator(&mut self, model: &mut Model, value: u128) -> Result<()> {
        let value = cancel(model, value, &mut self.numerator)?;
        if value > 1 {
            self.denominator.push(value);
        }
        Ok(())
    }
    fn numerator(&mut self, model: &mut Model, value: u128) -> Result<()> {
        let value = cancel(model, value, &mut self.denominator)?;
        if value > 1 {
            self.numerator.push(value);
        }
        Ok(())
    }
    fn square_multiplier(&mut self, model: &mut Model, value: u128) -> Result<()> {
        if value > 1 {
            self.numerator(model, value)?;
            self.numerator(model, value)?;
        }
        Ok(())
    }
    fn finish(mut self, model: &mut Model, mut polynomial: Vec<i128>) -> Result<Value> {
        let factor = content(&polynomial);
        if factor == 0 {
            return Ok(zero());
        }
        self.numerator(model, factor)?;
        let denominator = self
            .denominator
            .into_iter()
            .try_fold(1u128, |a, b| a.checked_mul(b).ok_or_else(arithmetic))?;
        model.spend(polynomial.len() * (1 + self.numerator.len()))?;
        for x in &mut polynomial {
            *x = divide_signed(*x, factor);
            for &scalar in &self.numerator {
                *x = multiply_signed(*x, scalar)?;
            }
        }
        persist(model, polynomial, denominator)
    }
}

fn primitive_state(model: &mut Model, state: Support, factors: &mut Factors) -> Result<Support> {
    let mut scalar = 0;
    for (_, polynomial) in state.iter() {
        model.spend(polynomial.len())?;
        scalar = gcd(scalar, content(polynomial));
        if scalar == 1 {
            return Ok(state);
        }
    }
    if scalar == 0 {
        return Ok(state);
    }
    factors.square_multiplier(model, scalar)?;
    Ok(state
        .iter()
        .map(|(weight, polynomial)| {
            let normalized = polynomial
                .iter()
                .map(|&x| divide_signed(x, scalar))
                .collect();
            (weight.clone(), normalized)
        })
        .collect())
}

pub(crate) fn word(model: &mut Model, grade: Grade, mut steps: &[Step]) -> Result<Value> {
    if steps.len() > MAX_WORD {
        return Err(Error::plain(
            "WZW gluing word exceeds the limit of 256 letters",
        ));
    }
    model.spend(steps.len())?;
    // T fixes the vacuum. A terminal S² is a scalar times charge,
    // and charge fixes the vacuum. Its squared scalar cancels A².
    while let [prefix @ .., Step::T(_)] | [prefix @ .., Step::S, Step::S] = steps {
        steps = prefix;
    }
    let mut state: Support = [(vec![0i128; grade.fusion.rank].into(), [1].into())].into();
    let mut factors = Factors::default();
    for &letter in steps.iter().rev() {
        state = match letter {
            Step::T(power) => rv_wzw::t_support(model, grade, &state, power)?.into(),
            Step::S => {
                // Factor A=nL^(n-1) without ever constructing A^m.
                factors.denominator(model, grade.fusion.rank as u128)?;
                for _ in 1..grade.fusion.rank {
                    factors.denominator(model, grade.period() as u128)?;
                }
                let raw = rv_wzw::s_support(model, grade, &state, false)?;
                primitive_state(model, raw.into(), &mut factors)?
            }
        };
    }
    let Some((weight, amp)) = state.first() else {
        return Ok(zero());
    };
    if weight[grade.fusion.rank - 1] != 0 {
        return Ok(zero());
    }
    let content = content(amp);
    factors.square_multiplier(model, content)?;
    model.spend(amp.len())?;
    let primitive: Vec<_> = amp.iter().map(|&x| divide_signed(x, content)).collect();
    let conjugate =
        rv_cyclotomic::substitute_root(model, grade.order, &primitive, grade.order - 1)?;
    let squared = rv_cyclotomic::product_coefficients(model, grade.order, &primitive, &conjugate)?;
    factors.finish(model, squared)
}

fn gluing(model: &mut Model, p: usize, q: i128) -> Result<Vec<Step>> {
    if p == 0 {
        return match q {
            1 => Ok(Vec::new()),
            -1 => Ok(vec![Step::S, Step::S]),
            _ => Err(Error::plain("lens p=0 requires q=1 or q=-1")),
        };
    }
    let mut y = p as i128;
    let mut x = q.rem_euclid(y);
    if gcd(x as u128, p as u128) != 1 {
        return Err(Error::plain("lens parameters p and q must be coprime"));
    }
    let mut result = Vec::new();
    while y != 0 {
        model.spend(1)?;
        let mut a = x.div_euclid(y);
        let mut remainder = x.rem_euclid(y);
        if 2 * remainder > y.abs() {
            a += y.signum();
            remainder -= y.abs();
        }
        // This is the inverse of (x,y) -> (a*x-y,x), the action of T^a S.
        (x, y) = (y, -remainder);
        result.extend([Step::T(a), Step::S]);
    }
    if x < 0 {
        result.extend([Step::S, Step::S]);
    }
    Ok(result)
}

pub(crate) fn lens(model: &mut Model, grade: Grade, p: usize, q: i128) -> Result<Value> {
    let steps = gluing(model, p, q)?;
    word(model, grade, &steps)
}
