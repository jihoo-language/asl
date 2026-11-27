use std::cmp::Ordering;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::{Model, Value};
use crate::rv_ore::gcd;
use crate::syntax::{Label, ParsedExpr};

use super::super::eigenpacket::Element;
use super::super::{multiply, Integer, Magnitude, INTEGER_WORDS};
use super::{binomial_table, eigenperiod::Eigenperiod};

pub(crate) enum Expression {
    Symbol(Label, [i128; 4]),
    Coefficient(usize, Box<ParsedExpr>),
    Twist(Label, i128, usize, usize),
}

#[derive(Clone, Copy, Eq, PartialEq)]
struct Cusp(i128, i128);

impl Cusp {
    const INFINITY: Self = Self(1, 0);
    const ZERO: Self = Self(0, 1);

    fn new(numerator: i128, denominator: i128) -> Result<Self> {
        if denominator == 0 {
            if numerator == 0 {
                return Err(Error::plain("a modular-symbol cusp cannot be 0/0"));
            }
            return Ok(Self::INFINITY);
        }
        let rational = Rational::ratio(numerator, denominator).map_err(|_| cusp_range())?;
        Ok(Self(rational.numerator, rational.denominator))
    }

    fn next(self, older: Self, quotient: i128) -> Result<Self> {
        let coordinate = |current: i128, previous| {
            current
                .checked_mul(quotient)
                .and_then(|product| product.checked_add(previous))
                .ok_or_else(cusp_range)
        };
        Ok(Self(
            coordinate(self.0, older.0)?,
            coordinate(self.1, older.1)?,
        ))
    }
}

pub(crate) struct Symbol(Vec<Element>);

impl Symbol {
    pub(crate) fn new(
        model: &mut Model,
        period: &Eigenperiod,
        [a, b, c, d]: [i128; 4],
    ) -> Result<Self> {
        let left = Cusp::new(a, b)?;
        let right = Cusp::new(c, d)?;
        let mut symbol = Self::zero(period);
        if left != right {
            symbol.accumulate_to_infinity(model, period, left, false)?;
            symbol.accumulate_to_infinity(model, period, right, true)?;
        }
        Ok(symbol)
    }

    pub(crate) fn zero(period: &Eigenperiod) -> Self {
        Self(vec![
            Element::zero(period.packet.degree());
            period.coefficients.len()
        ])
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(
            1 + self.0.len() * (1 + self.0[0].fraction_parts().0.len()) * INTEGER_WORDS,
        )?;
        Ok(Value::HeckeSymbol(Rc::new(self)))
    }

    pub(crate) fn elements(&self) -> &[Element] {
        &self.0
    }

    pub(crate) fn coefficient(&self, model: &mut Model, j: usize) -> Result<Value> {
        self.0
            .get(j)
            .cloned()
            .unwrap_or_else(|| Element::zero(self.0[0].fraction_parts().0.len()))
            .persist(model)
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let coefficients = self
            .0
            .iter()
            .zip(&rhs.0)
            .map(|(a, b)| a.add(model, b))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self(coefficients))
    }

    pub(crate) fn negate(&self, model: &mut Model) -> Result<Self> {
        let coefficients = self
            .0
            .iter()
            .map(|entry| entry.negate(model))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self(coefficients))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        for (a, b) in self.0.iter().zip(&rhs.0) {
            let cmp = a.compare(model, b)?;
            if !cmp.is_eq() {
                return Ok(cmp);
            }
        }
        Ok(Ordering::Equal)
    }

    fn accumulate_to_infinity(
        &mut self,
        model: &mut Model,
        period: &Eigenperiod,
        Cusp(mut numerator, mut denominator): Cusp,
        subtract: bool,
    ) -> Result<()> {
        let (mut older, mut previous) = (Cusp::ZERO, Cusp::INFINITY);
        let mut positive = true;
        while denominator != 0 {
            model.spend(1)?;
            let quotient = numerator.div_euclid(denominator);
            let remainder = numerator.rem_euclid(denominator);
            let current = previous.next(older, quotient)?;
            // The edge is current -> previous. The columns (previous,
            // current) have alternating determinants, starting at +1.
            // Negate the representative of current when necessary.
            let (b, d) = if positive {
                (current.0, current.1)
            } else {
                (negative(current.0)?, negative(current.1)?)
            };
            let inverse = [d, negative(b)?, negative(previous.1)?, previous.0];
            let edge = substitute(model, &period.coefficients, inverse)?;
            for (target, mut value) in self.0.iter_mut().zip(edge) {
                if !value.is_zero() {
                    if subtract {
                        value = value.negate(model)?;
                    }
                    *target = target.add(model, &value)?;
                }
            }
            (older, previous) = (previous, current);
            (numerator, denominator) = (denominator, remainder);
            positive = !positive;
        }
        Ok(())
    }

    fn shifted_coefficient(
        &self,
        model: &mut Model,
        index: usize,
        a: u64,
        conductor: u64,
        binomial: &[Vec<i128>],
    ) -> Result<Element> {
        let mut value = Element::zero(self.0[0].fraction_parts().0.len());
        let numerator = Integer::from_i128(i128::from(a));
        let denominator = Magnitude::from_small(u128::from(conductor));
        for (k, entry) in self.0.iter().enumerate().skip(index).rev() {
            value = value.scale(model, numerator, denominator)?;
            let term = integer_scale(model, entry, Integer::from_i128(binomial[k][index]))?;
            value = value.add(model, &term)?;
        }
        Ok(value)
    }
}

pub(super) fn scaled_moment(
    model: &mut Model,
    period: &Eigenperiod,
    modulus: u64,
    r: usize,
    residue: u64,
    binomial: &[Vec<i128>],
) -> Result<Element> {
    let w = period.grade.degree();
    let j = w - r;
    let cusp = Cusp::new(i128::from(residue), i128::from(modulus))?;
    let mut symbol = Symbol::zero(period);
    symbol.accumulate_to_infinity(model, period, cusp, false)?;

    // m^r times the shifted coefficient is the sum of
    // binom(k,j) c_k a^(k-j) m^(w-k), with m=p^n.
    // Evaluate by integral Horner to avoid large temporary denominators.
    let mut value = Element::zero(period.packet.degree());
    let mut power = Integer::from_i128(1);
    let a = Integer::from_i128(i128::from(residue));
    let m = Integer::from_i128(i128::from(modulus));
    for (k, entry) in symbol.0.iter().enumerate().skip(j).rev() {
        value = integer_scale(model, &value, a)?;
        if !entry.is_zero() {
            let scale = multiply(model, power, Integer::from_i128(binomial[k][j]))?;
            let term = integer_scale(model, entry, scale)?;
            value = value.add(model, &term)?;
        }
        if k != j {
            power = multiply(model, power, m)?;
        }
    }
    if j.is_multiple_of(2) {
        Ok(value)
    } else {
        value.negate(model)
    }
}

pub(crate) fn twist_ratio(
    model: &mut Model,
    period: &Eigenperiod,
    discriminant: i128,
    s: usize,
    t: usize,
) -> Result<Value> {
    let character = QuadraticCharacter::new(model, discriminant)?;
    let w = period.grade.degree();
    let top = w + 1;
    if !(1..=top).contains(&s) || !(1..=top).contains(&t) {
        return Err(Error::plain(
            "twisted critical L-value indices must be between 1 and k−1",
        ));
    }
    let (j, ell) = (top - s, top - t);
    let sign = period.grade.sign;
    if !sign.includes(ell) || sign.includes(j) == character.odd() {
        return Err(Error::plain(
            "twisted critical L-value indices and character must match the eigenperiod sign",
        ));
    }
    let reference = &period.coefficients[ell];
    if reference.is_zero() {
        return Err(Error::plain("the denominator critical L-value vanishes"));
    }
    // For square conductor twist of a level-one form the root sign is
    // (-1)^(k/2) * chi(-1).
    if s == period.grade.weight / 2 && (s % 2 == 1) != character.odd() {
        return Element::zero(period.packet.degree()).persist(model);
    }

    let conductor = character.conductor;
    let work = usize::try_from(conductor).map_err(|_| character_limit())?;
    model.spend(work).map_err(|_| character_limit())?;
    let binomial = binomial_table(model, w)?;
    let mut sum = Element::zero(period.packet.degree());
    for a in 0..conductor {
        let sign = character.value(a);
        if sign == 0 {
            continue;
        }
        let cusp = Cusp(i128::from(a), i128::from(conductor));
        // A nonzero primitive character value implies a/N is reduced.
        let mut symbol = Symbol::zero(period);
        symbol.accumulate_to_infinity(model, period, cusp, false)?;
        let mut term = symbol.shifted_coefficient(model, j, a, conductor, &binomial)?;
        if sign < 0 {
            term = term.negate(model)?;
        }
        sum = sum.add(model, &term)?;
    }
    let inverse = period.packet.inverse_element(model, reference)?;
    let ratio = period.packet.multiply_elements(model, &sum, &inverse)?;
    let numerator = if ell % 2 == j % 2 {
        binomial[w][ell]
    } else {
        -binomial[w][ell]
    };
    ratio
        .scale(
            model,
            Integer::from_i128(numerator),
            Magnitude::from_small(binomial[w][j] as u128),
        )?
        .persist(model)
}

fn substitute(model: &mut Model, coeffs: &[Element], matrix: [i128; 4]) -> Result<Vec<Element>> {
    let [a, b, c, d] = matrix.map(Integer::from_i128);
    let degree = coeffs[0].fraction_parts().0.len();
    let mut q = vec![coeffs.last().unwrap().clone()];
    let mut right = vec![Integer::from_i128(1)];
    // Homogeneous Horner in the two integral linear forms.
    for coefficient in coeffs[..coeffs.len() - 1].iter().rev() {
        model.spend(1 + q.len())?;
        let mut next = vec![Element::zero(degree); q.len() + 1];
        for (j, term) in q.iter().enumerate() {
            if !term.is_zero() {
                let constant = integer_scale(model, term, b)?;
                next[j] = next[j].add(model, &constant)?;
                next[j + 1] = integer_scale(model, term, a)?;
            }
        }
        let mut power = vec![Integer::from_i128(0); right.len() + 1];
        for (j, &entry) in right.iter().enumerate() {
            let constant = multiply(model, entry, d)?;
            power[j] = power[j].add(constant)?;
            power[j + 1] = multiply(model, entry, c)?;
        }
        if !coefficient.is_zero() {
            for (term, &entry) in next.iter_mut().zip(&power) {
                let extra = integer_scale(model, coefficient, entry)?;
                *term = term.add(model, &extra)?;
            }
        }
        (q, right) = (next, power);
    }
    Ok(q)
}

fn integer_scale(model: &mut Model, entry: &Element, scalar: Integer) -> Result<Element> {
    entry.scale(model, scalar, Magnitude::from_small(1))
}

fn negative(n: i128) -> Result<i128> {
    n.checked_neg().ok_or_else(cusp_range)
}

fn cusp_range() -> Error {
    Error::plain("modular-symbol cusp arithmetic exceeds the signed 128-bit range")
}

struct QuadraticCharacter {
    discriminant: i128,
    conductor: u64,
}

impl QuadraticCharacter {
    fn new(model: &mut Model, discriminant: i128) -> Result<Self> {
        let conductor =
            u64::try_from(discriminant.unsigned_abs()).map_err(|_| character_limit())?;
        let radicand = if discriminant.rem_euclid(4) == 1 {
            conductor
        } else if discriminant % 4 == 0 && matches!((discriminant / 4).rem_euclid(4), 2 | 3) {
            conductor / 4
        } else {
            return Err(not_fundamental());
        };
        let mut squarefree = radicand;
        let mut divisor = 2;
        while divisor <= squarefree / divisor {
            model.spend(1).map_err(|_| character_limit())?;
            if squarefree.is_multiple_of(divisor) {
                squarefree /= divisor;
                if squarefree.is_multiple_of(divisor) {
                    return Err(not_fundamental());
                }
            }
            divisor += if divisor == 2 { 1 } else { 2 };
        }
        Ok(Self {
            discriminant,
            conductor,
        })
    }

    fn odd(&self) -> bool {
        self.discriminant < 0
    }

    fn value(&self, a: u64) -> i8 {
        if self.conductor == 1 {
            return 1;
        }
        if gcd(u128::from(a), u128::from(self.conductor)) != 1 {
            return 0;
        }
        let power = a.trailing_zeros();
        let mut b = a >> power;
        let mut a = self.discriminant.rem_euclid(i128::from(b)) as u64;
        let mut sign = if power % 2 == 1 && matches!(self.discriminant.rem_euclid(8), 3 | 5) {
            -1
        } else {
            1
        };
        // Jacobi reciprocity on the odd denominator.
        while a != 0 {
            let exponent = a.trailing_zeros();
            a >>= exponent;
            if exponent % 2 == 1 && matches!(b % 8, 3 | 5) {
                sign = -sign;
            }
            if a % 4 == 3 && b % 4 == 3 {
                sign = -sign;
            }
            (a, b) = (b % a, a);
        }
        debug_assert_eq!(b, 1);
        sign
    }
}

fn character_limit() -> Error {
    Error::plain("quadratic twist conductor exceeds the exact computation limit")
}

fn not_fundamental() -> Error {
    Error::plain("quadratic twist requires a fundamental discriminant (including 1)")
}
