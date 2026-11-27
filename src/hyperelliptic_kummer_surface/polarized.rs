use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::jacobian_level::Basis;
use crate::model::{Model, Value};
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};

use super::nikulin;

const NODES: usize = 16;
pub(crate) const COEFFICIENTS: usize = NODES + 1;
type Coefficients = [i128; COEFFICIENTS];

pub(crate) enum Generator {
    Hyperplane,
    Node(usize),
    Trope(usize),
}

pub(crate) enum Expression {
    Generator(Label, Generator),
    Literal(Label, Box<Coefficients>),
    Scale(i128, Box<ParsedExpr>),
    Frobenius(Label, usize, Box<ParsedExpr>),
}

pub(crate) enum Command {
    Pair(Box<[ParsedExpr; 2]>),
    Fujiki(Fujiki, Box<[ParsedExpr; 2]>),
}

pub(crate) fn basis(model: &Model, name: &Label) -> Result<Rc<Basis>> {
    nikulin::basis_for(
        model,
        name,
        "marked Kummer divisors require a certified genus-two level-two basis",
    )
}

fn trope_origin(model: &mut Model, basis: &Basis) -> Result<u16> {
    let mut support = 0u16;
    for mask in 0..NODES {
        if basis.binary_combination(model, mask as u8)?.u.len() <= 2 {
            support |= 1 << mask;
        }
    }
    if support.count_ones() != 6 {
        return Err(Error::plain(
            "marked Kummer configuration must have six nodes on the zero trope",
        ));
    }
    Ok(support)
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Divisor(Coefficients);

impl Divisor {
    pub(crate) fn zero() -> Self {
        Self([0; COEFFICIENTS])
    }

    pub(crate) fn generator(
        model: &mut Model,
        basis: &Basis,
        generator: Generator,
    ) -> Result<Self> {
        let mut result = Self::zero();
        match generator {
            Generator::Hyperplane => result.0[0] = 2,
            Generator::Node(index) => result.0[index + 1] = 2,
            Generator::Trope(index) => {
                let origin = trope_origin(model, basis)?;
                result.0[0] = 1;
                for (node, coefficient) in result.0[1..].iter_mut().enumerate() {
                    if origin & (1 << (node ^ index)) != 0 {
                        *coefficient = -1;
                    }
                }
            }
        }
        Ok(result)
    }

    pub(crate) fn literal(
        model: &mut Model,
        basis: &Basis,
        coefficients: Coefficients,
    ) -> Result<Self> {
        model.spend(2 * NODES)?;
        let mut parity = 0;
        for (node, &coefficient) in coefficients[1..].iter().enumerate() {
            parity |= ((coefficient & 1) as u16) << node;
        }
        if coefficients[0] & 1 != 0 {
            parity ^= trope_origin(model, basis)?;
        }
        if !nikulin::affine_code().any(|table| table == parity) {
            return Err(Error::plain(
                "half-vector is outside the integral marked Kummer divisor lattice",
            ));
        }
        Ok(Self(coefficients))
    }

    pub(crate) fn coefficients(&self) -> &Coefficients {
        &self.0
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + 2 * COEFFICIENTS)?;
        Ok(Value::KummerDiv(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<std::cmp::Ordering> {
        model.spend(COEFFICIENTS)?;
        Ok(self.cmp(other))
    }

    pub(crate) fn add(&self, model: &mut Model, other: &Self) -> Result<Self> {
        model.spend(COEFFICIENTS)?;
        let mut result = Self::zero();
        for (i, coefficient) in result.0.iter_mut().enumerate() {
            *coefficient = self.0[i]
                .checked_add(other.0[i])
                .ok_or_else(coefficient_overflow)?;
        }
        Ok(result)
    }

    pub(crate) fn scale(&self, model: &mut Model, multiplier: i128) -> Result<Self> {
        model.spend(COEFFICIENTS)?;
        let mut result = Self::zero();
        for (coefficient, &original) in result.0.iter_mut().zip(&self.0) {
            *coefficient = original
                .checked_mul(multiplier)
                .ok_or_else(coefficient_overflow)?;
        }
        Ok(result)
    }

    pub(crate) fn frobenius(
        &self,
        model: &mut Model,
        basis: &Basis,
        base: &Curve,
        power: usize,
    ) -> Result<Self> {
        let permutation = nikulin::NodePermutation::new(model, basis, base, power)?;
        let mut result = Self::zero();
        result.0[0] = self.0[0];
        for (node, &coefficient) in self.0[1..].iter().enumerate() {
            result.0[usize::from(permutation.image(node)) + 1] = coefficient;
        }
        Ok(result)
    }

    fn intersection(&self, model: &mut Model, other: &Self) -> Result<Intersection> {
        model.spend(COEFFICIENTS * 12)?;
        let mut doubled = Intersection::product(self.0[0], other.0[0]).multiply(2)?;
        for (&left, &right) in self.0[1..].iter().zip(&other.0[1..]) {
            doubled = doubled.add(Intersection::product(left, right).negate())?;
        }
        let (magnitude, parity) = doubled.magnitude.divide_u32(2);
        if parity != 0 {
            return Err(Error::plain(
                "marked Kummer divisor pairing is not integral",
            ));
        }
        Ok(Intersection {
            magnitude,
            negative: doubled.negative,
        })
    }

    pub(crate) fn pair(&self, model: &mut Model, other: &Self) -> Result<i128> {
        self.intersection(model, other)?
            .signed()
            .ok_or_else(|| Error::plain("Kummer divisor intersection exceeds signed 128-bit range"))
    }
}

fn coefficient_overflow() -> Error {
    Error::plain("Kummer divisor coefficient exceeds signed 128-bit range")
}

// This width also holds 2(n−1) times the product of two i128 shifts.
type Intersection = Wide<12>;

const FUJIKI_LIMBS: usize = 32;
type FujikiInteger = Wide<FUJIKI_LIMBS>;

#[derive(Clone, Copy)]
pub(crate) struct Fujiki {
    pub(crate) dimension: usize,
    pub(crate) right_exponent: usize,
    pub(crate) left_shift: i128,
    pub(crate) right_shift: i128,
}

impl Fujiki {
    pub(crate) fn evaluate(self, model: &mut Model, x: &Divisor, y: &Divisor) -> Result<i128> {
        let Self {
            dimension: n,
            right_exponent: j,
            left_shift: k,
            right_shift: l,
        } = self;
        if n == 0 {
            return Err(Error::plain(
                "Kummer Fujiki dimension parameter n must be positive",
            ));
        }
        let top = n.checked_mul(2).ok_or_else(work_limit)?;
        let Some(left) = top.checked_sub(j) else {
            return Err(Error::plain(
                "Kummer Fujiki exponent j must lie between zero and 2n",
            ));
        };
        if n == 1 && (k != 0 || l != 0) {
            return Err(Error::plain("Kummer Fujiki shifts must be zero for n=1"));
        }
        let (a, b, x, y, k, l) = if left >= j {
            (left, j, x, y, k, l)
        } else {
            (j, left, y, x, l, k)
        };
        let width = b.checked_add(1).ok_or_else(work_limit)?;
        let work = a
            .checked_add(1)
            .and_then(|rows| rows.checked_mul(width))
            .and_then(|cells| cells.checked_mul(FUJIKI_LIMBS))
            .ok_or_else(work_limit)?;
        model.spend(work).map_err(|_| work_limit())?;
        crate::hw_algebra::storage(3 * width, FUJIKI_LIMBS + 2)?;
        let pairing =
            |model: &mut Model, x: &Divisor, y: &Divisor, k, l| -> Result<FujikiInteger> {
                let correction = Intersection::product(k, l).multiply(2 * (n as i128 - 1))?;
                x.intersection(model, y)?.add(correction.negate())?.resize()
            };
        let xx = pairing(model, x, x, k, k)?;
        let xy = pairing(model, x, y, k, l)?;
        let yy = pairing(model, y, y, l, l)?;
        pairings(a, b, xx, xy, yy)
            .and_then(FujikiInteger::signed)
            .ok_or_else(top_overflow)
    }
}

fn work_limit() -> Error {
    Error::plain("Kummer Fujiki pairing table exceeds the work limit")
}

fn top_overflow() -> Error {
    Error::plain("Kummer Fujiki intersection exceeds the supported exact integer range")
}

fn weighted(
    value: Option<FujikiInteger>,
    pairing: FujikiInteger,
    count: usize,
) -> Option<FujikiInteger> {
    if count == 0 || pairing.is_zero() {
        return Some(FujikiInteger::from_i128(0));
    }
    pairing
        .multiply_wide(value?)
        .and_then(|term| term.multiply(count as i128))
        .ok()
}

fn pairings(
    a: usize,
    b: usize,
    xx: FujikiInteger,
    xy: FujikiInteger,
    yy: FujikiInteger,
) -> Option<FujikiInteger> {
    let zero = Some(FujikiInteger::from_i128(0));
    let mut older = vec![zero; b + 1];
    let mut previous = vec![zero; b + 1];
    let mut current = vec![zero; b + 1];
    previous[0] = Some(FujikiInteger::from_i128(1));
    for j in (2..=b).step_by(2) {
        previous[j] = weighted(previous[j - 2], yy, j - 1);
    }
    for i in 1..=a {
        // Pair one specified left vector with any of the remaining vectors.
        for j in (i % 2..=b).step_by(2) {
            let same = weighted(older[j], xx, i - 1);
            let cross = if j > 0 {
                weighted(previous[j - 1], xy, j)
            } else {
                zero
            };
            current[j] = same
                .zip(cross)
                .and_then(|(same, cross)| same.add(cross).ok());
        }
        std::mem::swap(&mut older, &mut previous);
        std::mem::swap(&mut previous, &mut current);
    }
    previous[b]
}
