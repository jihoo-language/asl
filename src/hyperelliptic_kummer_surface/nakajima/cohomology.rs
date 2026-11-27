use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_ore::gcd;
use crate::syntax::{Label, ParsedExpr};
use crate::wide_integer::Natural;

use super::super::unimodular::Vector;
use super::{tautological, Accumulator, Factor, Integer, Mode, State, LIMBS, POINT};

pub(crate) mod schur;

type Denominator = Natural<LIMBS>;

pub(crate) enum Expression {
    Unit(Label, u128),
    Point(Label, u128),
    FromFock(i128, u128, Box<ParsedExpr>),
    Scale(i128, u128, Box<ParsedExpr>),
    ToFock(i128, Box<ParsedExpr>),
    Cup(Box<[ParsedExpr; 2]>),
    Chern(u128, u128, bool, Box<ParsedExpr>),
    Schur(u128, u128, Vec<u128>, Box<ParsedExpr>),
}

pub(crate) enum Report {
    Integral(Box<ParsedExpr>),
    Denominator(Box<ParsedExpr>),
    ExteriorEuler(u128, u128, Box<[ParsedExpr; 2]>),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Class {
    denominator: u128,
    numerator: State,
}

impl Class {
    pub(crate) fn zero() -> Self {
        Self {
            denominator: 1,
            numerator: State::zero(),
        }
    }

    pub(crate) fn unit(model: &mut Model, length: u128) -> Result<Self> {
        Self::degree_zero(model, length, length)
    }

    pub(crate) fn point(length: u128) -> Self {
        Self {
            numerator: pure(length, POINT),
            denominator: 1,
        }
    }

    fn degree_zero(model: &mut Model, length: u128, factorial_index: u128) -> Result<Self> {
        Ok(Self {
            numerator: pure(length, 0),
            denominator: factorial(model, factorial_index)?
                .small()
                .ok_or_else(denominator_range)?,
        })
    }

    pub(crate) fn denominator(&self) -> u128 {
        self.denominator
    }

    pub(crate) fn numerator(&self) -> &State {
        &self.numerator
    }

    pub(crate) fn top_integral(&self) -> String {
        let numerator = self.numerator.terms.values().next().copied().unwrap_or(0);
        if self.denominator == 1 {
            numerator.to_string()
        } else {
            format!("{numerator}/{}", self.denominator)
        }
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        match self.denominator.cmp(&other.denominator) {
            Ordering::Equal => self.numerator.compare(model, &other.numerator),
            unequal => Ok(unequal),
        }
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        super::check_cells(self.cells())?;
        model.allocate(self.cells())?;
        Ok(Value::HilbertCohomology(Rc::new(self)))
    }

    fn cells(&self) -> usize {
        4 + self.numerator.cells()
    }

    pub(crate) fn add(&self, model: &mut Model, other: &Self) -> Result<Self> {
        self.add_signed(model, other, false)
    }

    fn add_signed(&self, model: &mut Model, other: &Self, subtract: bool) -> Result<Self> {
        let common = gcd(self.denominator, other.denominator);
        let left_multiplier = other.denominator / common;
        let right_multiplier = self.denominator / common;
        let denominator = product(self.denominator, left_multiplier);
        let mut output = Accumulator::default();
        add_scaled(
            model,
            &mut output,
            &self.numerator,
            positive(left_multiplier),
        )?;
        let right_scale = positive(right_multiplier);
        let right_scale = if subtract {
            right_scale.negate()
        } else {
            right_scale
        };
        add_scaled(model, &mut output, &other.numerator, right_scale)?;
        normalize(model, output, denominator)
    }

    pub(crate) fn from_fock(
        model: &mut Model,
        state: &State,
        k: i128,
        denominator: u128,
    ) -> Result<Self> {
        Self::scaled_fock(model, state, k, Denominator::from_small(denominator))
    }

    pub(crate) fn scale(&self, model: &mut Model, k: i128, denominator: u128) -> Result<Self> {
        Self::scaled_fock(
            model,
            &self.numerator,
            k,
            product(self.denominator, denominator),
        )
    }

    fn scaled_fock(
        model: &mut Model,
        state: &State,
        k: i128,
        denominator: Denominator,
    ) -> Result<Self> {
        Self::scaled_wide(model, state, Integer::from_i128(k), denominator)
    }

    fn scaled_wide(
        model: &mut Model,
        state: &State,
        k: Integer,
        denominator: Denominator,
    ) -> Result<Self> {
        if denominator.is_zero() {
            return Err(Error::plain(
                "Hilbert cohomology denominator must be positive",
            ));
        }
        if k.is_zero() || state.terms.is_empty() {
            return Ok(Self::zero());
        }
        let mut output = Accumulator::default();
        add_scaled(model, &mut output, state, k)?;
        normalize(model, output, denominator)
    }

    pub(crate) fn to_fock(&self, model: &mut Model, scale: i128) -> Result<State> {
        let numerator = scale.unsigned_abs();
        if !numerator.is_multiple_of(self.denominator) {
            return Err(Error::plain(
                "Hilbert cohomology scale gives nonintegral Fock coefficients",
            ));
        }
        if numerator == 0 {
            return Ok(State::zero());
        }
        let multiplier = Integer {
            magnitude: Denominator::from_small(numerator / self.denominator),
            negative: scale < 0,
        };
        let mut output = Accumulator::default();
        add_scaled(model, &mut output, &self.numerator, multiplier)?;
        output.finish()
    }

    pub(crate) fn cup(&self, model: &mut Model, other: &Self, length: i128) -> Result<Self> {
        let output = self
            .numerator
            .accumulate_cup(model, &other.numerator, length)?;
        normalize(model, output, product(self.denominator, other.denominator))
    }

    pub(crate) fn chern(model: &mut Model, line: &Vector, n: u128, j: u128) -> Result<Self> {
        if j == 0 {
            return Self::unit(model, n);
        }
        if j > n {
            return Ok(Self::zero());
        }
        let denominator = tautological::denominator(model, line, n, j)?;
        let numerator = tautological::unsigned_chern(model, line, n, j, denominator)?;
        Ok(Self {
            numerator,
            denominator,
        })
    }

    pub(crate) fn character(model: &mut Model, line: &Vector, n: u128, j: u128) -> Result<Self> {
        if j == 0 {
            return if n == 0 {
                Ok(Self::zero())
            } else {
                Self::degree_zero(model, n, n - 1)
            };
        }
        if j.div_ceil(2) > n {
            return Ok(Self::zero());
        }
        let limit = usize::try_from(j).map_err(|_| character_limit())?;
        let power = Self::power_sums(model, line, n, limit)?
            .pop()
            .expect("positive character index");
        power.divide_by_factorial(model, j)
    }

    fn power_sums(model: &mut Model, line: &Vector, n: u128, limit: usize) -> Result<Vec<Self>> {
        let rank = n.min(limit as u128) as usize;
        let mut cells = 8usize
            .saturating_mul(limit)
            .saturating_add(8usize.saturating_mul(rank));
        character_storage(cells)?;
        model
            .spend(limit.saturating_mul(rank))
            .map_err(|_| character_limit())?;
        let length = i128::try_from(n).map_err(|_| super::grade_range())?;
        let mut chern = Vec::with_capacity(rank);
        for i in 1..=rank {
            let class = Self::chern(model, line, n, i as u128)?;
            cells += class.cells();
            character_storage(cells)?;
            chern.push(class);
        }
        let mut powers = Vec::with_capacity(limit);
        for degree in 1..=limit {
            let mut sum = Self::zero();
            for i in 1..=rank.min(degree - 1) {
                let product = chern[i - 1].cup(model, &powers[degree - i - 1], length)?;
                sum = sum.add_signed(model, &product, i.is_multiple_of(2))?;
            }
            if degree <= rank {
                let sign = if degree % 2 == 1 { 1i128 } else { -1 };
                let diagonal = chern[degree - 1].scale(model, sign * degree as i128, 1)?;
                sum = sum.add(model, &diagonal)?;
            }
            cells += sum.cells();
            character_storage(cells)?;
            powers.push(sum);
        }
        Ok(powers)
    }

    fn divide_by_factorial(&self, model: &mut Model, index: u128) -> Result<Self> {
        if self.numerator.terms.is_empty() {
            return Ok(Self::zero());
        }
        let denominator = factorial(model, index)?
            .multiply_u128(self.denominator)
            .ok_or_else(denominator_range)?;
        Self::scaled_fock(model, &self.numerator, 1, denominator)
    }
}

fn pure(length: u128, color: u8) -> State {
    if length == 0 {
        return State::vacuum();
    }
    let factor = Factor {
        mode: Mode { length: 1, color },
        exponent: length,
    };
    State {
        terms: BTreeMap::from([(vec![factor], 1)]),
    }
}

fn factorial(model: &mut Model, index: u128) -> Result<Denominator> {
    (2..=index).try_fold(Denominator::from_small(1), |product, i| {
        model.spend(LIMBS)?;
        product.multiply_u128(i).ok_or_else(denominator_range)
    })
}

fn positive(value: u128) -> Integer {
    Integer {
        magnitude: Denominator::from_small(value),
        negative: false,
    }
}

fn product(left: u128, right: u128) -> Denominator {
    Denominator::from_small(left)
        .multiply_u128(right)
        .expect("256-bit denominator product")
}

fn add_scaled(
    model: &mut Model,
    output: &mut Accumulator,
    state: &State,
    scale: Integer,
) -> Result<()> {
    for (key, &coefficient) in &state.terms {
        model.spend(1 + key.len() + LIMBS)?;
        let coefficient = scale
            .multiply(coefficient)
            .map_err(|_| coordinate_range())?;
        output.add(key.clone(), coefficient)?;
    }
    Ok(())
}

fn normalize(
    model: &mut Model,
    output: Accumulator,
    mut denominator: Denominator,
) -> Result<Class> {
    let one = Denominator::from_small(1);
    let mut content = denominator;
    let mut nonzero = false;
    for coefficient in output.terms.values() {
        if !coefficient.is_zero() {
            nonzero = true;
            if content == one {
                break;
            }
            content = content.gcd(model, coefficient.magnitude)?;
        }
    }
    if !nonzero {
        return Ok(Class::zero());
    }
    if content != one {
        denominator = divide(model, denominator, content)?;
    }
    let denominator = denominator.small().ok_or_else(denominator_range)?;
    let mut terms = BTreeMap::new();
    for (key, mut coefficient) in output.terms {
        if coefficient.is_zero() {
            continue;
        }
        if content != one {
            coefficient.magnitude = divide(model, coefficient.magnitude, content)?;
        }
        let value = coefficient.signed().ok_or_else(coordinate_range)?;
        terms.insert(key, value);
    }
    Ok(Class {
        denominator,
        numerator: State { terms },
    })
}

fn divide(model: &mut Model, value: Denominator, divisor: Denominator) -> Result<Denominator> {
    Ok(value
        .divide_exact(model, divisor)?
        .expect("the coefficient content divides each entry"))
}

fn character_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        Err(character_limit())
    } else {
        Ok(())
    }
}

fn character_limit() -> Error {
    Error::plain("Hilbert Chern character exceeds the Newton computation limit")
}

fn denominator_range() -> Error {
    Error::plain("Hilbert cohomology denominator exceeds the unsigned 128-bit range")
}

fn coordinate_range() -> Error {
    Error::plain("Hilbert cohomology coefficient exceeds the signed 128-bit range")
}
