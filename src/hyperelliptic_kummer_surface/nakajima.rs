use std::cmp::Ordering;
use std::collections::{btree_map::Entry, BTreeMap};
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};

use super::unimodular;

mod boundary;
mod cup;
pub(crate) mod primary;
pub(crate) mod tautological;
mod virasoro;
mod wick;

pub(crate) const COLORS: usize = unimodular::RANK + 2;
const POINT: u8 = COLORS as u8 - 1;
const LIMBS: usize = 16;
type Integer = Wide<LIMBS>;
type Monomial = Vec<Factor>;

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) basis: Rc<str>,
    pub(crate) length: i128,
    pub(crate) degree: i128,
}

impl Grade {
    pub(crate) fn supports_states(&self) -> bool {
        self.length >= 0
            && self.degree >= 0
            && self.degree % 2 == 0
            && (self.degree as u128).div_ceil(4) <= self.length as u128
    }

    pub(crate) fn complementary(&self, rhs: &Self) -> bool {
        if self.length != rhs.length || !self.supports_states() || !rhs.supports_states() {
            return false;
        }
        let total = self.degree as u128 + rhs.degree as u128;
        total.is_multiple_of(4) && total / 4 == self.length as u128
    }
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Mode {
    pub(crate) length: i128,
    pub(crate) color: u8,
}

impl Mode {
    pub(crate) fn new(length: u128, color: usize) -> Result<Self> {
        let length = i128::try_from(length).map_err(|_| grade_range())?;
        if length == 0 {
            return Err(Error::plain("Nakajima mode length must be positive"));
        }
        if color >= COLORS {
            return Err(Error::plain("Nakajima color must be between 0 and 23"));
        }
        Ok(Self {
            length,
            color: color as u8,
        })
    }

    fn color_degree(self) -> i128 {
        match self.color {
            0 => 0,
            POINT => 4,
            _ => 2,
        }
    }
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Factor {
    pub(crate) mode: Mode,
    pub(crate) exponent: u128,
}

pub(crate) enum Expression {
    Vacuum(Label),
    Monomial(Label, Vec<Factor>),
    Cup(Box<[ParsedExpr; 2]>),
    PrimaryGenerator(Label, i128, i128, usize),
    TautologicalChern(u128, u128, i128, Box<ParsedExpr>),
    Apply(Operation, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Create(Mode),
    Annihilate(Mode),
    Scale(i128),
    Virasoro(i128),
    Boundary,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Create(_) => "rvNakajimaCreate",
            Self::Annihilate(_) => "rvNakajimaAnn",
            Self::Scale(_) => "rvNakajimaScale",
            Self::Virasoro(_) => "rvVirasoro",
            Self::Boundary => "rvNakajimaBoundary",
        }
    }

    pub(crate) fn grades(self, source: &Grade) -> Result<(i128, i128)> {
        let (mode, sign) = match self {
            Self::Create(mode) => (mode, 1),
            Self::Annihilate(mode) => (mode, -1),
            Self::Scale(_) => return Ok((source.length, source.degree)),
            Self::Boundary => {
                let degree = source.degree.checked_add(2).ok_or_else(grade_range)?;
                return Ok((source.length, degree));
            }
            Self::Virasoro(index) => {
                let length = source.length.checked_sub(index).ok_or_else(grade_range)?;
                let degree = Integer::from_i128(source.degree)
                    .add(Integer::product(-2, index))?
                    .signed()
                    .ok_or_else(grade_range)?;
                return Ok((length, degree));
            }
        };
        let length = source
            .length
            .checked_add(sign * mode.length)
            .ok_or_else(grade_range)?;
        let degree = Integer::from_i128(source.degree)
            .add(Integer::product(2 * sign, mode.length))?
            .add(Integer::from_i128(mode.color_degree() - 2))?
            .signed()
            .ok_or_else(grade_range)?;
        Ok((length, degree))
    }
}

pub(crate) struct Gram([[i128; COLORS]; COLORS]);

impl Gram {
    pub(crate) fn new(model: &mut Model) -> Result<Self> {
        let middle = unimodular::integral_gram(model)?;
        model.allocate(COLORS * COLORS)?;
        let mut result = [[0; COLORS]; COLORS];
        result[0][POINT as usize] = 1;
        result[POINT as usize][0] = 1;
        for (i, row) in middle.iter().enumerate() {
            result[1 + i][1..=unimodular::RANK].copy_from_slice(row);
        }
        Ok(Self(result))
    }

    fn entry(&self, left: u8, right: u8) -> i128 {
        self.0[left as usize][right as usize]
    }
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct State {
    terms: BTreeMap<Monomial, i128>,
}

impl State {
    pub(crate) fn zero() -> Self {
        Self {
            terms: BTreeMap::new(),
        }
    }

    pub(crate) fn vacuum() -> Self {
        Self {
            terms: BTreeMap::from([(Vec::new(), 1)]),
        }
    }

    pub(crate) fn monomial(model: &mut Model, factors: Vec<Factor>) -> Result<(Self, i128, i128)> {
        let mut combined = BTreeMap::<Mode, u128>::new();
        let (mut length, mut degree) = (0u128, 0u128);
        for Factor { mode, exponent } in factors {
            model.spend(1 + combined.len().max(1).ilog2() as usize)?;
            if exponent == 0 {
                continue;
            }
            let m = mode.length as u128;
            let weight = (2 * (m - 1))
                .checked_add(mode.color_degree() as u128)
                .ok_or_else(grade_range)?;
            length = add_weight(length, m, exponent)?;
            degree = add_weight(degree, weight, exponent)?;
            let count = combined.entry(mode).or_default();
            *count = count.checked_add(exponent).ok_or_else(grade_range)?;
        }
        let monomial = combined
            .into_iter()
            .map(|(mode, exponent)| Factor { mode, exponent })
            .collect();
        let length = i128::try_from(length).map_err(|_| grade_range())?;
        let degree = i128::try_from(degree).map_err(|_| grade_range())?;
        Ok((
            Self {
                terms: BTreeMap::from([(monomial, 1)]),
            },
            length,
            degree,
        ))
    }

    pub(super) fn from_monomial_coefficients(
        model: &mut Model,
        terms: impl IntoIterator<Item = (Vec<Factor>, i128)>,
    ) -> Result<Self> {
        let mut state = Self::zero();
        for (factors, coefficient) in terms {
            if coefficient == 0 {
                continue;
            }
            let (monomial, _, _) = Self::monomial(model, factors)?;
            let key = monomial.terms.into_keys().next().expect("one monomial");
            model.spend(1 + state.terms.len().max(1).ilog2() as usize)?;
            state.add_term(key, coefficient)?;
        }
        Ok(state)
    }

    pub(crate) fn terms(&self) -> impl Iterator<Item = (&[Factor], i128)> {
        self.terms
            .iter()
            .map(|(key, &coefficient)| (key.as_slice(), coefficient))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(self.cells().min(rhs.cells()))?;
        Ok(self.cmp(rhs))
    }

    fn cells(&self) -> usize {
        1 + self
            .terms
            .keys()
            .map(|key| 6 + 5 * key.len())
            .sum::<usize>()
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        check_cells(self.cells())?;
        model.allocate(self.cells())?;
        Ok(Value::Nakajima(Rc::new(self)))
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        model.spend(self.cells() + rhs.cells())?;
        let mut result = Self {
            terms: self.terms.clone(),
        };
        for (key, &coefficient) in &rhs.terms {
            result.add_term(key.clone(), coefficient)?;
        }
        Ok(result)
    }

    fn add_term(&mut self, key: Monomial, coefficient: i128) -> Result<()> {
        match self.terms.entry(key) {
            Entry::Occupied(mut entry) => {
                let value = entry
                    .get()
                    .checked_add(coefficient)
                    .ok_or_else(coordinate_range)?;
                if value == 0 {
                    entry.remove();
                } else {
                    *entry.get_mut() = value;
                }
            }
            Entry::Vacant(entry) => {
                entry.insert(coefficient);
            }
        }
        Ok(())
    }

    pub(crate) fn scale(&self, model: &mut Model, scalar: i128) -> Result<Self> {
        if scalar == 0 {
            return Ok(Self::zero());
        }
        model.spend(self.cells())?;
        let terms = self
            .terms
            .iter()
            .map(|(key, &coefficient)| {
                let value = scalar
                    .checked_mul(coefficient)
                    .ok_or_else(coordinate_range)?;
                Ok((key.clone(), value))
            })
            .collect::<Result<_>>()?;
        Ok(Self { terms })
    }

    pub(crate) fn create(&self, model: &mut Model, mode: Mode) -> Result<Self> {
        model.spend(self.cells())?;
        let terms = self
            .terms
            .iter()
            .map(|(key, &coefficient)| {
                let mut key = key.clone();
                insert_factor(&mut key, mode)?;
                Ok((key, coefficient))
            })
            .collect::<Result<_>>()?;
        Ok(Self { terms })
    }

    pub(crate) fn annihilate(&self, model: &mut Model, gram: &Gram, mode: Mode) -> Result<Self> {
        model.spend(self.cells())?;
        let mut contracted = Accumulator::default();
        for (monomial, &coefficient) in &self.terms {
            let start = monomial.partition_point(|f| f.mode.length < mode.length);
            for (index, f) in monomial
                .iter()
                .enumerate()
                .skip(start)
                .take_while(|(_, f)| f.mode.length == mode.length)
            {
                let pair = gram.entry(mode.color, f.mode.color);
                if pair == 0 {
                    continue;
                }
                model.spend(1 + monomial.len() + LIMBS)?;
                let mut key = monomial.clone();
                remove_factor(&mut key, index);
                let contribution = Integer::product(coefficient, mode.length)
                    .multiply(f.exponent as i128)
                    .and_then(|value| value.multiply(pair))
                    .map_err(|_| coordinate_range())?;
                contracted.add(key, contribution)?;
            }
        }
        contracted.finish()
    }

    pub(crate) fn pair(&self, model: &mut Model, gram: &Gram, rhs: &Self) -> Result<i128> {
        wick::pair(model, gram, self, rhs)
    }
}

fn insert_factor(key: &mut Monomial, mode: Mode) -> Result<()> {
    match key.binary_search_by_key(&mode, |factor| factor.mode) {
        Ok(index) => {
            let count = &mut key[index].exponent;
            *count = count.checked_add(1).ok_or_else(grade_range)?;
        }
        Err(index) => key.insert(index, Factor { mode, exponent: 1 }),
    }
    Ok(())
}

fn length_weight(factor: &Factor, count: u128) -> Result<i128> {
    // The input length grade already bounds the sum of all occupied lengths.
    let count = i128::try_from(count).map_err(|_| coordinate_range())?;
    factor
        .mode
        .length
        .checked_mul(count)
        .ok_or_else(coordinate_range)
}

fn pair_weight(left: &Factor, right: &Factor, same: bool) -> Result<Integer> {
    let (mut a, mut b) = (left.exponent, right.exponent);
    if same {
        // Choose an unordered pair of copies of the same mode.
        b -= 1;
        if a.is_multiple_of(2) {
            a /= 2;
        } else {
            b /= 2;
        }
    }
    Ok(Integer::product(
        length_weight(left, a)?,
        length_weight(right, b)?,
    ))
}

fn remove_factor(key: &mut Monomial, index: usize) {
    if key[index].exponent == 1 {
        key.remove(index);
    } else {
        key[index].exponent -= 1;
    }
}

#[derive(Default)]
struct Accumulator {
    terms: BTreeMap<Monomial, Integer>,
    cells: usize,
}

impl Accumulator {
    fn add(&mut self, key: Monomial, contribution: Integer) -> Result<()> {
        match self.terms.entry(key) {
            Entry::Occupied(mut entry) => {
                *entry.get_mut() = entry
                    .get()
                    .add(contribution)
                    .map_err(|_| coordinate_range())?;
            }
            Entry::Vacant(entry) => {
                self.cells += LIMBS + 6 + 5 * entry.key().len();
                check_cells(self.cells)?;
                entry.insert(contribution);
            }
        }
        Ok(())
    }

    fn finish(self) -> Result<State> {
        let mut terms = BTreeMap::new();
        for (key, coefficient) in self.terms {
            if !coefficient.is_zero() {
                terms.insert(key, coefficient.signed().ok_or_else(coordinate_range)?);
            }
        }
        Ok(State { terms })
    }
}

fn add_weight(total: u128, weight: u128, exponent: u128) -> Result<u128> {
    weight
        .checked_mul(exponent)
        .and_then(|value| total.checked_add(value))
        .ok_or_else(grade_range)
}

fn check_cells(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "Nakajima polynomial exceeds the coefficient storage limit",
        ));
    }
    Ok(())
}

pub(crate) fn grade_range() -> Error {
    Error::plain("Nakajima grade exceeds signed 128-bit range")
}

fn coordinate_range() -> Error {
    Error::plain("Nakajima coordinate exceeds signed 128-bit range")
}
