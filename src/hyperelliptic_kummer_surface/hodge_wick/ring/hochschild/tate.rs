use std::cmp::Ordering;
use std::fmt::Write;

use crate::model::MAX_SHOWN_BYTES;

use super::{
    storage, BTreeMap, Entry, Error, Label, Model, ParsedExpr, Rational, Rc, Result, Value, LIMBS,
};

mod bar;
pub(crate) mod matrix_factorization;
mod resolution;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) enum Mode {
    Ordinary,
    Stable,
}

impl Mode {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Ordinary => "ordinary",
            Self::Stable => "stable",
        }
    }
}

pub(crate) enum Command {
    Declare(Label, Vec<usize>),
    Period(Label, usize, Mode, i128, usize, usize),
    TensorPeriod(Label, Vec<usize>, Vec<usize>, Vec<usize>),
    Check(Label, usize, Mode, i128),
    ComplexMap(Label, usize, Mode, i128, bool, usize),
}

#[derive(Clone, Copy)]
pub(crate) enum Bilinear {
    Cup,
    Bracket,
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Bv,
    Hodge(i128),
    Scale(Entry),
}

pub(crate) enum Expression {
    Term(Label, Mode, Vec<u128>, Vec<usize>, Vec<i128>),
    Bilinear(Bilinear, Box<[ParsedExpr; 2]>),
    Unary(Unary, Box<ParsedExpr>),
    Bar(Label, i128, Box<ParsedExpr>),
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) name: Rc<str>,
    pub(crate) mode: Mode,
    pub(crate) degree: i128,
    pub(crate) shift: i128,
}

impl Grade {
    pub(crate) fn target(&self, rhs: &Self, operation: Bilinear) -> Result<Self> {
        if self.name != rhs.name || self.mode != rhs.mode {
            return Err(Error::plain(
                "Tate operands require the same named complete intersection and mode",
            ));
        }
        let degree = match operation {
            Bilinear::Cup => self.degree.checked_add(rhs.degree),
            Bilinear::Bracket => self
                .degree
                .checked_sub(1)
                .and_then(|left| left.checked_add(rhs.degree))
                .or_else(|| {
                    rhs.degree
                        .checked_sub(1)
                        .and_then(|right| self.degree.checked_add(right))
                }),
        }
        .ok_or_else(grading_range)?;
        let shift = self
            .shift
            .checked_add(rhs.shift)
            .ok_or_else(grading_range)?;
        Ok(Self {
            degree,
            shift,
            ..self.clone()
        })
    }

    pub(crate) fn lower(&self) -> Result<Self> {
        let degree = self.degree.checked_sub(1).ok_or_else(grading_range)?;
        Ok(Self {
            degree,
            ..self.clone()
        })
    }

    pub(crate) fn bar_grade(&self, basis: Rc<str>, length: i128) -> Result<super::Grade> {
        if self.mode != Mode::Ordinary {
            return Err(Error::plain(
                "the normalized bar comparison requires ordinary Tate coordinates",
            ));
        }
        Ok(super::Grade {
            basis,
            length,
            arity: self.degree,
            shift: self.shift,
        })
    }
}

fn grading_range() -> Error {
    Error::plain("Verbitsky Tate degree or internal shift exceeds the signed integer range")
}

pub(crate) struct Definition {
    powers: Vec<usize>,
}

impl Definition {
    pub(crate) fn new(model: &mut Model, powers: Vec<usize>) -> Result<Self> {
        if powers.is_empty() || powers.iter().any(|&power| power < 2) {
            return Err(Error::plain(
                "Verbitsky Tate complete intersection needs relation exponents at least two",
            ));
        }
        storage(powers.len().saturating_add(12))?;
        model.allocate(powers.len() + 12)?;
        Ok(Self { powers })
    }

    fn factor(&self, index: usize) -> Result<usize> {
        self.powers
            .get(index)
            .copied()
            .ok_or_else(|| Error::plain("Verbitsky Tate factor index is out of range"))
    }

    pub(crate) fn validate_mode(&self, mode: Mode) -> Result<()> {
        if self.powers.len() != 1 && mode == Mode::Stable {
            Err(Error::plain(
                "the stable complete Tate resolution is supported in rank one",
            ))
        } else {
            Ok(())
        }
    }

    pub(crate) fn term(
        &self,
        model: &mut Model,
        name: Rc<str>,
        mode: Mode,
        monomial: Vec<u128>,
        exterior: Vec<usize>,
        period: Vec<i128>,
    ) -> Result<(Grade, State)> {
        self.validate_mode(mode)?;
        let count = self.powers.len();
        if monomial.len() != count || exterior.len() != count || period.len() != count {
            return Err(Error::plain(
                "Verbitsky Tate lists must each have one entry per factor",
            ));
        }
        let mut degree = 0i128;
        let mut shift = 0i128;
        let mut key = Vec::with_capacity(count);
        for (((exponent, exterior), periodic), &d) in monomial
            .into_iter()
            .zip(exterior)
            .zip(period)
            .zip(&self.powers)
        {
            if exterior > 1 {
                return Err(Error::plain(
                    "Verbitsky Tate exterior exponents must be zero or one",
                ));
            }
            if mode == Mode::Ordinary && periodic < 0 {
                return Err(Error::plain(
                    "ordinary Hochschild coordinates require nonnegative periodic exponents",
                ));
            }
            let value = periodic
                .checked_mul(2)
                .and_then(|v| v.checked_add(exterior as i128))
                .ok_or_else(grading_range)?;
            degree = degree.checked_add(value).ok_or_else(grading_range)?;
            let coordinate = i128::try_from(exponent).map_err(|_| grading_range())?;
            let weight = (d as i128)
                .checked_mul(periodic)
                .and_then(|v| coordinate.checked_sub(v))
                .ok_or_else(grading_range)?;
            shift = shift.checked_add(weight).ok_or_else(grading_range)?;
            key.push(Factor {
                exponent,
                exterior: exterior == 1,
                periodic,
            });
        }
        let mut state = State::default();
        if self.allowed(&key, mode) {
            state.insert(model, key, Rational::integer(1))?;
        }
        Ok((
            Grade {
                name,
                mode,
                degree,
                shift,
            },
            state,
        ))
    }

    fn allowed(&self, monomial: &[Factor], mode: Mode) -> bool {
        monomial.iter().zip(&self.powers).all(|(factor, &power)| {
            let bound = power as u128
                - u128::from(factor.exterior || factor.periodic != 0 || mode == Mode::Stable);
            factor.exponent < bound
        })
    }
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
struct Factor {
    exponent: u128,
    exterior: bool,
    periodic: i128,
}

type Monomial = Vec<Factor>;

#[derive(Default)]
pub(crate) struct State {
    terms: BTreeMap<Monomial, Rational>,
    cells: usize,
}

impl State {
    fn insert(&mut self, model: &mut Model, key: Monomial, coefficient: Rational) -> Result<()> {
        use std::collections::btree_map::Entry;
        if coefficient.is_zero() {
            return Ok(());
        }
        model.spend(key.len() + 1)?;
        let cells = 6 * key.len() + 2 * LIMBS + 12;
        match self.terms.entry(key) {
            Entry::Vacant(entry) => {
                storage(self.cells + cells)?;
                entry.insert(coefficient);
                self.cells += cells;
            }
            Entry::Occupied(mut entry) => {
                let sum = entry.get().add(model, coefficient)?;
                if sum.is_zero() {
                    entry.remove();
                    self.cells -= cells;
                } else {
                    entry.insert(sum);
                }
            }
        }
        Ok(())
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(self.cells + 1)?;
        Ok(Value::VerbitskyTate(Rc::new(self)))
    }

    pub(crate) fn add(&self, model: &mut Model, other: &Self) -> Result<Self> {
        let mut result = Self::default();
        result.add_scaled(model, self, Rational::integer(1))?;
        result.add_scaled(model, other, Rational::integer(1))?;
        Ok(result)
    }

    fn add_scaled(&mut self, model: &mut Model, other: &Self, scale: Rational) -> Result<()> {
        if !scale.is_zero() {
            for (key, &coefficient) in &other.terms {
                let value = coefficient.multiply(model, scale)?;
                self.insert(model, key.clone(), value)?;
            }
        }
        Ok(())
    }

    pub(crate) fn scale(&self, model: &mut Model, scalar: Entry) -> Result<Self> {
        let mut result = Self::default();
        let scalar = Rational::input(model, scalar)?;
        result.add_scaled(model, self, scalar)?;
        Ok(result)
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        for ((key, value), (other_key, other_value)) in self.terms.iter().zip(&other.terms) {
            model.spend(key.len() + 1)?;
            let order = key.cmp(other_key).then_with(|| value.compare(*other_value));
            if !order.is_eq() {
                return Ok(order);
            }
        }
        Ok(self.terms.len().cmp(&other.terms.len()))
    }

    pub(crate) fn text(&self, model: &mut Model) -> Result<String> {
        let mut text = "[".to_owned();
        for (i, (key, coefficient)) in self.terms.iter().enumerate() {
            if i > 0 {
                text.push_str(", ");
            }
            let factors = key
                .iter()
                .map(|value| (value.exponent, u8::from(value.exterior), value.periodic))
                .collect::<Vec<_>>();
            write!(text, "({factors:?}, {})", coefficient.text(model)?).unwrap();
            check_report(&text)?;
        }
        text.push(']');
        Ok(text)
    }

    pub(crate) fn hodge_project(&self, model: &mut Model, index: i128) -> Result<Self> {
        let mut result = Self::default();
        for (key, &coefficient) in &self.terms {
            let value = key.iter().try_fold(0i128, |sum, factor| {
                sum.checked_add(factor.periodic)
                    .and_then(|v| v.checked_add(i128::from(factor.exterior)))
                    .ok_or_else(grading_range)
            })?;
            if value == index {
                result.insert(model, key.clone(), coefficient)?;
            }
        }
        Ok(result)
    }
}

fn check_report(text: &str) -> Result<()> {
    if text.len() >= MAX_SHOWN_BYTES {
        Err(Error::plain(
            "Verbitsky Tate report exceeds the output limit",
        ))
    } else {
        Ok(())
    }
}

impl Definition {
    pub(crate) fn cup(
        &self,
        model: &mut Model,
        mode: Mode,
        left: &State,
        right: &State,
    ) -> Result<State> {
        let mut result = State::default();
        for (a, &coefficient) in &left.terms {
            for (b, &other) in &right.terms {
                model.spend(a.len() + 1)?;
                if a.iter().zip(b).any(|(a, b)| a.exterior && b.exterior) {
                    continue;
                }
                let mut odd = false;
                let mut right_parity = false;
                let mut output = Vec::with_capacity(a.len());
                for (a, b) in a.iter().zip(b) {
                    let periodic = a
                        .periodic
                        .checked_add(b.periodic)
                        .ok_or_else(grading_range)?;
                    output.push(Factor {
                        exponent: a.exponent + b.exponent,
                        exterior: a.exterior || b.exterior,
                        periodic,
                    });
                    odd ^= a.exterior && right_parity;
                    right_parity ^= b.exterior;
                }
                if self.allowed(&output, mode) {
                    let value = coefficient.multiply(model, other)?;
                    result.insert(model, output, if odd { value.negate() } else { value })?;
                }
            }
        }
        Ok(result)
    }

    pub(crate) fn bv(&self, model: &mut Model, value: &State) -> Result<State> {
        let mut output = State::default();
        for (key, &coefficient) in &value.terms {
            let mut odd = false;
            for (i, (factor, &d)) in key.iter().zip(&self.powers).enumerate() {
                if !factor.exterior {
                    continue;
                }
                let scalar = Rational::integer(factor.periodic)
                    .multiply(model, Rational::integer(d as i128))?;
                let scalar = scalar.add(
                    model,
                    Rational::integer(d as i128 - 1 - factor.exponent as i128),
                )?;
                let value = coefficient.multiply(model, scalar)?;
                let mut key = key.clone();
                key[i].exterior = false;
                output.insert(model, key, if odd { value.negate() } else { value })?;
                odd = !odd;
            }
        }
        Ok(output)
    }

    pub(crate) fn bracket(
        &self,
        model: &mut Model,
        mode: Mode,
        p: i128,
        q: i128,
        left: &State,
        right: &State,
    ) -> Result<State> {
        let product = self.cup(model, mode, left, right)?;
        let mut value = self.bv(model, &product)?;
        let dl = self.bv(model, left)?;
        let dr = self.bv(model, right)?;
        let first = self.cup(model, mode, &dl, right)?;
        let second = self.cup(model, mode, left, &dr)?;
        value.add_scaled(model, &first, Rational::integer(-1))?;
        value.add_scaled(
            model,
            &second,
            Rational::integer(if p.rem_euclid(2) == 0 { -1 } else { 1 }),
        )?;
        if p.rem_euclid(2) == 0 && q.rem_euclid(2) == 1 {
            Ok(value)
        } else {
            value.scale(model, (-1, 1))
        }
    }
}
