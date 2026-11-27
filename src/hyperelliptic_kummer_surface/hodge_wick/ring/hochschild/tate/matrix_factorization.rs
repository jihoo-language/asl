use std::cmp::Ordering;

use super::{
    check_report, storage, BTreeMap, Definition, Entry, Error, Label, Mode, Model, ParsedExpr,
    Rational, Rc, Result, State as TateState, Value, Write, LIMBS,
};

mod cone;
pub(crate) mod defect;
pub(crate) mod knorrer;

pub(crate) enum Command {
    Declare(Label, Label, Vec<usize>),
    Shift(Label, Label),
    Cone(Option<Label>, Box<ParsedExpr>),
    Object(Label),
    Hom(Label, Label, bool),
    Trace(Box<ParsedExpr>),
    Pair(Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Differential,
    Class,
    Representative,
    Primitive,
    Scale(Entry),
}

pub(crate) enum Expression {
    Term(Label, Label, bool, usize, usize, usize, u128, Entry),
    Basis(Label, Label, bool, usize, usize, u128),
    Identity(Label),
    Unary(Unary, Box<ParsedExpr>),
    Compose(Box<[ParsedExpr; 2]>),
    StableAction(Box<[ParsedExpr; 2]>),
    Shift(Label, Label, Box<ParsedExpr>),
}

pub(crate) struct Object {
    pub(crate) mark: Rc<str>,
    d: usize,
    factors: Vec<usize>,
}

impl Object {
    pub(crate) fn new(mark: Rc<str>, definition: &Definition, factors: Vec<usize>) -> Result<Self> {
        let [d] = definition.powers.as_slice() else {
            return Err(Error::plain(
                "matrix factorization requires a rank-one Tate marking",
            ));
        };
        if factors.iter().any(|&i| i == 0 || i >= *d) {
            return Err(Error::plain(
                "matrix factorization index must lie between 1 and d−1",
            ));
        }
        storage(factors.len().saturating_add(12))?;
        Ok(Self {
            mark,
            d: *d,
            factors,
        })
    }

    pub(crate) fn shifted(&self) -> Self {
        Self {
            mark: self.mark.clone(),
            d: self.d,
            factors: self.factors.iter().map(|i| self.d - i).collect(),
        }
    }

    pub(crate) fn is_shift_of(&self, other: &Self) -> bool {
        self.mark == other.mark
            && self.factors.len() == other.factors.len()
            && self
                .factors
                .iter()
                .zip(&other.factors)
                .all(|(&i, &j)| i == self.d - j)
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.factors.is_empty()
    }

    pub(crate) fn cells(&self) -> usize {
        self.factors.len() + 12
    }

    pub(crate) fn text(&self) -> Result<String> {
        let text = format!(
            "W=x^{} Tate={} factors={:?}",
            self.d, self.mark, self.factors
        );
        check_report(&text)?;
        Ok(text)
    }
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) source: Rc<str>,
    pub(crate) target: Rc<str>,
    pub(crate) odd: bool,
}

impl Grade {
    pub(crate) fn objects(&self, model: &Model) -> Result<[Rc<Object>; 2]> {
        Ok([
            model.matrix_factorization(&self.source)?,
            model.matrix_factorization(&self.target)?,
        ])
    }

    pub(crate) fn opposite(&self) -> Self {
        Self {
            odd: !self.odd,
            ..self.clone()
        }
    }

    pub(crate) fn compose(&self, right: &Self) -> Result<Self> {
        if self.source != right.target {
            return Err(Error::plain(
                "matrix factorization composition requires matching middle objects",
            ));
        }
        Ok(Self {
            source: right.source.clone(),
            target: self.target.clone(),
            odd: self.odd ^ right.odd,
        })
    }
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
struct Key {
    row: usize,
    column: usize,
    slot: bool,
    exponent: u128,
}

#[derive(Default)]
pub(crate) struct Morphism {
    terms: BTreeMap<Key, Rational>,
    cells: usize,
}

impl Morphism {
    fn insert(&mut self, model: &mut Model, key: Key, value: Rational) -> Result<()> {
        if value.is_zero() {
            return Ok(());
        }
        model.spend(1)?;
        let cells = 2 * LIMBS + 12;
        use std::collections::btree_map::Entry;
        match self.terms.entry(key) {
            Entry::Vacant(entry) => {
                storage(self.cells + cells)?;
                entry.insert(value);
                self.cells += cells;
            }
            Entry::Occupied(mut entry) => {
                let sum = entry.get().add(model, value)?;
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
        Ok(Value::MatrixFactorization(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        for ((key, value), (right, v)) in self.terms.iter().zip(&other.terms) {
            model.spend(1)?;
            let order = key.cmp(right).then_with(|| value.compare(*v));
            if !order.is_eq() {
                return Ok(order);
            }
        }
        Ok(self.terms.len().cmp(&other.terms.len()))
    }

    pub(crate) fn text(&self, model: &mut Model) -> Result<String> {
        let mut text = "[".to_owned();
        let mut sep = "";
        for (key, value) in &self.terms {
            write!(
                text,
                "{sep}({}, {}, {}, {}, {})",
                key.row,
                key.column,
                u8::from(key.slot),
                key.exponent,
                value.text(model)?
            )
            .unwrap();
            check_report(&text)?;
            sep = ", ";
        }
        text.push(']');
        Ok(text)
    }

    pub(crate) fn add(&self, model: &mut Model, other: &Self) -> Result<Self> {
        let mut result = Self::default();
        for (&key, &coefficient) in self.terms.iter().chain(&other.terms) {
            result.insert(model, key, coefficient)?;
        }
        Ok(result)
    }

    pub(crate) fn scale(&self, model: &mut Model, scale: Entry) -> Result<Self> {
        let scale = Rational::input(model, scale)?;
        let mut result = Self::default();
        for (&key, &value) in &self.terms {
            let scaled = scale.multiply(model, value)?;
            result.insert(model, key, scaled)?;
        }
        Ok(result)
    }

    pub(crate) fn shifted(&self, model: &mut Model) -> Result<Self> {
        let mut result = Self::default();
        for (&key, &value) in &self.terms {
            result.insert(
                model,
                Key {
                    slot: !key.slot,
                    ..key
                },
                value,
            )?;
        }
        Ok(result)
    }

    pub(crate) fn compose(&self, model: &mut Model, odd: bool, right: &Self) -> Result<Self> {
        let mut output = Self::default();
        for (&a, &f) in &self.terms {
            let start = Key {
                row: a.column,
                column: 0,
                slot: false,
                exponent: 0,
            };
            for (&b, &g) in right
                .terms
                .range(start..)
                .take_while(|(b, _)| b.row == a.column)
            {
                model.spend(1)?;
                if a.slot ^ odd != b.slot {
                    continue;
                }
                let exponent = add_exponents(a.exponent, b.exponent)?;
                let coefficient = f.multiply(model, g)?;
                output.insert(
                    model,
                    Key {
                        row: a.row,
                        column: b.column,
                        slot: a.slot,
                        exponent,
                    },
                    coefficient,
                )?;
            }
        }
        Ok(output)
    }
}

fn add_exponents(a: u128, b: u128) -> Result<u128> {
    a.checked_add(b).ok_or_else(|| {
        Error::plain("matrix factorization polynomial exponent exceeds the unsigned integer range")
    })
}

pub(crate) struct Hom<'a> {
    source: &'a Object,
    target: &'a Object,
    odd: bool,
}

impl<'a> Hom<'a> {
    pub(crate) fn new(source: &'a Object, target: &'a Object, odd: bool) -> Result<Self> {
        if source.mark != target.mark {
            return Err(Error::plain(
                "matrix factorization Hom requires one common Tate marking",
            ));
        }
        Ok(Self {
            source,
            target,
            odd,
        })
    }

    fn indices(&self, row: usize, column: usize) -> Result<(usize, usize)> {
        self.source
            .factors
            .get(column)
            .zip(self.target.factors.get(row))
            .map(|(&i, &j)| (i, j))
            .ok_or_else(|| Error::plain("matrix factorization Hom row or column is out of range"))
    }

    fn offsets(&self, i: usize, j: usize) -> [u128; 2] {
        if self.odd {
            let sum = i as u128 + j as u128;
            [
                (self.source.d as u128).saturating_sub(sum),
                sum.saturating_sub(self.source.d as u128),
            ]
        } else {
            [i.saturating_sub(j) as u128, j.saturating_sub(i) as u128]
        }
    }

    fn dimension(&self, i: usize, j: usize) -> u128 {
        i.min(j).min(self.source.d - i).min(self.source.d - j) as u128
    }

    pub(crate) fn term(
        &self,
        model: &mut Model,
        row: usize,
        column: usize,
        slot: usize,
        exponent: u128,
        coefficient: Entry,
    ) -> Result<Morphism> {
        self.indices(row, column)?;
        if slot > 1 {
            return Err(Error::plain(
                "matrix factorization block slot must be zero or one",
            ));
        }
        let value = Rational::input(model, coefficient)?;
        let mut output = Morphism::default();
        output.insert(
            model,
            Key {
                row,
                column,
                slot: slot == 1,
                exponent,
            },
            value,
        )?;
        Ok(output)
    }

    fn add_basis(
        &self,
        model: &mut Model,
        output: &mut Morphism,
        key: Key,
        power: u128,
        value: Rational,
    ) -> Result<()> {
        let (i, j) = self.indices(key.row, key.column)?;
        for (slot, offset) in self.offsets(i, j).into_iter().enumerate() {
            let exponent = add_exponents(offset, power)?;
            output.insert(
                model,
                Key {
                    slot: slot == 1,
                    exponent,
                    ..key
                },
                if slot == 1 && self.odd {
                    value.negate()
                } else {
                    value
                },
            )?;
        }
        Ok(())
    }

    pub(crate) fn basis(
        &self,
        model: &mut Model,
        row: usize,
        column: usize,
        power: u128,
    ) -> Result<Morphism> {
        let (i, j) = self.indices(row, column)?;
        if power >= self.dimension(i, j) {
            return Err(Error::plain(
                "matrix factorization Hom basis index is out of range",
            ));
        }
        let mut output = Morphism::default();
        let key = Key {
            row,
            column,
            slot: false,
            exponent: 0,
        };
        self.add_basis(model, &mut output, key, power, Rational::integer(1))?;
        Ok(output)
    }

    pub(crate) fn identity(&self, model: &mut Model) -> Result<Morphism> {
        let mut result = Morphism::default();
        for index in 0..self.source.factors.len() {
            let key = Key {
                row: index,
                column: index,
                slot: false,
                exponent: 0,
            };
            self.add_basis(model, &mut result, key, 0, Rational::integer(1))?;
        }
        Ok(result)
    }

    pub(crate) fn differential(&self, model: &mut Model, value: &Morphism) -> Result<Morphism> {
        let mut output = Morphism::default();
        let d = self.source.d;
        for (&key, &coefficient) in &value.terms {
            let (i, j) = self.indices(key.row, key.column)?;
            let left_shift = if key.slot { d - j } else { j };
            let right_shift = if key.slot ^ self.odd { i } else { d - i };
            let exponent = add_exponents(key.exponent, left_shift as u128)?;
            output.insert(
                model,
                Key {
                    slot: !key.slot,
                    exponent,
                    ..key
                },
                coefficient,
            )?;
            let exponent = add_exponents(key.exponent, right_shift as u128)?;
            output.insert(
                model,
                Key { exponent, ..key },
                if self.odd {
                    coefficient
                } else {
                    coefficient.negate()
                },
            )?;
        }
        Ok(output)
    }

    pub(crate) fn require_cycle(&self, model: &mut Model, value: &Morphism) -> Result<()> {
        if !self.differential(model, value)?.terms.is_empty() {
            return Err(Error::plain("matrix factorization cochain is not closed"));
        }
        Ok(())
    }

    pub(crate) fn class(&self, model: &mut Model, value: &Morphism) -> Result<Morphism> {
        self.require_cycle(model, value)?;
        let mut result = Morphism::default();
        for (&key, &coefficient) in value.terms.iter().filter(|(key, _)| !key.slot) {
            let (i, j) = self.indices(key.row, key.column)?;
            let power = self.reduce_power(key, i, j)?;
            if power < self.dimension(i, j) {
                self.add_basis(model, &mut result, key, power, coefficient)?;
            }
        }
        Ok(result)
    }

    fn reduce_power(&self, key: Key, i: usize, j: usize) -> Result<u128> {
        key.exponent
            .checked_sub(self.offsets(i, j)[0])
            .ok_or_else(|| {
                Error::plain("closed matrix factorization block has invalid leading exponent")
            })
    }

    pub(crate) fn primitive(&self, model: &mut Model, value: &Morphism) -> Result<Morphism> {
        self.require_cycle(model, value)?;
        let mut result = Morphism::default();
        let d = self.source.d;
        for (&key, &coefficient) in value.terms.iter().filter(|(key, _)| !key.slot) {
            let (i, j) = self.indices(key.row, key.column)?;
            let power = self.reduce_power(key, i, j)?;
            let m = self.dimension(i, j);
            let exponent = power
                .checked_sub(m)
                .ok_or_else(|| Error::plain("matrix factorization cycle is not a boundary"))?;
            // Use the term of smallest exponent in the top block of δ.
            let (slot, coefficient) = if self.odd {
                if j >= i {
                    (true, coefficient)
                } else {
                    (false, coefficient.negate())
                }
            } else {
                (d - j < i, coefficient)
            };
            result.insert(
                model,
                Key {
                    slot,
                    exponent,
                    ..key
                },
                coefficient,
            )?;
        }
        let opposite = Self {
            odd: !self.odd,
            ..*self
        };
        let check = opposite.differential(model, &result)?;
        if !check.compare(model, value)?.is_eq() {
            return Err(Error::plain(
                "matrix factorization boundary homotopy failed verification",
            ));
        }
        Ok(result)
    }

    pub(crate) fn report(&self, model: &mut Model) -> Result<String> {
        let mut dimension = 0u128;
        let mut text = " blocks=[".to_owned();
        let mut sep = "";
        for (r, &j) in self.target.factors.iter().enumerate() {
            for (s, &i) in self.source.factors.iter().enumerate() {
                model.spend(1)?;
                let m = self.dimension(i, j);
                dimension = dimension.checked_add(m).ok_or_else(|| {
                    Error::plain(
                        "matrix factorization Hom dimension exceeds the supported integer range",
                    )
                })?;
                let offset = self.offsets(i, j)[0] as i128;
                let start = 2 * offset
                    + if self.odd {
                        2 * (i as i128) - self.source.d as i128
                    } else {
                        0
                    };
                write!(text, "{sep}({r}, {s}, first_degree={start}, count={m})").unwrap();
                check_report(&text)?;
                sep = ", ";
            }
        }
        text.push(']');
        Ok(format!("dimension={dimension}{text}"))
    }

    pub(crate) fn trace(&self, model: &mut Model, value: &Morphism) -> Result<String> {
        self.require_cycle(model, value)?;
        let mut sum = Rational::zero();
        if self.odd {
            let d = self.source.d;
            for (&key, &coefficient) in &value.terms {
                if key.column != key.row {
                    continue;
                }
                let i = self.source.factors[key.column];
                let derivative = if key.slot { d - i } else { i };
                if key.exponent == (d - 1 - derivative) as u128 {
                    let scalar = Rational::integer(if key.slot {
                        -(derivative as i128)
                    } else {
                        derivative as i128
                    });
                    let term = coefficient.multiply(model, scalar)?;
                    sum = sum.add(model, term)?;
                }
            }
            let inverse = Rational::integer(d as i128).inverse()?;
            sum = sum.multiply(model, inverse)?;
        }
        sum.text(model)
    }

    pub(crate) fn stable_action(
        &self,
        model: &mut Model,
        value: &TateState,
        morphism: &Morphism,
    ) -> Result<Morphism> {
        let mut raw = Morphism::default();
        for (tate, &coefficient) in &value.terms {
            let [factor] = tate.as_slice() else {
                unreachable!()
            };
            debug_assert!(!factor.exterior);
            for (&key, &entry) in &morphism.terms {
                let exponent = add_exponents(key.exponent, factor.exponent)?;
                let value = coefficient.multiply(model, entry)?;
                raw.insert(model, Key { exponent, ..key }, value)?;
            }
        }
        self.class(model, &raw)
    }
}

pub(crate) fn check_stable(mark: &super::Grade, source: &Object) -> Result<()> {
    if mark.mode != Mode::Stable || mark.name != source.mark || mark.degree.rem_euclid(2) != 0 {
        return Err(Error::plain(
            "matrix factorization action requires an even stable Tate class of the same marking",
        ));
    }
    Ok(())
}
