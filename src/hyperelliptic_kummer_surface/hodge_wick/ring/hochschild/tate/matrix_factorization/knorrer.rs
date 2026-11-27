use std::cmp::Ordering;
use std::collections::BTreeSet;

use super::{
    add_exponents, check_report, storage, BTreeMap, Entry, Error, Grade as BaseGrade,
    Hom as BaseHom, Key as BaseKey, Label, Model, Morphism as BaseMorphism, Object, ParsedExpr,
    Rational, Rc, Result, Value, Write, LIMBS,
};

mod equivalence;
pub(crate) use equivalence::{include, project};

pub(crate) enum Command {
    Declare(Label, Label),
    Object(Label),
    Inspect(Query, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Weights,
    Cone,
    Trace,
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Diff,
    Homotopy,
    Primitive,
    Class,
    Representative,
    Project,
    Scale(Entry),
}

pub(crate) enum Expression {
    Term(Label, Label, bool, usize, usize, [u128; 3], Entry),
    Include(Label, Label, Box<ParsedExpr>),
    Shift(Label, Label, Box<ParsedExpr>),
    Compose(Box<[ParsedExpr; 2]>),
    Apply(Unary, Box<ParsedExpr>),
}

pub(crate) struct Doubled {
    pub(crate) base: Rc<str>,
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) source: Rc<str>,
    pub(crate) target: Rc<str>,
    pub(crate) odd: bool,
}

impl Grade {
    pub(crate) fn opposite(&self) -> Self {
        Self {
            odd: !self.odd,
            ..self.clone()
        }
    }

    pub(crate) fn base_grade(&self, model: &Model) -> Result<BaseGrade> {
        let source = model.knorrer_object(&self.source)?;
        let target = model.knorrer_object(&self.target)?;
        Ok(BaseGrade {
            source: source.base.clone(),
            target: target.base.clone(),
            odd: self.odd,
        })
    }

    pub(crate) fn compose(&self, right: &Self) -> Result<Self> {
        if self.source != right.target {
            return Err(Error::plain(
                "Knörrer composition requires matching middle objects",
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
    p_out: bool,
    q_out: bool,
    q_odd: bool,
    exponents: [u128; 3],
}

impl Key {
    fn p_in(self, odd: bool) -> bool {
        self.p_out ^ odd ^ self.q_odd
    }

    fn q_in(self) -> bool {
        self.q_out ^ self.q_odd
    }

    fn coordinates(self, odd: bool) -> (usize, usize) {
        (
            4 * self.row + 2 * usize::from(self.p_out) + usize::from(self.q_out),
            4 * self.column + 2 * usize::from(self.p_in(odd)) + usize::from(self.q_in()),
        )
    }

    fn physical_sign(self, odd: bool, value: Rational) -> Rational {
        if self.q_odd && self.p_in(odd) {
            value.negate()
        } else {
            value
        }
    }
}

#[derive(Default)]
pub(crate) struct Cochain {
    terms: BTreeMap<Key, Rational>,
    cells: usize,
}

impl Cochain {
    fn insert(&mut self, model: &mut Model, key: Key, value: Rational) -> Result<()> {
        if value.is_zero() {
            return Ok(());
        }
        model.spend(1)?;
        let cells = 2 * LIMBS + 16;
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
        Ok(Value::Knorrer(Rc::new(self)))
    }

    pub(crate) fn add(&self, model: &mut Model, other: &Self) -> Result<Self> {
        let mut value = Self::default();
        for (&key, &coefficient) in self.terms.iter().chain(&other.terms) {
            value.insert(model, key, coefficient)?;
        }
        Ok(value)
    }

    pub(crate) fn scale(&self, model: &mut Model, coefficient: Entry) -> Result<Self> {
        let scale = Rational::input(model, coefficient)?;
        let mut result = Self::default();
        for (&key, &coefficient) in &self.terms {
            let value = coefficient.multiply(model, scale)?;
            result.insert(model, key, value)?;
        }
        Ok(result)
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        for ((key, value), (right, coefficient)) in self.terms.iter().zip(&other.terms) {
            model.spend(1)?;
            let order = key.cmp(right).then_with(|| value.compare(*coefficient));
            if !order.is_eq() {
                return Ok(order);
            }
        }
        Ok(self.terms.len().cmp(&other.terms.len()))
    }

    pub(crate) fn text(&self, model: &mut Model, odd: bool) -> Result<String> {
        let mut report = "[".to_owned();
        let mut separator = "";
        for (&key, &coefficient) in &self.terms {
            let (row, column) = key.coordinates(odd);
            let value = key.physical_sign(odd, coefficient).text(model)?;
            write!(
                report,
                "{separator}({row}, {column}, {:?}, {value})",
                key.exponents
            )
            .unwrap();
            check_report(&report)?;
            separator = ", ";
        }
        report.push(']');
        Ok(report)
    }

    // Swapping P parity changes the physical sign of a Q-odd tensor term.
    // Conjugating by the Q parity sign compensates for that change.
    pub(crate) fn shifted(&self, model: &mut Model) -> Result<Self> {
        let mut output = Self::default();
        for (&key, &value) in &self.terms {
            output.insert(
                model,
                Key {
                    p_out: !key.p_out,
                    ..key
                },
                value,
            )?;
        }
        Ok(output)
    }

    /// Composition in the graded tensor basis, with the right map first.
    pub(crate) fn compose(
        &self,
        model: &mut Model,
        odd: bool,
        right: &Self,
        right_odd: bool,
    ) -> Result<Self> {
        let mut output = Self::default();
        for (&left, &f) in &self.terms {
            let start = Key {
                row: left.column,
                column: 0,
                p_out: false,
                q_out: false,
                q_odd: false,
                exponents: [0; 3],
            };
            for (&right, &g) in right
                .terms
                .range(start..)
                .take_while(|(key, _)| key.row == left.column)
            {
                model.spend(1)?;
                if left.p_in(odd) != right.p_out || left.q_in() != right.q_out {
                    continue;
                }
                let [a, b, c] = left.exponents;
                let [x, u, v] = right.exponents;
                let exponents = [
                    add_exponents(a, x)?,
                    add_exponents(b, u)?,
                    add_exponents(c, v)?,
                ];
                let negative = left.q_odd && (right_odd ^ right.q_odd);
                let coefficient = f.multiply(model, g)?;
                let key = Key {
                    column: right.column,
                    q_odd: left.q_odd ^ right.q_odd,
                    exponents,
                    ..left
                };
                output.insert(
                    model,
                    key,
                    if negative {
                        coefficient.negate()
                    } else {
                        coefficient
                    },
                )?;
            }
        }
        Ok(output)
    }
}

pub(crate) struct Hom<'a> {
    source: &'a Object,
    target: &'a Object,
    odd: bool,
}

impl<'a> Hom<'a> {
    pub(crate) fn new(source: &'a Object, target: &'a Object, odd: bool) -> Result<Self> {
        BaseHom::new(source, target, odd)?;
        Ok(Self {
            source,
            target,
            odd,
        })
    }

    fn base(&self) -> Result<BaseHom<'_>> {
        BaseHom::new(self.source, self.target, self.odd)
    }

    pub(crate) fn term(
        &self,
        model: &mut Model,
        row: usize,
        column: usize,
        exponents: [u128; 3],
        coefficient: Entry,
    ) -> Result<Cochain> {
        if row / 4 >= self.target.factors.len() || column / 4 >= self.source.factors.len() {
            return Err(Error::plain("Knörrer matrix row or column is out of range"));
        }
        let p_out = row % 4 >= 2;
        let p_in = column % 4 >= 2;
        let q_out = row % 2 == 1;
        let q_in = column % 2 == 1;
        let q_odd = q_out ^ q_in;
        if p_out ^ p_in ^ q_odd != self.odd {
            return Err(Error::plain("Knörrer term has the wrong matrix parity"));
        }
        let key = Key {
            row: row / 4,
            column: column / 4,
            p_out,
            q_out,
            q_odd,
            exponents,
        };
        let coefficient = key.physical_sign(self.odd, Rational::input(model, coefficient)?);
        let mut result = Cochain::default();
        result.insert(model, key, coefficient)?;
        Ok(result)
    }

    pub(crate) fn differential(&self, model: &mut Model, value: &Cochain) -> Result<Cochain> {
        if value.terms.is_empty() {
            return Ok(Cochain::default());
        }
        let left = differential_matrix(model, self.target)?;
        let right = differential_matrix(model, self.source)?;
        let mut output = left.compose(model, true, value, self.odd)?;
        let correction = value.compose(model, self.odd, &right, true)?;
        for (key, coefficient) in correction.terms {
            output.insert(
                model,
                key,
                if self.odd {
                    coefficient
                } else {
                    coefficient.negate()
                },
            )?;
        }
        Ok(output)
    }

    fn require_cycle(&self, model: &mut Model, value: &Cochain) -> Result<()> {
        if !self.differential(model, value)?.terms.is_empty() {
            return Err(Error::plain("Knörrer cochain is not closed"));
        }
        Ok(())
    }

    pub(crate) fn weights(&self, model: &mut Model, value: &Cochain) -> Result<String> {
        let d = self.source.d as i128;
        let mut weights = BTreeSet::new();
        for &key in value.terms.keys() {
            let [a, b, c] = key.exponents;
            let [Ok(a), Ok(b), Ok(c)] = [a, b, c].map(i128::try_from) else {
                return Err(weight_range());
            };
            let i = self.source.factors[key.column] as i128;
            let j = self.target.factors[key.row] as i128;
            let target = if key.p_out { d - 2 * j } else { 0 };
            let source = if key.p_in(self.odd) { d - 2 * i } else { 0 };
            let first = b
                .checked_add(c)
                .and_then(|sum| sum.checked_mul(d))
                .and_then(|q| a.checked_mul(2).and_then(|x| x.checked_add(q)))
                .and_then(|r| r.checked_add(target - source))
                .ok_or_else(weight_range)?;
            let second = (b - c)
                .checked_add(i128::from(key.q_in()) - i128::from(key.q_out))
                .ok_or_else(weight_range)?;
            model.spend(1)?;
            weights.insert((first, second));
        }
        let report = format!("{weights:?}");
        check_report(&report)?;
        Ok(report)
    }
}

fn weight_range() -> Error {
    Error::plain("Knörrer homogeneous weight exceeds the signed integer range")
}

fn differential_matrix(model: &mut Model, object: &Object) -> Result<Cochain> {
    let mut result = Cochain::default();
    for (i, &factor) in object.factors.iter().enumerate() {
        for p in [false, true] {
            for q in [false, true] {
                let exponent = if p { factor } else { object.d - factor } as u128;
                let key = Key {
                    row: i,
                    column: i,
                    p_out: p,
                    q_out: q,
                    q_odd: false,
                    exponents: [exponent, 0, 0],
                };
                result.insert(model, key, Rational::integer(1))?;
                let exponents = if q { [0, 1, 0] } else { [0, 0, 1] };
                result.insert(
                    model,
                    Key {
                        q_odd: true,
                        exponents,
                        ..key
                    },
                    Rational::integer(1),
                )?;
            }
        }
    }
    Ok(result)
}

pub(crate) fn object_report(model: &mut Model, object: &Object) -> Result<String> {
    let differential = differential_matrix(model, object)?;
    let square = differential.compose(model, true, &differential, true)?;
    let mut potential = Cochain::default();
    for i in 0..object.factors.len() {
        for p in [false, true] {
            for q in [false, true] {
                let key = Key {
                    row: i,
                    column: i,
                    p_out: p,
                    q_out: q,
                    q_odd: false,
                    exponents: [object.d as u128, 0, 0],
                };
                potential.insert(model, key, Rational::integer(1))?;
                potential.insert(
                    model,
                    Key {
                        exponents: [0, 1, 1],
                        ..key
                    },
                    Rational::integer(1),
                )?;
            }
        }
    }
    if !square.compare(model, &potential)?.is_eq() {
        return Err(Error::plain(
            "Knörrer differential failed to square to x^d+uv",
        ));
    }
    let report = format!(
        "W=x^{}+uv factors={:?} differential={} square=verified",
        object.d,
        object.factors,
        differential.text(model, true)?
    );
    check_report(&report)?;
    Ok(report)
}
