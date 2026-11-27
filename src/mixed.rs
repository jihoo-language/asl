use std::borrow::Cow;
use std::collections::{BTreeMap, HashMap};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::odd_ext;
use crate::prime_linear::{add_scaled, Coordinate};
use crate::syntax::Label;

pub(crate) type Terms = Vec<(Label, i64)>;
pub(crate) type Rules = Vec<(Label, Terms)>;
pub(crate) type HigherRules = Vec<(Label, Vec<(usize, Label, i64)>)>;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Id(pub(crate) usize);

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) complex: Id,
    pub(crate) b: i64,
}

impl Grade {
    pub(crate) fn shift(self, offset: i128) -> Result<Self> {
        let b = i64::try_from(i128::from(self.b) + offset)
            .map_err(|_| Error::plain("mixed complex degree is out of range"))?;
        Ok(Self { b, ..self })
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Differential {
    Down,
    Up,
    Total,
}

impl Differential {
    pub(crate) fn degree(self) -> i128 {
        if matches!(self, Self::Up) {
            1
        } else {
            -1
        }
    }
}

pub(crate) struct Complex {
    pub(crate) prime: Prime,
    names: Vec<String>,
    degrees: Vec<i64>,
    indices: HashMap<String, usize>,
    down: Vec<Vec<Coordinate>>,
    up: Vec<Vec<Coordinate>>,
    higher: Option<Vec<Vec<Coordinate>>>,
}

impl Complex {
    pub(crate) fn build(
        model: &mut Model,
        field: usize,
        generators: Vec<(Label, i64)>,
        down: Rules,
        up: Rules,
    ) -> Result<Self> {
        let mut complex = Self::declared(model, field, generators)?;
        complex.down = complex.columns(model, down, Differential::Down)?;
        complex.up = complex.columns(model, up, Differential::Up)?;
        complex.finish(model)
    }

    fn declared(model: &mut Model, field: usize, generators: Vec<(Label, i64)>) -> Result<Self> {
        let prime = match field {
            2 => Prime::validated(2),
            3.. if field % 2 == 1 => Prime::new(model, field)?,
            _ => return Err(Error::plain("mixedcomplex requires a prime field")),
        };
        let mut complex = Self {
            prime,
            names: Vec::new(),
            degrees: Vec::new(),
            indices: HashMap::new(),
            down: Vec::new(),
            up: Vec::new(),
            higher: None,
        };
        for (name, degree) in generators {
            model.spend(1 + name.text.len())?;
            let index = complex.names.len();
            if complex.indices.insert(name.text.clone(), index).is_some() {
                return Err(Error::new(name.span, "duplicate mixed complex generator"));
            }
            complex.names.push(name.text);
            complex.degrees.push(degree);
        }
        Ok(complex)
    }

    pub(crate) fn build_higher(
        model: &mut Model,
        field: usize,
        generators: Vec<(Label, i64)>,
        rules: HigherRules,
    ) -> Result<Self> {
        let complex = Self::declared(model, field, generators)?;
        model.spend(complex.len())?;
        let mut columns = vec![None; complex.len()];
        for (source, entries) in rules {
            let index = complex.index(&source)?;
            if columns[index].is_some() {
                return Err(Error::new(
                    source.span,
                    "duplicate higher mixed differential column",
                ));
            }
            let mut column = BTreeMap::new();
            for (power, target, coefficient) in entries {
                model.spend(1 + target.text.len())?;
                let i = complex.index(&target)?;
                if complex.delay(index, i) != 2 * power as i128 {
                    return Err(Error::new(
                        target.span,
                        "higher mixed term has the wrong homological degree for its power",
                    ));
                }
                let coefficient =
                    i128::from(coefficient).rem_euclid(i128::from(complex.prime.get())) as u32;
                let current = column.entry(i).or_insert(0);
                *current = complex.prime.add(*current, coefficient);
            }
            columns[index] = Some(column.into_iter().filter(|&(_, c)| c != 0).collect());
        }
        complex.finish_total(
            model,
            columns.into_iter().map(Option::unwrap_or_default).collect(),
        )
    }

    pub(crate) fn from_total(
        model: &mut Model,
        prime: Prime,
        degrees: Vec<i64>,
        columns: Vec<Vec<Coordinate>>,
    ) -> Result<Self> {
        Self::numbered(model, prime, degrees)?.finish_total(model, columns)
    }

    fn finish_total(mut self, model: &mut Model, columns: Vec<Vec<Coordinate>>) -> Result<Self> {
        if columns.len() != self.len() {
            return Err(Error::plain(
                "higher mixed differential has the wrong number of columns",
            ));
        }
        let mut extras = Vec::with_capacity(self.len());
        let mut found = false;
        for (source, column) in columns.into_iter().enumerate() {
            model.spend(1 + column.len())?;
            let mut down = Vec::new();
            let mut up = Vec::new();
            let mut higher = Vec::new();
            for term in column {
                if term.0 >= self.len() {
                    return Err(Error::plain(
                        "higher mixed differential target is out of range",
                    ));
                }
                match self.delay(source, term.0) {
                    0 => down.push(term),
                    2 => up.push(term),
                    delay if delay > 2 && delay % 2 == 0 => higher.push(term),
                    _ => {
                        return Err(Error::plain(
                            "higher mixed differential has a negative or nonintegral power of u",
                        ))
                    }
                }
            }
            found |= !higher.is_empty();
            self.down.push(down);
            self.up.push(up);
            extras.push(higher);
        }
        if found {
            self.higher = Some(extras);
        }
        self.finish(model)
    }

    fn delay(&self, source: usize, target: usize) -> i128 {
        i128::from(self.degrees[target]) - i128::from(self.degrees[source]) + 1
    }

    pub(crate) fn from_columns(
        model: &mut Model,
        prime: Prime,
        degrees: Vec<i64>,
        down: Vec<Vec<Coordinate>>,
        up: Vec<Vec<Coordinate>>,
    ) -> Result<Self> {
        if down.len() != degrees.len() || up.len() != degrees.len() {
            return Err(Error::plain(
                "imported mixed differential has the wrong number of columns",
            ));
        }
        for (columns, op) in [(&down, Differential::Down), (&up, Differential::Up)] {
            for (source, column) in columns.iter().enumerate() {
                let expected = i128::from(degrees[source]) + op.degree();
                model.spend(1 + column.len())?;
                if column.iter().any(|&(target, _)| {
                    degrees.get(target).map(|&d| i128::from(d)) != Some(expected)
                }) {
                    return Err(Error::plain(
                        "imported mixed differential has the wrong homological degree",
                    ));
                }
            }
        }
        let mut complex = Self::numbered(model, prime, degrees)?;
        complex.down = down;
        complex.up = up;
        complex.finish(model)
    }

    fn numbered(model: &mut Model, prime: Prime, degrees: Vec<i64>) -> Result<Self> {
        let mut names = Vec::with_capacity(degrees.len());
        let mut indices = HashMap::with_capacity(degrees.len());
        for index in 0..degrees.len() {
            let name = format!("w{index}");
            model.spend(1 + name.len())?;
            indices.insert(name.clone(), index);
            names.push(name);
        }
        Ok(Self {
            prime,
            names,
            degrees,
            indices,
            down: Vec::new(),
            up: Vec::new(),
            higher: None,
        })
    }

    fn finish(self, model: &mut Model) -> Result<Self> {
        self.validate(model)?;
        model.allocate(
            1 + 6 * self.len()
                + self
                    .down
                    .iter()
                    .chain(&self.up)
                    .map(|column| 2 * column.len())
                    .sum::<usize>()
                + 2 * self.names.iter().map(String::len).sum::<usize>()
                + self.higher.as_ref().map_or(0, |higher| {
                    3 * higher.len() + higher.iter().map(|col| 2 * col.len()).sum::<usize>()
                }),
        )?;
        Ok(self)
    }

    pub(crate) fn len(&self) -> usize {
        self.names.len()
    }

    pub(crate) fn column(&self, operation: Differential, index: usize) -> &[Coordinate] {
        match operation {
            Differential::Down => &self.down[index],
            Differential::Up => &self.up[index],
            Differential::Total => unreachable!(),
        }
    }

    pub(crate) fn is_strict(&self) -> bool {
        self.higher.is_none()
    }

    pub(crate) fn total_column(&self, model: &mut Model, index: usize) -> Result<Vec<Coordinate>> {
        self.apply(model, Differential::Total, &[(index, 1)])
    }

    pub(crate) fn operator_column(
        &self,
        model: &mut Model,
        op: Differential,
        index: usize,
    ) -> Result<Cow<'_, [Coordinate]>> {
        if matches!(op, Differential::Total) {
            Ok(Cow::Owned(self.total_column(model, index)?))
        } else {
            Ok(Cow::Borrowed(self.column(op, index)))
        }
    }

    pub(crate) fn delta(
        &self,
        model: &mut Model,
        power: usize,
        input: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        match power {
            0 => return self.apply(model, Differential::Down, input),
            1 => return self.apply(model, Differential::Up, input),
            _ => {}
        }
        let Some(columns) = &self.higher else {
            return Ok(Vec::new());
        };
        let mut result = Vec::new();
        for &(i, coefficient) in input {
            let column = &columns[i];
            model.spend(1 + column.len())?;
            let component: Vec<_> = column
                .iter()
                .copied()
                .filter(|&(j, _)| self.delay(i, j) == 2 * power as i128)
                .collect();
            result = add_scaled(model, self.prime, &result, &component, coefficient)?;
        }
        Ok(result)
    }

    pub(crate) fn raw_degree(&self, index: usize) -> i64 {
        self.degrees[index]
    }

    pub(crate) fn name(&self, index: usize) -> &str {
        &self.names[index]
    }

    pub(crate) fn index(&self, name: &Label) -> Result<usize> {
        self.indices.get(&name.text).copied().ok_or_else(|| {
            Error::new(
                name.span,
                format!("unknown mixed complex generator '{}'", name.text),
            )
        })
    }

    fn columns(
        &self,
        model: &mut Model,
        rules: Rules,
        op: Differential,
    ) -> Result<Vec<Vec<Coordinate>>> {
        model.spend(self.len())?;
        let mut result = vec![None; self.len()];
        for (source, entries) in rules {
            let index = self.index(&source)?;
            if result[index].is_some() {
                return Err(Error::new(
                    source.span,
                    "duplicate mixed complex differential column",
                ));
            }
            let degree = i128::from(self.degrees[index]) + op.degree();
            result[index] = Some(self.literal(model, degree, entries)?);
        }
        Ok(result.into_iter().map(Option::unwrap_or_default).collect())
    }

    pub(crate) fn literal(
        &self,
        model: &mut Model,
        degree: i128,
        entries: Terms,
    ) -> Result<Vec<Coordinate>> {
        let mut result = BTreeMap::new();
        for (target, coefficient) in entries {
            model.spend(1 + target.text.len())?;
            let i = self.index(&target)?;
            if i128::from(self.degrees[i]) != degree {
                return Err(Error::new(
                    target.span,
                    "mixed complex term has the wrong homological degree",
                ));
            }
            let coefficient =
                i128::from(coefficient).rem_euclid(i128::from(self.prime.get())) as u32;
            let previous = result.entry(i).or_insert(0);
            *previous = self.prime.add(*previous, coefficient);
        }
        Ok(result
            .into_iter()
            .filter(|&(_, coefficient)| coefficient != 0)
            .collect())
    }

    fn validate(&self, model: &mut Model) -> Result<()> {
        if !self.is_strict() {
            // For fixed source and target, homogeneity forces a unique
            // power of u in D², so specializing u=1 tests every coefficient.
            for i in 0..self.len() {
                let image = self.total_column(model, i)?;
                if !self.apply(model, Differential::Total, &image)?.is_empty() {
                    return Err(Error::plain(format!(
                        "higher mixed D squared is nonzero on '{}'",
                        self.names[i]
                    )));
                }
            }
            return Ok(());
        }
        for i in 0..self.len() {
            let down_down = self.apply(model, Differential::Down, &self.down[i])?;
            if !down_down.is_empty() {
                return Err(Error::plain(format!(
                    "b squared is nonzero on '{}'",
                    self.names[i]
                )));
            }
            let up_up = self.apply(model, Differential::Up, &self.up[i])?;
            if !up_up.is_empty() {
                return Err(Error::plain(format!(
                    "B squared is nonzero on '{}'",
                    self.names[i]
                )));
            }
            let down_up = self.apply(model, Differential::Down, &self.up[i])?;
            let up_down = self.apply(model, Differential::Up, &self.down[i])?;
            if !add_scaled(model, self.prime, &down_up, &up_down, 1)?.is_empty() {
                return Err(Error::plain(format!(
                    "bB + Bb is nonzero on '{}'",
                    self.names[i]
                )));
            }
        }
        Ok(())
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        op: Differential,
        input: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        model.spend(1 + input.len())?;
        let mut result = Vec::new();
        for &(index, coefficient) in input {
            if !matches!(op, Differential::Up) {
                result = add_scaled(model, self.prime, &result, &self.down[index], coefficient)?;
            }
            if !matches!(op, Differential::Down) {
                result = add_scaled(model, self.prime, &result, &self.up[index], coefficient)?;
            }
            if let (Differential::Total, Some(higher)) = (op, &self.higher) {
                result = add_scaled(model, self.prime, &result, &higher[index], coefficient)?;
            }
        }
        Ok(result)
    }

    pub(crate) fn filtration(&self, index: usize, b: i128) -> Option<i128> {
        let difference = i128::from(self.degrees[index]) - b;
        (difference >= 0 && difference % 2 == 0).then_some(difference / 2)
    }

    pub(crate) fn basis(&self, model: &mut Model, b: i128, total: bool) -> Result<Vec<usize>> {
        model.spend(self.len())?;
        Ok((0..self.len())
            .filter(|&i| {
                if total {
                    self.filtration(i, b).is_some()
                } else {
                    i128::from(self.degrees[i]) == b
                }
            })
            .collect())
    }

    pub(crate) fn elements(
        &self,
        model: &mut Model,
        grade: Grade,
        total: bool,
        limit: usize,
    ) -> Result<Vec<Value>> {
        let basis = self.basis(model, grade.b.into(), total)?;
        let mut values = odd_ext::elements(model, self.prime, basis.len(), limit)?;
        for value in &mut values {
            let Value::OddClass(support) = value else {
                unreachable!()
            };
            for (index, _) in Rc::make_mut(support) {
                *index = basis[*index];
            }
        }
        Ok(values)
    }
}
