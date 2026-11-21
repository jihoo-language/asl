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
}

impl Complex {
    pub(crate) fn build(
        model: &mut Model,
        field: usize,
        generators: Vec<(Label, i64)>,
        down: Rules,
        up: Rules,
    ) -> Result<Self> {
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
        complex.down = complex.columns(model, down, Differential::Down)?;
        complex.up = complex.columns(model, up, Differential::Up)?;
        complex.validate(model)?;
        model.allocate(
            1 + 6 * complex.len()
                + complex
                    .down
                    .iter()
                    .chain(&complex.up)
                    .map(|column| 2 * column.len())
                    .sum::<usize>()
                + 2 * complex.names.iter().map(String::len).sum::<usize>(),
        )?;
        Ok(complex)
    }

    pub(crate) fn len(&self) -> usize {
        self.names.len()
    }

    pub(crate) fn name(&self, index: usize) -> &str {
        &self.names[index]
    }

    fn index(&self, name: &Label) -> Result<usize> {
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
