use std::rc::Rc;

use crate::bp::{self, Shape};
use crate::bp_classes::{self, Grade, Support};
use crate::bp_ext;
use crate::bp_polynomial::Polynomial;
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result, Span};
use crate::model::{Model, TypeId, Value};

type Cochain = Rc<Polynomial>;
pub(crate) type Witness = (Span, Cochain);

pub(crate) fn interval_count(arity: usize) -> Result<usize> {
    if arity < 3 {
        return Err(Error::plain(
            "higher BP Massey operations require at least three inputs",
        ));
    }
    arity
        .checked_add(1)
        .and_then(|n| arity.checked_mul(n))
        .map(|twice| twice / 2)
        .filter(|&size| size <= MAX_PIVOT_CELLS / 8)
        .ok_or_else(|| Error::plain("BP defining system exceeds the interval storage limit"))
}

pub(crate) struct Layout {
    arity: usize,
    // Intervals are stored in increasing length, then increasing start.
    offsets: Vec<usize>,
    excess: Vec<usize>,
    grades: Vec<Grade>,
    pub(crate) target: Grade,
}

impl Layout {
    pub(crate) fn new(model: &mut Model, inputs: Vec<Grade>) -> Result<Self> {
        let arity = inputs.len();
        let count = interval_count(arity)?;
        model.spend(1 + count * 3)?;
        let base = Shape {
            arrows: 0,
            ..inputs[0].shape
        };
        let mut excess = vec![0usize];
        let mut weights = vec![0usize];
        for input in inputs {
            if (Shape {
                arrows: 0,
                ..input.shape
            }) != base
            {
                return Err(Error::plain(
                    "higher BP inputs must have the same prime and height",
                ));
            }
            let degree =
                input.shape.arrows.checked_sub(1).ok_or_else(|| {
                    Error::plain("higher BP inputs must have positive cobar degrees")
                })?;
            let new_excess =
                excess.last().unwrap().checked_add(degree).ok_or_else(|| {
                    Error::plain("BP defining-system cobar degree is out of range")
                })?;
            let new_weight = weights
                .last()
                .unwrap()
                .checked_add(input.t)
                .ok_or_else(|| {
                    Error::plain("BP defining-system internal degree is out of range")
                })?;
            excess.push(new_excess);
            weights.push(new_weight);
        }
        let mut offsets = vec![0];
        let mut grades = Vec::with_capacity(count);
        for length in 1..=arity {
            offsets.push(grades.len());
            for i in 0..=arity - length {
                let j = i + length;
                let arrows = 1 + (excess[j] - excess[i]);
                let shape = base.with_arrows(arrows)?;
                grades.push(Grade {
                    shape,
                    t: weights[j] - weights[i],
                });
            }
        }
        let interval = *grades.last().unwrap();
        let target = Grade {
            shape: interval.shape.next()?,
            ..interval
        };
        Ok(Self {
            arity,
            offsets,
            excess,
            grades,
            target,
        })
    }

    pub(crate) fn witness_grades(&self) -> &[Grade] {
        &self.grades[self.arity..self.grades.len() - 1]
    }

    fn index(&self, i: usize, j: usize) -> usize {
        self.offsets[j - i] + i
    }

    fn obstruction(
        &self,
        model: &mut Model,
        entries: &[Cochain],
        i: usize,
        j: usize,
    ) -> Result<Polynomial> {
        model.spend(j - i - 1)?;
        let mut output = Polynomial::default();
        for k in (i + 1)..j {
            let left = self.index(i, k);
            let right = self.index(k, j);
            if entries[left].is_zero() || entries[right].is_zero() {
                continue;
            }
            let product = bp::cup(
                model,
                self.grades[left].shape,
                self.grades[right].shape,
                &entries[left],
                &entries[right],
            )?;
            let negative = (self.excess[k - 1] - self.excess[i]) % 2 == 1;
            output.add_assign(model, &product, negative)?;
        }
        Ok(output)
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        leaves: Vec<(TypeId, Support)>,
        witnesses: Option<Vec<Witness>>,
        target: TypeId,
    ) -> Result<Value> {
        let mut entries = Vec::with_capacity(self.grades.len() - 1);
        let mut cells = self.grades.len();
        for (ty, support) in leaves {
            let representative = if support.is_empty() {
                Polynomial::default()
            } else {
                model.bp_ext(ty)?.representative(model, &support)?
            };
            keep(&mut entries, &mut cells, Rc::new(representative))?;
        }
        let mut witnesses = witnesses.map(Vec::into_iter);
        for length in 2..self.arity {
            for i in 0..=self.arity - length {
                let j = i + length;
                let theta = self.obstruction(model, &entries, i, j)?;
                let grade = self.grades[self.index(i, j)];
                let entry = if let Some(supplied) = &mut witnesses {
                    let (span, value) = supplied.next().unwrap();
                    let image = if value.is_zero() {
                        Polynomial::default()
                    } else {
                        bp_ext::validate_homogeneous(model, grade.shape, grade.t, &value)
                            .map_err(|error| error.at(span))?;
                        bp::differential(model, grade.shape, &value)
                            .map_err(|error| error.at(span))?
                    };
                    if !theta.compare(model, &image)?.is_eq() {
                        return Err(Error::new(
                            span,
                            format!(
                                "bpHigherSystem has the wrong differential on interval [{i},{j})"
                            ),
                        ));
                    }
                    value
                } else {
                    Rc::new(self.choose(model, grade, &theta, i, j)?)
                };
                keep(&mut entries, &mut cells, entry)?;
            }
        }
        let terminal = self.obstruction(model, &entries, 0, self.arity)?;
        bp_classes::class(model, target, &terminal)
    }

    fn choose(
        &self,
        model: &mut Model,
        grade: Grade,
        theta: &Polynomial,
        i: usize,
        j: usize,
    ) -> Result<Polynomial> {
        if theta.is_zero() {
            return Ok(Polynomial::default());
        }
        let shape = grade.shape.next()?;
        let ty = model.bp_ext_type(Grade { shape, ..grade })?;
        let group = model.bp_ext(ty)?;
        if !group.project(model, theta)?.is_empty() {
            return Err(Error::plain(format!(
                "bpChosenHigher has an obstruction in its chosen system on interval [{i},{j})"
            )));
        }
        bp_ext::null(model, shape, grade.t, theta)
    }
}

fn keep(entries: &mut Vec<Cochain>, cells: &mut usize, cochain: Cochain) -> Result<()> {
    *cells += cochain.cells();
    if *cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "BP defining system exceeds the cochain storage limit",
        ));
    }
    entries.push(cochain);
    Ok(())
}
