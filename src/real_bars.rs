use std::cmp::Ordering;
use std::collections::HashMap;

use crate::cobar;
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::linear::Basis;
use crate::model::{Model, TypeId, Value};
use crate::motivic::Grade;
use crate::real_ext;
use crate::real_motivic::{self, Cochain, Term};
use crate::rho_localized;

const MAX_BASIS: usize = 32_768;

#[derive(Default)]
pub(crate) struct Group {
    bars: Vec<Bar>,
}

pub(crate) struct Bar {
    birth: usize,
    representative: Cochain,
    finite: Option<Finite>,
}

struct Finite {
    length: usize,
    witness: Cochain,
}

impl Bar {
    pub(crate) fn exponent_at(&self, t: i128) -> Result<usize> {
        if !self.contains_shift(t, t) {
            return Err(Error::plain(
                "rho bar does not live in the requested degree",
            ));
        }
        let exponent = (self.birth as i128)
            .checked_sub(t)
            .and_then(|a| usize::try_from(a).ok());
        exponent.ok_or_else(|| Error::plain("rho bar exponent is outside the supported range"))
    }

    pub(crate) fn lift_at(&self, model: &mut Model, t: i128) -> Result<Vec<Term>> {
        real_motivic::shift_support(model, &self.representative, self.exponent_at(t)?)
    }

    pub(crate) fn grade(&self, s: usize, k: i128) -> Result<Grade> {
        rho_localized::section_grade(s, self.birth as i128, classical_degree(k)?)
    }

    pub(crate) fn witness_grade(&self, s: usize, k: i128) -> Result<Grade> {
        let length = self.finite.as_ref().ok_or_else(free_bar)?.length;
        let death = self.birth - length;
        rho_localized::section_grade(s - 1, death as i128, classical_degree(k)?)
    }

    pub(crate) fn representative(&self) -> Value {
        Value::RMCobar(self.representative.clone())
    }

    pub(crate) fn generator(&self, model: &mut Model, ty: TypeId) -> Result<Value> {
        let value = real_ext::class(model, ty, &self.representative)?;
        let Value::Class(coords) = &value else {
            unreachable!();
        };
        if coords.is_empty() {
            return Err(Error::plain(
                "rho bar representative is a boundary at its birth",
            ));
        }
        Ok(value)
    }

    pub(crate) fn witness(&self, model: &mut Model) -> Result<Value> {
        let finite = self.finite.as_ref().ok_or_else(free_bar)?;
        let image = real_motivic::differential_support(model, &finite.witness)?;
        let target = real_motivic::shift_support(model, &self.representative, finite.length)?;
        if real_motivic::compare(model, &image, &target)? != Ordering::Equal {
            return Err(Error::plain(
                "rho bar witness does not have the required differential",
            ));
        }
        Ok(Value::RMCobar(finite.witness.clone()))
    }

    fn cells(&self) -> usize {
        1 + self.representative.iter().map(Term::cells).sum::<usize>()
            + self.finite.as_ref().map_or(0, |finite| {
                1 + finite.witness.iter().map(Term::cells).sum::<usize>()
            })
    }

    fn contains_shift(&self, t: i128, target: i128) -> bool {
        t <= self.birth as i128
            && self
                .finite
                .as_ref()
                .is_none_or(|finite| target > (self.birth - finite.length) as i128)
    }
}

fn classical_degree(k: i128) -> Result<usize> {
    usize::try_from(k)
        .map_err(|_| Error::plain("rho bar degree t-w is outside the supported range"))
}

fn free_bar() -> Error {
    Error::plain("rmbarNull requires a finite rho bar")
}

impl Group {
    pub(crate) fn build(model: &mut Model, s: usize, k: i128) -> Result<Self> {
        if k < 0 {
            return Ok(Self::default());
        }
        let grade = rho_localized::section_grade(s, 0, classical_degree(k)?)?;
        let middle = stable_basis(model, grade)?;
        if middle.is_empty() {
            return Ok(Self::default());
        }
        let mut group = Self::default();
        let mut cells = 0;
        let incoming = if let Some(previous) = grade.previous() {
            let source = stable_basis(model, previous)?;
            let target = middle
                .iter()
                .enumerate()
                .map(|(index, term)| {
                    model.spend(term.cells())?;
                    Ok((term, index))
                })
                .collect::<Result<HashMap<_, _>>>()?;
            reduce_columns(
                model,
                &source,
                None,
                |model, image| coordinates(model, &target, image),
                |model, index, image, witness| {
                    let Some(&pivot) = image.last() else {
                        return Ok(());
                    };
                    let birth = middle[pivot].rho;
                    let death = source[index].rho;
                    let length = birth.checked_sub(death).ok_or_else(|| {
                        Error::plain("rho bar differential violates its filtration")
                    })?;
                    if length != 0 {
                        let representative = lift(model, &middle, image, birth)?;
                        let witness = lift(model, &source, witness, death)?;
                        let finite = Finite { length, witness };
                        group.record(
                            &mut cells,
                            Bar {
                                birth,
                                representative,
                                finite: Some(finite),
                            },
                        )?;
                    }
                    Ok(())
                },
            )?
        } else {
            Basis::default()
        };

        // Each incoming pivot already has a cycle with that leading
        // coordinate. Its outgoing image is in the earlier span.
        // Only dependence among the ordered source columns is needed here,
        // so the outgoing rows can be numbered in encounter order.
        let mut rows = real_motivic::Rows::default();
        reduce_columns(
            model,
            &middle,
            Some(&incoming),
            |model, image| rows.column(model, image),
            |model, index, image, witness| {
                if image.is_empty() {
                    let birth = middle[index].rho;
                    let representative = lift(model, &middle, witness, birth)?;
                    group.record(
                        &mut cells,
                        Bar {
                            birth,
                            representative,
                            finite: None,
                        },
                    )?;
                }
                Ok(())
            },
        )?;

        let count = group.bars.len();
        model.spend(count * (1 + count.max(1).ilog2() as usize))?;
        // Preserve the reduction order among bars with equal endpoints.
        group.bars.sort_by_key(|bar| {
            (
                bar.birth,
                bar.finite.as_ref().map_or(usize::MAX, |f| f.length),
            )
        });
        model.allocate(cells)?;
        Ok(group)
    }

    fn record(&mut self, cells: &mut usize, bar: Bar) -> Result<()> {
        *cells += bar.cells();
        if *cells > MAX_PIVOT_CELLS {
            return Err(Error::plain("rho bar witnesses exceed the storage limit"));
        }
        self.bars.push(bar);
        Ok(())
    }

    pub(crate) fn get(&self, index: usize) -> Result<&Bar> {
        self.bars.get(index).ok_or_else(|| {
            Error::plain(format!(
                "rho bar index {index} is out of range for {} bars",
                self.bars.len(),
            ))
        })
    }

    pub(crate) fn living(&self, model: &mut Model, t: i128) -> Result<Vec<usize>> {
        model.spend(self.bars.len())?;
        Ok(self
            .bars
            .iter()
            .enumerate()
            .filter_map(|(index, bar)| bar.contains_shift(t, t).then_some(index))
            .collect())
    }

    pub(crate) fn describe(&self, model: &mut Model) -> Result<String> {
        model.spend(self.bars.len())?;
        let bars = self
            .bars
            .iter()
            .map(|bar| {
                let birth = bar.birth;
                match &bar.finite {
                    Some(finite) => format!("({birth}, {})", finite.length),
                    None => format!("({birth}, infinite)"),
                }
            })
            .collect::<Vec<_>>()
            .join(", ");
        Ok(format!("[{bars}]"))
    }

    pub(crate) fn shift_rank(&self, model: &mut Model, t: i128, power: usize) -> Result<usize> {
        let target = t
            .checked_sub(power as i128)
            .ok_or_else(|| Error::plain("rmbarRank target degree is out of range"))?;
        model.spend(self.bars.len())?;
        Ok(self
            .bars
            .iter()
            .filter(|bar| bar.contains_shift(t, target))
            .count())
    }
}

fn stable_basis(model: &mut Model, grade: Grade) -> Result<Vec<Term>> {
    let mut terms = real_motivic::basis(model, grade, MAX_BASIS)?;
    model.spend(terms.len())?;
    terms.reverse();
    Ok(terms)
}

fn reduce_columns(
    model: &mut Model,
    source: &[Term],
    paired: Option<&Basis>,
    mut encode: impl FnMut(&mut Model, &[Term]) -> Result<Vec<usize>>,
    mut use_column: impl FnMut(&mut Model, usize, &[usize], &[usize]) -> Result<()>,
) -> Result<Basis> {
    let mut pivots = Basis::default();
    for (index, term) in source.iter().enumerate() {
        model.spend(1)?;
        if paired.is_some_and(|incoming| incoming.contains_pivot(index)) {
            continue;
        }
        let image = real_motivic::differential_support(model, std::slice::from_ref(term))?;
        let column = encode(model, &image)?;
        pivots.insert_traced(model, column, vec![index], |model, image, witness| {
            use_column(model, index, image, witness)
        })?;
    }
    Ok(pivots)
}

fn coordinates(
    model: &mut Model,
    basis: &HashMap<&Term, usize>,
    terms: &[Term],
) -> Result<Vec<usize>> {
    let mut result = Vec::with_capacity(terms.len());
    for term in terms {
        model.spend(2 * term.cells())?;
        let &coordinate = basis.get(term).ok_or_else(|| {
            Error::plain("rho differential term is outside the stable cobar basis")
        })?;
        result.push(coordinate);
    }
    // The cochain terms are ascending, while the persistence basis is reversed.
    result.reverse();
    Ok(result)
}

fn lift(model: &mut Model, basis: &[Term], support: &[usize], degree: usize) -> Result<Cochain> {
    let mut terms = Vec::with_capacity(support.len());
    let mut cells = 0;
    for &coordinate in support.iter().rev() {
        let term = &basis[coordinate];
        let rho = term
            .rho
            .checked_sub(degree)
            .ok_or_else(|| Error::plain("rho bar lift has a negative coefficient exponent"))?;
        cells += term.cells();
        cobar::check_storage(cells)?;
        model.spend(term.cells())?;
        terms.push(Term {
            rho,
            ..term.clone()
        });
    }
    Ok(terms.into())
}
