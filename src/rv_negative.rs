use std::cmp::Ordering;
use std::rc::Rc;

use crate::cohomology;
use crate::complex::{xor_sorted, MAX_PIVOT_CELLS};
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::rv_hh::{self, Differential, Grade as HochschildGrade, Word};
use crate::rv_hh_homology::{Slice, MAX_BASIS};
use crate::rv_poly;

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) ordinary: HochschildGrade,
    pub(crate) length: usize,
}

impl Grade {
    pub(crate) fn new(ordinary: HochschildGrade, length: usize) -> Result<Self> {
        if length == 0 {
            return Err(Error::plain(
                "negative cyclic truncation length must be positive",
            ));
        }
        Ok(Self { ordinary, length })
    }

    pub(crate) fn term(ordinary: HochschildGrade, length: usize, power: usize) -> Result<Self> {
        let mut grade = Self::new(ordinary, length)?;
        grade.check_power(power)?;
        grade.ordinary.b = degree(i128::from(grade.ordinary.b) - 2 * power as i128)?;
        Ok(grade)
    }

    pub(crate) fn part(&self, power: usize) -> Result<HochschildGrade> {
        self.check_power(power)?;
        let b = degree(i128::from(self.ordinary.b) + 2 * power as i128)?;
        Ok(HochschildGrade {
            b,
            ..self.ordinary.clone()
        })
    }

    pub(crate) fn differential(&self) -> Result<Self> {
        let b = degree(i128::from(self.ordinary.b) - 1)?;
        Ok(Self {
            ordinary: HochschildGrade {
                b,
                ..self.ordinary.clone()
            },
            length: self.length,
        })
    }

    fn check_power(&self, power: usize) -> Result<()> {
        if power >= self.length {
            return Err(Error::plain(
                "negative cyclic power must be below the truncation length",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Transfer {
    Shift(usize),
    Truncate(usize),
}

impl Transfer {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Shift(_) => "rvNCU",
            Self::Truncate(_) => "rvNCTruncate",
        }
    }

    pub(crate) fn target(self, grade: &Grade) -> Result<Grade> {
        match self {
            Self::Shift(power) => {
                let b = degree(i128::from(grade.ordinary.b) - 2 * power as i128)?;
                Ok(Grade {
                    ordinary: HochschildGrade {
                        b,
                        ..grade.ordinary.clone()
                    },
                    length: grade.length,
                })
            }
            Self::Truncate(length) => {
                if length > grade.length {
                    return Err(Error::plain(
                        "rvNCTruncate cannot increase the truncation length",
                    ));
                }
                Grade::new(grade.ordinary.clone(), length)
            }
        }
    }

    pub(crate) fn terms(
        self,
        model: &mut Model,
        target: &Grade,
        support: &[Term],
    ) -> Result<Vec<Term>> {
        model.spend(cells(support))?;
        let result = match self {
            Self::Shift(power) => support
                .iter()
                .take_while(|term| term.power < target.length.saturating_sub(power))
                .map(|term| Term {
                    power: term.power + power,
                    word: term.word.clone(),
                })
                .collect(),
            Self::Truncate(length) => support
                .iter()
                .take_while(|term| term.power < length)
                .cloned()
                .collect(),
        };
        Ok(result)
    }
}

fn degree(b: i128) -> Result<i64> {
    i64::try_from(b).map_err(|_| Error::plain("negative cyclic degree is out of range"))
}

#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Term {
    pub(crate) power: usize,
    pub(crate) word: Word,
}

pub(crate) type Support = Rc<[Term]>;

pub(crate) fn term_cells(term: &Term) -> usize {
    2 + rv_hh::word_cells(&term.word)
}

fn cells(support: &[Term]) -> usize {
    1 + support.iter().map(term_cells).sum::<usize>()
}

fn check_storage(size: usize) -> Result<()> {
    if size > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "negative cyclic calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

pub(crate) fn persist(model: &mut Model, support: Vec<Term>) -> Result<Support> {
    model.allocate(cells(&support))?;
    Ok(support.into())
}

pub(crate) fn compare(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Ordering> {
    model.spend(cells(left) + cells(right))?;
    Ok(left.iter().rev().cmp(right.iter().rev()))
}

pub(crate) fn add(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Support> {
    model.spend(cells(left) + cells(right))?;
    persist(model, xor_sorted(left, right))
}

pub(crate) fn term(model: &mut Model, power: usize, support: &[Word]) -> Result<Support> {
    model.spend(support.iter().map(|word| rv_hh::word_cells(word)).sum())?;
    let support = support
        .iter()
        .map(|word| Term {
            power,
            word: word.clone(),
        })
        .collect();
    persist(model, support)
}

pub(crate) fn part(model: &mut Model, power: usize, support: &[Term]) -> Result<rv_hh::Support> {
    let start = support.partition_point(|term| term.power < power);
    let end = support.partition_point(|term| term.power <= power);
    model.spend(1 + 2 * support.len().max(1).ilog2() as usize + end - start)?;
    rv_hh::persist(
        model,
        support[start..end]
            .iter()
            .map(|term| term.word.clone())
            .collect(),
    )
}

fn differential_support(
    model: &mut Model,
    grade: &Grade,
    differential: &mut Differential<'_>,
    terms: &[Term],
) -> Result<Vec<Term>> {
    let mut support = Vec::new();
    let mut pending: Option<(usize, Vec<Word>)> = None;
    let mut retained = 0;
    for chunk in terms.chunk_by(|a, b| a.power == b.power) {
        let power = chunk[0].power;
        model.spend(1 + chunk.len())?;
        let words: Vec<_> = chunk.iter().map(|term| term.word.clone()).collect();
        let mut image = differential.apply(model, &words)?;
        if let Some((p, part)) = pending.take() {
            if p == power {
                model.spend(
                    image
                        .iter()
                        .chain(&part)
                        .map(|word| rv_hh::word_cells(word))
                        .sum(),
                )?;
                image = xor_sorted(&image, &part);
            } else {
                append(&mut support, &mut retained, p, part)?;
            }
        }
        append(&mut support, &mut retained, power, image)?;
        if power + 1 < grade.length {
            pending = Some((power + 1, rv_hh::connes(model, &grade.ordinary, &words)?));
        }
    }
    if let Some((p, part)) = pending.take() {
        append(&mut support, &mut retained, p, part)?;
    }
    Ok(support)
}

fn append(
    output: &mut Vec<Term>,
    retained: &mut usize,
    power: usize,
    words: Vec<Word>,
) -> Result<()> {
    for word in words {
        let term = Term { power, word };
        *retained += term_cells(&term);
        check_storage(*retained)?;
        output.push(term);
    }
    Ok(())
}

pub(crate) fn differential(model: &mut Model, grade: &Grade, support: &[Term]) -> Result<Support> {
    if support.is_empty() {
        return Ok([].into());
    }
    let ideal = model.rv_ideal_named(&grade.ordinary.ideal).unwrap();
    let mut differential = Differential::new(model, &grade.ordinary, &ideal)?;
    let result = differential_support(model, grade, &mut differential, support)?;
    persist(model, result)
}

pub(crate) fn basis(
    model: &mut Model,
    grade: &Grade,
    total: i128,
    limit: usize,
) -> Result<Vec<Term>> {
    let lo = -total.div_euclid(2);
    let hi = (2 * rv_poly::degree(&grade.ordinary.weight)? as i128 - total).div_euclid(2);
    let first = lo.max(0);
    let last = hi.min(grade.length as i128 - 1);
    if first > last {
        return Ok(Vec::new());
    }
    let mut terms = Vec::new();
    let mut retained = 0;
    for power in first as usize..=last as usize {
        model.spend(1)?;
        let original = total + 2 * power as i128;
        let words = rv_hh::basis_at(model, &grade.ordinary, original, limit - terms.len())?;
        append(&mut terms, &mut retained, power, words)?;
    }
    Ok(terms)
}

pub(crate) fn build(model: &mut Model, grade: &Grade) -> Result<Slice<Term>> {
    let total = i128::from(grade.ordinary.b);
    let terms = basis(model, grade, total, MAX_BASIS)?;
    let ideal = model.rv_ideal_named(&grade.ordinary.ideal).unwrap();
    let mut differential = Differential::new(model, &grade.ordinary, &ideal)?;
    Slice::build(
        model,
        terms,
        |model| basis(model, grade, total + 1, MAX_BASIS),
        term_cells,
        "finite negative cyclic homology",
        |model, term| {
            differential_support(model, grade, &mut differential, std::slice::from_ref(term))
        },
    )
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(cohomology::zero());
    }
    let support = model
        .rv_negative_homology(ty)?
        .project(model, terms, "rvNCClass")?;
    model.allocate(support.len())?;
    Ok(Value::Class(support.into()))
}

pub(crate) fn representative(model: &mut Model, ty: TypeId, class: &[usize]) -> Result<Support> {
    if class.is_empty() {
        return Ok([].into());
    }
    let terms = model
        .rv_negative_homology(ty)?
        .representative(model, class)?;
    persist(model, terms)
}

pub(crate) fn transfer_class(
    model: &mut Model,
    operation: Transfer,
    input: TypeId,
    output: TypeId,
    support: &[usize],
) -> Result<Value> {
    if support.is_empty() {
        return Ok(cohomology::zero());
    }
    let Kind::RvNegativeCyclic(target) = model.kind(output).clone() else {
        unreachable!()
    };
    let representative = model
        .rv_negative_homology(input)?
        .representative(model, support)?;
    let mapped = operation.terms(model, &target, &representative)?;
    class(model, output, &mapped)
}
