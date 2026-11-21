use crate::cohomology;
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::rv_hh::{self, Grade as HochschildGrade};
use crate::rv_hh_homology::{self, Slice};
use crate::rv_negative::{self, Powers, Support, Term};

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) enum Mode {
    Negative,
    Periodic,
    Ordinary,
}

impl Mode {
    pub(crate) fn type_name(self, homology: bool) -> &'static str {
        match (self, homology) {
            (Self::Negative, false) => "RVHNChain",
            (Self::Periodic, false) => "RVHPChain",
            (Self::Ordinary, false) => "RVHCChain",
            (Self::Negative, true) => "RVHN",
            (Self::Periodic, true) => "RVHP",
            (Self::Ordinary, true) => "RVHC",
        }
    }

    pub(crate) fn dimension_name(self) -> &'static str {
        match self {
            Self::Negative => "rvHNdim",
            Self::Periodic => "rvHPdim",
            Self::Ordinary => "rvHCdim",
        }
    }

    pub(crate) fn term_name(self) -> &'static str {
        match self {
            Self::Negative => "rvHNTerm",
            Self::Periodic => "rvHPTerm",
            Self::Ordinary => "rvHCTerm",
        }
    }

    fn powers(self) -> Powers {
        Powers {
            first: if self == Self::Negative { 0 } else { i128::MIN },
            last: if self == Self::Ordinary { 0 } else { i128::MAX },
        }
    }

    fn check_power(self, power: i64) -> Result<()> {
        if !self.powers().includes(power.into()) {
            return Err(Error::plain(match self {
                Self::Negative => "negative cyclic powers must be nonnegative",
                Self::Ordinary => "ordinary cyclic powers must be nonpositive",
                Self::Periodic => unreachable!(),
            }));
        }
        Ok(())
    }

    fn label(self) -> &'static str {
        match self {
            Self::Negative => "negative cyclic homology",
            Self::Periodic => "periodic cyclic homology",
            Self::Ordinary => "cyclic homology",
        }
    }
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) ordinary: HochschildGrade,
    pub(crate) mode: Mode,
}

impl Grade {
    pub(crate) fn term(ordinary: HochschildGrade, mode: Mode, power: i64) -> Result<Self> {
        mode.check_power(power)?;
        Self { ordinary, mode }.shifted(-2 * i128::from(power))
    }

    pub(crate) fn part(&self, power: i64) -> Result<HochschildGrade> {
        self.mode.check_power(power)?;
        self.shifted(2 * i128::from(power))
            .map(|grade| grade.ordinary)
    }

    pub(crate) fn differential(&self) -> Result<Self> {
        self.shifted(-1)
    }

    pub(crate) fn connes_output(&self) -> Result<HochschildGrade> {
        if self.mode != Mode::Ordinary {
            return Err(Error::plain(
                "rvCyclicB requires an ordinary cyclic chain or class",
            ));
        }
        self.shifted(1).map(|grade| grade.ordinary)
    }

    fn shifted(&self, offset: i128) -> Result<Self> {
        let b = i64::try_from(i128::from(self.ordinary.b) + offset)
            .map_err(|_| Error::plain("cyclic degree is out of range"))?;
        Ok(Self {
            ordinary: HochschildGrade {
                b,
                ..self.ordinary.clone()
            },
            mode: self.mode,
        })
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Shift {
    Periodic(i64),
    Suspension,
}

impl Shift {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Periodic(_) => "rvCyclicU",
            Self::Suspension => "rvCyclicS",
        }
    }

    fn power(self) -> i128 {
        match self {
            Self::Periodic(power) => power.into(),
            Self::Suspension => 1,
        }
    }

    pub(crate) fn target(self, grade: &Grade) -> Result<Grade> {
        match (self, grade.mode) {
            (Self::Periodic(_), Mode::Periodic) | (Self::Suspension, Mode::Ordinary) => {
                grade.shifted(-2 * self.power())
            }
            (Self::Periodic(_), _) => Err(Error::plain(
                "rvCyclicU requires a periodic cyclic chain or class",
            )),
            (Self::Suspension, _) => Err(Error::plain(
                "rvCyclicS requires an ordinary cyclic chain or class",
            )),
        }
    }

    pub(crate) fn terms(self, model: &mut Model, mode: Mode, terms: &[Term]) -> Result<Vec<Term>> {
        model.spend(rv_negative::cells(terms))?;
        let powers = mode.powers();
        let power = self.power();
        Ok(terms
            .iter()
            .filter(|term| powers.includes(term.power + power))
            .map(|term| Term {
                power: term.power + power,
                word: term.word.clone(),
            })
            .collect())
    }
}

pub(crate) fn basis(model: &mut Model, grade: &Grade, limit: usize) -> Result<Vec<Term>> {
    rv_negative::basis_in(
        model,
        &grade.ordinary,
        grade.mode.powers(),
        grade.ordinary.b.into(),
        limit,
    )
}

pub(crate) fn build(model: &mut Model, grade: &Grade) -> Result<Slice<Term>> {
    rv_negative::reduce_in(
        model,
        &grade.ordinary,
        grade.mode.powers(),
        grade.mode.label(),
    )
}

pub(crate) fn differential(model: &mut Model, grade: &Grade, support: &[Term]) -> Result<Support> {
    rv_negative::differential_in(model, &grade.ordinary, grade.mode.powers(), support)
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(cohomology::zero());
    }
    let support = model
        .rv_cyclic_homology(ty)?
        .project(model, terms, "rvCyclicClass")?;
    model.allocate(support.len())?;
    Ok(Value::Class(support.into()))
}

pub(crate) fn representative(model: &mut Model, ty: TypeId, class: &[usize]) -> Result<Support> {
    if class.is_empty() {
        return Ok([].into());
    }
    let terms = model.rv_cyclic_homology(ty)?.representative(model, class)?;
    rv_negative::persist(model, terms)
}

pub(crate) fn include_class(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    support: &[usize],
) -> Result<Value> {
    if support.is_empty() {
        return Ok(cohomology::zero());
    }
    let lift = model
        .rv_hh_homology(input)?
        .representative(model, support)?;
    let terms = rv_negative::term(model, 0, &lift)?;
    class(model, output, &terms)
}

fn connes_words(
    model: &mut Model,
    ordinary: &HochschildGrade,
    terms: &[Term],
) -> Result<Vec<rv_hh::Word>> {
    let coefficient = rv_negative::part(model, 0, terms)?;
    rv_hh::connes(model, ordinary, &coefficient)
}

pub(crate) fn connes(model: &mut Model, grade: &Grade, terms: &[Term]) -> Result<rv_hh::Support> {
    let words = connes_words(model, &grade.ordinary, terms)?;
    rv_hh::persist(model, words)
}

pub(crate) fn connes_class(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    support: &[usize],
) -> Result<Value> {
    if support.is_empty() {
        return Ok(cohomology::zero());
    }
    let Kind::RvCyclicHomology(grade) = model.kind(input).clone() else {
        unreachable!()
    };
    let lift = model
        .rv_cyclic_homology(input)?
        .representative(model, support)?;
    let image = connes_words(model, &grade.ordinary, &lift)?;
    rv_hh_homology::class(model, output, &image)
}

pub(crate) fn shift_class(
    model: &mut Model,
    shift: Shift,
    input: TypeId,
    output: TypeId,
    support: &[usize],
) -> Result<Value> {
    if support.is_empty() {
        return Ok(cohomology::zero());
    }
    let Kind::RvCyclicHomology(grade) = model.kind(output).clone() else {
        unreachable!()
    };
    let lift = model
        .rv_cyclic_homology(input)?
        .representative(model, support)?;
    let mapped = shift.terms(model, grade.mode, &lift)?;
    class(model, output, &mapped)
}
