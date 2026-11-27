use crate::cobar::{self, Word};
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::motivic::{self, Grade as ChainGrade};
use crate::real_ext;
use crate::real_motivic::{self, Term};

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) exponent: usize,
    pub(crate) chain: ChainGrade,
}

impl Grade {
    pub(crate) fn new(exponent: usize, chain: ChainGrade) -> Result<Self> {
        if exponent == 0 {
            return Err(Error::plain("rho truncation exponent must be positive"));
        }
        Ok(Self { exponent, chain })
    }

    pub(crate) fn differential(self) -> Result<Self> {
        Ok(Self {
            chain: self.chain.differential()?,
            ..self
        })
    }

    pub(crate) fn connecting(self) -> Result<ChainGrade> {
        let mut chain = self.chain.differential()?;
        chain.t = chain
            .t
            .checked_add(self.exponent as i128)
            .ok_or_else(grade_limit)?;
        chain.w = i64::try_from(i128::from(chain.w) + self.exponent as i128)
            .map_err(|_| grade_limit())?;
        Ok(chain)
    }

    pub(crate) fn product(self, right: Self) -> Result<Self> {
        if self.exponent != right.exponent {
            return Err(Error::plain(
                "rho-truncated products require the same truncation exponent",
            ));
        }
        Ok(Self {
            chain: self.chain.product(right.chain)?,
            ..self
        })
    }
}

fn grade_limit() -> Error {
    Error::plain("rho connecting tridegree exceeds the supported range")
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Reduce(usize),
    Lift,
    Connect,
    ConnectRep,
    Beta,
    BetaRep,
    BetaNull,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Reduce(_) => "rhoReduce",
            Self::Lift => "rhoTruncLift",
            Self::Connect => "rhoConnect",
            Self::ConnectRep => "rhoConnectRep",
            Self::Beta => "rhoBeta",
            Self::BetaRep => "rhoBetaRep",
            Self::BetaNull => "rhoBetaNull",
        }
    }
}

pub(crate) fn basis(model: &mut Model, grade: Grade, limit: usize) -> Result<Vec<Term>> {
    real_motivic::basis_window(model, grade.chain, 0..=grade.exponent - 1, limit)
}

fn cut(terms: &[Term], exponent: usize) -> &[Term] {
    &terms[..terms.partition_point(|term| term.rho < exponent)]
}

pub(crate) fn reduce(model: &mut Model, terms: &[Term], exponent: usize) -> Result<Value> {
    let terms = cut(terms, exponent);
    model.spend(terms.iter().map(Term::cells).sum())?;
    real_motivic::persist(model, terms.to_vec())
}

pub(crate) fn reduce_class(
    model: &mut Model,
    input: TypeId,
    coords: &[usize],
    output: TypeId,
    exponent: usize,
) -> Result<Value> {
    if coords.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let terms = model.real_ext(input)?.representative(model, coords)?;
    real_ext::class(model, output, cut(&terms, exponent))
}

pub(crate) fn from_complex(
    model: &mut Model,
    grade: ChainGrade,
    words: &[Word],
) -> Result<Vec<Term>> {
    if words.is_empty() {
        return Ok(Vec::new());
    }
    let t = grade.classical_degree()?;
    let mut terms = Vec::with_capacity(words.len());
    let mut cells = 0;
    for word in words {
        let weight = motivic::intrinsic_weight(word, t) as i128;
        let tau = usize::try_from(weight - i128::from(grade.w)).map_err(|_| grade_limit())?;
        let term = Term {
            rho: 0,
            tau,
            word: word.clone(),
        };
        cells += term.cells();
        cobar::check_storage(cells)?;
        terms.push(term);
    }
    model.spend(cells * (1 + terms.len().max(1).ilog2() as usize))?;
    terms.sort_unstable();
    Ok(terms)
}

pub(crate) fn shift(
    model: &mut Model,
    terms: &[Term],
    exponent: usize,
    power: usize,
) -> Result<Vec<Term>> {
    if power >= exponent {
        return Ok(Vec::new());
    }
    real_motivic::shift_support(model, cut(terms, exponent - power), power)
}

pub(crate) fn shift_class(
    model: &mut Model,
    input: TypeId,
    coords: &[usize],
    output: TypeId,
    exponent: usize,
    power: usize,
) -> Result<Value> {
    if coords.is_empty() || power >= exponent {
        return Ok(Value::Class([].into()));
    }
    let terms = model.real_ext(input)?.representative(model, coords)?;
    let terms = shift(model, &terms, exponent, power)?;
    real_ext::class(model, output, &terms)
}

pub(crate) fn differential(model: &mut Model, grade: Grade, terms: &[Term]) -> Result<Value> {
    let terms = real_motivic::differential_below(model, terms, grade.exponent as u128)?;
    real_motivic::persist(model, terms)
}

pub(crate) fn concat(model: &mut Model, grade: Grade, a: &[Term], b: &[Term]) -> Result<Value> {
    let terms = real_motivic::concat_below(model, a, b, grade.exponent as u128)?;
    real_motivic::persist(model, terms)
}

pub(crate) fn yoneda(
    model: &mut Model,
    a: (TypeId, &[usize]),
    b: (TypeId, &[usize]),
    target: TypeId,
    grade: Grade,
) -> Result<Value> {
    if a.1.is_empty() || b.1.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let left = model.real_ext(a.0)?.representative(model, a.1)?;
    let right = model.real_ext(b.0)?.representative(model, b.1)?;
    let terms = real_motivic::concat_below(model, &left, &right, grade.exponent as u128)?;
    real_ext::class(model, target, &terms)
}

fn boundary_witness(
    model: &mut Model,
    grade: Grade,
    terms: &[Term],
    alternate: impl FnOnce(&mut Model) -> Result<Option<Vec<Term>>>,
) -> Result<Option<Vec<Term>>> {
    let Some(witness) =
        real_ext::truncated_preimage(model, grade.chain, grade.exponent, terms, alternate)?
    else {
        return Ok(None);
    };
    let image = real_motivic::differential_below(model, &witness, grade.exponent as u128)?;
    if !real_motivic::compare(model, &image, terms)?.is_eq() {
        return Err(Error::plain(
            "incorrect rho-truncated cobar boundary witness",
        ));
    }
    Ok(Some(witness))
}

pub(crate) fn primitive(
    model: &mut Model,
    grade: Grade,
    terms: &[Term],
    operation: &str,
) -> Result<Value> {
    let witness = boundary_witness(model, grade, terms, |_| Ok(None))?.ok_or_else(|| {
        Error::plain(format!(
            "{operation} requires a rho-truncated cobar boundary"
        ))
    })?;
    real_motivic::persist(model, witness)
}

pub(crate) fn connecting(
    model: &mut Model,
    op: Operation,
    input: TypeId,
    coords: &[usize],
    output: TypeId,
) -> Result<Value> {
    let Kind::RMRhoExt(grade) = *model.kind(input) else {
        unreachable!()
    };
    let lift = if coords.is_empty() {
        Vec::new()
    } else {
        model.real_ext(input)?.representative(model, coords)?
    };
    let mut terms = real_motivic::differential_support(model, &lift)?;
    for term in &mut terms {
        term.rho = term.rho.checked_sub(grade.exponent).ok_or_else(|| {
            Error::plain(format!(
                "connecting lift is not closed modulo rho^{}",
                grade.exponent
            ))
        })?;
    }
    match op {
        Operation::ConnectRep => real_motivic::persist(model, terms),
        Operation::Connect => real_ext::class(model, output, &terms),
        Operation::BetaRep => reduce(model, &terms, grade.exponent),
        Operation::Beta | Operation::BetaNull => {
            let terms = cut(&terms, grade.exponent);
            if matches!(op, Operation::Beta) && model.has_real_ext(output) {
                return real_ext::class(model, output, terms);
            }
            let target = Grade {
                chain: grade.connecting()?,
                ..grade
            };
            let witness = boundary_witness(model, target, terms, |model| {
                real_ext::degree_one_beta_witness(model, grade.chain, grade.exponent, &lift)
            })?;
            match (op, witness) {
                (Operation::BetaNull, Some(witness)) => real_motivic::persist(model, witness),
                (Operation::BetaNull, None) => Err(Error::plain(
                    "rhoBetaNull requires a rho-truncated cobar boundary",
                )),
                (_, Some(_)) => Ok(Value::Class([].into())),
                (_, None) => real_ext::class(model, output, terms),
            }
        }
        _ => unreachable!(),
    }
}
