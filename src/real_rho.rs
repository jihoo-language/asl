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
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Reduce(_) => "rhoReduce",
            Self::Lift => "rhoTruncLift",
            Self::Connect => "rhoConnect",
            Self::ConnectRep => "rhoConnectRep",
            Self::Beta => "rhoBeta",
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
    let mut terms = Vec::new();
    if !coords.is_empty() {
        let lift = model.real_ext(input)?.representative(model, coords)?;
        terms = real_motivic::differential_support(model, &lift)?;
        for term in &mut terms {
            term.rho = term.rho.checked_sub(grade.exponent).ok_or_else(|| {
                Error::plain(format!(
                    "connecting lift is not closed modulo rho^{}",
                    grade.exponent
                ))
            })?;
        }
    }
    match op {
        Operation::ConnectRep => real_motivic::persist(model, terms),
        Operation::Connect => real_ext::class(model, output, &terms),
        Operation::Beta => real_ext::class(model, output, cut(&terms, grade.exponent)),
        _ => unreachable!(),
    }
}
