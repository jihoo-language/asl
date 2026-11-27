use std::cmp::Ordering;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::local::LocalAlgebra;
use super::{Binary, Unary};

pub(crate) enum Statement {
    Declare(Label, Box<ParsedExpr>),
    Trace(Box<ParsedExpr>),
    Dimension(Label),
}

pub(crate) enum Expression {
    Map([Label; 2], Box<ParsedExpr>),
    Identity(Label),
    Compose(Box<[ParsedExpr; 2]>),
    Underlying(Box<ParsedExpr>),
}

pub(crate) struct Map {
    source: Rc<str>,
    target: Rc<str>,
    underlying: Rc<LocalAlgebra>,
}

pub(crate) fn certify_object(model: &mut Model, e: &LocalAlgebra) -> Result<()> {
    if e.binary(model, Binary::Compose, e)? != *e {
        return Err(Error::plain(
            "Karoubi object requires an idempotent diagram",
        ));
    }
    Ok(())
}

impl Map {
    pub(crate) fn corner(
        model: &mut Model,
        names: &[Label; 2],
        underlying: Rc<LocalAlgebra>,
    ) -> Result<Self> {
        let [source, target] = names;
        let a = model.karoubi_object(source)?;
        let b = model.karoubi_object(target)?;
        if underlying.rank() != a.rank() || a.rank() != b.rank() {
            return Err(Error::plain("Karoubi map rank must match both retracts"));
        }
        let right = underlying.binary(model, Binary::Compose, &a)?;
        if b.binary(model, Binary::Compose, &right)? != *underlying {
            return Err(Error::plain(
                "Karoubi map must satisfy f = e_target f e_source",
            ));
        }
        Ok(Self {
            source: source.text.as_str().into(),
            target: target.text.as_str().into(),
            underlying,
        })
    }

    pub(crate) fn identity(name: &Label, projector: Rc<LocalAlgebra>) -> Self {
        let source: Rc<str> = name.text.as_str().into();
        Self {
            source: source.clone(),
            target: source,
            underlying: projector,
        }
    }

    pub(crate) fn zero(source: Rc<str>, target: Rc<str>, rank: usize) -> Self {
        Self {
            source,
            target,
            underlying: Rc::new(LocalAlgebra::zero(rank)),
        }
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        self.underlying.charge(model)?;
        model.allocate(4 + self.source.len() + self.target.len())?;
        Ok(Value::RVKaroubiMap(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        let names = (&self.source, &self.target).cmp(&(&rhs.source, &rhs.target));
        if !names.is_eq() {
            return Ok(names);
        }
        self.underlying.compare(model, &rhs.underlying)
    }

    pub(crate) fn render(&self) -> String {
        format!(
            "RVKaroubiMap({}, {}){{{}}}",
            self.source,
            self.target,
            self.underlying.render()
        )
    }

    pub(crate) fn underlying(&self) -> &LocalAlgebra {
        &self.underlying
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        debug_assert!(self.source == rhs.source && self.target == rhs.target);
        Ok(Self {
            source: self.source.clone(),
            target: self.target.clone(),
            underlying: Rc::new(self.underlying.add(model, &rhs.underlying)?),
        })
    }

    pub(crate) fn negate(&self, model: &mut Model) -> Result<Self> {
        Ok(Self {
            source: self.source.clone(),
            target: self.target.clone(),
            underlying: Rc::new(self.underlying.unary(model, Unary::Scale(-1, 0))?),
        })
    }

    pub(crate) fn compose(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        if self.source != rhs.target {
            return Err(Error::plain(
                "Karoubi composition requires the same intermediate object",
            ));
        }
        Ok(Self {
            source: rhs.source.clone(),
            target: self.target.clone(),
            underlying: Rc::new(
                self.underlying
                    .binary(model, Binary::Compose, &rhs.underlying)?,
            ),
        })
    }

    pub(crate) fn trace(&self, model: &mut Model) -> Result<String> {
        if self.source != self.target {
            return Err(Error::plain("Karoubi trace requires an endomorphism"));
        }
        self.underlying.trace_report(model, "rvKaroubiTrace")
    }
}
