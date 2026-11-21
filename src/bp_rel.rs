use std::rc::Rc;

use crate::bp;
use crate::bp_ext::{self, Domain};
use crate::bp_polynomial::{Monomial, Polynomial};
use crate::error::{Error, Result};
use crate::model::{Model, Value};

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct Shape {
    pub(crate) upper: bp::Shape,
    pub(crate) lower: usize,
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct Grade {
    pub(crate) shape: Shape,
    pub(crate) t: usize,
}

impl Shape {
    pub(crate) fn new(
        model: &mut Model,
        prime: usize,
        lower: usize,
        upper: usize,
        arrows: usize,
    ) -> Result<Self> {
        Self::with_lower(bp::Shape::new(model, prime, upper, arrows)?, lower)
    }

    pub(crate) fn with_lower(upper: bp::Shape, lower: usize) -> Result<Self> {
        upper.require_integral("bpRel")?;
        if lower >= upper.height {
            return Err(Error::plain("BP relative heights must satisfy n < m"));
        }
        Ok(Self { upper, lower })
    }

    pub(crate) fn next(self) -> Result<Self> {
        Ok(Self {
            upper: self.upper.next()?,
            ..self
        })
    }

    pub(crate) fn previous(self) -> Result<Self> {
        Ok(Self {
            upper: bp_ext::previous_shape(self.upper)?,
            ..self
        })
    }

    pub(crate) fn lower_shape(self) -> bp::Shape {
        bp::Shape {
            height: self.lower,
            ..self.upper
        }
    }

    pub(crate) fn domain(self) -> Domain {
        Domain::relative(self.upper, self.lower)
    }
}

pub(crate) fn retains(monomial: &Monomial, lower: usize, upper: usize) -> bool {
    monomial
        .iter()
        .any(|&(variable, _)| variable % upper >= lower)
}

pub(crate) fn quotient(model: &mut Model, shape: Shape, input: &Polynomial) -> Result<Polynomial> {
    bp_ext::check_normalized(model, shape.upper, input)?;
    shape.domain().project(model, input)
}

pub(crate) fn connecting(
    model: &mut Model,
    shape: Shape,
    cycle: &Polynomial,
) -> Result<Polynomial> {
    let image = bp::differential(model, shape.upper, cycle)?;
    model.spend(image.cells())?;
    if image
        .terms()
        .any(|(monomial, _)| retains(monomial, shape.lower, shape.upper.height))
    {
        return Err(Error::plain(
            "bpHeightDelta requires a closed relative cochain",
        ));
    }
    bp::reindex_height(model, shape.upper.height, shape.lower, &image)
}

pub(crate) fn persist(model: &mut Model, cochain: Polynomial) -> Result<Value> {
    model.allocate(cochain.cells())?;
    Ok(Value::BPRel(Rc::new(cochain)))
}

pub(crate) fn zero() -> Value {
    Value::BPRel(Rc::new(Polynomial::default()))
}
