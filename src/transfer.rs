use std::rc::Rc;

use crate::adams;
use crate::cobar::{self, Word};
use crate::error::{Error, Result};
use crate::model::{Model, Value};

type Bidegree = (usize, usize);
type Input = (Bidegree, Rc<[usize]>);
type Cochain = Rc<[Word]>;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Product,
    Lift,
}

impl Operation {
    pub(crate) fn from_name(name: &str) -> Option<Self> {
        match name {
            "ainfty" => Some(Self::Product),
            "alift" => Some(Self::Lift),
            _ => None,
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Product => "ainfty",
            Self::Lift => "alift",
        }
    }

    pub(crate) fn require_arity(self, arity: usize) -> Result<()> {
        if arity < 2 {
            let name = self.name();
            return Err(Error::plain(format!(
                "{name} requires at least two Ext classes"
            )));
        }
        Ok(())
    }

    pub(crate) fn bidegree(
        self,
        inputs: impl IntoIterator<Item = Bidegree>,
        arity: usize,
    ) -> Result<Bidegree> {
        self.require_arity(arity)?;
        let name = self.name();
        let mut s = 0u128;
        let mut t = 0usize;
        for (length, degree) in inputs {
            s += length as u128;
            t = t
                .checked_add(degree)
                .ok_or_else(|| Error::plain(format!("{name} internal degree is out of range")))?;
        }
        let shift = match self {
            Self::Product => arity - 2,
            Self::Lift => arity - 1,
        };
        let s = s.checked_sub(shift as u128).ok_or_else(|| {
            Error::plain(format!("{name} produces a negative cohomological degree"))
        })?;
        let s = usize::try_from(s)
            .map_err(|_| Error::plain(format!("{name} cohomological degree is out of range")))?;
        Ok((s, t))
    }

    fn zero(self) -> Value {
        match self {
            Self::Product => Value::Class([].into()),
            Self::Lift => cobar::zero(),
        }
    }
}

pub(crate) fn evaluate(model: &mut Model, inputs: &[Input], operation: Operation) -> Result<Value> {
    let arity = inputs.len();
    let (result_s, result_t) =
        operation.bidegree(inputs.iter().map(|&(degree, _)| degree), arity)?;
    model.spend(arity)?;
    if inputs.iter().any(|(_, coordinates)| coordinates.is_empty()) {
        return Ok(operation.zero());
    }
    // Ext(0,0) is spanned by the strict unit, and other Ext(0,t) vanish.
    if let Some(unit) = inputs.iter().position(|((s, _), _)| *s == 0) {
        return Ok(match (operation, arity) {
            (Operation::Product, 2) => Value::Class(inputs[1 - unit].1.clone()),
            _ => operation.zero(),
        });
    }
    if let (Operation::Product, [((s, t), left), ((u, v), right)]) = (operation, inputs) {
        let left_ty = model.ext_type(*s, *t)?;
        let right_ty = model.ext_type(*u, *v)?;
        let ty = model.ext_type(result_s, result_t)?;
        return adams::yoneda(model, left_ty, left, right_ty, right, ty);
    }

    let prefix_count = arity.saturating_add(1);
    let interval_count = arity.saturating_mul(prefix_count) / 2;
    model.spend(interval_count)?;
    model.allocate(interval_count.saturating_add(prefix_count))?;
    let mut prefixes: Vec<Bidegree> = Vec::with_capacity(prefix_count);
    prefixes.push((0, 0));
    let mut singletons: Vec<Cochain> = Vec::with_capacity(arity);
    for &((s, t), ref coordinates) in inputs {
        let &(prev_s, prev_t) = prefixes.last().unwrap();
        let sum_s = prev_s
            .checked_add(s - 1)
            .ok_or_else(intermediate_overflow)?;
        let sum_t = prev_t.checked_add(t).ok_or_else(intermediate_overflow)?;
        prefixes.push((sum_s, sum_t));
        let ty = model.ext_type(s, t)?;
        let adams = model.adams(ty)?;
        let words = adams.representative(model, coordinates)?;
        model.allocate(words.len())?;
        singletons.push(Rc::from(words));
    }

    // Row k stores F on the consecutive intervals of length k + 1.
    let mut lifts = Vec::with_capacity(arity - 1);
    lifts.push(singletons);
    for length in 2..=arity {
        let mut row = Vec::with_capacity(arity - length + 1);
        for start in 0..=arity - length {
            let end = start + length;
            let s = (prefixes[end].0 - prefixes[start].0)
                .checked_add(2)
                .ok_or_else(intermediate_overflow)?;
            let bidegree = (s, prefixes[end].1 - prefixes[start].1);
            let mut source: Cochain = [].into();
            for left_len in 1..length {
                model.spend(1)?;
                let right_len = length - left_len;
                let left = &lifts[left_len - 1][start];
                let right = &lifts[right_len - 1][start + left_len];
                if left.is_empty() || right.is_empty() {
                    continue;
                }
                let term = cobar::concat_support(model, left, right)?;
                source = if source.is_empty() {
                    term.into()
                } else {
                    add(model, &source, &term)?
                };
            }
            if length == arity {
                return contract(model, bidegree, source, operation);
            }
            let Value::Cobar(lift) = contract(model, bidegree, source, Operation::Lift)? else {
                unreachable!();
            };
            row.push(lift);
        }
        lifts.push(row);
    }
    unreachable!()
}

fn intermediate_overflow() -> Error {
    Error::plain("transfer intermediate cobar bidegree is out of range")
}

fn add(model: &mut Model, left: &Cochain, right: &[Word]) -> Result<Cochain> {
    if right.is_empty() {
        return Ok(left.clone());
    }
    let Value::Cobar(sum) = cobar::add(model, left, right)? else {
        unreachable!();
    };
    Ok(sum)
}

fn contract(
    model: &mut Model,
    (s, t): Bidegree,
    source: Cochain,
    operation: Operation,
) -> Result<Value> {
    if source.is_empty() {
        return Ok(operation.zero());
    }

    // q = source + P(d source) is closed even when source is not.
    let differential = cobar::differential_support(model, &source)?;
    let closed = if differential.is_empty() {
        source
    } else {
        let next_s = s
            .checked_add(1)
            .ok_or_else(|| Error::plain("cobar differential degree is out of range"))?;
        let ty = model.ext_type(next_s, t)?;
        let next = model.adams(ty)?;
        let primitive = next
            .primitive(model, &differential)?
            .ok_or_else(|| Error::plain("cobar differential has no boundary preimage"))?;
        add(model, &source, &primitive)?
    };
    if closed.is_empty() {
        return Ok(operation.zero());
    }

    let ty = model.ext_type(s, t)?;
    let adams = model.adams(ty)?;
    let coordinates = adams.project(model, &closed)?;
    match operation {
        Operation::Product => {
            model.allocate(coordinates.len())?;
            Ok(Value::Class(coordinates.into()))
        }
        Operation::Lift => {
            let representative = adams.representative(model, &coordinates)?;
            let boundary = add(model, &closed, &representative)?;
            let primitive = adams.primitive(model, &boundary)?.ok_or_else(|| {
                Error::plain("cobar cycle differs from its Ext representative by a nonboundary")
            })?;
            model.allocate(primitive.len())?;
            Ok(Value::Cobar(primitive.into()))
        }
    }
}
