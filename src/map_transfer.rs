use std::rc::Rc;

use crate::error::{Error, Result};
use crate::map_cochain::{self, Support};
use crate::map_ext::{self, Term};
use crate::model::{Kind, Model, TypeId, Value};

type Input = (TypeId, Rc<[usize]>);

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Product,
    Lift,
}

impl Operation {
    pub(crate) fn from_name(name: &str) -> Option<Self> {
        match name {
            "eHigher" => Some(Self::Product),
            "eLift" => Some(Self::Lift),
            _ => None,
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Product => "eHigher",
            Self::Lift => "eLift",
        }
    }

    pub(crate) fn require_arity(self, arity: usize) -> Result<()> {
        if arity < 2 {
            return Err(Error::plain(format!(
                "{} requires at least two MapExt classes",
                self.name()
            )));
        }
        Ok(())
    }

    pub(crate) fn result_type(
        self,
        model: &mut Model,
        inputs: impl IntoIterator<Item = TypeId>,
        arity: usize,
    ) -> Result<TypeId> {
        self.require_arity(arity)?;
        let name = self.name();
        let mut inputs = inputs.into_iter();
        let outer = Hom::read(model, inputs.next().unwrap(), name)?;
        let mut start = outer.source;
        let mut s = outer.s as i128 - 1;
        let mut t = i128::from(outer.t);
        for ty in inputs {
            let inner = Hom::read(model, ty, name)?;
            if start != inner.target {
                return Err(Error::plain(format!(
                    "{name} requires composable MapExt endpoints"
                )));
            }
            start = inner.source;
            s += inner.s as i128 - 1;
            t += i128::from(inner.t);
        }
        s += match self {
            Self::Product => 2,
            Self::Lift => 1,
        };
        if s < 0 {
            return Err(Error::plain(format!(
                "{name} produces a negative cohomological degree"
            )));
        }
        let s = usize::try_from(s)
            .map_err(|_| Error::plain(format!("{name} cohomological degree is out of range")))?;
        let t = i64::try_from(t)
            .map_err(|_| Error::plain(format!("{name} internal degree is out of range")))?;
        match self {
            Self::Product => model.map_ext_type(start, outer.target, s, t),
            Self::Lift => model.map_cobar_type(start, outer.target, s, t),
        }
    }

    fn zero(self) -> Value {
        match self {
            Self::Product => Value::Class([].into()),
            Self::Lift => Value::MapCobar([].into()),
        }
    }
}

#[derive(Clone, Copy)]
struct Hom {
    source: TypeId,
    target: TypeId,
    s: usize,
    t: i64,
}

impl Hom {
    fn read(model: &Model, ty: TypeId, name: &str) -> Result<Self> {
        let Kind::MapExt(source, target, s, t) = *model.kind(ty) else {
            return Err(Error::plain(format!(
                "{name} requires MapExt classes, found {}",
                model.type_name(ty)
            )));
        };
        Ok(Self {
            source,
            target,
            s,
            t,
        })
    }
}

pub(crate) fn evaluate(model: &mut Model, inputs: &[Input], operation: Operation) -> Result<Value> {
    let arity = inputs.len();
    let ty = operation.result_type(model, inputs.iter().map(|&(ty, _)| ty), arity)?;
    model.spend(arity)?;
    if inputs.iter().any(|(_, class)| class.is_empty()) {
        return Ok(operation.zero());
    }
    if let (Operation::Product, [(outer_ty, outer), (inner_ty, inner)]) = (operation, inputs) {
        return map_ext::compose(model, (*outer_ty, outer), (*inner_ty, inner), ty);
    }

    let prefix_count = arity + 1;
    let interval_count = arity.saturating_mul(prefix_count) / 2;
    model.spend(interval_count)?;
    model.allocate(interval_count.saturating_add(prefix_count))?;
    let mut prefixes = Vec::<(i128, i128)>::with_capacity(prefix_count);
    prefixes.push((0, 0));
    let mut endpoints = Vec::with_capacity(prefix_count);
    let mut singletons = Vec::with_capacity(arity);
    for (index, &(ty, ref class)) in inputs.iter().enumerate() {
        let map = Hom::read(model, ty, operation.name())?;
        let (s, t) = prefixes[index];
        prefixes.push((s + map.s as i128 - 1, t + i128::from(map.t)));
        if index == 0 {
            endpoints.push(map.target);
        }
        endpoints.push(map.source);
        singletons.push(map_cochain::representative(model, ty, class)?);
    }

    // Row k stores F on every consecutive interval of length k + 1.
    // Inputs and interval endpoints run from the outer map to the inner map.
    let mut lifts = Vec::<Vec<Support>>::with_capacity(arity - 1);
    lifts.push(singletons);
    for length in 2..=arity {
        let mut row = Vec::with_capacity(arity - length + 1);
        for start in 0..=arity - length {
            let end = start + length;
            let s = prefixes[end].0 - prefixes[start].0 + 2;
            let t = prefixes[end].1 - prefixes[start].1;
            if length < arity && s < 1 {
                // This lift would be in negative cobar length.
                row.push([].into());
                continue;
            }
            let mut sum = Vec::new();
            for outer_len in 1..length {
                model.spend(1)?;
                let outer = &lifts[outer_len - 1][start];
                let inner = &lifts[length - outer_len - 1][start + outer_len];
                if outer.is_empty() || inner.is_empty() {
                    continue;
                }
                let product = map_ext::compose_terms(model, outer, inner)?;
                sum = if sum.is_empty() {
                    product
                } else {
                    map_ext::add_terms(model, &sum, &product)?
                };
            }
            let endpoints = (endpoints[end], endpoints[start]);
            if length == arity {
                return contract(model, endpoints, (s, t), sum, operation);
            }
            let Value::MapCobar(lift) = contract(model, endpoints, (s, t), sum, Operation::Lift)?
            else {
                unreachable!()
            };
            row.push(lift);
        }
        lifts.push(row);
    }
    unreachable!()
}

fn contract(
    model: &mut Model,
    (domain, range): (TypeId, TypeId),
    (s, t): (i128, i128),
    cochain: Vec<Term>,
    operation: Operation,
) -> Result<Value> {
    if cochain.is_empty() {
        return Ok(operation.zero());
    }
    let s = usize::try_from(s).map_err(|_| intermediate_overflow())?;
    let t = i64::try_from(t).map_err(|_| intermediate_overflow())?;

    // First remove P(d cochain) to obtain a cycle.
    let differential = map_ext::differential_support(model, domain, range, &cochain)?;
    let cycle = if differential.is_empty() {
        cochain
    } else {
        let next_s = s.checked_add(1).ok_or_else(intermediate_overflow)?;
        let ty = model.map_ext_type(domain, range, next_s, t)?;
        let primitive = model
            .map_ext(ty)?
            .primitive(model, &differential)?
            .ok_or_else(|| Error::plain("map cobar differential has no boundary preimage"))?;
        map_ext::add_terms(model, &cochain, &primitive)?
    };
    if cycle.is_empty() {
        return Ok(operation.zero());
    }

    let ty = model.map_ext_type(domain, range, s, t)?;
    let group = model.map_ext(ty)?;
    let class = group.project(model, &cycle)?;
    match operation {
        Operation::Product => {
            model.allocate(class.len())?;
            Ok(Value::Class(class.into()))
        }
        Operation::Lift => {
            let representative = group.representative(model, &class)?;
            let boundary = map_ext::add_terms(model, &cycle, &representative)?;
            let primitive = group.primitive(model, &boundary)?.ok_or_else(|| {
                Error::plain(
                    "map cobar cycle differs from its MapExt representative by a nonboundary",
                )
            })?;
            Ok(Value::MapCobar(map_cochain::persist(model, primitive)?))
        }
    }
}

fn intermediate_overflow() -> Error {
    Error::plain("map transfer intermediate bidegree is out of range")
}
