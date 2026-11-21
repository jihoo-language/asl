use std::rc::Rc;

use crate::error::{Error, Result};
use crate::linear::{self, Basis};
use crate::map_cochain;
use crate::map_ext::{self, Term};
use crate::model::{Kind, Model, TypeId, Value};

type Input = (TypeId, Rc<[usize]>);

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Predicate,
    Membership,
}

impl Operation {
    pub(crate) fn from_name(name: &str) -> Option<Self> {
        match name {
            "eBracket" => Some(Self::Predicate),
            "eInBracket" => Some(Self::Membership),
            _ => None,
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Predicate => "eBracket",
            Self::Membership => "eInBracket",
        }
    }

    pub(crate) fn require_arity(self, arity: usize) -> Result<()> {
        let (expected, message) = match self {
            Self::Predicate => (3, "eBracket requires three MapExt arguments"),
            Self::Membership => (4, "eInBracket requires four MapExt arguments"),
        };
        if arity == expected {
            Ok(())
        } else {
            Err(Error::plain(message))
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

pub(crate) fn result_type(model: &mut Model, [h, g, f]: [TypeId; 3], name: &str) -> Result<TypeId> {
    let h = Hom::read(model, h, name)?;
    let g = Hom::read(model, g, name)?;
    let f = Hom::read(model, f, name)?;
    if g.target != h.source || f.target != g.source {
        return Err(Error::plain(format!(
            "{name} requires composable MapExt endpoints"
        )));
    }
    let filtration = f.s as u128 + g.s as u128 + h.s as u128;
    let s = filtration
        .checked_sub(1)
        .ok_or_else(|| Error::plain(format!("{name} produces a negative cohomological degree")))?;
    let s = usize::try_from(s)
        .map_err(|_| Error::plain(format!("{name} cohomological degree is out of range")))?;
    let internal = i128::from(f.t) + i128::from(g.t) + i128::from(h.t);
    let t = i64::try_from(internal)
        .map_err(|_| Error::plain(format!("{name} internal degree is out of range")))?;
    model.map_ext_type(f.source, h.target, s, t)
}

pub(crate) fn predicate(model: &mut Model, target: TypeId, inputs: &[Input; 3]) -> Result<Value> {
    let points = model.elements(target)?;
    model.allocate(1 + points.len())?;
    let Some(system) = DefiningSystem::new(model, target, inputs)? else {
        return Ok(Value::Map(vec![Value::Point(0); points.len()].into()));
    };
    let directions = system.indeterminacy(model, target)?;
    let rows = points
        .iter()
        .map(|point| {
            let Value::Class(query) = point else {
                unreachable!();
            };
            affine_member(model, &system.base, &directions, query)
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(Value::Map(rows.into()))
}

pub(crate) fn member(
    model: &mut Model,
    target: TypeId,
    inputs: &[Input; 3],
    query: &[usize],
) -> Result<Value> {
    let Some(system) = DefiningSystem::new(model, target, inputs)? else {
        return Ok(Value::Point(0));
    };
    model.spend(system.base.len() + query.len())?;
    if query == system.base.as_slice() {
        return Ok(Value::Point(1));
    }
    let directions = system.indeterminacy(model, target)?;
    affine_member(model, &system.base, &directions, query)
}

fn affine_member(
    model: &mut Model,
    base: &[usize],
    directions: &Basis,
    query: &[usize],
) -> Result<Value> {
    let difference = linear::xor(model, base, query)?;
    let inside = directions.reduce(model, difference)?.is_empty();
    Ok(Value::Point(usize::from(inside)))
}

#[derive(Clone, Copy)]
struct Pair {
    source: TypeId,
    target: TypeId,
    s: u128,
    t: i128,
}

impl Pair {
    fn new(outer: Hom, inner: Hom) -> Self {
        Self {
            source: inner.source,
            target: outer.target,
            s: inner.s as u128 + outer.s as u128,
            t: i128::from(inner.t) + i128::from(outer.t),
        }
    }

    fn at(self, model: &mut Model, filtration: u128) -> Result<TypeId> {
        let s = usize::try_from(filtration).map_err(|_| intermediate_overflow())?;
        let t = i64::try_from(self.t).map_err(|_| intermediate_overflow())?;
        model.map_ext_type(self.source, self.target, s, t)
    }

    fn primitive(
        self,
        model: &mut Model,
        outer: &[Term],
        inner: &[Term],
        need_lift: bool,
    ) -> Result<Option<Vec<Term>>> {
        let composition = map_ext::compose_terms(model, outer, inner)?;
        if composition.is_empty() {
            return Ok(Some(Vec::new()));
        }
        if self.s == 0 {
            // There are no incoming cobar boundaries at length zero.
            return Ok(None);
        }
        let ty = self.at(model, self.s)?;
        let group = model.map_ext(ty)?;
        if need_lift {
            group.primitive(model, &composition)
        } else {
            Ok(group.is_boundary(model, &composition)?.then(Vec::new))
        }
    }

    fn add_directions(
        self,
        model: &mut Model,
        target: TypeId,
        directions: &mut Basis,
        mut product: impl FnMut(&mut Model, &[Term]) -> Result<Vec<Term>>,
    ) -> Result<()> {
        let Some(s) = self.s.checked_sub(1) else {
            return Ok(());
        };
        let ty = self.at(model, s)?;
        let group = model.map_ext(ty)?;
        for index in 0..group.rank() {
            model.spend(1)?;
            let cycle = group.representative(model, &[index])?;
            let terms = product(model, &cycle)?;
            if !terms.is_empty() {
                let coordinates = model.map_ext(target)?.project(model, &terms)?;
                directions.insert(model, coordinates)?;
            }
        }
        Ok(())
    }
}

fn intermediate_overflow() -> Error {
    Error::plain("triple Massey intermediate bidegree is out of range")
}

struct DefiningSystem {
    base: Vec<usize>,
    h: Vec<Term>,
    f: Vec<Term>,
    hg: Pair,
    gf: Pair,
}

impl DefiningSystem {
    fn new(model: &mut Model, target: TypeId, [h, g, f]: &[Input; 3]) -> Result<Option<Self>> {
        let h_type = Hom::read(model, h.0, "triple Massey product")?;
        let g_type = Hom::read(model, g.0, "triple Massey product")?;
        let f_type = Hom::read(model, f.0, "triple Massey product")?;
        let hg = Pair::new(h_type, g_type);
        let gf = Pair::new(g_type, f_type);

        let h = representative(model, h)?;
        let f = representative(model, f)?;
        let g = if h.is_empty() && f.is_empty() {
            Vec::new()
        } else {
            representative(model, g)?
        };
        let Some(u) = gf.primitive(model, &g, &f, !h.is_empty())? else {
            return Ok(None);
        };
        let Some(v) = hg.primitive(model, &h, &g, !f.is_empty())? else {
            return Ok(None);
        };
        let hu = map_ext::compose_terms(model, &h, &u)?;
        let vf = map_ext::compose_terms(model, &v, &f)?;
        let cycle = map_cochain::add(model, &hu, &vf)?;
        let base = if cycle.is_empty() {
            Vec::new()
        } else {
            model.map_ext(target)?.project(model, &cycle)?
        };
        Ok(Some(Self { base, h, f, hg, gf }))
    }

    fn indeterminacy(&self, model: &mut Model, target: TypeId) -> Result<Basis> {
        let mut directions = Basis::default();
        if !self.h.is_empty() {
            self.gf
                .add_directions(model, target, &mut directions, |model, cycle| {
                    map_ext::compose_terms(model, &self.h, cycle)
                })?;
        }
        if !self.f.is_empty() {
            self.hg
                .add_directions(model, target, &mut directions, |model, cycle| {
                    map_ext::compose_terms(model, cycle, &self.f)
                })?;
        }
        Ok(directions)
    }
}

fn representative(model: &mut Model, (ty, class): &Input) -> Result<Vec<Term>> {
    if class.is_empty() {
        return Ok(Vec::new());
    }
    model.map_ext(*ty)?.representative(model, class)
}
