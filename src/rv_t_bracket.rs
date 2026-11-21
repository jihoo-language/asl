use std::rc::Rc;

use crate::error::{Error, Result};
use crate::linear::{self, Basis};
use crate::model::{Kind, Model, TypeId, Value};
use crate::rv_ideal::Ideal;
use crate::rv_koszul::{self, Grade, Term};

type Input = (TypeId, Rc<[usize]>);

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Predicate,
    Membership,
}

impl Operation {
    pub(crate) fn from_name(name: &str) -> Option<Self> {
        match name {
            "rvTBracket" => Some(Self::Predicate),
            "rvTInBracket" => Some(Self::Membership),
            _ => None,
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Predicate => "rvTBracket",
            Self::Membership => "rvTInBracket",
        }
    }

    pub(crate) fn require_arity(self, arity: usize) -> Result<()> {
        match (self, arity) {
            (Self::Predicate, 3) | (Self::Membership, 4) => Ok(()),
            (Self::Predicate, _) => Err(Error::plain(
                "rvTBracket requires exactly three RVTor arguments",
            )),
            (Self::Membership, _) => Err(Error::plain(
                "rvTInBracket requires exactly four RVTor arguments",
            )),
        }
    }
}

fn grade(model: &Model, ty: TypeId, name: &str) -> Result<Grade> {
    let Kind::RVTor(grade) = model.kind(ty) else {
        return Err(Error::plain(format!(
            "{name} requires RVTor classes, found {}",
            model.type_name(ty)
        )));
    };
    Ok(grade.clone())
}

pub(crate) fn result_type(model: &mut Model, [a, b, c]: [TypeId; 3], name: &str) -> Result<TypeId> {
    let a = grade(model, a, name)?;
    let b = grade(model, b, name)?;
    let c = grade(model, c, name)?;
    if a.ideal != b.ideal || b.ideal != c.ideal {
        return Err(Error::plain(format!("{name} requires the same rank ideal")));
    }
    let target = a.product(&b)?.product(&c)?.primitive()?;
    model.rv_complex_type(target, true)
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
    if query == system.base {
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
    Ok(Value::Point(usize::from(
        directions.reduce(model, difference)?.is_empty(),
    )))
}

struct DefiningSystem {
    base: Vec<usize>,
    ideal: Rc<Ideal>,
    a: Vec<Term>,
    c: Vec<Term>,
    ab: Grade,
    bc: Grade,
}

impl DefiningSystem {
    fn new(model: &mut Model, target: TypeId, [a, b, c]: &[Input; 3]) -> Result<Option<Self>> {
        let grade_a = grade(model, a.0, "Koszul triple Massey product")?;
        let grade_b = grade(model, b.0, "Koszul triple Massey product")?;
        let grade_c = grade(model, c.0, "Koszul triple Massey product")?;
        let ab = grade_a.product(&grade_b)?;
        let bc = grade_b.product(&grade_c)?;
        let ideal = model.rv_ideal_named(&ab.ideal).unwrap();

        let a = representative(model, a)?;
        let c = representative(model, c)?;
        let b = if a.is_empty() && c.is_empty() {
            Vec::new()
        } else {
            representative(model, b)?
        };
        let Some(u) = solve(model, &ideal, &ab, &a, &b, !c.is_empty())? else {
            return Ok(None);
        };
        let Some(v) = solve(model, &ideal, &bc, &b, &c, !a.is_empty())? else {
            return Ok(None);
        };

        let av = rv_koszul::product_terms(model, &ideal, &a, &v)?;
        let uc = rv_koszul::product_terms(model, &ideal, &u, &c)?;
        let cycle = rv_koszul::sum_terms(model, &av, &uc)?;
        let base = if cycle.is_empty() {
            Vec::new()
        } else {
            model.rv_tor(target)?.project(model, &ideal, &cycle)?
        };
        Ok(Some(Self {
            base,
            ideal,
            a,
            c,
            ab,
            bc,
        }))
    }

    fn indeterminacy(&self, model: &mut Model, target: TypeId) -> Result<Basis> {
        let mut directions = Basis::default();
        for (pair, outside) in [(&self.ab, &self.c), (&self.bc, &self.a)] {
            if outside.is_empty() {
                continue;
            }
            let ambient = model.rv_tor(target)?;
            if ambient.rank() == 0 {
                break;
            }
            let ty = model.rv_complex_type(pair.primitive()?, true)?;
            let group = model.rv_tor(ty)?;
            for index in 0..group.rank() {
                model.spend(1)?;
                let cycle = group.representative_terms(model, &[index])?;
                // Characteristic two makes the two multiplication orders agree.
                let product = rv_koszul::product_terms(model, &self.ideal, &cycle, outside)?;
                if !product.is_empty() {
                    let class = ambient.project(model, &self.ideal, &product)?;
                    directions.insert(model, class)?;
                }
            }
        }
        Ok(directions)
    }
}

fn representative(model: &mut Model, (ty, class): &Input) -> Result<Vec<Term>> {
    if class.is_empty() {
        return Ok(Vec::new());
    }
    model.rv_tor(*ty)?.representative_terms(model, class)
}

fn solve(
    model: &mut Model,
    ideal: &Ideal,
    grade: &Grade,
    left: &[Term],
    right: &[Term],
    need_primitive: bool,
) -> Result<Option<Vec<Term>>> {
    let product = rv_koszul::product_terms(model, ideal, left, right)?;
    if product.is_empty() {
        return Ok(Some(Vec::new()));
    }
    let ty = model.rv_complex_type(grade.clone(), true)?;
    let group = model.rv_tor(ty)?;
    if need_primitive {
        group.primitive(model, ideal, grade, &product)
    } else {
        Ok(group.is_boundary(model, ideal, &product)?.then(Vec::new))
    }
}
