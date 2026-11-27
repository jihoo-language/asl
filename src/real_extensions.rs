use std::rc::Rc;

use crate::error::{Error, Result};
use crate::linear::Basis;
use crate::model::{Kind, Model, TypeId, Value};
use crate::motivic::Grade;
use crate::real_bars;
use crate::real_ext;
use crate::real_motivic::{self, Term};
use crate::rho_bock;

pub(crate) type Input<'a> = (TypeId, &'a [usize]);

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Coordinates,
    Class,
    Lift,
    Divide(usize),
    Edge(Option<usize>, usize),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Coordinates => "rmbarCoords",
            Self::Class => "rmbarClass",
            Self::Lift => "rmbarLift",
            Self::Divide(_) => "rhoDivide",
            Self::Edge(..) => "rhoEinf",
        }
    }
}

pub(crate) fn bars(model: &mut Model, grade: Grade) -> Result<Rc<real_bars::Group>> {
    let k = grade
        .t
        .checked_sub(i128::from(grade.w))
        .ok_or_else(grade_overflow)?;
    model.real_bars(grade.s, k)
}

fn grade_overflow() -> Error {
    Error::plain("rho extension grade exceeds the supported range")
}

fn grade(model: &Model, ty: TypeId) -> Grade {
    let Kind::RMExt(grade) = *model.kind(ty) else {
        unreachable!()
    };
    grade
}

fn coordinates(model: &mut Model, (ty, class): Input<'_>) -> Result<Vec<usize>> {
    if class.is_empty() {
        return Ok(Vec::new());
    }
    let grade = grade(model, ty);
    let bars = bars(model, grade)?;
    let ordinary = model.real_ext(ty)?;
    let mut change = Basis::default();
    for index in bars.living(model, grade.t)? {
        let z = bars.get(index)?.lift_at(model, grade.t)?;
        let column = ordinary.project(model, &z)?;
        if change
            .insert_with_witness(model, column, vec![index])?
            .is_some()
        {
            return Err(Error::plain("living rho bars are linearly dependent"));
        }
    }
    if change.rank() != ordinary.rank() {
        return Err(Error::plain(
            "living rho bars do not span the RMExt component",
        ));
    }
    let (remainder, support) = change.reduce_with_witness(model, class.to_vec(), Vec::new())?;
    if !remainder.is_empty() {
        return Err(Error::plain("RMExt class is outside its rho bar span"));
    }
    Ok(support)
}

pub(crate) fn decompose(model: &mut Model, input: Input<'_>) -> Result<Value> {
    let support = coordinates(model, input)?;
    model.allocate(support.len())?;
    Ok(Value::Class(support.into()))
}

fn raw_lift(model: &mut Model, grade: Grade, support: &[usize]) -> Result<Vec<Term>> {
    if support.is_empty() {
        return Ok(Vec::new());
    }
    let bars = bars(model, grade)?;
    let mut terms = Vec::new();
    for &i in support {
        let lift = bars.get(i)?.lift_at(model, grade.t)?;
        terms = if terms.is_empty() {
            lift
        } else {
            real_motivic::add(model, &terms, &lift)?
        };
    }
    Ok(terms)
}

pub(crate) fn lift(model: &mut Model, grade: Grade, support: &[usize]) -> Result<Value> {
    let terms = raw_lift(model, grade, support)?;
    real_motivic::persist(model, terms)
}

pub(crate) fn reconstruct(
    model: &mut Model,
    grade: Grade,
    support: &[usize],
    ty: TypeId,
) -> Result<Value> {
    let terms = raw_lift(model, grade, support)?;
    real_ext::class(model, ty, &terms)
}

pub(crate) fn depth_of_support(
    model: &mut Model,
    grade: Grade,
    support: &[usize],
) -> Result<Option<usize>> {
    let Some(&first) = support.first() else {
        return Ok(None);
    };
    // Bar indices are sorted first by birth. For a living bar of birth B
    // in degree t the coefficient is precisely rho^(B-t).
    Ok(Some(bars(model, grade)?.get(first)?.exponent_at(grade.t)?))
}

pub(crate) fn depth(model: &mut Model, input: Input<'_>) -> Result<Option<usize>> {
    let support = coordinates(model, input)?;
    depth_of_support(model, grade(model, input.0), &support)
}

pub(crate) fn quotient_grade(grade: Grade, power: usize) -> Result<Grade> {
    let t = grade
        .t
        .checked_add(power as i128)
        .ok_or_else(grade_overflow)?;
    let w = i64::try_from(i128::from(grade.w) + power as i128).map_err(|_| grade_overflow())?;
    Ok(Grade { t, w, ..grade })
}

fn requested_quotient(
    model: &mut Model,
    input: Input<'_>,
    power: usize,
) -> Result<(Grade, Vec<usize>)> {
    let grade = grade(model, input.0);
    let quotient = quotient_grade(grade, power)?;
    let support = coordinates(model, input)?;
    if let Some(depth) = depth_of_support(model, grade, &support)? {
        if power > depth {
            return Err(Error::plain(format!(
                "class has rho depth {depth} and cannot be divided by rho^{power}"
            )));
        }
    }
    Ok((quotient, support))
}

pub(crate) fn divide(
    model: &mut Model,
    input: Input<'_>,
    power: usize,
    output: TypeId,
) -> Result<Value> {
    if power == 0 || input.1.is_empty() {
        model.allocate(input.1.len())?;
        return Ok(Value::Class(input.1.into()));
    }
    let (quotient, support) = requested_quotient(model, input, power)?;
    let Value::Class(coords) = reconstruct(model, quotient, &support, output)? else {
        unreachable!()
    };
    let Value::Class(shifted) = real_ext::shift(model, output, &coords, input.0, power)? else {
        unreachable!()
    };
    if shifted.as_ref() != input.1 {
        return Err(Error::plain(
            "chosen rho quotient does not shift to the original class",
        ));
    }
    Ok(Value::Class(coords))
}

pub(crate) fn edge_grade(
    grade: Grade,
    page: Option<usize>,
    filtration: usize,
) -> Result<rho_bock::Grade> {
    let bound = if grade.t < i128::from(grade.w) {
        0
    } else {
        grade
            .t
            .checked_sub(2 * i128::from(grade.w))
            .ok_or_else(grade_overflow)?
            .max(0)
    };
    let stable = usize::try_from(bound)
        .ok()
        .and_then(|f| f.checked_add(1))
        .ok_or_else(|| Error::plain("rho-Bockstein stable page is outside the supported range"))?;
    let r = page.unwrap_or(stable);
    if r < stable {
        return Err(Error::plain(format!(
            "rhoEinf requires page {stable} or later for this tridegree"
        )));
    }
    Ok(rho_bock::Grade {
        r,
        chain: grade,
        f: filtration,
    })
}

pub(crate) fn edge(
    model: &mut Model,
    input: Input<'_>,
    output: TypeId,
    filtration: usize,
) -> Result<Value> {
    let support = coordinates(model, input)?;
    let grade = grade(model, input.0);
    let Some(depth) = depth_of_support(model, grade, &support)? else {
        return Ok(Value::Class([].into()));
    };
    if depth < filtration {
        return Err(Error::plain(format!(
            "rhoEinf requires depth at least {filtration}, found {depth}"
        )));
    }
    if depth > filtration {
        return Ok(Value::Class([].into()));
    }
    let terms = raw_lift(model, grade, &support)?;
    let Value::Class(edge) = rho_bock::class(model, output, &terms)? else {
        unreachable!()
    };
    if edge.is_empty() {
        return Err(Error::plain(
            "leading rho bar vanished on the stable Bockstein page",
        ));
    }
    Ok(Value::Class(edge))
}

pub(crate) fn product_null(
    model: &mut Model,
    left: Input<'_>,
    right: Input<'_>,
    power: usize,
) -> Result<Value> {
    if left.1.is_empty() || right.1.is_empty() {
        return Ok(real_motivic::zero());
    }
    let g = grade(model, left.0).product(grade(model, right.0))?;
    let a = model.real_ext(left.0)?.representative(model, left.1)?;
    let b = model.real_ext(right.0)?.representative(model, right.1)?;
    let product = real_motivic::concat_support(model, &a, &b)?;
    let ty = model.rmext_type(g)?;
    let Value::Class(class) = real_ext::class(model, ty, &product)? else {
        unreachable!()
    };
    let (quotient, support) = requested_quotient(model, (ty, &class), power)?;
    let rep = raw_lift(model, quotient, &support)?;
    let shifted = real_motivic::shift_support(model, &rep, power)?;
    let difference = real_motivic::add(model, &product, &shifted)?;
    real_ext::primitive(model, g, &difference)
}

pub(crate) fn hidden(model: &mut Model, left: Input<'_>, right: Input<'_>) -> Result<String> {
    let g = grade(model, left.0).product(grade(model, right.0))?;
    let ty = model.rmext_type(g)?;
    let Value::Class(class) = real_ext::yoneda(model, left, right, ty)? else {
        unreachable!()
    };
    let Some(actual) = depth(model, (ty, &class))? else {
        return Ok("rhohidden = zero".to_owned());
    };
    let (Some(a), Some(b)) = (depth(model, left)?, depth(model, right)?) else {
        return Err(Error::plain(
            "nonzero Yoneda product has a vanishing factor",
        ));
    };
    let expected = a.checked_add(b).ok_or_else(grade_overflow)?;
    if actual < expected {
        return Err(Error::plain(
            "Yoneda product does not preserve rho filtration",
        ));
    }
    Ok(if actual > expected {
        format!("rhohidden = hidden ({a} + {b} < {actual})")
    } else {
        format!("rhohidden = visible ({a} + {b} = {actual})")
    })
}
