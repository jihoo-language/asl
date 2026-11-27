use crate::adams;
use crate::cobar;
use crate::cobar_squares;
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::motivic::Grade;
use crate::profile::Profile;
use crate::real_ext;
use crate::real_motivic::{self, Term};
use crate::real_rho;
use crate::rho_bock;

fn grade_limit() -> Error {
    Error::plain("rhoFrob doubled grade is out of range")
}

fn doubled(grade: Grade) -> Result<Grade> {
    Ok(Grade {
        t: grade.t.checked_mul(2).ok_or_else(grade_limit)?,
        w: grade.w.checked_mul(2).ok_or_else(grade_limit)?,
        ..grade
    })
}

pub(crate) fn truncation(grade: real_rho::Grade) -> Result<real_rho::Grade> {
    Ok(real_rho::Grade {
        chain: doubled(grade.chain)?,
        exponent: grade.exponent.checked_mul(2).ok_or_else(grade_limit)?,
    })
}

pub(crate) fn target(model: &mut Model, input: TypeId) -> Result<Option<(TypeId, Grade)>> {
    let kind = model.kind(input);
    let (ty, chain) = match *kind {
        Kind::RMCobar(chain) => {
            let chain = doubled(chain)?;
            (model.real_cobar_type(chain)?, chain)
        }
        Kind::RMExt(chain) => {
            let chain = doubled(chain)?;
            (model.rmext_type(chain)?, chain)
        }
        Kind::MCobar(chain) => {
            let chain = doubled(chain)?;
            (model.mcobar_type(chain)?, chain)
        }
        Kind::MExt(chain) => {
            let chain = doubled(chain)?;
            (model.mext_type(chain)?, chain)
        }
        Kind::RMRhoCobar(grade) => {
            let grade = truncation(grade)?;
            (model.real_rho_cobar_type(grade)?, grade.chain)
        }
        Kind::RMRhoExt(grade) => {
            let grade = truncation(grade)?;
            (model.real_rho_ext_type(grade)?, grade.chain)
        }
        Kind::RhoBock(grade) => {
            let chain = doubled(grade.chain)?;
            let grade = rho_bock::Grade {
                chain,
                r: grade.r.checked_mul(2).ok_or_else(grade_limit)?,
                f: grade.f.checked_mul(2).ok_or_else(grade_limit)?,
            };
            (model.rho_page_type(grade)?, chain)
        }
        _ => return Ok(None),
    };
    Ok(Some((ty, chain)))
}

pub(crate) fn apply(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    value: Value,
) -> Result<Value> {
    let cutoff = match *model.kind(output) {
        Kind::RMRhoExt(grade) | Kind::RMRhoCobar(grade) => grade.exponent as u128,
        _ => u128::MAX,
    };
    match value {
        Value::RMCobar(terms) => {
            let terms = real_motivic::frobenius_below(model, &terms, cutoff)?;
            real_motivic::persist(model, terms)
        }
        Value::Cobar(words) => {
            let words = cobar_squares::frobenius(model, Profile::Full, &words)?;
            model.allocate(words.iter().map(cobar::word_cells).sum())?;
            Ok(Value::Cobar(words.into()))
        }
        Value::Class(coords) if coords.is_empty() => Ok(Value::Class(coords)),
        Value::Class(coords) => match model.kind(input) {
            Kind::RhoBock(..) => rho_bock::frobenius(model, input, output, &coords),
            Kind::MExt(..) => {
                let lift = model.adams(input)?.representative(model, &coords)?;
                let words = cobar_squares::frobenius(model, Profile::Full, &lift)?;
                adams::class(model, output, &words)
            }
            Kind::RMExt(..) | Kind::RMRhoExt(..) => {
                let lift = model.real_ext(input)?.representative(model, &coords)?;
                let terms = real_motivic::frobenius_below(model, &lift, cutoff)?;
                real_ext::class(model, output, &terms)
            }
            _ => unreachable!(),
        },
        _ => unreachable!(),
    }
}

pub(crate) fn connecting_null(model: &mut Model, input: TypeId, value: Value) -> Result<Value> {
    let grade = match *model.kind(input) {
        Kind::RMRhoExt(grade) | Kind::RMRhoCobar(grade) => grade,
        _ => unreachable!(),
    };
    let target = truncation(grade)?;
    let difference = match value {
        Value::RMCobar(lift) => {
            if !real_motivic::differential_below(model, &lift, grade.exponent as u128)?.is_empty() {
                return Err(Error::plain(
                    "rhoFrobConnectNull requires a closed rho-truncated cochain",
                ));
            }
            let mut image = real_motivic::frobenius_below(model, &lift, u128::MAX)?;
            image.split_off(real_rho::cut(&image, target.exponent).len())
        }
        Value::Class(coords) if coords.is_empty() => return Ok(real_motivic::zero()),
        Value::Class(coords) => {
            let lift = model.real_ext(input)?.representative(model, &coords)?;
            let image = real_motivic::frobenius_below(model, &lift, u128::MAX)?;
            canonical_difference(model, target, image)?
        }
        _ => unreachable!(),
    };
    let mut divided = difference;
    real_rho::divide_rho(model, &mut divided, target.exponent, || {
        Error::plain("rhoFrob connecting lift difference is not divisible by the target cutoff")
    })?;
    real_motivic::persist(model, divided)
}

fn canonical_difference(
    model: &mut Model,
    target: real_rho::Grade,
    image: Vec<Term>,
) -> Result<Vec<Term>> {
    let visible = real_rho::cut(&image, target.exponent);
    if visible.is_empty() {
        return Ok(image);
    }
    let ty = model.real_rho_ext_type(target)?;
    let group = model.real_ext(ty)?;
    let coords = group.project(model, visible)?;
    let chosen = group.representative(model, &coords)?;
    let mismatch = real_motivic::add(model, visible, &chosen)?;
    // Fresh representatives may differ by a truncated boundary. Correct
    // that boundary in the real lift before dividing the difference.
    let primitive = real_rho::boundary_witness(model, target, &mismatch, |_| Ok(None))?
        .ok_or_else(|| Error::plain("rhoFrob canonical lifts disagree in truncated cohomology"))?;
    let correction = real_motivic::differential_support(model, &primitive)?;
    let difference = real_motivic::add(model, &image, &chosen)?;
    real_motivic::add(model, &difference, &correction)
}
