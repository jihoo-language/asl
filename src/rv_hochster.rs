use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_ideal;
use crate::rv_koszul::{self, Grade, Support, Term};
use crate::syntax::Label;

fn squarefree_error() -> Error {
    Error::plain("Hochster operations require a squarefree monomial rank ideal")
}

fn indicator(width: usize, mask: usize) -> Vec<usize> {
    (0..width).map(|index| (mask >> index) & 1).collect()
}

pub(crate) fn face(
    model: &mut Model,
    name: Label,
    vertices: &[usize],
    face: &[usize],
) -> Result<(Grade, Support)> {
    let ideal = model.rv_ideal(&name)?;
    let initial = ideal
        .squarefree_masks(model)?
        .ok_or_else(squarefree_error)?;
    model.spend(1 + vertices.len() + face.len())?;
    let support = rv_koszul::mask(&ideal, vertices)?;
    let face_mask = rv_koszul::mask(&ideal, face)?;
    if face_mask & !support != 0 {
        return Err(Error::plain(
            "rvhface requires a face contained in the induced vertex set",
        ));
    }
    let width = ideal.selected.len();
    let powers = indicator(width, face_mask);
    if !rv_ideal::standard_support(model, initial, face_mask)? {
        return Err(Error::plain(
            "rvhface subset is not a face of the Stanley–Reisner complex",
        ));
    }
    let grade = Grade::local(
        model,
        &ideal,
        name.text.into(),
        vertices.len() - face.len(),
        indicator(width, support),
    )?;
    let terms = rv_koszul::persist(
        model,
        vec![Term {
            wedge: support ^ face_mask,
            powers: powers.into(),
        }],
    )?;
    Ok((grade, terms))
}

pub(crate) fn rank(model: &mut Model, name: &Label, vertices: &[usize], q: i64) -> Result<usize> {
    let ideal = model.rv_ideal(name)?;
    ideal
        .squarefree_masks(model)?
        .ok_or_else(squarefree_error)?;
    model.spend(1 + vertices.len())?;
    let support = rv_koszul::mask(&ideal, vertices)?;
    // The augmented simplicial complex begins at q = -1, with the empty face.
    let n = vertices.len() as i128;
    let s = n - i128::from(q) - 1;
    if s < 0 || s > n {
        return Ok(0);
    }
    let grade = Grade::local(
        model,
        &ideal,
        name.text.as_str().into(),
        s as usize,
        indicator(ideal.selected.len(), support),
    )?;
    let ty = model.rv_complex_type(grade, true)?;
    Ok(model.rv_tor(ty)?.rank())
}
