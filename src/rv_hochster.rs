use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::linear::Basis;
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

#[derive(Clone, Copy)]
pub(crate) struct Summand {
    pub(crate) s: usize,
    pub(crate) rank: usize,
}

pub(crate) fn squarefree_ranks(
    model: &mut Model,
    masks: &[usize],
    vertices: usize,
    vertex_limit: u32,
) -> Result<Option<Vec<Summand>>> {
    if !rv_ideal::standard_support(model, masks, 0)? {
        return Ok(Some(Vec::new()));
    }
    let size = vertices.count_ones() as usize;
    let mut reduced = vertices;
    for &mask in masks {
        model.spend(1)?;
        if mask.count_ones() == 1 {
            reduced &= !mask;
        }
    }
    if masks.iter().all(|mask| mask.count_ones() <= 2) {
        reduced = collapse_flag(model, masks, reduced)?;
    }
    let reduced_size = reduced.count_ones();
    if reduced_size > vertex_limit {
        return Ok(None);
    }
    let mut faces = vec![Vec::new(); reduced_size as usize + 1];
    let mut count = 0;
    collect_faces(model, masks, reduced, 0, 0, &mut faces, &mut count)?;
    let mut incoming = 0;
    let mut result = Vec::new();
    for (q, layer) in faces.into_iter().enumerate() {
        let mut differential = Basis::default();
        let mut outgoing = 0;
        for face in &layer {
            let mut image = Vec::new();
            let mut possible = reduced & !face;
            while possible != 0 {
                let bit = 1usize << possible.trailing_zeros();
                possible ^= bit;
                let extension = face | bit;
                if rv_ideal::standard_support(model, masks, extension)? {
                    image.push(extension);
                }
            }
            outgoing += usize::from(differential.insert(model, image)?);
        }
        let rank = layer
            .len()
            .checked_sub(incoming + outgoing)
            .ok_or_else(|| Error::plain("invalid augmented simplicial ranks"))?;
        if rank != 0 {
            result.push(Summand { s: size - q, rank });
        }
        incoming = outgoing;
    }
    result.reverse();
    Ok(Some(result))
}

fn collapse_flag(model: &mut Model, masks: &[usize], mut vertices: usize) -> Result<usize> {
    let mut neighborhoods = vec![vertices; usize::BITS as usize];
    for &mask in masks {
        model.spend(1)?;
        if mask.count_ones() == 2 && (mask & vertices) == mask {
            let lower = mask.trailing_zeros() as usize;
            let upper = mask.ilog2() as usize;
            neighborhoods[lower] &= !(1usize << upper);
            neighborhoods[upper] &= !(1usize << lower);
        }
    }
    loop {
        let mut candidates = vertices;
        let mut removed = 0;
        'search: while candidates != 0 {
            model.spend(1)?;
            let index = candidates.trailing_zeros() as usize;
            let vertex = 1usize << index;
            candidates ^= vertex;
            let closed = neighborhoods[index] & vertices;
            let mut neighbors = closed & !vertex;
            while neighbors != 0 {
                model.spend(1)?;
                let neighbor = neighbors.trailing_zeros() as usize;
                neighbors ^= 1usize << neighbor;
                if closed & !neighborhoods[neighbor] == 0 {
                    // In the clique complex the vertex is dominated by this
                    // neighbor. Deleting it preserves reduced cohomology.
                    removed = vertex;
                    break 'search;
                }
            }
        }
        if removed == 0 {
            return Ok(vertices);
        }
        vertices ^= removed;
    }
}

fn collect_faces(
    model: &mut Model,
    masks: &[usize],
    mut options: usize,
    face: usize,
    q: usize,
    layers: &mut [Vec<usize>],
    count: &mut usize,
) -> Result<()> {
    model.spend(1)?;
    *count += 1;
    if *count > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "Hochster face calculation exceeds the storage limit",
        ));
    }
    layers[q].push(face);
    while options != 0 {
        let bit = 1usize << options.trailing_zeros();
        options ^= bit;
        let extended = face | bit;
        if rv_ideal::standard_support(model, masks, extended)? {
            collect_faces(model, masks, options, extended, q + 1, layers, count)?;
        }
    }
    Ok(())
}
