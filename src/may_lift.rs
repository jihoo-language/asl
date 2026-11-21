use crate::cobar;
use crate::error::{Error, Result};
use crate::linear::Basis;
use crate::map_cochain;
use crate::map_ext;
use crate::map_may::MappingComplex;
use crate::may::{self, FilteredComplex, Rows, Sphere};
use crate::model::{Kind, Model, TypeId, Value};

pub(crate) fn evaluate(
    model: &mut Model,
    ty: TypeId,
    r: usize,
    f: usize,
    input: Value,
) -> Result<Value> {
    match input {
        Value::Cobar(words) => {
            let (profile, s, t) = model.kind(ty).cobar_grade().unwrap();
            let mut complex = Sphere::new(profile, t);
            let Some(correction) = correction(model, &mut complex, r, s, f, &words)? else {
                return Ok(Value::Cobar(words));
            };
            let words = cobar::add_support(model, &words, &correction)?;
            model.allocate(words.iter().map(cobar::word_cells).sum())?;
            Ok(Value::Cobar(words.into()))
        }
        Value::MapCobar(terms) => {
            if terms.is_empty() {
                return Ok(Value::MapCobar(terms));
            }
            let Kind::MapCobar(source, target, s, t) = *model.kind(ty) else {
                unreachable!()
            };
            let mut complex = MappingComplex::new(model, source, target, t)?;
            let Some(correction) = correction(model, &mut complex, r, s, f, &terms)? else {
                return Ok(Value::MapCobar(terms));
            };
            let terms = map_ext::add_terms(model, &terms, &correction)?;
            Ok(Value::MapCobar(map_cochain::persist(model, terms)?))
        }
        _ => unreachable!(),
    }
}

fn correction<C: FilteredComplex>(
    model: &mut Model,
    complex: &mut C,
    r: usize,
    s: usize,
    f: usize,
    input: &[C::Term],
) -> Result<Option<Vec<C::Term>>> {
    let cutoff = f.saturating_add(r);
    let mut below_cutoff = Vec::new();
    for term in input {
        model.spend(C::cells(term))?;
        let filtration = complex.filtration(term)?;
        if filtration < f {
            return Err(Error::plain(format!(
                "mayLiftTo input lies below requested filtration {f}"
            )));
        }
        // Terms at or above the cutoff cannot affect the projected differential.
        if filtration < cutoff {
            below_cutoff.push(term.clone());
        }
    }
    if below_cutoff.is_empty() {
        return Ok(None);
    }
    let image = complex.differential_support(model, &below_cutoff)?;
    let mut rows = Rows::new("mayLiftTo");
    let (required, lowest) = projected_rows(model, complex, &mut rows, cutoff, image)?;
    if required.is_empty() {
        return Ok(None);
    }
    let obstruction = || {
        Error::plain(format!(
            "mayLiftTo has a nonzero obstruction below filtration {cutoff} with filtration {f} terms fixed"
        ))
    };
    if lowest <= f {
        return Err(obstruction());
    }

    let mut basis = Basis::default();
    let mut candidates = Vec::new();
    let mut candidate_cells = 0;
    // The nonzero projected differential forces cutoff>f+1 at this point.
    // A correction at or above cutoff has no projected image.
    let band = f + 1..cutoff;
    for term in complex.basis(model, s, band.start)? {
        let cells = C::cells(&term);
        model.spend(cells)?;
        if !band.contains(&complex.filtration(&term)?) {
            continue;
        }
        let differential = complex.differential(model, &term)?;
        let (column, _) = projected_rows(model, complex, &mut rows, cutoff, differential)?;
        if !column.is_empty()
            && basis
                .insert_with_witness(model, column, vec![candidates.len()])?
                .is_none()
        {
            candidate_cells += cells;
            may::check_storage(candidate_cells + basis.cells(), "mayLiftTo")?;
            candidates.push(term);
        }
    }
    let (remainder, witness) = basis.reduce_with_witness(model, required, Vec::new())?;
    if !remainder.is_empty() {
        return Err(obstruction());
    }
    model.spend(witness.len())?;
    Ok(Some(
        witness
            .into_iter()
            .map(|index| candidates[index].clone())
            .collect(),
    ))
}

fn projected_rows<C: FilteredComplex>(
    model: &mut Model,
    complex: &C,
    rows: &mut Rows<C::Term>,
    cutoff: usize,
    terms: Vec<C::Term>,
) -> Result<(Vec<usize>, usize)> {
    let mut support = Vec::new();
    let mut lowest = usize::MAX;
    for term in terms {
        let cells = C::cells(&term);
        model.spend(cells)?;
        let filtration = complex.filtration(&term)?;
        if filtration < cutoff {
            lowest = lowest.min(filtration);
            support.push(rows.coordinate(model, cells, term)?);
        }
    }
    may::sort_support(model, &mut support)?;
    Ok((support, lowest))
}
