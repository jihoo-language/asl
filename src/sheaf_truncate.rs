use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::sheaf::{self, Cochains, Matrix};
use crate::sheaf_hyper::{Complex, Map, NamedSheaf};
use crate::sheaf_image::Site;
use crate::sheaf_leray::{component_name, Input, InputForm, Leray};
use crate::sheaf_stalk::{self, StalkQuotient};
use crate::syntax::Label;

pub(crate) fn truncate(
    model: &mut Model,
    name: &Label,
    complex: Rc<Complex>,
    cutoff: usize,
) -> Result<(Rc<Complex>, Option<NamedSheaf>)> {
    let Some(outgoing) = complex.maps().get(cutoff) else {
        return Ok((complex, None));
    };
    let original = complex.component(cutoff)?;
    let terminal_name = component_name(model, name, cutoff)?;
    let prime = complex.prime();
    let kernels = original
        .sheaf
        .dims()
        .iter()
        .enumerate()
        .map(|(p, &rank)| StalkQuotient::build(model, prime, rank, None, Some(outgoing.matrix(p))))
        .collect::<Result<Vec<_>>>()?;

    let sheaf = sheaf_stalk::induced_sheaf(model, &original.sheaf, &kernels)?;
    let terminal = NamedSheaf {
        name: terminal_name,
        sheaf: Rc::new(sheaf),
    };

    let mut components = complex.components()[..cutoff].to_vec();
    let mut maps = Vec::new();
    for (degree, map) in complex
        .maps()
        .iter()
        .take(cutoff.saturating_sub(1))
        .enumerate()
    {
        maps.push((differential_label(name, degree), map.clone()));
    }
    if cutoff > 0 {
        // The previous differential already lands in the kernel by d²=0.
        let before = &complex.maps()[cutoff - 1];
        let dims = before.source.sheaf.dims();
        let mut matrices = Vec::new();
        let mut cells = 0;
        for (point, kernel) in kernels.iter().enumerate() {
            let original = before.matrix(point);
            let matrix = if original.is_zero() {
                Matrix::default()
            } else {
                Matrix::from_action(model, dims[point], |model, index| {
                    kernel.coordinates(model, prime, original.column(index).to_vec())
                })?
            };
            cells += matrix.cells();
            sheaf::check_cells(cells)?;
            matrices.push(Rc::new(matrix));
        }
        let incoming =
            Map::from_matrices(model, before.source.clone(), terminal.clone(), matrices)?;
        maps.push((differential_label(name, cutoff - 1), Rc::new(incoming)));
    }
    components.push(terminal.clone());
    let result = Complex::build(model, components, maps)?;
    Ok((Rc::new(result), Some(terminal)))
}

pub(crate) fn deligne(
    model: &mut Model,
    name: &Label,
    source: Input,
    site: Site,
    cutoff: usize,
) -> Result<Rc<Complex>> {
    let form = source.form();
    if source.base() != site.source {
        let message = match form {
            InputForm::Sheaf => "shdeligne requires a sheaf on the open embedding's named source",
            InputForm::Complex => {
                "shdelignecomplex requires a complex on the open embedding's named source"
            }
        };
        return Err(Error::plain(message));
    }
    site.require_open_embedding(model, form.deligne_name())?;
    // The top kernel needs one additional Roos degree.
    let (_, _, direct) = Leray::build_through(model, name, source, site, cutoff.saturating_add(1))?;
    let (result, _) = truncate(model, name, Rc::new(direct), cutoff)?;
    Ok(result)
}

fn differential_label(name: &Label, degree: usize) -> Label {
    Label {
        span: name.span,
        text: format!("{} differential at {}", name.text, degree),
    }
}
