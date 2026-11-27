use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::{Basis, Coordinate, QuotientBasis};
use crate::sheaf::{self, Cochains, Matrix, Sheaf};
use crate::sheaf_hyper::{Complex, Map, NamedSheaf};
use crate::sheaf_image::Site;
use crate::sheaf_leray::{component_name, Leray};
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
        .map(|(p, &rank)| Kernel::build(model, prime, rank, outgoing.matrix(p)))
        .collect::<Result<Vec<_>>>()?;

    let dims = kernels.iter().map(Kernel::rank).collect();
    let sheaf =
        Sheaf::from_cover_maps(model, original.sheaf.base(), prime, dims, |model, a, b| {
            let map = original.sheaf.restriction(a, b);
            let restricted = if map.is_zero() {
                Matrix::default()
            } else {
                Matrix::from_action(model, kernels[a].rank(), |model, index| {
                    let lift = kernels[a].representative(model, prime, index)?;
                    let image = map.apply(model, prime, &lift)?;
                    kernels[b].coordinates(model, prime, image)
                })?
            };
            Ok(Rc::new(restricted))
        })?;
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
    source: NamedSheaf,
    site: Site,
    cutoff: usize,
) -> Result<Rc<Complex>> {
    if source.sheaf.base() != site.source {
        return Err(Error::plain(
            "shdeligne requires a sheaf on the open embedding's named source",
        ));
    }
    site.require_open_embedding(model)?;
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

enum Kernel {
    Full(usize),
    Computed(QuotientBasis),
}

impl Kernel {
    fn build(model: &mut Model, prime: Prime, rank: usize, differential: &Matrix) -> Result<Self> {
        if differential.is_zero() {
            return Ok(Self::Full(rank));
        }
        if rank > MAX_PIVOT_CELLS {
            return Err(Error::plain(
                "smart truncation stalk exceeds the storage limit",
            ));
        }
        let basis = QuotientBasis::build(model, prime, Basis::default(), rank, |model, i| {
            let column = differential.column(i);
            model.spend(1 + column.len())?;
            Ok(column.to_vec())
        })?;
        Ok(Self::Computed(basis))
    }

    fn rank(&self) -> usize {
        match self {
            Self::Full(rank) => *rank,
            Self::Computed(basis) => basis.rank(),
        }
    }

    fn representative(
        &self,
        model: &mut Model,
        prime: Prime,
        index: usize,
    ) -> Result<Vec<Coordinate>> {
        match self {
            Self::Full(_) => Ok(vec![(index, 1)]),
            Self::Computed(basis) => basis.representative(model, prime, &[(index, 1)]),
        }
    }

    fn coordinates(
        &self,
        model: &mut Model,
        prime: Prime,
        raw: Vec<Coordinate>,
    ) -> Result<Vec<Coordinate>> {
        match self {
            Self::Full(_) => Ok(raw),
            Self::Computed(basis) => basis.project(model, prime, raw)?.ok_or_else(|| {
                Error::plain(
                    "smart truncation encountered a map that does not preserve the stalk kernel",
                )
            }),
        }
    }
}
