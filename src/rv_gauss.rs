use crate::error::{Error, Result};
use crate::model::Model;
use crate::q_koszul::{self, Indices};
use crate::rv_ideal::{self, Ideal};
use crate::rv_poly::{self, Monomial, Polynomial};
use crate::{rv_rees, rv_scheme};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Graph,
    Dual,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Graph => "rvgauss",
            Self::Dual => "rvdual",
        }
    }
}

pub(crate) fn build(
    model: &mut Model,
    operation: Operation,
    selected: Indices,
    degree: usize,
    hypersurface: &Polynomial,
    dual_variables: Vec<usize>,
) -> Result<Ideal> {
    let name = operation.name();
    if degree == 0 || hypersurface.is_empty() {
        return Err(Error::plain(format!(
            "{name} requires a nonzero homogeneous polynomial of positive degree"
        )));
    }
    let width = selected.len();
    if dual_variables.len() != width {
        return Err(Error::plain(format!(
            "{name} requires one dual coordinate per original variable"
        )));
    }
    let dual = q_koszul::indices(model, dual_variables)?;
    if selected
        .iter()
        .any(|index| dual.binary_search(index).is_ok())
    {
        return Err(Error::plain(format!(
            "{name} requires fresh dual coordinates"
        )));
    }
    let full_width = width + dual.len();

    let mut gradient = Vec::with_capacity(width);
    let mut orders = vec![0; width];
    let mut cells = 0;
    for i in 0..width {
        orders[i] = 1;
        let partial = rv_poly::hasse(model, hypersurface, &orders)?;
        orders[i] = 0;
        cells += rv_poly::cells(&partial);
        rv_ideal::check_storage(cells)?;
        gradient.push(partial);
    }

    let graph = if gradient.iter().all(|partial| partial.is_empty()) {
        // There is no point where the hypersurface has a defined Gauss map.
        vec![rv_poly::monomial(vec![0; full_width])]
    } else {
        let raw = rv_rees::graph(model, std::iter::once(hypersurface), width, &gradient)?;
        let over_base = rv_scheme::saturate_coordinates(model, &raw, full_width, 0..width)?;
        rv_scheme::saturate_coordinates(model, &over_base, full_width, width..full_width)?
    };

    match operation {
        Operation::Graph => graph_ideal(model, &selected, &dual, graph),
        Operation::Dual => {
            // Contract to the dual coordinates. Saturation then takes its
            // projective scheme-theoretic image in the dual projective space.
            let image = rv_ideal::eliminate_prefix(model, graph, width)?;
            let image = rv_scheme::saturate_coordinates(model, &image, width, 0..width)?;
            Ideal::derived(model, dual, image)
        }
    }
}

fn graph_ideal(
    model: &mut Model,
    base: &[usize],
    dual: &[usize],
    mut graph: Vec<Polynomial>,
) -> Result<Ideal> {
    // Elimination keeps base coordinates before dual coordinates. Restore the
    // ordinary RVPoly coordinate order even when fresh labels interleave.
    let original: Vec<_> = base.iter().chain(dual).copied().collect();
    let mut permuted: Vec<_> = (0..original.len()).collect();
    permuted.sort_unstable_by_key(|&i| original[i]);
    let selected = q_koszul::indices(model, permuted.iter().map(|&i| original[i]).collect())?;
    if !permuted.iter().copied().eq(0..original.len()) {
        for polynomial in &mut graph {
            let cost = rv_poly::cells(polynomial);
            model.spend(cost * (1 + polynomial.len().max(1).ilog2() as usize))?;
            let mut terms: Vec<Monomial> = polynomial
                .iter()
                .map(|term| permuted.iter().map(|&i| term[i]).collect())
                .collect();
            terms.sort_unstable();
            *polynomial = terms.into();
        }
    }
    Ideal::derived(model, selected, graph)
}
