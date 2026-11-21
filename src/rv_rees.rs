use std::iter::{once, repeat_n};

use crate::error::{Error, Result};
use crate::model::Model;
use crate::q_koszul;
use crate::rv_ideal::{self, Ideal};
use crate::rv_poly::{self, Monomial, Polynomial};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Rees,
    NormalCone,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Rees => "rvrees",
            Self::NormalCone => "rvnormalcone",
        }
    }
}

pub(crate) fn build(
    model: &mut Model,
    operation: Operation,
    quotient: &Ideal,
    variables: Vec<usize>,
    generators: Vec<(usize, Polynomial)>,
) -> Result<Ideal> {
    let name = operation.name();
    if variables.len() != generators.len() {
        return Err(Error::plain(format!(
            "{name} requires one new variable per center generator"
        )));
    }
    if let Some((degree, _)) = generators.first() {
        if *degree == 0 || generators.iter().any(|(other, _)| other != degree) {
            return Err(Error::plain(format!(
                "{name} center generators must have a common positive degree"
            )));
        }
    }
    let extras = variables.len();
    let mut selected = quotient.selected.to_vec();
    selected.extend(variables);
    let selected = q_koszul::indices(model, selected)?;
    if extras == 0 {
        return Ideal::derived(model, selected, quotient.polynomials().cloned().collect());
    }

    let mut normals = Vec::with_capacity(extras);
    let mut cells = 0;
    for (_, generator) in generators {
        // Replacing h_j by its residue leaves the graph ideal unchanged.
        let normal = quotient.normal(model, generator)?;
        cells += rv_poly::cells(&normal);
        rv_ideal::check_storage(cells)?;
        normals.push(normal);
    }
    let mut relations = graph(
        model,
        quotient.polynomials(),
        quotient.selected.len(),
        &normals,
    )?;
    // gr_H(R) is Rees_R(H) modulo the degree-zero image of H.
    if matches!(operation, Operation::NormalCone) {
        let mut cells = relations.iter().map(|p| rv_poly::cells(p)).sum();
        for normal in normals.iter().filter(|p| !p.is_empty()) {
            charge(model, &mut cells, selected.len(), normal.len())?;
            relations.push(
                normal
                    .iter()
                    .map(|term| embed(term, None, extras))
                    .collect(),
            );
        }
    }
    Ideal::derived(model, selected, relations)
}

// Contract the graph of U_j -> T*h_j from the polynomial ring with T.
pub(crate) fn graph<'a>(
    model: &mut Model,
    base: impl IntoIterator<Item = &'a Polynomial>,
    width: usize,
    images: &[Polynomial],
) -> Result<Vec<Polynomial>> {
    let extras = images.len();
    let total_width = 1 + width + extras;
    let mut equations = Vec::new();
    let mut cells = 0;
    for polynomial in base {
        charge(model, &mut cells, total_width, polynomial.len())?;
        equations.push(
            polynomial
                .iter()
                .map(|term| embed(term, Some(0), extras))
                .collect(),
        );
    }
    for (index, image) in images.iter().enumerate() {
        charge(model, &mut cells, total_width, 1 + image.len())?;
        let mut variable = vec![0; total_width];
        variable[1 + width + index] = 1;
        // T comes first in lex order. The U_j term precedes the T*h_j terms.
        equations.push(
            once(variable.into())
                .chain(image.iter().map(|term| embed(term, Some(1), extras)))
                .collect(),
        );
    }
    rv_ideal::eliminate(model, equations)
}

fn embed(term: &Monomial, parameter: Option<usize>, extras: usize) -> Monomial {
    parameter
        .into_iter()
        .chain(term.iter().copied())
        .chain(repeat_n(0, extras))
        .collect()
}

fn charge(model: &mut Model, current: &mut usize, width: usize, terms: usize) -> Result<()> {
    let cost = 1 + (1 + width) * terms;
    *current += cost;
    rv_ideal::check_storage(*current)?;
    model.spend(cost)
}
