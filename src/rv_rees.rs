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

    let width = selected.len();
    let mut graph = Vec::new();
    let mut center = Vec::new();
    let mut cells = 0;
    for polynomial in quotient.polynomials() {
        charge(model, &mut cells, width + 1, polynomial.len())?;
        graph.push(
            polynomial
                .iter()
                .map(|term| embed(term, Some(0), extras))
                .collect(),
        );
    }
    for (index, (_, generator)) in generators.into_iter().enumerate() {
        // Replacing h_j by its residue leaves the graph ideal unchanged.
        let normal = quotient.normal(model, generator)?;
        charge(model, &mut cells, width + 1, 1 + normal.len())?;
        let mut variable = vec![0; width + 1];
        variable[1 + quotient.selected.len() + index] = 1;
        // T comes first in lex order. The U_j term precedes the T*h_j terms.
        graph.push(
            once(variable.into())
                .chain(normal.iter().map(|term| embed(term, Some(1), extras)))
                .collect(),
        );
        if matches!(operation, Operation::NormalCone) && !normal.is_empty() {
            charge(model, &mut cells, width, normal.len())?;
            center.push(
                normal
                    .iter()
                    .map(|term| embed(term, None, extras))
                    .collect(),
            );
        }
    }

    // Contract the graph of U_j -> T*h_j from the polynomial ring with T.
    // For the normal cone, gr_H(R) is Rees_R(H) modulo the degree-zero image of H.
    let mut relations = rv_ideal::eliminate(model, graph)?;
    rv_ideal::check_storage(
        relations
            .iter()
            .chain(&center)
            .map(|p| rv_poly::cells(p))
            .sum(),
    )?;
    relations.extend(center);
    Ideal::derived(model, selected, relations)
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
