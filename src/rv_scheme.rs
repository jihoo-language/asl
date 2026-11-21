use std::iter::once;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_ideal::{self, Ideal};
use crate::rv_poly::{self, Monomial, Polynomial};

pub(crate) fn intersect(model: &mut Model, left: &Ideal, right: &Ideal) -> Result<Ideal> {
    if left.selected != right.selected {
        return Err(Error::plain(
            "rvintersect requires ideals with the same primitive selection",
        ));
    }
    let left_generators: Vec<_> = left.polynomials().cloned().collect();
    let right_generators: Vec<_> = right.polynomials().cloned().collect();
    let result = intersection(model, &left_generators, &right_generators)?;
    Ideal::derived(model, left.selected.clone(), result)
}

pub(crate) fn colon(model: &mut Model, ideal: &Ideal, polynomial: Polynomial) -> Result<Ideal> {
    let width = ideal.selected.len();
    let result = if polynomial.is_empty() {
        unit(width)
    } else {
        let generators: Vec<_> = ideal.polynomials().cloned().collect();
        let multiples = intersection(model, &generators, std::slice::from_ref(&polynomial))?;
        multiples
            .into_iter()
            .map(|p| rv_poly::exact_divide(model, &p, &polynomial))
            .collect::<Result<_>>()?
    };
    Ideal::derived(model, ideal.selected.clone(), result)
}

pub(crate) fn saturate(model: &mut Model, ideal: &Ideal, polynomial: &Polynomial) -> Result<Ideal> {
    let generators: Vec<_> = ideal.polynomials().cloned().collect();
    let result = saturation(model, &generators, polynomial, ideal.selected.len())?;
    Ideal::derived(model, ideal.selected.clone(), result)
}

pub(crate) fn projective(model: &mut Model, ideal: &Ideal) -> Result<Ideal> {
    let width = ideal.selected.len();
    let generators: Vec<_> = ideal.polynomials().cloned().collect();
    let mut result = unit(width);
    for variable in 0..width {
        model.spend(1 + width)?;
        let mut powers = vec![0; width];
        powers[variable] = 1;
        let component = saturation(model, &generators, &rv_poly::monomial(powers), width)?;
        result = intersection(model, &result, &component)?;
        if result.is_empty() {
            break;
        }
    }
    Ideal::derived(model, ideal.selected.clone(), result)
}

fn unit(width: usize) -> Vec<Polynomial> {
    vec![rv_poly::monomial(vec![0; width])]
}

fn is_one(polynomial: &Polynomial) -> bool {
    polynomial.len() == 1 && polynomial[0].iter().all(|&power| power == 0)
}

fn intersection(
    model: &mut Model,
    left: &[Polynomial],
    right: &[Polynomial],
) -> Result<Vec<Polynomial>> {
    if left.is_empty() || right.is_empty() {
        return Ok(Vec::new());
    }
    if left.iter().any(is_one) {
        return Ok(right.to_vec());
    }
    if right.iter().any(is_one) {
        return Ok(left.to_vec());
    }

    let mut extended = Auxiliary::default();
    // Contract u*I + (u+1)*J to the original ring. Evaluating at u=1
    // and u=0 proves containment in both ideals; every common element
    // g belongs because u*g + (u+1)*g = g.
    for polynomial in left {
        extended.push(model, polynomial, &[1])?;
    }
    for polynomial in right {
        extended.push(model, polynomial, &[0, 1])?;
    }
    rv_ideal::eliminate(model, extended.generators)
}

fn saturation(
    model: &mut Model,
    generators: &[Polynomial],
    polynomial: &Polynomial,
    width: usize,
) -> Result<Vec<Polynomial>> {
    if polynomial.is_empty() {
        return Ok(unit(width));
    }
    if is_one(polynomial) || generators.is_empty() || generators.iter().any(is_one) {
        return Ok(generators.to_vec());
    }
    if polynomial.len() == 1 && generators.iter().all(|p| p.len() == 1) {
        // Inverting a monomial inverts every variable dividing it. For a
        // monomial ideal, contraction removes these factors from its generators.
        model.spend(1 + generators.len() * (1 + width))?;
        let mut result = Vec::with_capacity(generators.len());
        for generator in generators {
            let powers: Vec<_> = generator[0]
                .iter()
                .zip(polynomial[0].iter())
                .map(|(&power, &inverted)| if inverted == 0 { power } else { 0 })
                .collect();
            if powers.iter().all(|&power| power == 0) {
                return Ok(unit(width));
            }
            result.push(rv_poly::monomial(powers));
        }
        return Ok(result);
    }

    let mut extended = Auxiliary::default();
    for generator in generators {
        extended.push(model, generator, &[0])?;
    }
    // Adjoining u = 1/f gives the localization. Its contraction is J:f^∞.
    extended.charge(model, 2 + (1 + polynomial.len()) * (2 + width))?;
    let relation: Polynomial = once(vec![0; width + 1].into())
        .chain(lift(polynomial, 1))
        .collect();
    extended.generators.push(relation);
    rv_ideal::eliminate(model, extended.generators)
}

#[derive(Default)]
struct Auxiliary {
    generators: Vec<Polynomial>,
    cells: usize,
}

impl Auxiliary {
    fn charge(&mut self, model: &mut Model, cells: usize) -> Result<()> {
        self.cells += cells;
        rv_ideal::check_storage(self.cells)?;
        model.spend(cells)
    }

    fn push(
        &mut self,
        model: &mut Model,
        polynomial: &Polynomial,
        exponents: &[usize],
    ) -> Result<()> {
        if let Some(term) = polynomial.first() {
            self.charge(
                model,
                2 + polynomial.len() * exponents.len() * (2 + term.len()),
            )?;
            // Exponents are increasing, so these blocks stay sorted in pure lex
            // with the new variable first and the original variables following.
            let extended = exponents
                .iter()
                .flat_map(|&power| lift(polynomial, power))
                .collect();
            self.generators.push(extended);
        }
        Ok(())
    }
}

fn lift(polynomial: &Polynomial, power: usize) -> impl Iterator<Item = Monomial> + '_ {
    polynomial
        .iter()
        .map(move |term| once(power).chain(term.iter().copied()).collect())
}
