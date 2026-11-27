use std::rc::Rc;

use crate::ec_morphism::Rational;
use crate::ec_poly::{self as poly, Field, Polynomial};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;

use super::super::selmer::ResidueSquares;
use super::{
    constant_nonsquare, field, primitive_coordinates, read_polynomials, representatives, Cover,
    Descent,
};

mod local;

pub(crate) struct AuxiliaryInput {
    pub(crate) omitted: usize,
    pub(crate) coordinates: [Vec<Vec<i128>>; 3],
}

pub(crate) struct Pairing {
    pub(crate) cover: Rc<Cover>,
    tangents: [Tangent; 3],
}

struct Tangent {
    indices: [usize; 2],
    coefficients: [Polynomial; 3],
    tangent_x: Option<Rational>,
}

impl Pairing {
    pub(crate) fn new(
        model: &mut Model,
        cover: Rc<Cover>,
        inputs: Vec<AuxiliaryInput>,
    ) -> Result<Self> {
        if !cover
            .pencil
            .selmer_space(model)?
            .contains(model, &cover.b)?
        {
            return Err(Error::plain(
                "Cassels first covering must be everywhere locally soluble",
            ));
        }
        let field = field(&cover.pencil.curve);
        let roots = cover.pencil.split_roots(model)?;
        let mut tangents: [Option<Tangent>; 3] = std::array::from_fn(|_| None);
        for input in inputs {
            let slot = tangents
                .get_mut(input.omitted)
                .filter(|slot| slot.is_none())
                .ok_or_else(witness_indices)?;
            let q = read_polynomials(model, field, input.coordinates)?;
            *slot = Some(Tangent::new(
                model,
                field,
                &cover.representatives,
                roots,
                input.omitted,
                q,
            )?);
        }
        let [Some(t0), Some(t1), Some(t2)] = tangents else {
            return Err(witness_indices());
        };
        let tangents = [t0, t1, t2];
        model.allocate(
            12 + tangents
                .iter()
                .map(|tangent| tangent.cells(field))
                .sum::<usize>(),
        )?;
        Ok(Self { cover, tangents })
    }

    pub(crate) fn tangents(&self, name: &str) -> String {
        let coefficients = self.tangents.each_ref().map(|t| &t.coefficients);
        format!("rvK3CasselsTangents[{name}] = coefficients in remaining U indices then W, omitted indices [0, 1, 2]={coefficients:?}")
    }

    pub(crate) fn pair(&self, model: &mut Model, name: &str, right: &Descent) -> Result<String> {
        let selmer = self.cover.pencil.selmer_space(model)?;
        if !selmer.contains(model, right)? {
            return Err(Error::plain(
                "Cassels second class must be everywhere locally soluble",
            ));
        }
        let field = field(&self.cover.pencil.curve);
        let squares = ResidueSquares::new(model, field)?;
        let eta = constant_nonsquare(model, field)?;
        let right = representatives(model, field, right)?;
        let mut sign = 1;
        let mut factors = Vec::new();
        // The Selmer support includes every finite root collision and both
        // squareclass supports. Outside it, smooth conics and primitive
        // projective data give even tangent valuations, so all symbols are 1.
        for parameter in selmer.finite_parameters().map(Some).chain([None]) {
            let evidence = local::evaluate(model, field, &squares, &eta, self, &right, parameter)?;
            sign *= evidence.sign();
            let label = parameter.map_or_else(|| "infinity".into(), |a| format!("t={a:?}"));
            factors.push(format!("{label}: {}", evidence.report()));
        }
        let value = if sign > 0 { "0" } else { "1/2" };
        Ok(format!(
            "rvK3CasselsPair[{name}] = sign={sign}, value={value} in Q/Z, local factors=[{}]",
            factors.join("; ")
        ))
    }
}

impl Tangent {
    fn new(
        model: &mut Model,
        field: Field<'_>,
        representatives: &[Polynomial; 3],
        roots: [&Polynomial; 3],
        omitted: usize,
        mut q: [Polynomial; 3],
    ) -> Result<Self> {
        primitive_coordinates(
            model,
            field,
            &mut q,
            "Cassels conic witness must be a nonzero projective point",
        )?;
        let indices = match omitted {
            0 => [1, 2],
            1 => [0, 2],
            2 => [0, 1],
            _ => return Err(witness_indices()),
        };
        let [j, k] = indices;
        let diagonal = [
            representatives[j].clone(),
            poly::scale(model, field, &representatives[k], &field.signed(-1))?,
            poly::sum(model, field, roots[j], roots[k], true)?,
        ];
        let mut coefficients = std::array::from_fn(|_| Vec::new());
        let mut evaluation = Vec::new();
        for ((coefficient, value), alpha) in coefficients.iter_mut().zip(&q).zip(&diagonal) {
            *coefficient = poly::multiply(model, field, alpha, value)?;
            let term = poly::multiply(model, field, coefficient, value)?;
            evaluation = poly::sum(model, field, &evaluation, &term, false)?;
        }
        if !evaluation.is_empty() {
            return Err(Error::plain(format!(
                "Cassels witness for omitted index {omitted} does not lie on its conic"
            )));
        }
        let tangent_x = if q[2].is_empty() {
            None
        } else {
            let denominator = poly::multiply(model, field, &q[2], &q[2])?;
            let root = poly::multiply(model, field, roots[j], &denominator)?;
            let term = poly::multiply(model, field, &coefficients[0], &q[0])?;
            let numerator = poly::sum(model, field, &root, &term, false)?;
            Some(Rational::new(model, field, numerator, denominator)?)
        };
        Ok(Self {
            indices,
            coefficients,
            tangent_x,
        })
    }

    fn cells(&self, field: Field<'_>) -> usize {
        let coefficients = self.coefficients.iter().map(Vec::len).sum::<usize>();
        let x = self
            .tangent_x
            .as_ref()
            .map_or(0, |x| x.numerator.len() + x.denominator.len());
        5 + (coefficients + x) * field.cells()
    }
}

fn witness_indices() -> Error {
    Error::plain("Cassels witnesses require each omitted index 0, 1 and 2 exactly once")
}
