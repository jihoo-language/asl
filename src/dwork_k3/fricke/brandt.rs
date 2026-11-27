use std::cell::OnceCell;
use std::rc::Rc;

use crate::ec_poly::Scalar;
use crate::error::{Error, Result, Span};
use crate::model::{Model, Value};
use crate::rv_unram::Ring;
use crate::syntax::{Label, ParsedExpr};

mod geometry;
mod matrix;

pub(crate) enum Command {
    Define(Label, Label),
    Describe(Label),
    Matrix(Label, usize),
    Relations(Label),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<Option<String>> {
        match self {
            Self::Define(name, field) => {
                let ring = model.rv_unram(&field)?;
                let basis = Basis::geometry(model, field.text.into(), ring)
                    .map_err(|error| error.at(name.span))?;
                model.define_fricke_ss(name, basis)?;
                Ok(None)
            }
            Self::Describe(name) => Ok(Some(model.fricke_ss(&name)?.describe())),
            Self::Matrix(name, level) => {
                let basis = model.fricke_ss(&name)?;
                let matrix = basis.operator(model, level)?;
                matrix.spectra(model, &basis.frobenius, level).map(Some)
            }
            Self::Relations(name) => {
                let basis = model.fricke_ss(&name)?;
                let three = basis.operator(model, 3)?;
                let nine = basis.operator(model, 9)?;
                let five = if basis.p == 5 {
                    None
                } else {
                    Some(basis.operator(model, 5)?)
                };
                matrix::relations(model, &basis.frobenius, three, five, nine).map(Some)
            }
        }
    }
}

pub(crate) enum Expression {
    Divisor(Label, Vec<i128>),
    Basis(Label, usize),
    Scalar(Label, ScalarEntry),
    Action(Operation, ParsedExpr),
}

pub(crate) enum ScalarEntry {
    Vertex(usize),
    Endpoint(usize, usize),
    Polynomial { fricke: bool, index: usize },
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Hecke(usize),
    Frobenius,
    Negate,
}

pub(crate) struct Basis {
    field_name: Rc<str>,
    ring: Rc<Ring>,
    p: u64,
    js: Vec<Scalar>,
    j_polynomial: Vec<u64>,
    h_polynomial: Vec<u64>,
    vertices: Vec<Scalar>,
    endpoints: Vec<[Scalar; 2]>,
    frobenius: Vec<usize>,
    operators: [OnceCell<matrix::Matrix>; 3],
}

impl Basis {
    pub(crate) fn field_label(&self, span: Span) -> Label {
        Label {
            text: self.field_name.to_string(),
            span,
        }
    }

    pub(crate) fn rank(&self) -> usize {
        self.vertices.len()
    }

    pub(crate) fn div_basis(&self, index: usize) -> Result<Vec<i128>> {
        let mut result = vec![0; self.rank()];
        *result.get_mut(index).ok_or_else(vertex_range)? = 1;
        Ok(result)
    }

    pub(crate) fn scalar(&self, entry: ScalarEntry) -> Result<Scalar> {
        match entry {
            ScalarEntry::Vertex(i) => self.vertices.get(i).cloned().ok_or_else(vertex_range),
            ScalarEntry::Endpoint(i, j) => self
                .endpoints
                .get(i)
                .ok_or_else(vertex_range)?
                .get(j)
                .cloned()
                .ok_or_else(|| Error::plain("Fricke endpoint index must be 0 or 1")),
            ScalarEntry::Polynomial { fricke, index } => {
                let coefficients = if fricke {
                    &self.h_polynomial
                } else {
                    &self.j_polynomial
                };
                Ok(self
                    .ring
                    .scalar(coefficients.get(index).copied().unwrap_or(0)))
            }
        }
    }

    pub(crate) fn check_operation(&self, operation: Operation) -> Result<()> {
        if let Operation::Hecke(level) = operation {
            matrix::check_level(self.p, level)?;
        }
        Ok(())
    }

    fn operator(&self, model: &mut Model, level: usize) -> Result<&matrix::Matrix> {
        matrix::check_level(self.p, level)?;
        let cell = &self.operators[match level {
            3 => 0,
            5 => 1,
            9 => 2,
            _ => return Err(certificate("supported Brandt level")),
        }];
        if let Some(matrix) = cell.get() {
            return Ok(matrix);
        }
        let matrix = matrix::Matrix::hecke(model, self, &self.ring, level)?;
        cell.set(matrix)
            .map_err(|_| certificate("unique cached operator"))?;
        cell.get()
            .ok_or_else(|| certificate("cached operator present"))
    }

    fn describe(&self) -> String {
        format!(
            "geometric supersingular j values (p = {}) = {:?}\n             j polynomial (ascending) = {:?}\n             Fricke vertices = {:?}; unordered endpoints = {:?}\n             Fricke polynomial (ascending) = {:?}\n             Frobenius permutation = {:?}",
            self.p, self.js, self.j_polynomial, self.vertices, self.endpoints,
            self.h_polynomial, self.frobenius,
        )
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        operation: Operation,
        value: &[i128],
    ) -> Result<Value> {
        let result = match operation {
            Operation::Hecke(level) => self.operator(model, level)?.push(model, value)?,
            Operation::Frobenius => self.frobenius.iter().map(|&j| value[j]).collect(),
            Operation::Negate => value
                .iter()
                .map(|c| c.checked_neg().ok_or_else(arithmetic))
                .collect::<Result<_>>()?,
        };
        persist(model, result)
    }
}

pub(crate) fn persist(model: &mut Model, value: Vec<i128>) -> Result<Value> {
    model.spend(value.len())?;
    model.allocate(2 + 2 * value.len())?;
    Ok(Value::FrickeSSDiv(value.into()))
}

pub(crate) fn add(model: &mut Model, left: &[i128], right: &[i128]) -> Result<Value> {
    let sum = left
        .iter()
        .zip(right)
        .map(|(a, b)| a.checked_add(*b).ok_or_else(arithmetic))
        .collect::<Result<_>>()?;
    persist(model, sum)
}

fn vertex_range() -> Error {
    Error::plain("supersingular Fricke vertex index is out of range")
}

fn arithmetic() -> Error {
    Error::plain("supersingular Fricke integral arithmetic exceeds the signed 128-bit range")
}

fn certificate(identity: &str) -> Error {
    Error::plain(format!(
        "supersingular Fricke certificate failed: {identity}"
    ))
}
