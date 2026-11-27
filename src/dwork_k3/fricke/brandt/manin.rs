use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::lee_linear::Rational as Q;
use crate::model::{Model, Value};
use crate::rv_lens::gcd;
use crate::syntax::{Label, ParsedExpr};

use super::eigen::Fraction;
use geometry::Geometry;
use linear::{self as lin, Frame, Matrix, Vector, ONE, ZERO};

mod comparison;
pub(crate) mod frobenius;
mod geometry;
mod linear;

pub(crate) use comparison::Comparison;

pub(crate) enum Command {
    Define(Label, usize),
    Data(Label),
    Spectrum(Label, usize, usize),
    Compare(Label, Label, Label),
    ComparisonData(Label),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<Option<String>> {
        match self {
            Self::Define(name, p) => {
                let space = Space::new(model, p).map_err(|error| error.at(name.span))?;
                model.define_manin_space(name, space)?;
                Ok(None)
            }
            Self::Data(name) => Ok(Some(model.manin_space(&name)?.describe())),
            Self::Spectrum(name, n, cut) => {
                let space = model.manin_space(&name)?;
                space.spectrum(model, n, cut).map(Some)
            }
            Self::Compare(name, manin, brandt) => {
                let space = model.manin_space(&manin)?;
                let basis = model.fricke_ss(&brandt)?;
                let comparison = Comparison::new(model, manin, brandt, space, &basis)
                    .map_err(|error| error.at(name.span))?;
                model.define_manin_comparison(name, comparison)?;
                Ok(None)
            }
            Self::ComparisonData(name) => Ok(Some(model.manin_comparison(&name)?.describe())),
        }
    }
}

pub(crate) enum Expression {
    Symbol(Label, i128, i128),
    Path(Label, [[i128; 2]; 2]),
    Basis(Label, usize, usize),
    Action(Operator, Box<ParsedExpr>),
    Transfer(Label, bool, Box<ParsedExpr>),
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) enum Operator {
    J,
    Atkin(usize),
    Hecke(usize),
    Scale(Q),
}

pub(crate) struct Space {
    geometry: Geometry,
    cusp: Frame,
    plus: Frame,
    joint: Frame,
    operators: RefCell<BTreeMap<Operator, Rc<Matrix>>>,
}

impl Space {
    fn new(model: &mut Model, p: usize) -> Result<Self> {
        let geometry = Geometry::new(model, p)?;
        let d = geometry.dimension();
        let j = geometry.operator(model, Operator::J)?;
        let w2 = geometry.operator(model, Operator::Atkin(2))?;
        verify_involution(model, &j)?;
        verify_involution(model, &w2)?;
        verify_commutation(model, &j, &w2)?;
        let cusp = Frame::nullspace(model, d, geometry.boundary.clone())?;
        if d - cusp.vectors.len() != 3 {
            return Err(certificate(
                "the four cusp classes have boundary rank three",
            ));
        }
        let mut equations = geometry.boundary.clone();
        append_fixed(&mut equations, &j)?;
        let plus = Frame::nullspace(model, d, equations.clone())?;
        append_fixed(&mut equations, &w2)?;
        let joint = Frame::nullspace(model, d, equations)?;
        model.allocate(
            32 + d * (8 * d + 4 * (cusp.vectors.len() + plus.vectors.len() + joint.vectors.len())),
        )?;
        let operators = RefCell::new(BTreeMap::from([
            (Operator::J, Rc::new(j)),
            (Operator::Atkin(2), Rc::new(w2)),
        ]));
        Ok(Self {
            geometry,
            cusp,
            plus,
            joint,
            operators,
        })
    }

    pub(crate) fn dimension(&self) -> usize {
        self.geometry.dimension()
    }

    pub(crate) fn symbol(&self, c: i128, d: i128) -> Result<Vector> {
        self.geometry.symbol(c, d)
    }

    pub(crate) fn path(&self, model: &mut Model, [start, end]: [[i128; 2]; 2]) -> Result<Vector> {
        if [start, end]
            .iter()
            .flatten()
            .any(|x| x.unsigned_abs() > geometry::MAX_CUSP_ENTRY)
        {
            return Err(Error::plain(
                "Manin cusp coordinates must have absolute value at most 10^12",
            ));
        }
        self.geometry.path(model, start, end)
    }

    pub(crate) fn basis(&self, cut: usize, index: usize) -> Result<Vector> {
        let value = if let Some(frame) = self.frame(cut)? {
            frame.vectors.get(index).cloned()
        } else if index < self.dimension() {
            let mut vector = vec![ZERO; self.dimension()];
            vector[index] = ONE;
            Some(vector)
        } else {
            None
        };
        value.ok_or_else(|| Error::plain("Manin basis index is out of range"))
    }

    fn frame(&self, cut: usize) -> Result<Option<&Frame>> {
        match cut {
            0 => Ok(None),
            1 => Ok(Some(&self.cusp)),
            2 => Ok(Some(&self.plus)),
            3 => Ok(Some(&self.joint)),
            _ => Err(Error::plain("Manin sector must be 0 (relative), 1 (cuspidal), 2 (cuspidal J+), or 3 (cuspidal J+ W2+)")),
        }
    }

    pub(crate) fn check_operation(&self, operation: Operator) -> Result<()> {
        if !matches!(operation, Operator::Scale(_)) {
            self.geometry.transforms(operation)?;
        }
        Ok(())
    }

    fn operator(&self, model: &mut Model, operation: Operator) -> Result<Rc<Matrix>> {
        if let Some(matrix) = self.operators.borrow().get(&operation) {
            return Ok(Rc::clone(matrix));
        }
        let matrix = self.geometry.operator(model, operation)?;
        if matches!(operation, Operator::J | Operator::Atkin(_)) {
            verify_involution(model, &matrix)?;
        }
        for other in self.operators.borrow().values() {
            verify_commutation(model, &matrix, other)?;
        }
        model.allocate(4 * matrix.len().pow(2) + 4)?;
        let matrix = Rc::new(matrix);
        self.operators
            .borrow_mut()
            .insert(operation, Rc::clone(&matrix));
        Ok(matrix)
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        operation: Operator,
        value: &Fraction,
    ) -> Result<Value> {
        let source = rationals(value)?;
        let result = if let Operator::Scale(k) = operation {
            source
                .iter()
                .map(|&x| x.multiply(k))
                .collect::<Result<_>>()?
        } else {
            let matrix = self.operator(model, operation)?;
            lin::apply(model, &matrix, &source)?
        };
        persist(model, &result)
    }

    pub(crate) fn boundary(&self, model: &mut Model, value: &Fraction) -> Result<String> {
        let result = lin::apply(model, &self.geometry.boundary, &rationals(value)?)?;
        Ok(format!(
            "cusp denominators [1, 2, {}, {}]; boundary = {}",
            self.geometry.prime,
            self.geometry.level,
            lin::describe(&result)
        ))
    }

    fn describe(&self) -> String {
        format!("level {}: {} Manin symbols; relative dimension {}; cuspidal dimension {}; cuspidal J+ dimension {}; cuspidal J+ W2+ dimension {}",
            self.geometry.level,self.geometry.count(),self.dimension(),self.cusp.vectors.len(),self.plus.vectors.len(),self.joint.vectors.len())
    }

    fn spectrum(&self, model: &mut Model, n: usize, cut: usize) -> Result<String> {
        let matrix = self.operator(model, Operator::Hecke(n))?;
        let characteristic = if let Some(frame) = self.frame(cut)? {
            let restriction = frame.restrict(model, &matrix)?;
            lin::characteristic(model, &restriction)?
        } else {
            lin::characteristic(model, &matrix)?
        };
        Ok(format!(
            "T{n} sector {cut}: characteristic (descending) = {}",
            lin::describe(&characteristic)
        ))
    }
}

pub(crate) fn persist(model: &mut Model, vector: &[Q]) -> Result<Value> {
    let value = fraction(vector)?;
    model.allocate(6 + 2 * vector.len())?;
    Ok(Value::FrickeManin(Rc::new(value)))
}

fn fraction(vector: &[Q]) -> Result<Fraction> {
    let denominator = vector.iter().try_fold(1i128, |d, v| {
        let common = gcd(d as u128, v.denominator as u128) as i128;
        (d / common)
            .checked_mul(v.denominator)
            .ok_or_else(arithmetic)
    })?;
    let numerator = vector
        .iter()
        .map(|v| {
            v.numerator
                .checked_mul(denominator / v.denominator)
                .ok_or_else(arithmetic)
        })
        .collect::<Result<_>>()?;
    Ok(Fraction::new(numerator, denominator))
}

fn rationals(fraction: &Fraction) -> Result<Vector> {
    fraction
        .numerator
        .iter()
        .map(|&n| Q::new(n, fraction.denominator))
        .collect()
}

pub(crate) fn add(model: &mut Model, left: &Fraction, right: &Fraction) -> Result<Value> {
    let left = rationals(left)?;
    let right = rationals(right)?;
    let result = left
        .iter()
        .zip(&right)
        .map(|(&a, &b)| a.add(b))
        .collect::<Result<Vector>>()?;
    persist(model, &result)
}

fn append_fixed(equations: &mut Matrix, involution: &Matrix) -> Result<()> {
    for (i, row) in involution.iter().enumerate() {
        let mut row = row.clone();
        row[i] = lin::sub(row[i], ONE)?;
        equations.push(row);
    }
    Ok(())
}

fn verify_involution(model: &mut Model, matrix: &Matrix) -> Result<()> {
    if lin::multiply(model, matrix, matrix)? == lin::identity(matrix.len()) {
        Ok(())
    } else {
        Err(certificate("the cusp involution squares to identity"))
    }
}

fn verify_commutation(model: &mut Model, a: &Matrix, b: &Matrix) -> Result<()> {
    if lin::multiply(model, a, b)? == lin::multiply(model, b, a)? {
        Ok(())
    } else {
        Err(certificate("the exact Hecke and cusp involutions commute"))
    }
}

fn arithmetic() -> Error {
    Error::plain("exact Manin symbol arithmetic exceeds the signed 128-bit range")
}

fn certificate(identity: &str) -> Error {
    Error::plain(format!("Manin symbol certificate failed: {identity}"))
}
