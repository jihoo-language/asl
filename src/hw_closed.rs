use std::cell::{OnceCell, RefCell};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing};
use crate::hyperelliptic;
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::rv_unram::{self, Ring};
use crate::syntax::{Label, ParsedExpr};

type Coefficient = Vec<u64>;
pub(crate) type Matrix = hw_algebra::Matrix<Coefficient>;

pub(crate) enum Command {
    Define(Label, Label, Vec<ParsedExpr>),
    Report(Label, Query, Option<Box<ParsedExpr>>),
}

pub(crate) enum Expression {
    Curve(Label, CurveOp),
    Unary(Unary, Box<ParsedExpr>),
    Reduce(Label, Box<ParsedExpr>),
    Compose(Box<[ParsedExpr; 2]>),
}

impl CoefficientRing for Ring {
    type Element = Coefficient;

    fn cells(&self) -> usize {
        self.degree() + 1
    }

    fn zero(&self) -> Coefficient {
        self.scalar(0)
    }

    fn one(&self) -> Coefficient {
        self.scalar(1)
    }

    fn is_zero(&self, value: &Coefficient) -> bool {
        value.iter().all(|&c| c == 0)
    }

    fn is_unit(&self, value: &Coefficient) -> bool {
        value.iter().any(|c| !c.is_multiple_of(self.prime()))
    }

    fn add(&self, model: &mut Model, a: &Coefficient, b: &Coefficient) -> Result<Coefficient> {
        Ring::add(self, model, a, b)
    }

    fn subtract(&self, model: &mut Model, a: &Coefficient, b: &Coefficient) -> Result<Coefficient> {
        let negative = self.negate(model, b)?;
        Ring::add(self, model, a, &negative)
    }

    fn multiply(&self, model: &mut Model, a: &Coefficient, b: &Coefficient) -> Result<Coefficient> {
        Ring::multiply(self, model, a, b)
    }

    fn inverse_unit(&self, model: &mut Model, value: &Coefficient) -> Result<Coefficient> {
        self.inverse(model, value)
    }
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) ring: Rc<str>,
    pub(crate) genus: usize,
}

pub(crate) struct Curve {
    pub(crate) grade: Grade,
    pub(crate) ring: Rc<Ring>,
    residue: Ring,
    polynomial: Vec<Coefficient>,
    alpha: RefCell<HashMap<usize, Rc<Matrix>>>,
    ratio: OnceCell<Rc<Matrix>>,
    product: OnceCell<Rc<Matrix>>,
    weil: OnceCell<Box<[i128]>>,
}

impl Curve {
    pub(crate) fn new(
        model: &mut Model,
        ring_name: Rc<str>,
        ring: Rc<Ring>,
        mut polynomial: Vec<Coefficient>,
    ) -> Result<Self> {
        require_odd(&ring)?;
        hw_algebra::storage(polynomial.len(), ring.cells())?;
        model.spend(polynomial.len() * ring.cells())?;
        while polynomial
            .last()
            .is_some_and(|coefficient| ring.is_zero(coefficient))
        {
            polynomial.pop();
        }
        if polynomial.len() < 2 {
            return Err(Error::plain(
                "hyperelliptic polynomial degree must be positive",
            ));
        }
        if !ring.is_unit(polynomial.last().unwrap()) {
            return Err(Error::plain(
                "hyperelliptic leading coefficient must be a unit modulo p",
            ));
        }
        let genus = (polynomial.len() - 2) / 2;
        hw_algebra::storage(genus.saturating_mul(genus), ring.cells())?;
        let residue = ring.residue_ring(model)?;
        require_squarefree(model, &residue, &polynomial)?;
        model.allocate(8 + polynomial.len() * ring.cells())?;
        Ok(Self {
            grade: Grade {
                ring: ring_name,
                genus,
            },
            ring,
            residue,
            polynomial,
            alpha: RefCell::default(),
            ratio: OnceCell::new(),
            product: OnceCell::new(),
            weil: OnceCell::new(),
        })
    }

    fn alpha(&self, model: &mut Model, level: usize) -> Result<Rc<Matrix>> {
        if let Some(matrix) = self.alpha.borrow().get(&level) {
            return Ok(matrix.clone());
        }
        let q = u32::try_from(level)
            .ok()
            .and_then(|level| self.ring.prime().checked_pow(level))
            .ok_or_else(|| Error::plain("Hasse–Witt p^level exceeds u64"))?;
        let value = hw_algebra::hasse_witt(model, self.ring.as_ref(), &self.polynomial, q)?;
        let value = keep(model, &self.ring, value)?;
        self.alpha.borrow_mut().insert(level, value.clone());
        Ok(value)
    }

    fn semilinear_ratio(&self, model: &mut Model) -> Result<Rc<Matrix>> {
        if let Some(value) = self.ratio.get() {
            return Ok(value.clone());
        }
        let first = self.alpha(model, 1)?;
        if hw_algebra::inverse(model, self.ring.as_ref(), &first)?.is_none() {
            return Err(nonordinary());
        }
        let level = self.ring.precision();
        let upper = self.alpha(model, level)?;
        let lower = self.alpha(model, level - 1)?;
        let twisted = frobenius(model, &self.ring, &lower, 1)?;
        let inverse = hw_algebra::inverse(model, self.ring.as_ref(), &twisted)?
            .ok_or_else(|| Error::plain("ordinary Hasse–Witt denominator is not invertible"))?;
        let ratio = hw_algebra::multiply(model, self.ring.as_ref(), &upper, &inverse)?;
        let ratio = keep(model, &self.ring, ratio)?;
        self.ratio.get_or_init(|| ratio.clone());
        Ok(ratio)
    }

    fn ordered_product(&self, model: &mut Model) -> Result<Rc<Matrix>> {
        if let Some(value) = self.product.get() {
            return Ok(value.clone());
        }
        let ratio = self.semilinear_ratio(model)?;
        if self.ring.degree() == 1 {
            return Ok(ratio);
        }
        let mut product = ratio.as_ref().clone();
        model.spend(product.entries().len() * self.ring.cells())?;
        for conjugate in 1..self.ring.degree() {
            let factor = frobenius(model, &self.ring, &ratio, conjugate)?;
            product = hw_algebra::multiply(model, self.ring.as_ref(), &product, &factor)?;
        }
        let product = keep(model, &self.ring, product)?;
        self.product.get_or_init(|| product.clone());
        Ok(product)
    }

    fn count_trace(&self, model: &mut Model, extension_degree: usize) -> Result<(u64, i128)> {
        if extension_degree == 0 {
            return Err(Error::plain(
                "Hasse–Witt point extension degree must be positive",
            ));
        }
        let residue_degree = self
            .ring
            .degree()
            .checked_mul(extension_degree)
            .ok_or_else(|| Error::plain("Hasse–Witt residue extension degree is out of range"))?;
        let extension = if extension_degree == 1 {
            None
        } else {
            let prime = Prime::validated(self.ring.prime() as u32);
            Some(Ring::field_of_degree(model, prime, residue_degree)?)
        };
        let field = extension.as_ref().unwrap_or(&self.residue);
        let q = field.size_of_residue_field()?;
        let coefficients = self
            .ring
            .embed_residue_coefficients(model, field, &self.polynomial)?;
        let mut trace = i128::from(field.residue_quadratic_trace(model, |model, x| {
            field.evaluate_polynomial(model, &coefficients, x)
        })?);
        if !self.polynomial.len().is_multiple_of(2) {
            // Even degree has 1+chi(lc) rational points at infinity.
            let character = field.power(model, coefficients.last().unwrap(), (q - 1) / 2)?;
            trace -= match rv_unram::require_scalar(&character)? {
                1 => 1,
                c if c == field.prime() - 1 => -1,
                _ => {
                    return Err(Error::plain(
                        "invalid quadratic character of the hyperelliptic leading coefficient",
                    ))
                }
            };
        }
        Ok((q, trace))
    }

    fn weil(&self, model: &mut Model) -> Result<&[i128]> {
        if let Some(polynomial) = self.weil.get() {
            return Ok(polynomial);
        }
        let q = self.ring.size_of_residue_field()?;
        let traces = (1..=self.grade.genus)
            .map(|r| self.count_trace(model, r).map(|(_, trace)| trace))
            .collect::<Result<Vec<_>>>()?;
        let polynomial = hyperelliptic::weil_polynomial(model, q, &traces)?;
        model.allocate(polynomial.len() * 2)?;
        Ok(self.weil.get_or_init(|| polynomial.into_boxed_slice()))
    }

    pub(crate) fn evaluate(&self, model: &mut Model, operation: CurveOp) -> Result<Value> {
        let matrix = match operation {
            CurveOp::Alpha(level) => self.alpha(model, level)?,
            CurveOp::Ratio => self.semilinear_ratio(model)?,
            CurveOp::Product => self.ordered_product(model)?,
        };
        Ok(Value::HWWMatrix(matrix))
    }

    pub(crate) fn report(
        &self,
        model: &mut Model,
        name: &str,
        query: Query,
        candidate: Option<&Matrix>,
    ) -> Result<String> {
        match query {
            Query::Weil => return Ok(format!("rvHWWWeil[{name}] = {:?}", self.weil(model)?)),
            Query::Points(r) => {
                let (q, trace) = self.count_trace(model, r)?;
                let points = i128::from(q) + 1 - trace;
                return Ok(format!(
                    "rvHWWPoints[{name}, {r}] = {points}, trace={trace}"
                ));
            }
            Query::Verify => {}
        }
        let ring = self.ring.as_ref();
        let product;
        let matrix = match candidate {
            Some(matrix) => matrix,
            None => {
                product = self.ordered_product(model)?;
                &product
            }
        };
        let characteristic = hw_algebra::characteristic(model, ring, matrix)?;
        if !ring.is_unit(characteristic.last().unwrap()) {
            return Err(Error::plain(
                "Hasse–Witt unit-factor candidate must have unit determinant",
            ));
        }
        let unit = scalar_coefficients(&characteristic)?;
        let modulus = i128::from(ring.modulus());
        let weil = self.weil(model)?;
        let dividend: Vec<_> = weil
            .iter()
            .map(|a| ring.scalar(a.rem_euclid(modulus) as u64))
            .collect();
        let (quotient, remainder) =
            hw_algebra::divide_monic(model, ring, &dividend, &characteristic)?;
        if remainder.iter().any(|c| !ring.is_zero(c)) {
            return Err(Error::plain(
                "ordered Hasse–Witt factor does not divide the point-count Weil polynomial",
            ));
        }
        let quotient = scalar_coefficients(&quotient)?;
        Ok(format!(
            "rvHWWVerify[{name}] mod {modulus} = unit={unit:?}, quotient={quotient:?}"
        ))
    }
}

fn scalar_coefficients(values: &[Coefficient]) -> Result<Vec<u64>> {
    values.iter().map(|x| rv_unram::require_scalar(x)).collect()
}

fn nonordinary() -> Error {
    Error::plain(
        "hyperelliptic curve is nonordinary: the first Hasse–Witt matrix is singular modulo p",
    )
}

fn require_odd(ring: &Ring) -> Result<()> {
    if ring.prime() == 2 {
        return Err(Error::plain(
            "hyperelliptic Hasse–Witt requires odd residue characteristic",
        ));
    }
    Ok(())
}

fn strip_leading(field: &Ring, mut polynomial: Vec<Coefficient>) -> Vec<Coefficient> {
    let start = polynomial
        .iter()
        .position(|x| !field.is_zero(x))
        .unwrap_or(polynomial.len());
    polynomial.drain(..start);
    polynomial
}

fn require_squarefree(model: &mut Model, field: &Ring, polynomial: &[Coefficient]) -> Result<()> {
    let mut a: Vec<_> = polynomial.iter().rev().map(|c| field.residue(c)).collect();
    let derivative = polynomial
        .iter()
        .enumerate()
        .skip(1)
        .rev()
        .map(|(i, c)| field.scale(model, &field.residue(c), i as u64 % field.prime()))
        .collect::<Result<Vec<_>>>()?;
    let mut b = strip_leading(field, derivative);
    while !b.is_empty() {
        let inverse = field.inverse(model, &b[0])?;
        for coefficient in &mut b {
            *coefficient = field.multiply(model, coefficient, &inverse)?;
        }
        let (_, next) = hw_algebra::divide_monic(model, field, &a, &b)?;
        a = b;
        b = strip_leading(field, next);
    }
    if a.len() != 1 {
        return Err(Error::plain("hyperelliptic reduction must be squarefree"));
    }
    Ok(())
}

#[derive(Clone, Copy)]
pub(crate) enum CurveOp {
    Alpha(usize),
    Ratio,
    Product,
}

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Weil,
    Points(usize),
    Verify,
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Frobenius(usize),
    Inverse,
    Entry(usize, usize),
    Characteristic(usize),
    Descend,
    Negate,
}

impl Unary {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Frobenius(_) => "rvHWWFrob",
            Self::Inverse => "rvHWWInv",
            Self::Entry(_, _) => "rvHWWEntry",
            Self::Characteristic(_) => "rvHWWChar",
            Self::Descend => "rvHWWDescend",
            Self::Negate => "neg",
        }
    }
}

pub(crate) fn frobenius(
    model: &mut Model,
    ring: &Ring,
    matrix: &Matrix,
    power: usize,
) -> Result<Matrix> {
    let entries = matrix
        .entries()
        .iter()
        .map(|x| ring.frobenius(model, x, power))
        .collect::<Result<Vec<_>>>()?;
    Matrix::from_entries(matrix.size(), entries)
}

pub(crate) fn reduce(model: &mut Model, ring: &Ring, matrix: &Matrix) -> Result<Value> {
    model.spend(matrix.entries().len() * ring.cells())?;
    let entries = matrix
        .entries()
        .iter()
        .map(|x| x.iter().map(|a| a % ring.modulus()).collect())
        .collect();
    persist(model, ring, Matrix::from_entries(matrix.size(), entries)?)
}

pub(crate) fn zero(model: &mut Model, ring: &Ring, size: usize) -> Result<Value> {
    let count = matrix_size(ring, size)?;
    model.spend(count * ring.cells())?;
    persist(
        model,
        ring,
        Matrix::from_entries(size, vec![ring.zero(); count])?,
    )
}

pub(crate) fn matrix_size(ring: &Ring, size: usize) -> Result<usize> {
    require_odd(ring)?;
    let count = size
        .checked_mul(size)
        .ok_or_else(|| Error::plain("Hasse–Witt matrix size is out of range"))?;
    hw_algebra::storage(count, ring.cells())?;
    Ok(count)
}

pub(crate) fn count(ring: &Ring, size: usize) -> usize {
    ring.count()
        .saturating_pow(u32::try_from(size.saturating_mul(size)).unwrap_or(u32::MAX))
}

pub(crate) fn enumerate(
    model: &mut Model,
    ring: &Ring,
    size: usize,
    count: usize,
) -> Result<Vec<Value>> {
    let entries = matrix_size(ring, size)?;
    let modulus = ring.modulus();
    (0..count)
        .map(|index| {
            let mut index = index as u64;
            let mut values = vec![ring.zero(); entries];
            model.spend(entries * ring.cells())?;
            for value in values.iter_mut().rev() {
                for digit in value.iter_mut().rev() {
                    *digit = index % modulus;
                    index /= modulus;
                }
            }
            persist(model, ring, Matrix::from_entries(size, values)?)
        })
        .collect()
}

pub(crate) fn add(model: &mut Model, ring: &Ring, left: &Matrix, right: &Matrix) -> Result<Value> {
    let entries = left
        .entries()
        .iter()
        .zip(right.entries())
        .map(|(a, b)| ring.add(model, a, b))
        .collect::<Result<Vec<_>>>()?;
    persist(model, ring, Matrix::from_entries(left.size(), entries)?)
}

pub(crate) fn compare(model: &mut Model, left: &Matrix, right: &Matrix) -> Result<Ordering> {
    for (a, b) in left.entries().iter().zip(right.entries()) {
        let ordering = rv_unram::compare(model, a, b)?;
        if !ordering.is_eq() {
            return Ok(ordering);
        }
    }
    Ok(left.size().cmp(&right.size()))
}

pub(crate) fn compose(model: &mut Model, ring: &Ring, a: &Matrix, b: &Matrix) -> Result<Value> {
    let product = hw_algebra::multiply(model, ring, a, b)?;
    persist(model, ring, product)
}

pub(crate) fn apply(
    model: &mut Model,
    ring: &Ring,
    operation: Unary,
    matrix: &Matrix,
) -> Result<Value> {
    let output = match operation {
        Unary::Frobenius(power) => frobenius(model, ring, matrix, power)?,
        Unary::Inverse => hw_algebra::inverse(model, ring, matrix)?
            .ok_or_else(|| Error::plain("Hasse–Witt matrix is not invertible modulo p"))?,
        Unary::Entry(i, j) => return rv_unram::persist(model, matrix.entry(i - 1, j - 1).clone()),
        Unary::Characteristic(i) => {
            let coeff = hw_algebra::characteristic(model, ring, matrix)?;
            let value = coeff
                .into_iter()
                .nth(i)
                .expect("validated coefficient index");
            return rv_unram::persist(model, value);
        }
        Unary::Descend => {
            model.spend(matrix.entries().len() * ring.cells())?;
            let entries = scalar_coefficients(matrix.entries())?;
            let scalar = hw_algebra::Matrix::from_entries(matrix.size(), entries)?;
            return hyperelliptic::persist_matrix(model, scalar);
        }
        Unary::Negate => {
            let entries = matrix
                .entries()
                .iter()
                .map(|entry| ring.negate(model, entry))
                .collect::<Result<_>>()?;
            Matrix::from_entries(matrix.size(), entries)?
        }
    };
    persist(model, ring, output)
}

pub(crate) fn persist(model: &mut Model, ring: &Ring, matrix: Matrix) -> Result<Value> {
    Ok(Value::HWWMatrix(keep(model, ring, matrix)?))
}

fn keep(model: &mut Model, ring: &Ring, matrix: Matrix) -> Result<Rc<Matrix>> {
    model.allocate(1 + matrix.entries().len() * ring.cells())?;
    Ok(Rc::new(matrix))
}
