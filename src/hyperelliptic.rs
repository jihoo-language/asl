use std::cell::{OnceCell, RefCell};
use std::collections::HashMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;
use crate::rv_unram::{self, add_mod, mul_mod, reduce, sub_mod};
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Command {
    Define(Label, usize, Vec<i128>),
    Report(Label, Report),
    Valuation(Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Report {
    Status,
    Points(usize),
    Weil,
    Convergence(usize),
    Certificate(usize),
}

pub(crate) enum Expression {
    Raw(Label, usize, usize),
    Unit(Label, usize),
    Matrix(usize, usize, usize, Vec<i128>),
    Scalar(usize, usize, i128),
    Unary(Unary, Box<ParsedExpr>),
    Compose(Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Entry(usize, usize),
    Coefficient(usize),
    Reduce(usize),
    Invert,
    Negate,
}

impl Unary {
    pub(crate) fn apply(self, model: &mut Model, ring: Residues, input: Value) -> Result<Value> {
        match (self, input) {
            (Self::Entry(i, j), Value::HWMatrix(matrix)) => {
                Ok(Value::HWScalar(*matrix.entry(i - 1, j - 1)))
            }
            (Self::Coefficient(k), Value::HWMatrix(matrix)) => {
                let coefficients = hw_algebra::characteristic(model, &ring, &matrix)?;
                Ok(Value::HWScalar(coefficients[k]))
            }
            (Self::Reduce(s), Value::HWMatrix(matrix)) => {
                ring.reduce_to(s)?.reduce_matrix(model, &matrix)
            }
            (Self::Reduce(s), Value::HWScalar(value)) => {
                Ok(Value::HWScalar(value % ring.reduce_to(s)?.modulus))
            }
            (Self::Invert, Value::HWMatrix(matrix)) => {
                let inverse = hw_algebra::inverse(model, &ring, &matrix)?.ok_or_else(|| {
                    Error::plain("Hasse–Witt matrix has nonunit determinant modulo p")
                })?;
                persist_matrix(model, inverse)
            }
            (Self::Invert, Value::HWScalar(value)) => {
                if !ring.is_unit(&value) {
                    return Err(Error::plain(
                        "Hasse–Witt scalar inverse requires a unit modulo p",
                    ));
                }
                Ok(Value::HWScalar(ring.inverse_unit(model, &value)?))
            }
            (Self::Negate, Value::HWMatrix(matrix)) => ring.negate_matrix(model, &matrix),
            (Self::Negate, Value::HWScalar(value)) => {
                Ok(Value::HWScalar(sub_mod(0, value, ring.modulus)))
            }
            _ => unreachable!(),
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Entry(..) => "rvHWEntry",
            Self::Coefficient(_) => "rvHWChar",
            Self::Reduce(_) => "rvHWReduce",
            Self::Invert => "rvHWInv",
            Self::Negate => "neg",
        }
    }
}

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Residues {
    pub(crate) prime: Prime,
    pub(crate) precision: usize,
    pub(crate) modulus: u64,
}

impl Residues {
    pub(crate) fn checked(model: &mut Model, p: usize, precision: usize) -> Result<Self> {
        let p = odd_prime(model, p)?;
        Self::new(p, precision)
    }

    pub(crate) fn new(prime: Prime, precision: usize) -> Result<Self> {
        if precision == 0 {
            return Err(Error::plain("Hasse–Witt precision must be positive"));
        }
        let modulus = p_power(prime, precision)
            .ok_or_else(|| Error::plain("Hasse–Witt modulus p^s exceeds u64"))?;
        Ok(Self {
            prime,
            precision,
            modulus,
        })
    }

    pub(crate) fn coefficient(self, n: i128) -> u64 {
        reduce(n, self.modulus)
    }

    pub(crate) fn matrix_size(self, size: usize) -> Result<usize> {
        let entries = size
            .checked_mul(size)
            .ok_or_else(|| Error::plain("Hasse–Witt matrix size is out of range"))?;
        hw_algebra::storage(entries, self.cells())?;
        Ok(entries)
    }

    pub(crate) fn matrix_literal(
        self,
        model: &mut Model,
        size: usize,
        entries: &[i128],
    ) -> Result<Value> {
        if self.matrix_size(size)? != entries.len() {
            return Err(Error::plain(
                "Hasse–Witt matrix has the wrong number of entries",
            ));
        }
        model.spend(entries.len())?;
        let values = entries.iter().map(|&c| self.coefficient(c)).collect();
        persist_matrix(model, Matrix::from_entries(size, values)?)
    }

    pub(crate) fn zero_matrix(self, model: &mut Model, size: usize) -> Result<Value> {
        let entries = self.matrix_size(size)?;
        model.spend(entries)?;
        persist_matrix(model, Matrix::from_entries(size, vec![0; entries])?)
    }

    pub(crate) fn reduce_to(self, precision: usize) -> Result<Self> {
        if precision > self.precision {
            return Err(Error::plain("rvHWReduce requires a no greater precision"));
        }
        Self::new(self.prime, precision)
    }

    pub(crate) fn reduce_matrix(self, model: &mut Model, matrix: &Matrix<u64>) -> Result<Value> {
        model.spend(matrix.entries().len())?;
        let entries = matrix.entries().iter().map(|c| c % self.modulus).collect();
        persist_matrix(model, Matrix::from_entries(matrix.size(), entries)?)
    }

    pub(crate) fn add_matrices(
        self,
        model: &mut Model,
        a: &Matrix<u64>,
        b: &Matrix<u64>,
    ) -> Result<Value> {
        model.spend(a.entries().len())?;
        let entries = a
            .entries()
            .iter()
            .zip(b.entries())
            .map(|(&x, &y)| add_mod(x, y, self.modulus))
            .collect();
        persist_matrix(model, Matrix::from_entries(a.size(), entries)?)
    }

    pub(crate) fn negate_matrix(self, model: &mut Model, input: &Matrix<u64>) -> Result<Value> {
        model.spend(input.entries().len())?;
        let entries = input
            .entries()
            .iter()
            .map(|&x| sub_mod(0, x, self.modulus))
            .collect();
        persist_matrix(model, Matrix::from_entries(input.size(), entries)?)
    }

    pub(crate) fn valuation(self, model: &mut Model, mut coefficient: u64) -> Result<usize> {
        if coefficient == 0 {
            return Ok(self.precision);
        }
        let p = u64::from(self.prime.get());
        let mut order = 0;
        while coefficient.is_multiple_of(p) {
            model.spend(1)?;
            coefficient /= p;
            order += 1;
        }
        Ok(order)
    }

    pub(crate) fn count(self, width: usize) -> usize {
        let radix = usize::try_from(self.modulus).unwrap_or(usize::MAX);
        radix.saturating_pow(width as u32)
    }

    pub(crate) fn matrix_elements(
        self,
        model: &mut Model,
        size: usize,
        limit: usize,
    ) -> Result<Vec<Value>> {
        let width = self.matrix_size(size)?;
        let count = self.count(width);
        if count > limit {
            return Err(Error::plain(
                "Hasse–Witt matrix type exceeds the finite element limit",
            ));
        }
        model.spend(count * (width + 1))?;
        model.allocate(count)?;
        let mut output = Vec::with_capacity(count);
        for index in 0..count {
            let mut index = index as u64;
            let mut entries = vec![0; width];
            for entry in entries.iter_mut().rev() {
                *entry = index % self.modulus;
                index /= self.modulus;
            }
            output.push(persist_matrix(model, Matrix::from_entries(size, entries)?)?);
        }
        Ok(output)
    }
}

impl CoefficientRing for Residues {
    type Element = u64;

    fn cells(&self) -> usize {
        1
    }

    fn zero(&self) -> u64 {
        0
    }

    fn one(&self) -> u64 {
        1
    }

    fn is_zero(&self, value: &u64) -> bool {
        *value == 0
    }

    fn is_unit(&self, value: &u64) -> bool {
        !value.is_multiple_of(u64::from(self.prime.get()))
    }

    fn add(&self, _: &mut Model, a: &u64, b: &u64) -> Result<u64> {
        Ok(add_mod(*a, *b, self.modulus))
    }

    fn subtract(&self, _: &mut Model, a: &u64, b: &u64) -> Result<u64> {
        Ok(sub_mod(*a, *b, self.modulus))
    }

    fn multiply(&self, _: &mut Model, a: &u64, b: &u64) -> Result<u64> {
        Ok(mul_mod(*a, *b, self.modulus))
    }

    fn inverse_unit(&self, model: &mut Model, a: &u64) -> Result<u64> {
        let p = u64::from(self.prime.get());
        let mut inverse = u64::from(invert(model, self.prime, (a % p) as u32)?);
        let mut precision = 1;
        while precision < self.precision {
            model.spend(2)?;
            let product = mul_mod(*a, inverse, self.modulus);
            inverse = mul_mod(inverse, sub_mod(2, product, self.modulus), self.modulus);
            precision *= 2;
        }
        Ok(inverse)
    }
}

fn odd_prime(model: &mut Model, p: usize) -> Result<Prime> {
    if p < 3 || p.is_multiple_of(2) {
        return Err(Error::plain(
            "hyperelliptic Hasse–Witt calculations require an odd prime",
        ));
    }
    Prime::new(model, p)
}

fn p_power(prime: Prime, exponent: usize) -> Option<u64> {
    u32::try_from(exponent)
        .ok()
        .and_then(|e| u64::from(prime.get()).checked_pow(e))
}

pub(crate) fn persist_matrix(model: &mut Model, matrix: Matrix<u64>) -> Result<Value> {
    Ok(Value::HWMatrix(keep_matrix(model, matrix)?))
}

fn keep_matrix(model: &mut Model, matrix: Matrix<u64>) -> Result<Rc<Matrix<u64>>> {
    model.allocate(3 + matrix.entries().len())?;
    Ok(Rc::new(matrix))
}

type MatrixCache = RefCell<HashMap<(usize, usize), Rc<Matrix<u64>>>>;

pub(crate) struct Curve {
    pub(crate) prime: Prime,
    pub(crate) genus: usize,
    lift: Box<[i128]>,
    ordinary: OnceCell<bool>,
    raw: MatrixCache,
    units: RefCell<HashMap<usize, Rc<Matrix<u64>>>>,
    traces: RefCell<HashMap<usize, (u64, i64)>>,
    weil: OnceCell<Rc<[i128]>>,
}

impl Curve {
    pub(crate) fn new(model: &mut Model, prime: usize, lift: Vec<i128>) -> Result<Self> {
        let prime = odd_prime(model, prime)?;
        if lift.len() < 2 || !lift.len().is_multiple_of(2) || lift.last() != Some(&1) {
            return Err(Error::plain(
                "hyperelliptic equation requires a monic odd-degree polynomial",
            ));
        }
        hw_algebra::storage(lift.len(), 3)?;
        model.spend(lift.len())?;
        let p = u64::from(prime.get());
        let polynomial: Vec<_> = lift.iter().map(|&c| reduce(c, p) as u32).collect();
        let derivative = polynomial
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, &c)| prime.mul(prime.reduce(i), c))
            .collect();
        if !rv_unram::coprime(model, prime, polynomial, derivative)? {
            return Err(Error::plain(
                "hyperelliptic equation has repeated roots modulo p",
            ));
        }
        model.allocate(8 + 2 * lift.len())?;
        Ok(Self {
            prime,
            genus: (lift.len() - 2) / 2,
            lift: lift.into(),
            ordinary: OnceCell::new(),
            raw: RefCell::default(),
            units: RefCell::default(),
            traces: RefCell::default(),
            weil: OnceCell::new(),
        })
    }

    pub(crate) fn raw(
        &self,
        model: &mut Model,
        level: usize,
        ring: Residues,
    ) -> Result<Rc<Matrix<u64>>> {
        if ring.prime != self.prime {
            return Err(Error::plain(
                "Hasse–Witt matrix ring has the wrong residue characteristic",
            ));
        }
        let key = (level, ring.precision);
        if let Some(value) = self.raw.borrow().get(&key) {
            return Ok(value.clone());
        }
        let q = p_power(self.prime, level)
            .ok_or_else(|| Error::plain("Hasse–Witt truncation p^r exceeds u64"))?;
        model.spend(self.lift.len())?;
        let polynomial: Vec<_> = self.lift.iter().map(|&c| ring.coefficient(c)).collect();
        let matrix = hw_algebra::hasse_witt(model, &ring, &polynomial, q)?;
        let value = keep_matrix(model, matrix)?;
        self.raw.borrow_mut().insert(key, value.clone());
        Ok(value)
    }

    fn ordinary(&self, model: &mut Model) -> Result<bool> {
        if let Some(value) = self.ordinary.get() {
            return Ok(*value);
        }
        let ring = Residues::new(self.prime, 1)?;
        let a1 = self.raw(model, 1, ring)?;
        let value = hw_algebra::inverse(model, &ring, &a1)?.is_some();
        Ok(*self.ordinary.get_or_init(|| value))
    }

    pub(crate) fn unit(&self, model: &mut Model, ring: Residues) -> Result<Rc<Matrix<u64>>> {
        if ring.prime != self.prime {
            return Err(Error::plain(
                "Hasse–Witt matrix ring has the wrong residue characteristic",
            ));
        }
        let s = ring.precision;
        if let Some(unit) = self.units.borrow().get(&s) {
            return Ok(unit.clone());
        }
        if !self.ordinary(model)? {
            return Err(Error::plain(
                "smooth nonordinary curve has singular Hasse–Witt matrix modulo p",
            ));
        }
        let a = self.raw(model, s, ring)?;
        let b = self.raw(model, s - 1, ring)?;
        let inverse = hw_algebra::inverse(model, &ring, &b)?
            .ok_or_else(|| Error::plain("higher Hasse–Witt denominator has nonunit determinant"))?;
        let unit = hw_algebra::multiply(model, &ring, &a, &inverse)?;
        let unit = keep_matrix(model, unit)?;
        self.units.borrow_mut().insert(s, unit.clone());
        Ok(unit)
    }

    fn trace(&self, model: &mut Model, degree: usize) -> Result<(u64, i64)> {
        if let Some(trace) = self.traces.borrow().get(&degree) {
            return Ok(*trace);
        }
        let field = rv_unram::Ring::field_of_degree(model, self.prime, degree)?;
        let q = field.size_of_residue_field()?;
        let p = field.prime();
        model.spend(self.lift.len())?;
        let coefficients: Vec<_> = self.lift.iter().map(|&c| reduce(c, p)).collect();
        let trace = field.residue_polynomial_trace(model, &coefficients)?;
        model.allocate(4)?;
        self.traces.borrow_mut().insert(degree, (q, trace));
        Ok((q, trace))
    }

    fn weil(&self, model: &mut Model) -> Result<Rc<[i128]>> {
        if let Some(weil) = self.weil.get() {
            return Ok(weil.clone());
        }
        let traces = (1..=self.genus)
            .map(|r| self.trace(model, r).map(|(_, trace)| i128::from(trace)))
            .collect::<Result<Vec<_>>>()?;
        let p = u64::from(self.prime.get());
        let coefficients = weil_polynomial(model, p, &traces)?;
        model.allocate(4 + 2 * coefficients.len())?;
        Ok(self.weil.get_or_init(|| coefficients.into()).clone())
    }

    fn convergence(&self, model: &mut Model, ring: Residues) -> Result<u64> {
        if ring.precision < 2 {
            return Err(Error::plain(
                "Hasse–Witt convergence check requires precision at least two",
            ));
        }
        let lower = ring.reduce_to(ring.precision - 1)?;
        let observed = self.unit(model, ring)?;
        let previous = self.unit(model, lower)?;
        model.spend(observed.entries().len())?;
        if observed
            .entries()
            .iter()
            .zip(previous.entries())
            .any(|(&a, &b)| a % lower.modulus != b)
        {
            return Err(Error::plain(
                "Hasse–Witt matrices fail the requested precision congruence",
            ));
        }
        Ok(lower.modulus)
    }

    fn certificate(&self, model: &mut Model, ring: Residues) -> Result<Certificate> {
        let unit = self.unit(model, ring)?;
        let factor = hw_algebra::characteristic(model, &ring, &unit)?;
        if !ring.is_unit(factor.last().expect("monic characteristic polynomial")) {
            return Err(Error::plain(
                "Hasse–Witt characteristic polynomial has nonunit constant coefficient",
            ));
        }
        let weil = self.weil(model)?;
        model.spend(weil.len())?;
        let reduced: Vec<_> = weil.iter().map(|&c| ring.coefficient(c)).collect();
        let (complement, remainder) = hw_algebra::divide_monic(model, &ring, &reduced, &factor)?;
        if remainder.iter().any(|&c| c != 0) {
            return Err(Error::plain(
                "Hasse–Witt unit factor does not divide the point-count Weil polynomial",
            ));
        }
        model.spend(complement.len())?;
        let prime = u64::from(ring.prime.get());
        if complement.iter().skip(1).any(|c| !c.is_multiple_of(prime)) {
            return Err(Error::plain(
                "Hasse–Witt complementary factor does not reduce to T^g modulo p",
            ));
        }
        Ok(Certificate {
            factor,
            complement,
            remainder,
        })
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, report: Report) -> Result<String> {
        Ok(match report {
            Report::Status => format!(
                "rvHWStatus[{name}] = genus={}, ordinary={}",
                self.genus,
                if self.ordinary(model)? { "yes" } else { "no" }
            ),
            Report::Points(r) => {
                let (q, trace) = self.trace(model, r)?;
                let points = i128::from(q) + 1 - i128::from(trace);
                format!("rvHWPoints[{name}, {r}] = {points}, trace={trace}")
            }
            Report::Weil => format!("rvHWWeil[{name}] = {:?}", self.weil(model)?),
            Report::Convergence(s) => {
                let ring = Residues::new(self.prime, s)?;
                let modulus = self.convergence(model, ring)?;
                format!("rvHWConverges[{name}, {s}] mod {modulus} = yes")
            }
            Report::Certificate(s) => {
                let ring = Residues::new(self.prime, s)?;
                let cert = self.certificate(model, ring)?;
                format!(
                    "rvHWCertify[{name}, {s}] mod {} = unit={:?}, complement={:?}, remainder={:?}",
                    ring.modulus, cert.factor, cert.complement, cert.remainder
                )
            }
        })
    }
}

struct Certificate {
    factor: Vec<u64>,
    complement: Vec<u64>,
    remainder: Vec<u64>,
}

fn overflow() -> Error {
    Error::plain("hyperelliptic Weil coefficient exceeds signed 128-bit range")
}

pub(crate) fn weil_polynomial(model: &mut Model, q: u64, traces: &[i128]) -> Result<Vec<i128>> {
    let genus = traces.len();
    hw_algebra::storage(2 * genus + 1, 2)?;
    let mut coefficients = Vec::with_capacity(2 * genus + 1);
    coefficients.push(1i128);
    for i in 1..=genus {
        model.spend(i)?;
        let mut weighted = 0i128;
        for (&a, &t) in coefficients.iter().zip(traces[..i].iter().rev()) {
            weighted = weighted
                .checked_add(a.checked_mul(t).ok_or_else(overflow)?)
                .ok_or_else(overflow)?;
        }
        let denominator = i as i128;
        if weighted % denominator != 0 {
            return Err(Error::plain("point-count Newton division was not exact"));
        }
        let next = (weighted / denominator)
            .checked_neg()
            .ok_or_else(overflow)?;
        coefficients.push(next);
    }
    let mut scale = 1i128;
    for i in (0..genus).rev() {
        model.spend(1)?;
        scale = scale.checked_mul(i128::from(q)).ok_or_else(overflow)?;
        coefficients.push(coefficients[i].checked_mul(scale).ok_or_else(overflow)?);
    }
    Ok(coefficients)
}
