use std::cell::OnceCell;
use std::cmp::Ordering;
use std::mem::size_of;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};

use super::hecke::{self, compare_integers, integer, Operator};
use super::quadratic_hecke;
use super::{decimal, multiply as product, Integer, Magnitude, INTEGER_WORDS};

type DoubleInteger = Wide<33>;

pub(crate) enum Statement {
    Declare(Label, Box<ParsedExpr>),
    Data(Label),
    Report(Report, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Report {
    Trace,
    Norm,
    Coordinate(u128),
}

impl Report {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Trace => "rvHeckeEigenTrace",
            Self::Norm => "rvHeckeEigenNorm",
            Self::Coordinate(_) => "rvHeckeEigenCoord",
        }
    }
}

pub(crate) enum Expression {
    Image(Label, Box<ParsedExpr>),
    Coefficient(Label, u128),
    Generator(Label),
    Scalar(Label, i128, u128),
    Multiply(Box<[ParsedExpr; 2]>),
}

pub(crate) struct Packet {
    pub(crate) weight: u128,
    generator: Rc<Operator>,
    // Descending coefficients of the monic characteristic polynomial.
    polynomial: Vec<Integer>,
    moments: Vec<Integer>,
    gram: Elimination,
    quadratic_order: OnceCell<quadratic_hecke::Order>,
}

impl Packet {
    pub(crate) fn new(model: &mut Model, weight: u128, generator: Rc<Operator>) -> Result<Self> {
        let dimension = hecke::dimension(weight);
        if dimension == 0 {
            return Err(Error::plain(
                "Hecke eigenpacket requires a nonzero cusp space",
            ));
        }
        if dimension > 1 && generator.scalar().is_some() {
            return Err(nonseparating());
        }
        let degree = usize::try_from(dimension)
            .ok()
            .filter(|&d| d.saturating_mul(d) <= MAX_PIVOT_CELLS / INTEGER_WORDS)
            .ok_or_else(storage)?;
        let polynomial = generator.characteristic(model, weight, dimension)?;
        let moments = power_sums(model, &polynomial)?;
        let matrix = (0..degree)
            .map(|i| moments[i..i + degree].to_vec())
            .collect();
        // H = V^t V, where V is the Vandermonde matrix of the real eigenvalues.
        // H is nonsingular exactly when the characteristic polynomial is squarefree.
        let gram = Elimination::new(model, matrix)?.ok_or_else(nonseparating)?;
        model.allocate(INTEGER_WORDS * (degree * degree + moments.len() + polynomial.len()))?;
        Ok(Self {
            weight,
            generator,
            polynomial,
            moments,
            gram,
            quadratic_order: OnceCell::new(),
        })
    }

    pub(super) fn polynomial(&self) -> &[Integer] {
        &self.polynomial
    }

    pub(crate) fn quadratic_order(&self, model: &mut Model) -> Result<&quadratic_hecke::Order> {
        if let Some(order) = self.quadratic_order.get() {
            return Ok(order);
        }
        let order = quadratic_hecke::Order::new(model, self)?;
        Ok(self.quadratic_order.get_or_init(|| order))
    }

    pub(crate) fn degree(&self) -> usize {
        self.polynomial.len() - 1
    }

    pub(crate) fn data(&self, name: &str) -> String {
        let coefficients: Vec<_> = self.polynomial.iter().copied().map(decimal).collect();
        format!(
            "rvHeckeEigenData[{name}] = weight {}, [{}]",
            self.weight,
            coefficients.join(", ")
        )
    }

    pub(crate) fn image(&self, model: &mut Model, operator: &Operator) -> Result<Value> {
        self.image_element(model, operator)?.persist(model)
    }

    pub(super) fn image_element(&self, model: &mut Model, operator: &Operator) -> Result<Element> {
        if let Some(scalar) = operator.scalar() {
            let mut value = Element::zero(self.degree());
            value.numerator[0] = canonical(scalar);
            return Ok(value);
        }
        if operator.compare(model, &self.generator)?.is_eq() {
            return Ok(self.generator_element());
        }
        let mut trace_cache = hecke::Support::new();
        let mut mixed = operator.clone();
        let mut traces = Vec::with_capacity(self.degree());
        for i in 0..self.degree() {
            traces.push(mixed.trace_cached(model, self.weight, &mut trace_cache)?);
            if i + 1 < self.degree() {
                mixed = mixed.multiply(model, self.weight, &self.generator)?;
            }
        }
        self.gram.solve(model, traces)
    }

    pub(crate) fn element_generator(&self, model: &mut Model) -> Result<Value> {
        self.generator_element().persist(model)
    }

    fn generator_element(&self) -> Element {
        let mut value = Element::zero(self.degree());
        if self.degree() == 1 {
            value.numerator[0] = canonical(self.polynomial[1].negate());
        } else {
            value.numerator[1] = integer(1);
        }
        value
    }

    pub(crate) fn scalar(&self, model: &mut Model, n: i128, d: u128) -> Result<Value> {
        self.constant(model, integer(n), Magnitude::from_small(d))
    }

    fn constant(
        &self,
        model: &mut Model,
        numerator: Integer,
        denominator: Magnitude,
    ) -> Result<Value> {
        let mut value = Element::zero(self.degree());
        value.numerator[0] = numerator;
        Element::new(model, value.numerator, denominator)?.persist(model)
    }

    pub(crate) fn multiply(&self, model: &mut Model, a: &Element, b: &Element) -> Result<Value> {
        let degree = self.degree();
        let mut coefficients = vec![integer(0); 2 * degree - 1];
        for (i, &left) in a.numerator.iter().enumerate() {
            for (j, &right) in b.numerator.iter().enumerate() {
                coefficients[i + j] = add_product(model, coefficients[i + j], left, right)?;
            }
        }
        for i in (degree..coefficients.len()).rev() {
            let leading = coefficients[i];
            for (j, &f) in self.polynomial[1..].iter().enumerate() {
                coefficients[i - j - 1] =
                    add_product(model, coefficients[i - j - 1], leading, f.negate())?;
            }
        }
        coefficients.truncate(degree);
        let denominator =
            product(model, positive(a.denominator), positive(b.denominator))?.magnitude;
        Element::new(model, coefficients, denominator)?.persist(model)
    }

    pub(crate) fn report(
        &self,
        model: &mut Model,
        report: Report,
        value: &Element,
    ) -> Result<String> {
        let name = report.name();
        match report {
            Report::Trace => {
                let mut numerator = integer(0);
                for (&coefficient, &trace) in value.numerator.iter().zip(&self.moments) {
                    numerator = add_product(model, numerator, coefficient, trace)?;
                }
                Ok(format!(
                    "{name} = {}",
                    fraction(model, numerator, value.denominator)?
                ))
            }
            Report::Norm => {
                let d = self.degree();
                // Columns of multiplication by the integral numerator, in the power basis.
                let mut matrix = vec![vec![integer(0); d]; d];
                let mut column = value.numerator.clone();
                for j in 0..d {
                    for (row, &coefficient) in matrix.iter_mut().zip(&column) {
                        row[j] = coefficient;
                    }
                    if j + 1 < d {
                        let leading = column[d - 1];
                        column.rotate_right(1);
                        column[0] = integer(0);
                        for (coefficient, &f) in
                            column.iter_mut().zip(self.polynomial[1..].iter().rev())
                        {
                            *coefficient = add_product(model, *coefficient, leading, f.negate())?;
                        }
                    }
                }
                let Some(factors) = Elimination::new(model, matrix)? else {
                    return Ok(format!("{name} = 0"));
                };
                let norm = factors.determinant();
                let mut denominator = Magnitude::from_small(1);
                let mut numerator = norm;
                // Cancel one copy of the denominator at a time to avoid forming
                // an oversized denominator for an integral norm.
                for _ in 0..d {
                    let common = numerator.magnitude.gcd(model, value.denominator)?;
                    numerator = quotient(model, numerator, positive(common))?;
                    let rest = quotient(model, positive(value.denominator), positive(common))?;
                    denominator = product(model, positive(denominator), rest)?.magnitude;
                }
                Ok(format!(
                    "{name} = {}",
                    write_fraction(numerator, denominator)
                ))
            }
            Report::Coordinate(i) => {
                let numerator = usize::try_from(i)
                    .ok()
                    .and_then(|j| value.numerator.get(j).copied())
                    .unwrap_or_else(|| integer(0));
                Ok(format!(
                    "{name}[{i}] = {}",
                    fraction(model, numerator, value.denominator)?
                ))
            }
        }
    }
}

#[derive(Clone)]
pub(crate) struct Element {
    numerator: Vec<Integer>,
    denominator: Magnitude,
}

impl Element {
    fn new(model: &mut Model, mut numerator: Vec<Integer>, denominator: Magnitude) -> Result<Self> {
        if denominator.is_zero() {
            return Err(Error::plain(
                "Hecke eigenpacket scalar denominator must be positive",
            ));
        }
        let mut common = denominator;
        let one = Magnitude::from_small(1);
        for coefficient in &mut numerator {
            *coefficient = canonical(*coefficient);
            if common != one {
                common = common.gcd(model, coefficient.magnitude)?;
            }
        }
        let denominator = if common == one {
            denominator
        } else {
            for coefficient in &mut numerator {
                *coefficient = quotient(model, *coefficient, positive(common))?;
            }
            quotient(model, positive(denominator), positive(common))?.magnitude
        };
        Ok(Self {
            numerator,
            denominator,
        })
    }

    pub(crate) fn zero(degree: usize) -> Self {
        Self {
            numerator: vec![integer(0); degree],
            denominator: Magnitude::from_small(1),
        }
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(INTEGER_WORDS * (self.numerator.len() + 1))?;
        Ok(Value::HeckeEigen(Rc::new(self)))
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let common = positive(self.denominator.gcd(model, rhs.denominator)?);
        let left = quotient(model, positive(rhs.denominator), common)?;
        let right = quotient(model, positive(self.denominator), common)?;
        let mut numerator = Vec::with_capacity(self.numerator.len());
        for (&a, &b) in self.numerator.iter().zip(&rhs.numerator) {
            let first = product(model, a, left)?;
            numerator.push(add_product(model, first, b, right)?);
        }
        let denominator = product(model, positive(self.denominator), left)?.magnitude;
        Self::new(model, numerator, denominator)
    }

    pub(crate) fn negate(&self, model: &mut Model) -> Result<Self> {
        model.spend(self.numerator.len() * INTEGER_WORDS)?;
        Ok(Self {
            numerator: self
                .numerator
                .iter()
                .map(|value| canonical(value.negate()))
                .collect(),
            denominator: self.denominator,
        })
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend((self.numerator.len() + 1) * INTEGER_WORDS)?;
        let denominator = self.denominator.cmp(rhs.denominator);
        if !denominator.is_eq() {
            return Ok(denominator);
        }
        Ok(self
            .numerator
            .iter()
            .zip(&rhs.numerator)
            .map(|(&a, &b)| compare_integers(a, b))
            .find(|order| !order.is_eq())
            .unwrap_or(Ordering::Equal))
    }

    pub(super) fn fraction_parts(&self) -> (&[Integer], Magnitude) {
        (&self.numerator, self.denominator)
    }

    pub(crate) fn coordinates(&self) -> impl Iterator<Item = String> + '_ {
        self.numerator.iter().copied().map(decimal)
    }

    pub(crate) fn denominator(&self) -> String {
        decimal(positive(self.denominator))
    }
}

fn power_sums(model: &mut Model, polynomial: &[Integer]) -> Result<Vec<Integer>> {
    let degree = polynomial.len() - 1;
    let mut sums = vec![integer(degree as i128)];
    for index in 1..(2 * degree - 1) {
        let mut sum = if index <= degree {
            product(model, polynomial[index], integer(index as i128))?
        } else {
            integer(0)
        };
        for j in 1..=degree.min(index - 1) {
            sum = add_product(model, sum, polynomial[j], sums[index - j])?;
        }
        sums.push(canonical(sum.negate()));
    }
    Ok(sums)
}

// Bareiss elimination keeps all stored entries integral. The lower triangle
// records the exact row operations used again for each mixed moment vector.
struct Elimination {
    matrix: Vec<Vec<Integer>>,
    swaps: Vec<usize>,
}

impl Elimination {
    fn new(model: &mut Model, mut matrix: Vec<Vec<Integer>>) -> Result<Option<Self>> {
        let d = matrix.len();
        let mut previous = integer(1);
        let mut swaps = Vec::with_capacity(d);
        for k in 0..d {
            let Some(pivot_row) = (k..d).find(|&row| !matrix[row][k].is_zero()) else {
                return Ok(None);
            };
            swaps.push(pivot_row);
            // Keep the factors from previous steps on their original rows so
            // that right-hand sides can replay the swaps in the same order.
            if pivot_row != k {
                let (before, from_pivot) = matrix.split_at_mut(pivot_row);
                before[k][k..].swap_with_slice(&mut from_pivot[0][k..]);
            }
            let pivot = matrix[k][k];
            for i in k + 1..d {
                for j in k + 1..d {
                    matrix[i][j] = difference_quotient(
                        model,
                        matrix[i][j],
                        pivot,
                        matrix[k][j],
                        matrix[i][k],
                        previous,
                    )?;
                }
            }
            previous = pivot;
        }
        Ok(Some(Self { matrix, swaps }))
    }

    fn determinant(&self) -> Integer {
        let d = self.matrix.len();
        let pivot = self.matrix[d - 1][d - 1];
        let odd = self
            .swaps
            .iter()
            .enumerate()
            .filter(|&(i, row)| i != *row)
            .count();
        if odd.is_multiple_of(2) {
            pivot
        } else {
            pivot.negate()
        }
    }

    fn solve(&self, model: &mut Model, mut rhs: Vec<Integer>) -> Result<Element> {
        let d = rhs.len();
        let mut previous = integer(1);
        for (k, &row) in self.swaps.iter().enumerate() {
            rhs.swap(k, row);
            let pivot = self.matrix[k][k];
            for i in k + 1..d {
                rhs[i] =
                    difference_quotient(model, rhs[i], pivot, rhs[k], self.matrix[i][k], previous)?;
            }
            previous = pivot;
        }
        // Solve for det*G instead of introducing a denominator at each row.
        let mut result = vec![integer(0); d];
        for i in (0..d).rev() {
            let pivot = self.matrix[i][i];
            let mut numerator = double_product(model, previous, rhs[i])?;
            for (j, &known) in result.iter().enumerate().skip(i + 1) {
                numerator =
                    numerator.add(double_product(model, self.matrix[i][j], known.negate())?)?;
            }
            result[i] = double_quotient(model, numerator, pivot)?;
        }
        if previous.negative {
            for value in &mut result {
                *value = value.negate();
            }
        }
        Element::new(model, result, previous.magnitude)
    }
}

fn canonical(value: Integer) -> Integer {
    if value.is_zero() {
        integer(0)
    } else {
        value
    }
}

fn positive(magnitude: Magnitude) -> Integer {
    Integer {
        magnitude,
        negative: false,
    }
}

fn add_product(model: &mut Model, sum: Integer, a: Integer, b: Integer) -> Result<Integer> {
    sum.add(product(model, a, b)?).map_err(|_| overflow())
}

fn quotient(model: &mut Model, value: Integer, divisor: Integer) -> Result<Integer> {
    value.divide_exact(model, divisor)?.ok_or_else(inexact)
}

fn difference_quotient(
    model: &mut Model,
    a: Integer,
    b: Integer,
    c: Integer,
    d: Integer,
    divisor: Integer,
) -> Result<Integer> {
    let left = double_product(model, a, b)?;
    let right = double_product(model, c, d)?;
    double_quotient(model, left.add(right.negate())?, divisor)
}

// Keep room for a sum of products before the exact division of a Bareiss entry.
fn double_product(model: &mut Model, a: Integer, b: Integer) -> Result<DoubleInteger> {
    model.spend(1 + a.magnitude.bits().div_ceil(32) * b.magnitude.bits().div_ceil(32))?;
    a.resize::<33>()?.multiply_wide(b.resize()?)
}

fn double_quotient(model: &mut Model, value: DoubleInteger, divisor: Integer) -> Result<Integer> {
    value
        .divide_exact(model, divisor.resize()?)?
        .ok_or_else(inexact)?
        .resize()
        .map_err(|_| overflow())
}

fn fraction(model: &mut Model, numerator: Integer, denominator: Magnitude) -> Result<String> {
    let common = positive(numerator.magnitude.gcd(model, denominator)?);
    let numerator = quotient(model, numerator, common)?;
    let denominator = quotient(model, positive(denominator), common)?.magnitude;
    Ok(write_fraction(numerator, denominator))
}

fn write_fraction(numerator: Integer, denominator: Magnitude) -> String {
    if denominator == Magnitude::from_small(1) {
        decimal(numerator)
    } else {
        format!("{}/{}", decimal(numerator), decimal(positive(denominator)))
    }
}

fn nonseparating() -> Error {
    Error::plain("Hecke eigenpacket generator does not have a squarefree characteristic polynomial")
}
fn storage() -> Error {
    Error::plain("Hecke eigenpacket exceeds the exact matrix storage limit")
}
fn overflow() -> Error {
    Error::plain(format!(
        "Hecke eigenpacket exact intermediate exceeds {} bits",
        8 * size_of::<Magnitude>()
    ))
}
fn inexact() -> Error {
    Error::plain("Hecke eigenpacket integer division is inexact")
}
