use std::cmp::Ordering;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::lee_linear::{Echelon, Rational};
use crate::model::{Model, Value};
use crate::syntax::ParsedExpr;

use super::hecke::{self, Operator};

pub(crate) mod amice;
pub(crate) mod eigenperiod;
pub(crate) mod haberland;
pub(crate) mod mazur_tate;
pub(crate) mod modular_symbol;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum Sign {
    Plus,
    Minus,
}

impl Sign {
    pub(crate) fn text(self) -> &'static str {
        match self {
            Self::Plus => "plus",
            Self::Minus => "minus",
        }
    }

    fn includes(self, degree: usize) -> bool {
        degree.is_multiple_of(2) == (self == Self::Plus)
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) weight: usize,
    pub(crate) sign: Sign,
}

impl Grade {
    pub(crate) fn new(weight: u128, sign: Sign) -> Result<Self> {
        hecke::validate_weight(weight)?;
        let weight = usize::try_from(weight)
            .ok()
            .filter(|&k| k.saturating_mul(k) <= MAX_PIVOT_CELLS / 4)
            .ok_or_else(|| {
                Error::plain("period-polynomial weight exceeds the exact matrix storage limit")
            })?;
        Ok(Self { weight, sign })
    }

    pub(crate) fn is_trivial(self) -> bool {
        hecke::dimension(self.weight as u128) == 0
    }

    fn degree(self) -> usize {
        self.weight - 2
    }
}

pub(crate) enum Expression {
    Basis(Grade, usize),
    IntegralBasis(Grade, usize),
    Polynomial(Grade, Vec<Rational>),
    Scale(Rational, Box<ParsedExpr>),
    Hecke(Box<[ParsedExpr; 2]>),
}

pub(crate) enum Report {
    Dimension(Grade),
    Coefficient(u128, Box<ParsedExpr>),
    Matrix(Grade, Box<ParsedExpr>),
    Pair(Box<[ParsedExpr; 2]>),
    Pairing(Grade),
}

pub(crate) struct Basis {
    free: Vec<usize>,
    vectors: Vec<Polynomial>,
}

impl Basis {
    pub(crate) fn new(model: &mut Model, grade: Grade) -> Result<Self> {
        let w = grade.degree();
        let binomial = binomial_table(model, w)?;
        let mut echelon = Echelon::default();
        let mut vectors = Vec::new();
        let mut free = Vec::new();
        for j in (1..=w / 2).filter(|&j| grade.sign.includes(j)) {
            let mirror = w - j;
            if j == mirror && j.is_multiple_of(2) {
                // Such a middle coefficient must vanish under 1 + S.
                continue;
            }
            let reflected = if j.is_multiple_of(2) { -1 } else { 1 };
            let mut relation = Vec::new();
            for r in 0..=w {
                let mut entry = manin_column(&binomial, w, j, r)?;
                if j != mirror {
                    let extra = manin_column(&binomial, w, mirror, r)?;
                    entry = entry
                        .checked_add(extra.checked_mul(reflected).ok_or_else(arithmetic)?)
                        .ok_or_else(arithmetic)?;
                }
                if entry != 0 {
                    relation.push((r, entry));
                }
            }
            model.spend(1 + w)?;
            let Some(kernel) = echelon.insert(model, relation, vec![(j, 1)])? else {
                continue;
            };
            // Each dependency uses only previous independent columns, so the
            // coefficients at the free indices are our basis coordinates.
            let denominator = kernel.last().expect("dependency has its new column").1;
            let mut polynomial = Polynomial::zero(grade);
            for (j, numerator) in kernel {
                let coefficient =
                    Rational::ratio(numerator, denominator).map_err(|_| arithmetic())?;
                polynomial.0[j] = coefficient;
                if j != w - j {
                    polynomial.0[w - j] = if j.is_multiple_of(2) {
                        negate(coefficient)?
                    } else {
                        coefficient
                    };
                }
            }
            free.push(j);
            vectors.push(polynomial);
        }
        if vectors.len() as u128 != hecke::dimension(grade.weight as u128) {
            return Err(Error::plain(
                "Manin period-polynomial dimension disagrees with the cusp space",
            ));
        }
        model.allocate(1 + vectors.len() * (1 + 4 * (w + 1)))?;
        Ok(Self { free, vectors })
    }

    pub(crate) fn dimension(&self) -> usize {
        self.vectors.len()
    }

    pub(crate) fn vector(&self, index: usize) -> Result<Polynomial> {
        self.vectors
            .get(index)
            .cloned()
            .ok_or_else(|| Error::plain("period-polynomial basis index is outside this cusp space"))
    }

    pub(crate) fn hecke_matrix(
        &self,
        model: &mut Model,
        grade: Grade,
        operator: &Operator,
    ) -> Result<String> {
        let columns: Vec<_> = self
            .vectors
            .iter()
            .map(|p| p.hecke(model, grade, operator))
            .collect::<Result<_>>()?;
        let rows: Vec<_> = self
            .free
            .iter()
            .map(|&j| {
                let entries: Vec<_> = columns
                    .iter()
                    .map(|column| column.0[j].to_string())
                    .collect();
                format!("[{}]", entries.join(", "))
            })
            .collect();
        Ok(format!("rvPeriodHeckeMatrix = [{}]", rows.join(", ")))
    }
}

#[derive(Clone)]
pub(crate) struct Polynomial(Vec<Rational>);

impl Polynomial {
    pub(crate) fn new(model: &mut Model, grade: Grade, coordinates: Vec<Rational>) -> Result<Self> {
        if coordinates.len() != grade.degree() + 1 {
            return Err(Error::plain(
                "period polynomial requires exactly k−1 ascending coefficients",
            ));
        }
        let mut result = Self(coordinates);
        result.check(model, grade)?;
        result.normalize(grade)?;
        Ok(result)
    }

    pub(crate) fn zero(grade: Grade) -> Self {
        Self(vec![Rational::ZERO; grade.degree() + 1])
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + 4 * self.0.len())?;
        Ok(Value::PeriodPolynomial(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(self.0.len())?;
        Ok(self.0.cmp(&rhs.0))
    }

    pub(crate) fn coordinates(&self) -> &[Rational] {
        &self.0
    }

    pub(crate) fn coefficient(&self, i: u128) -> Rational {
        usize::try_from(i)
            .ok()
            .and_then(|j| self.0.get(j).copied())
            .unwrap_or(Rational::ZERO)
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        model.spend(self.0.len())?;
        let sum = self
            .0
            .iter()
            .zip(&rhs.0)
            .map(|(&a, &b)| plus(a, b))
            .collect::<Result<_>>()?;
        Ok(Self(sum))
    }

    pub(crate) fn scale(&self, model: &mut Model, scalar: Rational) -> Result<Self> {
        model.spend(self.0.len())?;
        let entries = self
            .0
            .iter()
            .map(|&value| multiply(value, scalar))
            .collect::<Result<_>>()?;
        Ok(Self(entries))
    }

    pub(crate) fn hecke(
        &self,
        model: &mut Model,
        grade: Grade,
        operator: &Operator,
    ) -> Result<Self> {
        let mut result = self.hecke_full(model, grade, operator)?;
        result.normalize(grade)?;
        Ok(result)
    }

    fn hecke_full(&self, model: &mut Model, grade: Grade, operator: &Operator) -> Result<Self> {
        let mut result = Self::zero(grade);
        if self.0.iter().all(|r| r.numerator == 0) {
            return Ok(result);
        }
        if let Some(scalar) = operator.scalar() {
            return self.scale(
                model,
                Rational::integer(scalar.signed().ok_or_else(arithmetic)?),
            );
        }
        for (index, coefficient) in operator.terms() {
            let scalar = Rational::integer(coefficient.signed().ok_or_else(arithmetic)?);
            let mut image = Self::zero(grade);
            if index == 1 {
                image = self.clone();
            } else {
                merel_matrices(model, index, |model, matrix| {
                    let transformed = self.transform(model, matrix)?;
                    image = image.add(model, &transformed)?;
                    Ok(())
                })?;
            }
            let scaled = image.scale(model, scalar)?;
            result = result.add(model, &scaled)?;
        }
        result.check(model, grade)?;
        Ok(result)
    }

    fn check(&self, model: &mut Model, grade: Grade) -> Result<()> {
        let w = grade.degree();
        for (j, &value) in self.0.iter().enumerate() {
            if value.numerator != 0 && !grade.sign.includes(j) {
                return Err(Error::plain("period polynomial has the wrong X-parity"));
            }
            let reflected = if j.is_multiple_of(2) {
                self.0[w - j]
            } else {
                negate(self.0[w - j])?
            };
            if plus(value, reflected)?.numerator != 0 {
                return Err(Error::plain(
                    "period polynomial violates the Manin relation 1+S",
                ));
            }
        }
        if self.0.iter().any(|value| value.numerator != 0) {
            let u = self.transform(model, [0, -1, 1, 1])?;
            let u_squared = self.transform(model, [-1, -1, 1, 0])?;
            if self
                .add(model, &u)?
                .add(model, &u_squared)?
                .0
                .iter()
                .any(|r| r.numerator != 0)
            {
                return Err(Error::plain(
                    "period polynomial violates the Manin relation 1+U+U²",
                ));
            }
        }
        Ok(())
    }

    fn normalize(&mut self, grade: Grade) -> Result<()> {
        let w = grade.degree();
        if w > 0 && grade.sign == Sign::Plus {
            // P − p_w(X^w − Y^w).
            self.0[0] = plus(self.0[0], self.0[w])?;
            self.0[w] = Rational::ZERO;
        }
        Ok(())
    }

    fn transform(&self, model: &mut Model, [a, b, c, d]: [i128; 4]) -> Result<Self> {
        // Homogeneous Horner expansion by the forms aX+bY and cX+dY.
        let mut coefficients = self.0.iter().rev().copied();
        let mut result = vec![coefficients.next().expect("nonempty polynomial")];
        let mut right_power = vec![Rational::integer(1)];
        for value in coefficients {
            model.spend(1 + 4 * result.len())?;
            result = times_linear(&result, a, b)?;
            right_power = times_linear(&right_power, c, d)?;
            if value.numerator != 0 {
                for (coefficient, &power) in result.iter_mut().zip(&right_power) {
                    *coefficient = plus(*coefficient, multiply(value, power)?)?;
                }
            }
        }
        Ok(Self(result))
    }
}

pub(crate) fn ratio(numerator: i128, denominator: i128) -> Result<Rational> {
    if denominator <= 0 {
        return Err(Error::plain(
            "period-polynomial rational denominator must be positive",
        ));
    }
    Rational::new(numerator, denominator).map_err(|_| arithmetic())
}

fn plus(a: Rational, b: Rational) -> Result<Rational> {
    a.add(b).map_err(|_| arithmetic())
}

fn multiply(a: Rational, b: Rational) -> Result<Rational> {
    a.multiply(b).map_err(|_| arithmetic())
}

fn negate(a: Rational) -> Result<Rational> {
    multiply(a, Rational::integer(-1))
}

fn times_linear(poly: &[Rational], x: i128, y: i128) -> Result<Vec<Rational>> {
    let mut result = vec![Rational::ZERO; poly.len() + 1];
    let x = Rational::integer(x);
    let y = Rational::integer(y);
    for (j, &coefficient) in poly.iter().enumerate() {
        result[j] = plus(result[j], multiply(coefficient, y)?)?;
        result[j + 1] = multiply(coefficient, x)?;
    }
    Ok(result)
}

fn binomial_table(model: &mut Model, w: usize) -> Result<Vec<Vec<i128>>> {
    model.spend((w + 1) * (w + 1))?;
    let mut table = vec![vec![1_i128]];
    for k in 1..=w {
        let previous = &table[k - 1];
        let mut row = vec![1];
        for pair in previous.windows(2) {
            row.push(pair[0].checked_add(pair[1]).ok_or_else(arithmetic)?);
        }
        row.push(1);
        table.push(row);
    }
    Ok(table)
}

fn manin_column(binomial: &[Vec<i128>], w: usize, j: usize, r: usize) -> Result<i128> {
    // U sends (X,Y) to (−Y,X+Y), and U² to (−X−Y,X).
    let first = binomial[w - j].get(r).copied().unwrap_or(0);
    let second = r
        .checked_sub(w - j)
        .and_then(|index| binomial[j].get(index))
        .copied()
        .unwrap_or(0);
    let combined = first.checked_add(second).ok_or_else(arithmetic)?;
    i128::from(j == r)
        .checked_add(if j.is_multiple_of(2) {
            combined
        } else {
            -combined
        })
        .ok_or_else(arithmetic)
}

// In this right substitution convention use the transpose of the common
// Merel set: a>c>=0, d>b>=0. Put h=a-c and l=d-b, so n=a*l+h*b.
fn merel_matrices(
    model: &mut Model,
    n: u128,
    mut visit: impl FnMut(&mut Model, [i128; 4]) -> Result<()>,
) -> Result<()> {
    let n = i128::try_from(n).map_err(|_| {
        Error::plain("period-polynomial Hecke index exceeds the exact matrix-entry range")
    })?;
    // Even before the divisor loops, the a-range alone has n entries.
    let outer_steps = usize::try_from(n)
        .map_err(|_| Error::plain("period-polynomial Hecke index exceeds the computation limit"))?;
    model.spend(outer_steps)?;
    for a in 1..=n {
        for l in 1..=n / a {
            model.spend(1)?;
            let remainder = n - a * l;
            if remainder == 0 {
                for h in 1..=a {
                    model.spend(1)?;
                    visit(model, [a, 0, a - h, l])?;
                }
            } else {
                for h in 1..=a.min(remainder.isqrt()) {
                    model.spend(1)?;
                    if (remainder as u128).is_multiple_of(h as u128) {
                        let other = remainder / h;
                        visit(model, [a, other, a - h, l + other])?;
                        if other != h && other <= a {
                            visit(model, [a, h, a - other, l + h])?;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

fn arithmetic() -> Error {
    Error::plain("period polynomial exact arithmetic exceeds the signed 128-bit range")
}
