use std::cmp::Ordering;
use std::rc::Rc;

use super::unimodular;
use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::lee_linear::Rational;
use crate::model::{Model, Value};
use crate::rv_lens::gcd;
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod llv;

const RANK: usize = 1 + unimodular::RANK;
const MIDDLE: usize = RANK * (RANK + 1) / 2;

pub(crate) enum Expression {
    Basis(Label, i128, usize),
    Divisor(i128, Box<ParsedExpr>),
    Cup(Box<[ParsedExpr; 2]>),
    Scale(i128, i128, Box<ParsedExpr>),
    Chern(Label, bool),
    Llv(llv::Operation, Box<[ParsedExpr; 2]>),
    Rotate(Box<[ParsedExpr; 3]>),
}

pub(crate) fn dimension(degree: i128) -> usize {
    match degree {
        0 | 8 => 1,
        2 | 6 => RANK,
        4 => MIDDLE,
        _ => 0,
    }
}

pub(crate) fn cup_degree(left: i128, right: i128) -> Result<i128> {
    left.checked_add(right).ok_or_else(|| {
        Error::plain("Hilbert-square cup product degree exceeds the signed 128-bit range")
    })
}

pub(crate) fn ratio(numerator: i128, denominator: i128) -> Result<Rational> {
    if denominator == 0 {
        return Err(Error::plain(
            "Hilbert-square rational denominator cannot vanish",
        ));
    }
    checked(Rational::ratio(numerator, denominator))
}

fn checked<T>(result: Result<T>) -> Result<T> {
    result.map_err(|_| {
        Error::plain(
            "Hilbert-square rational arithmetic exceeds the signed 128-bit coefficient range",
        )
    })
}

pub(crate) struct Class {
    degree: i128,
    coordinates: Vec<Rational>,
}

impl Class {
    pub(crate) fn zero(degree: i128) -> Self {
        Self {
            degree,
            coordinates: vec![Rational::ZERO; dimension(degree)],
        }
    }

    pub(crate) fn basis(degree: i128, index: usize) -> Result<Self> {
        let mut value = Self::zero(degree);
        let coordinate = value.coordinates.get_mut(index).ok_or_else(|| {
            Error::plain(format!(
            "rvHilbSquareBasis index {index} is out of range for degree {degree} of dimension {}",
            dimension(degree)
        ))
        })?;
        *coordinate = Rational::integer(1);
        Ok(value)
    }

    pub(crate) fn divisor(
        model: &mut Model,
        class: &unimodular::Vector,
        delta: i128,
    ) -> Result<Self> {
        let coordinates = class.basis_coordinates(model)?.into_iter().chain([delta]);
        Ok(Self {
            degree: 2,
            coordinates: coordinates.map(Rational::integer).collect(),
        })
    }

    pub(crate) fn coordinates(&self) -> &[Rational] {
        &self.coordinates
    }

    pub(crate) fn integral(&self) -> Rational {
        debug_assert_eq!(self.degree, 8);
        self.coordinates[0]
    }

    pub(crate) fn c4() -> Self {
        Self {
            degree: 8,
            coordinates: vec![Rational::integer(324)],
        }
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(2 + 2 * self.coordinates.len())?;
        Ok(Value::HilbertSquare(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        model.spend(1 + self.coordinates.len())?;
        Ok(self.coordinates.cmp(&other.coordinates))
    }

    pub(crate) fn add(&self, model: &mut Model, other: &Self) -> Result<Self> {
        model.spend(self.coordinates.len())?;
        let coordinates = self
            .coordinates
            .iter()
            .zip(&other.coordinates)
            .map(|(&a, &b)| checked(a.add(b)))
            .collect::<Result<_>>()?;
        Ok(Self {
            degree: self.degree,
            coordinates,
        })
    }

    pub(crate) fn scale(&self, model: &mut Model, scalar: Rational) -> Result<Self> {
        model.spend(self.coordinates.len())?;
        let coordinates = self
            .coordinates
            .iter()
            .map(|&value| checked(value.multiply(scalar)))
            .collect::<Result<_>>()?;
        Ok(Self {
            degree: self.degree,
            coordinates,
        })
    }

    pub(crate) fn cup(&self, model: &mut Model, form: &Form, other: &Self) -> Result<Self> {
        let (left, right) = if self.degree <= other.degree {
            (self, other)
        } else {
            (other, self)
        };
        let degree = cup_degree(left.degree, right.degree)?;
        let mut result = Self::zero(degree);
        model.spend(left.coordinates.len() + right.coordinates.len())?;
        if result.coordinates.is_empty()
            || nonzero(&left.coordinates).next().is_none()
            || nonzero(&right.coordinates).next().is_none()
        {
            return Ok(result);
        }
        match (left.degree, right.degree) {
            (0, _) => return right.scale(model, left.coordinates[0]),
            (2, 2) => {
                model.spend(RANK * RANK)?;
                for (i, a) in nonzero(&left.coordinates) {
                    for (j, b) in nonzero(&right.coordinates) {
                        let position = monomial_index(i.min(j), i.max(j));
                        accumulate(&mut result.coordinates[position], a, b)?;
                    }
                }
            }
            (2, 4) => {
                let contractions = form.contract(model, &left.coordinates)?;
                let mut trace = Rational::ZERO;
                model.spend(MIDDLE * 3 + RANK)?;
                for (index, coefficient) in nonzero(&right.coordinates) {
                    let (i, j) = form.monomials[index];
                    accumulate(&mut trace, coefficient, Rational::integer(form.gram[i][j]))?;
                    accumulate(&mut result.coordinates[j], coefficient, contractions[i])?;
                    accumulate(&mut result.coordinates[i], coefficient, contractions[j])?;
                }
                for (out, &coordinate) in result.coordinates.iter_mut().zip(&left.coordinates) {
                    accumulate(out, trace, coordinate)?;
                }
            }
            (2, 6) => {
                let contractions = form.contract(model, &left.coordinates)?;
                model.spend(RANK)?;
                for (index, coefficient) in nonzero(&right.coordinates) {
                    accumulate(&mut result.coordinates[0], contractions[index], coefficient)?;
                }
            }
            (4, 4) => {
                model.spend(MIDDLE * MIDDLE)?;
                for (first, a) in nonzero(&left.coordinates) {
                    let (i, j) = form.monomials[first];
                    for (second, b) in nonzero(&right.coordinates) {
                        let (k, l) = form.monomials[second];
                        let q = &form.gram;
                        let pairing = q[i][j] * q[k][l] + q[i][k] * q[j][l] + q[i][l] * q[j][k];
                        if pairing != 0 {
                            let term = multiply_weighted(a, b, pairing)?;
                            result.coordinates[0] = checked(result.coordinates[0].add(term))?;
                        }
                    }
                }
            }
            _ => unreachable!("only supported degrees remain"),
        }
        Ok(result)
    }
}

fn nonzero(coordinates: &[Rational]) -> impl Iterator<Item = (usize, Rational)> + '_ {
    coordinates
        .iter()
        .copied()
        .enumerate()
        .filter(|(_, value)| value.numerator != 0)
}

fn accumulate(out: &mut Rational, a: Rational, b: Rational) -> Result<()> {
    if a.numerator != 0 && b.numerator != 0 {
        *out = checked(out.add(checked(a.multiply(b))?))?;
    }
    Ok(())
}

fn multiply_weighted(mut a: Rational, mut b: Rational, mut weight: i128) -> Result<Rational> {
    // Cancel against both denominators before multiplying a weighted term.
    for factor in [&mut a, &mut b] {
        let divisor = gcd(weight.unsigned_abs(), factor.denominator as u128) as i128;
        weight /= divisor;
        factor.denominator /= divisor;
    }
    let product = checked(a.multiply(b))?;
    checked(product.multiply(Rational::integer(weight)))
}

fn monomial_index(i: usize, j: usize) -> usize {
    i * RANK - i * i.saturating_sub(1) / 2 + (j - i)
}

pub(crate) struct Form {
    gram: [[i128; RANK]; RANK],
    monomials: Vec<(usize, usize)>,
}

impl Form {
    pub(crate) fn new(model: &mut Model) -> Result<Self> {
        let mut gram = [[0; RANK]; RANK];
        let surface = unimodular::integral_gram(model)?;
        for (row, surface_row) in gram.iter_mut().zip(surface) {
            row[..unimodular::RANK].copy_from_slice(&surface_row);
        }
        gram[unimodular::RANK][unimodular::RANK] = -2;
        model.allocate(RANK * RANK + 2 * MIDDLE)?;
        let monomials = (0..RANK)
            .flat_map(|i| (i..RANK).map(move |j| (i, j)))
            .collect();
        Ok(Self { gram, monomials })
    }

    fn contract(&self, model: &mut Model, vector: &[Rational]) -> Result<[Rational; RANK]> {
        model.spend(RANK * RANK)?;
        let mut pairings = [Rational::ZERO; RANK];
        for (out, row) in pairings.iter_mut().zip(&self.gram) {
            for (&pairing, &coordinate) in row.iter().zip(vector) {
                accumulate(out, Rational::integer(pairing), coordinate)?;
            }
        }
        Ok(pairings)
    }

    pub(crate) fn c2(&self, model: &mut Model) -> Result<Class> {
        let entries = self
            .gram
            .iter()
            .flatten()
            .copied()
            .map(Rational::integer)
            .collect();
        let matrix = Matrix::from_entries(RANK, entries)?;
        let inverse = hw_algebra::inverse(model, &RationalField, &matrix)?
            .ok_or_else(|| Error::plain("Hilbert-square Beauville form is degenerate"))?;
        model.spend(MIDDLE)?;
        let coordinates = self
            .monomials
            .iter()
            .map(|&(i, j)| {
                let numerator = if i == j { 6 } else { 12 };
                let scale = ratio(numerator, 5)?;
                checked(inverse.entry(i, j).multiply(scale))
            })
            .collect::<Result<_>>()?;
        Ok(Class {
            degree: 4,
            coordinates,
        })
    }
}

struct RationalField;

impl CoefficientRing for RationalField {
    type Element = Rational;

    fn cells(&self) -> usize {
        2
    }
    fn zero(&self) -> Rational {
        Rational::ZERO
    }
    fn one(&self) -> Rational {
        Rational::integer(1)
    }
    fn is_zero(&self, value: &Rational) -> bool {
        value.numerator == 0
    }
    fn is_unit(&self, value: &Rational) -> bool {
        !self.is_zero(value)
    }

    fn add(&self, _: &mut Model, a: &Rational, b: &Rational) -> Result<Rational> {
        checked(a.add(*b))
    }

    fn subtract(&self, _: &mut Model, a: &Rational, b: &Rational) -> Result<Rational> {
        let negative = checked(b.multiply(Rational::integer(-1)))?;
        checked(a.add(negative))
    }

    fn multiply(&self, _: &mut Model, a: &Rational, b: &Rational) -> Result<Rational> {
        checked(a.multiply(*b))
    }

    fn inverse_unit(&self, _: &mut Model, value: &Rational) -> Result<Rational> {
        ratio(value.denominator, value.numerator)
    }
}
