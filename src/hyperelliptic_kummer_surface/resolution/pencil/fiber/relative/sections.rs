use std::cmp::Ordering;
use std::rc::Rc;

use crate::ec_morphism::Rational;
use crate::ec_poly::{self as poly, Field, Polynomial};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::{field, Pencil, RelativeJacobian};

pub(crate) mod covering;
mod selmer;
mod square_classes;
mod torsion;

pub(crate) use selmer::Selmer;
pub(crate) use square_classes::Descent;

pub(crate) enum Expression {
    Literal(Label, Box<[Vec<Vec<i128>>; 4]>),
    KummerLiteral(Label, Box<[Vec<Vec<i128>>; 3]>),
    CoverPoint(Label, Box<[Vec<Vec<i128>>; 4]>),
    CoverMap(Box<ParsedExpr>),
    CoverLift(Label, Box<ParsedExpr>),
    ShaClass(Label, Label),
    Torsion(Label, usize),
    SelmerBasis(Label, usize),
    SelmerContains(Box<ParsedExpr>),
    Descent(Box<ParsedExpr>, bool),
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) enum Section {
    Infinity,
    Affine(Rational, Rational),
}

impl Section {
    pub(crate) fn literal(
        model: &mut Model,
        pencil: &Pencil,
        coordinates: [Vec<Vec<i128>>; 4],
    ) -> Result<Self> {
        let field = field(&pencil.curve);
        let [xn, xd, yn, yd] = read_polynomials(model, field, coordinates)?;
        if xd.is_empty() || yd.is_empty() {
            return Err(Error::plain(
                "K3 pencil section coordinate has zero denominator",
            ));
        }
        let x = Rational::new(model, field, xn, xd)?;
        let y = Rational::new(model, field, yn, yd)?;
        let relative = pencil.relative(model)?;
        if !on_curve(model, field, relative, &x, &y)? {
            return Err(Error::plain(
                "point does not lie on the generic K3 pencil Jacobian",
            ));
        }
        Ok(Self::Affine(x, y))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(self.cells())?;
        Ok(Value::K3PencilSection(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(self.cells().min(rhs.cells()))?;
        Ok(self.cmp(rhs))
    }

    fn cells(&self) -> usize {
        match self {
            Self::Infinity => 1,
            Self::Affine(x, y) => {
                5 + [&x.numerator, &x.denominator, &y.numerator, &y.denominator]
                    .into_iter()
                    .flatten()
                    .map(|coordinate| 3 + coordinate.len())
                    .sum::<usize>()
            }
        }
    }

    pub(crate) fn negate(&self, model: &mut Model, pencil: &Pencil) -> Result<Self> {
        let Self::Affine(x, y) = self else {
            return Ok(Self::Infinity);
        };
        let field = field(&pencil.curve);
        model.spend(self.cells())?;
        Ok(Self::Affine(
            x.clone(),
            y.scale(model, field, &field.signed(-1))?,
        ))
    }

    pub(crate) fn multiple(
        &self,
        model: &mut Model,
        pencil: &Pencil,
        mut multiple: u128,
    ) -> Result<Self> {
        let mut sum = Self::Infinity;
        let mut power = self.clone();
        while multiple != 0 {
            if multiple & 1 != 0 {
                sum = sum.add(model, pencil, &power)?;
            }
            multiple >>= 1;
            if multiple != 0 {
                power = power.add(model, pencil, &power)?;
            }
        }
        Ok(sum)
    }

    pub(crate) fn add(&self, model: &mut Model, pencil: &Pencil, rhs: &Self) -> Result<Self> {
        let (Self::Affine(x, y), Self::Affine(u, v)) = (self, rhs) else {
            let point = if matches!(self, Self::Infinity) {
                rhs
            } else {
                self
            };
            model.spend(point.cells())?;
            return Ok(point.clone());
        };
        model.spend(self.cells().min(rhs.cells()))?;
        let field = field(&pencil.curve);
        let relative = pencil.relative(model)?;
        let a2 = from_polynomial(field, relative.0.a2.clone());
        let slope = if x == u {
            if y != v || y.numerator.is_empty() {
                return Ok(Self::Infinity);
            }
            let square = x.multiply(model, field, x)?;
            let lead = square.scale(model, field, &field.signed(3))?;
            let term = a2.multiply(model, field, x)?;
            let term = term.scale(model, field, &field.signed(2))?;
            let numerator = sum(model, field, &lead, &term, false)?;
            let a4 = from_polynomial(field, relative.0.a4.clone());
            let numerator = sum(model, field, &numerator, &a4, false)?;
            let denominator = y.scale(model, field, &field.signed(2))?;
            quotient(model, field, &numerator, &denominator)?
        } else {
            let numerator = sum(model, field, v, y, true)?;
            let denominator = sum(model, field, u, x, true)?;
            quotient(model, field, &numerator, &denominator)?
        };
        let square = slope.multiply(model, field, &slope)?;
        let x3 = sum(model, field, &square, &a2, true)?;
        let x3 = sum(model, field, &x3, x, true)?;
        let x3 = sum(model, field, &x3, u, true)?;
        let difference = sum(model, field, x, &x3, true)?;
        let y3 = slope.multiply(model, field, &difference)?;
        let y3 = sum(model, field, &y3, y, true)?;
        Ok(Self::Affine(x3, y3))
    }
}

fn read_polynomials<const N: usize>(
    model: &mut Model,
    field: Field<'_>,
    rows: [Vec<Vec<i128>>; N],
) -> Result<[Polynomial; N]> {
    let mut result = std::array::from_fn(|_| Vec::new());
    for (polynomial, rows) in result.iter_mut().zip(rows) {
        field.check_polynomial(rows.len())?;
        let coefficients = rows
            .into_iter()
            .map(|row| field.ring.literal(model, row))
            .collect::<Result<_>>()?;
        *polynomial = poly::trim(coefficients);
    }
    Ok(result)
}

fn on_curve(
    model: &mut Model,
    field: Field<'_>,
    relative: &RelativeJacobian,
    x: &Rational,
    y: &Rational,
) -> Result<bool> {
    // The right-hand side in Horner form. Rational arithmetic keeps every
    // intermediate reduced, so equality is exact and independent of input.
    let mut curve = Rational::scalar(field, field.one());
    for coefficient in [&relative.0.a2, &relative.0.a4, &relative.0.a6] {
        let term = from_polynomial(field, coefficient.clone());
        curve = curve.multiply(model, field, x)?;
        curve = sum(model, field, &curve, &term, false)?;
    }
    Ok(y.multiply(model, field, y)? == curve)
}

fn from_polynomial(field: Field<'_>, polynomial: Polynomial) -> Rational {
    Rational {
        numerator: polynomial,
        denominator: vec![field.one()],
    }
}

fn sum(
    model: &mut Model,
    field: Field<'_>,
    left: &Rational,
    right: &Rational,
    subtract: bool,
) -> Result<Rational> {
    model.spend(left.denominator.len().min(right.denominator.len()) * field.cells())?;
    if left.denominator == right.denominator {
        let numerator = poly::sum(model, field, &left.numerator, &right.numerator, subtract)?;
        return Rational::new(model, field, numerator, left.denominator.clone());
    }
    let a = poly::multiply(model, field, &left.numerator, &right.denominator)?;
    let b = poly::multiply(model, field, &right.numerator, &left.denominator)?;
    let numerator = poly::sum(model, field, &a, &b, subtract)?;
    let denominator = poly::multiply(model, field, &left.denominator, &right.denominator)?;
    Rational::new(model, field, numerator, denominator)
}

fn quotient(
    model: &mut Model,
    field: Field<'_>,
    left: &Rational,
    right: &Rational,
) -> Result<Rational> {
    let numerator = poly::multiply(model, field, &left.numerator, &right.denominator)?;
    let denominator = poly::multiply(model, field, &left.denominator, &right.numerator)?;
    Rational::new(model, field, numerator, denominator)
}
