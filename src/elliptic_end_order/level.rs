use std::collections::BTreeMap;
use std::rc::Rc;

use super::{realization::Realization, Certificate, Integer};
use crate::ec_poly::Scalar;
use crate::elliptic::{Curve, Point};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::legendre_cm::Element;
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};
use crate::wide_integer::Natural;

mod weil;

pub(crate) type Vector = [u64; 2];

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Matrix(pub(crate) [Vector; 2]);

pub(crate) enum Statement {
    Define(Label, Label, usize, Box<[ParsedExpr; 2]>),
    Data(Label),
    Invariants(Box<ParsedExpr>),
}

pub(crate) enum Expression {
    Vector(Label, [i128; 2]),
    Matrix(Label, [[i128; 2]; 2]),
    Frobenius(Label),
    BasisPairing(Label),
    Unary(Unary, Box<ParsedExpr>),
    Binary(Binary, Box<[ParsedExpr; 2]>),
}

pub(crate) enum Unary {
    OnOrder(Label),
    Coordinates(Label),
    Project(Label),
    Point,
}

pub(crate) enum Binary {
    Multiply,
    Apply,
    Weil(Label),
}

pub(crate) struct Level {
    pub(crate) modulus: u64,
    pub(crate) action: Rc<Realization>,
    pub(crate) curve: Rc<Curve>,
    certificate: Rc<Certificate>,
    generators: [Rc<Point>; 2],
    coordinates: BTreeMap<Point, Vector>,
    pub(crate) frobenius: Matrix,
    pub(crate) pairing: Scalar,
}

impl Level {
    pub(crate) fn new(
        model: &mut Model,
        action: Rc<Realization>,
        modulus: usize,
        generators: [Rc<Point>; 2],
    ) -> Result<Self> {
        let modulus = u64::try_from(modulus).map_err(|_| limit())?;
        let curve = model.ec_curve(&action.curve_name)?;
        if modulus < 2 || modulus.is_multiple_of(curve.ring.prime()) {
            return Err(Error::plain(
                "elliptic torsion level must be at least two and prime to the characteristic",
            ));
        }
        let count = usize::try_from(u128::from(modulus).pow(2)).map_err(|_| limit())?;
        model.allocate(
            count
                .checked_mul(10 + 2 * curve.ring.degree())
                .ok_or_else(limit)?,
        )?;
        model.spend(
            count
                .checked_mul(1 + count.ilog2() as usize * curve.ring.degree())
                .ok_or_else(limit)?,
        )?;

        for input in &generators {
            if curve.multiple(model, input, i128::from(modulus))? != Point::Infinity {
                return Err(Error::plain(
                    "elliptic level generator is not killed by the level",
                ));
            }
        }
        let coordinates = basis_table(model, &curve, modulus, &generators)?;
        let pairing = weil::pairing(model, &curve, modulus, &generators[0], &generators[1])?;
        exact_order(model, &curve, modulus, &pairing)?;
        let certificate = model.ec_end_order(&action.certificate_name().text)?;
        let mut level = Self {
            modulus,
            action,
            curve,
            certificate,
            generators,
            coordinates,
            frobenius: Matrix([[0; 2]; 2]),
            pairing,
        };
        let [a, b] = &level.generators;
        let [p, q] = [
            level.action.frobenius_point(model, a)?,
            level.action.frobenius_point(model, b)?,
        ];
        level.frobenius = level.certify_matrix(model, level.action.frobenius(), &p, &q)?;
        Ok(level)
    }

    pub(crate) fn point(&self, model: &mut Model, [a, b]: Vector) -> Result<Point> {
        let p = self
            .curve
            .multiple(model, &self.generators[0], i128::from(a))?;
        let q = self
            .curve
            .multiple(model, &self.generators[1], i128::from(b))?;
        self.curve.add(model, &p, &q)
    }

    pub(crate) fn coord(&self, model: &mut Model, point: &Point) -> Result<Vector> {
        model.spend(1 + self.coordinates.len().ilog2() as usize * self.curve.ring.degree())?;
        self.coordinates
            .get(point)
            .copied()
            .ok_or_else(|| Error::plain("elliptic point is outside the certified torsion level"))
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        source: &Self,
        vector: Vector,
    ) -> Result<Vector> {
        let point = source.point(model, vector)?;
        let quotient = i128::from(source.modulus / self.modulus);
        let image = self.curve.multiple(model, &point, quotient)?;
        self.coord(model, &image)
    }

    pub(crate) fn matrix(&self, model: &mut Model, element: Element) -> Result<Matrix> {
        let p = self.action.apply(model, element, &self.generators[0])?;
        let q = self.action.apply(model, element, &self.generators[1])?;
        self.certify_matrix(model, element, &p, &q)
    }

    fn certify_matrix(
        &self,
        model: &mut Model,
        element: Element,
        p: &Point,
        q: &Point,
    ) -> Result<Matrix> {
        let matrix = Matrix([self.coord(model, p)?, self.coord(model, q)?]);
        let norm = self.certificate.norm_value(model, element)?;
        let trace = self.certificate.trace_value(model, element)?;
        if matrix.det(self.modulus) != residue(model, self.modulus, norm)?
            || matrix.trace(self.modulus) != residue(model, self.modulus, trace)?
        {
            return Err(Error::plain(
                "elliptic level matrix failed the CM trace or norm check",
            ));
        }
        let actual = weil::pairing(model, &self.curve, self.modulus, p, q)?;
        let expected = self
            .curve
            .ring
            .power(model, &self.pairing, matrix.det(self.modulus))?;
        if actual != expected {
            return Err(Error::plain(
                "elliptic level matrix failed the Weil determinant relation",
            ));
        }
        Ok(matrix)
    }

    pub(crate) fn weil(&self, model: &mut Model, p: &Point, q: &Point) -> Result<Scalar> {
        let [a, b] = [self.coord(model, p)?, self.coord(model, q)?];
        let actual = weil::pairing(model, &self.curve, self.modulus, p, q)?;
        let power = Matrix([a, b]).det(self.modulus);
        if actual != self.curve.ring.power(model, &self.pairing, power)? {
            return Err(Error::plain(
                "elliptic Weil value disagreed with its coordinate determinant",
            ));
        }
        Ok(actual)
    }

    pub(crate) fn report(&self, name: &str) -> String {
        format!(
            "rvECLevelData[{name}] = curve={}, level={}, basis_pairing={:?}, frobenius_columns={:?}",
            self.action.curve_name, self.modulus, self.pairing, self.frobenius.0,
        )
    }
}

impl Matrix {
    pub(crate) fn trace(self, modulus: u64) -> u64 {
        add_mod(self.0[0][0], self.0[1][1], modulus)
    }

    pub(crate) fn det(self, modulus: u64) -> u64 {
        let [[a, c], [b, d]] = self.0;
        add_mod(
            mul_mod(a, d, modulus),
            modulus - mul_mod(b, c, modulus),
            modulus,
        )
    }

    pub(crate) fn apply(self, [x, y]: Vector, modulus: u64) -> Vector {
        let [[a, c], [b, d]] = self.0;
        [
            add_mod(mul_mod(a, x, modulus), mul_mod(b, y, modulus), modulus),
            add_mod(mul_mod(c, x, modulus), mul_mod(d, y, modulus), modulus),
        ]
    }

    pub(crate) fn product(self, other: Self, modulus: u64) -> Self {
        Self(other.0.map(|column| self.apply(column, modulus)))
    }
}

pub(crate) fn add_vectors(a: Vector, b: Vector, modulus: u64) -> Vector {
    std::array::from_fn(|i| add_mod(a[i], b[i], modulus))
}

fn add_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((u128::from(a) + u128::from(b)) % u128::from(modulus)) as u64
}

fn mul_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((u128::from(a) * u128::from(b)) % u128::from(modulus)) as u64
}

fn residue(model: &mut Model, modulus: u64, integer: Integer) -> Result<u64> {
    let (_, remainder) = integer
        .magnitude
        .div_rem(model, Natural::from_small(u128::from(modulus)))?;
    let remainder = remainder.small().expect("remainder smaller than modulus") as u64;
    Ok(if integer.negative && remainder != 0 {
        modulus - remainder
    } else {
        remainder
    })
}

fn basis_table(
    model: &mut Model,
    curve: &Curve,
    modulus: u64,
    generators: &[Rc<Point>; 2],
) -> Result<BTreeMap<Point, Vector>> {
    let mut table = BTreeMap::new();
    let mut row = Point::Infinity;
    for a in 0..modulus {
        let mut point = row.clone();
        for b in 0..modulus {
            if table.insert(point.clone(), [a, b]).is_some() {
                return Err(Error::plain(
                    "elliptic level generators are dependent modulo the level",
                ));
            }
            point = curve.add(model, &point, &generators[1])?;
        }
        row = curve.add(model, &row, &generators[0])?;
    }
    Ok(table)
}

fn exact_order(model: &mut Model, curve: &Curve, order: u64, value: &Scalar) -> Result<()> {
    let field = curve.field();
    if curve.ring.power(model, value, order)? != field.one() {
        return Err(Error::plain(
            "elliptic Weil value is not a root of unity of the level",
        ));
    }
    let mut rest = order;
    let mut prime = 2;
    while prime <= rest / prime {
        model.spend(1)?;
        if rest.is_multiple_of(prime) {
            reject_factor(model, curve, value, order, prime)?;
            while rest.is_multiple_of(prime) {
                rest /= prime;
            }
        }
        prime += 1;
    }
    if rest > 1 {
        reject_factor(model, curve, value, order, rest)?;
    }
    Ok(())
}

fn reject_factor(
    model: &mut Model,
    curve: &Curve,
    value: &Scalar,
    order: u64,
    prime: u64,
) -> Result<()> {
    if curve.ring.power(model, value, order / prime)? == curve.field().one() {
        return Err(Error::plain(
            "elliptic level basis Weil pairing does not have exact order",
        ));
    }
    Ok(())
}

fn limit() -> Error {
    Error::plain("elliptic level basis table exceeds the supported size")
}
