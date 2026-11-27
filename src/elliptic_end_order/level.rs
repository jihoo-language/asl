use std::collections::BTreeMap;
use std::rc::Rc;

use super::{realization::Realization, Integer};
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
    Plain(Label, Label, Label, usize, Box<[ParsedExpr; 2]>),
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
    action: Option<Rc<Realization>>,
    pub(crate) base_name: Rc<str>,
    pub(crate) curve_name: Rc<str>,
    pub(crate) curve: Rc<Curve>,
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
        let certificate = model.ec_end_order(&action.certificate_name().text)?;
        Self::from_basis(
            model,
            certificate.curve_name().into(),
            Rc::clone(&action.curve_name),
            modulus,
            generators,
            Some(action),
        )
    }

    pub(crate) fn plain(
        model: &mut Model,
        base_name: Rc<str>,
        curve_name: Rc<str>,
        modulus: usize,
        generators: [Rc<Point>; 2],
    ) -> Result<Self> {
        let base = model.ec_curve(&base_name)?;
        let target = model.ec_curve(&curve_name)?;
        let coefficients = base.ring.embed_residue_coefficients(
            model,
            &target.ring,
            &[base.a.clone(), base.b.clone()],
        )?;
        if coefficients != [target.a.clone(), target.b.clone()] {
            return Err(Error::plain(
                "target elliptic equation is not the declared torsion base extension",
            ));
        }
        Self::from_basis(model, base_name, curve_name, modulus, generators, None)
    }

    fn from_basis(
        model: &mut Model,
        base_name: Rc<str>,
        curve_name: Rc<str>,
        modulus: usize,
        generators: [Rc<Point>; 2],
        action: Option<Rc<Realization>>,
    ) -> Result<Self> {
        let modulus = u64::try_from(modulus).map_err(|_| limit())?;
        let curve = model.ec_curve(&curve_name)?;
        let base = model.ec_curve(&base_name)?;
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
        let mut level = Self {
            modulus,
            action,
            base_name,
            curve_name,
            curve,
            generators,
            coordinates,
            frobenius: Matrix([[0; 2]; 2]),
            pairing,
        };
        let [a, b] = &level.generators;
        let power = base.ring.degree();
        let p = level.curve.frobenius_point(model, a, power)?;
        let q = level.curve.frobenius_point(model, b, power)?;
        let (order, trace) = base.trace(model)?;
        level.frobenius = level.certify_matrix(
            model,
            &p,
            &q,
            Integer::from_i128(i128::from(trace)),
            Integer::from_i128(i128::from(order)),
        )?;
        Ok(level)
    }

    pub(crate) fn action(&self) -> Result<&Realization> {
        self.action.as_deref().ok_or_else(|| {
            Error::plain("elliptic CM matrix requires a level attached to a witnessed CM action")
        })
    }

    pub(super) fn is_plain(&self) -> bool {
        self.action.is_none()
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
        let action = self.action()?;
        let p = action.apply(model, element, &self.generators[0])?;
        let q = action.apply(model, element, &self.generators[1])?;
        let certificate = model.ec_end_order(&action.certificate_name().text)?;
        let trace = certificate.trace_value(model, element)?;
        let norm = certificate.norm_value(model, element)?;
        self.certify_matrix(model, &p, &q, trace, norm)
    }

    fn certify_matrix(
        &self,
        model: &mut Model,
        p: &Point,
        q: &Point,
        trace: Integer,
        norm: Integer,
    ) -> Result<Matrix> {
        let matrix = Matrix([self.coord(model, p)?, self.coord(model, q)?]);
        if matrix.det(self.modulus) != residue(model, self.modulus, norm)?
            || matrix.trace(self.modulus) != residue(model, self.modulus, trace)?
        {
            return Err(Error::plain(
                "elliptic level matrix failed its trace or norm check",
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
            self.curve_name, self.modulus, self.pairing, self.frobenius.0,
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
