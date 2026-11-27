use std::cmp::Ordering;
use std::rc::Rc;

use crate::ec_morphism::Rational;
use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::finite_field_factor::factor_with_multiplicity;
use crate::hw_algebra::CoefficientRing;
use crate::hyperelliptic_kummer_surface::square_root;
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::square_classes::is_square;
use super::{
    field, from_polynomial, on_curve, quotient, read_polynomials, sum, Descent, Pencil, Section,
};

mod image;
pub(crate) use image::MordellImage;

pub(crate) enum Command {
    Define(Label, Label, Box<ParsedExpr>),
    Quadrics(Label),
    DefineImage(Label, Label, Vec<ParsedExpr>),
    Image(Label),
    Status(Label, Option<Label>),
}

pub(crate) struct Cover {
    pub(crate) pencil_name: Label,
    pub(crate) pencil: Rc<Pencil>,
    b: Rc<Descent>,
    representatives: [Polynomial; 3],
    product_root: Polynomial,
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Point([Polynomial; 4]);

impl Cover {
    pub(crate) fn new(
        model: &mut Model,
        pencil_name: Label,
        pencil: Rc<Pencil>,
        b: Rc<Descent>,
    ) -> Result<Self> {
        pencil.split_roots(model)?;
        let field = field(&pencil.curve);
        let mut representatives = std::array::from_fn(|_| vec![field.one()]);
        let nonsquare = if b.coordinates().any(|(odd, _)| odd) {
            constant_nonsquare(model, field)?
        } else {
            field.one()
        };
        let mut product = vec![field.one()];
        for (representative, (odd, factors)) in representatives.iter_mut().zip(b.coordinates()) {
            if odd {
                representative[0].clone_from(&nonsquare);
            }
            for factor in factors {
                *representative = poly::multiply(model, field, representative, factor)?;
            }
            product = poly::multiply(model, field, &product, representative)?;
        }
        let product_root = polynomial_square_root(model, field, &product)?;
        model.allocate(
            10 + pencil_name.text.len()
                + field.cells()
                    * (product_root.len() + representatives.iter().map(Vec::len).sum::<usize>()),
        )?;
        Ok(Self {
            pencil_name,
            pencil,
            b,
            representatives,
            product_root,
        })
    }

    pub(crate) fn quadrics(&self, model: &mut Model, name: &str) -> Result<String> {
        let field = field(&self.pencil.curve);
        let [e0, e1, e2] = self.pencil.split_roots(model)?;
        let [b0, b1, b2] = &self.representatives;
        let first = [
            b0.clone(),
            poly::scale(model, field, b1, &field.signed(-1))?,
            Vec::new(),
            poly::sum(model, field, e0, e1, true)?,
        ];
        let second = [
            b0.clone(),
            Vec::new(),
            poly::scale(model, field, b2, &field.signed(-1))?,
            poly::sum(model, field, e0, e2, true)?,
        ];
        Ok(format!(
            "rvK3PencilCoverQuadrics[{name}] = diagonal (U0²,U1²,U2²,W²)={:?}, product square root={:?}",
            [first, second], self.product_root,
        ))
    }

    pub(crate) fn point(
        &self,
        model: &mut Model,
        coordinates: [Vec<Vec<i128>>; 4],
    ) -> Result<Point> {
        let field = field(&self.pencil.curve);
        let coordinates = read_polynomials(model, field, coordinates)?;
        self.checked_point(model, field, coordinates)
    }

    fn checked_point(
        &self,
        model: &mut Model,
        field: Field<'_>,
        mut coordinates: [Polynomial; 4],
    ) -> Result<Point> {
        let mut nonzero = coordinates.iter().filter(|p| !p.is_empty());
        let first = nonzero
            .next()
            .ok_or_else(|| Error::plain("K3 pencil covering projective point must be nonzero"))?;
        let mut common = first.clone();
        for coordinate in nonzero {
            if common.len() <= 1 {
                break;
            }
            common = poly::gcd(model, field, common, coordinate.clone())?;
        }
        if common.len() > 1 {
            for c in &mut coordinates {
                *c = poly::exact_divide(model, field, c, &common)?;
            }
        }
        let leading = coordinates
            .iter()
            .find_map(|p| p.last())
            .expect("nonzero coordinate");
        if *leading != field.one() {
            let inverse = field.inverse_unit(model, leading)?;
            for c in &mut coordinates {
                *c = poly::scale(model, field, c, &inverse)?;
            }
        }
        let [u0, u1, u2, w] = &coordinates;
        let [e0, e1, e2] = self.pencil.split_roots(model)?;
        let a = poly::multiply(model, field, u0, u0)?;
        let a = poly::multiply(model, field, &self.representatives[0], &a)?;
        let w2 = poly::multiply(model, field, w, w)?;
        for (i, (u, e)) in [(u1, e1), (u2, e2)].into_iter().enumerate() {
            let delta = poly::sum(model, field, e, e0, true)?;
            let right = poly::multiply(model, field, &delta, &w2)?;
            let square = poly::multiply(model, field, u, u)?;
            let weighted = poly::multiply(model, field, &self.representatives[i + 1], &square)?;
            let left = poly::sum(model, field, &a, &weighted, true)?;
            if left != right {
                return Err(Error::plain(
                    "point does not lie on both K3 pencil covering quadrics",
                ));
            }
        }
        Ok(Point(coordinates))
    }

    pub(crate) fn map(&self, model: &mut Model, point: &Point) -> Result<Section> {
        let [a, b, c, w] = &point.0;
        let field = field(&self.pencil.curve);
        let section = if w.is_empty() {
            Section::Infinity
        } else {
            let [e0, ..] = self.pencil.split_roots(model)?;
            let relative = self.pencil.relative(model)?;
            let w2 = poly::multiply(model, field, w, w)?;
            let w3 = poly::multiply(model, field, &w2, w)?;
            let a2 = poly::multiply(model, field, a, a)?;
            let term = poly::multiply(model, field, &self.representatives[0], &a2)?;
            let root_term = poly::multiply(model, field, e0, &w2)?;
            let numerator = poly::sum(model, field, &term, &root_term, false)?;
            let x = Rational::new(model, field, numerator, w2)?;
            let ab = poly::multiply(model, field, a, b)?;
            let abc = poly::multiply(model, field, &ab, c)?;
            let numerator = poly::multiply(model, field, &self.product_root, &abc)?;
            let y = Rational::new(model, field, numerator, w3)?;
            if !on_curve(model, field, relative, &x, &y)? {
                return Err(Error::plain(
                    "K3 pencil covering map failed the Jacobian equation",
                ));
            }
            Section::Affine(x, y)
        };
        if section.descent(model, &self.pencil)? != *self.b {
            return Err(Error::plain(
                "K3 pencil covering map failed the prescribed Kummer class",
            ));
        }
        Ok(section)
    }

    pub(crate) fn lift(&self, model: &mut Model, section: &Section) -> Result<Point> {
        if section.descent(model, &self.pencil)? != *self.b {
            return Err(Error::plain(
                "K3 pencil section is outside the covering image",
            ));
        }
        let field = field(&self.pencil.curve);
        let Section::Affine(x, y) = section else {
            let one = vec![field.one()];
            return self.checked_point(model, field, [one.clone(), one.clone(), one, Vec::new()]);
        };
        let mut roots = std::array::from_fn::<_, 3, _>(|_| Rational::scalar(field, field.zero()));
        let mut w = vec![field.one()];
        for ((root, e), b) in roots
            .iter_mut()
            .zip(self.pencil.split_roots(model)?)
            .zip(&self.representatives)
        {
            let diff = sum(model, field, x, &from_polynomial(field, e.clone()), true)?;
            let square = quotient(model, field, &diff, &from_polynomial(field, b.clone()))?;
            *root = rational_square_root(model, field, &square)?;
            w = poly::multiply(model, field, &w, &root.denominator)?;
        }
        let mut coordinates = std::array::from_fn(|_| Vec::new());
        for (u, root) in coordinates[..3].iter_mut().zip(roots) {
            let scale = poly::exact_divide(model, field, &w, &root.denominator)?;
            *u = poly::multiply(model, field, &root.numerator, &scale)?;
        }
        coordinates[3] = w;
        let point = self.checked_point(model, field, coordinates)?;
        let image = self.map(model, &point)?;
        if image == *section {
            return Ok(point);
        }
        // Negating one root changes only the sign of y. At a branch point
        // y=0 and the first checked lift already maps to the right point.
        if image != Section::Affine(x.clone(), y.scale(model, field, &field.signed(-1))?) {
            return Err(Error::plain(
                "K3 pencil covering lift failed to recover the section",
            ));
        }
        let mut coordinates = point.0;
        coordinates[0] = poly::scale(model, field, &coordinates[0], &field.signed(-1))?;
        let point = self.checked_point(model, field, coordinates)?;
        if self.map(model, &point)? != *section {
            return Err(Error::plain(
                "K3 pencil covering lift failed its sign check",
            ));
        }
        Ok(point)
    }
}

impl Point {
    pub(crate) fn coordinates(&self) -> &[Polynomial; 4] {
        &self.0
    }

    fn cells(&self) -> usize {
        5 + self
            .0
            .iter()
            .flatten()
            .map(|row| row.len() + 3)
            .sum::<usize>()
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(self.cells())?;
        Ok(Value::K3PencilCoverPoint(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(self.cells() + rhs.cells())?;
        Ok(self.cmp(rhs))
    }
}

impl Pencil {
    fn split_roots(&self, model: &mut Model) -> Result<[&Polynomial; 3]> {
        let [e0, e1, e2] = self.torsion_roots(model)? else {
            return Err(Error::plain(
                "K3 pencil two-coverings require full split two-torsion",
            ));
        };
        Ok([e0, e1, e2])
    }

    pub(crate) fn kummer_literal(
        &self,
        model: &mut Model,
        representatives: [Vec<Vec<i128>>; 3],
    ) -> Result<Descent> {
        self.split_roots(model)?;
        let field = field(&self.curve);
        let representatives = read_polynomials(model, field, representatives)?;
        Descent::from_representatives(model, field, representatives)
    }
}

fn constant_nonsquare(model: &mut Model, field: Field<'_>) -> Result<Scalar> {
    for index in 1..field.ring.count() {
        model.spend(field.ring.degree())?;
        let candidate = field.ring.coordinate_row(index);
        if !is_square(model, field, &candidate)? {
            return Ok(candidate);
        }
    }
    Err(Error::plain(
        "K3 pencil covering field has no nonsquare constant",
    ))
}

fn polynomial_square_root(model: &mut Model, field: Field<'_>, p: &[Scalar]) -> Result<Polynomial> {
    let Some(leading) = p.last() else {
        return Ok(Vec::new());
    };
    let constant = if *leading == field.one() {
        field.one()
    } else {
        square_root(model, field, leading)?.ok_or_else(not_square)?
    };
    let mut root = vec![constant];
    for (_, order, factor) in factor_with_multiplicity(model, field, p)? {
        if !order.is_multiple_of(2) {
            return Err(not_square());
        }
        for _ in 0..order / 2 {
            root = poly::multiply(model, field, &root, &factor)?;
        }
    }
    if poly::multiply(model, field, &root, &root)? != p {
        return Err(not_square());
    }
    Ok(root)
}

fn rational_square_root(
    model: &mut Model,
    field: Field<'_>,
    square: &Rational,
) -> Result<Rational> {
    let numerator = polynomial_square_root(model, field, &square.numerator)?;
    let denominator = polynomial_square_root(model, field, &square.denominator)?;
    Rational::new(model, field, numerator, denominator)
}

fn not_square() -> Error {
    Error::plain("K3 pencil covering square root does not exist in the rational function field")
}
