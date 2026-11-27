use std::collections::BTreeMap;
use std::rc::Rc;

use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::elliptic::{self, Curve, DivisionPolynomials, Isogeny, Kernel, Point};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

#[derive(Clone, Copy)]
pub(crate) enum Report {
    Degree,
    Certify,
}

pub(crate) enum Expression {
    Identity(Label),
    Multiple(Label, i128),
    Chart(Label, Label, Box<ParsedExpr>),
    Unary(Unary, Box<ParsedExpr>),
    Binary(Binary, Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    FromVelu,
    Dual,
    Coefficient(Coordinate, usize),
    KernelCoefficient(usize),
}

impl Unary {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::FromVelu => "rvECMap",
            Self::Dual => "rvECDual",
            Self::Coefficient(Coordinate::XNumerator, _) => "rvECMapXNum",
            Self::Coefficient(Coordinate::XDenominator, _) => "rvECMapXDen",
            Self::Coefficient(Coordinate::YNumerator, _) => "rvECMapYNum",
            Self::Coefficient(Coordinate::YDenominator, _) => "rvECMapYDen",
            Self::KernelCoefficient(_) => "rvECMapKCoeff",
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Coordinate {
    XNumerator,
    XDenominator,
    YNumerator,
    YDenominator,
}

impl Coordinate {
    fn polynomial(self, map: &AffineMap) -> &[Scalar] {
        match self {
            Self::XNumerator => &map.x.numerator,
            Self::XDenominator => &map.x.denominator,
            Self::YNumerator => &map.y_multiplier.numerator,
            Self::YDenominator => &map.y_multiplier.denominator,
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Binary {
    Compose,
    Apply,
    Push,
}

impl Binary {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Compose => "rvECMapCompose",
            Self::Apply => "rvECMapApply",
            Self::Push => "rvECPush",
        }
    }
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Rational {
    pub(crate) numerator: Polynomial,
    pub(crate) denominator: Polynomial,
}

impl Rational {
    pub(crate) fn new(
        model: &mut Model,
        field: Field<'_>,
        numerator: Polynomial,
        denominator: Polynomial,
    ) -> Result<Self> {
        let numerator = poly::trim(numerator);
        let denominator = poly::trim(denominator);
        if denominator.is_empty() {
            return Err(Error::plain("rational map has a zero denominator"));
        }
        if numerator.is_empty() {
            return Ok(Self::scalar(field, field.zero()));
        }
        model.spend(numerator.len() + denominator.len())?;
        let divisor = poly::gcd(model, field, numerator.clone(), denominator.clone())?;
        let (numerator, denominator) = if divisor.len() > 1 {
            (
                poly::exact_divide(model, field, &numerator, &divisor)?,
                poly::exact_divide(model, field, &denominator, &divisor)?,
            )
        } else {
            (numerator, denominator)
        };
        let monic = denominator.last().expect("nonzero denominator");
        if *monic == field.one() {
            return Ok(Self {
                numerator,
                denominator,
            });
        }
        let unit = field.inverse_unit(model, monic)?;
        let numerator = poly::scale(model, field, &numerator, &unit)?;
        let denominator = poly::scale(model, field, &denominator, &unit)?;
        Ok(Self {
            numerator,
            denominator,
        })
    }

    pub(crate) fn x(field: Field<'_>) -> Self {
        Self {
            numerator: field.polynomial(&[0, 1]),
            denominator: vec![field.one()],
        }
    }

    pub(crate) fn scalar(field: Field<'_>, value: Scalar) -> Self {
        Self {
            numerator: if poly::zero(&value) {
                Vec::new()
            } else {
                vec![value]
            },
            denominator: vec![field.one()],
        }
    }

    pub(crate) fn scale(
        &self,
        model: &mut Model,
        field: Field<'_>,
        value: &Scalar,
    ) -> Result<Self> {
        if poly::zero(value) {
            return Ok(Self::scalar(field, field.zero()));
        }
        model.spend(self.denominator.len())?;
        let numerator = poly::scale(model, field, &self.numerator, value)?;
        Ok(Self {
            numerator,
            denominator: self.denominator.clone(),
        })
    }

    pub(crate) fn multiply(
        &self,
        model: &mut Model,
        field: Field<'_>,
        other: &Self,
    ) -> Result<Self> {
        let numerator = poly::multiply(model, field, &self.numerator, &other.numerator)?;
        let denominator = poly::multiply(model, field, &self.denominator, &other.denominator)?;
        Self::new(model, field, numerator, denominator)
    }

    fn compose(&self, model: &mut Model, field: Field<'_>, inner: &Self) -> Result<Self> {
        let degree = self.degree();
        field.check_polynomial(degree.saturating_mul(inner.degree()).saturating_add(1))?;
        let numerator = inner.homogenize(model, field, &self.numerator, degree)?;
        let denominator = inner.homogenize(model, field, &self.denominator, degree)?;
        Self::new(model, field, numerator, denominator)
    }

    fn homogenize(
        &self,
        model: &mut Model,
        field: Field<'_>,
        coefficient: &[Scalar],
        degree: usize,
    ) -> Result<Polynomial> {
        let mut value = Vec::new();
        let mut power = vec![field.one()];
        for k in (0..=degree).rev() {
            model.spend(1)?;
            value = poly::multiply(model, field, &value, &self.numerator)?;
            if let Some(c) = coefficient.get(k) {
                let term = poly::scale(model, field, &power, c)?;
                value = poly::sum(model, field, &value, &term, false)?;
            }
            if k != 0 {
                power = poly::multiply(model, field, &power, &self.denominator)?;
            }
        }
        Ok(value)
    }

    fn derivative(&self, model: &mut Model, field: Field<'_>) -> Result<Self> {
        let a = poly::derivative(model, field, &self.numerator)?;
        let b = poly::derivative(model, field, &self.denominator)?;
        let left = poly::multiply(model, field, &a, &self.denominator)?;
        let right = poly::multiply(model, field, &self.numerator, &b)?;
        let numerator = poly::sum(model, field, &left, &right, true)?;
        let denominator = poly::multiply(model, field, &self.denominator, &self.denominator)?;
        Self::new(model, field, numerator, denominator)
    }

    fn residue(
        &self,
        model: &mut Model,
        field: Field<'_>,
        modulus: &[Scalar],
    ) -> Result<Polynomial> {
        let inverse = poly::inverse_mod(model, field, &self.denominator, modulus)?;
        poly::multiply_mod(model, field, &self.numerator, &inverse, modulus)
    }

    fn finite_value(
        &self,
        model: &mut Model,
        field: Field<'_>,
        x: &[u64],
    ) -> Result<Option<Scalar>> {
        let denominator = poly::evaluate(model, field, &self.denominator, x)?;
        if poly::zero(&denominator) {
            return Ok(None);
        }
        let numerator = poly::evaluate(model, field, &self.numerator, x)?;
        let inverse = field.inverse_unit(model, &denominator)?;
        Ok(Some(field.multiply(model, &numerator, &inverse)?))
    }

    pub(crate) fn evaluate(
        &self,
        model: &mut Model,
        field: Field<'_>,
        x: &[u64],
    ) -> Result<Scalar> {
        self.finite_value(model, field, x)?
            .ok_or_else(|| Error::plain("rational coordinate has an affine pole"))
    }

    fn degree(&self) -> usize {
        self.numerator.len().max(self.denominator.len()) - 1
    }

    fn cells(&self) -> usize {
        6 + self
            .numerator
            .iter()
            .chain(&self.denominator)
            .map(|c| 3 + c.len())
            .sum::<usize>()
    }
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct AffineMap {
    pub(crate) x: Rational,
    pub(crate) y_multiplier: Rational,
}

#[derive(Clone)]
struct PrimeCertificate {
    kernel: Rc<Kernel>,
    differential: Scalar,
}

#[derive(Clone)]
pub(crate) struct Morphism {
    pub(crate) source: Rc<str>,
    pub(crate) target: Rc<str>,
    pub(crate) coordinates: Option<AffineMap>,
    prime: Option<PrimeCertificate>,
}

impl Morphism {
    pub(crate) fn zero(source: Rc<str>, target: Rc<str>) -> Self {
        Self {
            source,
            target,
            coordinates: None,
            prime: None,
        }
    }

    pub(crate) fn from_affine(
        model: &mut Model,
        source: Rc<str>,
        target: Rc<str>,
        domain: &Curve,
        range: &Curve,
        coordinates: AffineMap,
    ) -> Result<Self> {
        let map = Self {
            source,
            target,
            coordinates: Some(coordinates),
            prime: None,
        };
        map.certify(model, domain, range)?;
        Ok(map)
    }

    pub(crate) fn from_velu(model: &mut Model, domain: &Curve, isogeny: &Isogeny) -> Result<Self> {
        let field = domain.field();
        let coordinates = AffineMap {
            x: Rational::new(
                model,
                field,
                isogeny.numerator.clone(),
                isogeny.denominator.clone(),
            )?,
            y_multiplier: Rational::new(
                model,
                field,
                isogeny.derivative_numerator.clone(),
                isogeny.derivative_denominator.clone(),
            )?,
        };
        Ok(Self {
            source: isogeny.source.clone(),
            target: isogeny.target.clone(),
            coordinates: Some(coordinates),
            prime: Some(PrimeCertificate {
                kernel: isogeny.kernel.clone(),
                differential: field.one(),
            }),
        })
    }

    pub(crate) fn identity(model: &mut Model, name: Rc<str>, curve: &Curve) -> Result<Self> {
        Self::chart(
            model,
            name.clone(),
            name,
            curve,
            curve,
            &curve.field().one(),
        )
    }

    pub(crate) fn chart(
        model: &mut Model,
        source: Rc<str>,
        target: Rc<str>,
        domain: &Curve,
        range: &Curve,
        scale: &Scalar,
    ) -> Result<Self> {
        let (u2, u3) = chart_scalars(model, domain, range, scale)?;
        let field = domain.field();
        Ok(Self {
            source,
            target,
            coordinates: Some(AffineMap {
                x: Rational::x(field).scale(model, field, &u2)?,
                y_multiplier: Rational::scalar(field, u3),
            }),
            prime: None,
        })
    }

    pub(crate) fn degree(&self) -> usize {
        self.coordinates
            .as_ref()
            .map_or(0, |coordinates| coordinates.x.degree())
    }

    pub(crate) fn multiple(
        model: &mut Model,
        name: Rc<str>,
        curve: &Curve,
        signed: i128,
    ) -> Result<Self> {
        let n = u64::try_from(signed.unsigned_abs())
            .map_err(|_| Error::plain("elliptic multiplication degree is out of range"))?;
        let field = curve.field();
        if n == 0 {
            return Ok(Self::zero(name.clone(), name));
        }
        if n == 1 {
            let sign = field.signed(signed);
            return Self::chart(model, name.clone(), name, curve, curve, &sign);
        }

        let mut division = DivisionPolynomials::new(model, curve, None)?;
        let (x_num, x_den) = division.x_image(model, n)?;
        let (mut y_num, y_den) = division.y_multiplier(model, n)?;
        if signed < 0 {
            y_num = poly::scale(model, field, &y_num, &field.signed(-1))?;
        }
        let x = Rational::new(model, field, x_num, x_den)?;
        let y_multiplier = Rational::new(model, field, y_num, y_den)?;
        let map = Self::from_affine(
            model,
            name.clone(),
            name,
            curve,
            curve,
            AffineMap { x, y_multiplier },
        )?;
        if map.degree() as u128 != u128::from(n).pow(2) {
            return Err(Error::plain("multiplication map has the wrong degree"));
        }
        Ok(map)
    }

    pub(crate) fn compose(&self, model: &mut Model, domain: &Curve, inner: &Self) -> Result<Self> {
        if inner.target != self.source {
            return Err(Error::plain("elliptic morphisms have mismatched endpoints"));
        }
        let (Some(outer), Some(coords)) = (&self.coordinates, &inner.coordinates) else {
            return Ok(Self::zero(inner.source.clone(), self.target.clone()));
        };
        let degree = self
            .degree()
            .checked_mul(inner.degree())
            .ok_or_else(|| Error::plain("elliptic morphism composition degree is out of range"))?;
        let field = domain.field();
        field.check_polynomial(degree.saturating_add(1))?;
        let x = outer.x.compose(model, field, &coords.x)?;
        let y_multiplier = outer.y_multiplier.compose(model, field, &coords.x)?;
        let y_multiplier = y_multiplier.multiply(model, field, &coords.y_multiplier)?;
        let map = Self {
            source: inner.source.clone(),
            target: self.target.clone(),
            coordinates: Some(AffineMap { x, y_multiplier }),
            prime: None,
        };
        if map.degree() != degree {
            return Err(Error::plain(
                "elliptic morphism composition has the wrong degree",
            ));
        }
        Ok(map)
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        domain: &Curve,
        range: &Curve,
        point: &Point,
    ) -> Result<Point> {
        let (Some(coords), Point::Affine(x, y)) = (&self.coordinates, point) else {
            return Ok(Point::Infinity);
        };
        let field = domain.field();
        let Some(next_x) = coords.x.finite_value(model, field, x)? else {
            return Ok(Point::Infinity);
        };
        let multiplier = coords.y_multiplier.evaluate(model, field, x)?;
        let next_y = field.multiply(model, &multiplier, y)?;
        range.point(model, next_x, next_y)
    }

    pub(crate) fn certify(&self, model: &mut Model, source: &Curve, target: &Curve) -> Result<()> {
        same_field(source, target)?;
        let Some(coords) = &self.coordinates else {
            return Ok(());
        };
        let field = source.field();
        let a = &coords.x.numerator;
        let b = &coords.x.denominator;
        let c = &coords.y_multiplier.numerator;
        let d = &coords.y_multiplier.denominator;
        if a.len() <= b.len() {
            return Err(Error::plain(
                "an affine elliptic morphism must carry infinity to infinity",
            ));
        }
        // For f C² B³ = D²(A³ + a' A B² + b' B³), coprimality of each
        // coordinate and squarefreeness of f imply J = D²/B³ divides f.
        // Check the divisions before forming the much smaller chart identity.
        let b2 = poly::multiply(model, field, b, b)?;
        let b3 = poly::multiply(model, field, &b2, b)?;
        let d2 = poly::multiply(model, field, d, d)?;
        let (j, remainder) = poly::divide(model, field, d2, &b3)?;
        let chart_error =
            || Error::plain("elliptic morphism fails the exact cleared chart equation");
        if j.is_empty() || j.len() > source.polynomial().len() || !remainder.is_empty() {
            return Err(chart_error());
        }
        let (f, remainder) = poly::divide(model, field, source.polynomial().to_vec(), &j)?;
        if !remainder.is_empty() {
            return Err(chart_error());
        }
        let c2 = poly::multiply(model, field, c, c)?;
        let left = poly::multiply(model, field, &c2, &f)?;
        let a2 = poly::multiply(model, field, a, a)?;
        let a3 = poly::multiply(model, field, &a2, a)?;
        let linear = poly::multiply(model, field, a, &b2)?;
        let linear = poly::scale(model, field, &linear, &target.a)?;
        let constant = poly::scale(model, field, &b3, &target.b)?;
        let right = poly::sum(model, field, &a3, &linear, false)?;
        let right = poly::sum(model, field, &right, &constant, false)?;
        model.spend(left.len().min(right.len()) * field.cells())?;
        if left != right {
            return Err(chart_error());
        }
        Ok(())
    }

    pub(crate) fn dual(&self, model: &mut Model, domain: &Curve, range: &Curve) -> Result<Self> {
        let prime = self.prime.as_ref().ok_or_else(|| {
            Error::plain("rvECDual requires a morphism carrying a certified prime kernel")
        })?;
        let ell = prime.kernel.ell;
        let field = domain.field();
        let x = &self
            .coordinates
            .as_ref()
            .expect("prime isogeny is affine")
            .x;
        let psi = elliptic::division_polynomial(model, domain, ell)?;
        let other = poly::exact_divide(model, field, &psi, &prime.kernel.polynomial)?;
        let dual = image_kernel(model, x, range, ell, &other)?;
        let (velu, middle) = Isogeny::new(
            model,
            range,
            self.target.clone(),
            self.source.clone(),
            dual.clone(),
        )?;
        // If phi*omega is lambda omega, post-compose normalized Vélu
        // with the chart of scale lambda/ell. The dual multiplier is ell/lambda.
        let ell_scalar = field.signed(i128::from(ell));
        let ell_inverse = field.inverse_unit(model, &ell_scalar)?;
        let scale = field.multiply(model, &prime.differential, &ell_inverse)?;
        let (u2, u3) = chart_scalars(model, &middle, domain, &scale)?;
        let x = Rational::new(model, field, velu.numerator, velu.denominator)?
            .scale(model, field, &u2)?;
        let y_multiplier = Rational::new(
            model,
            field,
            velu.derivative_numerator,
            velu.derivative_denominator,
        )?
        .scale(model, field, &u3)?;
        let mut result = Self::from_affine(
            model,
            self.target.clone(),
            self.source.clone(),
            range,
            domain,
            AffineMap { x, y_multiplier },
        )?;
        let differential = field.inverse_unit(model, &scale)?;
        let coords = result.coordinates.as_ref().expect("dual is affine");
        let observed = coords.x.derivative(model, field)?;
        let expected = coords.y_multiplier.scale(model, field, &differential)?;
        model.spend(observed.cells().min(expected.cells()))?;
        if observed != expected {
            return Err(Error::plain(
                "dual isogeny has the wrong invariant differential",
            ));
        }
        result.prime = Some(PrimeCertificate {
            kernel: dual,
            differential,
        });
        Ok(result)
    }

    pub(crate) fn coefficient(
        &self,
        field: Field<'_>,
        coordinate: Coordinate,
        index: usize,
    ) -> Result<Scalar> {
        let affine = self.coordinates.as_ref().ok_or_else(|| {
            Error::plain("the constant infinity morphism has no affine rational coordinates")
        })?;
        Ok(coordinate
            .polynomial(affine)
            .get(index)
            .cloned()
            .unwrap_or_else(|| field.zero()))
    }

    pub(crate) fn kernel(&self) -> Option<&Kernel> {
        self.prime.as_ref().map(|tag| tag.kernel.as_ref())
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<std::cmp::Ordering> {
        let cells = |map: &Self| {
            map.coordinates
                .as_ref()
                .map_or(0, |c| c.x.cells() + c.y_multiplier.cells())
        };
        model.spend(1 + cells(self).min(cells(other)))?;
        Ok((&self.source, &self.target, &self.coordinates).cmp(&(
            &other.source,
            &other.target,
            &other.coordinates,
        )))
    }
}

fn chart_scalars(
    model: &mut Model,
    domain: &Curve,
    range: &Curve,
    scale: &Scalar,
) -> Result<(Scalar, Scalar)> {
    same_field(domain, range)?;
    if poly::zero(scale) {
        return Err(Error::plain("elliptic chart scale must be nonzero"));
    }
    let field = domain.field();
    let u2 = field.multiply(model, scale, scale)?;
    let u3 = field.multiply(model, &u2, scale)?;
    let u4 = field.multiply(model, &u2, &u2)?;
    let u6 = field.multiply(model, &u3, &u3)?;
    if field.multiply(model, &domain.a, &u4)? != range.a
        || field.multiply(model, &domain.b, &u6)? != range.b
    {
        return Err(Error::plain(
            "elliptic chart scale does not give the target equation",
        ));
    }
    Ok((u2, u3))
}

fn same_field(a: &Curve, b: &Curve) -> Result<()> {
    if a.ring_name != b.ring_name {
        return Err(Error::plain(
            "elliptic morphisms require the same named residue field",
        ));
    }
    Ok(())
}

pub(crate) fn persist(model: &mut Model, map: Morphism) -> Result<Value> {
    let cells = map
        .coordinates
        .as_ref()
        .map_or(0, |a| a.x.cells() + a.y_multiplier.cells());
    model.allocate(8 + cells)?;
    Ok(Value::ECMap(Rc::new(map)))
}

pub(crate) fn push(
    model: &mut Model,
    domain: &Curve,
    range: &Curve,
    map: &Isogeny,
    kernel: &Kernel,
) -> Result<Rc<Kernel>> {
    if map.kernel.ell == kernel.ell {
        return Err(Error::plain(
            "kernel push requires coprime prime isogeny and kernel degrees",
        ));
    }
    let x = Rational::new(
        model,
        domain.field(),
        map.numerator.clone(),
        map.denominator.clone(),
    )?;
    image_kernel(model, &x, range, kernel.ell, &kernel.polynomial)
}

fn image_kernel(
    model: &mut Model,
    x: &Rational,
    target: &Curve,
    ell: u64,
    modulus: &[Scalar],
) -> Result<Rc<Kernel>> {
    let field = target.field();
    let value = x.residue(model, field, modulus)?;
    let expected = (ell as usize - 1) / 2;
    let minimum = minimal_polynomial(model, field, &value, modulus, expected)?;
    if minimum.len() != expected + 1 {
        return Err(Error::plain(
            "image kernel minimal polynomial has the wrong degree",
        ));
    }
    Ok(Rc::new(Kernel::new(model, target, ell as usize, minimum)?))
}

struct Dependency {
    image: Polynomial,
    witness: Polynomial,
}

fn minimal_polynomial(
    model: &mut Model,
    field: Field<'_>,
    value: &[Scalar],
    modulus: &[Scalar],
    bound: usize,
) -> Result<Polynomial> {
    let mut rows = BTreeMap::<usize, Dependency>::new();
    let mut power = vec![field.one()];
    let mut cells = 0;
    // The first dependence of 1,x,x²,... gives the monic minimal
    // polynomial in the finite kernel algebra, even with repeated images.
    for degree in 0..=bound {
        model.spend(1 + degree + power.len())?;
        let mut image = power.clone();
        let mut witness = vec![field.zero(); degree];
        witness.push(field.one());
        while let Some(pivot) = image.len().checked_sub(1) {
            let Some(row) = rows.get(&pivot) else {
                break;
            };
            let scalar = image.last().expect("nonzero pivot");
            let change = poly::scale(model, field, &row.image, scalar)?;
            let witness_change = poly::scale(model, field, &row.witness, scalar)?;
            image = poly::sum(model, field, &image, &change, true)?;
            witness = poly::sum(model, field, &witness, &witness_change, true)?;
        }
        if image.is_empty() {
            return Ok(witness);
        }
        if degree == bound {
            break;
        }
        let inverse = field.inverse_unit(model, image.last().expect("nonzero pivot"))?;
        let pivot = image.len() - 1;
        let image = poly::scale(model, field, &image, &inverse)?;
        let witness = poly::scale(model, field, &witness, &inverse)?;
        cells += image.len() + witness.len() + 1;
        field.check_polynomial(cells)?;
        rows.insert(pivot, Dependency { image, witness });
        power = poly::multiply_mod(model, field, &power, value, modulus)?;
    }
    Err(Error::plain(
        "image kernel has no minimal polynomial within the predicted degree",
    ))
}
