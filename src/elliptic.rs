use std::cell::OnceCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::rv_unram::Ring;
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Command {
    Curve(Label, Label, Box<[ParsedExpr; 2]>),
    Isogeny(Label, Label, Label, Box<ParsedExpr>),
    PointCount(Label),
    Kernels(Label, usize),
    Certify(Label),
}

pub(crate) enum Expression {
    Curve(Label),
    CurveCoefficient(Label, bool),
    Isogeny(Label),
    Kernel(Label, usize, Vec<Vec<i128>>),
    Point(Label, Box<[ParsedExpr; 2]>),
    Infinity(Label),
    Apply(Box<[ParsedExpr; 2]>),
    Unary(Unary, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    PointCoordinate(bool),
    Multiple(i128),
    KernelCoefficient(usize),
    IsogenyKernel,
    IsogenyCoefficient(MapPolynomial, usize),
}

#[derive(Clone, Copy)]
pub(crate) enum MapPolynomial {
    XNumerator,
    XDenominator,
    DerivativeNumerator,
    DerivativeDenominator,
}

impl Unary {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::PointCoordinate(false) => "rvECX",
            Self::PointCoordinate(true) => "rvECY",
            Self::Multiple(_) => "rvECScale",
            Self::KernelCoefficient(_) => "rvECKCoeff",
            Self::IsogenyKernel => "rvECGetKernel",
            Self::IsogenyCoefficient(MapPolynomial::XNumerator, _) => "rvECXNum",
            Self::IsogenyCoefficient(MapPolynomial::XDenominator, _) => "rvECXDen",
            Self::IsogenyCoefficient(MapPolynomial::DerivativeNumerator, _) => "rvECDNum",
            Self::IsogenyCoefficient(MapPolynomial::DerivativeDenominator, _) => "rvECDDen",
        }
    }
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) enum Point {
    Infinity,
    Affine(Scalar, Scalar),
}

pub(crate) struct Curve {
    pub(crate) ring_name: Rc<str>,
    pub(crate) ring: Rc<Ring>,
    pub(crate) a: Scalar,
    pub(crate) b: Scalar,
    polynomial: Polynomial,
    trace: OnceCell<(u64, i64)>,
}

impl Curve {
    pub(crate) fn new(
        model: &mut Model,
        ring_name: Rc<str>,
        ring: Rc<Ring>,
        a: Scalar,
        b: Scalar,
    ) -> Result<Self> {
        if ring.precision() != 1 || ring.prime() < 5 {
            return Err(Error::plain(
                "short Weierstrass curves require a named residue field at precision one and p>=5",
            ));
        }
        let field = Field { ring: &ring };
        if a.len() != ring.degree() || b.len() != ring.degree() {
            return Err(Error::plain(
                "Weierstrass coefficients have the wrong residue field",
            ));
        }
        let a2 = field.multiply(model, &a, &a)?;
        let a3 = field.multiply(model, &a2, &a)?;
        let b2 = field.multiply(model, &b, &b)?;
        let t = field.multiply(model, &a3, &field.signed(4))?;
        let u = field.multiply(model, &b2, &field.signed(27))?;
        if poly::zero(&field.add(model, &t, &u)?) {
            return Err(Error::plain("short Weierstrass curve is singular"));
        }
        let polynomial = vec![b.clone(), a.clone(), field.zero(), field.one()];
        model.allocate(8 + 6 * field.cells())?;
        Ok(Self {
            ring_name,
            ring,
            a,
            b,
            polynomial,
            trace: OnceCell::new(),
        })
    }

    pub(crate) fn field(&self) -> Field<'_> {
        Field { ring: &self.ring }
    }

    pub(crate) fn polynomial(&self) -> &[Scalar] {
        &self.polynomial
    }

    pub(crate) fn trace(&self, model: &mut Model) -> Result<(u64, i64)> {
        if let Some(value) = self.trace.get() {
            return Ok(*value);
        }
        let q = self.ring.size_of_residue_field()?;
        let trace = self.ring.residue_quadratic_trace(model, |model, x| {
            poly::evaluate(model, self.field(), &self.polynomial, x)
        })?;
        Ok(*self.trace.get_or_init(|| (q, trace)))
    }

    pub(crate) fn has_kernel(&self, model: &mut Model, ell: u64) -> Result<bool> {
        // A rational cyclic kernel is an F_ell-line fixed by Frobenius.
        // The characteristic polynomial on E[ell] is T² - trace*T + q.
        let (q, trace) = self.trace(model)?;
        let modulus = i128::from(ell);
        let trace = i128::from(trace).rem_euclid(modulus);
        let discriminant = (trace * trace - 4 * i128::from(q)).rem_euclid(modulus) as u64;
        Ok(discriminant == 0 || integer_power_mod(model, discriminant, (ell - 1) / 2, ell)? == 1)
    }

    pub(crate) fn elements(&self, model: &mut Model, limit: usize) -> Result<Vec<Value>> {
        let q = self.ring.size_of_residue_field()?;
        if q > 4 * limit as u64 {
            return Err(Error::plain(
                "elliptic point type exceeds the finite element limit",
            ));
        }
        // The Hasse bound gives |E(F_q)| >= q/4 for q>=4.
        let field = self.field();
        field.check_polynomial(3 * q as usize)?;
        model.spend(q as usize)?;
        let coordinates = (0..q)
            .map(|mut code| {
                let mut value = field.zero();
                for digit in &mut value {
                    *digit = code % self.ring.prime();
                    code /= self.ring.prime();
                }
                value
            })
            .collect::<Vec<_>>();
        let mut roots = BTreeMap::<Scalar, Vec<&Scalar>>::new();
        for y in &coordinates {
            let square = field.multiply(model, y, y)?;
            model.spend(1 + roots.len().max(1).ilog2() as usize)?;
            roots.entry(square).or_default().push(y);
        }
        let mut points = vec![Point::Infinity];
        for x in &coordinates {
            let value = poly::evaluate(model, field, &self.polynomial, x)?;
            model.spend(1 + roots.len().max(1).ilog2() as usize)?;
            if let Some(ys) = roots.get(&value) {
                if points.len() + ys.len() > limit {
                    return Err(Error::plain(
                        "elliptic point type exceeds the finite element limit",
                    ));
                }
                points.extend(ys.iter().map(|y| Point::Affine(x.clone(), (*y).clone())));
            }
        }
        model.spend(points.len() * (1 + points.len().ilog2() as usize) * field.cells())?;
        points.sort_unstable();
        points
            .into_iter()
            .map(|point| persist_point(model, point))
            .collect()
    }

    pub(crate) fn point(&self, model: &mut Model, x: Scalar, y: Scalar) -> Result<Point> {
        let field = self.field();
        if field.multiply(model, &y, &y)? != poly::evaluate(model, field, &self.polynomial, &x)? {
            return Err(Error::plain(
                "point does not satisfy the named Weierstrass equation",
            ));
        }
        Ok(Point::Affine(x, y))
    }

    pub(crate) fn add(&self, model: &mut Model, a: &Point, b: &Point) -> Result<Point> {
        let (Point::Affine(x, y), Point::Affine(z, w)) = (a, b) else {
            return Ok(if *a == Point::Infinity {
                b.clone()
            } else {
                a.clone()
            });
        };
        let field = self.field();
        let slope = if x == z {
            if poly::zero(&field.add(model, y, w)?) {
                return Ok(Point::Infinity);
            }
            let x2 = field.multiply(model, x, x)?;
            let t = field.multiply(model, &x2, &field.signed(3))?;
            let numerator = field.add(model, &t, &self.a)?;
            let denom = field.add(model, y, y)?;
            divide_scalar(model, field, &numerator, &denom)?
        } else {
            let numerator = field.subtract(model, w, y)?;
            let denom = field.subtract(model, z, x)?;
            divide_scalar(model, field, &numerator, &denom)?
        };
        let square = field.multiply(model, &slope, &slope)?;
        let subtract_x = field.subtract(model, &square, x)?;
        let target_x = field.subtract(model, &subtract_x, z)?;
        let difference = field.subtract(model, x, &target_x)?;
        let product = field.multiply(model, &slope, &difference)?;
        let target_y = field.subtract(model, &product, y)?;
        Ok(Point::Affine(target_x, target_y))
    }

    pub(crate) fn multiple(&self, model: &mut Model, point: &Point, n: i128) -> Result<Point> {
        let mut power = point.clone();
        if n < 0 {
            if let Point::Affine(_, y) = &mut power {
                *y = self.ring.negate(model, y)?;
            }
        }
        let mut magnitude = n.unsigned_abs();
        let mut result = Point::Infinity;
        while magnitude != 0 {
            model.spend(1)?;
            if !magnitude.is_multiple_of(2) {
                result = self.add(model, &result, &power)?;
            }
            magnitude /= 2;
            if magnitude != 0 {
                power = self.add(model, &power, &power)?;
            }
        }
        Ok(result)
    }
}

pub(crate) fn persist_point(model: &mut Model, point: Point) -> Result<Value> {
    if let Point::Affine(x, y) = &point {
        model.allocate(6 + x.len() + y.len())?;
    }
    Ok(Value::ECPoint(Rc::new(point)))
}

pub(crate) fn point_coordinate(point: &Point, is_y: bool) -> Result<&Scalar> {
    match point {
        Point::Affine(x, y) => Ok(if is_y { y } else { x }),
        Point::Infinity => Err(Error::plain(
            "the point at infinity has no affine coordinate",
        )),
    }
}

fn divide_scalar(model: &mut Model, field: Field<'_>, a: &Scalar, b: &Scalar) -> Result<Scalar> {
    let inverse = field.inverse_unit(model, b)?;
    field.multiply(model, a, &inverse)
}

pub(crate) struct Kernel {
    pub(crate) ell: u64,
    pub(crate) polynomial: Rc<[Scalar]>,
}

impl Kernel {
    pub(crate) fn new(
        model: &mut Model,
        curve: &Curve,
        ell: usize,
        polynomial: Polynomial,
    ) -> Result<Self> {
        let ell = validate_ell(model, curve, ell)?;
        let size = (ell - 1) as usize / 2;
        let field = curve.field();
        if polynomial.len() != size + 1 || polynomial.last() != Some(&field.one()) {
            return Err(Error::plain(
                "kernel polynomial must be monic of degree (ell-1)/2",
            ));
        }
        field.check_polynomial(polynomial.len())?;
        let derivative = poly::derivative(model, field, &polynomial)?;
        if !poly::coprime(model, field, polynomial.clone(), derivative)? {
            return Err(Error::plain("kernel polynomial has repeated roots"));
        }
        let mut division = DivisionPolynomials::new(model, curve, Some(&polynomial))?;
        if !division.get(model, ell)?.is_empty() {
            return Err(Error::plain(
                "candidate kernel does not divide the ell-division polynomial",
            ));
        }
        let g = generator_mod_sign(model, ell)?;
        if !division.stable(model, g, &polynomial)? {
            return Err(Error::plain(
                "candidate is not stable under the cyclic generator modulo sign",
            ));
        }
        Self::from_certified(model, curve, ell, polynomial)
    }

    pub(crate) fn from_certified(
        model: &mut Model,
        curve: &Curve,
        ell: u64,
        polynomial: Polynomial,
    ) -> Result<Self> {
        model.allocate(4 + polynomial.len() * curve.field().cells())?;
        Ok(Self {
            ell,
            polynomial: polynomial.into(),
        })
    }
}

pub(crate) fn validate_ell(model: &mut Model, curve: &Curve, ell: usize) -> Result<u64> {
    if ell < 3 || ell.is_multiple_of(2) {
        return Err(Error::plain("cyclic isogeny degree must be an odd prime"));
    }
    let ell = u64::from(Prime::new(model, ell)?.get());
    if ell == curve.ring.prime() {
        return Err(Error::plain(
            "isogeny degree ell equals the residue characteristic",
        ));
    }
    Ok(ell)
}

fn generator_mod_sign(model: &mut Model, ell: u64) -> Result<u64> {
    let n = ell - 1;
    let mut remaining = n;
    let mut factors = Vec::new();
    let mut q = 2;
    while q <= remaining / q {
        model.spend(1)?;
        if remaining.is_multiple_of(q) {
            factors.push(q);
            while remaining.is_multiple_of(q) {
                model.spend(1)?;
                remaining /= q;
            }
        }
        q += 1;
    }
    if remaining != 1 {
        factors.push(remaining);
    }
    for g in 2..ell {
        model.spend(1)?;
        let mut okay = true;
        for &prime in &factors {
            if integer_power_mod(model, g, n / prime, ell)? == 1 {
                okay = false;
                break;
            }
        }
        if okay {
            // Either sign represents the same generator of F_ell^*/{±1}.
            return Ok(g.min(ell - g));
        }
    }
    Err(Error::plain("unable to find a cyclic generator modulo ell"))
}

fn integer_power_mod(model: &mut Model, mut a: u64, mut n: u64, m: u64) -> Result<u64> {
    let mut value = 1;
    while n != 0 {
        model.spend(1)?;
        if !n.is_multiple_of(2) {
            value = value * a % m;
        }
        n /= 2;
        a = a * a % m;
    }
    Ok(value)
}

pub(crate) fn division_polynomial(model: &mut Model, curve: &Curve, m: u64) -> Result<Polynomial> {
    if m == 0 || m.is_multiple_of(2) {
        return Err(Error::plain(
            "x-division polynomial requires a positive odd index",
        ));
    }
    let degree = (u128::from(m).pow(2) - 1) / 2;
    let coefficients = usize::try_from(degree + 1)
        .map_err(|_| Error::plain("division polynomial degree is out of range"))?;
    curve.field().check_polynomial(coefficients)?;
    let mut division = DivisionPolynomials::new(model, curve, None)?;
    let psi = division.get(model, m)?;
    model.spend(psi.len())?;
    Ok(psi.to_vec())
}

pub(crate) type Coefficients = Rc<[Scalar]>;

pub(crate) struct DivisionPolynomials<'a> {
    field: Field<'a>,
    modulus: Option<&'a [Scalar]>,
    z: Polynomial,
    z_squared: Polynomial,
    cache: BTreeMap<u64, Coefficients>,
    cells: usize,
}

impl<'a> DivisionPolynomials<'a> {
    pub(crate) fn new(
        model: &mut Model,
        curve: &'a Curve,
        modulus: Option<&'a [Scalar]>,
    ) -> Result<Self> {
        let field = curve.field();
        let z = poly::scale(model, field, curve.polynomial(), &field.signed(4))?;
        let mut division = Self {
            field,
            modulus,
            z: Vec::new(),
            z_squared: Vec::new(),
            cache: BTreeMap::new(),
            cells: 0,
        };
        division.z = division.reduce(model, z)?;
        division.z_squared = division.mul(model, &division.z, &division.z)?;
        let [e3, e4] = base_polynomials(model, curve)?;
        for (i, e) in [Vec::new(), vec![field.one()], vec![field.one()], e3, e4]
            .into_iter()
            .enumerate()
        {
            let e = division.reduce(model, e)?;
            division.insert(i as u64, e)?;
        }
        Ok(division)
    }

    fn insert(&mut self, i: u64, value: Polynomial) -> Result<Coefficients> {
        self.cells += 1 + value.len();
        self.field.check_polynomial(self.cells)?;
        let entry: Coefficients = value.into();
        self.cache.insert(i, entry.clone());
        Ok(entry)
    }

    fn reduce(&self, model: &mut Model, value: Polynomial) -> Result<Polynomial> {
        if let Some(modulus) = self.modulus {
            poly::remainder(model, self.field, value, modulus)
        } else {
            Ok(value)
        }
    }

    fn mul(&self, model: &mut Model, a: &[Scalar], b: &[Scalar]) -> Result<Polynomial> {
        let product = poly::multiply(model, self.field, a, b)?;
        self.reduce(model, product)
    }

    fn weighted(
        &self,
        model: &mut Model,
        a: &[Scalar],
        b: &[Scalar],
        power: usize,
    ) -> Result<Polynomial> {
        let mut value = self.mul(model, b, b)?;
        if power == 3 {
            value = self.mul(model, &value, b)?;
        }
        self.mul(model, a, &value)
    }

    pub(crate) fn get(&mut self, model: &mut Model, n: u64) -> Result<Coefficients> {
        model.spend(1)?;
        if let Some(value) = self.cache.get(&n) {
            return Ok(value.clone());
        }
        let m = n / 2;
        // e_n=psi_n for n odd and psi_n/(2y) for n even. This writes
        // every recurrence over base-field polynomials without y or division.
        let left = self.get(model, m + 2)?;
        let right = self.get(model, m + 1)?;
        let before = self.get(model, m - 1)?;
        let center = self.get(model, m)?;
        let value = if n.is_multiple_of(2) {
            let first = self.weighted(model, &left, &before, 2)?;
            let least = self.get(model, m - 2)?;
            let second = self.weighted(model, &least, &right, 2)?;
            let difference = poly::sum(model, self.field, &first, &second, true)?;
            self.mul(model, &center, &difference)?
        } else {
            let mut first = self.weighted(model, &left, &center, 3)?;
            let mut second = self.weighted(model, &before, &right, 3)?;
            if m.is_multiple_of(2) {
                first = self.mul(model, &first, &self.z_squared)?;
            } else {
                second = self.mul(model, &second, &self.z_squared)?;
            }
            poly::sum(model, self.field, &first, &second, true)?
        };
        self.insert(n, value)
    }

    pub(crate) fn x_image(
        &mut self,
        model: &mut Model,
        n: u64,
    ) -> Result<(Polynomial, Polynomial)> {
        if n == 0 {
            return Err(Error::plain(
                "the zero multiple has no affine x-rational function",
            ));
        }
        let next = n
            .checked_add(1)
            .ok_or_else(|| Error::plain("elliptic multiplication index is out of range"))?;
        if self.modulus.is_none() {
            let size = usize::try_from(u128::from(n).pow(2) + 1)
                .map_err(|_| Error::plain("elliptic multiplication degree is out of range"))?;
            self.field.check_polynomial(size)?;
        }
        let en = self.get(model, n)?;
        let lower = self.get(model, n - 1)?;
        let upper = self.get(model, next)?;
        let mut denominator = self.mul(model, &en, &en)?;
        let mut sub = self.mul(model, &lower, &upper)?;
        if n.is_multiple_of(2) {
            denominator = self.mul(model, &denominator, &self.z)?;
        } else {
            sub = self.mul(model, &sub, &self.z)?;
        }
        let x = self.field.polynomial(&[0, 1]);
        let x_times_denominator = self.mul(model, &x, &denominator)?;
        let numerator = poly::sum(model, self.field, &x_times_denominator, &sub, true)?;
        Ok((numerator, denominator))
    }

    fn stable(&mut self, model: &mut Model, g: u64, kernel: &[Scalar]) -> Result<bool> {
        if g == 1 {
            return Ok(true);
        }
        let (numerator, denominator) = self.x_image(model, g)?;
        if !poly::coprime(model, self.field, kernel.to_vec(), denominator.clone())? {
            return Ok(false);
        }
        // Homogeneous Horner evaluates denominator^degree K(numerator/denominator).
        let mut power = vec![self.field.one()];
        let mut value = power.clone();
        for coefficient in kernel[..kernel.len() - 1].iter().rev() {
            power = self.mul(model, &power, &denominator)?;
            let next = self.mul(model, &value, &numerator)?;
            let term = poly::scale(model, self.field, &power, coefficient)?;
            value = poly::sum(model, self.field, &next, &term, false)?;
        }
        Ok(value.is_empty())
    }
}

fn base_polynomials(model: &mut Model, curve: &Curve) -> Result<[Polynomial; 2]> {
    let f = curve.field();
    let a = &curve.a;
    let b = &curve.b;
    let a2 = f.multiply(model, a, a)?;
    let a3 = f.multiply(model, &a2, a)?;
    let ab = f.multiply(model, a, b)?;
    let b2 = f.multiply(model, b, b)?;
    let t = f.multiply(model, &b2, &f.signed(-16))?;
    let u = f.multiply(model, &a3, &f.signed(-2))?;
    let e3 = poly::trim(vec![
        f.ring.negate(model, &a2)?,
        f.multiply(model, b, &f.signed(12))?,
        f.multiply(model, a, &f.signed(6))?,
        f.zero(),
        f.signed(3),
    ]);
    let e4 = poly::trim(vec![
        f.add(model, &t, &u)?,
        f.multiply(model, &ab, &f.signed(-8))?,
        f.multiply(model, &a2, &f.signed(-10))?,
        f.multiply(model, b, &f.signed(40))?,
        f.multiply(model, a, &f.signed(10))?,
        f.zero(),
        f.signed(2),
    ]);
    Ok([e3, e4])
}

pub(crate) fn enumerate_kernels(
    model: &mut Model,
    curve: &Curve,
    ell: usize,
) -> Result<Vec<Polynomial>> {
    let ell = validate_ell(model, curve, ell)?;
    let degree = (ell - 1) as usize / 2;
    let q = curve.ring.size_of_residue_field()?;
    let candidates = q
        .checked_pow(degree as u32)
        .and_then(|n| usize::try_from(n).ok())
        .ok_or_else(|| Error::plain("cyclic kernel enumeration exceeds the work limit"))?;
    model.spend(candidates)?;
    let psi = division_polynomial(model, curve, ell)?;
    let field = curve.field();
    let generator = generator_mod_sign(model, ell)?;
    let p = curve.ring.prime();
    let mut result = Vec::new();
    let mut cells = 0;
    for mut index in 0..candidates as u64 {
        model.spend(degree * curve.ring.degree() + psi.len())?;
        let mut k = Vec::with_capacity(degree + 1);
        for _ in 0..degree {
            let mut code = index % q;
            index /= q;
            let mut coefficient = field.zero();
            for digit in &mut coefficient {
                *digit = code % p;
                code /= p;
            }
            k.push(coefficient);
        }
        k.push(field.one());
        if !poly::remainder(model, field, psi.clone(), &k)?.is_empty() {
            continue;
        }
        // For prime ell different from p, psi_ell is squarefree. Its
        // monic divisors have distinct ell-torsion x-coordinates.
        if generator != 1
            && !DivisionPolynomials::new(model, curve, Some(&k))?.stable(model, generator, &k)?
        {
            continue;
        }
        cells += k.len() + 1;
        field.check_polynomial(cells)?;
        result.push(k);
    }
    Ok(result)
}

pub(crate) struct Isogeny {
    pub(crate) source: Rc<str>,
    pub(crate) target: Rc<str>,
    pub(crate) kernel: Rc<Kernel>,
    pub(crate) numerator: Polynomial,
    pub(crate) denominator: Polynomial,
    pub(crate) derivative_numerator: Polynomial,
    pub(crate) derivative_denominator: Polynomial,
}

impl Isogeny {
    pub(crate) fn polynomial(&self, which: MapPolynomial) -> &[Scalar] {
        match which {
            MapPolynomial::XNumerator => &self.numerator,
            MapPolynomial::XDenominator => &self.denominator,
            MapPolynomial::DerivativeNumerator => &self.derivative_numerator,
            MapPolynomial::DerivativeDenominator => &self.derivative_denominator,
        }
    }

    pub(crate) fn new(
        model: &mut Model,
        curve: &Curve,
        source: Rc<str>,
        target: Rc<str>,
        kernel: Rc<Kernel>,
    ) -> Result<(Self, Curve)> {
        let field = curve.field();
        let k = &kernel.polynomial;
        let sums = poly::power_sums(model, field, k)?;
        let [power_zero, power_one, power_two, power_three] = &sums;
        let a = &curve.a;
        let b = &curve.b;

        let t0 = field.multiply(model, &field.signed(6), power_two)?;
        let t1 = field.multiply(model, &field.signed(2), a)?;
        let t2 = field.multiply(model, &t1, power_zero)?;
        let t = field.add(model, &t0, &t2)?;
        let w0 = field.multiply(model, &field.signed(10), power_three)?;
        let w1 = field.multiply(model, &field.signed(6), a)?;
        let w2 = field.multiply(model, &w1, power_one)?;
        let w3 = field.multiply(model, &field.signed(4), b)?;
        let w4 = field.multiply(model, &w3, power_zero)?;
        let w = field.add(model, &w0, &w2)?;
        let w = field.add(model, &w, &w4)?;
        let a_shift = field.multiply(model, &field.signed(5), &t)?;
        let b_shift = field.multiply(model, &field.signed(7), &w)?;
        let a_prime = field.subtract(model, a, &a_shift)?;
        let b_prime = field.subtract(model, b, &b_shift)?;
        let quotient = Curve::new(
            model,
            curve.ring_name.clone(),
            curve.ring.clone(),
            a_prime,
            b_prime,
        )?;

        let k_prime = poly::derivative(model, field, k)?;
        let u = vec![t1, field.zero(), field.signed(6)];
        let v = poly::scale(model, field, curve.polynomial(), &field.signed(4))?;
        let positive = weighted_numerator(model, field, k, &k_prime, &u, &sums)?;
        let negative = weighted_numerator(model, field, k, &k_prime, &v, &sums)?;
        let negative_prime = poly::derivative(model, field, &negative)?;
        let denominator = poly::multiply(model, field, k, k)?;
        let derivative_denominator = poly::multiply(model, field, k, &denominator)?;
        let x = field.polynomial(&[0, 1]);
        let linear = poly::multiply(model, field, &x, &denominator)?;
        let positive = poly::multiply(model, field, &positive, k)?;
        let negative_left = poly::multiply(model, field, &negative_prime, k)?;
        let negative_right = poly::multiply(model, field, &negative, &k_prime)?;
        let first = poly::sum(model, field, &linear, &positive, false)?;
        let first = poly::sum(model, field, &first, &negative_left, true)?;
        let numerator = poly::sum(model, field, &first, &negative_right, false)?;

        let slope = poly::derivative(model, field, &numerator)?;
        let preimage = poly::multiply(model, field, &slope, k)?;
        let denominator_slope = poly::multiply(model, field, &numerator, &k_prime)?;
        let denominator_slope = poly::scale(model, field, &denominator_slope, &field.signed(2))?;
        let derivative_numerator = poly::sum(model, field, &preimage, &denominator_slope, true)?;

        if numerator.len() != kernel.ell as usize + 1
            || numerator.last() != Some(&field.one())
            || !poly::coprime(model, field, numerator.clone(), k.to_vec())?
        {
            return Err(Error::plain(
                "normalized Vélu numerator failed degree or pole certification",
            ));
        }
        let map = Self {
            source,
            target,
            kernel,
            numerator,
            denominator,
            derivative_numerator,
            derivative_denominator,
        };
        map.certify(model, curve, &quotient)?;
        model.allocate(
            12 + field.cells()
                * (map.numerator.len()
                    + map.denominator.len()
                    + map.derivative_numerator.len()
                    + map.derivative_denominator.len()),
        )?;
        Ok((map, quotient))
    }

    pub(crate) fn certify(&self, model: &mut Model, curve: &Curve, target: &Curve) -> Result<()> {
        let field = curve.field();
        let derivative_square = poly::multiply(
            model,
            field,
            &self.derivative_numerator,
            &self.derivative_numerator,
        )?;
        let left = poly::multiply(model, field, curve.polynomial(), &derivative_square)?;
        let n2 = poly::multiply(model, field, &self.numerator, &self.numerator)?;
        let n3 = poly::multiply(model, field, &n2, &self.numerator)?;
        let k4 = poly::multiply(model, field, &self.denominator, &self.denominator)?;
        let k6 = poly::multiply(
            model,
            field,
            &self.derivative_denominator,
            &self.derivative_denominator,
        )?;
        let linear = poly::multiply(model, field, &self.numerator, &k4)?;
        let linear = poly::scale(model, field, &linear, &target.a)?;
        let constant = poly::scale(model, field, &k6, &target.b)?;
        let right = poly::sum(model, field, &n3, &linear, false)?;
        let right = poly::sum(model, field, &right, &constant, false)?;
        model.spend(left.len().min(right.len()) * field.cells())?;
        if left != right {
            return Err(Error::plain(
                "Vélu map failed the cleared Weierstrass equation",
            ));
        }
        Ok(())
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        source: &Curve,
        target: &Curve,
        input: &Point,
    ) -> Result<Point> {
        let Point::Affine(x, y) = input else {
            return Ok(Point::Infinity);
        };
        let field = source.field();
        let kernel = poly::evaluate(model, field, &self.kernel.polynomial, x)?;
        if poly::zero(&kernel) {
            return Ok(Point::Infinity);
        }
        let inverse = field.inverse_unit(model, &kernel)?;
        let inverse2 = field.multiply(model, &inverse, &inverse)?;
        let inverse3 = field.multiply(model, &inverse2, &inverse)?;
        let u = poly::evaluate(model, field, &self.numerator, x)?;
        let v = poly::evaluate(model, field, &self.derivative_numerator, x)?;
        let u = field.multiply(model, &u, &inverse2)?;
        let v = field.multiply(model, &v, &inverse3)?;
        let v = field.multiply(model, &v, y)?;
        target.point(model, u, v)
    }
}

// For a polynomial U, form a numerator of sum_r U(r)/(x-r):
// U(x) K'(x) - K(x) sum_r (U(x)-U(r))/(x-r). The second term
// uses only Newton sums and does not split the kernel polynomial.
fn weighted_numerator(
    model: &mut Model,
    field: Field<'_>,
    k: &[Scalar],
    derivative: &[Scalar],
    weights: &[Scalar],
    sums: &[Scalar; 4],
) -> Result<Polynomial> {
    let mut difference = vec![field.zero(); weights.len().saturating_sub(1)];
    for (i, coefficient) in weights.iter().enumerate().skip(1) {
        for (j, sum) in sums.iter().enumerate().take(i) {
            let term = field.multiply(model, coefficient, sum)?;
            let output = &mut difference[i - 1 - j];
            *output = field.add(model, output, &term)?;
        }
    }
    let left = poly::multiply(model, field, weights, derivative)?;
    let right = poly::multiply(model, field, k, &poly::trim(difference))?;
    poly::sum(model, field, &left, &right, true)
}
