use std::cmp::Ordering;
use std::rc::Rc;

use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::hyperelliptic_chart::Direction;
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::{Basis, Coordinate};
use crate::rv_unram::Ring;
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod miller;
pub(crate) mod point_halving;
pub(crate) mod quadratic_trace;

pub(crate) enum Command {
    Base(Label, Label, Label),
    Chart(Label, Label, Box<ParsedExpr>),
    Count(Label, bool),
    Order(Box<ParsedExpr>),
    Compose(Box<[ParsedExpr; 2]>),
    Kummer(Box<ParsedExpr>),
    HalfFiber(Label, Box<[ParsedExpr; 2]>),
}

pub(crate) enum Expression {
    Literal(Label, Vec<Vec<i128>>, Vec<Vec<i128>>),
    Point(Label, Box<[ParsedExpr; 2]>, Option<Vec<Vec<i128>>>),
    HalfAt(Label, String, Box<[ParsedExpr; 2]>),
    RichelotPoint(Label, Box<[ParsedExpr; 2]>),
    Transport(Label, Direction, Box<ParsedExpr>),
    QuadraticTrace(Label, Box<ParsedExpr>, Box<[Vec<Vec<i128>>; 4]>),
    Apply(Operation, Box<ParsedExpr>),
    Pair(Pairing, Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Copy)]
pub(crate) enum Pairing {
    Weil,
    Tate,
    MillerTate(u64),
    MillerWeil(u64),
}

impl Pairing {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Weil => "rvJacWeil2",
            Self::Tate => "rvJacTate2",
            Self::MillerTate(_) => "rvJacTate",
            Self::MillerWeil(_) => "rvJacWeil",
        }
    }

    pub(crate) fn evaluate(
        self,
        model: &mut Model,
        curve: &Curve,
        left: &Divisor,
        right: &Divisor,
    ) -> Result<Scalar> {
        match self {
            Self::MillerTate(n) => return miller::tate(model, curve, n, left, right),
            Self::MillerWeil(n) => return miller::weil(model, curve, n, left, right),
            _ => {}
        }
        // In reduced Mumford form, negation fixes u and sends v to −v.
        if !left.v.is_empty() {
            return Err(Error::plain(format!(
                "{} first argument must have order dividing 2",
                self.name()
            )));
        }
        if matches!(self, Self::Weil) && !right.v.is_empty() {
            return Err(Error::plain(
                "rvJacWeil2 second argument must have order dividing 2",
            ));
        }
        let field = field(curve);
        if left.is_zero() || right.is_zero() {
            return Ok(field.one());
        }
        model.spend((left.u.len() + right.u.len()) * field.cells())?;
        let common = poly::gcd(model, field, left.u.clone(), right.u.clone())?;
        match self {
            Self::Weil => {
                // Odd-degree branch subsets also contain the point at infinity.
                let intersection = (left.u.len() - 1) * (right.u.len() - 1) + common.len() - 1;
                Ok(if intersection.is_multiple_of(2) {
                    field.one()
                } else {
                    field.signed(-1)
                })
            }
            Self::Tate => left.tate_pairing(model, curve, right, common),
            Self::MillerTate(_) | Self::MillerWeil(_) => unreachable!(),
        }
    }
}

pub(crate) enum Operation {
    Scale(i128),
    Coefficient(bool, usize),
    Extend(Label),
    Descend(Label),
    Frobenius(Label, usize),
    Weil(Label),
    Divisible2,
    KummerSupport(usize),
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Divisor {
    pub(crate) u: Polynomial,
    pub(crate) v: Polynomial,
}

enum CantorFactor<'a> {
    Composition(&'a [Scalar]),
    Reduction {
        interpolation: &'a [Scalar],
        complement: &'a [Scalar],
        norm: &'a [Scalar],
    },
}

impl Divisor {
    pub(crate) fn zero(curve: &Curve) -> Self {
        Self {
            u: vec![curve.ring.scalar(1)],
            v: Vec::new(),
        }
    }

    fn is_zero(&self) -> bool {
        self.u.len() == 1
    }

    pub(crate) fn new(
        model: &mut Model,
        curve: &Curve,
        u: Polynomial,
        v: Polynomial,
    ) -> Result<Self> {
        let field = field(curve);
        let u = poly::trim(u);
        if u.last() != Some(&field.one()) {
            return Err(Error::plain("Mumford u must be nonzero and monic"));
        }
        let v = poly::remainder(model, field, poly::trim(v), &u)?;
        let norm = norm(model, curve, &v)?;
        if !poly::remainder(model, field, norm, &u)?.is_empty() {
            return Err(Error::plain("Mumford u must divide f−v²"));
        }
        Self { u, v }.reduce(model, curve)
    }

    pub(crate) fn norm_quotient(&self, model: &mut Model, curve: &Curve) -> Result<Polynomial> {
        let norm = norm(model, curve, &self.v)?;
        poly::exact_divide(model, field(curve), &norm, &self.u)
    }

    pub(crate) fn from_point(
        model: &mut Model,
        curve: &Curve,
        x: &[u64],
        y: Scalar,
    ) -> Result<Self> {
        let constant = curve.ring.negate(model, x)?;
        Self::new(model, curve, vec![constant, curve.ring.scalar(1)], vec![y])
    }

    fn reduce(self, model: &mut Model, curve: &Curve) -> Result<Self> {
        self.reduce_with(model, curve, |_, _| Ok(()))
    }

    fn reduce_with(
        mut self,
        model: &mut Model,
        curve: &Curve,
        mut on_factor: impl FnMut(&mut Model, CantorFactor<'_>) -> Result<()>,
    ) -> Result<Self> {
        let field = field(curve);
        while self.u.len() - 1 > curve.grade.genus {
            // The conjugate complementary effective divisor is obtained from
            // the zeros of y−v after removing the zeros prescribed by u.
            let norm = norm(model, curve, &self.v)?;
            let quotient = poly::exact_divide(model, field, &norm, &self.u)?;
            if quotient.is_empty() || quotient.len() >= self.u.len() {
                return Err(Error::plain(
                    "Cantor reduction did not decrease the divisor degree",
                ));
            }
            let inverse = field.inverse_unit(model, quotient.last().unwrap())?;
            let u = poly::scale(model, field, &quotient, &inverse)?;
            on_factor(
                model,
                CantorFactor::Reduction {
                    interpolation: &self.v,
                    complement: &u,
                    norm: &norm,
                },
            )?;
            let v = poly::scale(model, field, &self.v, &field.signed(-1))?;
            let v = poly::remainder(model, field, v, &u)?;
            self = Self { u, v };
        }
        Ok(self)
    }

    pub(crate) fn negate(&self, model: &mut Model, curve: &Curve) -> Result<Self> {
        let field = field(curve);
        model.spend(self.u.len() * field.cells())?;
        let v = poly::scale(model, field, &self.v, &field.signed(-1))?;
        Ok(Self {
            u: self.u.clone(),
            v,
        })
    }

    pub(crate) fn add(&self, model: &mut Model, curve: &Curve, rhs: &Self) -> Result<Self> {
        self.add_with(model, curve, rhs, |_, _| Ok(()))
    }

    fn add_with(
        &self,
        model: &mut Model,
        curve: &Curve,
        rhs: &Self,
        mut on_factor: impl FnMut(&mut Model, CantorFactor<'_>) -> Result<()>,
    ) -> Result<Self> {
        if self.is_zero() || rhs.is_zero() {
            let result = if self.is_zero() { rhs } else { self };
            model.spend(result.cells(curve))?;
            return Ok(result.clone());
        }
        let (semireduced, gcd) = self.compose_with_gcd(model, curve, rhs)?;
        on_factor(model, CantorFactor::Composition(&gcd))?;
        semireduced.reduce_with(model, curve, on_factor)
    }

    pub(crate) fn compose(&self, model: &mut Model, curve: &Curve, rhs: &Self) -> Result<Self> {
        self.compose_with_gcd(model, curve, rhs).map(|(sum, _)| sum)
    }

    fn compose_with_gcd(
        &self,
        model: &mut Model,
        curve: &Curve,
        rhs: &Self,
    ) -> Result<(Self, Polynomial)> {
        let field = field(curve);
        let v_sum = poly::sum(model, field, &self.v, &rhs.v, false)?;
        let (d0, e1, e2) = poly::bezout(model, field, &self.u, &rhs.u)?;
        let (d, c1, c2) = poly::bezout(model, field, &d0, &v_sum)?;
        let left = poly::multiply(model, field, &e1, &self.u)?;
        let left = poly::multiply(model, field, &left, &rhs.v)?;
        let right = poly::multiply(model, field, &e2, &rhs.u)?;
        let right = poly::multiply(model, field, &right, &self.v)?;
        let first = poly::sum(model, field, &left, &right, false)?;
        let first = poly::multiply(model, field, &c1, &first)?;
        let second = poly::multiply(model, field, &self.v, &rhs.v)?;
        let second = poly::sum(model, field, curve.polynomial(), &second, false)?;
        let second = poly::multiply(model, field, &c2, &second)?;
        let numerator = poly::sum(model, field, &first, &second, false)?;
        let v = poly::exact_divide(model, field, &numerator, &d)?;
        let product = poly::multiply(model, field, &self.u, &rhs.u)?;
        let square = poly::multiply(model, field, &d, &d)?;
        let u = poly::exact_divide(model, field, &product, &square)?;
        let v = poly::remainder(model, field, v, &u)?;
        Ok((Self { u, v }, d))
    }

    pub(crate) fn scale(&self, model: &mut Model, curve: &Curve, n: i128) -> Result<Self> {
        let mut result = Self::zero(curve);
        let mut n_abs = n.unsigned_abs();
        if self.is_zero() || n_abs == 0 {
            return Ok(result);
        }
        let mut current = if n < 0 {
            self.negate(model, curve)?
        } else {
            model.spend(self.cells(curve))?;
            self.clone()
        };
        while n_abs != 0 {
            if n_abs & 1 != 0 {
                result = result.add(model, curve, &current)?;
            }
            n_abs >>= 1;
            if n_abs != 0 {
                current = current.add(model, curve, &current)?;
            }
        }
        Ok(result)
    }

    pub(crate) fn order(&self, model: &mut Model, curve: &Curve) -> Result<i128> {
        if self.is_zero() {
            return Ok(1);
        }
        let mut order = cardinality(model, curve)?;
        let mut rest = order;
        let mut prime = 2;
        while prime <= rest / prime {
            model.spend(1)?;
            if rest % prime == 0 {
                order = self.strip(model, curve, order, prime)?;
                while rest % prime == 0 {
                    rest /= prime;
                }
            }
            prime += if prime == 2 { 1 } else { 2 };
        }
        if rest > 1 {
            order = self.strip(model, curve, order, rest)?;
        }
        Ok(order)
    }

    fn strip(
        &self,
        model: &mut Model,
        curve: &Curve,
        mut order: i128,
        prime: i128,
    ) -> Result<i128> {
        while order % prime == 0 && self.scale(model, curve, order / prime)?.is_zero() {
            order /= prime;
        }
        Ok(order)
    }

    pub(crate) fn extend(&self, model: &mut Model, base: &Curve, target: &Curve) -> Result<Self> {
        require_extension(model, base, target)?;
        let u = base
            .ring
            .embed_residue_coefficients(model, &target.ring, &self.u)?;
        let v = base
            .ring
            .embed_residue_coefficients(model, &target.ring, &self.v)?;
        Ok(Self { u, v })
    }

    pub(crate) fn descend(&self, model: &mut Model, source: &Curve, base: &Curve) -> Result<Self> {
        require_extension(model, base, source)?;
        let width = base.ring.degree();
        let mut powers = vec![base.ring.scalar(0); width];
        for (index, element) in powers.iter_mut().enumerate() {
            element[index] = 1;
        }
        let images = base
            .ring
            .embed_residue_coefficients(model, &source.ring, &powers)?;
        let prime = Prime::validated(base.ring.prime() as u32);
        let mut image = Basis::default();
        for (index, value) in images.iter().enumerate() {
            image.insert_with_witness(model, prime, coordinates(value), vec![(index, 1)])?;
        }
        let descend = |model: &mut Model, polynomial: &[Scalar]| -> Result<Polynomial> {
            polynomial
                .iter()
                .map(|coefficient| {
                    let (remainder, witness) =
                        image.reduce_with_witness(model, prime, coordinates(coefficient))?;
                    if !remainder.is_empty() {
                        return Err(Error::plain(
                            "Mumford coefficient lies outside the requested subfield",
                        ));
                    }
                    let mut value = base.ring.scalar(0);
                    for (index, coefficient) in witness {
                        value[index] = u64::from(coefficient);
                    }
                    Ok(value)
                })
                .collect()
        };
        Ok(Self {
            u: descend(model, &self.u)?,
            v: descend(model, &self.v)?,
        })
    }

    fn frobenius_coefficients(
        &self,
        model: &mut Model,
        target: &Curve,
        p_power: usize,
    ) -> Result<Self> {
        let transform = |model: &mut Model, coefficients: &[Scalar]| -> Result<Polynomial> {
            coefficients
                .iter()
                .map(|c| target.ring.frobenius(model, c, p_power))
                .collect()
        };
        let u = transform(model, &self.u)?;
        let v = transform(model, &self.v)?;
        Ok(Self { u, v })
    }

    pub(crate) fn frobenius(
        &self,
        model: &mut Model,
        base: &Curve,
        target: &Curve,
        power: usize,
    ) -> Result<Self> {
        require_extension(model, base, target)?;
        let base_degree = base.ring.degree();
        let extension_degree = target.ring.degree() / base_degree;
        self.frobenius_coefficients(model, target, (power % extension_degree) * base_degree)
    }

    pub(crate) fn weil(&self, model: &mut Model, base: &Curve, target: &Curve) -> Result<Self> {
        require_extension(model, base, target)?;
        let mut result = Self::zero(target);
        for &coefficient in base.weil(model)? {
            let image = result.frobenius_coefficients(model, target, base.ring.degree())?;
            let scaled = self.scale(model, target, coefficient)?;
            result = image.add(model, target, &scaled)?;
        }
        Ok(result)
    }

    fn tate_pairing(
        &self,
        model: &mut Model,
        curve: &Curve,
        rhs: &Self,
        common: Polynomial,
    ) -> Result<Scalar> {
        let field = field(curve);
        let h = &self.u;
        let value = if common.len() == 1 {
            model.spend(rhs.u.len() * field.cells())?;
            poly::resultant(model, field, rhs.u.clone(), h)?
        } else {
            let leading = common.last().expect("nonzero gcd of monic polynomials");
            let common = if leading == &field.one() {
                common
            } else {
                let inverse = field.inverse_unit(model, leading)?;
                poly::scale(model, field, &common, &inverse)?
            };
            let disjoint = poly::exact_divide(model, field, &rhs.u, &common)?;
            let ordinary = poly::resultant(model, field, disjoint, h)?;
            // At a common branch point y is a uniformizer and h/y² = 1/(f/h).
            // The inverse has the same square class. At infinity h has leading unit 1.
            let complement = poly::exact_divide(model, field, curve.polynomial(), h)?;
            let ramified = poly::resultant(model, field, common, &complement)?;
            field.multiply(model, &ordinary, &ramified)?
        };
        // The field norm accounts for (q−1)/(p−1), even if q does not fit in u64.
        let norm = curve.ring.conjugate_sum_or_product(model, &value, true)?;
        curve.ring.power(model, &norm, (curve.ring.prime() - 1) / 2)
    }

    pub(crate) fn coefficient(&self, curve: &Curve, v: bool, k: usize) -> Scalar {
        let polynomial = if v { &self.v } else { &self.u };
        polynomial
            .get(k)
            .cloned()
            .unwrap_or_else(|| curve.ring.scalar(0))
    }

    fn cells(&self, curve: &Curve) -> usize {
        2 + (self.u.len() + self.v.len()) * field(curve).cells()
    }

    pub(crate) fn persist(self, model: &mut Model, curve: &Curve) -> Result<Value> {
        model.allocate(self.cells(curve))?;
        Ok(Value::Jacobian(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        let degree = self.u[0].len();
        model.spend(
            2 + (self.u.len().min(rhs.u.len()) + self.v.len().min(rhs.v.len())) * (1 + degree),
        )?;
        Ok(self.cmp(rhs))
    }
}

fn coordinates(value: &[u64]) -> Vec<Coordinate> {
    value
        .iter()
        .enumerate()
        .filter_map(|(i, &v)| (v != 0).then_some((i, v as u32)))
        .collect()
}

fn field(curve: &Curve) -> Field<'_> {
    Field { ring: &curve.ring }
}

pub(crate) fn require_curve(curve: &Curve) -> Result<()> {
    if curve.ring.precision() != 1 {
        return Err(Error::plain(
            "Jacobians require a named residue field of precision one",
        ));
    }
    let f = curve.polynomial();
    if !f.len().is_multiple_of(2) || f.last() != Some(&curve.ring.scalar(1)) {
        return Err(Error::plain(
            "Jacobians require a smooth monic odd-degree hyperelliptic model",
        ));
    }
    Ok(())
}

pub(crate) fn coefficients_over(
    model: &mut Model,
    base: &Curve,
    target: &Ring,
) -> Result<Polynomial> {
    require_curve(base)?;
    base.ring
        .embed_residue_coefficients(model, target, base.polynomial())
}

fn require_extension(model: &mut Model, base: &Curve, target: &Curve) -> Result<()> {
    require_curve(target)?;
    let polynomial = coefficients_over(model, base, &target.ring)?;
    if polynomial != target.polynomial() {
        return Err(Error::plain(
            "target hyperelliptic equation is not the declared base extension",
        ));
    }
    Ok(())
}

fn norm(model: &mut Model, curve: &Curve, v: &[Scalar]) -> Result<Polynomial> {
    let field = field(curve);
    let square = poly::multiply(model, field, v, v)?;
    poly::sum(model, field, curve.polynomial(), &square, true)
}

pub(crate) fn cardinality(model: &mut Model, curve: &Curve) -> Result<i128> {
    let order = curve
        .weil(model)?
        .iter()
        .try_fold(0i128, |sum, &value| sum.checked_add(value));
    order.filter(|&order| order > 0).ok_or_else(|| {
        Error::plain("Jacobian cardinality exceeds the supported exact integer range")
    })
}

fn digits(ring: &Ring, q: usize, length: usize, mut index: usize) -> Polynomial {
    let p = ring.prime() as usize;
    (0..length)
        .map(|_| {
            let mut digit = index % q;
            index /= q;
            (0..ring.degree())
                .map(|_| {
                    let coefficient = digit % p;
                    digit /= p;
                    coefficient as u64
                })
                .collect()
        })
        .collect()
}

fn enumerate(
    model: &mut Model,
    curve: &Curve,
    mut visit: impl FnMut(Divisor) -> Result<()>,
) -> Result<()> {
    let field = field(curve);
    visit(Divisor::zero(curve))?;
    if curve.grade.genus == 0 {
        return Ok(());
    }
    let overflow = || Error::plain("Mumford pair enumeration exceeds the work limit");
    let q = usize::try_from(curve.ring.size_of_residue_field()?).map_err(|_| overflow())?;
    let mut width = 1usize;
    for degree in 1..=curve.grade.genus {
        width = width.checked_mul(q).ok_or_else(overflow)?;
        let lower_bound = width.checked_mul(width).ok_or_else(overflow)?;
        model.spend(lower_bound)?;
        for index in 0..width {
            model.spend(degree * field.cells())?;
            let mut u = digits(&curve.ring, q, degree, index);
            u.push(field.one());
            let f = poly::remainder(model, field, curve.polynomial().to_vec(), &u)?;
            for value in 0..width {
                model.spend(degree * field.cells())?;
                let v = poly::trim(digits(&curve.ring, q, degree, value));
                if poly::multiply_mod(model, field, &v, &v, &u)? == f {
                    model.spend(u.len() * field.cells())?;
                    visit(Divisor { u: u.clone(), v })?;
                }
            }
        }
    }
    Ok(())
}

pub(crate) fn count_pairs(model: &mut Model, curve: &Curve) -> Result<usize> {
    let mut count = 0;
    enumerate(model, curve, |_| {
        count += 1;
        Ok(())
    })?;
    Ok(count)
}

pub(crate) fn elements(model: &mut Model, curve: &Curve, limit: usize) -> Result<Vec<Value>> {
    let mut output = Vec::new();
    enumerate(model, curve, |divisor| {
        if output.len() == limit {
            return Err(Error::plain("Jacobian exceeds the finite element limit"));
        }
        output.push(divisor);
        Ok(())
    })?;
    let max_cells = 2 + (2 * curve.grade.genus + 1) * field(curve).cells();
    model.spend(
        output
            .len()
            .saturating_mul(1 + output.len().max(1).ilog2() as usize)
            .saturating_mul(max_cells),
    )?;
    output.sort_unstable();
    output
        .into_iter()
        .map(|divisor| divisor.persist(model, curve))
        .collect()
}
