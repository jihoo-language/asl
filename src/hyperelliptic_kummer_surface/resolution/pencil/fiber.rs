use crate::ec_poly::{self as poly, Field, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::rv_unram;
use crate::syntax::{Label, ParsedExpr};

use super::super::super::{field, k3};
use super::super::quartic;
use super::{Parameter, Pencil};

pub(super) mod relative;

pub(crate) struct Command(
    pub(crate) Label,
    pub(crate) Query,
    pub(crate) Box<ParsedExpr>,
);

pub(crate) enum Query {
    Quartic,
    Jacobian,
    Points(u128),
    Register(Label),
}

type Quadratic<T> = [T; 3];

trait FiberRing: CoefficientRing {
    fn scale(&self, model: &mut Model, value: &Self::Element, factor: u64)
        -> Result<Self::Element>;
}

impl FiberRing for Field<'_> {
    fn scale(&self, model: &mut Model, value: &Scalar, factor: u64) -> Result<Scalar> {
        self.ring.scale(model, value, factor)
    }
}

pub(crate) struct Fiber {
    a: Quadratic<Scalar>,
    b: Quadratic<Scalar>,
    c: Quadratic<Scalar>,
    binary: [Scalar; 5],
    jacobian: Jacobian<Scalar>,
}

struct Jacobian<T> {
    a2: T,
    a4: T,
    a6: T,
    c4: T,
    c6: T,
    discriminant: T,
}

impl Pencil {
    fn transverse_axes(&self, field: Field<'_>) -> Result<[usize; 2]> {
        let [first, second] = self.forms.each_ref().map(|form| {
            form.iter()
                .position(|entry| !field.is_zero(entry))
                .ok_or_else(|| Error::plain("K3 pencil lost a plane coordinate"))
        });
        Ok([first?, second?])
    }

    pub(crate) fn equation(&self, model: &mut Model, parameter: &Parameter) -> Result<Fiber> {
        let field = field(&self.curve);
        let mut transverse = std::array::from_fn(|_| field.zero());
        for (pivot, value) in self
            .transverse_axes(field)?
            .into_iter()
            .zip(parameter.coordinates())
        {
            // The leading columns of the two forms are an identity block.
            transverse[pivot] = value.clone();
        }
        let [i, j] = &self.nodes;
        let plane = [j.coordinates(), &transverse, i.coordinates()];
        let terms = quartic(model, &self.curve)?;
        let [a, b, c] = plane_quadratics(model, &field, &terms, plane)?;
        let binary = binary_quartic(model, &field, &a, &b, &c)?;
        let jacobian = Jacobian::new(model, &field, &binary)?;
        Ok(Fiber {
            a,
            b,
            c,
            binary,
            jacobian,
        })
    }
}

fn plane_quadratics<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    terms: &[([usize; 4], R::Element)],
    plane: [&[R::Element; 4]; 3],
) -> Result<[Quadratic<R::Element>; 3]> {
    let coefficients = restrict(model, ring, terms, plane)?;
    // Indexed by (Y exponent, Z exponent), with X exponent implicit.
    // The two double nodes exclude every monomial with X^3 or Z^3.
    for (y, row) in coefficients.iter().enumerate() {
        for (z, coefficient) in row.iter().enumerate().take(5 - y) {
            if (y + z < 2 || z > 2) && !ring.is_zero(coefficient) {
                return Err(Error::plain(
                    "K3 fiber plane quartic lacks the two prescribed double nodes",
                ));
            }
        }
    }
    let a = std::array::from_fn(|i| coefficients[i][2].clone());
    let b = std::array::from_fn(|i| coefficients[i + 1][1].clone());
    let c = std::array::from_fn(|i| coefficients[i + 2][0].clone());
    if ring.is_zero(&a[0]) {
        return Err(Error::plain(
            "K3 fiber plane quartic lost the secant coefficient",
        ));
    }
    Ok([a, b, c])
}

fn binary_quartic<R: FiberRing>(
    model: &mut Model,
    ring: &R,
    a: &Quadratic<R::Element>,
    b: &Quadratic<R::Element>,
    c: &Quadratic<R::Element>,
) -> Result<[R::Element; 5]> {
    let mut binary = std::array::from_fn(|_| ring.zero());
    for i in 0..3 {
        for j in 0..3 {
            let square = ring.multiply(model, &b[i], &b[j])?;
            let cross = ring.multiply(model, &a[i], &c[j])?;
            let cross = ring.scale(model, &cross, 4)?;
            let term = ring.subtract(model, &square, &cross)?;
            binary[i + j] = ring.add(model, &binary[i + j], &term)?;
        }
    }
    Ok(binary)
}

impl<T: Clone> Jacobian<T> {
    fn new<R: FiberRing<Element = T>>(
        model: &mut Model,
        ring: &R,
        binary: &[T; 5],
    ) -> Result<Self> {
        let [a, b, c, d, e] = binary;
        let four_ae = ring.multiply(model, a, e)?;
        let four_ae = ring.scale(model, &four_ae, 4)?;
        let bd = ring.multiply(model, b, d)?;
        let a4 = ring.subtract(model, &bd, &four_ae)?;
        let d2 = ring.multiply(model, d, d)?;
        let ad2 = ring.multiply(model, a, &d2)?;
        let b2 = ring.multiply(model, b, b)?;
        let b2e = ring.multiply(model, &b2, e)?;
        let sum = ring.add(model, &ad2, &b2e)?;
        let four_ace = ring.multiply(model, &four_ae, c)?;
        let a6 = ring.subtract(model, &sum, &four_ace)?;
        let a2 = c.clone();
        let (c4, c6, discriminant) = invariants(model, ring, &a2, &a4, &a6)?;
        Ok(Self {
            a2,
            a4,
            a6,
            c4,
            c6,
            discriminant,
        })
    }
}

impl Fiber {
    fn status(&self) -> &'static str {
        if !poly::zero(&self.jacobian.discriminant) {
            "smooth"
        } else if self.binary.iter().all(|entry| poly::zero(entry)) {
            "singular (zero binary quartic)"
        } else {
            "singular (repeated branch root)"
        }
    }

    fn require_smooth(&self) -> Result<()> {
        if poly::zero(&self.jacobian.discriminant) {
            return Err(Error::plain("K3 pencil fiber is singular"));
        }
        Ok(())
    }

    pub(crate) fn quartic_report(&self, name: &str) -> String {
        format!(
            "rvK3FiberQuartic[{name}] = A2={:?}, B2={:?}, C2={:?}, F4={:?}, status={}",
            self.a,
            self.b,
            self.c,
            self.binary,
            self.status(),
        )
    }

    pub(crate) fn jacobian_report(&self, name: &str) -> String {
        format!(
            "rvK3FiberJacobian[{name}] = a2={:?}, a4={:?}, a6={:?}, discriminant={:?}, status={}",
            self.jacobian.a2,
            self.jacobian.a4,
            self.jacobian.a6,
            self.jacobian.discriminant,
            self.status(),
        )
    }

    pub(crate) fn j(&self, model: &mut Model, pencil: &Pencil) -> Result<Scalar> {
        self.require_smooth()?;
        let field = field(&pencil.curve);
        let cube = field.ring.power(model, &self.jacobian.c4, 3)?;
        let inverse = field.inverse_unit(model, &self.jacobian.discriminant)?;
        field.multiply(model, &cube, &inverse)
    }

    pub(crate) fn register(&self, model: &mut Model, pencil: &Pencil, name: Label) -> Result<()> {
        self.require_smooth()?;
        let field = field(&pencil.curve);
        if field.ring.prime() < 5 {
            return Err(Error::new(name.span, "rvk3fiberjac requires characteristic at least five for a named short Weierstrass curve"));
        }
        let third = field.inverse_unit(model, &field.signed(3))?;
        let shift = field.multiply(model, &self.jacobian.a2, &third)?;
        let cross = field.multiply(model, &self.jacobian.a2, &shift)?;
        let a = field.subtract(model, &self.jacobian.a4, &cross)?;
        let linear = field.multiply(model, &self.jacobian.a4, &shift)?;
        let b = field.subtract(model, &self.jacobian.a6, &linear)?;
        let cubic = field.ring.power(model, &shift, 3)?;
        let cubic = field.ring.scale(model, &cubic, 2)?;
        let b = field.add(model, &b, &cubic)?;
        let ring = Label {
            text: pencil.curve.grade.ring.to_string(),
            span: name.span,
        };
        model.define_ec_curve(name, ring, a, b)
    }

    pub(crate) fn points(&self, model: &mut Model, pencil: &Pencil, degree: u128) -> Result<u128> {
        self.require_smooth()?;
        let field = field(&pencil.curve);
        let q = field.ring.size_of_residue_field()?;
        // Even with maximal cancellation in the Hasse bound, a field size
        // >=2^129 cannot give a count fitting u128. Keep smaller cases wide.
        let degree = usize::try_from(degree)
            .ok()
            .filter(|&r| (1..=128).contains(&r))
            .ok_or_else(count_overflow)?;
        let mut order = Wide::<8>::from_i128(1);
        for _ in 0..degree {
            order = order
                .multiply(i128::from(q))
                .map_err(|_| count_overflow())?;
            if order.magnitude.bits() > 129 {
                return Err(count_overflow());
            }
        }
        let coefficients: Vec<_> = self.binary.iter().rev().cloned().collect();
        let mut trace = i128::from(field.ring.residue_quadratic_trace(model, |model, x| {
            poly::evaluate(model, field, &coefficients, x)
        })?);
        if !field.is_zero(&self.binary[0]) {
            // The point at infinity of the base contributes chi(F4(1,0)).
            let character = field.ring.power(model, &self.binary[0], (q - 1) / 2)?;
            trace -= match rv_unram::require_scalar(&character)? {
                1 => 1,
                c if c == field.ring.prime() - 1 => -1,
                _ => {
                    return Err(Error::plain(
                        "invalid quadratic character of K3 fiber quartic",
                    ))
                }
            };
        }
        let traces =
            k3::h1_power_sums::<8>(model, &[1, -trace, i128::from(q)], degree, count_overflow)?;
        let count = order
            .add(Wide::from_i128(1))
            .and_then(|total| total.add(traces[degree].negate()))
            .map_err(|_| count_overflow())?;
        if count.negative && !count.is_zero() {
            return Err(Error::plain("K3 fiber point count became negative"));
        }
        count.magnitude.small().ok_or_else(count_overflow)
    }
}

fn invariants<R: FiberRing>(
    model: &mut Model,
    ring: &R,
    a2: &R::Element,
    a4: &R::Element,
    a6: &R::Element,
) -> Result<(R::Element, R::Element, R::Element)> {
    let b2 = ring.scale(model, a2, 4)?;
    let b4 = ring.scale(model, a4, 2)?;
    let b6 = ring.scale(model, a6, 4)?;
    let product = ring.multiply(model, a2, &b6)?;
    let square = ring.multiply(model, a4, a4)?;
    let b8 = ring.subtract(model, &product, &square)?;
    let b2_square = ring.multiply(model, &b2, &b2)?;
    let scaled = ring.scale(model, &b4, 24)?;
    let c4 = ring.subtract(model, &b2_square, &scaled)?;
    let term = ring.multiply(model, &b2_square, &b8)?;
    let mut discriminant = ring.subtract(model, &ring.zero(), &term)?;
    let b4_cube = cube(model, ring, &b4)?;
    let term = ring.scale(model, &b4_cube, 8)?;
    discriminant = ring.subtract(model, &discriminant, &term)?;
    let b6_square = ring.multiply(model, &b6, &b6)?;
    let term = ring.scale(model, &b6_square, 27)?;
    discriminant = ring.subtract(model, &discriminant, &term)?;
    let term = ring.multiply(model, &b2, &b4)?;
    let term = ring.multiply(model, &term, &b6)?;
    let term = ring.scale(model, &term, 9)?;
    discriminant = ring.add(model, &discriminant, &term)?;
    let b2_cube = ring.multiply(model, &b2_square, &b2)?;
    let b2b4 = ring.multiply(model, &b2, &b4)?;
    let middle = ring.scale(model, &b2b4, 36)?;
    let last = ring.scale(model, &b6, 216)?;
    let c6 = ring.subtract(model, &middle, &b2_cube)?;
    let c6 = ring.subtract(model, &c6, &last)?;
    Ok((c4, c6, discriminant))
}

fn cube<R: CoefficientRing>(model: &mut Model, ring: &R, value: &R::Element) -> Result<R::Element> {
    let square = ring.multiply(model, value, value)?;
    ring.multiply(model, &square, value)
}

type PlaneCoefficients<T> = [[T; 5]; 5];

fn restrict<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    terms: &[([usize; 4], R::Element)],
    plane: [&[R::Element; 4]; 3],
) -> Result<PlaneCoefficients<R::Element>> {
    let zero = || std::array::from_fn(|_| std::array::from_fn(|_| ring.zero()));
    let mut total = zero();
    for (slots, coefficient) in terms {
        let mut current = zero();
        current[0][0] = coefficient.clone();
        for (degree, &slot) in slots.iter().enumerate() {
            let mut next = zero();
            for (y, row) in current.iter().enumerate().take(degree + 1) {
                for (z, term) in row.iter().enumerate().take(degree + 1 - y) {
                    if ring.is_zero(term) {
                        continue;
                    }
                    for (axis, vector) in plane.iter().enumerate() {
                        if ring.is_zero(&vector[slot]) {
                            continue;
                        }
                        let coefficient =
                            &mut next[y + usize::from(axis == 1)][z + usize::from(axis == 2)];
                        let product = ring.multiply(model, term, &vector[slot])?;
                        *coefficient = ring.add(model, coefficient, &product)?;
                    }
                }
            }
            current = next;
        }
        for (y, row) in current.iter().enumerate() {
            for (z, term) in row.iter().enumerate().take(5 - y) {
                total[y][z] = ring.add(model, &total[y][z], term)?;
            }
        }
    }
    Ok(total)
}

fn count_overflow() -> Error {
    Error::plain("K3 smooth fiber point count exceeds unsigned 128-bit range")
}
