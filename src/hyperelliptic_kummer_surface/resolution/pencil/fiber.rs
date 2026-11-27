use crate::ec_poly::{self as poly, Field, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::rv_unram;
use crate::syntax::{Label, ParsedExpr};

use super::super::super::{field, k3};
use super::super::{quartic, QuarticTerm};
use super::{Parameter, Pencil};

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

type Quadratic = [Scalar; 3];

pub(crate) struct Fiber {
    a: Quadratic,
    b: Quadratic,
    c: Quadratic,
    binary: [Scalar; 5],
    a2: Scalar,
    a4: Scalar,
    a6: Scalar,
    c4: Scalar,
    discriminant: Scalar,
}

impl Pencil {
    pub(crate) fn equation(&self, model: &mut Model, parameter: &Parameter) -> Result<Fiber> {
        let field = field(&self.curve);
        let mut transverse = std::array::from_fn(|_| field.zero());
        for (form, value) in self.forms.iter().zip(parameter.coordinates()) {
            // The leading columns of the two forms are an identity block.
            let pivot = form
                .iter()
                .position(|entry| !field.is_zero(entry))
                .ok_or_else(|| Error::plain("K3 pencil lost a plane coordinate"))?;
            transverse[pivot] = value.clone();
        }
        let [i, j] = &self.nodes;
        let plane = [j.coordinates(), &transverse, i.coordinates()];
        let terms = quartic(model, &self.curve)?;
        let coefficients = restrict(model, field, &terms, plane)?;
        // Indexed by (Y exponent, Z exponent), with X exponent implicit.
        // The two double nodes exclude every monomial with X^3 or Z^3.
        for (y, row) in coefficients.iter().enumerate() {
            for (z, coefficient) in row.iter().enumerate().take(5 - y) {
                if (y + z < 2 || z > 2) && !field.is_zero(coefficient) {
                    return Err(Error::plain(
                        "K3 fiber plane quartic lacks the two prescribed double nodes",
                    ));
                }
            }
        }
        let a = std::array::from_fn(|i| coefficients[i][2].clone());
        let b = std::array::from_fn(|i| coefficients[i + 1][1].clone());
        let c = std::array::from_fn(|i| coefficients[i + 2][0].clone());
        if field.is_zero(&a[0]) {
            return Err(Error::plain(
                "K3 fiber plane quartic lost the secant coefficient",
            ));
        }
        Fiber::new(model, field, a, b, c)
    }
}

impl Fiber {
    fn new(
        model: &mut Model,
        field: Field<'_>,
        a: Quadratic,
        b: Quadratic,
        c: Quadratic,
    ) -> Result<Self> {
        let mut binary = std::array::from_fn(|_| field.zero());
        for i in 0..3 {
            for j in 0..3 {
                let square = field.multiply(model, &b[i], &b[j])?;
                let cross = field.multiply(model, &a[i], &c[j])?;
                let cross = field.ring.scale(model, &cross, 4)?;
                let term = field.subtract(model, &square, &cross)?;
                binary[i + j] = field.add(model, &binary[i + j], &term)?;
            }
        }
        let [a0, b0, c0, d0, e0] = &binary;
        let four_ae = field.multiply(model, a0, e0)?;
        let four_ae = field.ring.scale(model, &four_ae, 4)?;
        let bd = field.multiply(model, b0, d0)?;
        let a4 = field.subtract(model, &bd, &four_ae)?;
        let d2 = field.multiply(model, d0, d0)?;
        let ad2 = field.multiply(model, a0, &d2)?;
        let b2 = field.multiply(model, b0, b0)?;
        let b2e = field.multiply(model, &b2, e0)?;
        let sum = field.add(model, &ad2, &b2e)?;
        let four_ace = field.multiply(model, &four_ae, c0)?;
        let a6 = field.subtract(model, &sum, &four_ace)?;
        let a2 = c0.clone();
        let (c4, discriminant) = invariants(model, field, &a2, &a4, &a6)?;
        Ok(Self {
            a,
            b,
            c,
            binary,
            a2,
            a4,
            a6,
            c4,
            discriminant,
        })
    }

    fn status(&self) -> &'static str {
        if !poly::zero(&self.discriminant) {
            "smooth"
        } else if self.binary.iter().all(|entry| poly::zero(entry)) {
            "singular (zero binary quartic)"
        } else {
            "singular (repeated branch root)"
        }
    }

    fn require_smooth(&self) -> Result<()> {
        if poly::zero(&self.discriminant) {
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
            self.a2,
            self.a4,
            self.a6,
            self.discriminant,
            self.status(),
        )
    }

    pub(crate) fn j(&self, model: &mut Model, pencil: &Pencil) -> Result<Scalar> {
        self.require_smooth()?;
        let field = field(&pencil.curve);
        let cube = field.ring.power(model, &self.c4, 3)?;
        let inverse = field.inverse_unit(model, &self.discriminant)?;
        field.multiply(model, &cube, &inverse)
    }

    pub(crate) fn register(&self, model: &mut Model, pencil: &Pencil, name: Label) -> Result<()> {
        self.require_smooth()?;
        let field = field(&pencil.curve);
        if field.ring.prime() < 5 {
            return Err(Error::new(name.span, "rvk3fiberjac requires characteristic at least five for a named short Weierstrass curve"));
        }
        let third = field.inverse_unit(model, &field.signed(3))?;
        let shift = field.multiply(model, &self.a2, &third)?;
        let cross = field.multiply(model, &self.a2, &shift)?;
        let a = field.subtract(model, &self.a4, &cross)?;
        let linear = field.multiply(model, &self.a4, &shift)?;
        let b = field.subtract(model, &self.a6, &linear)?;
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

fn invariants(
    model: &mut Model,
    field: Field<'_>,
    a2: &Scalar,
    a4: &Scalar,
    a6: &Scalar,
) -> Result<(Scalar, Scalar)> {
    let b2 = field.ring.scale(model, a2, 4)?;
    let b4 = field.ring.scale(model, a4, 2)?;
    let b6 = field.ring.scale(model, a6, 4)?;
    let product = field.multiply(model, a2, &b6)?;
    let square = field.multiply(model, a4, a4)?;
    let b8 = field.subtract(model, &product, &square)?;
    let b2_square = field.multiply(model, &b2, &b2)?;
    let scaled = field.ring.scale(model, &b4, 24)?;
    let c4 = field.subtract(model, &b2_square, &scaled)?;
    let term = field.multiply(model, &b2_square, &b8)?;
    let mut discriminant = field.ring.negate(model, &term)?;
    let b4_cube = field.ring.power(model, &b4, 3)?;
    let term = field.ring.scale(model, &b4_cube, 8)?;
    discriminant = field.subtract(model, &discriminant, &term)?;
    let b6_square = field.multiply(model, &b6, &b6)?;
    let term = field.ring.scale(model, &b6_square, 27)?;
    discriminant = field.subtract(model, &discriminant, &term)?;
    let term = field.multiply(model, &b2, &b4)?;
    let term = field.multiply(model, &term, &b6)?;
    let term = field.ring.scale(model, &term, 9)?;
    discriminant = field.add(model, &discriminant, &term)?;
    Ok((c4, discriminant))
}

type PlaneCoefficients = [[Scalar; 5]; 5];

fn restrict(
    model: &mut Model,
    field: Field<'_>,
    terms: &[QuarticTerm],
    plane: [&[Scalar; 4]; 3],
) -> Result<PlaneCoefficients> {
    let zero = || std::array::from_fn(|_| std::array::from_fn(|_| field.zero()));
    let mut total = zero();
    for (slots, coefficient) in terms {
        let mut current = zero();
        current[0][0] = coefficient.clone();
        for (degree, &slot) in slots.iter().enumerate() {
            let mut next = zero();
            for (y, row) in current.iter().enumerate().take(degree + 1) {
                for (z, term) in row.iter().enumerate().take(degree + 1 - y) {
                    if field.is_zero(term) {
                        continue;
                    }
                    for (axis, vector) in plane.iter().enumerate() {
                        if field.is_zero(&vector[slot]) {
                            continue;
                        }
                        let coefficient =
                            &mut next[y + usize::from(axis == 1)][z + usize::from(axis == 2)];
                        let product = field.multiply(model, term, &vector[slot])?;
                        *coefficient = field.add(model, coefficient, &product)?;
                    }
                }
            }
            current = next;
        }
        for (y, row) in current.iter().enumerate() {
            for (z, term) in row.iter().enumerate().take(5 - y) {
                total[y][z] = field.add(model, &total[y][z], term)?;
            }
        }
    }
    Ok(total)
}

fn count_overflow() -> Error {
    Error::plain("K3 smooth fiber point count exceeds unsigned 128-bit range")
}
