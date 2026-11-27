use std::rc::Rc;

use super::clifford;
use super::rational_sparse::decimal;
use super::spinor;
use super::unimodular::{self, Vector, RANK};
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};
use crate::wide_integer::Natural;

// The doubled i128 coordinates give Gram entries below 2^260 and
// determinants below 2^521. This also covers squared saturation indices.
const LIMBS: usize = 18;
type Magnitude = Natural<LIMBS>;

pub(crate) enum Statement {
    Define(Label, Box<[ParsedExpr; 2]>),
    Report(Label),
    Is11(Label, Box<ParsedExpr>),
}

pub(crate) enum Expression {
    Blade(Label),
    CM(Label, Box<ParsedExpr>),
}

pub(crate) struct Period {
    marking: Rc<str>,
    u: Rc<Vector>,
    v: Rc<Vector>,
    discriminant: Magnitude,
    index: Magnitude,
    primitive_discriminant: Magnitude,
}

impl Period {
    pub(crate) fn new(
        model: &mut Model,
        marking: Rc<str>,
        u: Rc<Vector>,
        v: Rc<Vector>,
    ) -> Result<Self> {
        let uu = u.wide_pair(model, &u)?;
        if uu.negative || uu.is_zero() {
            return Err(Error::plain(
                "K3 period requires the first lattice vector to have positive square",
            ));
        }
        let vv = v.wide_pair(model, &v)?;
        let uv = u.wide_pair(model, &v)?.resize::<LIMBS>()?;
        let product: Wide<LIMBS> = uu.resize()?.multiply_wide(vv.resize()?)?;
        let determinant = product.add(uv.multiply_wide(uv)?.negate())?;
        if determinant.negative || determinant.is_zero() {
            return Err(Error::plain(
                "K3 period requires a positive Gram determinant",
            ));
        }
        let discriminant = determinant.magnitude;
        let index = saturation_index(model, &u, &v)?;
        let divisor = index.checked_product(index).ok_or_else(arithmetic)?;
        let primitive_discriminant =
            discriminant.divide_exact(model, divisor)?.ok_or_else(|| {
                Error::plain(
                    "K3 period Gram determinant is not divisible by the squared saturation index",
                )
            })?;
        model.allocate(4 + 3 * LIMBS)?;
        Ok(Self {
            marking,
            u,
            v,
            discriminant,
            index,
            primitive_discriminant,
        })
    }

    pub(crate) fn marking(&self) -> &str {
        &self.marking
    }

    pub(crate) fn report(&self, name: &str) -> String {
        let discriminant = decimal(self.discriminant);
        let index = decimal(self.index);
        let primitive = decimal(self.primitive_discriminant);
        let half_spin = spinor::DIMENSION;
        format!(
            "rvK3PeriodReport[{name}] = discriminant={discriminant}, saturation_index={index}, primitive_discriminant={primitive}, NS_rank=20, NS_signature=(1, 19), NS_determinant=-{primitive}, half_spin_dimension={half_spin}, half_spin_CM_rank={}, even_H1_dimension={}",
            half_spin / 2,
            1 << (RANK - 1),
        )
    }

    pub(crate) fn blade(&self, model: &mut Model) -> Result<clifford::Element> {
        clifford::Element::blade(model, &self.u, &self.v)
    }

    pub(crate) fn is_11(&self, model: &mut Model, vector: &Vector) -> Result<bool> {
        Ok(
            self.u.wide_pair(model, vector)?.is_zero()
                && self.v.wide_pair(model, vector)?.is_zero(),
        )
    }
}

fn saturation_index(model: &mut Model, u: &Vector, v: &Vector) -> Result<Magnitude> {
    let u = u.wide_basis_coordinates(model)?;
    let v = v.wide_basis_coordinates(model)?;
    model.spend(RANK * (RANK - 1) * LIMBS * 2)?;
    let mut index = Natural::default();
    for i in 0..RANK {
        for j in (i + 1)..RANK {
            let minor: unimodular::Integer = u[i]
                .multiply_wide(v[j])?
                .add(u[j].multiply_wide(v[i])?.negate())?;
            index = index.gcd(model, minor.magnitude)?;
            if index.small() == Some(1) {
                return Ok(Magnitude::from_small(1));
            }
        }
    }
    if index.is_zero() {
        return Err(Error::plain(
            "K3 period requires independent lattice vectors",
        ));
    }
    index.resize().ok_or_else(arithmetic)
}

fn arithmetic() -> Error {
    Error::plain("K3 period lattice arithmetic exceeds the supported exact range")
}
