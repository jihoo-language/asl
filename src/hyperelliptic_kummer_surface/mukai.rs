use std::cmp::Ordering;
use std::rc::Rc;

use super::polarized::{self, Divisor};
use super::unimodular::{self, Integer, Vector as K3Vector};
use crate::error::{Error, Result};
use crate::jacobian_level::Basis;
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};
use crate::wide_integer::Natural;

pub(crate) mod abelian;
pub(crate) mod moduli;

pub(crate) enum Expression {
    FromDivisor(Box<ParsedExpr>),
    ToDivisor(Box<ParsedExpr>),
    Vector(i128, i128, Box<ParsedExpr>),
    Center(Box<ParsedExpr>),
    Line(Box<ParsedExpr>),
    Ideal(Label, u128),
    Tensor(Box<[ParsedExpr; 2]>),
    Twist(Box<[ParsedExpr; 2]>),
}

pub(crate) enum Command {
    Pair(Box<[ParsedExpr; 2]>),
    Fujiki(usize, usize, Box<[ParsedExpr; 2]>),
}

pub(crate) fn from_divisor(model: &mut Model, divisor: &Divisor) -> Result<K3Vector> {
    let (&h, exceptional) = divisor
        .coefficients()
        .split_first()
        .expect("divisor has a hyperplane coefficient");
    model.spend(unimodular::RANK)?;
    // In wedge coordinates the Weil polarization is e02 + e13.
    K3Vector::new(
        [0, h, 0, 0, h, 0],
        exceptional.try_into().expect("sixteen nodes"),
    )
}

pub(crate) fn to_divisor(model: &mut Model, basis: &Basis, vector: &K3Vector) -> Result<Divisor> {
    let (wedge, exceptional) = vector.coordinates();
    let h = wedge[1];
    if *wedge != [0, h, 0, 0, h, 0] {
        return Err(Error::plain(
            "K3 lattice vector lies outside the marked Kummer divisor sublattice",
        ));
    }
    let mut coordinates = [0; polarized::COEFFICIENTS];
    coordinates[0] = h;
    coordinates[1..].copy_from_slice(exceptional);
    Divisor::literal(model, basis, coordinates)
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Vector {
    pub(crate) r: i128,
    pub(crate) c: Rc<K3Vector>,
    pub(crate) s: i128,
}

impl Vector {
    pub(crate) fn zero() -> Self {
        Self {
            r: 0,
            c: Rc::new(K3Vector::zero()),
            s: 0,
        }
    }

    pub(crate) fn ideal(length: u128) -> Result<Self> {
        let s = if length == 0 {
            1
        } else {
            Natural::<13>::from_small(length - 1)
                .signed(true)
                .ok_or_else(coordinate_range)?
        };
        Ok(Self {
            r: 1,
            c: Rc::new(K3Vector::zero()),
            s,
        })
    }

    pub(crate) fn line(model: &mut Model, divisor: &Divisor) -> Result<Self> {
        let c = from_divisor(model, divisor)?;
        let s = Integer::from_i128(1)
            .add(half_square(model, &c)?)?
            .signed()
            .ok_or_else(coordinate_range)?;
        Ok(Self {
            r: 1,
            c: Rc::new(c),
            s,
        })
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(2 + unimodular::RANK)?;
        Ok(self.cmp(rhs))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(2 * (2 + unimodular::RANK) + 1)?;
        Ok(Value::Mukai(Rc::new(self)))
    }

    pub(crate) fn persist_quotient(self, model: &mut Model) -> Result<Value> {
        model.allocate(2 * (2 + unimodular::RANK) + 1)?;
        Ok(Value::MukaiQuotient(Rc::new(self)))
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let r = self.r.checked_add(rhs.r).ok_or_else(coordinate_range)?;
        let s = self.s.checked_add(rhs.s).ok_or_else(coordinate_range)?;
        Ok(Self {
            r,
            c: Rc::new(self.c.add(model, &rhs.c)?),
            s,
        })
    }

    pub(crate) fn negate(&self, model: &mut Model) -> Result<Self> {
        let r = self.r.checked_neg().ok_or_else(coordinate_range)?;
        let s = self.s.checked_neg().ok_or_else(coordinate_range)?;
        Ok(Self {
            r,
            c: Rc::new(self.c.negate(model)?),
            s,
        })
    }

    fn wide_pair(&self, model: &mut Model, rhs: &Self) -> Result<Integer> {
        model.spend(60)?;
        self.c
            .wide_pair(model, &rhs.c)?
            .add(Integer::product(self.r, rhs.s).negate())?
            .add(Integer::product(rhs.r, self.s).negate())
    }

    pub(crate) fn pair(&self, model: &mut Model, rhs: &Self) -> Result<i128> {
        self.wide_pair(model, rhs)?
            .signed()
            .ok_or_else(|| Error::plain("Mukai pairing exceeds the signed i128 range"))
    }

    pub(crate) fn tensor(&self, model: &mut Model, divisor: &Divisor) -> Result<Self> {
        let d = from_divisor(model, divisor)?;
        let correction = self
            .c
            .wide_pair(model, &d)?
            .add(half_square(model, &d)?.multiply(self.r)?)?;
        let s = Integer::from_i128(self.s)
            .add(correction)?
            .signed()
            .ok_or_else(coordinate_range)?;
        let c = self.c.add_scaled(model, &d, Integer::from_i128(self.r))?;
        Ok(Self {
            r: self.r,
            c: Rc::new(c),
            s,
        })
    }

    pub(crate) fn twist(&self, model: &mut Model, root: &Self) -> Result<Self> {
        if root.wide_pair(model, root)?.signed() != Some(-2) {
            return Err(Error::plain("rvMukaiTwist requires a vector of square −2"));
        }
        let pairing = self.wide_pair(model, root)?;
        let combine = |value, root| -> Result<i128> {
            Integer::from_i128(value)
                .add(pairing.multiply(root)?)?
                .signed()
                .ok_or_else(coordinate_range)
        };
        let r = combine(self.r, root.r)?;
        let s = combine(self.s, root.s)?;
        let c = self.c.add_scaled(model, &root.c, pairing)?;
        Ok(Self {
            r,
            c: Rc::new(c),
            s,
        })
    }

    pub(crate) fn hilbert_fujiki(
        &self,
        model: &mut Model,
        rhs: &Self,
        n: usize,
        j: usize,
    ) -> Result<i128> {
        if n < 2 {
            return Err(Error::plain("rvMukaiHilbFujiki requires n at least two"));
        }
        let left = polarized::fujiki_left_exponent(n, j)?;
        let hilbert_vector = Self::ideal(n as u128)?;
        if !self.wide_pair(model, &hilbert_vector)?.is_zero()
            || !rhs.wide_pair(model, &hilbert_vector)?.is_zero()
        {
            return Err(Error::plain(
                "rvMukaiHilbFujiki requires both vectors orthogonal to the ideal sheaf vector",
            ));
        }
        let xx = self.wide_pair(model, self)?;
        let xy = self.wide_pair(model, rhs)?;
        let yy = rhs.wide_pair(model, rhs)?;
        polarized::fujiki_pairings(model, left, j, xx, xy, yy)
    }
}

fn half_square(model: &mut Model, divisor: &K3Vector) -> Result<Integer> {
    let square = divisor.wide_pair(model, divisor)?;
    let (magnitude, remainder) = square.magnitude.divide_u32(2);
    if remainder != 0 {
        return Err(Error::plain("K3 lattice must have even squares"));
    }
    Ok(Integer {
        negative: square.negative,
        magnitude,
    })
}

fn coordinate_range() -> Error {
    Error::plain("Mukai coordinate exceeds the signed i128 range")
}
