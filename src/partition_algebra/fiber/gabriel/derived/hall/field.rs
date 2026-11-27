use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::syntax::Label;

#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) enum FieldName {
    Prime(u32),
    Residue(Rc<str>),
}

impl FieldName {
    pub(super) fn render(&self) -> String {
        match self {
            Self::Prime(p) => p.to_string(),
            Self::Residue(name) => name.to_string(),
        }
    }
}

pub(crate) enum FieldSpec {
    Prime(u128),
    Residue(Label),
}

impl FieldSpec {
    pub(super) fn resolve(&self, model: &mut Model) -> Result<FieldName> {
        match self {
            Self::Prime(p) => {
                let p = u32::try_from(*p)
                    .map_err(|_| Error::plain("partition Hall prime must fit u32"))?;
                if p != 2 {
                    Prime::new(model, p as usize).map_err(|_| {
                        Error::plain("partition Hall coefficient characteristic must be prime")
                    })?;
                }
                Ok(FieldName::Prime(p))
            }
            Self::Residue(name) => {
                model.rv_unram(name)?;
                Ok(FieldName::Residue(name.text.as_str().into()))
            }
        }
    }
}

pub(super) enum Field {
    Prime(Prime),
    Table(Table),
}

pub(super) struct Table {
    q: usize,
    add: Vec<u32>,
    multiply: Vec<u32>,
    negative: Vec<u32>,
    inverse: Vec<u32>,
}

impl Field {
    pub(super) fn new(model: &mut Model, name: &FieldName) -> Result<Self> {
        let FieldName::Residue(name) = name else {
            let FieldName::Prime(p) = name else {
                unreachable!()
            };
            return Ok(Self::Prime(Prime::validated(*p)));
        };
        let ring = model
            .rv_unram_named(name)
            .ok_or_else(|| Error::plain("unknown partition Hall residue field"))?;
        let q = ring.size_of_residue_field()?;
        if q > 64 {
            return Err(Error::plain(
                "named partition Hall residue fields are limited to size 64",
            ));
        }
        let q = q as usize;
        let ring = ring.residue_ring(model)?;
        let p = ring.prime();
        let encode = |value: Vec<u64>| -> u32 {
            value
                .into_iter()
                .rev()
                .fold(0, |row, digit| row * p + digit) as u32
        };
        let scalars: Vec<_> = (0..q)
            .map(|i| {
                let mut digits = ring.coordinate_row(i);
                digits.reverse();
                digits
            })
            .collect();
        let mut add = Vec::with_capacity(q * q);
        let mut multiply = Vec::with_capacity(q * q);
        let mut negative = Vec::with_capacity(q);
        let mut inverse = vec![0];
        for (i, x) in scalars.iter().enumerate() {
            negative.push(encode(ring.negate(model, x)?));
            if i != 0 {
                inverse.push(encode(ring.inverse(model, x)?));
            }
            for y in &scalars {
                add.push(encode(ring.add(model, x, y)?));
                multiply.push(encode(ring.multiply(model, x, y)?));
            }
        }
        // Base-p enumeration puts the scalar one in the least significant place.
        if encode(ring.scalar(1)) != 1 || (1..q).any(|i| multiply[i * q + inverse[i] as usize] != 1)
        {
            return Err(Error::plain(
                "partition Hall field table failed inverse certification",
            ));
        }
        model.allocate(2 * q * q + 2 * q)?;
        Ok(Self::Table(Table {
            q,
            add,
            multiply,
            negative,
            inverse,
        }))
    }

    pub(super) fn order(&self) -> u32 {
        match self {
            Self::Prime(p) => p.get(),
            Self::Table(t) => t.q as u32,
        }
    }

    pub(super) fn add(&self, a: u32, b: u32) -> u32 {
        match self {
            Self::Prime(p) => p.add(a, b),
            Self::Table(t) => t.add[a as usize * t.q + b as usize],
        }
    }

    pub(super) fn mul(&self, a: u32, b: u32) -> u32 {
        match self {
            Self::Prime(p) => p.mul(a, b),
            Self::Table(t) => t.multiply[a as usize * t.q + b as usize],
        }
    }

    pub(super) fn sub(&self, a: u32, b: u32) -> u32 {
        match self {
            Self::Prime(p) => p.add(a, p.negate(b)),
            Self::Table(t) => t.add[a as usize * t.q + t.negative[b as usize] as usize],
        }
    }

    pub(super) fn inverse(&self, model: &mut Model, x: u32) -> Result<u32> {
        debug_assert_ne!(x, 0);
        if x == 1 {
            return Ok(1);
        }
        match self {
            Self::Table(t) => Ok(t.inverse[x as usize]),
            Self::Prime(p) => {
                let mut exponent = p.get() - 2;
                let mut base = x;
                let mut value = 1;
                model.spend(2 * u32::BITS as usize)?;
                while exponent != 0 {
                    if exponent & 1 == 1 {
                        value = p.mul(value, base);
                    }
                    base = p.mul(base, base);
                    exponent >>= 1;
                }
                Ok(value)
            }
        }
    }
}
