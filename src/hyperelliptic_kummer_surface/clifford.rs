use std::collections::BTreeMap;
use std::rc::Rc;

use super::rational_sparse::{Coefficients, Domain, Integer, Magnitude, WORK_LIMBS};
use super::unimodular::{Vector, RANK};
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

pub(crate) const VOLUME_MASK: u32 = (1 << RANK) - 1;
// The first six slots are the U,V pairs for (01,23), (02,13), (03,12).
const WEDGE_MASK: u32 = 0x3f;
pub(super) const NEGATIVE_SQUARES: u32 = (VOLUME_MASK ^ WEDGE_MASK) | 0b011001;

pub(crate) enum Expression {
    One(Label),
    Volume(Label),
    Vector(Box<ParsedExpr>),
    Reverse(Box<ParsedExpr>),
    Scale(i128, u128, Box<ParsedExpr>),
    Multiply(Box<[ParsedExpr; 2]>),
    Act(Box<[ParsedExpr; 2]>),
}

pub(crate) struct Coefficient(pub(crate) u32, pub(crate) Box<ParsedExpr>);

pub(crate) enum Clifford {}

impl Domain for Clifford {
    const NAME: &'static str = "Clifford";
}

pub(crate) type Element = Coefficients<Clifford>;

impl Element {
    pub(crate) fn one() -> Self {
        Self::singleton(0, 1)
    }

    pub(crate) fn volume() -> Self {
        Self::singleton(VOLUME_MASK, 8)
    }

    pub(crate) fn from_vector(model: &mut Model, vector: &Vector) -> Result<Self> {
        model.spend(RANK * WORK_LIMBS)?;
        let terms = vector_numerators(vector)?;
        Self::normalized(
            model,
            terms.iter().map(|(mask, entry)| (*mask, entry)),
            Magnitude::from_small(4),
        )
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        self.charge_storage(model)?;
        Ok(Value::K3Clifford(Rc::new(self)))
    }

    pub(crate) fn reverse(&self, model: &mut Model) -> Result<Self> {
        self.signed_map(model, reverse_sign)
    }

    pub(crate) fn multiply(&self, model: &mut Model, other: &Self) -> Result<Self> {
        let left = self.widen(model, false)?;
        let right = other.widen(model, false)?;
        let terms = product(model, &left, &right)?;
        let denominator = self.product_denominator(other)?;
        Self::normalized(
            model,
            terms.iter().map(|(&mask, entry)| (mask, entry)),
            denominator,
        )
    }

    pub(crate) fn act(&self, model: &mut Model, vector: &Vector, parity: u8) -> Result<Vector> {
        // Cancel the common denominator of g against that of its reverse.
        let left = self.widen(model, false)?;
        let reverse = self.widen(model, true)?;
        let norm = product(model, &left, &reverse)?;
        let Some(norm) = norm.get(&0).filter(|_| norm.len() == 1) else {
            return Err(Error::plain(
                "rvCliffordAct requires g times its reverse to be a nonzero scalar",
            ));
        };
        let first = product(model, &left, &vector_numerators(vector)?)?;
        let intermediate: Vec<_> = first.into_iter().collect();
        let image = product(model, &intermediate, &reverse)?;
        if image.keys().any(|mask| mask.count_ones() != 1) {
            return Err(Error::plain(
                "rvCliffordAct twisted conjugation did not return a vector",
            ));
        }
        // The embedded vector had denominator four. Recovering doubled lattice
        // coordinates from the diagonal Clifford generators multiplies by two.
        let denominator = norm
            .magnitude
            .multiply_u128(2)
            .ok_or_else(Self::arithmetic)?;
        let numerator = |i| {
            image
                .get(&(1u32 << i))
                .copied()
                .unwrap_or_else(|| Integer::from_i128(0))
        };
        let coordinate = |model: &mut Model, value: Integer| -> Result<i128> {
            let magnitude = value
                .magnitude
                .divide_exact(model, denominator)?
                .ok_or_else(|| {
                    Error::plain("rvCliffordAct image has a nonintegral K3 lattice coordinate")
                })?;
            magnitude
                .signed(value.negative ^ norm.negative ^ (parity == 1))
                .ok_or_else(|| {
                    Error::plain("rvCliffordAct image coordinate exceeds signed i128 range")
                })
        };
        let mut w = [0; 6];
        let mut b = [0; 16];
        for i in 0..3 {
            let u = numerator(2 * i);
            let v = numerator(2 * i + 1);
            w[i] = coordinate(model, u.add(v)?)?;
            w[5 - i] = coordinate(model, u.add(v.negate())?)?;
        }
        for (i, entry) in b.iter_mut().enumerate() {
            *entry = coordinate(model, numerator(6 + i))?;
        }
        Vector::new(w, b)
    }

    fn widen(&self, model: &mut Model, reverse: bool) -> Result<Vec<(u32, Integer)>> {
        model.spend(WORK_LIMBS * self.terms().len())?;
        self.terms()
            .iter()
            .map(|&(mask, value)| {
                let value = value.resize()?;
                Ok((
                    mask,
                    if reverse && reverse_sign(mask) {
                        value.negate()
                    } else {
                        value
                    },
                ))
            })
            .collect()
    }
}

fn vector_numerators(vector: &Vector) -> Result<Vec<(u32, Integer)>> {
    let (w, b) = vector.coordinates();
    let mut terms = Vec::with_capacity(RANK);
    for i in 0..3 {
        let left = Integer::from_i128(w[i]);
        let right = Integer::from_i128(w[5 - i]);
        terms.push((1 << (2 * i), left.add(right)?));
        terms.push((1 << (2 * i + 1), left.add(right.negate())?));
    }
    for (i, &entry) in b.iter().enumerate() {
        terms.push((1 << (6 + i), Integer::from_i128(entry).multiply(2)?));
    }
    terms.retain(|(_, value)| !value.is_zero());
    Ok(terms)
}

fn reverse_sign(mask: u32) -> bool {
    (mask.count_ones() & 2) != 0
}

fn product(
    model: &mut Model,
    left: &[(u32, Integer)],
    right: &[(u32, Integer)],
) -> Result<BTreeMap<u32, Integer>> {
    model
        .spend(
            left.len()
                .saturating_mul(right.len())
                .saturating_mul(4 * WORK_LIMBS),
        )
        .map_err(|_| Element::work_limit())?;
    let mut terms = BTreeMap::new();
    for &(a, x) in left {
        for &(b, y) in right {
            let overlap = a & b;
            // Move the right factors to their ordered slots, then square repeats.
            let mut odd = (overlap & NEGATIVE_SQUARES).count_ones() & 1;
            let mut pending = a;
            while pending != 0 {
                let slot = pending.trailing_zeros();
                odd ^= (b & ((1u32 << slot) - 1)).count_ones() & 1;
                pending &= pending - 1;
            }
            let scale = 1i128 << (overlap & WEDGE_MASK).count_ones();
            let coefficient =
                x.multiply_wide(y)?
                    .multiply(if odd == 0 { scale } else { -scale })?;
            Element::accumulate(&mut terms, a ^ b, coefficient)?;
        }
    }
    Ok(terms)
}
