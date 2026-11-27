use std::cmp::Ordering;
use std::rc::Rc;

use super::nikulin::{self, AFFINE, PAIRS, QUADRATIC};
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};

pub(crate) const RANK: usize = 22;
const GLUE_DIMENSION: usize = QUADRATIC.len() + AFFINE.len();
// A pairing sums at most 22 products of two i128 coordinates. A reflection
// then multiplies its half by one more i128 coordinate, below 2^386.
pub(crate) type Integer = Wide<13>;

pub(crate) enum Expression {
    Literal(Label, Box<([i128; 6], [i128; 16])>),
    Basis(Label, usize),
    Reflect(Box<[ParsedExpr; 2]>),
    Transvection(u8, i128, Box<ParsedExpr>),
}

pub(crate) enum Command {
    Pair(Box<[ParsedExpr; 2]>),
    Gram(Label),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Vector {
    w: [i128; 6],
    b: [i128; 16],
}

impl Vector {
    pub(crate) fn zero() -> Self {
        Self {
            w: [0; 6],
            b: [0; 16],
        }
    }

    pub(crate) fn new(w: [i128; 6], b: [i128; 16]) -> Result<Self> {
        let w_mask = parity_bits(&w) as u8;
        let difference = parity_bits(&b) ^ nikulin::quadratic_table(w_mask);
        if !nikulin::affine_code().any(|table| table == difference) {
            return Err(Error::plain(
                "K3 lattice coordinates do not satisfy the affine gluing condition",
            ));
        }
        Ok(Self { w, b })
    }

    pub(crate) fn coordinates(&self) -> (&[i128; 6], &[i128; 16]) {
        (&self.w, &self.b)
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(RANK)?;
        Ok(self.cmp(rhs))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + RANK)?;
        Ok(Value::K3Lattice(Rc::new(self)))
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        model.spend(RANK)?;
        self.map(|i, value| {
            value
                .checked_add(rhs.coordinate(i))
                .ok_or_else(coordinate_range)
        })
    }

    pub(crate) fn negate(&self, model: &mut Model) -> Result<Self> {
        model.spend(RANK)?;
        self.map(|_, value| value.checked_neg().ok_or_else(coordinate_range))
    }

    pub(crate) fn pair(&self, model: &mut Model, rhs: &Self) -> Result<i128> {
        self.wide_pair(model, rhs)?
            .signed()
            .ok_or_else(|| Error::plain("K3 lattice pairing exceeds the signed i128 range"))
    }

    pub(crate) fn wide_pair(&self, model: &mut Model, rhs: &Self) -> Result<Integer> {
        model.spend(30 * RANK)?;
        let mut sum = Integer::from_i128(0);
        // Complementary wedge pairs, with orientation making e02+e13 positive.
        for i in 0..6 {
            let term = Integer::product(self.w[i], rhs.w[5 - i]);
            sum = sum.add(if i == 1 || i == 4 {
                term
            } else {
                term.negate()
            })?;
        }
        for (&left, &right) in self.b.iter().zip(&rhs.b) {
            sum = sum.add(Integer::product(left, right).negate())?;
        }
        let (magnitude, remainder) = sum.magnitude.divide_u32(2);
        if remainder != 0 {
            return Err(Error::plain("K3 lattice pairing must be integral"));
        }
        Ok(Integer {
            magnitude,
            negative: sum.negative,
        })
    }

    pub(crate) fn reflect(&self, model: &mut Model, root: &Self) -> Result<Self> {
        if root.wide_pair(model, root)?.signed() != Some(-2) {
            return Err(Error::plain("rvK3Reflect requires a root of square −2"));
        }
        let pairing = self.wide_pair(model, root)?;
        self.add_scaled(model, root, pairing)
    }

    pub(crate) fn add_scaled(&self, model: &mut Model, rhs: &Self, scale: Integer) -> Result<Self> {
        model.spend(70 * RANK)?;
        self.map(|i, value| {
            Integer::from_i128(value)
                .add(scale.multiply(rhs.coordinate(i))?)?
                .signed()
                .ok_or_else(coordinate_range)
        })
    }

    pub(crate) fn transvection(&self, model: &mut Model, mask: u8, amount: i128) -> Result<Self> {
        model.spend(45 * RANK)?;
        let u: [i128; 4] = std::array::from_fn(|i| i128::from((mask >> i) & 1));
        let g = [u[2], u[3], -u[0], -u[1]]; // g=Ju
        let mut contractions = [Integer::from_i128(0); 4];
        for (&(i, j), &value) in PAIRS.iter().zip(&self.w) {
            let (i, j) = (usize::from(i), usize::from(j));
            contractions[j] = contractions[j].add(Integer::product(u[i], value))?;
            contractions[i] = contractions[i].add(Integer::product(u[j], value).negate())?;
        }
        // M^(-t)=I−amount*g*u^t. The degree-two term on the wedge vanishes
        // since the perturbation has rank one. contractions[j]=(u^t w)_j.
        let mut w = self.w;
        for (&(i, j), value) in PAIRS.iter().zip(&mut w) {
            let (i, j) = (usize::from(i), usize::from(j));
            let change = contractions[j]
                .multiply(g[i])?
                .add(contractions[i].multiply(g[j])?.negate())?;
            *value = Integer::from_i128(*value)
                .add(change.multiply(amount)?.negate())?
                .signed()
                .ok_or_else(coordinate_range)?;
        }
        let mut b = self.b;
        if amount & 1 != 0 {
            let dual = (mask >> 2) | ((mask & 3) << 2);
            for node in 0..16u8 {
                let image = if (dual & node).count_ones() & 1 == 1 {
                    node ^ mask
                } else {
                    node
                };
                b[usize::from(image)] = self.b[usize::from(node)];
            }
        }
        Self::new(w, b)
    }

    fn map(&self, mut apply: impl FnMut(usize, i128) -> Result<i128>) -> Result<Self> {
        let mut result = Self::zero();
        for (i, entry) in result.w.iter_mut().chain(&mut result.b).enumerate() {
            *entry = apply(i, self.coordinate(i))?;
        }
        Ok(result)
    }

    fn coordinate(&self, i: usize) -> i128 {
        if i < self.w.len() {
            self.w[i]
        } else {
            self.b[i - self.w.len()]
        }
    }
}

fn coordinate_range() -> Error {
    Error::plain("K3 lattice coordinate exceeds the signed i128 range")
}

fn parity_bits(coordinates: &[i128]) -> u16 {
    coordinates
        .iter()
        .enumerate()
        .fold(0, |mask, (i, &entry)| mask | (((entry & 1) as u16) << i))
}

fn integral_basis(model: &mut Model) -> Result<Vec<Vector>> {
    model.spend(RANK * GLUE_DIMENSION * 2)?;
    let mut rows = [0u32; GLUE_DIMENSION];
    for (i, &table) in QUADRATIC.iter().enumerate() {
        rows[i] = (1 << i) | (u32::from(table) << 6);
    }
    for (row, &table) in rows[6..].iter_mut().zip(&AFFINE) {
        *row = u32::from(table) << 6;
    }
    let pivots = nikulin::binary_rref(&mut rows, RANK);
    if pivots.len() != GLUE_DIMENSION {
        return Err(Error::plain(
            "K3 lattice gluing code has the wrong dimension",
        ));
    }
    // In doubled coordinates each RREF pivot is 1. Doubling its row recovers
    // the original 2e_i using the 2e_j at nonpivots.
    let mut result = Vec::with_capacity(RANK);
    for row in rows {
        let coordinate = |i| i128::from((row >> i) & 1u32);
        result.push(Vector::new(
            std::array::from_fn(coordinate),
            std::array::from_fn(|i| coordinate(i + 6)),
        )?);
    }
    for column in (0..RANK).filter(|column| !pivots.contains(column)) {
        result.push(Vector::zero().map(|i, _| Ok(if i == column { 2 } else { 0 }))?);
    }
    Ok(result)
}

pub(crate) fn basis_vector(model: &mut Model, index: usize) -> Result<Vector> {
    integral_basis(model)?
        .into_iter()
        .nth(index)
        .ok_or_else(|| Error::plain("rvK3LatticeBasis index must be between 0 and 21"))
}

pub(crate) fn gram_report(model: &mut Model, name: &Label) -> Result<String> {
    nikulin::basis_for(
        model,
        name,
        "K3 lattices require a certified genus-two level-2 basis",
    )?;
    let gram = integral_gram(model)?;
    let index = 1i128 << GLUE_DIMENSION;
    // The unglued U(2)^3 ⊕ A_1^16 has determinant −2^22. Passing to an
    // overlattice of this index divides its determinant by the square index.
    let determinant = -(1i128 << RANK) / (index * index);
    Ok(format!("rvK3LatticeGram[{}] = rank={RANK}, signature=(3, 19), index={index}, determinant={determinant}, gram={gram:?}", name.text))
}

pub(crate) fn integral_gram(model: &mut Model) -> Result<[[i128; RANK]; RANK]> {
    let basis = integral_basis(model)?;
    let mut gram = [[0i128; RANK]; RANK];
    for (i, left) in basis.iter().enumerate() {
        for (j, right) in basis[..=i].iter().enumerate() {
            let pair = left.pair(model, right)?;
            if i == j && pair % 2 != 0 {
                return Err(Error::plain("K3 lattice basis has an odd square"));
            }
            gram[i][j] = pair;
            gram[j][i] = pair;
        }
    }
    Ok(gram)
}
