use std::fmt;

use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::hilbert::{
    elliptic_genus::{decimal, half_k3_coefficient},
    Integer, Magnitude,
};
use crate::model::Model;

/// Coordinates doubled in the basis (a, b, r) of monomials P^a q^b y^r.
#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Weight(pub(crate) [i128; 3]);

const SIMPLE: [Weight; 3] = [Weight([0, 0, -2]), Weight([2, 0, 2]), Weight([0, 2, 2])];

impl fmt::Display for Weight {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        let [a, b, r] = self.0;
        write!(out, "{a}, {b}, {r}")
    }
}

impl Weight {
    pub(crate) const ZERO: Self = Self([0; 3]);

    pub(crate) fn add(self, rhs: Self) -> Result<Self> {
        let mut sum = self;
        for (left, right) in sum.0.iter_mut().zip(rhs.0) {
            *left = left.checked_add(right).ok_or_else(coordinate_overflow)?;
        }
        Ok(sum)
    }

    pub(crate) fn negate(self) -> Result<Self> {
        let mut opposite = self;
        for component in &mut opposite.0 {
            *component = component.checked_neg().ok_or_else(coordinate_overflow)?;
        }
        Ok(opposite)
    }

    pub(crate) fn reflect(self, index: usize) -> Result<Self> {
        let root = SIMPLE[index];
        let pairing = self.twice_pair(root)?;
        let mut image = self;
        for (coordinate, root_coordinate) in image.0.iter_mut().zip(root.0) {
            // The root coordinates are doubled. Cancel that two with the
            // factor in twice_pair before rounding to the coordinate range.
            let correction = pairing.multiply(root_coordinate / 2)?;
            *coordinate = Integer::from_i128(*coordinate)
                .add(correction.negate())?
                .signed()
                .ok_or_else(coordinate_overflow)?;
        }
        Ok(image)
    }

    pub(crate) fn pair_report(self, other: Self) -> Result<String> {
        let value = self.twice_pair(other)?;
        let (magnitude, remainder) = value.magnitude.divide_u32(2);
        let pairing = if remainder == 0 {
            decimal(Integer { magnitude, ..value })
        } else {
            format!("{}/2", decimal(value))
        };
        Ok(format!("rvSiegelPair = {pairing}"))
    }

    fn twice_pair(self, other: Self) -> Result<Integer> {
        let [a, b, r] = self.0;
        let [c, d, s] = other.0;
        let cross = Integer::from_i128(a)
            .multiply(d)?
            .add(Integer::from_i128(b).multiply(c)?)?;
        Integer::from_i128(r).multiply(s)?.add(cross.multiply(-2)?)
    }

    pub(crate) fn root_report(self, model: &mut Model) -> Result<String> {
        if self.0.iter().any(|coordinate| coordinate % 2 != 0) {
            return Err(Error::plain(
                "rvSiegelRoot requires integral weight coordinates",
            ));
        }
        let [a, b, r] = self.0.map(|coordinate| coordinate / 2);
        let c0 = Integer::from_i128(a)
            .add(Integer::from_i128(b))?
            .add(Integer::from_i128(r).negate())?;
        let simple = [c0, Integer::from_i128(a), Integer::from_i128(b)];
        let orientation = if self == Self::ZERO {
            "neither"
        } else if simple
            .iter()
            .all(|coordinate| !coordinate.negative || coordinate.is_zero())
        {
            "positive"
        } else if simple
            .iter()
            .all(|coordinate| coordinate.negative || coordinate.is_zero())
        {
            "negative"
        } else {
            "neither"
        };
        let square = self.twice_pair(self)?;
        let (magnitude, _) = square.magnitude.divide_u32(2);
        let norm = Integer {
            magnitude,
            ..square
        };
        let multiplicity = if orientation == "neither" {
            Integer::from_i128(0)
        } else {
            root_multiplicity(model, norm, r)?
        };
        let kind = if multiplicity.is_zero() {
            "absent"
        } else if norm.is_zero() {
            "null"
        } else if norm.negative {
            "imaginary"
        } else {
            "real"
        };
        Ok(format!(
            "rvSiegelRoot = {kind}, orientation={orientation}, simple=({}, {a}, {b}), norm={}, supermultiplicity={}",
            decimal(c0),
            decimal(norm),
            decimal(multiplicity)
        ))
    }
}

fn root_multiplicity(model: &mut Model, norm: Integer, r: i128) -> Result<Integer> {
    // Elliptic translation brings r to ε ∈ {0,1} and ab to (4ab-r²+ε)/4.
    // The weak Jacobi form has no coefficients below discriminant -1.
    if !norm.negative && norm.magnitude.small().is_none_or(|value| value > 2) {
        return Ok(Integer::from_i128(0));
    }
    let parity = r.rem_euclid(2);
    let (magnitude, _) = norm.magnitude.divide_u32(2);
    let discriminant = Integer {
        magnitude,
        ..norm.negate()
    };
    let degree = reduced_jacobi_degree(discriminant, parity)?.ok_or_else(|| {
        Error::plain("Siegel root coefficient exceeds the elliptic genus coefficient storage limit")
    })?;
    half_k3_coefficient(model, degree, parity)
}

pub(crate) fn discriminant(a: u128, b: u128, r: i128) -> Result<Integer> {
    let product = Integer {
        magnitude: Magnitude::from_small(a),
        negative: false,
    }
    .multiply_wide(Integer {
        magnitude: Magnitude::from_small(b),
        negative: false,
    })?;
    product.multiply(4)?.add(Integer::product(r, r).negate())
}

pub(crate) fn reduced_jacobi_degree(discriminant: Integer, parity: i128) -> Result<Option<u128>> {
    let shifted = discriminant.add(Integer::from_i128(parity))?;
    let (degree, remainder) = shifted.magnitude.divide_u32(4);
    if remainder != 0 || (shifted.negative && !shifted.is_zero()) {
        return Err(Error::plain(
            "invalid index-one Jacobi discriminant and parity",
        ));
    }
    Ok(degree.small())
}

fn coordinate_overflow() -> Error {
    Error::plain("Siegel weight doubled coordinate exceeds signed 128-bit range")
}
