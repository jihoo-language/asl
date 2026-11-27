use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_ore::gcd;

use super::super::sections::Section;
use super::{
    closed_places, field, Field, Kodaira, MinimalReduction, Pencil, RelativeJacobian,
    DISCRIMINANT_DEGREE,
};

pub(super) struct HeightData {
    pub(super) multiplier: u128,
    finite: Vec<(usize, String, usize)>,
    infinity: String,
}

impl Pencil {
    pub(crate) fn section_height_report(
        &self,
        model: &mut Model,
        name: &str,
        section: &Section,
        other: Option<&Section>,
    ) -> Result<String> {
        let heights = self.relative(model)?.heights(model, field(&self.curve))?;
        let left = heights.numerator(model, self, section)?;
        let square = heights.multiplier * heights.multiplier;
        let (operation, value) = if let Some(other) = other {
            let right = heights.numerator(model, self, other)?;
            let pair = heights.twice_pair(model, self, section, other, left, right)?;
            ("Pair", signed_ratio(pair, 2 * square))
        } else {
            ("Height", signed_ratio(left, square))
        };
        Ok(format!(
            "rvK3PencilSection{operation}[{name}] = {value}, {}",
            heights.provenance(),
        ))
    }
}

impl RelativeJacobian {
    pub(super) fn heights(&self, model: &mut Model, field: Field<'_>) -> Result<HeightData> {
        if field.ring.prime() < 5 {
            return Err(Error::plain(
                "K3 pencil section heights require characteristic at least five",
            ));
        }
        let discriminant = &self.0.discriminant;
        let infinity = self.minimal_at_factor(model, field, None)?;
        let mut degree = infinity.delta;
        let mut multiplier = exponent(&infinity)?;
        let mut finite = Vec::new();
        for (closed_degree, order, polynomial) in closed_places(model, field, discriminant)? {
            let reduction = self.minimal_at_factor(model, field, Some(&polynomial))?;
            let exponent = exponent(&reduction)?;
            if reduction.delta != order {
                return Err(inconsistent_model());
            }
            degree += closed_degree * order;
            multiplier = lcm(multiplier, exponent);
            finite.push((closed_degree, reduction.kodaira.symbol(), order));
        }
        if degree != DISCRIMINANT_DEGREE {
            return Err(inconsistent_model());
        }
        Ok(HeightData {
            multiplier,
            finite,
            infinity: infinity.kodaira.symbol(),
        })
    }
}

impl HeightData {
    pub(super) fn provenance(&self) -> String {
        format!(
            "component multiplier={}, minimal discriminant degree=24, finite fibers (degree, type, discriminant order)={:?}, infinity={}",
            self.multiplier, self.finite, self.infinity,
        )
    }

    pub(super) fn numerator(
        &self,
        model: &mut Model,
        pencil: &Pencil,
        point: &Section,
    ) -> Result<i128> {
        let multiple = point.multiple(model, pencil, self.multiplier)?;
        let Section::Affine(x, _) = multiple else {
            return Ok(0);
        };
        let numerator_degree = x.numerator.len().saturating_sub(1);
        let denominator_degree = x.denominator.len() - 1;
        // A finite x pole of order 2k contributes k to (mP.O).
        // At infinity use s^4 x, so the excess above four gives the
        // remaining contribution. mP meets only identity components.
        let height = numerator_degree.max(4 + denominator_degree) as i128;
        if !height.unsigned_abs().is_multiple_of(self.multiplier) {
            return Err(Error::plain(
                "K3 pencil section height failed its local denominator bound",
            ));
        }
        Ok(height)
    }

    pub(super) fn twice_pair(
        &self,
        model: &mut Model,
        pencil: &Pencil,
        left: &Section,
        right: &Section,
        left_height: i128,
        right_height: i128,
    ) -> Result<i128> {
        let sum = left.add(model, pencil, right)?;
        let pair = self.numerator(model, pencil, &sum)? - left_height - right_height;
        if !pair.unsigned_abs().is_multiple_of(2 * self.multiplier) {
            return Err(Error::plain(
                "K3 pencil section pairing failed its local denominator bound",
            ));
        }
        Ok(pair)
    }
}

fn exponent(reduction: &MinimalReduction) -> Result<u128> {
    if reduction.scales != 0 {
        return Err(Error::plain(
            "K3 pencil section heights require a globally minimal Weierstrass model",
        ));
    }
    component_exponent(reduction.kodaira)
}

pub(super) fn component_exponent(kind: Kodaira) -> Result<u128> {
    match kind {
        Kodaira::Good => Ok(1),
        Kodaira::Multiplicative(n) => Ok(n as u128),
        Kodaira::Star(0) => Ok(2),
        _ => Err(Error::plain(
            "K3 pencil section heights support geometric I_n and I0* fibers",
        )),
    }
}

pub(super) fn lcm(a: u128, b: u128) -> u128 {
    // The component numbers are bounded by the degree-24 discriminant.
    a / gcd(a, b) * b
}

pub(super) fn signed_ratio(numerator: i128, denominator: u128) -> String {
    let divisor = gcd(numerator.unsigned_abs(), denominator);
    let residue = numerator / divisor as i128;
    let denominator = denominator / divisor;
    if denominator == 1 {
        residue.to_string()
    } else {
        format!("{residue}/{denominator}")
    }
}

fn inconsistent_model() -> Error {
    Error::plain("K3 pencil section height model failed its minimal discriminant check")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn component_exponents_require_minimal_supported_fibers() {
        let minimal = |c4, c6, d| {
            MinimalReduction::new(Some(c4), Some(c6), Some(d))
                .unwrap_or_else(|error| panic!("{error}"))
        };
        assert!(matches!(exponent(&minimal(0, 0, 4)), Ok(4)));
        assert!(matches!(exponent(&minimal(2, 3, 6)), Ok(2)));
        assert!(exponent(&minimal(4, 6, 12)).is_err());
        assert!(exponent(&minimal(2, 3, 7)).is_err());
    }
}
