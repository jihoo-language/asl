use crate::error::{Error, Result};
use crate::model::Model;

use super::super::multiply;
use super::table::{
    factorial, from_power_basis, power_basis, CharacterOracle, Coefficients, MAX_WEIGHT,
};
use super::Character;

impl Character {
    pub(crate) fn kronecker(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        debug_assert_eq!(self.weight, rhs.weight);
        if self.terms.is_empty() || rhs.terms.is_empty() {
            return Ok(Self::zero(self.weight));
        }
        if let Some(scalar) = self.trivial_scalar() {
            return rhs.scale(model, scalar);
        }
        if let Some(scalar) = rhs.trivial_scalar() {
            return self.scale(model, scalar);
        }
        if self.weight > u128::from(MAX_WEIGHT) {
            return Err(Error::plain(
                "ordinary Kronecker product exceeds the weight limit of 16 boxes",
            ));
        }
        let degree = self.weight as u8;
        let mut oracle = CharacterOracle::default();
        let left = self.class_values(model, &mut oracle, degree)?;
        let right = rhs.class_values(model, &mut oracle, degree)?;
        let product = left
            .into_iter()
            .filter_map(|(class, a)| {
                right
                    .get(&class)
                    .map(|&b| multiply(model, a, b).map(|value| (class, value)))
            })
            .collect::<Result<Coefficients>>()?;
        let basis = power_basis(model, degree, product)?;
        from_power_basis(model, degree, factorial(degree), &basis, &mut oracle)
    }

    fn trivial_scalar(&self) -> Option<i128> {
        let (shape, &coefficient) = self.terms.first_key_value()?;
        (self.terms.len() == 1 && shape.rows().len() <= 1).then_some(coefficient)
    }
}
