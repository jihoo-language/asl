use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::rc::Rc;

use crate::ec_morphism::Rational;
use crate::ec_poly::{Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::finite_field_factor::factor_with_multiplicity;
use crate::hw_algebra::CoefficientRing;
use crate::model::{Model, Value};

use super::{field, from_polynomial, sum, Pencil, Section};

#[derive(Clone, Default, Eq, Ord, PartialEq, PartialOrd)]
struct SquareClass {
    nonsquare: bool,
    support: Vec<Polynomial>,
}

#[derive(Default, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Descent([SquareClass; 3]);

impl Descent {
    pub(super) fn from_coordinates(coordinates: [(bool, Vec<Polynomial>); 3]) -> Self {
        Self(coordinates.map(|(nonsquare, support)| SquareClass { nonsquare, support }))
    }

    pub(super) fn from_representatives(
        model: &mut Model,
        field: Field<'_>,
        representatives: [Polynomial; 3],
    ) -> Result<Self> {
        let mut result = Self::default();
        let mut parity = false;
        let mut odd_factors = BTreeSet::new();
        for (coordinate, polynomial) in result.0.iter_mut().zip(representatives) {
            if polynomial.is_empty() {
                return Err(Error::plain(
                    "K3 pencil squareclass representative must be nonzero",
                ));
            }
            *coordinate = SquareClass::of(model, field, &from_polynomial(field, polynomial))?;
            parity ^= coordinate.nonsquare;
            for factor in &coordinate.support {
                if !odd_factors.insert(factor.clone()) {
                    odd_factors.remove(factor);
                }
            }
        }
        if parity || !odd_factors.is_empty() {
            return Err(Error::plain(
                "K3 pencil Kummer triple must have square product",
            ));
        }
        Ok(result)
    }

    pub(crate) fn coordinates(&self) -> impl Iterator<Item = (bool, &[Polynomial])> {
        self.0
            .iter()
            .map(|coordinate| (coordinate.nonsquare, coordinate.support.as_slice()))
    }

    pub(crate) fn is_trivial(&self) -> bool {
        self.0
            .iter()
            .all(|coordinate| !coordinate.nonsquare && coordinate.support.is_empty())
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(self.cells())?;
        Ok(Value::K3PencilKummer2(Rc::new(self)))
    }

    fn cells(&self) -> usize {
        4 + self.0.iter().map(SquareClass::cells).sum::<usize>()
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(self.cells() + rhs.cells())?;
        Ok(self.cmp(rhs))
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let mut result = Self::default();
        for ((out, left), right) in result.0.iter_mut().zip(&self.0).zip(&rhs.0) {
            let mut factors = BTreeSet::new();
            model.spend(
                (left.cells() + right.cells())
                    * (1 + (left.support.len() + right.support.len()).max(1).ilog2() as usize),
            )?;
            for factor in left.support.iter().chain(&right.support) {
                if !factors.insert(factor) {
                    factors.remove(factor);
                }
            }
            out.nonsquare = left.nonsquare ^ right.nonsquare;
            out.support = factors.into_iter().cloned().collect();
        }
        Ok(result)
    }
}

impl Section {
    pub(crate) fn descent(&self, model: &mut Model, pencil: &Pencil) -> Result<Descent> {
        let field = field(&pencil.curve);
        let roots = pencil.torsion_roots(model)?;
        let [a, b, c] = roots else {
            return Err(Error::plain(
                "K3 pencil two-descent requires all three affine two-torsion points",
            ));
        };
        let Self::Affine(x, _) = self else {
            return Ok(Descent::default());
        };
        let roots = [a, b, c].map(|root| from_polynomial(field, root.clone()));
        let mut result = Descent::default();
        for (i, coordinate) in result.0.iter_mut().enumerate() {
            let difference = sum(model, field, x, &roots[i], true)?;
            let representative = if difference.numerator.is_empty() {
                let j = (i + 1) % 3;
                let k = (i + 2) % 3;
                let a = sum(model, field, x, &roots[j], true)?;
                let b = sum(model, field, x, &roots[k], true)?;
                a.multiply(model, field, &b)?
            } else {
                difference
            };
            *coordinate = SquareClass::of(model, field, &representative)?;
        }
        Ok(result)
    }
}

impl SquareClass {
    fn cells(&self) -> usize {
        2 + self
            .support
            .iter()
            .map(|factor| {
                1 + factor
                    .iter()
                    .map(|coordinate| 3 + coordinate.len())
                    .sum::<usize>()
            })
            .sum::<usize>()
    }

    fn of(model: &mut Model, field: Field<'_>, function: &Rational) -> Result<Self> {
        let leading = function
            .numerator
            .last()
            .ok_or_else(|| Error::plain("zero has no K3 pencil square class"))?;
        let mut support = Vec::new();
        // The reduced denominator is monic. Its odd valuations also
        // contribute at finite primes, since −1=1 in the parity group.
        for part in [&function.numerator, &function.denominator] {
            support.extend(
                factor_with_multiplicity(model, field, part)?
                    .into_iter()
                    .filter_map(|(_, order, factor)| (!order.is_multiple_of(2)).then_some(factor)),
            );
        }
        support.sort_unstable();
        Ok(Self {
            nonsquare: !is_square(model, field, leading)?,
            support,
        })
    }
}

pub(super) fn is_square(model: &mut Model, field: Field<'_>, value: &Scalar) -> Result<bool> {
    // The norm has the same square character over the prime field, even
    // when q does not fit a machine word.
    let mut norm = value.clone();
    let mut conjugate = value.clone();
    for _ in 1..field.ring.degree() {
        conjugate = field.ring.frobenius(model, &conjugate, 1)?;
        norm = field.multiply(model, &norm, &conjugate)?;
    }
    let character = field
        .ring
        .power(model, &norm, (field.ring.prime() - 1) / 2)?;
    match crate::rv_unram::require_scalar(&character)? {
        1 => Ok(true),
        c if c == field.ring.prime() - 1 => Ok(false),
        _ => Err(Error::plain(
            "K3 pencil square class has an invalid quadratic character",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rv_unram::Ring;

    #[test]
    fn norm_gives_the_constant_square_character_over_an_extension() {
        let mut model = Model::new();
        let ring =
            Ring::new(&mut model, 3, 1, vec![1, 0, 1]).unwrap_or_else(|error| panic!("{error}"));
        let field = Field { ring: &ring };
        assert!(is_square(&mut model, field, &vec![0, 1]).unwrap_or_else(|error| panic!("{error}")));
        assert!(
            !is_square(&mut model, field, &vec![1, 1]).unwrap_or_else(|error| panic!("{error}"))
        );
    }
}
