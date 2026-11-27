use std::cmp::Ordering;
use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::model::Model;

use super::orthogonal::{self, DIMENSION};
use super::{arithmetic, decimal, guard, reduced, Magnitude, WEYL_LIMBS};

struct Restriction {
    target: usize,
    memo: [HashMap<Vec<u128>, Magnitude>; DIMENSION + 1],
    cells: usize,
}

impl Restriction {
    fn new(rho: u128) -> Result<Self> {
        if !(1..=(DIMENSION - 2) as u128).contains(&rho) {
            return Err(Error::plain("Hodge polarization rank must be in 1..=21"));
        }
        Ok(Self {
            target: DIMENSION - rho as usize,
            memo: Default::default(),
            cells: 0,
        })
    }

    fn invariant(
        &mut self,
        model: &mut Model,
        dimension: usize,
        shape: &[u128],
    ) -> Result<Magnitude> {
        let steps = dimension - self.target;
        if shape.is_empty() {
            return Ok(Magnitude::from_small(1));
        }
        // One interlacing can remove at most the last nonzero coordinate.
        if shape.len() > steps {
            return Ok(Magnitude::default());
        }
        if steps == 1 {
            return Ok(Magnitude::from_small(1));
        }
        if let Some(&count) = self.memo[dimension].get(shape) {
            return Ok(count);
        }
        let mut count = Magnitude::default();
        orthogonal::branch(
            model,
            dimension,
            shape,
            steps - 1,
            &mut |model, child, signs| {
                let value = self.invariant(model, dimension - 1, child)?;
                model.spend(WEYL_LIMBS)?;
                // The two D chiralities are distinct, but their restrictions to
                // the next odd orthogonal group depend on the same absolute
                // coordinates. At SO(2) only charge zero survives as an invariant.
                count.add_multiple(value, signs).ok_or_else(arithmetic)?;
                Ok(())
            },
        )?;
        self.cells += WEYL_LIMBS + 4 * shape.len() + 8;
        guard(self.cells)?;
        self.memo[dimension].insert(shape.to_vec(), count);
        Ok(count)
    }
}

pub(super) fn multiplicity(model: &mut Model, rho: u128, shape: Vec<u128>) -> Result<String> {
    let shape = reduced::validate_b11(shape)?;
    let count = Restriction::new(rho)?.invariant(model, DIMENSION, &shape)?;
    Ok(format!(
        "rvHilbHodgeInvariantMultiplicity[{rho}]{shape:?} = {}",
        decimal(model, count)?
    ))
}

pub(super) fn dimension(
    model: &mut Model,
    n: u128,
    d: i128,
    rho: u128,
    primitive: bool,
) -> Result<String> {
    let mut restriction = Restriction::new(rho)?;
    let mut positive = Magnitude::default();
    let mut negative = Magnitude::default();
    reduced::visit_summands(model, n, d, primitive, |model, shape, multiplicity| {
        let invariant = restriction.invariant(model, DIMENSION, shape)?;
        let contribution = invariant
            .multiply_u128(multiplicity.unsigned_abs())
            .ok_or_else(arithmetic)?;
        let sum = if multiplicity >= 0 {
            &mut positive
        } else {
            &mut negative
        };
        sum.add_multiple(contribution, 1).ok_or_else(arithmetic)?;
        Ok(())
    })?;
    if positive.cmp(negative) == Ordering::Less {
        return Err(Error::plain(
            "primitive LLV Hodge invariant dimension is negative",
        ));
    }
    let name = if primitive {
        "rvHilbPrimitiveHodgeInvariant"
    } else {
        "rvHilbHodgeInvariant"
    };
    let count = decimal(model, positive.difference(negative))?;
    Ok(format!("{name}[{n}, {d}, {rho}] = {count}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_last_interlacing_has_pairs_of_nonzero_torus_weights() {
        let mut model = Model::new();
        let mut weights = Vec::new();
        orthogonal::branch(&mut model, 3, &[2], 1, &mut |_, shape, signs| {
            weights.push((shape.to_vec(), signs));
            Ok(())
        })
        .unwrap_or_else(|err| panic!("{err}"));
        assert_eq!(weights, [(vec![], 1), (vec![1], 2), (vec![2], 2)]);

        let mut restriction = Restriction::new(21).unwrap_or_else(|err| panic!("{err}"));
        assert_eq!(
            restriction
                .invariant(&mut model, 2, &[1])
                .unwrap_or_else(|err| panic!("{err}"))
                .small(),
            Some(0)
        );
        assert_eq!(
            restriction
                .invariant(&mut model, 3, &[2])
                .unwrap_or_else(|err| panic!("{err}"))
                .small(),
            Some(1)
        );
    }
}
