use crate::elliptic_modular::decimal;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_witt_spectrum::{character::StableCharacter, SchurShape};

use super::{intermediate_overflow, Integer};

const CHARACTER_CHECK_RANK: usize = 8;

pub(crate) enum Report {
    Dimensions(usize, u128),
    Centralizer(usize, u128),
}

impl Report {
    pub(crate) fn execute(self, model: &mut Model) -> Result<String> {
        match self {
            Self::Dimensions(rank, n) => {
                let (generic, finite) = dimensions(model, 2 * rank, n)?;
                Ok(format!(
                    "rvPartitionDimensions[{rank}, {n}] = generic={}, finite={}",
                    decimal(generic),
                    decimal(finite)
                ))
            }
            Self::Centralizer(rank, n) => {
                if rank > CHARACTER_CHECK_RANK {
                    return Err(Error::plain(
                        "partition centralizer character check exceeds rank 8",
                    ));
                }
                let (generic, finite) = dimensions(model, 2 * rank, n)?;
                let unit = StableCharacter::seed(SchurShape::new(Vec::new())?);
                let standard = StableCharacter::seed(SchurShape::new(vec![1])?);
                let permutation = unit.add(model, &standard)?;
                let mut tensor = unit;
                for _ in 0..rank {
                    tensor = tensor.product(model, &permutation)?;
                }
                let specialized = tensor.specialize(model, n)?;
                let sum = specialized.squared_multiplicities(model)?;
                if sum.magnitude != finite.magnitude {
                    return Err(Error::plain("partition centralizer dimension disagrees with the character multiplicities"));
                }
                Ok(format!("rvPartitionCentralizerCheck[{rank}, {n}] = generic={}, finite={}, multiplicity squares={}",decimal(generic),decimal(finite),decimal(sum)))
            }
        }
    }
}

pub(super) fn dimensions(
    model: &mut Model,
    vertices: usize,
    n: u128,
) -> Result<(Integer, Integer)> {
    let mut row = vec![Integer::from_i128(0); vertices + 1];
    row[0] = Integer::from_i128(1);
    for length in 1..=vertices {
        model.spend(length + 1)?;
        for blocks in (1..=length).rev() {
            row[blocks] = row[blocks]
                .multiply(blocks as i128)
                .and_then(|value| value.add(row[blocks - 1]))
                .map_err(|_| intermediate_overflow())?;
        }
        row[0] = Integer::from_i128(0);
    }
    let sum = |end: usize| {
        row[..=end].iter().try_fold(Integer::from_i128(0), |a, &b| {
            a.add(b).map_err(|_| intermediate_overflow())
        })
    };
    Ok((sum(vertices)?, sum(n.min(vertices as u128) as usize)?))
}
