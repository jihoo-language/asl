use std::collections::BTreeMap;
use std::mem::size_of;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::model::Model;
use crate::syntax::Label;

use super::{
    hilbert::{checked_product, graded_product, intermediate_overflow, Integer},
    stratum::Partition,
};

const SURFACE_RANK: usize = 24;

pub(crate) struct Declaration(pub Label, pub Label, pub Vec<(u128, Vec<u128>)>);

pub(crate) struct Decoration {
    pub(crate) shape: Rc<Partition>,
    diagrams: BTreeMap<u128, Box<[u128]>>,
}

impl Declaration {
    pub(crate) fn execute(self, model: &mut Model) -> Result<()> {
        let Self(name, shape, diagrams) = self;
        let shape = model.hilbert_stratum(&shape)?;
        let mut assigned = BTreeMap::new();
        for (size, rows) in diagrams {
            model.spend(1 + rows.len() + shape.parts().len().max(1).ilog2() as usize)?;
            let index = shape
                .parts()
                .binary_search_by_key(&size, |&(m, _)| m)
                .map_err(|_| {
                    Error::plain("Specht diagram size must occur in the Hilbert–Chow partition")
                })?;
            if rows.is_empty() || rows.contains(&0) || rows.windows(2).any(|pair| pair[0] < pair[1])
            {
                return Err(Error::plain(
                    "Specht diagram rows must be positive and weakly decreasing",
                ));
            }
            let count = rows
                .iter()
                .try_fold(0u128, |sum, row| sum.checked_add(*row));
            if count != Some(shape.parts()[index].1) {
                return Err(Error::plain(
                    "Specht diagram size must equal the part multiplicity",
                ));
            }
            if assigned.insert(size, rows.into_boxed_slice()).is_some() {
                return Err(Error::plain(
                    "Specht diagram for this part size is already specified",
                ));
            }
        }
        // A one-row diagram is the default, so only store nontrivial ones.
        assigned.retain(|_, rows| rows.len() > 1);
        model.define_specht_stratum(
            name,
            Decoration {
                shape,
                diagrams: assigned,
            },
        )
    }
}

impl Decoration {
    pub(crate) fn storage(&self) -> usize {
        2 + self
            .diagrams
            .values()
            .map(|rows| 3 + 2 * rows.len())
            .sum::<usize>()
    }

    pub(super) fn supports_degree(&self, model: &mut Model, mut degree: u128) -> Result<bool> {
        // Every box below the first row has half-degree at least one.
        // In a full-height column the last box has half-degree two.
        for rows in self.diagrams.values() {
            model.spend(rows.len())?;
            if rows.len() > SURFACE_RANK {
                return Ok(false);
            }
            let extra =
                rows[1..].iter().sum::<u128>() + rows.get(SURFACE_RANK - 1).copied().unwrap_or(0);
            let Some(left) = degree.checked_sub(extra) else {
                return Ok(false);
            };
            degree = left;
        }
        Ok(true)
    }

    pub(super) fn coefficient(
        &self,
        model: &mut Model,
        symmetric: &[Vec<Integer>],
    ) -> Result<Integer> {
        let limit = (symmetric.len() - 1) as u128;
        let width = symmetric[0].len();
        let mut product = vec![Integer::from_i128(1)];
        for &(size, count) in self.shape.parts() {
            product = if let Some(rows) = self.diagrams.get(&size) {
                let schur = jacobi_trudi(model, rows, symmetric, width)?;
                graded_product(model, &product, &schur, width)?
            } else {
                graded_product(
                    model,
                    &product,
                    &symmetric[count.min(limit) as usize],
                    width,
                )?
            };
        }
        Ok(product
            .get(width - 1)
            .copied()
            .unwrap_or_else(|| Integer::from_i128(0)))
    }
}

fn jacobi_trudi(
    model: &mut Model,
    diagram: &[u128],
    symmetric: &[Vec<Integer>],
    width: usize,
) -> Result<Vec<Integer>> {
    let size = diagram.len();
    let ring = GradedRing(width);
    if diagram[0] == 1 {
        // A single column is exterior power. Use E(t)H(−t)=1, the
        // triangular Jacobi–Trudi relation, even for a full height column.
        let limit = symmetric.len() - 1;
        hw_algebra::storage(2 * (size + 1), ring.cells()).map_err(|_| {
            Error::plain("Specht determinant exceeds the coefficient storage limit")
        })?;
        let mut exterior = vec![ring.one()];
        for n in 1..=size {
            let mut value = ring.zero();
            for j in 1..=n {
                let term = Series::new(graded_product(
                    model,
                    &symmetric[j.min(limit)],
                    &exterior[n - j].0,
                    width,
                )?);
                value = ring.combine(model, &value, &term, j.is_multiple_of(2))?;
            }
            exterior.push(value);
        }
        return Ok(exterior.pop().unwrap().0);
    }
    hw_algebra::storage(size * size + 5 * (size + 1), ring.cells())
        .map_err(|_| Error::plain("Specht determinant exceeds the coefficient storage limit"))?;
    model.spend(size * size * width)?;
    let limit = (symmetric.len() - 1) as u128;
    let mut entries = Vec::with_capacity(size * size);
    for (i, &row) in diagram.iter().enumerate() {
        for j in 0..size {
            // A Jacobi–Trudi index is at most the size of the diagram.
            let index = if j >= i {
                Some(row + (j - i) as u128)
            } else {
                row.checked_sub((i - j) as u128)
            };
            entries.push(match index {
                Some(index) => Series::new(symmetric[index.min(limit) as usize].clone()),
                None => ring.zero(),
            });
        }
    }
    let matrix = Matrix::from_entries(size, entries)?;
    let mut determinant = hw_algebra::characteristic(model, &ring, &matrix)?
        .pop()
        .unwrap()
        .0;
    if !size.is_multiple_of(2) {
        for value in &mut determinant {
            *value = value.negate();
        }
    }
    Ok(determinant)
}

#[derive(Clone)]
struct Series(Vec<Integer>);

impl Series {
    fn new(mut values: Vec<Integer>) -> Self {
        while values.last().is_some_and(|value| value.is_zero()) {
            values.pop();
        }
        Self(values)
    }
}

impl PartialEq for Series {
    fn eq(&self, other: &Self) -> bool {
        self.0.len() == other.0.len()
            && self.0.iter().zip(&other.0).all(|(a, b)| {
                a.magnitude == b.magnitude && (a.is_zero() || a.negative == b.negative)
            })
    }
}

impl Eq for Series {}

struct GradedRing(usize);

impl GradedRing {
    fn combine(&self, model: &mut Model, a: &Series, b: &Series, subtract: bool) -> Result<Series> {
        let len = a.0.len().max(b.0.len());
        model.spend(len)?;
        let mut values = vec![Integer::from_i128(0); len];
        values[..a.0.len()].copy_from_slice(&a.0);
        for (value, &term) in values.iter_mut().zip(&b.0) {
            let term = if subtract { term.negate() } else { term };
            *value = value.add(term).map_err(|_| intermediate_overflow())?;
        }
        Ok(Series::new(values))
    }
}

impl CoefficientRing for GradedRing {
    type Element = Series;

    fn cells(&self) -> usize {
        self.0 * size_of::<Integer>().div_ceil(size_of::<usize>()) + 3
    }

    fn zero(&self) -> Series {
        Series(Vec::new())
    }

    fn one(&self) -> Series {
        Series(vec![Integer::from_i128(1)])
    }

    fn is_zero(&self, value: &Series) -> bool {
        value.0.is_empty()
    }

    fn is_unit(&self, value: &Series) -> bool {
        matches!(
            value.0.first().and_then(|coefficient| coefficient.signed()),
            Some(1 | -1)
        )
    }

    fn add(&self, model: &mut Model, a: &Series, b: &Series) -> Result<Series> {
        self.combine(model, a, b, false)
    }

    fn subtract(&self, model: &mut Model, a: &Series, b: &Series) -> Result<Series> {
        self.combine(model, a, b, true)
    }

    fn multiply(&self, model: &mut Model, a: &Series, b: &Series) -> Result<Series> {
        Ok(Series::new(graded_product(model, &a.0, &b.0, self.0)?))
    }

    fn inverse_unit(&self, model: &mut Model, value: &Series) -> Result<Series> {
        let Some(sign @ (1 | -1)) = value.0.first().and_then(|coefficient| coefficient.signed())
        else {
            return Err(Error::plain("Specht graded series has no integral inverse"));
        };
        let mut reciprocal = vec![Integer::from_i128(sign)];
        for _ in 1..self.0 {
            let mut sum = Integer::from_i128(0);
            for (&a, &b) in value.0.iter().skip(1).zip(reciprocal.iter().rev()) {
                let term = checked_product(model, a, b)?;
                sum = sum.add(term).map_err(|_| intermediate_overflow())?;
            }
            reciprocal.push(if sign == 1 { sum.negate() } else { sum });
        }
        Ok(Series::new(reciprocal))
    }
}
