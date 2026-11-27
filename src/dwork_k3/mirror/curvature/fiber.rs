use std::fmt::Debug;

use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::model::Model;

use super::{integrity, Data};

impl Data {
    pub(super) fn specialize<R: CoefficientRing>(
        &self,
        model: &mut Model,
        rank: usize,
        ring: &R,
        point: &R::Element,
        lift: impl Fn(u32) -> R::Element,
    ) -> Result<Matrix<R::Element>> {
        if ring.is_zero(point) {
            return Err(Error::plain(
                "Dwork p-curvature requires a regular fiber with z nonzero",
            ));
        }
        let d = evaluate(
            model,
            ring,
            &[1, self.prime.negate(self.prime.reduce(256))],
            point,
            &lift,
        )?;
        if ring.is_zero(&d) {
            return Err(Error::plain(
                "Dwork p-curvature excludes the singular fiber 1-256z=0",
            ));
        }
        let inverse = ring.inverse_unit(model, &d)?;
        let scale = scalar_power(model, ring, inverse, self.prime.get())?;
        let entries = self
            .numerator(rank)
            .entries()
            .iter()
            .map(|polynomial| {
                let value = evaluate(model, ring, polynomial, point, &lift)?;
                ring.multiply(model, &value, &scale)
            })
            .collect::<Result<_>>()?;
        Matrix::from_entries(rank, entries)
    }

    pub(super) fn fiber<R: CoefficientRing>(
        &self,
        model: &mut Model,
        rank: usize,
        ring: &R,
        point: &R::Element,
        lift: impl Fn(u32) -> R::Element,
    ) -> Result<String>
    where
        R::Element: Debug,
    {
        let matrix = self.specialize(model, rank, ring, point, &lift)?;
        let rows = matrix.entries().chunks(rank).collect::<Vec<_>>();
        let mut ranks = vec![rank];
        let mut kernels = Vec::new();
        let mut current = matrix.clone();
        for exponent in 1..=rank {
            let basis = kernel(model, ring, &current)?;
            ranks.push(rank - basis.len());
            kernels.push(basis);
            if exponent < rank {
                current = hw_algebra::multiply(model, ring, &current, &matrix)?;
            }
        }
        if ranks[rank] != 0 {
            return Err(integrity());
        }
        let index = ranks.iter().position(|&r| r == 0).ok_or_else(integrity)?;
        let jet = self
            .hasse
            .iter()
            .map(|coefficient| evaluate(model, ring, coefficient, point, &lift))
            .collect::<Result<Vec<_>>>()?;
        let order = jet
            .iter()
            .position(|c| !ring.is_zero(c))
            .ok_or_else(integrity)?;
        let mut description = format!(
            "matrix={rows:?}, power ranks={ranks:?}, nilpotency index={index}, kernel bases of Psi^1 through Psi^{rank}={kernels:?}, Hasse jet={jet:?}, Hasse order={order}"
        );
        if rank == 3 {
            if !annihilates(model, ring, &matrix, &jet)? {
                return Err(integrity());
            }
            let line = if kernels[0].len() == 1 {
                "ker(Psi)"
            } else {
                "a line in ker(Psi)"
            };
            let contact = match order {
                0 => "nonzero first coordinate",
                1 => "first coordinate vanishes and second is nonzero",
                _ => "third coordinate line",
            };
            description.push_str(&format!(
                ", Hasse line={line}, coordinate contact={contact}"
            ));
        }
        Ok(description)
    }
}

fn evaluate<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    polynomial: &[u32],
    point: &R::Element,
    lift: &impl Fn(u32) -> R::Element,
) -> Result<R::Element> {
    let mut value = ring.zero();
    for &coefficient in polynomial.iter().rev() {
        value = ring.multiply(model, &value, point)?;
        if coefficient != 0 {
            value = ring.add(model, &value, &lift(coefficient))?;
        }
    }
    Ok(value)
}

fn scalar_power<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    mut value: R::Element,
    mut exponent: u32,
) -> Result<R::Element> {
    let mut result = ring.one();
    while exponent != 0 {
        if exponent & 1 == 1 {
            result = ring.multiply(model, &result, &value)?;
        }
        exponent >>= 1;
        if exponent != 0 {
            value = ring.multiply(model, &value, &value)?;
        }
    }
    Ok(result)
}

pub(super) fn iterate<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    matrix: Matrix<R::Element>,
    exponent: usize,
) -> Result<Matrix<R::Element>> {
    let rank = matrix.size();
    match exponent {
        0 => hw_algebra::identity(model, ring, rank),
        1 => Ok(matrix),
        k if k >= rank => Matrix::from_entries(rank, vec![ring.zero(); rank * rank]),
        _ => hw_algebra::multiply(model, ring, &matrix, &matrix),
    }
}

fn kernel<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    matrix: &Matrix<R::Element>,
) -> Result<Vec<Vec<R::Element>>> {
    let rank = matrix.size();
    let mut rows = matrix
        .entries()
        .chunks(rank)
        .map(<[_]>::to_vec)
        .collect::<Vec<_>>();
    let mut pivots = Vec::new();
    for col in 0..rank {
        let pivot_row = pivots.len();
        let Some(row) = (pivot_row..rank).find(|&i| !ring.is_zero(&rows[i][col])) else {
            continue;
        };
        rows.swap(row, pivot_row);
        let scale = ring.inverse_unit(model, &rows[pivot_row][col])?;
        for coefficient in &mut rows[pivot_row][col..] {
            *coefficient = ring.multiply(model, coefficient, &scale)?;
        }
        let pivot = rows[pivot_row].clone();
        for (i, row) in rows.iter_mut().enumerate() {
            if i != pivot_row && !ring.is_zero(&row[col]) {
                let scale = row[col].clone();
                for (coefficient, value) in row[col..].iter_mut().zip(&pivot[col..]) {
                    let term = ring.multiply(model, &scale, value)?;
                    *coefficient = ring.subtract(model, coefficient, &term)?;
                }
            }
        }
        pivots.push(col);
    }
    let mut basis = Vec::new();
    for free in (0..rank).filter(|col| !pivots.contains(col)) {
        let mut vector = vec![ring.zero(); rank];
        vector[free] = ring.one();
        for (row, &pivot) in pivots.iter().enumerate() {
            vector[pivot] = ring.subtract(model, &ring.zero(), &rows[row][free])?;
        }
        if !annihilates(model, ring, matrix, &vector)? {
            return Err(integrity());
        }
        basis.push(vector);
    }
    Ok(basis)
}

fn annihilates<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    matrix: &Matrix<R::Element>,
    vector: &[R::Element],
) -> Result<bool> {
    for row in matrix.entries().chunks(matrix.size()) {
        let mut value = ring.zero();
        for (a, b) in row.iter().zip(vector) {
            let term = ring.multiply(model, a, b)?;
            value = ring.add(model, &value, &term)?;
        }
        if !ring.is_zero(&value) {
            return Ok(false);
        }
    }
    Ok(true)
}
