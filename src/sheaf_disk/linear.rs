use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::{self, Basis, Coordinate};
use crate::sheaf::{self, Matrix};

pub(super) fn literal(
    model: &mut Model,
    prime: Prime,
    height: usize,
    width: usize,
    values: Vec<i64>,
) -> Result<Matrix> {
    let cells = height.saturating_mul(width);
    sheaf::check_cells(cells)?;
    if values.len() != cells {
        return Err(Error::plain(format!(
            "disk matrix requires {height} rows and {width} columns ({cells} entries), found {}",
            values.len()
        )));
    }
    Matrix::literal(model, prime, width, values)
}

pub(super) fn transpose(model: &mut Model, matrix: &Matrix) -> Result<Matrix> {
    model.spend(matrix.cells())?;
    let mut rows = BTreeMap::<usize, Vec<Coordinate>>::new();
    for (column, support) in matrix.columns() {
        for &(row, value) in support {
            rows.entry(row).or_default().push((column, value));
        }
    }
    Matrix::from_columns(
        rows.into_iter()
            .map(|(i, entries)| (i, entries.into()))
            .collect(),
    )
}

pub(super) fn negate(model: &mut Model, p: Prime, matrix: &Matrix) -> Result<Matrix> {
    model.spend(matrix.cells())?;
    Matrix::from_columns(
        matrix
            .columns()
            .map(|(index, values)| {
                (
                    index,
                    values.iter().map(|&(i, a)| (i, p.negate(a))).collect(),
                )
            })
            .collect(),
    )
}

pub(super) fn one_minus(
    model: &mut Model,
    p: Prime,
    rank: usize,
    product: &Matrix,
) -> Result<Matrix> {
    Matrix::from_action(model, rank, |model, j| {
        prime_linear::add_scaled(model, p, &[(j, 1)], product.column(j), p.negate(1))
    })
}

pub(super) fn inverse(model: &mut Model, p: Prime, rank: usize, matrix: &Matrix) -> Result<Matrix> {
    let mut basis = Basis::default();
    for j in 0..rank {
        if basis
            .insert_with_witness(model, p, matrix.column(j).to_vec(), vec![(j, 1)])?
            .is_some()
        {
            return Err(Error::plain("nearby monodromy I - v u is singular"));
        }
    }
    coordinates(model, p, &basis, rank, |j| vec![(j, 1)])
}

pub(super) fn rank(model: &mut Model, p: Prime, matrix: &Matrix) -> Result<usize> {
    let mut basis = Basis::default();
    for (_, col) in matrix.columns() {
        basis.insert(model, p, col.to_vec())?;
    }
    Ok(basis.rank())
}

pub(super) struct Image {
    pub(super) projection: Matrix,
    pub(super) inclusion: Matrix,
    pub(super) basis: Basis,
}

pub(super) fn image_factorization(
    model: &mut Model,
    p: Prime,
    width: usize,
    matrix: &Matrix,
) -> Result<(Matrix, Matrix, usize)> {
    let image = image(model, p, width, matrix)?;
    Ok((image.projection, image.inclusion, image.basis.rank()))
}

pub(super) fn image(model: &mut Model, p: Prime, width: usize, matrix: &Matrix) -> Result<Image> {
    let mut basis = Basis::default();
    let mut included = Vec::new();
    for (_, column) in matrix.columns() {
        if basis
            .insert_with_witness(model, p, column.to_vec(), vec![(included.len(), 1)])?
            .is_none()
        {
            included.push((included.len(), column.to_vec().into()));
        }
    }
    let quotient = coordinates(model, p, &basis, width, |j| matrix.column(j).to_vec())?;
    let inclusion = Matrix::from_columns(included)?;
    Ok(Image {
        projection: quotient,
        inclusion,
        basis,
    })
}

pub(super) fn factor_through(
    model: &mut Model,
    p: Prime,
    inclusion: &Matrix,
    width: usize,
    map: &Matrix,
) -> Result<Matrix> {
    let mut basis = Basis::default();
    for (j, column) in inclusion.columns() {
        let dependent = basis.insert_with_witness(model, p, column.to_vec(), vec![(j, 1)])?;
        debug_assert!(dependent.is_none());
    }
    coordinates(model, p, &basis, width, |j| map.column(j).to_vec())
}

fn coordinates(
    model: &mut Model,
    p: Prime,
    basis: &Basis,
    width: usize,
    mut column: impl FnMut(usize) -> Vec<Coordinate>,
) -> Result<Matrix> {
    Matrix::from_action(model, width, |model, j| {
        let (remainder, entries) = basis.reduce_with_witness(model, p, column(j))?;
        if !remainder.is_empty() {
            return Err(Error::plain("disk map is outside the required image"));
        }
        Ok(entries)
    })
}

pub(super) fn report(
    model: &mut Model,
    matrix: &Matrix,
    height: usize,
    width: usize,
) -> Result<String> {
    let size = height.saturating_mul(width);
    sheaf::check_cells(size + height)?;
    model.spend(size + height)?;
    let mut rows = vec![vec![0u32; width]; height];
    for (column, entries) in matrix.columns() {
        for &(row, value) in entries {
            rows[row][column] = value;
        }
    }
    Ok(format!("{rows:?}"))
}
