use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;

pub(crate) struct Matrix {
    rows: usize,
    columns: usize,
    entries: Vec<i128>,
}

impl Matrix {
    pub(crate) fn new(rows: usize, columns: usize) -> Result<Self> {
        let entries = rows
            .checked_mul(columns)
            .filter(|&size| size <= MAX_PIVOT_CELLS / 2)
            .ok_or_else(|| Error::plain("BP integral matrix exceeds the storage limit"))?;
        Ok(Self {
            rows,
            columns,
            entries: vec![0; entries],
        })
    }

    pub(crate) fn row(&self, row: usize) -> &[i128] {
        let start = row * self.columns;
        &self.entries[start..start + self.columns]
    }

    pub(crate) fn set(&mut self, row: usize, column: usize, value: i128) {
        self.entries[row * self.columns + column] = value;
    }

    pub(crate) fn augmented(&self, model: &mut Model, target: &[i128]) -> Result<Self> {
        debug_assert_eq!(target.len(), self.rows);
        let mut augmented = Self::new(self.rows, self.columns + 1)?;
        model.spend(self.entries.len() + self.rows)?;
        for (row, &coefficient) in target.iter().enumerate() {
            let new = row * augmented.columns;
            augmented.entries[new..new + self.columns].copy_from_slice(self.row(row));
            augmented.entries[new + self.columns] = coefficient;
        }
        Ok(augmented)
    }

    pub(crate) fn apply(&self, model: &mut Model, vector: &[i128]) -> Result<Vec<i128>> {
        debug_assert_eq!(vector.len(), self.columns);
        model.spend(self.entries.len() + self.rows)?;
        let mut result = vec![0i128; self.rows];
        for (row, target) in result.iter_mut().enumerate() {
            for (&entry, &coefficient) in self.row(row).iter().zip(vector) {
                let value = entry
                    .checked_mul(coefficient)
                    .ok_or_else(witness_overflow)?;
                *target = target.checked_add(value).ok_or_else(witness_overflow)?;
            }
        }
        Ok(result)
    }
}

pub(crate) fn witness_overflow() -> Error {
    Error::plain("BP integral witness arithmetic exceeds the signed i128 range")
}

fn incompatible_rank() -> Error {
    Error::plain("BP primary reduction disagrees with the rational rank")
}

// The incoming maps have p-primary cokernel torsion. Finding every pivot
// of the known rational rank certifies their invariant factors over Z.
pub(crate) fn primary_factors(
    model: &mut Model,
    matrix: &Matrix,
    prime: u32,
    rank: usize,
) -> Result<Vec<i128>> {
    if rank > matrix.rows.min(matrix.columns) {
        return Err(incompatible_rank());
    }
    if rank == 0 {
        return Ok(Vec::new());
    }
    let prime = u64::from(prime);
    let mut modulus = prime * prime;
    loop {
        if let Some(factors) = primary_diagonal(model, matrix, prime, modulus, rank)? {
            return Ok(factors);
        }
        let previous = modulus;
        // Double the precision, or use the largest power fitting a u64.
        modulus = if let Some(square) = modulus.checked_mul(modulus) {
            square
        } else {
            while let Some(power) = modulus.checked_mul(prime) {
                modulus = power;
            }
            modulus
        };
        if modulus == previous {
            return Err(Error::plain(
                "BP primary Smith precision exceeds the u64 modulus limit",
            ));
        }
    }
}

struct ModularMatrix {
    entries: Vec<u64>,
    rows: usize,
    columns: usize,
    prime: u64,
    modulus: u64,
    rhs: Option<Vec<u64>>,
    permutation: Option<Vec<usize>>,
}

fn primary_diagonal(
    model: &mut Model,
    matrix: &Matrix,
    prime: u64,
    modulus: u64,
    rank: usize,
) -> Result<Option<Vec<i128>>> {
    let mut reduced = ModularMatrix::new(model, matrix, prime, modulus)?;
    let mut factors = Vec::with_capacity(rank);
    for pivot in 0..rank {
        let Some((row, column, divisor)) = reduced.find_pivot(model, pivot)? else {
            return Ok(None);
        };
        reduced.swap_pivot(model, pivot, row, column)?;
        reduced.clear_column(model, pivot, divisor)?;
        factors.push(i128::from(divisor));
    }
    Ok(Some(factors))
}

pub(crate) fn modular_preimage(
    model: &mut Model,
    matrix: &Matrix,
    prime: u32,
    modulus: u64,
    target: &[i128],
    factor: u64,
    rank: usize,
) -> Result<Vec<i128>> {
    if rank > matrix.rows.min(matrix.columns) {
        return Err(incompatible_rank());
    }
    let mut reduced = ModularMatrix::new(model, matrix, u64::from(prime), modulus)?;
    model.spend(target.len() + matrix.columns)?;
    reduced.rhs = Some(
        target
            .iter()
            .map(|&coefficient| {
                let residue = coefficient.rem_euclid(i128::from(modulus)) as u64;
                multiply_modulo(residue, factor, modulus)
            })
            .collect(),
    );
    reduced.permutation = Some((0..matrix.columns).collect());
    let mut diagonal = Vec::with_capacity(rank);
    for pivot in 0..rank {
        let (row, column, divisor) = reduced
            .find_pivot(model, pivot)?
            .ok_or_else(incompatible_rank)?;
        reduced.swap_pivot(model, pivot, row, column)?;
        let inverse = reduced.clear_column(model, pivot, divisor)?;
        diagonal.push((divisor, inverse));
    }
    reduced.solution(model, &diagonal)
}

impl ModularMatrix {
    fn new(model: &mut Model, matrix: &Matrix, prime: u64, modulus: u64) -> Result<Self> {
        model.spend(matrix.entries.len())?;
        Ok(Self {
            entries: matrix
                .entries
                .iter()
                .map(|&value| value.rem_euclid(i128::from(modulus)) as u64)
                .collect(),
            rows: matrix.rows,
            columns: matrix.columns,
            prime,
            modulus,
            rhs: None,
            permutation: None,
        })
    }

    fn solution(self, model: &mut Model, diagonal: &[(u64, u64)]) -> Result<Vec<i128>> {
        let rhs = self.rhs.as_ref().unwrap();
        let permutation = self.permutation.unwrap();
        let rank = diagonal.len();
        model.spend(rhs.len() + self.columns)?;
        if rhs[rank..].iter().any(|&residue| residue != 0) {
            return Err(no_modular_solution());
        }
        let mut coordinates = vec![0; self.columns];
        for (row, &(divisor, inverse)) in diagonal.iter().enumerate().rev() {
            let mut residual = rhs[row];
            let base = row * self.columns;
            model.spend(rank - row)?;
            for (column, &coefficient) in coordinates.iter().enumerate().take(rank).skip(row + 1) {
                let product =
                    multiply_modulo(coefficient, self.entries[base + column], self.modulus);
                residual = subtract_modulo(residual, product, self.modulus);
            }
            if !residual.is_multiple_of(divisor) {
                return Err(no_modular_solution());
            }
            coordinates[row] = multiply_modulo(residual / divisor, inverse, self.modulus / divisor);
        }
        let mut result = vec![0; self.columns];
        for (position, value) in permutation.into_iter().zip(coordinates) {
            result[position] = if value <= self.modulus / 2 {
                i128::from(value)
            } else {
                -(i128::from(self.modulus - value))
            };
        }
        Ok(result)
    }

    fn find_pivot(&self, model: &mut Model, pivot: usize) -> Result<Option<(usize, usize, u64)>> {
        let mut best = None;
        let mut best_divisor = self.modulus;
        for row in pivot..self.rows {
            model.spend(self.columns - pivot)?;
            for column in pivot..self.columns {
                let mut value = self.entries[row * self.columns + column];
                if value == 0 {
                    continue;
                }
                let mut divisor = 1;
                while value.is_multiple_of(self.prime) {
                    model.spend(1)?;
                    value /= self.prime;
                    divisor *= self.prime;
                }
                if divisor == 1 {
                    return Ok(Some((row, column, 1)));
                }
                if divisor < best_divisor {
                    best = Some((row, column, divisor));
                    best_divisor = divisor;
                }
            }
        }
        Ok(best)
    }

    fn swap_pivot(
        &mut self,
        model: &mut Model,
        pivot: usize,
        row: usize,
        column: usize,
    ) -> Result<()> {
        if row != pivot {
            model.spend(self.columns - pivot)?;
            for next in pivot..self.columns {
                self.entries
                    .swap(pivot * self.columns + next, row * self.columns + next);
            }
            if let Some(rhs) = &mut self.rhs {
                rhs.swap(pivot, row);
            }
        }
        if column != pivot {
            let first = if self.permutation.is_some() { 0 } else { pivot };
            model.spend(self.rows - first)?;
            for next in first..self.rows {
                self.entries
                    .swap(next * self.columns + pivot, next * self.columns + column);
            }
            if let Some(columns) = &mut self.permutation {
                columns.swap(pivot, column);
            }
        }
        Ok(())
    }

    fn clear_column(&mut self, model: &mut Model, pivot: usize, divisor: u64) -> Result<u64> {
        let small_modulus = self.modulus / divisor;
        let diagonal = self.entries[pivot * self.columns + pivot] / divisor;
        let inverse = inverse_modulo(model, diagonal, small_modulus)?;
        model.spend(self.rows - pivot)?;
        for row in pivot + 1..self.rows {
            let base = row * self.columns;
            let value = self.entries[base + pivot];
            if value == 0 {
                continue;
            }
            let factor = multiply_modulo(value / divisor, inverse, small_modulus);
            model.spend(self.columns - pivot)?;
            for col in pivot + 1..self.columns {
                let product = multiply_modulo(
                    factor,
                    self.entries[pivot * self.columns + col],
                    self.modulus,
                );
                let residue = self.entries[base + col];
                self.entries[base + col] = subtract_modulo(residue, product, self.modulus);
            }
            if let Some(rhs) = &mut self.rhs {
                let product = multiply_modulo(factor, rhs[pivot], self.modulus);
                rhs[row] = subtract_modulo(rhs[row], product, self.modulus);
            }
        }
        // Clearing the pivot row by column operations would leave the
        // trailing block unchanged, since the pivot column below it is zero.
        Ok(inverse)
    }
}

fn no_modular_solution() -> Error {
    Error::plain("BP integral witness congruence has no solution")
}

fn subtract_modulo(a: u64, b: u64, modulus: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        modulus - (b - a)
    }
}

fn multiply_modulo(a: u64, b: u64, modulus: u64) -> u64 {
    (u128::from(a) * u128::from(b) % u128::from(modulus)) as u64
}

fn inverse_modulo(model: &mut Model, a: u64, modulus: u64) -> Result<u64> {
    if a == 1 {
        return Ok(1);
    }
    let mut previous = i128::from(modulus);
    let mut current = i128::from(a);
    let mut x = 0;
    let mut y = 1;
    while current != 0 {
        model.spend(1)?;
        let quotient = previous / current;
        (previous, current) = (current, previous % current);
        (x, y) = (y, x - quotient * y);
    }
    // Remainders are u64 and all intermediate Bezout coefficients are
    // bounded by the modulus, so the signed arithmetic above fits i128.
    debug_assert_eq!(previous, 1);
    Ok(x.rem_euclid(i128::from(modulus)) as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primary_precision_expands_until_it_detects_all_pivots() {
        let mut model = Model::new();
        let matrix = Matrix {
            rows: 2,
            columns: 2,
            entries: vec![16, 16, 16, 48],
        };
        assert_eq!(
            primary_factors(&mut model, &matrix, 2, 2).unwrap_or_else(|error| panic!("{error}")),
            [16, 32]
        );
    }

    #[test]
    fn modular_preimage_tracks_columns_with_nonunit_pivots() {
        let mut model = Model::new();
        let matrix = Matrix {
            rows: 2,
            columns: 2,
            entries: vec![8, 4, 0, 12],
        };
        let witness = modular_preimage(&mut model, &matrix, 2, 16, &[4, 12], 1, 2)
            .unwrap_or_else(|error| panic!("{error}"));
        let image = matrix
            .apply(&mut model, &witness)
            .unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(image[0].rem_euclid(16), 4);
        assert_eq!(image[1].rem_euclid(16), 12);
    }
}
