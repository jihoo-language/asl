use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_cyclotomic;
use crate::rv_lens::{self, Rational};

const MAX_MATRIX_CELLS: usize = 65_536;

pub(super) struct Field<'a> {
    pub(super) model: &'a mut Model,
    pub(super) conductor: usize,
    roots: Vec<Option<Rational>>,
}

impl<'a> Field<'a> {
    pub(super) fn new(model: &'a mut Model, conductor: usize) -> Self {
        Self {
            model,
            conductor,
            roots: vec![None; conductor],
        }
    }

    pub(super) fn phase(&mut self, exponent: usize) -> Result<Rational> {
        let i = exponent % self.conductor;
        if let Some(root) = &self.roots[i] {
            return Ok(root.clone());
        }
        let mut polynomial = vec![0; i + 1];
        polynomial[i] = 1;
        let polynomial = rv_cyclotomic::reduced(self.model, self.conductor, polynomial)?;
        let root = rv_lens::fraction(self.model, &polynomial, 1)?;
        self.roots[i] = Some(root.clone());
        Ok(root)
    }

    pub(super) fn add(&mut self, left: &Rational, right: &Rational) -> Result<Rational> {
        rv_lens::sum(self.model, left, right)
    }

    pub(super) fn multiply(&mut self, left: &Rational, right: &Rational) -> Result<Rational> {
        rv_lens::multiply(self.model, self.conductor, left, right)
    }

    pub(super) fn conjugate(&mut self, value: &Rational) -> Result<Rational> {
        if value.numerator.len() <= 1 {
            return Ok(value.clone());
        }
        let polynomial = rv_cyclotomic::substitute_root(
            self.model,
            self.conductor,
            &value.numerator,
            self.conductor - 1,
        )?;
        rv_lens::fraction(self.model, &polynomial, value.denominator)
    }

    pub(super) fn shift(&mut self, value: &Rational, phase: usize) -> Result<Rational> {
        rv_lens::shifted(self.model, self.conductor, value, phase as i128)
    }
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(super) struct Matrix {
    pub(super) rows: usize,
    pub(super) cols: usize,
    values: Vec<Rational>,
}

impl Matrix {
    pub(super) fn zero(rows: usize, cols: usize) -> Result<Self> {
        let cells = rows
            .checked_mul(cols)
            .filter(|&n| n <= MAX_MATRIX_CELLS)
            .ok_or_else(|| Error::plain("double fusion matrix exceeds its storage limit"))?;
        Ok(Self {
            rows,
            cols,
            values: vec![Rational::zero(); cells],
        })
    }

    pub(super) fn identity(size: usize) -> Result<Self> {
        let mut matrix = Self::zero(size, size)?;
        for i in 0..size {
            matrix.set(i, i, Rational::integer(1));
        }
        Ok(matrix)
    }

    pub(super) fn cells(&self) -> usize {
        3 + self
            .values
            .iter()
            .map(|v| 3 + 2 * v.numerator.len())
            .sum::<usize>()
    }

    pub(super) fn get(&self, row: usize, col: usize) -> &Rational {
        &self.values[row * self.cols + col]
    }

    pub(super) fn set(&mut self, row: usize, col: usize, value: Rational) {
        self.values[row * self.cols + col] = value;
    }

    pub(super) fn entry(&self, [row, col]: [usize; 2]) -> Result<Rational> {
        if row >= self.rows || col >= self.cols {
            return Err(Error::plain("double fusion matrix index is out of range"));
        }
        Ok(self.get(row, col).clone())
    }

    pub(super) fn multiply(&self, field: &mut Field<'_>, right: &Self) -> Result<Self> {
        debug_assert_eq!(self.cols, right.rows);
        let mut out = Self::zero(self.rows, right.cols)?;
        field.model.spend(self.rows * self.cols * right.cols)?;
        for i in 0..self.rows {
            for k in 0..self.cols {
                let x = self.get(i, k);
                if x.is_zero() {
                    continue;
                }
                for j in 0..right.cols {
                    let y = right.get(k, j);
                    if y.is_zero() {
                        continue;
                    }
                    let term = field.multiply(x, y)?;
                    out.set(i, j, field.add(out.get(i, j), &term)?);
                }
            }
        }
        Ok(out)
    }

    pub(super) fn tensor(&self, field: &mut Field<'_>, right: &Self) -> Result<Self> {
        let mut out = Self::zero(self.rows * right.rows, self.cols * right.cols)?;
        field.model.spend(out.values.len())?;
        for i in 0..self.rows {
            for j in 0..self.cols {
                let x = self.get(i, j);
                if x.is_zero() {
                    continue;
                }
                for k in 0..right.rows {
                    for l in 0..right.cols {
                        let y = right.get(k, l);
                        if y.is_zero() {
                            continue;
                        }
                        out.set(
                            i * right.rows + k,
                            j * right.cols + l,
                            field.multiply(x, y)?,
                        );
                    }
                }
            }
        }
        Ok(out)
    }

    pub(super) fn transformed(
        &self,
        field: &mut Field<'_>,
        image: impl Fn(usize) -> (usize, usize),
    ) -> Result<Self> {
        field.model.spend(self.values.len())?;
        let mut out = Self::zero(self.rows, self.cols)?;
        for i in 0..self.rows {
            let (target, phase) = image(i);
            for j in 0..self.cols {
                out.set(target, j, field.shift(self.get(i, j), phase)?);
            }
        }
        Ok(out)
    }

    pub(super) fn inner(&self, field: &mut Field<'_>, other: &Self) -> Result<Rational> {
        debug_assert_eq!((self.rows, self.cols), (other.rows, other.cols));
        field.model.spend(self.values.len())?;
        let mut result = Rational::zero();
        for (x, y) in self.values.iter().zip(&other.values) {
            if x.is_zero() || y.is_zero() {
                continue;
            }
            let x = field.conjugate(x)?;
            let term = field.multiply(&x, y)?;
            result = field.add(&result, &term)?;
        }
        Ok(result)
    }
}

pub(super) fn coordinates(
    field: &mut Field<'_>,
    target: &[&Matrix],
    columns: &[&Matrix],
) -> Result<Matrix> {
    if target.len() != columns.len() {
        return Err(Error::plain("double fusion bases have unequal dimension"));
    }
    let mut matrix = Matrix::zero(target.len(), columns.len())?;
    for (i, right) in target.iter().enumerate() {
        for &other in &target[..i] {
            if !right.inner(field, other)?.is_zero() {
                return Err(Error::plain("double fusion target basis is not orthogonal"));
            }
        }
        let norm = right.inner(field, right)?;
        let &[norm_top] = norm.numerator.as_ref() else {
            return Err(Error::plain(
                "double fusion basis has nonrational or zero norm",
            ));
        };
        if norm_top <= 0 {
            return Err(Error::plain("double fusion basis norm is not positive"));
        }
        let norm_bottom = i128::try_from(norm.denominator).map_err(|_| {
            Error::plain("cyclotomic norm denominator exceeds signed 128-bit range")
        })?;
        let reciprocal = rv_lens::fraction(field.model, &[norm_bottom], norm_top as u128)?;
        for (j, col) in columns.iter().enumerate() {
            let numerator = right.inner(field, col)?;
            matrix.set(i, j, field.multiply(&numerator, &reciprocal)?);
        }
    }
    // Check the solved equation in the actual tensor coordinates.
    for (j, col) in columns.iter().enumerate() {
        field.model.spend(col.values.len() * target.len())?;
        for (position, desired) in col.values.iter().enumerate() {
            let mut actual = Rational::zero();
            for (i, right) in target.iter().enumerate() {
                let u = &right.values[position];
                let v = matrix.get(i, j);
                if u.is_zero() || v.is_zero() {
                    continue;
                }
                let term = field.multiply(u, v)?;
                actual = field.add(&actual, &term)?;
            }
            if actual != *desired {
                return Err(Error::plain(
                    "double fusion coordinates do not reconstruct the intertwiner",
                ));
            }
        }
    }
    Ok(matrix)
}
