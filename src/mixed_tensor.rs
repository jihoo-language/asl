use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::mixed::{Complex, Differential, Id};
use crate::mixed_bars::{Bar, Torsion};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::{add_scaled, Coordinate};

#[derive(Clone, Copy)]
pub(crate) enum Certificate {
    Product,
    ProductNull,
    Secondary,
    SecondaryNull,
}

pub(crate) struct Tensor {
    pub(crate) left: Id,
    pub(crate) right: Id,
    width: usize,
    prime: Prime,
}

impl Tensor {
    pub(crate) fn build(model: &mut Model, left: Id, right: Id) -> Result<(Self, Complex)> {
        let a = model.mixed_complex(left);
        let b = model.mixed_complex(right);
        if a.prime != b.prime {
            return Err(Error::plain(
                "mixtensor factors must have the same prime field",
            ));
        }
        let tensor = Self {
            left,
            right,
            width: b.len(),
            prime: a.prime,
        };
        let dimension = a.len().saturating_mul(b.len());
        check_storage(6usize.saturating_mul(dimension))?;
        model.spend(dimension)?;
        let mut degrees = Vec::with_capacity(dimension);
        for i in 0..a.len() {
            for j in 0..b.len() {
                degrees.push(degree(
                    i128::from(a.raw_degree(i)) + i128::from(b.raw_degree(j)),
                )?);
            }
        }
        let mut cells = 6 * dimension;
        let down = tensor.columns(model, &a, &b, Differential::Down, &mut cells)?;
        let up = tensor.columns(model, &a, &b, Differential::Up, &mut cells)?;
        let complex = Complex::from_columns(model, tensor.prime, degrees, down, up)?;
        Ok((tensor, complex))
    }

    fn columns(
        &self,
        model: &mut Model,
        left: &Complex,
        right: &Complex,
        operation: Differential,
        cells: &mut usize,
    ) -> Result<Vec<Vec<Coordinate>>> {
        let mut columns = Vec::with_capacity(left.len() * right.len());
        for i in 0..left.len() {
            let a = left.column(operation, i);
            let odd = left.raw_degree(i) % 2 != 0;
            for j in 0..right.len() {
                let b = right.column(operation, j);
                let count = a.len() + b.len();
                *cells += 2 * count;
                check_storage(*cells)?;
                model.spend(1 + count * (1 + count.max(1).ilog2() as usize))?;
                let mut image = Vec::with_capacity(count);
                image.extend(a.iter().map(|&(k, c)| (self.index(k, j), c)));
                image.extend(
                    b.iter().map(|&(k, c)| {
                        (self.index(i, k), if odd { self.prime.negate(c) } else { c })
                    }),
                );
                image.sort_unstable();
                columns.push(image);
            }
        }
        Ok(columns)
    }

    fn index(&self, left: usize, right: usize) -> usize {
        left * self.width + right
    }

    pub(crate) fn product(
        &self,
        model: &mut Model,
        left: &[Coordinate],
        right: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let count = left.len().saturating_mul(right.len());
        check_storage(2usize.saturating_mul(count))?;
        model.spend(count)?;
        let mut result = Vec::with_capacity(count);
        for &(i, a) in left {
            for &(j, b) in right {
                result.push((self.index(i, j), self.prime.mul(a, b)));
            }
        }
        Ok(result)
    }

    pub(crate) fn transpose(
        &self,
        model: &mut Model,
        target: &Self,
        terms: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let left = model.mixed_complex(self.left);
        let right = model.mixed_complex(self.right);
        model.spend(terms.len() * (1 + terms.len().max(1).ilog2() as usize))?;
        let mut result = Vec::with_capacity(terms.len());
        for &(index, coefficient) in terms {
            let (i, j) = (index / self.width, index % self.width);
            let odd = left.raw_degree(i) % 2 != 0 && right.raw_degree(j) % 2 != 0;
            result.push((
                target.index(j, i),
                if odd {
                    self.prime.negate(coefficient)
                } else {
                    coefficient
                },
            ));
        }
        result.sort_unstable();
        Ok(result)
    }

    pub(crate) fn certificate(
        &self,
        model: &mut Model,
        kind: Certificate,
        left: &Bar,
        right: &Bar,
    ) -> Result<(i64, Vec<Coordinate>)> {
        let birth = i128::from(left.birth) + i128::from(right.birth);
        let a = left.torsion.as_ref();
        let b = right.torsion.as_ref();
        let odd = left.birth % 2 != 0;
        let sign = if odd { self.prime.negate(1) } else { 1 };
        match kind {
            Certificate::Product => {
                let support = self.product(model, &left.cycle, &right.cycle)?;
                Ok((degree(birth)?, support))
            }
            Certificate::ProductNull => {
                let result = match (a, b) {
                    (Some(a), Some(b)) if a.length > b.length => {
                        self.right_null(model, left, b, sign)?
                    }
                    (Some(a), _) => {
                        let terms = self.product(model, &a.witness, &right.cycle)?;
                        (i128::from(a.null_degree) + i128::from(right.birth), terms)
                    }
                    (None, Some(b)) => self.right_null(model, left, b, sign)?,
                    (None, None) => {
                        return Err(Error::plain(
                            "mixTensorNull requires at least one torsion bar",
                        ))
                    }
                };
                Ok((degree(result.0)?, result.1))
            }
            Certificate::Secondary | Certificate::SecondaryNull => {
                let (Some(a), Some(b)) = (a, b) else {
                    return Err(Error::plain(
                        "the secondary tensor certificate requires two torsion bars",
                    ));
                };
                if matches!(kind, Certificate::SecondaryNull) {
                    let total = i128::from(a.null_degree) + i128::from(b.null_degree);
                    let support = self.product(model, &a.witness, &b.witness)?;
                    return Ok((degree(total)?, support));
                }

                // The result degree incorporates the nonnegative powers
                // l-min(l,m) and m-min(l,m) of the two summands.
                let degree = degree(birth - 2 * i128::from(a.length.max(b.length)) + 1)?;
                let first = self.product(model, &left.cycle, &b.witness)?;
                let second = self.product(model, &a.witness, &right.cycle)?;
                let support =
                    add_scaled(model, self.prime, &first, &second, self.prime.negate(sign))?;
                Ok((degree, support))
            }
        }
    }

    fn right_null(
        &self,
        model: &mut Model,
        left: &Bar,
        right: &Torsion,
        sign: u32,
    ) -> Result<(i128, Vec<Coordinate>)> {
        let mut support = self.product(model, &left.cycle, &right.witness)?;
        if sign != 1 {
            model.spend(support.len())?;
            for (_, coefficient) in &mut support {
                *coefficient = self.prime.mul(*coefficient, sign);
            }
        }
        Ok((
            i128::from(left.birth) + i128::from(right.null_degree),
            support,
        ))
    }
}

fn degree(value: i128) -> Result<i64> {
    i64::try_from(value).map_err(|_| Error::plain("mixed tensor degree is out of range"))
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "mixed tensor complex exceeds the storage limit",
        ));
    }
    Ok(())
}
