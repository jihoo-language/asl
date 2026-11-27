use crate::error::Result;
use crate::integer_smith::{self, Track};
use crate::lee_linear::Rational as Q;
use crate::model::Model;
use crate::rv_qp_quad::Wide;

use super::linear::{self as lin, Matrix, ONE, ZERO};
use super::{arithmetic, certificate};

pub(super) struct Lattice {
    pub(super) gram: Matrix,
    pub(super) inverse: Matrix,
    pub(super) smith: Vec<i128>,
    pub(super) order: usize,
    quotient: Vec<i128>,
    lift: Vec<i128>,
    linking: Q,
}

impl Lattice {
    pub(super) fn new(model: &mut Model, widths: &[i128]) -> Result<Self> {
        let width = *widths
            .last()
            .ok_or_else(|| certificate("at least one geometric edge"))?;
        let d = widths.len() - 1;
        let gram: Matrix = widths[..d]
            .iter()
            .enumerate()
            .map(|(i, &w)| {
                (0..d)
                    .map(|j| Q::integer(width + if i == j { w } else { 0 }))
                    .collect()
            })
            .collect();
        let (inverse, determinant) = lin::inverse(model, &gram)?
            .ok_or_else(|| certificate("monodromy Gram is nonsingular over Q"))?;
        let mut integral = integer_smith::Matrix::new(model, d, d)?;
        for (i, row) in gram.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                integral.write(model, i, j, c.numerator)?;
            }
        }
        let frame = integer_smith::reduce(model, integral, Track::Both, Track::None)?;
        let order = frame.factors.iter().product::<i128>();
        let weight_product = widths.iter().product::<i128>();
        let expected = widths.iter().map(|w| weight_product / w).sum::<i128>();
        if order != expected
            || determinant != Q::integer(order)
            || frame.factors.len() != d
            || frame
                .factors
                .iter()
                .take(d.saturating_sub(1))
                .any(|&c| c != 1)
        {
            return Err(certificate(
                "cyclic Smith cokernel and weighted graph determinant",
            ));
        }
        let order = usize::try_from(order).map_err(|_| arithmetic())?;
        let mut quotient = vec![0; d];
        let mut lift = vec![0; d];
        if order > 1 {
            for (&i, &c) in frame.rows.coordinate(d - 1) {
                quotient[i] = c;
            }
            for (&i, &c) in frame.rows.basis(d - 1) {
                lift[i] = c;
            }
        }
        let image = lin::apply(
            model,
            &inverse,
            &lift.iter().copied().map(Q::integer).collect::<Vec<_>>(),
        )?;
        let linking = image
            .iter()
            .zip(&lift)
            .try_fold(ZERO, |sum, (&c, &a)| sum.add(c.multiply(Q::integer(a))?))?;
        let linking = Q::new(
            linking.numerator.rem_euclid(linking.denominator),
            linking.denominator,
        )?;
        if linking.denominator != order as i128 {
            return Err(certificate("nondegenerate component linking form"));
        }
        let lattice = Self {
            gram,
            inverse,
            smith: frame.factors,
            order,
            quotient,
            lift,
            linking,
        };
        for j in 0..d {
            if lattice.class(
                &lattice
                    .gram
                    .iter()
                    .map(|row| row[j].numerator)
                    .collect::<Vec<_>>(),
            ) != 0
            {
                return Err(certificate("component classes annihilate the Gram image"));
            }
        }
        if order > 1 && lattice.class(&lattice.lift) != 1 {
            return Err(certificate(
                "the Smith generator lifts to an actual dual vector",
            ));
        }
        Ok(lattice)
    }

    pub(super) fn class(&self, value: &[i128]) -> usize {
        let d = self.order as i128;
        let residue = self.quotient.iter().zip(value).fold(0, |sum, (&a, &b)| {
            (sum + a.rem_euclid(d) * b.rem_euclid(d)).rem_euclid(d)
        });
        residue as usize
    }

    pub(super) fn lift(&self, residue: usize) -> Result<Vec<i128>> {
        self.lift
            .iter()
            .map(|a| a.checked_mul(residue as i128).ok_or_else(arithmetic))
            .collect()
    }

    pub(super) fn link(&self, x: usize, y: usize) -> Result<Q> {
        let a = self.linking;
        Q::new(
            (a.numerator * (x * y) as i128).rem_euclid(a.denominator),
            a.denominator,
        )
    }
}

pub(super) fn permutation(images: &[usize]) -> Matrix {
    let n = images.len();
    let mut matrix = vec![vec![ZERO; n]; n];
    for (i, &j) in images.iter().enumerate() {
        matrix[j][i] = ONE;
    }
    matrix
}

pub(super) fn restrict(matrix: &Matrix) -> Result<Matrix> {
    let n = matrix.len();
    let d = n - 1;
    let mut output = vec![vec![ZERO; d]; d];
    for j in 0..d {
        let mut sum = ZERO;
        for (i, row) in matrix.iter().enumerate() {
            let image = lin::sub(row[j], row[d])?;
            if i < d {
                output[i][j] = image;
            }
            sum = sum.add(image)?;
        }
        if sum != ZERO {
            return Err(certificate("edge action preserves the cycle lattice"));
        }
    }
    Ok(output)
}

pub(super) fn transpose(matrix: &Matrix) -> Matrix {
    let width = matrix.first().map_or(0, Vec::len);
    (0..width)
        .map(|j| matrix.iter().map(|row| row[j]).collect())
        .collect()
}

pub(super) fn negative(matrix: &Matrix) -> Result<Matrix> {
    matrix
        .iter()
        .map(|row| row.iter().copied().map(lin::neg).collect())
        .collect()
}

pub(super) fn integral_apply(
    model: &mut Model,
    matrix: &Matrix,
    vector: &[i128],
) -> Result<Vec<i128>> {
    model.spend(9 * matrix.len() * vector.len())?;
    matrix
        .iter()
        .map(|row| {
            let sum =
                row.iter()
                    .zip(vector)
                    .try_fold(Wide::<9>::from_i128(0), |sum, (a, &b)| {
                        if a.denominator != 1 {
                            return Err(certificate("dual and cycle operators are integral"));
                        }
                        sum.add(Wide::product(a.numerator, b))
                            .map_err(|_| arithmetic())
                    })?;
            sum.signed().ok_or_else(arithmetic)
        })
        .collect()
}

pub(super) fn intertwine(
    model: &mut Model,
    left: &Matrix,
    right: &Matrix,
    map: &Matrix,
) -> Result<()> {
    if lin::multiply(model, map, left)? != lin::multiply(model, right, map)? {
        return Err(certificate(
            "exact action commutes with the specified comparison map",
        ));
    }
    Ok(())
}
