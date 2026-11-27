use std::cell::RefCell;

use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::rv_lens::{self, Rational as Scalar};

use super::field::Field;
use super::fusion;
use super::{check_report, Defect, Error, Model, Rc, Result, Value};

mod split;
use split::Split;

type Square = Matrix<Scalar>;

pub(crate) struct Convolution {
    pub(crate) d: usize,
    mark: Rc<str>,
    supports: [Vec<usize>; 2],
    intervals: Option<[(usize, usize); 2]>,
    weights: Vec<usize>,
    powers: Vec<Square>,
    split: RefCell<Option<Rc<Split>>>,
}

impl Convolution {
    pub(crate) fn new(model: &mut Model, left: &Defect, right: &Defect) -> Result<Self> {
        if left.mark != right.mark || left.d != right.d {
            return Err(Error::plain(
                "LG convolution requires the same Tate marking",
            ));
        }
        let field = Field::new(model, left.d)?;
        let (r, s) = (left.indices.len(), right.indices.len());
        let n = r * s;
        let cells = n
            .saturating_mul(n)
            .saturating_mul(left.d + 3)
            .saturating_mul(field.cells());
        if cells > crate::complex::MAX_PIVOT_CELLS {
            return Err(Error::plain(
                "LG convolution exceeds the coefficient storage limit",
            ));
        }
        let weights: Vec<_> = (0..n).map(|i| i / s + i % s).collect();
        let x = action(model, &field, left, right)?;
        check_weights(&x, &weights, &weights, 1)?;
        let identity = hw_algebra::identity(model, &field, n)?;
        let mut powers = Vec::with_capacity(left.d);
        powers.push(identity);
        powers.push(x);
        for _ in 2..left.d {
            powers.push(hw_algebra::multiply(
                model,
                &field,
                powers.last().unwrap(),
                &powers[1],
            )?);
        }
        if hw_algebra::multiply(model, &field, powers.last().unwrap(), &powers[1])? != powers[0] {
            return Err(Error::plain(
                "LG convolution failed the relation X(z)^d=z^d I",
            ));
        }
        // The two displayed factorization matrices telescope to x^d I-X(z)^d.
        model.allocate(cells + n + r + s + 12)?;
        let intervals = fusion::interval(left)
            .zip(fusion::interval(right))
            .map(|(a, b)| [a, b]);
        Ok(Self {
            d: left.d,
            mark: left.mark.clone(),
            supports: [left.indices.clone(), right.indices.clone()],
            intervals,
            weights,
            powers,
            split: RefCell::new(None),
        })
    }

    pub(crate) fn report(&self) -> Result<String> {
        let report = format!(
            "W=x^{}−z^{} supports={:?} rank={} weights={:?} power=verified factorization=verified",
            self.d,
            self.d,
            self.supports,
            self.weights.len(),
            self.weights
        );
        check_report(&report)?;
        Ok(report)
    }

    fn split(&self, model: &mut Model) -> Result<Rc<Split>> {
        if let Some(split) = self.split.borrow().as_ref() {
            return Ok(split.clone());
        }
        let split = Rc::new(Split::compute(model, self)?);
        *self.split.borrow_mut() = Some(split.clone());
        Ok(split)
    }

    pub(crate) fn split_report(&self, model: &mut Model) -> Result<String> {
        self.split(model)?.report(self.d)
    }

    pub(crate) fn part(&self, model: &mut Model, i: usize) -> Result<Defect> {
        let split = self.split(model)?;
        let part = split
            .parts
            .get(i)
            .ok_or_else(|| Error::plain("LG convolution summand index is out of range"))?;
        if part.support.len() == self.d {
            return Err(Error::plain(
                "LG convolution summand is a free block with full root support",
            ));
        }
        Defect::new(model, self.mark.clone(), self.d, part.support.clone())
    }

    pub(crate) fn term(
        &self,
        model: &mut Model,
        [kind, row, col, x, z]: [usize; 5],
    ) -> Result<Value> {
        if row >= self.weights.len() || col >= self.weights.len() {
            return Err(Error::plain("LG convolution matrix index is out of range"));
        }
        let scalar = match kind {
            0 if x == 0 => homogeneous(
                &self.powers[1],
                &self.weights,
                &self.weights,
                1,
                row,
                col,
                z,
            ),
            1 if x == 1 && z == 0 && row == col => Scalar::integer(1),
            1 if x == 0 => rv_lens::negative(
                model,
                &homogeneous(
                    &self.powers[1],
                    &self.weights,
                    &self.weights,
                    1,
                    row,
                    col,
                    z,
                ),
            )?,
            2 if x < self.d => {
                let k = self.d - 1 - x;
                homogeneous(
                    &self.powers[k],
                    &self.weights,
                    &self.weights,
                    k,
                    row,
                    col,
                    z,
                )
            }
            3..=5 => {
                let split = self.split(model)?;
                match kind {
                    _ if x != 0 => Scalar::zero(),
                    3 => homogeneous(
                        &split.forward,
                        &self.weights,
                        &split.weights,
                        0,
                        row,
                        col,
                        z,
                    ),
                    4 => homogeneous(
                        &split.inverse,
                        &split.weights,
                        &self.weights,
                        0,
                        row,
                        col,
                        z,
                    ),
                    _ => homogeneous(
                        &split.action,
                        &split.weights,
                        &split.weights,
                        1,
                        row,
                        col,
                        z,
                    ),
                }
            }
            0..=2 => Scalar::zero(),
            _ => {
                return Err(Error::plain(
                    "LG convolution matrix kind must be between zero and five",
                ))
            }
        };
        rv_lens::rational_value(model, scalar)
    }
}

// Specialization z=1 records the scalar in each homogeneous matrix entry.
// Its unique z exponent is column degree + operator degree - row degree.
fn homogeneous(
    matrix: &Square,
    rows: &[usize],
    cols: &[usize],
    degree: usize,
    row: usize,
    col: usize,
    power: usize,
) -> Scalar {
    if (cols[col] + degree).checked_sub(rows[row]) == Some(power) {
        matrix.entry(row, col).clone()
    } else {
        Scalar::zero()
    }
}

fn check_weights(matrix: &Square, rows: &[usize], cols: &[usize], degree: usize) -> Result<()> {
    for (row, &r) in rows.iter().enumerate() {
        for (col, &c) in cols.iter().enumerate() {
            if r > c + degree && !matrix.entry(row, col).is_zero() {
                return Err(Error::plain(
                    "LG convolution matrix has a negative homogeneous exponent",
                ));
            }
        }
    }
    Ok(())
}

fn action(model: &mut Model, field: &Field, left: &Defect, right: &Defect) -> Result<Square> {
    let (r, s) = (left.indices.len(), right.indices.len());
    let n = r * s;
    let mut entries = vec![Scalar::zero(); n * n];
    let px: Vec<_> = left.p[..r]
        .iter()
        .map(|c| rv_lens::fraction(model, c, 1))
        .collect::<Result<_>>()?;
    let py: Vec<_> = right.p[..s]
        .iter()
        .map(|c| rv_lens::fraction(model, c, 1))
        .collect::<Result<_>>()?;
    for a in 0..r - 1 {
        for b in 0..s {
            let column = a * s + b;
            entries[((a + 1) * s + b) * n + column] = Scalar::integer(1);
        }
    }
    for b in 0..s {
        let column = (r - 1) * s + b;
        let mut image = vec![Scalar::zero(); s];
        image[b] = Scalar::integer(1);
        for (a, c) in px.iter().enumerate().rev() {
            // Multiplication by y, reducing its top power by p_J(y,1).
            let top = image.pop().unwrap();
            image.insert(0, Scalar::zero());
            model.spend(s)?;
            for (h, (value, p)) in image.iter_mut().zip(&py).enumerate() {
                let reduction = field.multiply(model, &top, p)?;
                *value = field.subtract(model, value, &reduction)?;
                let product = field.multiply(model, value, c)?;
                entries[(a * s + h) * n + column] = rv_lens::negative(model, &product)?;
            }
        }
    }
    Matrix::from_entries(n, entries)
}
