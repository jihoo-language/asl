use std::mem::size_of;

use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::model::Model;
use crate::wide_integer::Natural;

use super::norm::{half, power};
use super::{
    add, determinant_from_ghosts, exact_divide, integrity, multiply, one, power_sums, require_rank,
    trace_at_product, zero, Integer, Spectrum,
};

mod major;
pub(crate) use major::Branch;

#[derive(Clone)]
pub(crate) struct Shape {
    rows: Box<[u128]>,
    weight: u128,
}

impl Shape {
    pub(crate) fn new(rows: Vec<u128>) -> Result<Self> {
        if rows.last() == Some(&0) || rows.windows(2).any(|pair| pair[0] < pair[1]) {
            return Err(Error::plain(
                "Schur shape rows must be positive and weakly decreasing",
            ));
        }
        let weight = rows
            .iter()
            .try_fold(0u128, |weight, &row| weight.checked_add(row))
            .ok_or_else(|| Error::plain("Schur shape weight exceeds unsigned 128-bit range"))?;
        Ok(Self {
            rows: rows.into_boxed_slice(),
            weight,
        })
    }

    pub(crate) fn rows(&self) -> &[u128] {
        &self.rows
    }

    fn elementary_major(&self) -> Option<u128> {
        if self.rows.len() <= 1 {
            Some(0)
        } else if self.rows[0] == 1 {
            Some(if self.weight.is_multiple_of(2) {
                self.weight / 2
            } else {
                0
            })
        } else {
            None
        }
    }

    fn bounded_weight(&self) -> Result<usize> {
        usize::try_from(self.weight)
            .ok()
            .filter(|&n| n <= 256)
            .ok_or_else(|| Error::plain("Schur shape exceeds the weight limit of 256 boxes"))
    }
}

impl Spectrum {
    pub(crate) fn schur(&self, model: &mut Model, shape: &Shape) -> Result<Self> {
        if shape.weight == 0 {
            return Self::line(1);
        }
        if !self.in_hook(shape) {
            return Ok(Self::zero());
        }
        if shape.weight == 1 {
            return Ok(self.clone());
        }
        if let Some((value, odd)) = self.schur_line(model, shape.weight)? {
            let spectrum = Self::effective(vec![one(), value.negate()])?;
            return Ok(if odd { spectrum.negate() } else { spectrum });
        }
        let plan = JacobiTrudi::new(model, shape)?;
        let a = Integer::from_i128((self.numerator.len() - 1) as i128);
        let b = Integer::from_i128((self.denominator.len() - 1) as i128);
        let ranks = plan.parities(model, |_| (a, b))?;
        let (even, odd) = (require_rank(ranks[0])?, require_rank(ranks[1])?);
        let rank = even.max(odd);
        let limit = rank * plan.limit;
        hw_algebra::storage(2 * (limit + 1), IntegerRing.cells()).map_err(|_| {
            Error::plain("Schur spectrum power sums exceed the coefficient storage limit")
        })?;
        let a = power_sums(model, &self.numerator, limit)?;
        let b = power_sums(model, &self.denominator, limit)?;
        model.spend(even + odd)?;
        let mut p = Vec::with_capacity(even);
        let mut q = Vec::with_capacity(odd);
        for k in 1..=rank {
            let [positive, negative] = plan.parities(model, |i| (a[k * i], b[k * i]))?;
            if k <= even {
                p.push(positive);
            }
            if k <= odd {
                q.push(negative);
            }
        }
        let p = determinant_from_ghosts(model, &p)?;
        let q = determinant_from_ghosts(model, &q)?;
        if p[even].is_zero() || q[odd].is_zero() {
            return Err(integrity());
        }
        Self::canonical(model, p, q)
    }

    pub(super) fn schur_ghost(&self, model: &mut Model, shape: &Shape, k: u128) -> Result<Integer> {
        if k == 0 {
            return Err(Error::plain("Schur spectrum ghost index must be positive"));
        }
        if shape.weight == 0 {
            return Ok(one());
        }
        if !self.in_hook(shape) {
            return Ok(zero());
        }
        if shape.weight == 1 {
            return self.ghost(model, k);
        }
        if let Some((value, odd)) = self.schur_line(model, shape.weight)? {
            let value = power(model, value, k)?;
            return Ok(if odd { value.negate() } else { value });
        }
        let plan = JacobiTrudi::new(model, shape)?;
        plan.trace(model, |model, i| {
            let a = trace_at_product(model, &self.numerator, k, i as u128)?;
            let b = trace_at_product(model, &self.denominator, k, i as u128)?;
            add(a, b.negate())
        })
    }

    fn in_hook(&self, shape: &Shape) -> bool {
        // The Schur functor on an (a|b)-dimensional supermodule vanishes
        // outside the hook λ_(a+1) <= b. Check before bounding the weight.
        let a = self.numerator.len() - 1;
        let b = (self.denominator.len() - 1) as u128;
        shape.rows.get(a).is_none_or(|&row| row <= b)
    }

    fn schur_line(&self, model: &mut Model, weight: u128) -> Result<Option<(Integer, bool)>> {
        let (coefficient, odd) = match (&*self.numerator, &*self.denominator) {
            ([1, coefficient], [1]) => (coefficient, false),
            ([1], [1, coefficient]) => (coefficient, !weight.is_multiple_of(2)),
            _ => return Ok(None),
        };
        let value = power(model, Integer::from_i128(*coefficient).negate(), weight)?;
        Ok(Some((value, odd)))
    }
}

struct JacobiTrudi {
    rows: Vec<usize>,
    elementary: bool,
    limit: usize,
}

impl JacobiTrudi {
    fn new(model: &mut Model, shape: &Shape) -> Result<Self> {
        let weight = shape.bounded_weight()?;
        let width = shape.rows[0] as usize;
        let height = shape.rows.len();
        let elementary = width < height;
        model.spend(weight)?;
        let rows: Vec<_> = if elementary {
            let mut column_height = height;
            (1..=width)
                .map(|column| {
                    while shape.rows[column_height - 1] < column as u128 {
                        column_height -= 1;
                    }
                    column_height
                })
                .collect()
        } else {
            shape.rows.iter().map(|&row| row as usize).collect()
        };
        let order = rows.len();
        hw_algebra::storage(order * order, IntegerRing.cells())
            .map_err(|_| Error::plain("Schur determinant exceeds the coefficient storage limit"))?;
        let limit = rows[0] + order - 1;
        Ok(Self {
            rows,
            elementary,
            limit,
        })
    }

    fn parities(
        &self,
        model: &mut Model,
        parts: impl Fn(usize) -> (Integer, Integer),
    ) -> Result<[Integer; 2]> {
        let trace = |model: &mut Model, signed| {
            self.trace(model, |_, i| {
                let (a, b) = parts(i);
                if signed || i.is_multiple_of(2) {
                    add(a, b.negate())
                } else {
                    add(a, b)
                }
            })
        };
        let supertrace = trace(model, true)?;
        let trace = trace(model, false)?;
        Ok([
            half(add(trace, supertrace)?)?,
            half(add(trace, supertrace.negate())?)?,
        ])
    }

    fn trace(
        &self,
        model: &mut Model,
        mut ghost: impl FnMut(&mut Model, usize) -> Result<Integer>,
    ) -> Result<Integer> {
        let mut p = Vec::with_capacity(self.limit);
        for i in 1..=self.limit {
            let value = ghost(model, i)?;
            p.push(if self.elementary && i.is_multiple_of(2) {
                value.negate()
            } else {
                value
            });
        }
        let mut complete = vec![one()];
        for i in 1..=self.limit {
            let mut sum = zero();
            for (&a, &b) in complete.iter().rev().zip(&p) {
                sum = add(sum, multiply(model, a, b)?)?;
            }
            complete.push(exact_divide(model, sum, Natural::from_small(i as u128))?);
        }
        let order = self.rows.len();
        model.spend(order * order)?;
        let mut entries = Vec::with_capacity(order * order);
        for (i, &row) in self.rows.iter().enumerate() {
            for j in 0..order {
                entries.push(Scalar(
                    (row + j).checked_sub(i).map_or_else(zero, |d| complete[d]),
                ));
            }
        }
        let matrix = Matrix::from_entries(order, entries)?;
        let Scalar(value) = hw_algebra::characteristic(model, &IntegerRing, &matrix)?
            .pop()
            .unwrap();
        Ok(if order.is_multiple_of(2) {
            value
        } else {
            value.negate()
        })
    }
}

#[derive(Clone)]
struct Scalar(Integer);

impl PartialEq for Scalar {
    fn eq(&self, rhs: &Self) -> bool {
        self.0.magnitude == rhs.0.magnitude
            && (self.0.is_zero() || self.0.negative == rhs.0.negative)
    }
}

impl Eq for Scalar {}

struct IntegerRing;

impl CoefficientRing for IntegerRing {
    type Element = Scalar;

    fn cells(&self) -> usize {
        size_of::<Integer>().div_ceil(size_of::<usize>())
    }
    fn zero(&self) -> Scalar {
        Scalar(zero())
    }
    fn one(&self) -> Scalar {
        Scalar(one())
    }
    fn is_zero(&self, value: &Scalar) -> bool {
        value.0.is_zero()
    }
    fn is_unit(&self, value: &Scalar) -> bool {
        value.0.magnitude.small() == Some(1)
    }
    fn add(&self, _: &mut Model, a: &Scalar, b: &Scalar) -> Result<Scalar> {
        Ok(Scalar(add(a.0, b.0)?))
    }
    fn subtract(&self, _: &mut Model, a: &Scalar, b: &Scalar) -> Result<Scalar> {
        Ok(Scalar(add(a.0, b.0.negate())?))
    }
    fn multiply(&self, model: &mut Model, a: &Scalar, b: &Scalar) -> Result<Scalar> {
        Ok(Scalar(multiply(model, a.0, b.0)?))
    }
    fn inverse_unit(&self, _: &mut Model, value: &Scalar) -> Result<Scalar> {
        if !self.is_unit(value) {
            return Err(Error::plain("Schur scalar has no integral inverse"));
        }
        Ok(value.clone())
    }
}
