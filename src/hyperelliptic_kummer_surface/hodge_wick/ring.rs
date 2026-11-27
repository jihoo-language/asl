use std::fmt::Write;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, MAX_SHOWN_BYTES};
use crate::syntax::{Label, ParsedExpr};

use super::moments::Moments;
use super::rational::{Rational, LIMBS};
use super::{multiply_powers, storage, work_limit, Form, Grade, Power, Tensor};

pub(crate) mod koszul;
mod linear;
use linear::Echelon;

pub(crate) enum Command {
    Rank(Label, i128, i128),
    Relations(Label, i128, i128),
    Signature(i128, Box<ParsedExpr>),
}

type Monomial = Vec<Power>;

pub(super) struct Slice {
    generators: Vec<Monomial>,
    // With i <= n, Verbitsky injectivity makes the reduction the identity.
    reduction: Option<Echelon>,
}

impl Slice {
    fn rank(&self) -> usize {
        self.reduction
            .as_ref()
            .map_or(self.generators.len(), |r| r.pivots.len())
    }

    fn pivots(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.rank()).map(|i| self.reduction.as_ref().map_or(i, |r| r.pivots[i]))
    }

    fn cells(&self) -> usize {
        let words: usize = self.generators.iter().map(|m| 8 + 5 * m.len()).sum();
        8 + words
            + self.reduction.as_ref().map_or(0, |r| {
                r.rows.len() * (8 + self.generators.len() * (2 * LIMBS + 1))
            })
    }

    fn projected(&self, model: &mut Model, input: &Tensor) -> Result<Tensor> {
        let Some(reduction) = &self.reduction else {
            return Ok(input.clone());
        };
        let mut output = Tensor::default();
        for (powers, &coefficient) in &input.terms {
            let column = self
                .generators
                .binary_search(powers)
                .map_err(|_| grading_error())?;
            if reduction.pivots.binary_search(&column).is_ok() {
                output.insert(model, powers.clone(), coefficient)?;
                continue;
            }
            for (&pivot, row) in reduction.pivots.iter().zip(&reduction.rows) {
                let value = coefficient.multiply(model, row[column])?;
                output.insert(model, self.generators[pivot].clone(), value)?;
            }
        }
        Ok(output)
    }

    fn expression(&self, model: &mut Model, values: &[Rational]) -> Result<Tensor> {
        let mut tensor = Tensor::default();
        for (index, &value) in self.pivots().zip(values) {
            tensor.insert(model, self.generators[index].clone(), value)?;
        }
        Ok(tensor)
    }
}

impl Form {
    fn slice(&self, model: &mut Model, grade: &Grade) -> Result<Rc<Slice>> {
        if let Some(slice) = self
            .slices
            .borrow()
            .get(&(grade.length, grade.degree))
            .cloned()
        {
            return Ok(slice);
        }
        let slice = if !grade.supports_states() {
            Slice {
                generators: Vec::new(),
                reduction: None,
            }
        } else {
            let n = grade.length as u128;
            let weight = grade.degree as u128 / 2;
            let generators = monomials(model, self, weight)?;
            let reduction = if weight <= n || generators.is_empty() {
                None
            } else {
                // The row space of the complementary-degree pairing is the
                // dual of the quotient by its radical. RREF pivot columns
                // identify monomials and give coordinates for every column.
                let dual_weight = 2 * n - weight;
                let rank = monomial_count(self.gram.size(), self.field.degree(), dual_weight)?;
                if rank == generators.len() {
                    None
                } else {
                    let dual = monomials(model, self, dual_weight)?;
                    let mut gram = linear::matrix(dual.len(), generators.len())?;
                    let mut moments = Moments::new(self);
                    for (row, complement) in gram.iter_mut().zip(&dual) {
                        for (value, powers) in row.iter_mut().zip(&generators) {
                            let product = multiply_powers(complement, powers)?;
                            *value = moments.monomial(model, &product)?;
                        }
                    }
                    let reduced = linear::rref(model, gram, generators.len())?;
                    if reduced.pivots.len() != dual.len() {
                        return Err(Error::plain(
                            "Verbitsky complementary pairing has unexpected rank",
                        ));
                    }
                    Some(reduced)
                }
            };
            Slice {
                generators,
                reduction,
            }
        };
        storage(slice.cells())?;
        model.allocate(slice.cells())?;
        let slice = Rc::new(slice);
        self.slices
            .borrow_mut()
            .insert((grade.length, grade.degree), slice.clone());
        Ok(slice)
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        grade: &Grade,
        input: &Tensor,
    ) -> Result<Tensor> {
        if !grade.supports_states() || input.terms.is_empty() {
            return Ok(Tensor::default());
        }
        if grade.degree as u128 / 2 <= grade.length as u128 {
            return Ok(input.clone());
        }
        self.slice(model, grade)?.projected(model, input)
    }

    pub(crate) fn basis(&self, model: &mut Model, grade: &Grade, index: usize) -> Result<Tensor> {
        let slice = self.slice(model, grade)?;
        let pivot = slice
            .pivots()
            .nth(index)
            .ok_or_else(|| Error::plain("Verbitsky basis index is outside the homogeneous rank"))?;
        Ok(Tensor::monomial(slice.generators[pivot].clone()))
    }

    pub(crate) fn rank(&self, grade: &Grade) -> Result<usize> {
        if !grade.supports_states() {
            return Ok(0);
        }
        let weight = grade.degree as u128 / 2;
        let dual = 2 * grade.length as u128 - weight;
        monomial_count(self.gram.size(), self.field.degree(), weight.min(dual))
    }

    pub(crate) fn relations(&self, model: &mut Model, grade: &Grade) -> Result<String> {
        if !grade.supports_states() || grade.degree as u128 / 2 <= grade.length as u128 {
            return Ok("[]".to_owned());
        }
        let slice = self.slice(model, grade)?;
        let Some(reduction) = &slice.reduction else {
            return Ok("[]".to_owned());
        };
        let mut result = "[".to_owned();
        let mut separator = "";
        for (column, powers) in slice.generators.iter().enumerate() {
            if reduction.pivots.binary_search(&column).is_ok() {
                continue;
            }
            let mut expression = Tensor::monomial(powers.clone());
            for (&pivot, row) in reduction.pivots.iter().zip(&reduction.rows) {
                expression.insert(model, slice.generators[pivot].clone(), row[column].negate())?;
            }
            write!(result, "{separator}{}", expression.text(model, self)?).unwrap();
            if result.len() >= MAX_SHOWN_BYTES {
                return Err(Error::plain(
                    "Verbitsky relation report exceeds the output limit",
                ));
            }
            separator = ", ";
        }
        result.push(']');
        Ok(result)
    }

    pub(crate) fn lefschetz(
        &self,
        model: &mut Model,
        grade: &Grade,
        h: &Tensor,
        primitive_index: Option<usize>,
    ) -> Result<Lefschetz> {
        self.validate_grade(grade)?;
        if grade.degree as u128 / 2 > grade.length as u128 {
            return Err(Error::plain(
                "Verbitsky primitive degree must be at most the middle cohomological degree",
            ));
        }
        if self.divisor_square(model, h)?.is_zero() {
            return Err(Error::plain(
                "Verbitsky Lefschetz divisor must be nonisotropic",
            ));
        }
        let slice = self.slice(model, grade)?;
        let weight = grade.degree as u128 / 2;
        let exponent = 2 * (grade.length as u128 - weight);
        let h_power = tensor_power(model, h, exponent)?;
        let mut form = linear::matrix(slice.rank(), slice.rank())?;
        let mut moments = Moments::new(self);
        for (i, a) in slice.generators.iter().enumerate() {
            for (j, b) in slice.generators.iter().enumerate().skip(i) {
                let powers = multiply_powers(a, b)?;
                let mut sum = Rational::zero();
                for (other, &coefficient) in &h_power.terms {
                    let product = multiply_powers(&powers, other)?;
                    let moment = moments.monomial(model, &product)?;
                    let value = coefficient.multiply(model, moment)?;
                    sum = sum.add(model, value)?;
                }
                form[i][j] = sum;
                form[j][i] = sum;
            }
        }
        let lower = if weight == 0 {
            Vec::new()
        } else {
            self.slice(
                model,
                &Grade {
                    degree: grade.degree - 2,
                    ..grade.clone()
                },
            )?
            .generators
            .clone()
        };
        let mut map = linear::matrix(lower.len(), slice.rank())?;
        // hA^(i-1) is tested by the form with h^(2n-2i). Its
        // orthogonal complement is the kernel of h^(2n-2i+1).
        for (row, base) in map.iter_mut().zip(&lower) {
            for (term, &coefficient) in &h.terms {
                let lifted = multiply_powers(base, term)?;
                let index = slice
                    .generators
                    .binary_search(&lifted)
                    .map_err(|_| grading_error())?;
                linear::add_multiple(model, row, &form[index], coefficient)?;
            }
        }
        let reduction = linear::rref(model, map, slice.rank())?;
        let kernel = reduction.kernel(slice.rank())?;
        if let Some(index) = primitive_index {
            let coordinates = kernel.get(index).ok_or_else(|| {
                Error::plain("Verbitsky primitive index is outside the kernel dimension")
            })?;
            return Ok(Lefschetz::Primitive(slice.expression(model, coordinates)?));
        }
        let restricted = linear::restrict(model, &form, &kernel)?;
        Ok(Lefschetz::Signature {
            full: linear::inertia(model, form)?,
            primitive: linear::inertia(model, restricted)?,
        })
    }

    fn divisor_square(&self, model: &mut Model, h: &Tensor) -> Result<Rational> {
        let mut sum = Rational::zero();
        for (left, &c) in &h.terms {
            for (right, &d) in &h.terms {
                let (Some(a), Some(b)) = (left.first(), right.first()) else {
                    return Err(grading_error());
                };
                let q = *self.gram.entry(a.index.into(), b.index.into());
                let product = c.multiply(model, d)?.multiply(model, q)?;
                sum = sum.add(model, product)?;
            }
        }
        Ok(sum)
    }
}

pub(crate) enum Lefschetz {
    Signature {
        full: [usize; 3],
        primitive: [usize; 3],
    },
    Primitive(Tensor),
}

fn tensor_power(model: &mut Model, tensor: &Tensor, mut exponent: u128) -> Result<Tensor> {
    let mut result = Tensor::unit();
    let mut base = tensor.clone();
    while exponent > 0 {
        if !exponent.is_multiple_of(2) {
            result = result.cup(model, &base)?;
        }
        exponent /= 2;
        if exponent > 0 {
            base = base.cup(model, &base)?;
        }
    }
    Ok(result)
}

fn monomial_count(divisors: usize, quadratics: usize, weight: u128) -> Result<usize> {
    if divisors == 0 {
        if !weight.is_multiple_of(2) {
            return Ok(0);
        }
        return symmetric_count(quadratics, weight / 2);
    }
    if weight / 2 >= MAX_PIVOT_CELLS as u128 {
        return Err(work_limit());
    }
    let mut sum = 0usize;
    for pairs in 0..=weight / 2 {
        let a = symmetric_count(divisors, weight - 2 * pairs)?;
        let b = symmetric_count(quadratics, pairs)?;
        sum = sum.saturating_add(a.saturating_mul(b));
        storage(sum)?;
    }
    Ok(sum)
}

fn symmetric_count(variables: usize, weight: u128) -> Result<usize> {
    let mut count = 1u128;
    for i in 1..variables as u128 {
        count = weight
            .checked_add(i)
            .and_then(|j| count.checked_mul(j))
            .ok_or_else(work_limit)?
            / i;
        if count > MAX_PIVOT_CELLS as u128 {
            return Err(work_limit());
        }
    }
    Ok(count as usize)
}

fn monomials(model: &mut Model, form: &Form, weight: u128) -> Result<Vec<Monomial>> {
    let divisors = form.gram.size();
    let quadratics = form.field.degree();
    let count = monomial_count(divisors, quadratics, weight)?;
    storage(count.saturating_mul(8 + 5 * (divisors + quadratics)))?;
    model.spend(count).map_err(|_| work_limit())?;
    let mut result = Vec::with_capacity(count);
    if count > 0 {
        let mut weights = vec![1; divisors];
        weights.extend(std::iter::repeat_n(2, quadratics));
        enumerate(model, &weights, 0, weight, &mut Vec::new(), &mut result)?;
        result.sort_unstable();
    }
    Ok(result)
}

fn enumerate(
    model: &mut Model,
    weights: &[u128],
    index: usize,
    left: u128,
    powers: &mut Monomial,
    output: &mut Vec<Monomial>,
) -> Result<()> {
    model.spend(1)?;
    if left == 0 {
        output.push(powers.clone());
        return Ok(());
    }
    let weight = weights[index];
    if index + 1 == weights.len() {
        if left.is_multiple_of(weight) {
            powers.push(Power {
                index: index as u8,
                exponent: left / weight,
            });
            output.push(powers.clone());
            powers.pop();
        }
        return Ok(());
    }
    if weight == 2 && !left.is_multiple_of(2) {
        return Ok(());
    }
    for exponent in 0..=left / weight {
        if exponent > 0 {
            powers.push(Power {
                index: index as u8,
                exponent,
            });
        }
        enumerate(
            model,
            weights,
            index + 1,
            left - exponent * weight,
            powers,
            output,
        )?;
        if exponent > 0 {
            powers.pop();
        }
    }
    Ok(())
}

fn grading_error() -> Error {
    Error::plain("inconsistent monomial grading in Verbitsky Hodge ring")
}
