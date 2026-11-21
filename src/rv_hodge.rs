use std::collections::BTreeMap;

use crate::binary_field::BinaryField;
use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId};
use crate::rv_descent;
use crate::rv_hom::Grade;
use crate::rv_laurent::{self, add_polynomial_product, add_powers, Polynomial, Term};
use crate::rv_newton;

#[derive(Clone, Copy)]
pub(crate) enum Boundary {
    Exact,
    Closure,
}

impl Boundary {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Exact => "rvInADLV",
            Self::Closure => "rvInADLVClosure",
        }
    }
}

pub(crate) fn square(model: &Model, ty: TypeId, name: &str) -> Result<Grade> {
    let grade = rv_laurent::grade(model, ty, name)?;
    if grade.source != grade.target {
        return Err(Error::plain(format!(
            "{name} requires a Laurent Margolis endomorphism"
        )));
    }
    Ok(grade)
}

pub(crate) fn compatible(model: &Model, a: TypeId, x: TypeId, name: &str) -> Result<Grade> {
    let grade = square(model, a, name)?;
    square(model, x, name)?;
    if a != x {
        return Err(Error::plain(format!(
            "{name} requires both Laurent frames on the same endpoint and affine field point"
        )));
    }
    Ok(grade)
}

pub(crate) fn relative_step(
    model: &mut Model,
    a: TypeId,
    x: TypeId,
    modulus: usize,
    root: usize,
    name: &str,
) -> Result<usize> {
    let grade = compatible(model, a, x, name)?;
    let (embedding, _) = rv_descent::rational_point(model, &grade, modulus, root, name)?;
    Ok(embedding.small().degree())
}

pub(crate) fn check_coweight(model: &mut Model, ty: TypeId, mu: &[i64]) -> Result<()> {
    model.spend(mu.len())?;
    if mu.windows(2).any(|pair| pair[0] > pair[1]) {
        return Err(Error::plain("the ADLV coweight must be weakly increasing"));
    }
    if mu.len() != rv_laurent::dimensions(model, ty)?.0 {
        return Err(Error::plain(
            "the ADLV coweight must have one entry per Margolis dimension",
        ));
    }
    Ok(())
}

impl Boundary {
    pub(crate) fn contains(self, model: &mut Model, position: &[i128], mu: &[i64]) -> Result<bool> {
        model.spend(mu.len())?;
        match self {
            Self::Exact => Ok(position
                .iter()
                .copied()
                .eq(mu.iter().map(|&x| i128::from(x)))),
            Self::Closure => {
                let mut actual = 0;
                let mut bound = 0;
                for (&h, &m) in position.iter().zip(mu) {
                    actual += h;
                    bound += i128::from(m);
                    if actual < bound {
                        return Ok(false);
                    }
                }
                Ok(actual == bound)
            }
        }
    }
}

pub(crate) fn hodge(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Vec<i128>> {
    let grade = square(model, ty, "rvHodge")?;
    let size = rv_laurent::dimensions(model, ty)?.0;
    determinantal_exponents(model, grade.field, size, terms, "rvHodge")
}

pub(crate) fn relative(model: &mut Model, ty: TypeId, x: &[Term], y: &[Term]) -> Result<Vec<i128>> {
    let grade = square(model, ty, "rvPosition")?;
    let size = rv_laurent::dimensions(model, ty)?.0;
    adjugate_position(model, grade.field, size, x, y, "rvPosition")
}

pub(crate) fn mu_at(
    model: &mut Model,
    ty: TypeId,
    a: &[Term],
    x: &[Term],
    step: usize,
    name: &str,
) -> Result<Vec<i128>> {
    let grade = square(model, ty, name)?;
    let size = rv_laurent::dimensions(model, ty)?.0;
    let conjugate = rv_laurent::frobenius(model, grade.field, x, step)?;
    let image = rv_laurent::product(model, grade.field, size, size, a, &conjugate)?;
    adjugate_position(model, grade.field, size, x, &image, name)
}

fn scale_polynomial(
    model: &mut Model,
    field: BinaryField,
    terms: &[Term],
    polynomial: &Polynomial,
) -> Result<Vec<Term>> {
    let mut output = BTreeMap::new();
    for &term in terms {
        model.spend(1)?;
        for &(power, coefficient) in polynomial {
            let power = add_powers(term.power, power)?;
            model.spend(field.degree())?;
            rv_laurent::accumulate(
                model,
                &mut output,
                (term.coordinate, power),
                field.multiply(term.coefficient, coefficient),
            )?;
        }
    }
    Ok(output
        .into_iter()
        .map(|((coordinate, power), coefficient)| Term {
            coordinate,
            power,
            coefficient,
        })
        .collect())
}

fn adjugate_position(
    model: &mut Model,
    field: BinaryField,
    size: usize,
    x: &[Term],
    y: &[Term],
    name: &str,
) -> Result<Vec<i128>> {
    let (transformed, valuation) = adjugate_action(model, field, size, x, y, name)?;
    let mut result = determinantal_exponents(model, field, size, &transformed, name)?;
    for exponent in &mut result {
        *exponent -= i128::from(valuation);
    }
    Ok(result)
}

pub(crate) fn adjugate_action(
    model: &mut Model,
    field: BinaryField,
    size: usize,
    x: &[Term],
    y: &[Term],
    name: &str,
) -> Result<(Vec<Term>, i64)> {
    if size == 0 {
        return Ok((Vec::new(), 0));
    }
    let characteristic = rv_newton::characteristic(model, field, size, x)?;
    let &(valuation, _) = characteristic[0]
        .first()
        .ok_or_else(|| not_invertible(name))?;
    // Cayley-Hamilton in characteristic two gives
    // adj(X)Y = (X^(n-1) + c_(n-1) X^(n-2) + ... + c_1)Y.
    model.spend(y.len())?;
    let mut transformed = y.to_vec();
    for coefficient in characteristic[1..size].iter().rev() {
        let product = rv_laurent::product(model, field, size, size, x, &transformed)?;
        let lower = scale_polynomial(model, field, y, coefficient)?;
        transformed = rv_laurent::add_support(model, &product, &lower)?;
    }
    Ok((transformed, valuation))
}

fn not_invertible(name: &str) -> Error {
    Error::plain(format!(
        "{name} requires invertible Laurent endomorphisms over L((u))"
    ))
}

fn subtract_powers(left: i64, right: i64) -> Result<i64> {
    left.checked_sub(right)
        .ok_or_else(|| Error::plain("Laurent exponent is out of the signed 64-bit range"))
}

fn exact_quotient(
    model: &mut Model,
    field: BinaryField,
    mut numerator: BTreeMap<i64, usize>,
    divisor: &Polynomial,
) -> Result<Polynomial> {
    if numerator.is_empty() || divisor == &[(0, 1)] {
        return Ok(numerator.into_iter().collect());
    }
    let minimum = subtract_powers(*numerator.first_key_value().unwrap().0, divisor[0].0)?;
    let &(leading_power, leading_coefficient) = divisor.last().unwrap();
    let inverse = field.inverse(model, leading_coefficient)?;
    let mut quotient = Vec::new();
    while let Some((&power, &coefficient)) = numerator.last_key_value() {
        let offset = subtract_powers(power, leading_power)?;
        if offset < minimum {
            return Err(Error::plain(
                "Hodge reduction encountered a nonexact Laurent quotient",
            ));
        }
        model.spend(field.degree() * (1 + divisor.len()))?;
        let scalar = field.multiply(coefficient, inverse);
        for &(power, value) in divisor {
            let power = add_powers(power, offset)?;
            rv_laurent::accumulate(model, &mut numerator, power, field.multiply(scalar, value))?;
        }
        check_storage(2 * (quotient.len() + numerator.len() + 1))?;
        quotient.push((offset, scalar));
    }
    quotient.reverse();
    Ok(quotient)
}

fn determinantal_exponents(
    model: &mut Model,
    field: BinaryField,
    size: usize,
    entries: &[Term],
    name: &str,
) -> Result<Vec<i128>> {
    Ok(
        pivots(model, field, size, entries, name, PivotOrder::Valuation)?
            .into_iter()
            .map(|pivot| pivot.exponent)
            .collect(),
    )
}

#[derive(Clone, Copy)]
pub(crate) enum PivotOrder {
    Valuation,
    Iwahori,
}

pub(crate) struct Pivot {
    pub(crate) row: usize,
    pub(crate) column: usize,
    pub(crate) exponent: i128,
}

pub(crate) fn pivots(
    model: &mut Model,
    field: BinaryField,
    size: usize,
    entries: &[Term],
    name: &str,
    order: PivotOrder,
) -> Result<Vec<Pivot>> {
    let (mut matrix, mut cells) = rv_laurent::polynomial_matrix(model, size, entries)?;
    let mut divisor = vec![(0, 1)];
    let mut previous = 0;
    let mut result = Vec::with_capacity(size);
    let mut rows: Vec<_> = (0..size).collect();
    let mut columns = rows.clone();
    for pivot in 0..size {
        model.spend((size - pivot) * (size - pivot))?;
        let mut chosen = None;
        for (row, entries) in matrix.iter().enumerate().skip(pivot) {
            for (column, polynomial) in entries.iter().enumerate().skip(pivot) {
                if let Some(&(valuation, _)) = polynomial.first() {
                    let weight = match order {
                        PivotOrder::Valuation => i128::from(valuation),
                        PivotOrder::Iwahori => {
                            (size as i128) * i128::from(valuation) + rows[row] as i128
                                - columns[column] as i128
                        }
                    };
                    let candidate = (weight, rows[row], columns[column], row, column, valuation);
                    if chosen.is_none_or(|previous| candidate < previous) {
                        chosen = Some(candidate);
                    }
                }
            }
        }
        let (_, original_row, original_column, row, column, valuation) =
            chosen.ok_or_else(|| not_invertible(name))?;
        matrix.swap(pivot, row);
        rows.swap(pivot, row);
        if column != pivot {
            model.spend(size - pivot)?;
            columns.swap(pivot, column);
            for row in &mut matrix[pivot..] {
                row.swap(pivot, column);
            }
        }

        // The current Bareiss block is the Schur complement times the
        // previous pivot, which contributes a common shift to its valuations.
        result.push(Pivot {
            row: original_row,
            column: original_column,
            exponent: i128::from(valuation) - i128::from(previous),
        });
        previous = valuation;
        for row in pivot + 1..size {
            for column in pivot + 1..size {
                let mut numerator = BTreeMap::new();
                model.spend(1)?;
                add_polynomial_product(
                    model,
                    field,
                    &mut numerator,
                    &matrix[pivot][pivot],
                    &matrix[row][column],
                )?;
                add_polynomial_product(
                    model,
                    field,
                    &mut numerator,
                    &matrix[row][pivot],
                    &matrix[pivot][column],
                )?;
                let polynomial = exact_quotient(model, field, numerator, &divisor)?;
                cells += 2 * polynomial.len();
                cells -= 2 * matrix[row][column].len();
                check_storage(cells + 2 * divisor.len())?;
                matrix[row][column] = polynomial;
            }
        }

        divisor = std::mem::take(&mut matrix[pivot][pivot]);
        cells -= 2 * divisor.len();
        for entry in &mut matrix[pivot][pivot + 1..] {
            cells -= 2 * std::mem::take(entry).len();
        }
        for row in &mut matrix[pivot + 1..] {
            cells -= 2 * std::mem::take(&mut row[pivot]).len();
        }
    }
    Ok(result)
}
