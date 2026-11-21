use std::collections::BTreeMap;

use crate::binary_field::{BinaryField, Embedding};
use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::rv_descent;
use crate::rv_hom::Grade;
use crate::rv_laurent::{self, add_polynomial_product, Polynomial, Term};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Norm,
    Coefficient(usize),
    Slopes,
}

impl Operation {
    fn name(self) -> &'static str {
        match self {
            Self::Norm => "rvUNorm",
            Self::Coefficient(_) => "rvNewtonCoeff",
            Self::Slopes => "rvNewton",
        }
    }
}

pub(crate) fn plan(
    model: &mut Model,
    ty: TypeId,
    modulus: usize,
    root: usize,
    operation: Operation,
) -> Result<(TypeId, Embedding)> {
    let grade = rv_laurent::grade(model, ty, operation.name())?;
    if grade.source != grade.target {
        return Err(Error::plain(format!(
            "{} requires a Laurent Margolis endomorphism",
            operation.name()
        )));
    }
    let (embedding, point) =
        rv_descent::rational_point(model, &grade, modulus, root, operation.name())?;
    let output = if matches!(operation, Operation::Coefficient(_)) {
        let one = model.named("One").unwrap();
        model.allocate(1 + point.len())?;
        model.rv_laurent_type(Grade {
            source: one,
            target: one,
            field: embedding.small(),
            point: point.into(),
            ..grade
        })?
    } else {
        ty
    };
    Ok((output, embedding))
}

fn norm_and_characteristic(
    model: &mut Model,
    ty: TypeId,
    embedding: &Embedding,
    support: &[Term],
    operation: Operation,
) -> Result<(Vec<Term>, Vec<Polynomial>)> {
    let field = embedding.large();
    let step = embedding.small().degree();
    let order = field.degree() / step;
    let size = rv_laurent::dimensions(model, ty)?.0;
    model.spend(support.len())?;
    let mut norm = support.to_vec();
    for _ in 1..order {
        let conjugate = rv_laurent::frobenius(model, field, &norm, step)?;
        norm = rv_laurent::product(model, field, size, size, support, &conjugate)?;
    }
    let characteristic = characteristic(model, field, size, &norm)?;
    if characteristic[0].is_empty() {
        return Err(Error::plain(format!(
            "{} requires an invertible Laurent endomorphism over L((u))",
            operation.name()
        )));
    }
    Ok((norm, characteristic))
}

pub(crate) fn evaluate(
    model: &mut Model,
    ty: TypeId,
    embedding: &Embedding,
    support: &[Term],
    operation: Operation,
) -> Result<Value> {
    let (norm, characteristic) = norm_and_characteristic(model, ty, embedding, support, operation)?;
    match operation {
        Operation::Norm => rv_laurent::persist(model, norm),
        Operation::Coefficient(degree) => {
            let coefficient = characteristic.get(degree).map_or(&[][..], Vec::as_slice);
            let mut terms = Vec::new();
            for &(power, value) in coefficient {
                let coefficient = embedding.preimage(model, value)?.ok_or_else(|| {
                    Error::plain(
                        "Newton characteristic coefficient does not lie in the specified subfield",
                    )
                })?;
                terms.push(Term {
                    coordinate: 0,
                    power,
                    coefficient,
                });
            }
            rv_laurent::persist(model, terms)
        }
        Operation::Slopes => unreachable!(),
    }
}

fn push_polynomial(
    output: &mut Vec<Polynomial>,
    cells: &mut usize,
    value: Polynomial,
) -> Result<()> {
    *cells += 1 + 2 * value.len();
    check_storage(*cells)?;
    output.push(value);
    Ok(())
}

pub(crate) fn characteristic(
    model: &mut Model,
    field: BinaryField,
    size: usize,
    norm: &[Term],
) -> Result<Vec<Polynomial>> {
    let (matrix, _) = rv_laurent::polynomial_matrix(model, size, norm)?;

    // Division-free Berkowitz recurrence on trailing principal blocks.
    // For a block [a R; S M], the descending coefficient vector is the
    // truncated convolution of char(M) with [1, a, RS, RMS, ...] in char 2.
    let one = vec![(0, 1)];
    let mut coefficients = vec![one.clone()];
    for start in (0..size).rev() {
        let degree = size - start;
        let mut series = vec![one.clone(), matrix[start][start].clone()];
        let mut cells = 4 + 2 * series[1].len();
        let tail = start + 1;
        if degree > 1 {
            let mut image: Vec<_> = matrix[tail..]
                .iter()
                .map(|row| row[start].clone())
                .collect();
            for power in 2..=degree {
                let mut coefficient = BTreeMap::new();
                model.spend(image.len())?;
                for (a, b) in matrix[start][tail..].iter().zip(&image) {
                    add_polynomial_product(model, field, &mut coefficient, a, b)?;
                }
                push_polynomial(&mut series, &mut cells, coefficient.into_iter().collect())?;
                if power < degree {
                    let mut next = Vec::with_capacity(image.len());
                    let mut image_cells = 0;
                    for row in &matrix[tail..] {
                        let mut output = BTreeMap::new();
                        model.spend(image.len())?;
                        for (a, b) in row[tail..].iter().zip(&image) {
                            add_polynomial_product(model, field, &mut output, a, b)?;
                        }
                        push_polynomial(&mut next, &mut image_cells, output.into_iter().collect())?;
                    }
                    image = next;
                }
            }
        }
        let mut next = Vec::with_capacity(degree + 1);
        let mut result_cells = 0;
        for i in 0..=degree {
            let mut coefficient = BTreeMap::new();
            let count = (i + 1).min(coefficients.len());
            model.spend(count)?;
            for (j, item) in coefficients.iter().enumerate().take(count) {
                add_polynomial_product(model, field, &mut coefficient, &series[i - j], item)?;
            }
            push_polynomial(
                &mut next,
                &mut result_cells,
                coefficient.into_iter().collect(),
            )?;
        }
        coefficients = next;
    }
    coefficients.reverse();
    Ok(coefficients)
}

#[derive(Clone, Copy)]
struct Vertex {
    degree: usize,
    valuation: i64,
}

fn at_least_as_steep(a: Vertex, b: Vertex, c: Vertex) -> bool {
    let first = i128::from(b.valuation) - i128::from(a.valuation);
    let second = i128::from(c.valuation) - i128::from(b.valuation);
    first * (c.degree - b.degree) as i128 >= second * (b.degree - a.degree) as i128
}

fn gcd(mut left: u128, mut right: u128) -> u128 {
    while right != 0 {
        (left, right) = (right, left % right);
    }
    left
}

fn slope_list(
    model: &mut Model,
    characteristic: &[Polynomial],
    extension: usize,
) -> Result<Vec<String>> {
    let mut vertices = Vec::new();
    for (degree, coefficient) in characteristic.iter().enumerate() {
        model.spend(1)?;
        let Some(&(valuation, _)) = coefficient.first() else {
            continue;
        };
        let next = Vertex { degree, valuation };
        while let [.., a, b] = &vertices[..] {
            if !at_least_as_steep(*a, *b, next) {
                break;
            }
            model.spend(1)?;
            vertices.pop();
        }
        vertices.push(next);
    }
    let mut slopes = Vec::new();
    for edge in vertices.windows(2).rev() {
        let [a, b] = edge else { unreachable!() };
        let width = b.degree - a.degree;
        let numerator = i128::from(a.valuation) - i128::from(b.valuation);
        let denominator = width as i128 * extension as i128;
        let divisor = gcd(numerator.unsigned_abs(), denominator as u128) as i128;
        let (numerator, denominator) = (numerator / divisor, denominator / divisor);
        let slope = if denominator == 1 {
            format!("{numerator}")
        } else {
            format!("{numerator}/{denominator}")
        };
        slopes.push(format!("{slope}: {width}"));
    }
    Ok(slopes)
}

pub(crate) fn slopes(
    model: &mut Model,
    ty: TypeId,
    embedding: &Embedding,
    support: &[Term],
) -> Result<String> {
    let (_, characteristic) =
        norm_and_characteristic(model, ty, embedding, support, Operation::Slopes)?;
    let order = embedding.large().degree() / embedding.small().degree();
    Ok(slope_list(model, &characteristic, order)?.join(", "))
}
