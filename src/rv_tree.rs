use std::collections::BTreeMap;

use crate::binary_field::BinaryField;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId};
use crate::rv_descent;
use crate::rv_hodge::{self, Boundary};
use crate::rv_laurent::{self, add_polynomial_product, Polynomial, Term};

const MAX_VERTICES: usize = 4096;

pub(crate) enum Report {
    Shell,
    Ball,
    Count(Boundary, Vec<i64>),
}

impl Report {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            Self::Shell => "rvTreeShell",
            Self::Ball => "rvTreeBall",
            Self::Count(Boundary::Exact, _) => "rvTreeCount",
            Self::Count(Boundary::Closure, _) => "rvTreeClosureCount",
        }
    }
}

pub(crate) struct Plan {
    field: BinaryField,
    step: usize,
    shells: Vec<Shell>,
}

struct Shell {
    radius: i64,
    first: usize,
    second: usize,
}

fn limit() -> Error {
    Error::plain(format!(
        "Bruhat–Tits tree query exceeds the limit of {MAX_VERTICES} vertices"
    ))
}

fn shells(order: usize, radius: usize, ball: bool) -> Result<Vec<Shell>> {
    let mut shells = Vec::new();
    let mut total = 0;
    if ball || radius == 0 {
        total = 1;
        shells.push(Shell {
            radius: 0,
            first: 0,
            second: 0,
        });
    }
    let mut power = 1usize;
    for r in 1..=radius {
        let previous = power;
        power = power
            .checked_mul(order)
            .filter(|&count| count <= MAX_VERTICES)
            .ok_or_else(limit)?;
        if ball || r == radius {
            total += power + previous;
            if total > MAX_VERTICES {
                return Err(limit());
            }
            shells.push(Shell {
                radius: r as i64,
                first: power,
                second: previous,
            });
        }
    }
    Ok(shells)
}

pub(crate) fn plan(
    model: &mut Model,
    ty: TypeId,
    modulus: usize,
    root: usize,
    radius: usize,
    report: &Report,
) -> Result<Plan> {
    let name = report.name();
    let grade = rv_hodge::square(model, ty, name)?;
    if rv_laurent::dimensions(model, ty)?.0 != 2 {
        return Err(Error::plain(format!(
            "{name} requires a rank-two Laurent Margolis endomorphism"
        )));
    }
    if let Report::Count(_, mu) = report {
        rv_hodge::check_coweight(model, ty, mu)?;
    }
    let (embedding, _) = rv_descent::rational_point(model, &grade, modulus, root, name)?;
    let range = shells(
        grade.field.order(),
        radius,
        !matches!(report, Report::Shell),
    )?;
    Ok(Plan {
        field: grade.field,
        step: embedding.small().degree(),
        shells: range,
    })
}

#[derive(Clone, Copy)]
struct Matrix<'a> {
    a: &'a Polynomial,
    b: &'a Polynomial,
    c: &'a Polynomial,
    d: &'a Polynomial,
}

impl<'a> Matrix<'a> {
    fn swap(self) -> Self {
        Self {
            a: self.d,
            b: self.c,
            c: self.b,
            d: self.a,
        }
    }
}

fn polynomial_map(model: &mut Model, polynomial: &Polynomial) -> Result<BTreeMap<i64, usize>> {
    model.spend(polynomial.len())?;
    Ok(polynomial.iter().copied().collect())
}

fn residue(
    model: &mut Model,
    field: BinaryField,
    step: usize,
    mut index: usize,
    mut power: i64,
) -> Result<(Polynomial, Polynomial)> {
    let mut value = Vec::new();
    let mut conjugate = Vec::new();
    while index != 0 {
        model.spend(1)?;
        let coefficient = index % field.order();
        index /= field.order();
        if coefficient != 0 {
            value.push((power, coefficient));
            conjugate.push((power, field.frobenius(model, coefficient, step)?));
        }
        power += 1;
    }
    Ok((value, conjugate))
}

fn lower_position(
    model: &mut Model,
    field: BinaryField,
    matrix: Matrix<'_>,
    determinant: i128,
    radius: i64,
    alpha: &Polynomial,
    conjugate: &Polynomial,
) -> Result<[i128; 2]> {
    // For X = [1 0; alpha u^r], adj(X) A sigma(X) has entries
    // u^r(a+b sigma(alpha)), u^(2r)b,
    // c+d sigma(alpha)+alpha(a+b sigma(alpha)), u^r(d+alpha b).
    let mut first = polynomial_map(model, matrix.a)?;
    add_polynomial_product(model, field, &mut first, matrix.b, conjugate)?;

    let mut lower = polynomial_map(model, matrix.c)?;
    add_polynomial_product(model, field, &mut lower, matrix.d, conjugate)?;
    let shift = i128::from(radius);
    let upper_left = first
        .first_key_value()
        .map(|(&power, _)| i128::from(power) + shift);
    let first = first.into_iter().collect();
    add_polynomial_product(model, field, &mut lower, alpha, &first)?;

    let mut last = polynomial_map(model, matrix.d)?;
    add_polynomial_product(model, field, &mut last, alpha, matrix.b)?;
    let least = [
        upper_left,
        matrix
            .b
            .first()
            .map(|&(power, _)| i128::from(power) + 2 * shift),
        lower.first_key_value().map(|(&power, _)| i128::from(power)),
        last.first_key_value()
            .map(|(&power, _)| i128::from(power) + shift),
    ]
    .into_iter()
    .flatten()
    .min()
    .unwrap();
    let first = least - shift;
    Ok([first, determinant - first])
}

fn add_vertex(
    model: &mut Model,
    histogram: &mut BTreeMap<[i128; 2], usize>,
    position: [i128; 2],
) -> Result<()> {
    model.spend(1 + histogram.len().max(1).ilog2() as usize)?;
    *histogram.entry(position).or_insert(0) += 1;
    Ok(())
}

pub(crate) fn evaluate(
    model: &mut Model,
    plan: Plan,
    support: &[Term],
    report: &Report,
    modulus: usize,
    root: usize,
    radius: usize,
) -> Result<String> {
    let Plan {
        field,
        step,
        shells,
    } = plan;
    let (table, _) = rv_laurent::polynomial_matrix(model, 2, support)?;
    let matrix = Matrix {
        a: &table[0][0],
        b: &table[0][1],
        c: &table[1][0],
        d: &table[1][1],
    };
    let mut determinant = BTreeMap::new();
    add_polynomial_product(model, field, &mut determinant, matrix.a, matrix.d)?;
    add_polynomial_product(model, field, &mut determinant, matrix.b, matrix.c)?;
    let (&valuation, _) = determinant.first_key_value().ok_or_else(|| {
        Error::plain(format!(
            "{} requires an invertible rank-two Laurent endomorphism",
            report.name()
        ))
    })?;
    let determinant = i128::from(valuation);
    let mut histogram = BTreeMap::new();
    for shell in shells {
        if shell.radius == 0 {
            let first = i128::from(support.iter().map(|term| term.power).min().unwrap());
            add_vertex(model, &mut histogram, [first, determinant - first])?;
            continue;
        }
        for index in 0..shell.first {
            let (alpha, conjugate) = residue(model, field, step, index, 0)?;
            let position = lower_position(
                model,
                field,
                matrix,
                determinant,
                shell.radius,
                &alpha,
                &conjugate,
            )?;
            add_vertex(model, &mut histogram, position)?;
        }
        // The second chart interchanges the two coordinates and restricts its
        // free residue to u times a polynomial of degree less than r-1.
        let swapped = matrix.swap();
        for index in 0..shell.second {
            let (alpha, conjugate) = residue(model, field, step, index, 1)?;
            let position = lower_position(
                model,
                field,
                swapped,
                determinant,
                shell.radius,
                &alpha,
                &conjugate,
            )?;
            add_vertex(model, &mut histogram, position)?;
        }
    }
    let name = report.name();
    if let Report::Count(boundary, mu) = report {
        let mut count = 0;
        for (position, multiplicity) in &histogram {
            if boundary.contains(model, position, mu)? {
                count += multiplicity;
            }
        }
        Ok(format!(
            "{name}[{modulus}, {root}, {radius}, {mu:?}] = {count}"
        ))
    } else {
        Ok(format!(
            "{name}[{modulus}, {root}, {radius}] = {histogram:?}"
        ))
    }
}
