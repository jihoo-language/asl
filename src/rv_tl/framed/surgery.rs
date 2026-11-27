use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_cyclotomic;
use crate::rv_lens::{self, Rational};

use super::{crossing, trace, twist_exponent, Context, Grade, Link, MAX_WORD};
use crate::rv_tl::power;

#[derive(Clone, Copy)]
pub(crate) enum Surgery {
    Wrt,
    Doubled,
}

fn linking_arithmetic() -> Error {
    Error::plain("framed linking matrix integer arithmetic is out of range")
}

fn divide_positive(n: i128, d: u128) -> i128 {
    match i128::try_from(d) {
        Ok(d) => n / d,
        Err(_) => -i128::from(n != 0),
    }
}

fn pivot_swap(matrix: &mut [Vec<i128>], from: usize, to: usize) {
    if from != to {
        matrix.swap(from, to);
        for row in matrix.iter_mut() {
            row.swap(from, to);
        }
    }
}

fn schur(
    model: &mut Model,
    matrix: &mut [Vec<i128>],
    start: usize,
    width: usize,
    pivot: i128,
) -> Result<()> {
    let next = start + width;
    let n = matrix.len();
    let remaining = n - next;
    model.spend(width * remaining)?;
    let across = matrix[start..next].iter().flat_map(|row| &row[next..]);
    if across.clone().all(|&value| value == 0) {
        return Ok(());
    }
    // Divide a common positive content out of the pivot and its row
    // before computing a positive integral multiple of the Schur form.
    let g = across.fold(pivot.unsigned_abs(), |g, x| {
        rv_lens::gcd(g, x.unsigned_abs())
    });
    let scale = divide_positive(pivot, g);
    let mut content = 0;
    model.spend(remaining * (remaining + 1))?;
    for i in next..n {
        let (through, following) = matrix.split_at_mut(i + 1);
        let (row, pivots) = through.split_last_mut().unwrap();
        let left = divide_positive(pivots[start][i], g);
        let right = if width == 2 {
            divide_positive(pivots[start + 1][i], g)
        } else {
            0
        };
        for (j, entry) in row.iter_mut().enumerate().skip(i) {
            let (main, other) = if width == 1 {
                (pivots[start][j], 0)
            } else {
                (pivots[start + 1][j], pivots[start][j])
            };
            let value = scale
                .checked_mul(*entry)
                .and_then(|v| left.checked_mul(main).and_then(|term| v.checked_sub(term)))
                .and_then(|v| {
                    right
                        .checked_mul(other)
                        .and_then(|term| v.checked_sub(term))
                })
                .and_then(|v| if pivot > 0 { Some(v) } else { v.checked_neg() })
                .ok_or_else(linking_arithmetic)?;
            content = rv_lens::gcd(content, value.unsigned_abs());
            *entry = value;
            if j > i {
                following[j - i - 1][i] = value;
            }
        }
    }
    if content > 1 {
        model.spend(remaining * remaining)?;
        for row in &mut matrix[next..] {
            for value in &mut row[next..] {
                *value = divide_positive(*value, content);
            }
        }
    }
    Ok(())
}

fn inertia(model: &mut Model, mut matrix: Vec<Vec<i128>>) -> Result<(usize, usize)> {
    let n = matrix.len();
    let (mut positive, mut negative) = (0, 0);
    let mut start = 0;
    while start < n {
        model.spend(n - start)?;
        let diagonal = (start..n)
            .filter(|&i| matrix[i][i] != 0)
            .min_by_key(|&i| matrix[i][i].unsigned_abs());
        if let Some(diagonal) = diagonal {
            pivot_swap(&mut matrix, diagonal, start);
            let pivot = matrix[start][start];
            positive += usize::from(pivot > 0);
            negative += usize::from(pivot < 0);
            schur(model, &mut matrix, start, 1, pivot)?;
            start += 1;
        } else {
            let mut off_diagonal = None;
            for (i, row) in matrix.iter().enumerate().skip(start) {
                model.spend(n - i - 1)?;
                if let Some(j) = (i + 1..n).find(|&j| row[j] != 0) {
                    off_diagonal = Some((i, j));
                    break;
                }
            }
            let Some((i, j)) = off_diagonal else {
                break;
            };
            pivot_swap(&mut matrix, i, start);
            pivot_swap(&mut matrix, j, start + 1);
            let pivot = matrix[start][start + 1];
            // The form [[0,pivot],[pivot,0]] has one sign of each
            // kind. Its inverse swaps the two entries and divides.
            positive += 1;
            negative += 1;
            schur(model, &mut matrix, start, 2, pivot)?;
            start += 2;
        }
    }
    Ok((positive, negative))
}

fn partition_components(model: &mut Model, link: &Link) -> Result<(Vec<usize>, Vec<usize>)> {
    model.spend(link.component_of.len() + link.braid.len())?;
    let mut moving = link.component_of.clone();
    let mut touched = vec![false; link.components()];
    for &letter in link.braid.iter().rev() {
        let (index, _) = crossing(letter);
        touched[moving[index - 1]] = true;
        touched[moving[index]] = true;
        moving.swap(index - 1, index);
    }
    Ok((0..link.components()).partition(|&component| touched[component]))
}

fn coloring_count(model: &mut Model, grade: Grade, link: &Link, active: &[usize]) -> Result<usize> {
    let largest_cable = grade.level * grade.level * link.braid.len();
    if largest_cable > MAX_WORD {
        return Err(Error::plain(
            "Kirby colors include a cabled braid exceeding the limit of 4096 elementary crossings",
        ));
    }
    let per_coloring = 1 + active.len() + link.component_of.len() + link.braid.len();
    let count = active.iter().try_fold(1usize, |total, _| {
        total
            .checked_mul(grade.level + 1)
            .ok_or_else(|| Error::plain("Kirby colorings exceed the computation limit"))
    })?;
    let work = count
        .checked_mul(per_coloring)
        .ok_or_else(|| Error::plain("Kirby colorings exceed the computation limit"))?;
    model.spend(work)?;
    Ok(count)
}

fn raw_sum(
    model: &mut Model,
    grade: Grade,
    link: &Link,
    active: &[usize],
    framings: &[i128],
    context: &mut Context,
) -> Result<Rational> {
    if active.is_empty() {
        return Ok(Rational::integer(1));
    }
    let count = coloring_count(model, grade, link, active)?;
    let mut result = Rational::zero();
    let mut colors = vec![0; link.components()];
    for _ in 0..count {
        let value = trace(model, grade, link, &colors, framings, context)?;
        if !value.is_zero() {
            let mut weighted = value;
            for &component in active {
                let color = colors[component];
                if color > 0 {
                    let dim = context.dimension(model, grade, color, false)?;
                    weighted = rv_lens::multiply(model, grade.order, &weighted, &dim)?;
                }
            }
            result = rv_lens::sum(model, &result, &weighted)?;
        }
        for &component in active {
            if colors[component] < grade.level {
                colors[component] += 1;
                break;
            }
            colors[component] = 0;
        }
    }
    Ok(result)
}

fn total_dimension_inverse(model: &mut Model, grade: Grade) -> Result<Rational> {
    // D^2=N/(2 sin(pi/N)^2), with N=k+2 and A=e^(pi i/(2N)).
    let mut coefficients = vec![0; grade.order];
    coefficients[0] = 2;
    coefficients[4] -= 1;
    coefficients[grade.order - 4] -= 1;
    let numerator = rv_cyclotomic::reduced(model, grade.order, coefficients)?;
    rv_lens::fraction(model, &numerator, (2 * (grade.level + 2)) as u128)
}

fn gauss(model: &mut Model, grade: Grade, squares: &[Rational], framing: i128) -> Result<Rational> {
    let mut result = Rational::zero();
    for (color, dim_square) in squares.iter().enumerate() {
        let weighted = rv_lens::shifted(
            model,
            grade.order,
            dim_square,
            twist_exponent(grade, color, framing),
        )?;
        result = rv_lens::sum(model, &result, &weighted)?;
    }
    Ok(result)
}

fn hermitian_square(model: &mut Model, grade: Grade, value: &Rational) -> Result<Rational> {
    let polynomial =
        rv_cyclotomic::substitute_root(model, grade.order, &value.numerator, grade.order - 1)?;
    let conjugate = rv_lens::fraction(model, &polynomial, value.denominator)?;
    rv_lens::multiply(model, grade.order, value, &conjugate)
}

pub(crate) fn evaluate(
    model: &mut Model,
    grade: Grade,
    word: &[i128],
    framings: &[i128],
    kind: Surgery,
) -> Result<Value> {
    let link = Link::new(model, grade.strands, word)?;
    link.check_framings(framings)?;
    if link.components() == 0 {
        let value = match kind {
            Surgery::Wrt => Rational::integer(1),
            Surgery::Doubled => total_dimension_inverse(model, grade)?,
        };
        return rv_lens::rational_value(model, value);
    }
    let (active, split) = partition_components(model, &link)?;
    let signature = match kind {
        Surgery::Wrt => {
            let matrix = link.matrix(model, framings, &active)?;
            let (positive, negative) = inertia(model, matrix)?;
            if positive + negative != active.len() || split.iter().any(|&i| framings[i] == 0) {
                return Err(Error::plain(
                    "rvTLSurgery requires a nonsingular linking matrix",
                ));
            }
            Some((positive, negative))
        }
        Surgery::Doubled => None,
    };
    let mut context = Context::default();
    let mut sum = raw_sum(model, grade, &link, &active, framings, &mut context)?;
    if sum.is_zero() {
        return rv_lens::rational_value(model, sum);
    }
    let inverse_d2 = total_dimension_inverse(model, grade)?;
    let mut split_framings = BTreeMap::<_, usize>::new();
    for &i in &split {
        model.spend(1 + split_framings.len().max(1).ilog2() as usize)?;
        let framing = framings[i];
        let positive = signature.is_some() && framing > 0;
        let phase = framing.rem_euclid(grade.order as i128);
        *split_framings.entry((phase, positive)).or_default() += 1;
    }
    let mut squares = Vec::with_capacity(grade.level + 1);
    for color in 0..=grade.level {
        let dim = context.dimension(model, grade, color, false)?;
        squares.push(rv_lens::multiply(model, grade.order, &dim, &dim)?);
    }
    if let Some((positive, negative)) = signature {
        // Gauss normalizers obey p_+ p_-=D^2 in this pivotal convention.
        let gauss_minus = gauss(model, grade, &squares, -1)?;
        let gauss_plus = gauss(model, grade, &squares, 1)?;
        let recip_plus = rv_lens::multiply(model, grade.order, &gauss_minus, &inverse_d2)?;
        let recip_minus = rv_lens::multiply(model, grade.order, &gauss_plus, &inverse_d2)?;
        for inverse in std::iter::repeat_n(&recip_plus, positive)
            .chain(std::iter::repeat_n(&recip_minus, negative))
        {
            sum = rv_lens::multiply(model, grade.order, &sum, inverse)?;
        }
        for &component in &split {
            let framing = framings[component];
            let inverse = if framing > 0 {
                &recip_plus
            } else {
                &recip_minus
            };
            let value = gauss(model, grade, &squares, framing)?;
            let factor = rv_lens::multiply(model, grade.order, &value, inverse)?;
            sum = rv_lens::multiply(model, grade.order, &sum, &factor)?;
        }
    } else {
        sum = hermitian_square(model, grade, &sum)?;
        for _ in 0..=active.len() {
            sum = rv_lens::multiply(model, grade.order, &sum, &inverse_d2)?;
        }
        for (&(phase, _), &multiplicity) in &split_framings {
            let value = gauss(model, grade, &squares, phase)?;
            let norm = hermitian_square(model, grade, &value)?;
            let factor = rv_lens::multiply(model, grade.order, &norm, &inverse_d2)?;
            let factor = power(model, grade.order, factor, multiplicity)?;
            sum = rv_lens::multiply(model, grade.order, &sum, &factor)?;
        }
    }
    rv_lens::rational_value(model, sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integral_inertia_handles_hyperbolic_pivots_and_lattice_slides() {
        let model = &mut Model::new();
        let mut check = |matrix, expected| {
            model.reset_work();
            let actual = inertia(model, matrix).unwrap_or_else(|error| panic!("{error}"));
            assert_eq!(actual, expected);
        };
        check(vec![vec![0, 2, 3], vec![2, 0, 4], vec![3, 4, 0]], (1, 2));
        check(vec![vec![0, -2, 3], vec![-2, 0, 4], vec![3, 4, 0]], (2, 1));
        check(vec![vec![0, 1, 0], vec![1, 0, 1], vec![0, 1, 0]], (1, 1));

        let base = vec![
            vec![0, 2, 0, 0, 0],
            vec![2, 0, 0, 0, 0],
            vec![0, 0, -3, 0, 0],
            vec![0, 0, 0, 1, 0],
            vec![0, 0, 0, 0, 0],
        ];
        for i in 0..base.len() {
            for j in 0..base.len() {
                if i == j {
                    continue;
                }
                let mut transform = base.clone();
                let source = transform[j].clone();
                for (target, value) in transform[i].iter_mut().zip(source) {
                    *target += 5 * value;
                }
                for row in &mut transform {
                    row[i] += 5 * row[j];
                }
                check(transform, (2, 2));
            }
        }
    }
}
