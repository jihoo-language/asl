use crate::error::Result;
use crate::integer_smith::{self, Smith, Track};
use crate::lee_linear::Rational;
use crate::model::Model;
use crate::rv_lens::gcd;
use crate::rv_qp_quad::Wide;

use super::{arithmetic, certificate};

pub(super) type Matrix = Vec<Vec<i128>>;

pub(super) fn identity(n: usize) -> Matrix {
    (0..n)
        .map(|i| (0..n).map(|j| i128::from(i == j)).collect())
        .collect()
}

pub(super) fn transpose(matrix: &Matrix) -> Matrix {
    let n = matrix.len();
    (0..n)
        .map(|j| matrix.iter().map(|row| row[j]).collect())
        .collect()
}

pub(super) fn flat(matrix: &Matrix) -> Vec<i128> {
    matrix.iter().flatten().copied().collect()
}

pub(super) fn powers(model: &mut Model, matrix: &Matrix, degree: usize) -> Result<Vec<Matrix>> {
    let mut powers = vec![identity(matrix.len())];
    for _ in 0..degree {
        powers.push(multiply(model, powers.last().unwrap(), matrix)?);
    }
    Ok(powers)
}

pub(super) fn multiply(model: &mut Model, left: &Matrix, right: &Matrix) -> Result<Matrix> {
    let n = left.len();
    model.spend(n * n * n)?;
    let mut product = vec![vec![0; n]; n];
    for (row, u) in product.iter_mut().zip(left) {
        for (&a, v) in u.iter().zip(right) {
            for (entry, &b) in row.iter_mut().zip(v) {
                *entry = add(*entry, mul(a, b)?)?;
            }
        }
    }
    Ok(product)
}

pub(super) fn apply(model: &mut Model, matrix: &Matrix, vector: &[i128]) -> Result<Vec<i128>> {
    model.spend(matrix.len() * vector.len())?;
    matrix.iter().map(|row| dot(row, vector)).collect()
}

pub(super) fn dot(row: &[i128], vector: &[i128]) -> Result<i128> {
    row.iter()
        .zip(vector)
        .try_fold(0, |sum, (&a, &b)| add(sum, mul(a, b)?))
}

pub(super) fn combine(
    model: &mut Model,
    coefficients: &[i128],
    powers: &[Matrix],
) -> Result<Matrix> {
    let n = powers[0].len();
    model.spend(coefficients.len() * n * n)?;
    let mut matrix = vec![vec![0; n]; n];
    for (&c, power) in coefficients.iter().zip(powers) {
        for (entry, &a) in matrix.iter_mut().flatten().zip(power.iter().flatten()) {
            *entry = add(*entry, mul(c, a)?)?;
        }
    }
    Ok(matrix)
}

pub(super) fn smith(model: &mut Model, matrix: &Matrix, track: Track) -> Result<Smith> {
    let mut sparse = integer_smith::Matrix::new(model, matrix.len(), matrix.len())?;
    for (i, row) in matrix.iter().enumerate() {
        for (j, &value) in row.iter().enumerate().filter(|(_, value)| **value != 0) {
            sparse.write(model, i, j, value)?;
        }
    }
    integer_smith::reduce(model, sparse, Track::None, track)
}

// Only the rational row span of a kernel matrix matters. Work with an
// equivalent primitive echelon matrix before asking Smith for an integral
// frame. This controls growth without changing the saturated kernel.
pub(super) fn kernel_frame(model: &mut Model, matrix: &Matrix) -> Result<Smith> {
    let n = matrix.len();
    let mut rows = matrix.clone();
    for row in &mut rows {
        let content = row.iter().fold(0, |g, c| gcd(g, c.unsigned_abs()));
        if content > 1 {
            let content = i128::try_from(content).map_err(|_| arithmetic())?;
            for c in row {
                *c /= content;
            }
        }
    }
    let mut rank = 0;
    for j in 0..n {
        let Some(index) = (rank..n)
            .filter(|&i| rows[i][j] != 0)
            .min_by_key(|&i| rows[i][j].unsigned_abs())
        else {
            continue;
        };
        rows.swap(rank, index);
        let pivot = rows[rank].clone();
        for (i, row) in rows
            .iter_mut()
            .enumerate()
            .filter(|(i, row)| *i != rank && row[j] != 0)
        {
            let divisor = gcd(pivot[j].unsigned_abs(), row[j].unsigned_abs());
            let divisor = i128::try_from(divisor).map_err(|_| arithmetic())?;
            let a = pivot[j] / divisor;
            let b = row[j] / divisor;
            let start = if i < rank { 0 } else { j };
            let mut entries = Vec::with_capacity(n - start);
            let mut content = Wide::<9>::from_i128(0).magnitude;
            for (&r, &s) in row[start..].iter().zip(&pivot[start..]) {
                model.spend(18)?;
                let c: Wide<9> = Wide::product(a, r)
                    .add(Wide::product(b, s).negate())
                    .map_err(|_| arithmetic())?;
                content = content.gcd(model, c.magnitude)?;
                entries.push(c);
            }
            if content.is_zero() {
                row[start..].fill(0);
                continue;
            }
            let content = Wide {
                magnitude: content,
                negative: false,
            };
            for (r, c) in row[start..].iter_mut().zip(entries) {
                *r = c
                    .divide_exact(model, content)?
                    .and_then(Wide::signed)
                    .ok_or_else(arithmetic)?;
            }
        }
        rank += 1;
    }
    let frame = smith(model, &rows, Track::Both)?;
    if frame.factors.len() != rank {
        return Err(certificate("primitive echelon and Smith ranks agree"));
    }
    Ok(frame)
}

pub(super) fn kernel(frame: &Smith, n: usize) -> Vec<Vec<i128>> {
    (frame.factors.len()..n)
        .map(|j| dense(frame.columns.basis(j), n))
        .collect()
}

pub(super) fn restrict(
    model: &mut Model,
    operator: &Matrix,
    frame: &Smith,
    kernel: &[Vec<i128>],
) -> Result<Matrix> {
    let rank = frame.factors.len();
    let n = operator.len();
    let m = kernel.len();
    let duals: Vec<_> = (0..n)
        .map(|i| dense(frame.columns.coordinate(i), n))
        .collect();
    let mut result = vec![vec![0; m]; m];
    for (j, v) in kernel.iter().enumerate() {
        let image = apply(model, operator, v)?;
        for (i, dual) in duals.iter().enumerate() {
            let coordinate = dot(dual, &image)?;
            if i < rank {
                if coordinate != 0 {
                    return Err(certificate("integral block is invariant"));
                }
            } else {
                result[i - rank][j] = coordinate;
            }
        }
        let mut recovered = vec![0; n];
        for (v, row) in kernel.iter().zip(&result) {
            for (target, &value) in recovered.iter_mut().zip(v) {
                *target = add(*target, mul(value, row[j])?)?;
            }
        }
        if recovered != image {
            return Err(certificate("coordinates in saturated block basis"));
        }
    }
    Ok(result)
}

fn dense(sparse: &integer_smith::Vector, n: usize) -> Vec<i128> {
    let mut result = vec![0; n];
    for (&i, &value) in sparse {
        result[i] = value;
    }
    result
}

// Fraction-free elimination keeps the pivots and augmented entries as
// signed minors. Back substitution uses the last nonzero pivot as a
// common denominator, so every division stays exact until final reduction.
type Minor = Wide<80>;

pub(super) fn solve(
    model: &mut Model,
    columns: &[Vec<i128>],
    target: &[i128],
) -> Result<Option<Vec<Rational>>> {
    let m = columns.len();
    let n = target.len();
    let mut augmented: Vec<Vec<Minor>> = target
        .iter()
        .enumerate()
        .map(|(i, &entry)| {
            columns
                .iter()
                .map(|column| Minor::from_i128(column[i]))
                .chain([Minor::from_i128(entry)])
                .collect()
        })
        .collect();
    let mut previous = Minor::from_i128(1);
    for j in 0..m {
        let Some(pivot) = (j..n)
            .filter(|&i| !augmented[i][j].is_zero())
            .min_by_key(|&i| augmented[i][j].magnitude.bits())
        else {
            return Ok(None);
        };
        augmented.swap(j, pivot);
        let (head, tail) = augmented.split_at_mut(j + 1);
        let source = &head[j];
        let pivot = source[j];
        for row in tail {
            let multiplier = row[j];
            for (a, &b) in row[j + 1..].iter_mut().zip(&source[j + 1..]) {
                let numerator = product(model, *a, pivot)?
                    .add(product(model, multiplier, b)?.negate())
                    .map_err(|_| arithmetic())?;
                *a = exact(model, numerator, previous)?;
            }
            row[j] = Minor::from_i128(0);
        }
        previous = pivot;
    }
    if augmented[m..].iter().any(|row| !row[m].is_zero()) {
        return Ok(None);
    }
    if target.iter().all(|&c| c == 0) {
        return Ok(Some(vec![Rational::ZERO; m]));
    }
    let mut scaled = vec![Minor::from_i128(0); m];
    for j in (0..m).rev() {
        let row = &augmented[j];
        let mut numerator = product(model, previous, row[m])?;
        for (&value, &solution) in row[j + 1..m].iter().zip(&scaled[j + 1..]) {
            numerator = numerator
                .add(product(model, value, solution)?.negate())
                .map_err(|_| arithmetic())?;
        }
        scaled[j] = exact(model, numerator, row[j])?;
    }
    let coefficients = scaled
        .into_iter()
        .map(|numerator| lower_fraction(model, numerator, previous))
        .collect::<Result<_>>()?;
    Ok(Some(coefficients))
}

fn product(model: &mut Model, a: Minor, b: Minor) -> Result<Minor> {
    if a.is_zero() || b.is_zero() {
        return Ok(Minor::from_i128(0));
    }
    model.spend(80 + a.magnitude.bits().div_ceil(32) * b.magnitude.bits().div_ceil(32))?;
    a.multiply_wide(b).map_err(|_| arithmetic())
}

fn exact(model: &mut Model, numerator: Minor, denominator: Minor) -> Result<Minor> {
    if numerator.is_zero() {
        return Ok(numerator);
    }
    numerator
        .divide_exact(model, denominator)?
        .ok_or_else(|| certificate("fraction-free matrix elimination divides exactly"))
}

fn lower_fraction(model: &mut Model, numerator: Minor, denominator: Minor) -> Result<Rational> {
    if numerator.is_zero() {
        return Ok(Rational::ZERO);
    }
    let divisor = numerator.magnitude.gcd(model, denominator.magnitude)?;
    let divisor = Minor {
        magnitude: divisor,
        negative: denominator.negative,
    };
    let a = exact(model, numerator, divisor)?
        .signed()
        .ok_or_else(arithmetic)?;
    let b = exact(model, denominator, divisor)?
        .signed()
        .ok_or_else(arithmetic)?;
    Rational::new(a, b)
}

pub(super) fn lcm(a: i128, b: i128) -> Result<i128> {
    mul(a / gcd(a as u128, b as u128) as i128, b)
}

pub(super) fn clear(values: &[Rational]) -> Result<(Vec<i128>, i128)> {
    let denominator = values.iter().try_fold(1, |d, v| lcm(d, v.denominator))?;
    let numerators = values
        .iter()
        .map(|v| mul(v.numerator, denominator / v.denominator))
        .collect::<Result<_>>()?;
    Ok((numerators, denominator))
}

pub(super) fn add(a: i128, b: i128) -> Result<i128> {
    a.checked_add(b).ok_or_else(arithmetic)
}

pub(super) fn mul(a: i128, b: i128) -> Result<i128> {
    a.checked_mul(b).ok_or_else(arithmetic)
}

pub(super) fn negate(a: i128) -> Result<i128> {
    a.checked_neg().ok_or_else(arithmetic)
}
