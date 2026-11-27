use std::collections::{BTreeMap, HashMap};

use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_fusion::{self, Grade, Support, Weight};
use crate::rv_ideal::Natural;

const MAX_RANK: usize = 64;
const MAX_SUPPORT: usize = 4096;
const MAX_DIGITS: usize = 65_536;

type Matrix<T> = Vec<Vec<T>>;

#[derive(Clone, Copy)]
pub(crate) enum Report {
    Series,
    Numerator,
    Denominator,
}

struct Handle {
    weights: Vec<Weight>,
    matrix: Matrix<i128>,
}

fn arithmetic() -> Error {
    Error::plain("conformal series integer arithmetic is out of range")
}

fn normalize(
    model: &mut Model,
    grade: Grade,
    labels: Vec<Vec<i128>>,
) -> Result<Option<Vec<Weight>>> {
    let mut degree = 0;
    let mut weights = Vec::with_capacity(labels.len());
    for label in labels {
        let normalized = grade.normalize(model, label)?;
        degree = (degree + normalized.iter().sum::<i128>()) % grade.rank as i128;
        weights.push(normalized);
    }
    Ok((degree == 0).then_some(weights))
}

fn vacuum(grade: Grade) -> Weight {
    vec![0; grade.rank].into()
}

impl Handle {
    fn build(model: &mut Model, grade: Grade) -> Result<Self> {
        let handle = rv_fusion::handle_support(model, grade)?;
        let vacuum = vacuum(grade);
        let mut positions = HashMap::from([(vacuum.clone(), 0)]);
        let mut weights = vec![vacuum];
        let mut columns = Vec::new();
        while columns.len() < weights.len() {
            let input: Support = [(weights[columns.len()].clone(), 1)].into();
            let mapped = rv_fusion::multiply_support(model, grade, &handle, &input)?;
            let mut column = Vec::new();
            for (weight, scalar) in mapped.iter() {
                model.spend(1 + weight.len())?;
                let i = if let Some(&i) = positions.get(weight) {
                    i
                } else {
                    let i = weights.len();
                    if i == MAX_RANK {
                        return Err(Error::plain(
                            "conformal handle orbit exceeds the limit of 64 weights",
                        ));
                    }
                    positions.insert(weight.clone(), i);
                    weights.push(weight.clone());
                    i
                };
                column.push((i, *scalar));
            }
            columns.push(column);
        }
        let rank = weights.len();
        let mut matrix = vec![vec![0; rank]; rank];
        for (j, column) in columns.into_iter().enumerate() {
            for (i, coefficient) in column {
                matrix[i][j] = coefficient;
            }
        }
        Ok(Self { weights, matrix })
    }

    fn row(&self, model: &mut Model, product: &Support) -> Result<Vec<i128>> {
        self.weights
            .iter()
            .map(|weight| {
                model.spend(weight.len() * (1 + product.len().max(1).ilog2() as usize))?;
                let conjugate = rv_fusion::dual_weight(weight);
                Ok(product
                    .binary_search_by(|(key, _)| key.cmp(&conjugate))
                    .map(|i| product[i].1)
                    .unwrap_or(0))
            })
            .collect()
    }
}

fn checked_dot(a: &[i128], b: impl Iterator<Item = i128>) -> Result<i128> {
    a.iter().zip(b).try_fold(0i128, |total, (&x, y)| {
        x.checked_mul(y)
            .and_then(|term| total.checked_add(term))
            .ok_or_else(arithmetic)
    })
}

fn characteristic(model: &mut Model, matrix: &Matrix<i128>) -> Result<Vec<i128>> {
    let rank = matrix.len();
    let mut traces = Vec::with_capacity(rank);
    let mut power = matrix.clone();
    for degree in 1..=rank {
        model.spend(rank)?;
        let trace = power
            .iter()
            .enumerate()
            .try_fold(0i128, |trace, (i, row)| {
                trace.checked_add(row[i]).ok_or_else(arithmetic)
            })?;
        traces.push(trace);
        if degree < rank {
            model.spend(rank * rank * rank)?;
            for row in &mut power {
                let next = (0..rank)
                    .map(|j| checked_dot(row, matrix.iter().map(|x| x[j])))
                    .collect::<Result<_>>()?;
                *row = next;
            }
        }
    }
    // Integer Newton identities compute det(I-tM). The division by the
    // index is exact in Z and takes place only during the exact query.
    let mut denominator = vec![1];
    for degree in 1..=rank {
        model.spend(degree)?;
        let signed = checked_dot(&denominator, traces[..degree].iter().rev().copied())?;
        let coefficient = signed.checked_neg().ok_or_else(arithmetic)?;
        let divisor = degree as i128;
        if coefficient % divisor != 0 {
            return Err(Error::plain(
                "conformal characteristic division was not exact",
            ));
        }
        denominator.push(coefficient / divisor);
    }
    trim(&mut denominator);
    Ok(denominator)
}

fn trim(poly: &mut Vec<i128>) {
    poly.truncate(poly.iter().rposition(|&c| c != 0).map_or(0, |i| i + 1));
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn divide_unsigned(value: i128, divisor: u128) -> i128 {
    if let Ok(divisor) = i128::try_from(divisor) {
        value / divisor
    } else {
        -i128::from(value == i128::MIN)
    }
}

fn primitive(poly: &mut [i128]) {
    let content = poly.iter().map(|c| c.unsigned_abs()).fold(0, gcd);
    if content > 1 {
        for c in poly {
            *c = divide_unsigned(*c, content);
        }
    }
}

fn polynomial_gcd(model: &mut Model, mut a: Vec<i128>, mut b: Vec<i128>) -> Result<Vec<i128>> {
    primitive(&mut a);
    primitive(&mut b);
    while !b.is_empty() {
        while a.len() >= b.len() {
            model.spend(a.len())?;
            let x = *a.last().unwrap();
            let y = *b.last().unwrap();
            let factor = gcd(x.unsigned_abs(), y.unsigned_abs());
            let x = divide_unsigned(x, factor);
            let y = divide_unsigned(y, factor);
            let offset = a.len() - b.len();
            for (i, value) in a.iter_mut().enumerate() {
                let term = i.checked_sub(offset).map_or(0, |j| b[j]);
                *value = value
                    .checked_mul(y)
                    .and_then(|v| term.checked_mul(x).and_then(|w| v.checked_sub(w)))
                    .ok_or_else(arithmetic)?;
            }
            trim(&mut a);
            primitive(&mut a);
        }
        (a, b) = (b, a);
    }
    if a.first().is_some_and(|value| *value < 0) {
        for value in &mut a {
            *value = value.checked_neg().ok_or_else(arithmetic)?;
        }
    }
    Ok(a)
}

fn exact_divide(model: &mut Model, mut poly: Vec<i128>, divisor: &[i128]) -> Result<Vec<i128>> {
    let mut result = vec![0; poly.len().saturating_sub(divisor.len() - 1)];
    let leading = *divisor.last().unwrap();
    while poly.len() >= divisor.len() {
        model.spend(divisor.len())?;
        let top = *poly.last().unwrap();
        let coefficient = top.checked_div(leading).ok_or_else(arithmetic)?;
        let offset = poly.len() - divisor.len();
        result[offset] = coefficient;
        for (value, &subtract) in poly[offset..].iter_mut().zip(divisor) {
            *value = coefficient
                .checked_mul(subtract)
                .and_then(|term| value.checked_sub(term))
                .ok_or_else(arithmetic)?;
        }
        if poly.last() != Some(&0) {
            return Err(Error::plain(
                "conformal series polynomial division was not exact",
            ));
        }
        trim(&mut poly);
    }
    if !poly.is_empty() {
        return Err(Error::plain(
            "conformal series polynomial division was not exact",
        ));
    }
    trim(&mut result);
    Ok(result)
}

fn rational(
    model: &mut Model,
    grade: Grade,
    punctures: Vec<Vec<i128>>,
) -> Result<(Vec<i128>, Vec<i128>)> {
    let Some(punctures) = normalize(model, grade, punctures)? else {
        return Ok((vec![], vec![1]));
    };
    let product = rv_fusion::puncture_product(model, grade, punctures)?;
    let handle = Handle::build(model, grade)?;
    let row = handle.row(model, &product)?;
    if row.iter().all(|&c| c == 0) {
        return Ok((vec![], vec![1]));
    }
    let rank = row.len();
    let denominator = characteristic(model, &handle.matrix)?;
    let mut state = vec![0; rank];
    state[0] = 1;
    let mut moments = Vec::with_capacity(rank);
    for degree in 0..rank {
        model.spend(rank)?;
        moments.push(checked_dot(&row, state.iter().copied())?);
        if degree + 1 < rank {
            model.spend(rank * rank)?;
            state = handle
                .matrix
                .iter()
                .map(|row| checked_dot(row, state.iter().copied()))
                .collect::<Result<_>>()?;
        }
    }
    let mut numerator = Vec::with_capacity(rank);
    for j in 0..rank {
        let len = (j + 1).min(denominator.len());
        model.spend(len)?;
        numerator.push(checked_dot(
            &denominator[..len],
            moments[j + 1 - len..=j].iter().rev().copied(),
        )?);
    }
    trim(&mut numerator);
    if numerator.is_empty() {
        return Ok((vec![], vec![1]));
    }
    // Cayley–Hamilton proves this numerator accounts for every subsequent
    // coefficient, not just the first moments; exact gcd removes all poles
    // absent from this particular vacuum pairing.
    let divisor = polynomial_gcd(model, denominator.clone(), numerator.clone())?;
    if divisor == [1] {
        return Ok((numerator, denominator));
    }
    let numerator = exact_divide(model, numerator, &divisor)?;
    let denominator = exact_divide(model, denominator, &divisor)?;
    if denominator.first() != Some(&1) {
        return Err(Error::plain(
            "conformal series denominator lost its unit constant",
        ));
    }
    Ok((numerator, denominator))
}

pub(crate) fn report(
    model: &mut Model,
    grade: Grade,
    punctures: Vec<Vec<i128>>,
    report: Report,
) -> Result<String> {
    let (p, q) = rational(model, grade, punctures)?;
    Ok(match report {
        Report::Series => format!("rvConformalSeries = {p:?} / {q:?}"),
        Report::Numerator => format!("rvConformalNum = {p:?}"),
        Report::Denominator => format!("rvConformalDen = {q:?}"),
    })
}

fn mod_add(a: u128, b: u128, m: u128) -> u128 {
    if a >= m - b {
        a - (m - b)
    } else {
        a + b
    }
}

fn mod_mul(model: &mut Model, mut a: u128, mut b: u128, m: u128) -> Result<u128> {
    if let Some(product) = a.checked_mul(b) {
        return Ok(product % m);
    }
    let mut result = 0;
    while b != 0 {
        model.spend(1)?;
        if b & 1 != 0 {
            result = mod_add(result, a, m);
        }
        b >>= 1;
        if b != 0 {
            a = mod_add(a, a, m);
        }
    }
    Ok(result)
}

fn mod_dot(model: &mut Model, a: &[u128], b: impl Iterator<Item = u128>, m: u128) -> Result<u128> {
    model.spend(a.len())?;
    let mut total = 0;
    for (&x, y) in a.iter().zip(b) {
        if x != 0 && y != 0 {
            total = mod_add(total, mod_mul(model, x, y, m)?, m);
        }
    }
    Ok(total)
}

fn mod_punctures(
    model: &mut Model,
    grade: Grade,
    punctures: Vec<Weight>,
    m: u128,
) -> Result<BTreeMap<Weight, u128>> {
    let mut product = BTreeMap::from([(vacuum(grade), 1)]);
    for weight in punctures {
        let mut next = BTreeMap::new();
        for (state, coefficient) in product {
            let image = rv_fusion::basis_product_support(model, grade, &weight, &state)?;
            for (target, multiplicity) in image.iter() {
                model.spend(1 + target.len())?;
                let contribution = mod_mul(model, coefficient, *multiplicity as u128 % m, m)?;
                if contribution == 0 {
                    continue;
                }
                let entry = next.entry(target.clone()).or_insert(0);
                *entry = mod_add(*entry, contribution, m);
                if next.len() > MAX_SUPPORT {
                    return Err(Error::plain(
                        "modular conformal puncture support exceeds the limit of 4096 weights",
                    ));
                }
            }
        }
        next.retain(|_, c| *c != 0);
        product = next;
    }
    Ok(product)
}

pub(crate) fn modular(
    model: &mut Model,
    grade: Grade,
    mut genus: u128,
    m: u128,
    punctures: Vec<Vec<i128>>,
) -> Result<u128> {
    if m == 0 {
        return Err(Error::plain("conformal modulus must be positive"));
    }
    let Some(punctures) = normalize(model, grade, punctures)? else {
        return Ok(0);
    };
    if m == 1 {
        return Ok(0);
    }
    if grade.level == 0 {
        return Ok(1);
    }
    let product = mod_punctures(model, grade, punctures, m)?;
    if genus == 0 {
        return Ok(product.get(&vacuum(grade)).copied().unwrap_or(0));
    }
    if product.is_empty() {
        return Ok(0);
    }
    let handle = Handle::build(model, grade)?;
    let row: Vec<_> = handle
        .weights
        .iter()
        .map(|weight| {
            product
                .get(&rv_fusion::dual_weight(weight))
                .copied()
                .unwrap_or(0)
        })
        .collect();
    if row.iter().all(|&c| c == 0) {
        return Ok(0);
    }
    let rank = row.len();
    let mut state = vec![0; rank];
    state[0] = 1;
    let mut matrix: Matrix<u128> = handle
        .matrix
        .iter()
        .map(|row| row.iter().map(|&value| value as u128 % m).collect())
        .collect();
    loop {
        if genus & 1 != 0 {
            state = matrix
                .iter()
                .map(|row| mod_dot(model, row, state.iter().copied(), m))
                .collect::<Result<_>>()?;
        }
        genus >>= 1;
        if genus == 0 {
            break;
        }
        let mut square = Vec::with_capacity(rank);
        for row in &matrix {
            square.push(
                (0..rank)
                    .map(|j| mod_dot(model, row, matrix.iter().map(|r| r[j]), m))
                    .collect::<Result<_>>()?,
            );
        }
        matrix = square;
    }
    mod_dot(model, &row, state.into_iter(), m)
}

pub(crate) fn exact(
    model: &mut Model,
    grade: Grade,
    genus: u128,
    punctures: Vec<Vec<i128>>,
) -> Result<String> {
    let Some(punctures) = normalize(model, grade, punctures)? else {
        return Ok("0".into());
    };
    if grade.level == 0 {
        return Ok("1".into());
    }
    let product = rv_fusion::puncture_product(model, grade, punctures)?;
    let handle = Handle::build(model, grade)?;
    let row = handle.row(model, &product)?;
    if row.iter().all(|&c| c == 0) {
        return Ok("0".into());
    }
    let rank = row.len();
    // Every reachable weight has a positive diagonal entry at least two
    // (the handle contains all vacuum copies). A positive pairing then grows
    // at least as 2^(g-rank), bounding the possible shown decimal size.
    if genus > 4 * MAX_DIGITS as u128 + rank as u128 {
        return Err(Error::plain(
            "conformal dimension exceeds the decimal output limit",
        ));
    }
    let mut state: Vec<_> = (0..rank).map(|_| Natural::default()).collect();
    state[0] = Natural::unit();
    for _ in 0..genus {
        let mut next = Vec::with_capacity(rank);
        for row in &handle.matrix {
            model.spend(rank)?;
            let mut value = Natural::default();
            for (&coefficient, entry) in row.iter().zip(&state) {
                if !entry.is_zero() {
                    value.add_multiple(model, entry, coefficient as u128)?;
                }
            }
            if value.decimal_bound() > MAX_DIGITS + 8 {
                return Err(Error::plain(
                    "conformal dimension exceeds the decimal output limit",
                ));
            }
            next.push(value);
        }
        state = next;
    }
    let mut result = Natural::default();
    for (&c, entry) in row.iter().zip(&state) {
        if !entry.is_zero() {
            result.add_multiple(model, entry, c as u128)?;
        }
    }
    let result = result.decimal();
    if result.len() > MAX_DIGITS {
        return Err(Error::plain(
            "conformal dimension exceeds the decimal output limit",
        ));
    }
    Ok(result)
}
