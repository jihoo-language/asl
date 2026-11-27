use std::collections::{BTreeMap, HashMap};

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_ore::{self, arithmetic, gcd, Powers};
use crate::rv_qgr::{self, Grade, Schubert, Shape, Support};

const MAX_RANK: usize = 64;
const MAX_POLYNOMIAL_CELLS: usize = 300_000;
const MAX_PERIOD: usize = 4096;

mod integer_gcd;

pub(crate) type Polynomial = Vec<i128>;

pub(crate) fn storage(size: usize) -> Result<()> {
    if size > MAX_POLYNOMIAL_CELLS {
        return Err(Error::plain(
            "quantum differential calculation exceeds the polynomial storage limit",
        ));
    }
    Ok(())
}

fn derivative_power(q: usize, c: i128) -> Result<i128> {
    (q as i128).checked_mul(c).ok_or_else(arithmetic)
}

fn nabla_support(
    model: &mut Model,
    grade: Grade,
    input: &[(Schubert, i128)],
) -> Result<Vec<(Schubert, i128)>> {
    let mut sum = rv_qgr::Sum::default();
    for (term, c) in input.iter() {
        sum.insert(model, term.clone(), derivative_power(term.quantum, *c)?, 0)?;
        for mut next in rv_qgr::hyperplane(model, grade, &term.shape, sum.cells)? {
            next.quantum = next
                .quantum
                .checked_add(term.quantum)
                .ok_or_else(arithmetic)?;
            sum.insert(model, next, *c, 0)?;
        }
    }
    Ok(sum.terms.into_iter().collect())
}

pub(crate) fn nabla(model: &mut Model, grade: Grade, input: &Support) -> Result<Value> {
    let support = nabla_support(model, grade, input)?;
    rv_qgr::persist(model, support)
}

pub(crate) fn apply(
    model: &mut Model,
    grade: Grade,
    op: &rv_ore::Support,
    input: &Support,
) -> Result<Value> {
    if input.is_empty() || op.is_empty() {
        return Ok(rv_qgr::zero());
    }
    let mut by_order: BTreeMap<usize, Vec<(usize, i128)>> = BTreeMap::new();
    for &(Powers { q, theta }, c) in op.iter() {
        model.spend(1 + by_order.len().max(1).ilog2() as usize)?;
        by_order.entry(theta).or_default().push((q, c));
    }
    model.spend(*by_order.last_key_value().unwrap().0)?;
    let mut sum = rv_qgr::Sum::default();
    let mut current = input.clone();
    let mut previous = 0;
    for (order, factors) in by_order {
        for _ in previous..order {
            current = nabla_support(model, grade, &current)?.into();
            if current.is_empty() {
                return rv_qgr::persist(model, sum.terms.into_iter().collect());
            }
        }
        for (power, scalar) in factors {
            for (term, c) in current.iter() {
                let quantum = term.quantum.checked_add(power).ok_or_else(arithmetic)?;
                let term = Schubert {
                    quantum,
                    shape: term.shape.clone(),
                };
                sum.insert(
                    model,
                    term,
                    c.checked_mul(scalar).ok_or_else(arithmetic)?,
                    0,
                )?;
            }
        }
        previous = order;
    }
    rv_qgr::persist(model, sum.terms.into_iter().collect())
}

pub(crate) fn schubert_basis(model: &mut Model, grade: Grade, limit: usize) -> Result<Vec<Shape>> {
    let mut rank = 1u128;
    for index in 1..=grade.rank {
        rank = rank * (grade.width() as u128 + index as u128) / index as u128;
        if rank > limit as u128 {
            return Err(Error::plain(format!(
                "quantum connection Schubert rank exceeds the limit of {limit}"
            )));
        }
    }
    model.spend(rank as usize * grade.rank)?;
    let mut shapes: Vec<Shape> = Vec::with_capacity(rank as usize);
    let mut beta: Vec<_> = (0..grade.rank).collect();
    loop {
        let descending: Vec<_> = beta.iter().rev().copied().collect();
        shapes.push(rv_qgr::beta_shape(&descending));
        let Some(index) = (0..grade.rank).rev().find(|&i| beta[i] < grade.width() + i) else {
            break;
        };
        beta[index] += 1;
        for next in index + 1..grade.rank {
            beta[next] = beta[next - 1] + 1;
        }
    }
    Ok(shapes)
}

struct Connection {
    arrows: Vec<Vec<(usize, usize)>>,
    weights: Vec<usize>,
    high_to_low: Vec<usize>,
}

impl Connection {
    fn build(model: &mut Model, mut grade: Grade) -> Result<Self> {
        if grade.rank > grade.width() {
            grade = grade.complement();
        }
        let shapes = schubert_basis(model, grade, MAX_RANK)?;
        let positions: HashMap<_, _> = shapes
            .iter()
            .enumerate()
            .map(|(i, shape)| (shape.clone(), i))
            .collect();
        let mut arrows = Vec::with_capacity(shapes.len());
        for shape in &shapes {
            let action = rv_qgr::hyperplane(model, grade, shape, 0)?;
            let action = action
                .into_iter()
                .map(|term| (positions[&term.shape], term.quantum))
                .collect();
            arrows.push(action);
        }
        let weights: Vec<usize> = shapes.iter().map(|shape| shape.iter().sum()).collect();
        let mut high_to_low: Vec<_> = (0..shapes.len()).collect();
        high_to_low.sort_unstable_by_key(|&i| std::cmp::Reverse(weights[i]));
        Ok(Self {
            arrows,
            weights,
            high_to_low,
        })
    }

    fn next(&self, model: &mut Model, vector: &[Polynomial]) -> Result<Vec<Polynomial>> {
        let mut result = vec![Vec::new(); vector.len()];
        for (i, polynomial) in vector.iter().enumerate() {
            model.spend(polynomial.len() * (1 + self.arrows[i].len()))?;
            for (q, &coefficient) in polynomial.iter().enumerate() {
                add_at(&mut result[i], q, derivative_power(q, coefficient)?)?;
                for &(j, degree) in &self.arrows[i] {
                    add_at(&mut result[j], q + degree, coefficient)?;
                }
            }
        }
        for polynomial in &mut result {
            trim(polynomial);
        }
        Ok(result)
    }
}

fn trim(poly: &mut Polynomial) {
    poly.truncate(
        poly.iter()
            .rposition(|&entry| entry != 0)
            .map_or(0, |i| i + 1),
    );
}

fn add_at(poly: &mut Polynomial, exponent: usize, coefficient: i128) -> Result<()> {
    if coefficient != 0 {
        if exponent >= poly.len() {
            storage(exponent + 1)?;
            poly.resize(exponent + 1, 0);
        }
        poly[exponent] = poly[exponent]
            .checked_add(coefficient)
            .ok_or_else(arithmetic)?;
    }
    Ok(())
}

fn product(model: &mut Model, left: &[i128], right: &[i128]) -> Result<Polynomial> {
    if left.is_empty() || right.is_empty() {
        return Ok(Vec::new());
    }
    model.spend(left.len() * right.len())?;
    let size = left.len() + right.len() - 1;
    storage(size)?;
    let mut result = vec![0i128; size];
    for (i, &a) in left.iter().enumerate() {
        if a == 0 {
            continue;
        }
        for (j, &b) in right.iter().enumerate() {
            result[i + j] = a
                .checked_mul(b)
                .and_then(|ab| result[i + j].checked_add(ab))
                .ok_or_else(arithmetic)?;
        }
    }
    trim(&mut result);
    Ok(result)
}

fn subtract(model: &mut Model, mut left: Polynomial, right: &[i128]) -> Result<Polynomial> {
    model.spend(1 + left.len().max(right.len()))?;
    left.resize(left.len().max(right.len()), 0);
    for (entry, &value) in left.iter_mut().zip(right) {
        *entry = entry.checked_sub(value).ok_or_else(arithmetic)?;
    }
    trim(&mut left);
    Ok(left)
}

enum Division {
    Exact(Polynomial),
    NotExact,
    OutOfRange,
}

fn checked_divide(model: &mut Model, mut value: Polynomial, divisor: &[i128]) -> Result<Division> {
    if divisor == [1] || value.is_empty() {
        return Ok(Division::Exact(value));
    }
    let mut quotient = vec![0; value.len().saturating_sub(divisor.len() - 1)];
    let lead = *divisor.last().unwrap();
    while value.len() >= divisor.len() {
        let current = *value.last().unwrap();
        let Some(coefficient) = current.checked_div(lead) else {
            return Ok(Division::OutOfRange);
        };
        if current % lead != 0 {
            break;
        }
        let shift = value.len() - divisor.len();
        quotient[shift] = coefficient;
        model.spend(divisor.len())?;
        for (entry, &a) in value[shift..].iter_mut().zip(divisor) {
            let Some(next) = coefficient
                .checked_mul(a)
                .and_then(|term| entry.checked_sub(term))
            else {
                return Ok(Division::OutOfRange);
            };
            *entry = next;
        }
        trim(&mut value);
    }
    if !value.is_empty() {
        return Ok(Division::NotExact);
    }
    trim(&mut quotient);
    Ok(Division::Exact(quotient))
}

pub(crate) fn exact_divide(
    model: &mut Model,
    value: Polynomial,
    divisor: &[i128],
) -> Result<Polynomial> {
    if divisor == [1] || value.is_empty() {
        return Ok(value);
    }
    let result = match checked_divide(model, value.clone(), divisor)? {
        Division::OutOfRange => integer_gcd::wide_divide(model, &value, divisor)?,
        result => result,
    };
    match result {
        Division::Exact(quotient) => Ok(quotient),
        Division::NotExact => Err(Error::plain(
            "quantum differential fraction-free division was not exact",
        )),
        Division::OutOfRange => Err(arithmetic()),
    }
}

fn checked_cells_of(polys: impl IntoIterator<Item = usize>) -> Result<usize> {
    let mut cells = 0;
    for length in polys {
        cells += 1 + 2 * length;
        storage(cells)?;
    }
    Ok(cells)
}

fn pivot(
    model: &mut Model,
    matrix: &mut [Vec<Polynomial>],
    j: usize,
    previous: &[i128],
    cells: &mut usize,
) -> Result<()> {
    let (upper, lower) = matrix.split_at_mut(j + 1);
    let row = &upper[j];
    for target in lower {
        if target[j].is_empty() && row[j] == previous {
            continue;
        }
        for k in j + 1..row.len() {
            let first = product(model, &row[j], &target[k])?;
            let second = product(model, &target[j], &row[k])?;
            let numerator = subtract(model, first, &second)?;
            let next = exact_divide(model, numerator, previous)?;
            *cells = *cells - 2 * target[k].len() + 2 * next.len();
            storage(*cells)?;
            target[k] = next;
        }
        *cells -= 2 * target[j].len();
        target[j].clear();
    }
    Ok(())
}

fn null_relation(model: &mut Model, mut matrix: Vec<Vec<Polynomial>>) -> Result<Vec<Polynomial>> {
    let rank = matrix.len();
    let mut cells = checked_cells_of(matrix.iter().flatten().map(Vec::len))?;
    let mut lead = vec![1];
    let mut order = 0;
    while order < rank {
        let Some(row) = (order..rank).find(|&row| !matrix[row][order].is_empty()) else {
            break;
        };
        matrix.swap(order, row);
        pivot(model, &mut matrix, order, &lead, &mut cells)?;
        lead = matrix[order][order].clone();
        order += 1;
    }
    // The last pivot is the selected square minor. Cramer's rule gives
    // polynomial coefficients with that minor on the last column.
    let mut relation = vec![Vec::new(); order + 1];
    relation[order] = lead;
    for i in (0..order).rev() {
        let mut numerator = Vec::new();
        for j in i + 1..=order {
            let term = product(model, &matrix[i][j], &relation[j])?;
            numerator = subtract(model, numerator, &term)?;
        }
        relation[i] = exact_divide(model, numerator, &matrix[i][i])?;
    }
    primitive_relation(model, &mut relation)?;
    Ok(relation)
}

fn divide_content(poly: &mut Polynomial, content: u128) {
    if content <= 1 {
        return;
    }
    if let Ok(divisor) = i128::try_from(content) {
        for entry in poly {
            *entry /= divisor;
        }
    } else {
        for entry in poly {
            *entry = -i128::from(*entry == i128::MIN);
        }
    }
}

fn content(poly: &[i128]) -> u128 {
    poly.iter().map(|x| x.unsigned_abs()).fold(0, gcd)
}

fn primitive(poly: &mut Polynomial) {
    divide_content(poly, content(poly));
}

fn primitive_gcd(
    model: &mut Model,
    mut a: Polynomial,
    mut b: Polynomial,
) -> Result<Option<Polynomial>> {
    while !b.is_empty() {
        while a.len() >= b.len() {
            let x = *a.last().unwrap();
            let y = *b.last().unwrap();
            let common = gcd(x.unsigned_abs(), y.unsigned_abs());
            let scale = |value| {
                if let Ok(common) = i128::try_from(common) {
                    value / common
                } else {
                    -1
                }
            };
            let (x, y) = (scale(x), scale(y));
            let offset = a.len() - b.len();
            model.spend(a.len())?;
            for (i, value) in a.iter_mut().enumerate() {
                let term = i.checked_sub(offset).map_or(0, |j| b[j]);
                let Some(next) = value
                    .checked_mul(y)
                    .and_then(|l| term.checked_mul(x).and_then(|r| l.checked_sub(r)))
                else {
                    return Ok(None);
                };
                *value = next;
            }
            trim(&mut a);
            primitive(&mut a);
        }
        (a, b) = (b, a);
    }
    if a.last().is_some_and(|entry| *entry < 0) {
        for entry in &mut a {
            *entry = entry.checked_neg().ok_or_else(arithmetic)?;
        }
    }
    Ok(Some(a))
}

pub(crate) fn polynomial_gcd(
    model: &mut Model,
    mut a: Polynomial,
    mut b: Polynomial,
) -> Result<Polynomial> {
    model.spend(a.len() + b.len())?;
    primitive(&mut a);
    primitive(&mut b);
    if a.len() == 1 || b.len() == 1 {
        return Ok(vec![1]);
    }
    if let Some(result) = primitive_gcd(model, a.clone(), b.clone())? {
        Ok(result)
    } else {
        integer_gcd::gcd_fallback(model, &a, &b)
    }
}

fn primitive_relation(model: &mut Model, relation: &mut [Polynomial]) -> Result<()> {
    let mut divisor = Vec::new();
    let integer = relation.iter().map(|p| content(p)).fold(0, gcd);
    for polynomial in relation.iter_mut() {
        divide_content(polynomial, integer);
        if !polynomial.is_empty() && divisor != [1] {
            divisor = polynomial_gcd(model, divisor, polynomial.clone())?;
        }
    }
    for polynomial in relation.iter_mut() {
        *polynomial = exact_divide(model, std::mem::take(polynomial), &divisor)?;
    }
    if relation
        .last()
        .and_then(|polynomial| polynomial.last())
        .is_some_and(|entry| *entry < 0)
    {
        for entry in relation.iter_mut().flatten() {
            *entry = entry.checked_neg().ok_or_else(arithmetic)?;
        }
    }
    Ok(())
}

pub(crate) fn annihilator(model: &mut Model, grade: Grade) -> Result<Value> {
    let mut result = rv_ore::Sum::default();
    if grade.rank.min(grade.width()) == 1 {
        result.add(
            model,
            Powers {
                q: 0,
                theta: grade.ambient,
            },
            1,
        )?;
        result.add(model, Powers { q: 1, theta: 0 }, -1)?;
        return result.finish(model);
    }
    let connection = Connection::build(model, grade)?;
    let rank = connection.arrows.len();
    let mut matrix = vec![Vec::with_capacity(rank + 1); rank];
    let mut vector = vec![Vec::new(); rank];
    vector[0] = vec![1];
    let mut cells = 0;
    for j in 0..=rank {
        cells += checked_cells_of(vector.iter().map(Vec::len))?;
        storage(cells)?;
        let next = if j < rank {
            connection.next(model, &vector)?
        } else {
            Vec::new()
        };
        for (row, polynomial) in matrix.iter_mut().zip(vector) {
            row.push(polynomial);
        }
        vector = next;
    }
    let relation = null_relation(model, matrix)?;
    for (theta, polynomial) in relation.into_iter().enumerate() {
        for (q, coefficient) in polynomial.into_iter().enumerate() {
            result.add(model, Powers { q, theta }, coefficient)?;
        }
    }
    result.finish(model)
}

#[derive(Clone, Copy)]
pub(crate) struct Rational {
    numerator: i128,
    denominator: i128,
}

impl Rational {
    pub(crate) const ZERO: Self = Self {
        numerator: 0,
        denominator: 1,
    };
    const ONE: Self = Self {
        numerator: 1,
        denominator: 1,
    };

    pub(crate) fn add(self, other: Self) -> Result<Self> {
        let common = gcd(self.denominator as u128, other.denominator as u128) as i128;
        let left = self.denominator / common;
        let right = other.denominator / common;
        let numerator = self
            .numerator
            .checked_mul(right)
            .and_then(|a| {
                other
                    .numerator
                    .checked_mul(left)
                    .and_then(|b| a.checked_add(b))
            })
            .ok_or_else(arithmetic)?;
        let reduce = gcd(numerator.unsigned_abs(), common as u128) as i128;
        let denominator = left
            .checked_mul(other.denominator / reduce)
            .ok_or_else(arithmetic)?;
        Ok(Self {
            numerator: numerator / reduce,
            denominator,
        })
    }

    pub(crate) fn is_zero(self) -> bool {
        self.numerator == 0
    }

    pub(crate) fn scale(self, numerator: i128, denominator: i128) -> Result<Self> {
        if self.is_zero() || numerator == 0 {
            return Ok(Self::ZERO);
        }
        let divisor = gcd(numerator.unsigned_abs(), denominator as u128) as i128;
        let numerator = numerator / divisor;
        let denominator = denominator / divisor;
        let cross_num = gcd(numerator.unsigned_abs(), self.denominator as u128) as i128;
        let cross_den = gcd(self.numerator.unsigned_abs(), denominator as u128) as i128;
        let numerator = (self.numerator / cross_den)
            .checked_mul(numerator / cross_num)
            .ok_or_else(arithmetic)?;
        let denominator = (self.denominator / cross_num)
            .checked_mul(denominator / cross_den)
            .ok_or_else(arithmetic)?;
        Ok(Self {
            numerator,
            denominator,
        })
    }

    fn divide(self, divisor: i128) -> Result<Self> {
        let reduce = gcd(self.numerator.unsigned_abs(), divisor as u128) as i128;
        let denominator = self
            .denominator
            .checked_mul(divisor / reduce)
            .ok_or_else(arithmetic)?;
        Ok(Self {
            numerator: self.numerator / reduce,
            denominator,
        })
    }
}

impl std::fmt::Display for Rational {
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(out, "{}", self.numerator)?;
        if self.denominator != 1 {
            write!(out, "/{}", self.denominator)?;
        }
        Ok(())
    }
}

pub(crate) fn check_period_degree(degree: usize) -> Result<()> {
    if degree > MAX_PERIOD {
        return Err(Error::plain(
            "quantum period degree exceeds the limit of 4096",
        ));
    }
    Ok(())
}

#[derive(Clone, Copy)]
pub(crate) enum PeriodScale {
    Standard,
    Factorial,
}

pub(crate) fn period_coefficients(
    model: &mut Model,
    grade: Grade,
    degree: usize,
    scale: PeriodScale,
) -> Result<Vec<Rational>> {
    check_period_degree(degree)?;
    let factorial = matches!(scale, PeriodScale::Factorial);
    let mut values = vec![Rational::ONE];
    if grade.rank.min(grade.width()) == 1 {
        if factorial {
            model.spend(degree)?;
            values.resize(degree + 1, Rational::ONE);
            return Ok(values);
        }
        for d in 1..=degree {
            model.spend(1)?;
            let power = if d == 1 {
                1
            } else {
                let ambient = u32::try_from(grade.ambient).map_err(|_| arithmetic())?;
                (d as i128).checked_pow(ambient).ok_or_else(arithmetic)?
            };
            values.push(values.last().copied().unwrap().divide(power)?);
        }
    } else if degree > 0 {
        let connection = Connection::build(model, grade)?;
        let mut previous = vec![Rational::ZERO; connection.arrows.len()];
        previous[0] = Rational::ONE;
        let mut quantum_targets = vec![false; previous.len()];
        for &(target, quantum) in connection.arrows.iter().flatten() {
            quantum_targets[target] |= quantum != 0;
        }
        for d in 1..=degree {
            if factorial {
                // For d>0, U_d[i]=(d!)^n l_d[i]/d^|i|. A classical arrow
                // raises |i| by one, so its contribution to U_d has weight 1.
                // A quantum arrow lowers |i| by n-1. It has weight
                // ((d-1)/d)^|j| multiplying U_(d-1)[j] at target j.
                // In degree zero only the weight-zero unit is nonzero.
                for (i, &target) in quantum_targets.iter().enumerate() {
                    if target && !previous[i].is_zero() {
                        let weight = connection.weights[i];
                        model.spend(weight)?;
                        for _ in 0..weight {
                            previous[i] = previous[i].scale((d - 1) as i128, d as i128)?;
                        }
                    }
                }
            }
            let mut current = vec![Rational::ZERO; previous.len()];
            // Classical arrows raise codimension, so these entries of the
            // current row are already known when solving for its i-th entry.
            for &i in &connection.high_to_low {
                let mut value = Rational::ZERO;
                model.spend(1 + connection.arrows[i].len())?;
                for &(j, quantum) in &connection.arrows[i] {
                    value = value.add(if quantum == 0 {
                        current[j]
                    } else {
                        previous[j]
                    })?;
                }
                current[i] = if factorial {
                    value
                } else {
                    value.divide(d as i128)?
                };
            }
            values.push(current[0]);
            previous = current;
        }
    }
    Ok(values)
}

pub(crate) fn period(model: &mut Model, grade: Grade, degree: usize) -> Result<String> {
    let values = period_coefficients(model, grade, degree, PeriodScale::Standard)?;
    Ok(format_period(&values))
}

pub(crate) fn regularized_period(model: &mut Model, grade: Grade, degree: usize) -> Result<String> {
    let values = period_coefficients(model, grade, degree, PeriodScale::Factorial)?;
    Ok(format_period(&values))
}

fn format_period(values: &[Rational]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(", ")
    )
}
