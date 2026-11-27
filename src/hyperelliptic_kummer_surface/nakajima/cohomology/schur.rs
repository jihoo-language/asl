use std::cmp::Ordering;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;

use super::{divide, factorial, positive, product, pure, Class, Denominator, Integer, Vector};

type Series = Vec<Class>;

pub(crate) fn character(
    model: &mut Model,
    line: &Vector,
    n: u128,
    degree: u128,
    shape: &[u128],
) -> Result<Class> {
    if shape.len() as u128 > n || degree.div_ceil(2) > n {
        return Ok(Class::zero());
    }
    if shape.is_empty() {
        return if degree == 0 {
            Class::unit(model, n)
        } else {
            Ok(Class::zero())
        };
    }
    if degree == 0 {
        return rank(model, n, shape);
    }
    if shape == [1] {
        return Class::character(model, line, n, degree);
    }
    let degree = usize::try_from(degree).map_err(|_| work_limit())?;
    guard_cells(degree.saturating_add(1).saturating_mul(16))?;
    let length = i128::try_from(n).map_err(|_| super::super::grade_range())?;

    // Remove the determinant twist before planning Jacobi–Trudi. This also
    // handles all powers of a line and full-height rectangles without a table.
    let shift = if shape.len() as u128 == n {
        *shape.last().unwrap()
    } else {
        0
    };
    let reduced: Vec<_> = shape
        .iter()
        .map(|row| row - shift)
        .take_while(|&row| row != 0)
        .collect();
    if reduced.is_empty() {
        let first = Class::chern(model, line, n, 1)?;
        return exponential(model, first, n, degree, shift).map(|mut series| series.pop().unwrap());
    }

    let (rows, exterior) = if reduced[0] < reduced.len() as u128 {
        let width = reduced[0] as usize;
        let transpose = (1..=width)
            .map(|column| reduced.partition_point(|&row| row >= column as u128) as u128)
            .collect();
        (transpose, true)
    } else {
        (reduced, false)
    };
    let size = rows.len();
    let subsets = u32::try_from(size)
        .ok()
        .and_then(|shift| 1usize.checked_shl(shift))
        .ok_or_else(work_limit)?;
    let bound = rows[0].saturating_add(size as u128 - 1);
    let bound =
        usize::try_from(if exterior { bound.min(n) } else { bound }).map_err(|_| work_limit())?;
    let series_cells = (degree + 1).saturating_mul(10);
    let units = subsets.saturating_add(3usize.saturating_mul(bound.saturating_add(1)));
    let scratch = units.saturating_mul(series_cells);
    guard_cells(scratch)?;
    let operations = bound.saturating_mul(bound.saturating_add(1)) / 2;
    let operations = operations.saturating_add(subsets.saturating_mul(size));
    let convolution = (degree + 1).saturating_mul(degree + 2) / 2;
    model
        .spend(operations.saturating_mul(convolution))
        .map_err(|_| work_limit())?;

    let first = if shift == 0 {
        None
    } else {
        Some(Class::chern(model, line, n, 1)?)
    };
    let mut coeffs = Vec::with_capacity(degree + 1);
    coeffs.push(Class::character(model, line, n, 0)?);
    for (index, power) in Class::power_sums(model, line, n, degree)?
        .into_iter()
        .enumerate()
    {
        coeffs.push(power.divide_by_factorial(model, index as u128 + 1)?);
    }
    let mut cells = scratch;
    for class in &coeffs {
        charge(&mut cells, class.cells())?;
    }
    let table = symmetric_series(model, &coeffs, n, bound, exterior, &mut cells)?;
    let mut result = determinant(model, &rows, &table, length, degree, &mut cells)?;
    if let Some(first) = first {
        let twist = exponential(model, first, n, degree, shift)?;
        let mut top = Class::zero();
        for (a, x) in twist.iter().enumerate() {
            add_product(model, &mut top, x, &result[degree - a], length, false)?;
        }
        return Ok(top);
    }
    Ok(result.pop().unwrap())
}

fn rank(model: &mut Model, n: u128, shape: &[u128]) -> Result<Class> {
    let denominator = factorial(model, n)?;
    let size = usize::try_from(n).map_err(|_| work_limit())?;
    let mut numerator = Denominator::from_small(1);
    let mut divisor = Denominator::from_small(1);
    for right in 1..size {
        let right_row = shape.get(right).copied().unwrap_or(0);
        for (left, &left_row) in shape[..shape.len().min(right)].iter().enumerate() {
            model.spend(1)?;
            let offset = (right - left) as u128;
            if left_row == right_row {
                continue;
            }
            let mut factor = positive(left_row - right_row)
                .add(positive(offset))?
                .magnitude;
            let mut base = Denominator::from_small(offset);
            cancel(model, &mut factor, &mut base)?;
            cancel(model, &mut numerator, &mut base)?;
            cancel(model, &mut factor, &mut divisor)?;
            numerator = numerator
                .checked_product(factor)
                .ok_or_else(arithmetic_range)?;
            divisor = divisor.checked_product(base).ok_or_else(arithmetic_range)?;
        }
        debug_assert_eq!(divisor.small(), Some(1));
    }
    let rank = Integer {
        magnitude: numerator,
        negative: false,
    };
    Class::scaled_wide(model, &pure(n, 0), rank, denominator)
}

fn cancel(
    model: &mut Model,
    numerator: &mut Denominator,
    denominator: &mut Denominator,
) -> Result<()> {
    if denominator.small() == Some(1) {
        return Ok(());
    }
    let common = numerator.gcd(model, *denominator)?;
    if common.small() != Some(1) {
        *numerator = divide(model, *numerator, common)?;
        *denominator = divide(model, *denominator, common)?;
    }
    Ok(())
}

fn exponential(
    model: &mut Model,
    first: Class,
    n: u128,
    degree: usize,
    shift: u128,
) -> Result<Series> {
    let length = i128::try_from(n).map_err(|_| super::super::grade_range())?;
    let mut series = vec![Class::unit(model, n)?];
    let mut cells = (degree + 1) * 8;
    for index in 1..=degree {
        let power = series[index - 1].cup(model, &first, length)?;
        let term = Class::scaled_wide(
            model,
            &power.numerator,
            positive(shift),
            product(power.denominator, index as u128),
        )?;
        charge(&mut cells, term.cells())?;
        series.push(term);
    }
    Ok(series)
}

fn symmetric_series(
    model: &mut Model,
    ch: &[Class],
    n: u128,
    bound: usize,
    exterior: bool,
    cells: &mut usize,
) -> Result<Vec<Series>> {
    let length = i128::try_from(n).map_err(|_| super::super::grade_range())?;
    let degree = ch.len() - 1;
    let mut adams = Vec::with_capacity(bound);
    for power in 1..=bound {
        let mut row = Vec::with_capacity(ch.len());
        let mut scalar = Some(Denominator::from_small(1));
        for component in ch {
            let result = if is_zero(component) {
                Class::zero()
            } else {
                let weight = Integer {
                    magnitude: scalar.ok_or_else(arithmetic_range)?,
                    negative: false,
                };
                Class::scaled_wide(
                    model,
                    &component.numerator,
                    weight,
                    Denominator::from_small(component.denominator),
                )?
            };
            charge(cells, result.cells())?;
            row.push(result);
            scalar = scalar.and_then(|value| value.multiply_u128(power as u128));
        }
        adams.push(row);
    }
    let mut table = Vec::with_capacity(bound + 1);
    table.push(unit_series(model, n, degree)?);
    for index in 1..=bound {
        let mut sum = zeros(degree);
        for k in 1..=index {
            add_series_product(
                model,
                &mut sum,
                &adams[k - 1],
                &table[index - k],
                length,
                exterior && k.is_multiple_of(2),
            )?;
        }
        let row = sum
            .into_iter()
            .map(|coefficient| coefficient.scale(model, 1, index as u128))
            .collect::<Result<Vec<_>>>()?;
        charge(cells, row.iter().map(Class::cells).sum())?;
        table.push(row);
    }
    Ok(table)
}

fn determinant(
    model: &mut Model,
    rows: &[u128],
    table: &[Series],
    length: i128,
    degree: usize,
    cells: &mut usize,
) -> Result<Series> {
    let end = (1usize << rows.len()) - 1;
    let mut partials: Vec<Option<Series>> = (0..=end).map(|_| None).collect();
    let unit = unit_series(model, length as u128, degree)?;
    partials[0] = Some(unit);
    for mask in 0..end {
        let Some(partial) = partials[mask].take() else {
            continue;
        };
        let row = mask.count_ones() as usize;
        for column in 0..rows.len() {
            let bit = 1 << column;
            if mask & bit != 0 {
                continue;
            }
            let Some(index) = rows[row]
                .checked_add(column as u128)
                .and_then(|value| value.checked_sub(row as u128))
                .and_then(|value| usize::try_from(value).ok())
            else {
                continue;
            };
            let Some(entry) = table.get(index) else {
                continue;
            };
            let negative = !(mask >> (column + 1)).count_ones().is_multiple_of(2);
            let target = partials[mask | bit].get_or_insert_with(|| zeros(degree));
            let previous: usize = target.iter().map(Class::cells).sum();
            add_series_product(model, target, &partial, entry, length, negative)?;
            let current: usize = target.iter().map(Class::cells).sum();
            charge(cells, current.saturating_sub(previous))?;
        }
    }
    Ok(partials[end].take().unwrap_or_else(|| zeros(degree)))
}

fn add_series_product(
    model: &mut Model,
    target: &mut Series,
    left: &[Class],
    right: &[Class],
    length: i128,
    negative: bool,
) -> Result<()> {
    for (degree, sum) in target.iter_mut().enumerate() {
        for (i, x) in left.iter().take(degree + 1).enumerate() {
            add_product(model, sum, x, &right[degree - i], length, negative)?;
        }
    }
    Ok(())
}

fn add_product(
    model: &mut Model,
    sum: &mut Class,
    x: &Class,
    y: &Class,
    length: i128,
    negative: bool,
) -> Result<()> {
    if !is_zero(x) && !is_zero(y) {
        let value = x.cup(model, y, length)?;
        *sum = sum.add_signed(model, &value, negative)?;
    }
    Ok(())
}

fn is_zero(class: &Class) -> bool {
    class.numerator.terms.is_empty()
}

fn zeros(degree: usize) -> Series {
    (0..=degree).map(|_| Class::zero()).collect()
}

fn unit_series(model: &mut Model, n: u128, degree: usize) -> Result<Series> {
    let mut series = zeros(degree);
    series[0] = Class::unit(model, n)?;
    Ok(series)
}

pub(crate) fn exterior_euler(
    model: &mut Model,
    line: &Vector,
    twist: &Vector,
    n: u128,
    power: u128,
) -> Result<i128> {
    if power > n {
        return Ok(0);
    }
    if n == 0 {
        return Ok(1);
    }
    let base = half_square(model, twist)?.add(Integer::from_i128(2))?;
    let line_twist = base
        .add(half_square(model, line)?)?
        .add(line.wide_pair(model, twist)?.resize()?)?;
    let remaining = n - power;
    let symmetric_top = base.add(positive(remaining))?.add(Integer::from_i128(-1))?;
    // Detect zero factors before computing a potentially large other factor.
    let Some(exterior) = Binomial::new(line_twist, power)? else {
        return Ok(0);
    };
    let Some(symmetric) = Binomial::new(symmetric_top, remaining)? else {
        return Ok(0);
    };
    let abs = exterior
        .magnitude(model)?
        .checked_product(symmetric.magnitude(model)?)
        .ok_or_else(euler_range)?;
    abs.signed(exterior.negative != symmetric.negative)
        .ok_or_else(euler_range)
}

fn half_square(model: &mut Model, vector: &Vector) -> Result<Integer> {
    let pair = vector.wide_pair(model, vector)?;
    let (magnitude, remainder) = pair.magnitude.divide_u32(2);
    debug_assert_eq!(remainder, 0, "the K3 lattice is even");
    Ok(Integer {
        magnitude: magnitude.resize().expect("wide K3 lattice square"),
        negative: pair.negative,
    })
}

struct Binomial {
    high: Denominator,
    steps: u128,
    negative: bool,
}

impl Binomial {
    fn new(top: Integer, lower: u128) -> Result<Option<Self>> {
        if lower == 0 {
            return Ok(Some(Self {
                high: Denominator::from_small(1),
                steps: 0,
                negative: false,
            }));
        }
        let absolute = Denominator::from_small(lower);
        let negative = top.negative && !top.is_zero();
        let high = if negative {
            top.negate().add(positive(lower - 1))?.magnitude
        } else {
            top.magnitude
        };
        if high.cmp(absolute) == Ordering::Less {
            return Ok(None);
        }
        let complement = high.difference(absolute).small();
        let steps = complement.map_or(lower, |other| lower.min(other));
        Ok(Some(Self {
            high,
            steps,
            negative: negative && !lower.is_multiple_of(2),
        }))
    }

    fn magnitude(&self, model: &mut Model) -> Result<Denominator> {
        let one = Denominator::from_small(1);
        let mut value = one;
        let mut factor = self.high;
        for step in 1..=self.steps {
            model.spend(1)?;
            let product = value.checked_product(factor).ok_or_else(euler_range)?;
            let (quotient, remainder) = product.divide(model, step)?;
            debug_assert_eq!(remainder, 0);
            value = quotient;
            if value.bits() > 128 {
                return Err(euler_range());
            }
            factor = factor.difference(one);
        }
        Ok(value)
    }
}

fn charge(cells: &mut usize, additional: usize) -> Result<()> {
    *cells = cells.saturating_add(additional);
    guard_cells(*cells)
}

fn guard_cells(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        Err(work_limit())
    } else {
        Ok(())
    }
}

fn work_limit() -> Error {
    Error::plain("Hilbert Schur character shape and degree exceed the computation limit")
}

fn arithmetic_range() -> Error {
    Error::plain("Hilbert Schur character exceeds the exact arithmetic range")
}

fn euler_range() -> Error {
    Error::plain("Hilbert exterior Euler characteristic exceeds the signed 128-bit range")
}
