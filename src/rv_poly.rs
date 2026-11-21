use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use crate::binary_field::BinaryField;
use crate::comodule::Generator;
use crate::complex::xor_sorted;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId};
use crate::rank_variety::Endpoint;

const MAX_CELLS: usize = 300_000;

pub(crate) type Monomial = Rc<[usize]>;
pub(crate) type Polynomial = Rc<[Monomial]>;

fn degree_error() -> Error {
    Error::plain("rank polynomial degree is out of range")
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain(
            "symbolic rank calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

pub(crate) fn cells(polynomial: &[Monomial]) -> usize {
    polynomial
        .first()
        .map_or(1, |term| 1 + (1 + term.len()) * polynomial.len())
}

pub(crate) fn persist(model: &mut Model, polynomial: Polynomial) -> Result<Polynomial> {
    model.allocate(cells(&polynomial))?;
    Ok(polynomial)
}

pub(crate) fn degree(powers: &[usize]) -> Result<usize> {
    powers.iter().try_fold(0usize, |sum, &power| {
        sum.checked_add(power).ok_or_else(degree_error)
    })
}

pub(crate) fn monomial(powers: Vec<usize>) -> Polynomial {
    [powers.into()].into()
}

pub(crate) fn add(model: &mut Model, left: &[Monomial], right: &[Monomial]) -> Result<Polynomial> {
    model.spend(cells(left) + cells(right))?;
    let result = xor_sorted(left, right);
    check_storage(cells(&result))?;
    Ok(result.into())
}

pub(crate) fn compare(
    model: &mut Model,
    left: &[Monomial],
    right: &[Monomial],
) -> Result<Ordering> {
    for (a, b) in left.iter().rev().zip(right.iter().rev()) {
        model.spend(1 + a.len())?;
        let order = a.cmp(b);
        if order != Ordering::Equal {
            return Ok(order);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

pub(crate) fn basis(
    model: &mut Model,
    width: usize,
    degree: usize,
    limit: usize,
) -> Result<Vec<Monomial>> {
    fn enumerate(
        model: &mut Model,
        powers: &mut [usize],
        position: usize,
        degree: usize,
        output: &mut Vec<Monomial>,
        limit: usize,
    ) -> Result<()> {
        model.spend(1)?;
        let free = powers.len() - position;
        if free <= 1 {
            if free == 1 {
                powers[position] = degree;
            } else if degree != 0 {
                return Ok(());
            }
            if output.len() == limit {
                return Err(Error::plain(
                    "rank polynomial type exceeds the enumeration limit",
                ));
            }
            model.spend(1 + powers.len())?;
            output.push(Rc::from(&*powers));
            return Ok(());
        }
        for power in 0..=degree {
            powers[position] = power;
            enumerate(model, powers, position + 1, degree - power, output, limit)?;
        }
        Ok(())
    }

    let mut basis = Vec::new();
    enumerate(model, &mut vec![0; width], 0, degree, &mut basis, limit)?;
    model.allocate(cells(&basis))?;
    Ok(basis)
}

fn product_term(left: &[usize], right: &[usize]) -> Result<Monomial> {
    left.iter()
        .zip(right)
        .map(|(&a, &b)| a.checked_add(b).ok_or_else(degree_error))
        .collect()
}

fn toggle(terms: &mut BTreeSet<Monomial>, term: Monomial) {
    if !terms.remove(&term) {
        terms.insert(term);
    }
}

pub(crate) fn multiply(
    model: &mut Model,
    left: &[Monomial],
    right: &[Monomial],
) -> Result<Polynomial> {
    let Some(width) = left.first().map(|word| word.len()) else {
        return Ok([].into());
    };
    let mut terms = BTreeSet::new();
    let cost = 1 + width;
    for a in left {
        for b in right {
            model.spend(cost * (1 + terms.len().max(1).ilog2() as usize))?;
            toggle(&mut terms, product_term(a, b)?);
            check_storage(cost * terms.len())?;
        }
    }
    Ok(terms.into_iter().collect())
}

pub(crate) fn hasse(
    model: &mut Model,
    polynomial: &[Monomial],
    orders: &[usize],
) -> Result<Polynomial> {
    model.spend(cells(polynomial))?;
    // Lucas parity: each binary digit of a derivative order must occur in
    // the matching exponent. Subtracting a fixed order preserves sorting.
    Ok(polynomial
        .iter()
        .filter(|term| {
            term.iter()
                .zip(orders)
                .all(|(&power, &order)| power & order == order)
        })
        .map(|term| {
            term.iter()
                .zip(orders)
                .map(|(&power, &order)| power - order)
                .collect::<Monomial>()
        })
        .collect())
}

pub(crate) fn vanishes(
    model: &mut Model,
    polynomial: &[Monomial],
    field: BinaryField,
    point: &[usize],
) -> Result<bool> {
    let mut sum = 0;
    for term in polynomial {
        let mut product = 1;
        model.spend(1 + field.degree() * term.len())?;
        for (&coordinate, &exponent) in point.iter().zip(term.iter()) {
            product = field.multiply(product, power(model, field, coordinate, exponent)?);
            if product == 0 {
                break;
            }
        }
        sum ^= product;
    }
    Ok(sum == 0)
}

fn power(
    model: &mut Model,
    field: BinaryField,
    mut value: usize,
    mut exponent: usize,
) -> Result<usize> {
    if exponent == 0 || value == 1 {
        return Ok(1);
    }
    if value == 0 {
        return Ok(0);
    }
    exponent %= field.order() - 1;
    let mut output = 1;
    while exponent != 0 {
        model.spend(2 * field.degree())?;
        if exponent & 1 != 0 {
            output = field.multiply(output, value);
        }
        value = field.multiply(value, value);
        exponent >>= 1;
    }
    Ok(output)
}

fn exact_divide(
    model: &mut Model,
    numerator: &[Monomial],
    denominator: &[Monomial],
) -> Result<Polynomial> {
    let leading = denominator.last().expect("nonzero preceding pivot");
    if denominator.len() == 1 {
        model.spend(cells(numerator))?;
        return numerator
            .iter()
            .map(|term| quotient_term(term, leading))
            .collect();
    }
    model.spend(cells(numerator) * (1 + numerator.len().max(1).ilog2() as usize))?;
    let mut remainder: BTreeSet<Monomial> = numerator.iter().cloned().collect();
    let mut quotient = Vec::new();
    let cost = 1 + leading.len();
    while let Some(term) = remainder.pop_last() {
        let factor = quotient_term(&term, leading)?;
        for lower in &denominator[..denominator.len() - 1] {
            model.spend(cost * (1 + remainder.len().max(1).ilog2() as usize))?;
            toggle(&mut remainder, product_term(&factor, lower)?);
            check_storage(cost * (remainder.len() + quotient.len()))?;
        }
        quotient.push(factor);
    }
    quotient.reverse();
    Ok(quotient.into())
}

fn quotient_term(term: &[usize], divisor: &[usize]) -> Result<Monomial> {
    term.iter()
        .zip(divisor)
        .map(|(&a, &b)| {
            a.checked_sub(b).ok_or_else(|| {
                Error::plain("symbolic determinant encountered a nonexact pivot division")
            })
        })
        .collect()
}

pub(crate) fn minor(
    model: &mut Model,
    base: TypeId,
    selected: &[usize],
    targets: &[Generator],
    sources: &[Generator],
) -> Result<Polynomial> {
    model.check_q_endpoints(base, base, selected)?;
    if targets.len() != sources.len() {
        return Err(Error::plain(
            "rvMinor requires equally many target rows and source columns",
        ));
    }
    if targets.is_empty() {
        return Ok(monomial(vec![0; selected.len()]));
    }
    let endpoint = Endpoint::build(model, base, selected)?;
    let rows = minor_coordinates(model, &endpoint, targets, "row")?;
    let columns = minor_coordinates(model, &endpoint, sources, "column")?;
    determinant(model, &endpoint, selected.len(), &rows, &columns)
}

fn minor_coordinates(
    model: &mut Model,
    endpoint: &Endpoint,
    coordinates: &[Generator],
    axis: &str,
) -> Result<Vec<usize>> {
    let mut indices = BTreeSet::new();
    for &coordinate in coordinates {
        model.spend(1 + endpoint.columns.len().max(1).ilog2() as usize)?;
        let index = endpoint.coordinate(coordinate).ok_or_else(|| {
            Error::plain(format!(
                "rvMinor {axis} coordinate ({}, {}) is out of range",
                coordinate.degree, coordinate.index,
            ))
        })?;
        // All permutation signs are +1 in F2. Sorting also gives stable pivots.
        if !indices.insert(index) {
            return Err(Error::plain(format!("rvMinor repeats a {axis} coordinate")));
        }
    }
    Ok(indices.into_iter().collect())
}

fn determinant(
    model: &mut Model,
    endpoint: &Endpoint,
    width: usize,
    rows: &[usize],
    columns: &[usize],
) -> Result<Polynomial> {
    let size = rows.len();
    check_storage(2 * size)?;
    let mut row_edges = vec![BTreeMap::new(); size];
    let mut column_edges = vec![BTreeSet::new(); size];
    let mut entries = 0;
    for (index, &column) in columns.iter().enumerate() {
        for &(row, variable) in &endpoint.columns[column] {
            model.spend(1 + size.max(1).ilog2() as usize)?;
            if let Ok(row) = rows.binary_search(&row) {
                row_edges[row].insert(index, variable);
                column_edges[index].insert(row);
                entries += 3;
                check_storage(2 * size + entries)?;
            }
        }
    }
    let mut forced_rows = Vec::new();
    let mut forced_columns = Vec::new();
    for (i, (row, column)) in row_edges.iter().zip(&column_edges).enumerate() {
        if row.is_empty() || column.is_empty() {
            return Ok([].into());
        }
        if row.len() == 1 {
            forced_rows.push(i);
        }
        if column.len() == 1 {
            forced_columns.push(i);
        }
    }
    model.spend(1 + size)?;
    let mut kept_rows = vec![true; size];
    let mut kept_columns = vec![true; size];
    let mut forced = vec![0; width];

    // A row or column with one entry forces it in every determinant term.
    // Propagate those choices before fraction-free elimination.
    while !forced_rows.is_empty() || !forced_columns.is_empty() {
        let (row, column) = if let Some(row) = forced_rows.pop() {
            if !kept_rows[row] {
                continue;
            }
            (row, *row_edges[row].first_key_value().unwrap().0)
        } else {
            let column = forced_columns.pop().unwrap();
            if !kept_columns[column] {
                continue;
            }
            (*column_edges[column].first().unwrap(), column)
        };
        kept_rows[row] = false;
        kept_columns[column] = false;
        forced[row_edges[row][&column]] += 1;
        for other in std::mem::take(&mut column_edges[column]) {
            model.spend(1 + row_edges[other].len().max(1).ilog2() as usize)?;
            if other == row {
                continue;
            }
            row_edges[other].remove(&column);
            match row_edges[other].len() {
                0 => return Ok([].into()),
                1 => forced_rows.push(other),
                _ => {}
            }
        }
        for (other, _) in std::mem::take(&mut row_edges[row]) {
            model.spend(1 + column_edges[other].len().max(1).ilog2() as usize)?;
            if other == column {
                continue;
            }
            column_edges[other].remove(&row);
            match column_edges[other].len() {
                0 => return Ok([].into()),
                1 => forced_columns.push(other),
                _ => {}
            }
        }
    }

    let remaining = kept_rows.iter().filter(|&&keep| keep).count();
    if remaining == 0 {
        return Ok(monomial(forced));
    }
    let active_columns: Vec<_> = kept_columns
        .into_iter()
        .enumerate()
        .filter_map(|(column, keep)| keep.then_some(column))
        .collect();
    let mut storage = remaining.saturating_mul(remaining);
    check_storage(storage)?;
    model.spend(storage + width * width)?;
    let variables: Vec<_> = (0..width)
        .map(|variable| {
            let mut powers = vec![0; width];
            powers[variable] = 1;
            monomial(powers)
        })
        .collect();
    let mut matrix = Vec::with_capacity(remaining);
    for row in row_edges
        .iter()
        .zip(kept_rows)
        .filter_map(|(row, keep)| keep.then_some(row))
    {
        let mut entries = Vec::with_capacity(remaining);
        model.spend(remaining * (1 + row.len().max(1).ilog2() as usize))?;
        for column in &active_columns {
            let polynomial = row
                .get(column)
                .map_or_else(|| [].into(), |&i| variables[i].clone());
            storage += cells(&polynomial);
            check_storage(storage)?;
            entries.push(polynomial);
        }
        matrix.push(entries);
    }
    let determinant = bareiss(model, matrix, storage)?;
    multiply(model, &monomial(forced), &determinant)
}

fn bareiss(
    model: &mut Model,
    mut matrix: Vec<Vec<Polynomial>>,
    mut storage: usize,
) -> Result<Polynomial> {
    let size = matrix.len();
    let mut previous: Option<Polynomial> = None;
    for i in 0..size - 1 {
        model.spend(size - i)?;
        let Some(row) = (i..size).find(|&row| !matrix[row][i].is_empty()) else {
            return Ok([].into());
        };
        matrix.swap(i, row);
        let pivot = matrix[i][i].clone();
        for row in i + 1..size {
            let factor = matrix[row][i].clone();
            for column in i + 1..size {
                model.spend(1)?;
                let current = &matrix[row][column];
                let diagonal = multiply(model, &pivot, current)?;
                let cross = multiply(model, &factor, &matrix[i][column])?;
                let numerator = add(model, &diagonal, &cross)?;
                let value = match &previous {
                    Some(previous) => exact_divide(model, &numerator, previous)?,
                    None => numerator,
                };
                storage = storage - cells(current) + cells(&value);
                check_storage(storage)?;
                matrix[row][column] = value;
            }
            storage = storage - cells(&matrix[row][i]) + 1;
            matrix[row][i] = [].into();
        }
        previous = Some(pivot);
    }
    Ok(matrix[size - 1][size - 1].clone())
}

#[derive(Default)]
pub(crate) struct Minors {
    pivots: BTreeMap<Monomial, Polynomial>,
    storage: usize,
}

impl Minors {
    pub(crate) fn build(
        model: &mut Model,
        base: TypeId,
        selected: &[usize],
        degree: usize,
    ) -> Result<Self> {
        model.check_q_endpoints(base, base, selected)?;
        let mut result = Self::default();
        if degree == 0 {
            result.insert(model, monomial(vec![0; selected.len()]))?;
            return Ok(result);
        }
        let endpoint = Endpoint::build(model, base, selected)?;
        // The parameterized operator squares to zero over the fraction field.
        if degree > endpoint.columns.len() / 2 {
            return Ok(result);
        }
        if degree == 1 {
            let mut occurring = vec![false; selected.len()];
            for entries in &endpoint.columns {
                model.spend(1 + entries.len())?;
                for &(_, variable) in entries {
                    occurring[variable] = true;
                }
            }
            for (variable, present) in occurring.into_iter().enumerate() {
                if present {
                    let mut powers = vec![0; selected.len()];
                    powers[variable] = 1;
                    result.insert(model, monomial(powers))?;
                }
            }
            return Ok(result);
        }
        let mut rows = BTreeSet::new();
        let mut columns = Vec::new();
        for (column, entries) in endpoint.columns.iter().enumerate() {
            model
                .spend(1 + entries.len() * (1 + endpoint.columns.len().max(1).ilog2() as usize))?;
            if !entries.is_empty() {
                columns.push(column);
                rows.extend(entries.iter().map(|&(row, _)| row));
            }
        }
        if degree > rows.len().min(columns.len()) {
            return Ok(result);
        }
        let rows: Vec<_> = rows.into_iter().collect();
        let mut row_indices = (0..degree).collect::<Vec<_>>();
        let mut minor_rows = Vec::with_capacity(degree);
        let mut minor_columns = Vec::with_capacity(degree);
        let mut column_indices = Vec::with_capacity(degree);
        let dimension = small_dimension(selected.len(), degree);
        'rows: loop {
            minor_rows.clear();
            minor_rows.extend(row_indices.iter().map(|&index| rows[index]));
            column_indices.clear();
            column_indices.extend(0..degree);
            loop {
                model.spend(1 + 2 * degree)?;
                minor_columns.clear();
                minor_columns.extend(column_indices.iter().map(|&index| columns[index]));
                let value = determinant(
                    model,
                    &endpoint,
                    selected.len(),
                    &minor_rows,
                    &minor_columns,
                )?;
                result.insert(model, value)?;
                if Some(result.pivots.len()) == dimension {
                    break 'rows;
                }
                if !next_subset(&mut column_indices, columns.len()) {
                    break;
                }
            }
            if !next_subset(&mut row_indices, rows.len()) {
                break;
            }
        }
        result.reduce_rows(model)?;
        Ok(result)
    }

    fn insert(&mut self, model: &mut Model, mut polynomial: Polynomial) -> Result<()> {
        while let Some(leading) = polynomial.last() {
            model.spend(1 + leading.len() * (1 + self.pivots.len().max(1).ilog2() as usize))?;
            if let Some(row) = self.pivots.get(leading) {
                polynomial = add(model, &polynomial, row)?;
            } else {
                self.storage += cells(&polynomial) + leading.len() + 1;
                check_storage(self.storage)?;
                self.pivots.insert(leading.clone(), polynomial);
                break;
            }
        }
        Ok(())
    }

    pub(crate) fn contains(&self, model: &mut Model, mut polynomial: Polynomial) -> Result<bool> {
        while let Some(leading) = polynomial.last() {
            model.spend(1 + leading.len() * (1 + self.pivots.len().max(1).ilog2() as usize))?;
            let Some(row) = self.pivots.get(leading) else {
                return Ok(false);
            };
            polynomial = add(model, &polynomial, row)?;
        }
        Ok(true)
    }

    pub(crate) fn storage(&self) -> usize {
        self.storage
    }

    pub(crate) fn polynomials(&self) -> impl Iterator<Item = &Polynomial> {
        self.pivots.values()
    }

    fn reduce_rows(&mut self, model: &mut Model) -> Result<()> {
        // Remove every pivot except the leading one to obtain the unique
        // reduced row basis in lexicographic monomial order.
        let mut result = BTreeMap::<Monomial, Polynomial>::new();
        for (leading, row) in std::mem::take(&mut self.pivots) {
            let mut remaining: BTreeSet<_> = row.iter().cloned().collect();
            let mut reduced = Vec::new();
            model.spend(cells(&row) * (1 + row.len().max(1).ilog2() as usize))?;
            while let Some(term) = remaining.pop_last() {
                model.spend(1 + term.len() * (1 + result.len().max(1).ilog2() as usize))?;
                if let Some(pivot) = result.get(&term) {
                    for term in &pivot[..pivot.len() - 1] {
                        model.spend(
                            (1 + term.len()) * (1 + remaining.len().max(1).ilog2() as usize),
                        )?;
                        toggle(&mut remaining, term.clone());
                        check_storage((1 + term.len()) * (remaining.len() + reduced.len()))?;
                    }
                } else {
                    reduced.push(term);
                }
            }
            reduced.reverse();
            self.storage = self.storage - cells(&row) + cells(&reduced);
            check_storage(self.storage)?;
            result.insert(leading, reduced.into());
        }
        self.pivots = result;
        Ok(())
    }
}

fn small_dimension(width: usize, degree: usize) -> Option<usize> {
    let gaps = width.checked_sub(1)?;
    let positions = degree + gaps;
    let mut count = 1;
    for i in 1..=degree.min(gaps) {
        count = count * (positions + 1 - i) as u128 / i as u128;
        if count > MAX_CELLS as u128 {
            return None;
        }
    }
    Some(count as usize)
}

fn next_subset(indices: &mut [usize], bound: usize) -> bool {
    let size = indices.len();
    for position in (0..size).rev() {
        if indices[position] < bound - (size - position) {
            let first = indices[position] + 1;
            for (offset, index) in indices[position..].iter_mut().enumerate() {
                *index = first + offset;
            }
            return true;
        }
    }
    false
}
