use std::cmp::Ordering;
use std::collections::BTreeMap;

use crate::binary_field::BinaryField;
use crate::comodule::Generator;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId};
use crate::q_koszul::{self, Complex};

const MAX_CELLS: usize = 300_000;
const MAX_PROJECTIVE_POINTS: usize = 4096;

type Entry = (usize, usize);
pub(crate) type Column = Vec<Entry>;

// An entry records the row and the parameter multiplying it. The underlying
// Q_i matrices are binary, and weighting them produces a field-valued column.
pub(crate) struct Endpoint {
    pub(crate) columns: Vec<Column>,
    groups: BTreeMap<usize, (usize, usize)>,
}

impl Endpoint {
    pub(crate) fn build(model: &mut Model, base: TypeId, selected: &[usize]) -> Result<Self> {
        let groups = q_koszul::margolis_groups(model, base)?;
        let mut offsets = BTreeMap::new();
        let mut size = 0;
        for &(degree, rank) in &groups {
            model.spend(1)?;
            offsets.insert(degree, (size, rank));
            size += rank;
            check_storage(size)?;
        }
        let mut complex = Complex::new(model, base, base, selected)?;
        let mut columns = Vec::with_capacity(size);
        let mut cells = size + 3 * offsets.len();
        check_storage(cells)?;
        for (degree, rank) in groups {
            for index in 0..rank {
                model.spend(1 + selected.len())?;
                let generator = Generator { degree, index };
                let mut column = Vec::new();
                for (variable, &primitive) in selected.iter().enumerate() {
                    let Some(lower) = degree.checked_sub(q_koszul::weight(primitive)) else {
                        continue;
                    };
                    let Some(&(offset, _)) = offsets.get(&lower) else {
                        continue;
                    };
                    let image = complex.margolis_image(model, base, primitive, lower, generator)?;
                    cells += 2 * image.len();
                    check_storage(cells)?;
                    model.spend(image.len())?;
                    column.extend(image.iter().map(|&row| (offset + row, variable)));
                }
                // Distinct primitive degrees give disjoint row groups.
                model.spend(column.len() * (1 + column.len().max(1).ilog2() as usize))?;
                column.sort_unstable();
                columns.push(column);
            }
        }
        Ok(Self {
            columns,
            groups: offsets,
        })
    }

    pub(crate) fn coordinate(&self, generator: Generator) -> Option<usize> {
        let &(offset, rank) = self.groups.get(&generator.degree)?;
        if generator.index < rank {
            Some(offset + generator.index)
        } else {
            None
        }
    }

    fn dimension(&self, model: &mut Model, field: BinaryField, point: &[usize]) -> Result<usize> {
        let mut pivots = BTreeMap::<usize, Column>::new();
        let mut cells = 0;
        for entries in &self.columns {
            model.spend(1 + entries.len())?;
            let mut column: Column = entries
                .iter()
                .filter_map(|&(row, variable)| {
                    let coefficient = point[variable];
                    (coefficient != 0).then_some((row, coefficient))
                })
                .collect();
            while let Some(&(leading, coefficient)) = column.last() {
                model.spend(1)?;
                if let Some(pivot) = pivots.get(&leading) {
                    column = add_multiple(model, field, &column, pivot, coefficient)?;
                } else {
                    cells += 1 + 2 * column.len();
                    check_storage(cells)?;
                    if coefficient != 1 {
                        let inverse = field.inverse(model, coefficient)?;
                        model.spend(field.degree() * column.len())?;
                        for (_, value) in &mut column {
                            *value = field.multiply(*value, inverse);
                        }
                    }
                    pivots.insert(leading, column);
                    break;
                }
            }
        }
        self.columns
            .len()
            .checked_sub(2 * pivots.len())
            .ok_or_else(|| {
                Error::plain("Margolis differential rank is incompatible with square zero")
            })
    }
}

pub(crate) fn add_multiple(
    model: &mut Model,
    field: BinaryField,
    mut left: &[Entry],
    mut right: &[Entry],
    factor: usize,
) -> Result<Column> {
    model.spend(left.len() + right.len() * if factor == 1 { 1 } else { field.degree() })?;
    let mut result = Vec::new();
    let scale = |value| {
        if factor == 1 {
            value
        } else {
            field.multiply(factor, value)
        }
    };
    while let (Some((&(a, x), rest_left)), Some((&(b, y), rest_right))) =
        (left.split_first(), right.split_first())
    {
        match a.cmp(&b) {
            Ordering::Less => {
                result.push((a, x));
                left = rest_left;
            }
            Ordering::Equal => {
                let coefficient = x ^ scale(y);
                if coefficient != 0 {
                    result.push((a, coefficient));
                }
                left = rest_left;
                right = rest_right;
            }
            Ordering::Greater => {
                result.push((b, scale(y)));
                right = rest_right;
            }
        }
    }
    result.extend_from_slice(left);
    result.extend(right.iter().map(|&(row, value)| (row, scale(value))));
    check_storage(2 * result.len())?;
    Ok(result)
}

pub(crate) struct Endpoints {
    source: Endpoint,
    target: Option<Endpoint>,
}

impl Endpoints {
    pub(crate) fn build(
        model: &mut Model,
        source: TypeId,
        target: TypeId,
        selected: &[usize],
    ) -> Result<Self> {
        let from = Endpoint::build(model, source, selected)?;
        let to = if target == source {
            None
        } else {
            Some(Endpoint::build(model, target, selected)?)
        };
        Ok(Self {
            source: from,
            target: to,
        })
    }

    pub(crate) fn dimension(
        &self,
        model: &mut Model,
        field: BinaryField,
        point: &[usize],
    ) -> Result<usize> {
        let source = self.source.dimension(model, field, point)?;
        let target = match &self.target {
            Some(target) => target.dimension(model, field, point)?,
            None => source,
        };
        source
            .checked_mul(target)
            .ok_or_else(|| Error::plain("Margolis Hom dimension is out of range"))
    }
}

pub(crate) fn dimension(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    selected: &[usize],
    modulus: usize,
    point: &[usize],
) -> Result<usize> {
    let field = BinaryField::new(model, modulus)?;
    model.check_q_endpoints(source, target, selected)?;
    field.check_point(point, selected.len(), "rvdim")?;
    let (active, coefficients): (Vec<_>, Vec<_>) = selected
        .iter()
        .copied()
        .zip(point.iter().copied())
        .filter(|&(_, coefficient)| coefficient != 0)
        .unzip();
    Endpoints::build(model, source, target, &active)?.dimension(model, field, &coefficients)
}

pub(crate) fn spectrum(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    selected: &[usize],
    modulus: usize,
) -> Result<BTreeMap<usize, usize>> {
    let field = BinaryField::new(model, modulus)?;
    model.check_q_endpoints(source, target, selected)?;
    if selected.is_empty() {
        return Ok(BTreeMap::new());
    }
    let space = ProjectiveSpace::new(field, selected.len(), Some(MAX_PROJECTIVE_POINTS))?;
    let endpoints = Endpoints::build(model, source, target, selected)?;
    let mut counts = BTreeMap::new();
    space.visit(model, |model, point| {
        let dimension = endpoints.dimension(model, field, point)?;
        *counts.entry(dimension).or_insert(0) += 1;
        Ok(())
    })?;
    Ok(counts)
}

pub(crate) struct ProjectiveSpace {
    order: usize,
    width: usize,
    size: usize,
}

impl ProjectiveSpace {
    pub(crate) fn new(field: BinaryField, width: usize, limit: Option<usize>) -> Result<Self> {
        let order = field.order();
        let mut size = 0usize;
        for _ in 0..width {
            size = size
                .checked_mul(order)
                .and_then(|n| n.checked_add(1))
                .filter(|&n| limit.is_none_or(|limit| n <= limit))
                .ok_or_else(|| match limit {
                    Some(limit) => Error::plain(format!(
                        "projective spectrum exceeds the limit of {limit} points"
                    )),
                    None => Error::plain("projective point count is out of range"),
                })?;
        }
        Ok(Self { order, width, size })
    }

    pub(crate) fn visit(
        &self,
        model: &mut Model,
        mut visit: impl FnMut(&mut Model, &[usize]) -> Result<()>,
    ) -> Result<()> {
        model.spend(self.size)?;
        let mut point = vec![0; self.width];
        let mut chart_size = 1;
        // Every projective point has a unique vector with first nonzero digit one.
        for first in (0..self.width).rev() {
            for mut digits in 0..chart_size {
                model.spend(point.len())?;
                point[first] = 1;
                for digit in &mut point[first + 1..] {
                    *digit = digits % self.order;
                    digits /= self.order;
                }
                visit(model, &point)?;
            }
            if first > 0 {
                chart_size *= self.order;
            }
        }
        Ok(())
    }
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain(
            "rank variety calculation exceeds the storage limit",
        ));
    }
    Ok(())
}
