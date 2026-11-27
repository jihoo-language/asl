use std::collections::BTreeMap;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;

pub(crate) type Vector = BTreeMap<usize, i128>;

pub(crate) fn arithmetic() -> Error {
    Error::plain("integral calculation exceeds the signed 128-bit coefficient range")
}

pub(crate) fn storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "integral calculation exceeds the sparse storage limit",
        ));
    }
    Ok(())
}

fn cells(entries: usize) -> usize {
    6 * entries
}

fn lookup_work(len: usize) -> usize {
    1 + len.max(1).ilog2() as usize
}

pub(crate) fn add_scaled(
    model: &mut Model,
    target: &mut Vector,
    source: &Vector,
    factor: i128,
) -> Result<()> {
    if factor == 0 {
        return Ok(());
    }
    model.spend(source.len() * lookup_work(target.len() + source.len()))?;
    for (&i, &value) in source {
        let old = target.get(&i).copied().unwrap_or(0);
        let change = factor.checked_mul(value).ok_or_else(arithmetic)?;
        let new = old.checked_add(change).ok_or_else(arithmetic)?;
        if new == 0 {
            target.remove(&i);
        } else {
            if old == 0 {
                storage(cells(target.len() + 1))?;
            }
            target.insert(i, new);
        }
    }
    Ok(())
}

#[derive(Clone, Copy)]
enum Axis {
    Row,
    Column,
}

impl Axis {
    fn other(self) -> Self {
        match self {
            Self::Row => Self::Column,
            Self::Column => Self::Row,
        }
    }
}

pub(crate) struct Matrix {
    rows: Vec<Vector>,
    columns: Vec<Vector>,
    entries: usize,
}

impl Matrix {
    pub(crate) fn new(model: &mut Model, rows: usize, columns: usize) -> Result<Self> {
        let size = rows.checked_add(columns).ok_or_else(arithmetic)?;
        storage(size.checked_mul(3).ok_or_else(arithmetic)?)?;
        model.spend(size)?;
        Ok(Self {
            rows: vec![Vector::new(); rows],
            columns: vec![Vector::new(); columns],
            entries: 0,
        })
    }

    pub(crate) fn cells(&self) -> usize {
        3 * (self.rows.len() + self.columns.len()) + 2 * cells(self.entries)
    }

    fn axis(&self, axis: Axis) -> &[Vector] {
        match axis {
            Axis::Row => &self.rows,
            Axis::Column => &self.columns,
        }
    }

    fn axis_mut(&mut self, axis: Axis) -> &mut [Vector] {
        match axis {
            Axis::Row => &mut self.rows,
            Axis::Column => &mut self.columns,
        }
    }

    fn value(&self, axis: Axis, index: usize, opposite: usize) -> i128 {
        self.axis(axis)[index].get(&opposite).copied().unwrap_or(0)
    }

    pub(crate) fn write(
        &mut self,
        model: &mut Model,
        row: usize,
        column: usize,
        value: i128,
    ) -> Result<()> {
        self.store(model, Axis::Row, row, column, value)
    }

    fn store(
        &mut self,
        model: &mut Model,
        axis: Axis,
        index: usize,
        opposite: usize,
        value: i128,
    ) -> Result<()> {
        model.spend(
            lookup_work(self.axis(axis)[index].len())
                + lookup_work(self.axis(axis.other())[opposite].len()),
        )?;
        if value == 0 {
            if self.axis_mut(axis)[index].remove(&opposite).is_some() {
                self.axis_mut(axis.other())[opposite].remove(&index);
                self.entries -= 1;
            }
        } else {
            if !self.axis(axis)[index].contains_key(&opposite) {
                storage(self.cells() + 2 * cells(1))?;
                self.entries += 1;
            }
            self.axis_mut(axis)[index].insert(opposite, value);
            self.axis_mut(axis.other())[opposite].insert(index, value);
        }
        Ok(())
    }

    fn swap(&mut self, model: &mut Model, axis: Axis, a: usize, b: usize) -> Result<()> {
        if a == b {
            return Ok(());
        }
        let side = self.axis(axis);
        model.spend(side[a].len() + side[b].len())?;
        let mut affected: Vec<_> = side[a].keys().chain(side[b].keys()).copied().collect();
        model.spend(affected.len() * lookup_work(affected.len()))?;
        affected.sort_unstable();
        affected.dedup();
        self.axis_mut(axis).swap(a, b);
        for j in affected {
            let (x, y) = (self.value(axis, a, j), self.value(axis, b, j));
            let opposite = &mut self.axis_mut(axis.other())[j];
            model.spend(4 * lookup_work(opposite.len()))?;
            set_entry(opposite, a, x);
            set_entry(opposite, b, y);
        }
        Ok(())
    }

    fn add(
        &mut self,
        model: &mut Model,
        axis: Axis,
        target: usize,
        source: usize,
        factor: i128,
    ) -> Result<()> {
        if factor == 0 {
            return Ok(());
        }
        let source = &self.axis(axis)[source];
        model.spend(source.len())?;
        let source: Vec<_> = source.iter().map(|(&i, &x)| (i, x)).collect();
        for (opposite, value) in source {
            let old = self.value(axis, target, opposite);
            let change = value.checked_mul(factor).ok_or_else(arithmetic)?;
            let new = old.checked_add(change).ok_or_else(arithmetic)?;
            self.store(model, axis, target, opposite, new)?;
        }
        Ok(())
    }

    fn negate(&mut self, model: &mut Model, axis: Axis, index: usize) -> Result<()> {
        let side = &self.axis(axis)[index];
        model.spend(side.len())?;
        let updates: Vec<_> = side.iter().map(|(&i, &x)| (i, x)).collect();
        for (opposite, value) in updates {
            self.store(
                model,
                axis,
                index,
                opposite,
                value.checked_neg().ok_or_else(arithmetic)?,
            )?;
        }
        Ok(())
    }

    fn smallest(&self, model: &mut Model, start: usize) -> Result<Option<(usize, usize)>> {
        let mut smallest = None;
        for (row, entries) in self.rows.iter().enumerate().skip(start) {
            model.spend(lookup_work(entries.len()))?;
            for (&column, &value) in entries.range(start..) {
                model.spend(1)?;
                let candidate = (value.unsigned_abs(), row, column);
                if smallest.is_none_or(|previous| candidate < previous) {
                    if candidate.0 == 1 {
                        return Ok(Some((row, column)));
                    }
                    smallest = Some(candidate);
                }
            }
        }
        Ok(smallest.map(|(_, row, column)| (row, column)))
    }

    fn offending(&self, model: &mut Model, start: usize, pivot: i128) -> Result<Option<usize>> {
        if pivot == 1 {
            return Ok(None);
        }
        for (row, entries) in self.rows.iter().enumerate().skip(start + 1) {
            model.spend(lookup_work(entries.len()))?;
            for (_, &value) in entries.range(start + 1..) {
                model.spend(1)?;
                if value % pivot != 0 {
                    return Ok(Some(row));
                }
            }
        }
        Ok(None)
    }

    fn neighbor(&self, model: &mut Model, fixed: Axis, index: usize) -> Result<Option<usize>> {
        let line = &self.axis(fixed)[index];
        model.spend(lookup_work(line.len()))?;
        Ok(line
            .range(index + 1..)
            .next()
            .map(|(&neighbor, _)| neighbor))
    }
}

fn set_entry(line: &mut Vector, index: usize, value: i128) {
    if value == 0 {
        line.remove(&index);
    } else {
        line.insert(index, value);
    }
}

struct Frame {
    rows: Vec<Vector>,
    cells: usize,
}

impl Frame {
    fn identity(model: &mut Model, width: usize) -> Result<Self> {
        let count = (3 + cells(1)) * width;
        storage(count)?;
        model.spend(count)?;
        Ok(Self {
            rows: (0..width).map(|i| [(i, 1)].into()).collect(),
            cells: count,
        })
    }

    fn add(&mut self, model: &mut Model, target: usize, source: usize, factor: i128) -> Result<()> {
        let (target, source) = if target < source {
            let (head, tail) = self.rows.split_at_mut(source);
            (&mut head[target], &tail[0])
        } else {
            let (head, tail) = self.rows.split_at_mut(target);
            (&mut tail[0], &head[source])
        };
        let before = target.len();
        add_scaled(model, target, source, factor)?;
        self.cells += cells(target.len());
        self.cells -= cells(before);
        Ok(())
    }

    fn negate(&mut self, model: &mut Model, index: usize) -> Result<()> {
        model.spend(self.rows[index].len())?;
        for value in self.rows[index].values_mut() {
            *value = value.checked_neg().ok_or_else(arithmetic)?;
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Track {
    None,
    Basis,
    Both,
}

#[derive(Default)]
pub(crate) struct Coordinates {
    basis: Option<Frame>,
    coordinate: Option<Frame>,
}

impl Coordinates {
    fn new(model: &mut Model, width: usize, track: Track) -> Result<Self> {
        let basis = (!matches!(track, Track::None))
            .then(|| Frame::identity(model, width))
            .transpose()?;
        let coordinate = matches!(track, Track::Both)
            .then(|| Frame::identity(model, width))
            .transpose()?;
        Ok(Self { basis, coordinate })
    }

    pub(crate) fn basis(&self, index: usize) -> &Vector {
        &self.basis.as_ref().expect("tracked basis").rows[index]
    }

    pub(crate) fn coordinate(&self, index: usize) -> &Vector {
        &self.coordinate.as_ref().expect("tracked coordinates").rows[index]
    }

    pub(crate) fn cells(&self) -> usize {
        self.basis.as_ref().map_or(0, |x| x.cells) + self.coordinate.as_ref().map_or(0, |x| x.cells)
    }

    fn swap(&mut self, a: usize, b: usize) {
        for frame in [&mut self.basis, &mut self.coordinate]
            .into_iter()
            .flatten()
        {
            frame.rows.swap(a, b);
        }
    }

    fn add(
        &mut self,
        model: &mut Model,
        axis: Axis,
        target: usize,
        source: usize,
        factor: i128,
    ) -> Result<()> {
        // Frames are stored in original coordinates. A row operation acts
        // directly on the coordinate forms, and inversely on the basis.
        let (direct, inverse) = match axis {
            Axis::Row => (&mut self.coordinate, &mut self.basis),
            Axis::Column => (&mut self.basis, &mut self.coordinate),
        };
        if let Some(frame) = direct {
            frame.add(model, target, source, factor)?;
        }
        if let Some(frame) = inverse {
            frame.add(
                model,
                source,
                target,
                factor.checked_neg().ok_or_else(arithmetic)?,
            )?;
        }
        Ok(())
    }

    fn negate(&mut self, model: &mut Model, index: usize) -> Result<()> {
        for frame in [&mut self.basis, &mut self.coordinate]
            .into_iter()
            .flatten()
        {
            frame.negate(model, index)?;
        }
        Ok(())
    }
}

pub(crate) struct Smith {
    pub(crate) factors: Vec<i128>,
    pub(crate) rows: Coordinates,
    pub(crate) columns: Coordinates,
}

impl Smith {
    pub(crate) fn empty() -> Self {
        Self {
            factors: Vec::new(),
            rows: Coordinates::default(),
            columns: Coordinates::default(),
        }
    }
}

struct Reduction {
    matrix: Matrix,
    rows: Coordinates,
    columns: Coordinates,
}

impl Reduction {
    fn coordinates_mut(&mut self, axis: Axis) -> &mut Coordinates {
        match axis {
            Axis::Row => &mut self.rows,
            Axis::Column => &mut self.columns,
        }
    }

    fn check(&self) -> Result<()> {
        storage(self.matrix.cells() + self.rows.cells() + self.columns.cells())
    }

    fn swap(&mut self, model: &mut Model, axis: Axis, i: usize, j: usize) -> Result<()> {
        self.matrix.swap(model, axis, i, j)?;
        self.coordinates_mut(axis).swap(i, j);
        Ok(())
    }

    fn add(
        &mut self,
        model: &mut Model,
        axis: Axis,
        target: usize,
        source: usize,
        factor: i128,
    ) -> Result<()> {
        self.matrix.add(model, axis, target, source, factor)?;
        self.coordinates_mut(axis)
            .add(model, axis, target, source, factor)?;
        self.check()
    }

    fn pivot(&mut self, model: &mut Model, index: usize) -> Result<Option<i128>> {
        let Some((row, column)) = self.matrix.smallest(model, index)? else {
            return Ok(None);
        };
        self.swap(model, Axis::Row, index, row)?;
        self.swap(model, Axis::Column, index, column)?;
        loop {
            model.spend(1)?;
            let mut pivot = self.matrix.value(Axis::Row, index, index);
            if pivot < 0 {
                self.matrix.negate(model, Axis::Row, index)?;
                self.rows.negate(model, index)?;
                pivot = pivot.checked_neg().ok_or_else(arithmetic)?;
            }
            if let Some(row) = self.matrix.neighbor(model, Axis::Column, index)? {
                let value = self.matrix.value(Axis::Row, row, index);
                let factor = value
                    .div_euclid(pivot)
                    .checked_neg()
                    .ok_or_else(arithmetic)?;
                self.add(model, Axis::Row, row, index, factor)?;
                if value.rem_euclid(pivot) != 0 {
                    self.swap(model, Axis::Row, row, index)?;
                }
                continue;
            }
            if let Some(column) = self.matrix.neighbor(model, Axis::Row, index)? {
                let value = self.matrix.value(Axis::Row, index, column);
                let factor = value
                    .div_euclid(pivot)
                    .checked_neg()
                    .ok_or_else(arithmetic)?;
                self.add(model, Axis::Column, column, index, factor)?;
                if value.rem_euclid(pivot) != 0 {
                    self.swap(model, Axis::Column, column, index)?;
                }
                continue;
            }
            if let Some(row) = self.matrix.offending(model, index, pivot)? {
                // The pivot row and column are clear. Add a row with an
                // indivisible entry to expose a smaller Euclidean remainder.
                self.add(model, Axis::Row, index, row, 1)?;
                continue;
            }
            return Ok(Some(pivot));
        }
    }
}

pub(crate) fn reduce(
    model: &mut Model,
    matrix: Matrix,
    rows: Track,
    columns: Track,
) -> Result<Smith> {
    let limit = matrix.rows.len().min(matrix.columns.len());
    let mut reduction = Reduction {
        rows: Coordinates::new(model, matrix.rows.len(), rows)?,
        columns: Coordinates::new(model, matrix.columns.len(), columns)?,
        matrix,
    };
    reduction.check()?;
    let mut factors = Vec::new();
    for i in 0..limit {
        let Some(value) = reduction.pivot(model, i)? else {
            break;
        };
        factors.push(value);
    }
    Ok(Smith {
        factors,
        rows: reduction.rows,
        columns: reduction.columns,
    })
}
