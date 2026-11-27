use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;

use super::{arithmetic, Dims, Distribution, Field, Hall, Shape};

const SUBSPACE_LIMIT: u128 = 50_000;
type Row = Vec<u32>;

pub(super) fn count(
    model: &mut Model,
    hall: &Hall,
    module: &Shape,
    dimensions: &Dims,
) -> Result<Distribution> {
    let geometry = Geometry::new(hall, module);
    let mut frames = Vec::new();
    let mut total = 1u128;
    let q = u128::from(hall.field.order());
    for (&width, &dimension) in geometry.widths.iter().zip(dimensions) {
        let patterns = frames_for(width, dimension, q)?;
        let choices = patterns.iter().map(|f| f.count as u128).sum::<u128>();
        total = total
            .checked_mul(choices)
            .filter(|&x| x <= SUBSPACE_LIMIT)
            .ok_or_else(subspace_bound)?;
        frames.push(patterns);
    }
    model.spend(total as usize)?;
    let mut walk = Walk {
        hall,
        geometry,
        frames,
        output: BTreeMap::new(),
    };
    walk.visit(model, &mut Vec::new())?;
    Ok(walk.output)
}

struct Frame {
    pivots: Vec<usize>,
    free: Vec<(usize, usize)>,
    count: usize,
}

fn frames_for(width: usize, dim: usize, q: u128) -> Result<Vec<Frame>> {
    if dim > width {
        return Err(Error::plain(
            "partition Hall subspace dimension exceeds its vertex space",
        ));
    }
    let mut frames = Vec::new();
    for mask in 0..1usize << width {
        if mask.count_ones() as usize != dim {
            continue;
        }
        let pivots: Vec<_> = (0..width).filter(|&j| mask & (1 << j) != 0).collect();
        let free: Vec<_> = pivots
            .iter()
            .enumerate()
            .flat_map(|(row, &pivot)| {
                (pivot + 1..width)
                    .filter(move |&col| mask & (1 << col) == 0)
                    .map(move |col| (row, col))
            })
            .collect();
        let count = q
            .checked_pow(free.len() as u32)
            .filter(|&x| x <= SUBSPACE_LIMIT)
            .ok_or_else(subspace_bound)? as usize;
        frames.push(Frame {
            pivots,
            free,
            count,
        });
    }
    Ok(frames)
}

fn subspace_bound() -> Error {
    Error::plain("partition Hall subspace choices exceed 50,000 for one dimension vector")
}

struct Subspace {
    pivots: Vec<usize>,
    rows: Vec<Row>,
}

struct Walk<'a> {
    hall: &'a Hall,
    geometry: Geometry,
    frames: Vec<Vec<Frame>>,
    output: Distribution,
}

impl Walk<'_> {
    fn visit(&mut self, model: &mut Model, chosen: &mut Vec<Subspace>) -> Result<()> {
        let vertex = chosen.len();
        if vertex == self.frames.len() {
            return self.accept(model, chosen);
        }
        let width = self.geometry.widths[vertex];
        let q = self.hall.field.order() as usize;
        for frame in 0..self.frames[vertex].len() {
            for index in 0..self.frames[vertex][frame].count {
                let frame = &self.frames[vertex][frame];
                model.spend(1 + width * frame.pivots.len())?;
                let mut rows = vec![vec![0; width]; frame.pivots.len()];
                for (i, &pivot) in frame.pivots.iter().enumerate() {
                    rows[i][pivot] = 1;
                }
                let mut remaining = index;
                for &(r, c) in &frame.free {
                    rows[r][c] = (remaining % q) as u32;
                    remaining /= q;
                }
                chosen.push(Subspace {
                    pivots: frame.pivots.clone(),
                    rows,
                });
                self.visit(model, chosen)?;
                chosen.pop();
            }
        }
        Ok(())
    }

    fn accept(&mut self, model: &mut Model, chosen: &[Subspace]) -> Result<()> {
        let field = &self.hall.field;
        for (source, vertex_paths) in self.geometry.paths.iter().enumerate() {
            if let Some(map) = vertex_paths.get(1).and_then(Option::as_ref) {
                let (target, width) = (map.target, self.geometry.widths[map.target]);
                let span = Span::from_subspace(width, &chosen[target]);
                for row in &chosen[source].rows {
                    let image = map.image(row, width);
                    if span.reduce(model, field, image)?.iter().any(|&v| v != 0) {
                        return Ok(());
                    }
                }
            }
        }
        let limit = self.hall.max_length + 1;
        let mut sub = vec![vec![0; limit + 1]; chosen.len()];
        let mut quo = sub.clone();
        for (source, paths) in self.geometry.paths.iter().enumerate() {
            for (k, map) in paths.iter().enumerate() {
                let Some(map) = map else {
                    continue;
                };
                let t = map.target;
                let width = self.geometry.widths[t];
                if k == 0 {
                    sub[source][0] = chosen[source].rows.len();
                    quo[source][0] = width - sub[source][0];
                    continue;
                }
                let mut restricted = Span::new(width);
                for row in &chosen[source].rows {
                    restricted.insert(model, field, map.image(row, width))?;
                }
                sub[source][k] = restricted.rank;
                let mut quotient = Span::from_subspace(width, &chosen[t]);
                for &axis in map.axis.iter().flatten() {
                    let mut unit = vec![0; width];
                    unit[axis] = 1;
                    quotient.insert(model, field, unit)?;
                }
                quo[source][k] = quotient.rank - chosen[t].rows.len();
            }
        }
        let n = classify(self.hall, &sub)?;
        let m = classify(self.hall, &quo)?;
        let coefficient = self.output.entry((m, n)).or_default();
        *coefficient = coefficient.checked_add(1).ok_or_else(arithmetic)?;
        Ok(())
    }
}

fn classify(hall: &Hall, ranks: &[Vec<usize>]) -> Result<Shape> {
    let failure =
        || Error::plain("partition Hall path ranks did not certify a uniserial decomposition");
    let mut parts = Vec::new();
    for (i, row) in ranks.iter().enumerate() {
        let count = |k: usize| {
            row[k]
                .checked_sub(hall.predecessor[i].map_or(0, |j| ranks[j][k + 1]))
                .ok_or_else(failure)
        };
        for k in 0..hall.max_length {
            let multiplicity = count(k)?.checked_sub(count(k + 1)?).ok_or_else(failure)?;
            if multiplicity != 0 {
                parts.push(hall.part(i as u128, (k + 1) as u128, multiplicity as u128)?);
            }
        }
    }
    let shape = Shape(parts);
    let mut check = vec![vec![0; hall.max_length + 2]; hall.vertices];
    for &(i, l, m) in &shape.0 {
        let h = &hall.path_heads[usize::from(i)][..usize::from(l)];
        for (start, &vertex) in h.iter().enumerate() {
            for rank in &mut check[vertex][..h.len() - start] {
                *rank += usize::from(m);
            }
        }
    }
    if check != ranks {
        return Err(failure());
    }
    Ok(shape)
}

struct Geometry {
    widths: Vec<usize>,
    paths: Vec<Vec<Option<PathMap>>>,
}

struct PathMap {
    target: usize,
    axis: Vec<Option<usize>>,
}

impl PathMap {
    fn image(&self, row: &Row, width: usize) -> Row {
        let mut output = vec![0; width];
        for (&x, j) in row.iter().zip(&self.axis) {
            if let Some(j) = j {
                output[*j] = x;
            }
        }
        output
    }
}

impl Geometry {
    fn new(hall: &Hall, module: &Shape) -> Self {
        let mut axes = vec![Vec::new(); hall.vertices];
        let mut copy = 0;
        for &(i, l, m) in &module.0 {
            for _ in 0..m {
                for (layer, &vertex) in hall.path_heads[usize::from(i)][..usize::from(l)]
                    .iter()
                    .enumerate()
                {
                    axes[vertex].push((copy, layer));
                }
                copy += 1;
            }
        }
        let paths = axes
            .iter()
            .enumerate()
            .map(|(source, axis)| {
                let mut target = Some(source);
                (0..hall.max_length)
                    .map(|k| {
                        let end = target?;
                        target = hall.successor[end];
                        let images = axis
                            .iter()
                            .map(|&(i, l)| {
                                axes[end].iter().position(|&(j, m)| i == j && l + k == m)
                            })
                            .collect();
                        Some(PathMap {
                            target: end,
                            axis: images,
                        })
                    })
                    .collect()
            })
            .collect();
        Self {
            widths: axes.iter().map(Vec::len).collect(),
            paths,
        }
    }
}

struct Span {
    rows: Vec<Option<Row>>,
    rank: usize,
}

impl Span {
    fn new(width: usize) -> Self {
        Self {
            rows: vec![None; width],
            rank: 0,
        }
    }

    fn from_subspace(width: usize, sub: &Subspace) -> Self {
        let mut span = Self::new(width);
        for (&pivot, row) in sub.pivots.iter().zip(&sub.rows) {
            span.rows[pivot] = Some(row.clone());
        }
        span.rank = sub.rows.len();
        span
    }

    fn reduce(&self, model: &mut Model, field: &Field, mut row: Row) -> Result<Row> {
        for (pivot, old) in self.rows.iter().enumerate() {
            if let Some(old) = old {
                let scalar = row[pivot];
                if scalar == 0 {
                    continue;
                }
                model.spend(row.len() - pivot)?;
                for (x, &y) in row[pivot..].iter_mut().zip(&old[pivot..]) {
                    *x = field.sub(*x, field.mul(scalar, y));
                }
            }
        }
        Ok(row)
    }

    fn insert(&mut self, model: &mut Model, field: &Field, row: Row) -> Result<()> {
        let mut row = self.reduce(model, field, row)?;
        if let Some(pivot) = row.iter().position(|&x| x != 0) {
            let scale = field.inverse(model, row[pivot])?;
            model.spend(row.len() - pivot)?;
            for term in &mut row[pivot..] {
                *term = field.mul(*term, scale);
            }
            self.rows[pivot] = Some(row);
            self.rank += 1;
        }
        Ok(())
    }
}
