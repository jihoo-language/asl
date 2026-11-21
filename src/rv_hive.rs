use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};

const MAX_RANK: usize = 24;
const MAX_HORN_CELLS: usize = 300_000;

pub(crate) type Interior = Rc<[i128]>;

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Boundary(pub(crate) [Rc<[i128]>; 3]);

impl Boundary {
    pub(crate) fn checked(model: &mut Model, weights: [Vec<i128>; 3]) -> Result<Self> {
        let rank = weights[0].len();
        if weights.iter().any(|weight| weight.len() != rank) {
            return Err(Error::plain(
                "hive boundary weights must have the same rank",
            ));
        }
        if rank > MAX_RANK {
            return Err(Error::plain("hive rank exceeds the limit of 24"));
        }
        model.spend(3 * rank)?;
        if weights
            .iter()
            .any(|weight| weight.windows(2).any(|pair| pair[0] > pair[1]))
        {
            return Err(Error::plain(
                "hive boundary weights must be weakly increasing",
            ));
        }
        Ok(Self(weights.map(Rc::from)))
    }

    pub(crate) fn rank(&self) -> usize {
        self.0[0].len()
    }
}

fn arithmetic() -> Error {
    Error::plain("hive integer arithmetic is out of range")
}

fn add(left: i128, right: i128) -> Result<i128> {
    left.checked_add(right).ok_or_else(arithmetic)
}

fn sub(left: i128, right: i128) -> Result<i128> {
    left.checked_sub(right).ok_or_else(arithmetic)
}

fn sum(terms: impl IntoIterator<Item = i128>) -> Result<i128> {
    terms.into_iter().try_fold(0, add)
}

struct Triangle {
    rows: Vec<Vec<i128>>,
    expected: i128,
    actual: i128,
}

impl Triangle {
    fn new(model: &mut Model, boundary: &Boundary) -> Result<Self> {
        let [lambda, mu, nu] = &boundary.0;
        let rank = boundary.rank();
        model.spend((rank + 1) * (rank + 2) / 2)?;
        let mut rows: Vec<_> = (0..=rank).map(|j| vec![0; rank + 1 - j]).collect();
        for (i, &weight) in lambda.iter().rev().enumerate() {
            rows[0][i + 1] = add(rows[0][i], weight)?;
        }
        let mut right = rows[0][rank];
        let mut left = 0;
        for (j, (&m, &n)) in mu.iter().rev().zip(nu.iter().rev()).enumerate() {
            left = add(left, n)?;
            right = add(right, m)?;
            rows[j + 1][0] = left;
            rows[j + 1][rank - j - 1] = right;
        }
        Ok(Self {
            rows,
            expected: right,
            actual: left,
        })
    }

    fn balanced(&self) -> bool {
        self.expected == self.actual
    }

    fn interiors(&self) -> impl Iterator<Item = i128> + '_ {
        self.rows[1..]
            .iter()
            .flat_map(|row| row.iter().skip(1).take(row.len().saturating_sub(2)))
            .copied()
    }

    fn set_interior(&mut self, values: &[i128]) -> Result<()> {
        let rank = self.rows.len() - 1;
        let size = rank.saturating_sub(1) * rank.saturating_sub(2) / 2;
        if values.len() != size {
            return Err(Error::plain(format!(
                "rank-{rank} hive needs {size} interior coordinates, found {}",
                values.len()
            )));
        }
        let mut values = values.iter().copied();
        for row in &mut self.rows[1..] {
            let size = row.len().saturating_sub(2);
            for vertex in row.iter_mut().skip(1).take(size) {
                *vertex = values.next().unwrap();
            }
        }
        Ok(())
    }

    fn vertical_right(&self, i: usize, j: usize) -> Result<i128> {
        sub(self.rows[j - 1][i + 1], self.rows[j - 2][i + 1])
    }

    fn vertical_ok(&self, i: usize, j: usize) -> Result<bool> {
        Ok(j < 2 || sub(self.rows[j][i], self.rows[j - 1][i])? <= self.vertical_right(i, j)?)
    }

    fn horizontal_range(&self, i: usize, j: usize) -> Result<(i128, i128)> {
        let previous = &self.rows[j - 1];
        Ok((
            sub(previous[i + 1], previous[i])?,
            sub(previous[i], previous[i - 1])?,
        ))
    }

    fn check(&self, model: &mut Model) -> Result<()> {
        if !self.balanced() {
            return Err(Error::plain(
                "hive boundary totals must satisfy |lambda| + |mu| = |nu|",
            ));
        }
        for j in 1..self.rows.len() {
            model.spend(self.rows[j].len())?;
            for i in 0..self.rows[j].len() {
                if i > 0 {
                    let slope = sub(self.rows[j][i], self.rows[j][i - 1])?;
                    let (lower, upper) = self.horizontal_range(i, j)?;
                    if slope < lower || slope > upper {
                        return Err(Error::plain(format!(
                            "hive horizontal rhombus inequality fails at ({i}, {j})"
                        )));
                    }
                }
                if !self.vertical_ok(i, j)? {
                    return Err(Error::plain(format!(
                        "hive vertical rhombus inequality fails at ({i}, {j})"
                    )));
                }
            }
        }
        Ok(())
    }

    fn visit(
        &mut self,
        model: &mut Model,
        mut solution: impl FnMut(&mut Model, &Self) -> Result<bool>,
    ) -> Result<()> {
        if self.balanced() {
            self.visit_row(model, 1, &mut solution)?;
        }
        Ok(())
    }

    fn visit_row(
        &mut self,
        model: &mut Model,
        j: usize,
        solution: &mut impl FnMut(&mut Model, &Self) -> Result<bool>,
    ) -> Result<bool> {
        model.spend(1)?;
        if j == self.rows.len() {
            return solution(model, self);
        }
        let width = self.rows[j].len() - 1;
        if !self.vertical_ok(0, j)? || !self.vertical_ok(width, j)? {
            return Ok(false);
        }
        let previous = &self.rows[j - 1];
        let total = sub(self.rows[j][width], self.rows[j][0])?;
        if total < sub(previous[width + 1], previous[1])?
            || total > sub(previous[width], previous[0])?
        {
            return Ok(false);
        }
        if width <= 1 {
            return self.visit_row(model, j + 1, solution);
        }
        self.visit_vertex(model, 1, j, solution)
    }

    fn visit_vertex(
        &mut self,
        model: &mut Model,
        i: usize,
        j: usize,
        solution: &mut impl FnMut(&mut Model, &Self) -> Result<bool>,
    ) -> Result<bool> {
        model.spend(1)?;
        let row = &self.rows[j];
        let previous = &self.rows[j - 1];
        let width = row.len() - 1;
        let (low_step, high_step) = self.horizontal_range(i, j)?;
        // The as-yet-unassigned slopes must also reach the fixed right boundary.
        // Their interlacing intervals telescope along the preceding row.
        let mut lower =
            add(row[i - 1], low_step)?.max(sub(row[width], sub(previous[width], previous[i])?)?);
        let mut upper = add(row[i - 1], high_step)?
            .min(sub(row[width], sub(previous[width + 1], previous[i + 1])?)?);
        if j >= 2 {
            upper = upper.min(add(previous[i], self.vertical_right(i, j)?)?);
        }
        while lower <= upper {
            model.spend(1)?;
            self.rows[j][i] = lower;
            let stop = if i + 1 == width {
                self.visit_row(model, j + 1, solution)?
            } else {
                self.visit_vertex(model, i + 1, j, solution)?
            };
            if stop {
                return Ok(true);
            }
            if lower == upper {
                break;
            }
            lower += 1;
        }
        Ok(false)
    }

    fn witness(&mut self, model: &mut Model) -> Result<Option<Interior>> {
        let mut witness = None;
        self.visit(model, |_, triangle| {
            witness = Some(triangle.interiors().collect());
            Ok(true)
        })?;
        Ok(witness)
    }
}

pub(crate) fn witness(model: &mut Model, boundary: &Boundary) -> Result<Option<Interior>> {
    Triangle::new(model, boundary)?.witness(model)
}

pub(crate) fn count(model: &mut Model, boundary: &Boundary) -> Result<u128> {
    let mut count = 0u128;
    Triangle::new(model, boundary)?.visit(model, |_, _| {
        count = count.checked_add(1).ok_or_else(arithmetic)?;
        Ok(false)
    })?;
    Ok(count)
}

pub(crate) fn elements(model: &mut Model, boundary: &Boundary, limit: usize) -> Result<Vec<Value>> {
    let mut values = Vec::new();
    Triangle::new(model, boundary)?.visit(model, |model, triangle| {
        if values.len() == limit {
            return Err(Error::plain(format!(
                "type exceeds the enumeration limit of {limit} elements"
            )));
        }
        let interior: Interior = triangle.interiors().collect();
        model.allocate(1 + interior.len())?;
        values.push(Value::RVHive(interior));
        Ok(false)
    })?;
    Ok(values)
}

pub(crate) fn literal(model: &mut Model, boundary: &Boundary, values: Vec<i128>) -> Result<Value> {
    let mut triangle = Triangle::new(model, boundary)?;
    triangle.set_interior(&values)?;
    triangle.check(model)?;
    model.allocate(values.len())?;
    Ok(Value::RVHive(values.into()))
}

pub(crate) fn vertex(
    model: &mut Model,
    boundary: &Boundary,
    interior: &[i128],
    i: usize,
    j: usize,
) -> Result<i128> {
    let mut triangle = Triangle::new(model, boundary)?;
    triangle.set_interior(interior)?;
    triangle
        .rows
        .get(j)
        .and_then(|row| row.get(i))
        .copied()
        .ok_or_else(|| Error::plain("hive vertex lies outside the triangle"))
}

pub(crate) enum Horn {
    Admissible,
    Total {
        expected: i128,
        actual: i128,
    },
    Violation {
        indices: [Vec<usize>; 3],
        left: i128,
        right: i128,
    },
}

struct Subset {
    indices: Vec<usize>,
    weight: [i128; 3],
    shape: Rc<[i128]>,
}

fn subsets(model: &mut Model, boundary: &Boundary, size: usize) -> Result<Vec<Vec<Subset>>> {
    let rank = boundary.rank();
    let mut groups: Vec<Vec<Subset>> = (0..=size * (rank - size)).map(|_| Vec::new()).collect();
    let mut indices: Vec<_> = (0..size).collect();
    let mut cells = 0;
    loop {
        model.spend(size)?;
        cells += 2 * size + 4;
        if cells > MAX_HORN_CELLS {
            return Err(Error::plain("Horn subset search exceeds the storage limit"));
        }
        // The indices refer to descending weights. The associated partition,
        // stored in increasing order, has entries i_a-a in one-based notation.
        let shape: Rc<[_]> = indices
            .iter()
            .enumerate()
            .map(|(position, &index)| (index - position) as i128)
            .collect();
        let degree = shape.iter().sum::<i128>() as usize;
        let [lambda, mu, nu] = &boundary.0;
        let [lambda, mu, nu] = [lambda, mu, nu]
            .map(|weight| sum(indices.iter().map(|&index| weight[rank - index - 1])));
        groups[degree].push(Subset {
            indices: indices.iter().map(|&index| index + 1).collect(),
            weight: [lambda?, mu?, nu?],
            shape,
        });
        let Some(position) = (0..size)
            .rev()
            .find(|&position| indices[position] < rank - size + position)
        else {
            break;
        };
        indices[position] += 1;
        for next in position + 1..size {
            indices[next] = indices[next - 1] + 1;
        }
    }
    Ok(groups)
}

pub(crate) fn horn(model: &mut Model, boundary: &Boundary) -> Result<Horn> {
    let triangle = Triangle::new(model, boundary)?;
    if !triangle.balanced() {
        return Ok(Horn::Total {
            expected: triangle.expected,
            actual: triangle.actual,
        });
    }
    for size in 1..boundary.rank() {
        let groups = subsets(model, boundary, size)?;
        for (degree_i, group_i) in groups.iter().enumerate() {
            for (degree_j, group_j) in groups.iter().take(groups.len() - degree_i).enumerate() {
                let group_k = &groups[degree_i + degree_j];
                for i in group_i {
                    for j in group_j {
                        model.spend(1)?;
                        let left = add(i.weight[0], j.weight[1])?;
                        for k in group_k {
                            model.spend(1)?;
                            let right = k.weight[2];
                            if left >= right {
                                continue;
                            }
                            let test =
                                Boundary([i.shape.clone(), j.shape.clone(), k.shape.clone()]);
                            if witness(model, &test)?.is_some() {
                                return Ok(Horn::Violation {
                                    indices: [
                                        i.indices.clone(),
                                        j.indices.clone(),
                                        k.indices.clone(),
                                    ],
                                    left,
                                    right,
                                });
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(Horn::Admissible)
}
