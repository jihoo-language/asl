use std::collections::{btree_map::Entry, BTreeMap};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::sheaf::{
    self, Cochain, Cochains, Cohomology, Generator, Layer, Place, Sheaf, Sink, Term,
};
use crate::sheaf_hyper::Map;
use crate::syntax::Label;

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Ends {
    pub(crate) source: Rc<str>,
    pub(crate) target: Rc<str>,
}

impl Ends {
    pub(crate) fn endomorphisms(name: Rc<str>) -> Self {
        Self {
            source: name.clone(),
            target: name,
        }
    }
}

pub(crate) struct Hom {
    pub(crate) source: Rc<Sheaf>,
    pub(crate) target: Rc<Sheaf>,
    pub(crate) cohomology: Cohomology,
    longest: Vec<usize>,
    pub(crate) height: usize,
}

impl Hom {
    pub(crate) fn build(model: &mut Model, source: Rc<Sheaf>, target: Rc<Sheaf>) -> Result<Self> {
        if source.base() != target.base() || source.prime != target.prime {
            return Err(Error::plain(
                "derived sheaf Hom requires the same named base and prime field",
            ));
        }
        let size = source.dims().len();
        // Ignore chains whose first stalk is zero, without enumerating those
        // chains when their intermediate or terminal stalks are nonzero.
        let mut longest = vec![0; size];
        let mut height = 0;
        for &point in source.topological() {
            model.spend(1 + source.lower(point).len())?;
            let prefix = source
                .lower(point)
                .iter()
                .map(|&lower| longest[lower])
                .max()
                .unwrap_or(0);
            longest[point] = if prefix > 0 || source.dims()[point] != 0 {
                prefix + 1
            } else {
                0
            };
            if target.dims()[point] != 0 {
                height = height.max(longest[point]);
            }
        }
        model.allocate(1 + size)?;
        Ok(Self {
            source,
            target,
            cohomology: Cohomology::default(),
            longest,
            height,
        })
    }

    pub(crate) fn decode(&self, generator: &Generator) -> (usize, usize) {
        let rows = self.target.dims()[*generator.place.chain().last().unwrap()];
        (generator.index / rows, generator.index % rows)
    }

    fn dimension(&self, first: usize, last: usize) -> Result<usize> {
        self.source.dims()[first]
            .checked_mul(self.target.dims()[last])
            .ok_or_else(matrix_size)
    }

    fn chains(
        &self,
        model: &mut Model,
        remaining: usize,
        path: &mut Vec<usize>,
        layer: &mut Layer,
        limit: usize,
    ) -> Result<()> {
        let point = *path.last().unwrap();
        if remaining == 0 {
            let dimension = self.dimension(point, path[0])?;
            if dimension != 0 {
                layer.push_block(
                    model,
                    path.iter().rev().copied().collect(),
                    dimension,
                    limit,
                )?;
            }
            return Ok(());
        }
        model.spend(1 + self.source.lower(point).len())?;
        for &lower in self.source.lower(point) {
            if self.longest[lower] >= remaining {
                path.push(lower);
                self.chains(model, remaining - 1, path, layer, limit)?;
                path.pop();
            }
        }
        Ok(())
    }

    pub(crate) fn term(
        &self,
        model: &mut Model,
        chain: &[Label],
        column: usize,
        row: usize,
        coefficient: i64,
    ) -> Result<Value> {
        if chain.is_empty() {
            return Err(Error::plain("shHomTerm requires a nonempty strict chain"));
        }
        model.spend(chain.len())?;
        let points: Vec<_> = chain
            .iter()
            .map(|label| model.point(self.source.base(), label))
            .collect::<Result<_>>()?;
        for (pair, labels) in points.windows(2).zip(chain.windows(2)) {
            if !self.source.poset.strict_less(pair[0], pair[1]) {
                return Err(Error::new(
                    labels[1].span,
                    "shHomTerm requires a strict increasing chain",
                ));
            }
        }
        let rows = self.target.dims()[*points.last().unwrap()];
        if column >= self.source.dims()[points[0]] || row >= rows {
            return Err(Error::plain(
                "shHomTerm matrix coordinate is outside the endpoint stalks",
            ));
        }
        let scalar = sheaf::residue(self.prime(), coefficient);
        if scalar == 0 {
            return Ok(Value::SheafC([].into()));
        }
        sheaf::persist(
            model,
            vec![(coordinate(points.into(), rows, column, row)?, scalar)],
        )
    }
}

fn matrix_size() -> Error {
    Error::plain("derived sheaf Hom matrix size is out of range")
}

pub(crate) fn coordinate(
    chain: Rc<[usize]>,
    rows: usize,
    column: usize,
    row: usize,
) -> Result<Generator> {
    let index = column
        .checked_mul(rows)
        .and_then(|value| value.checked_add(row))
        .ok_or_else(matrix_size)?;
    Ok(Generator {
        vertical: 0,
        place: Place::Chain(chain),
        index,
    })
}

impl Cochains for Hom {
    fn prime(&self) -> Prime {
        self.source.prime
    }

    fn cofaces(
        &self,
        model: &mut Model,
        generator: &Generator,
        coefficient: u32,
        emit: &mut Sink<'_>,
    ) -> Result<()> {
        let chain = generator.place.chain();
        if chain.len() >= self.height {
            return Ok(());
        }
        let (column, row) = self.decode(generator);
        let rows = self.target.dims()[*chain.last().unwrap()];
        let prime = self.prime();
        for slot in 0..=chain.len() {
            model.spend(self.source.dims().len())?;
            let scalar = if slot % 2 == 0 {
                coefficient
            } else {
                prime.negate(coefficient)
            };
            for point in 0..self.source.dims().len() {
                if (slot > 0 && !self.source.poset.strict_less(chain[slot - 1], point))
                    || (slot < chain.len() && !self.source.poset.strict_less(point, chain[slot]))
                {
                    continue;
                }
                if slot == 0 {
                    let source = self.source.restriction(point, chain[0]);
                    let action = source.row(model, column)?;
                    if action.is_empty() {
                        continue;
                    }
                    let next = sheaf::inserted(model, chain, slot, point)?;
                    for (new_column, factor) in action {
                        emit(
                            model,
                            coordinate(next.clone(), rows, new_column, row)?,
                            prime.mul(scalar, factor),
                        )?;
                    }
                } else if slot == chain.len() {
                    let target = self.target.restriction(chain[slot - 1], point);
                    let action = target.column(row);
                    if action.is_empty() {
                        continue;
                    }
                    let next = sheaf::inserted(model, chain, slot, point)?;
                    for &(new_row, factor) in action {
                        emit(
                            model,
                            coordinate(next.clone(), self.target.dims()[point], column, new_row)?,
                            prime.mul(scalar, factor),
                        )?;
                    }
                } else {
                    let next = sheaf::inserted(model, chain, slot, point)?;
                    emit(
                        model,
                        Generator {
                            place: Place::Chain(next),
                            ..generator.clone()
                        },
                        scalar,
                    )?;
                }
            }
        }
        Ok(())
    }

    fn layer(&self, model: &mut Model, degree: usize, limit: usize) -> Result<Layer> {
        let mut output = Layer::default();
        if degree >= self.height {
            return Ok(output);
        }
        let mut path = Vec::new();
        for (end, &rows) in self.target.dims().iter().enumerate() {
            if rows != 0 && self.longest[end] > degree {
                path.push(end);
                self.chains(model, degree, &mut path, &mut output, limit)?;
                path.pop();
            }
        }
        output.order_blocks(model, degree, |chain| {
            self.source.dims()[chain[0]] * self.target.dims()[*chain.last().unwrap()]
        })?;
        Ok(output)
    }
}

pub(crate) fn identity(model: &mut Model, data: &Sheaf) -> Result<Value> {
    let mut result = Cochain::default();
    for (point, &dimension) in data.dims().iter().enumerate() {
        let chain = [point].into();
        model.spend(dimension)?;
        for index in 0..dimension {
            result.add(
                model,
                data.prime,
                coordinate(Rc::clone(&chain), dimension, index, index)?,
                1,
            )?;
        }
    }
    sheaf::persist(model, result.into_terms())
}

pub(crate) fn map_cochain(model: &mut Model, map: &Map) -> Result<Value> {
    let mut result = Cochain::default();
    let prime = map.source.sheaf.prime;
    for (point, &rows) in map.target.sheaf.dims().iter().enumerate() {
        let chain = [point].into();
        for (column, entries) in map.matrix(point).columns() {
            model.spend(1 + entries.len())?;
            for &(row, scalar) in entries {
                result.add(
                    model,
                    prime,
                    coordinate(Rc::clone(&chain), rows, column, row)?,
                    scalar,
                )?;
            }
        }
    }
    sheaf::persist(model, result.into_terms())
}

pub(crate) fn compose(
    model: &mut Model,
    outer: (&Hom, usize, &[Term]),
    inner: (&Hom, usize, &[Term]),
) -> Result<Vec<Term>> {
    let (a, m, after) = outer;
    let (b, n, before) = inner;
    if after.is_empty() || before.is_empty() {
        return Ok(Vec::new());
    }
    let prime = a.prime();
    let sign = if m & n & 1 == 1 { prime.negate(1) } else { 1 };
    let mut cells = 0;
    let mut matching = BTreeMap::new();
    for (generator, coefficient) in after {
        model.spend(1 + matching.len().max(1).ilog2() as usize)?;
        let (column, row) = a.decode(generator);
        let bucket = match matching.entry((generator.place.chain()[0], column)) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => {
                cells += 5;
                entry.insert(Vec::new())
            }
        };
        cells += 3;
        sheaf::check_cells(cells)?;
        bucket.push((generator, row, *coefficient));
    }
    let mut result = Cochain::default();
    for (generator, coefficient) in before {
        model.spend(1 + matching.len().max(1).ilog2() as usize)?;
        let (column, middle) = b.decode(generator);
        let junction = *generator.place.chain().last().unwrap();
        if let Some(candidates) = matching.get(&(junction, middle)) {
            for &(after, row, factor) in candidates {
                model.spend(generator.place.chain().len() + after.place.chain().len())?;
                let chain = generator
                    .place
                    .chain()
                    .iter()
                    .copied()
                    .chain(after.place.chain().iter().skip(1).copied())
                    .collect();
                let rows = a.target.dims()[*after.place.chain().last().unwrap()];
                let scalar = prime.mul(sign, prime.mul(*coefficient, factor));
                result.add(model, prime, coordinate(chain, rows, column, row)?, scalar)?;
            }
        }
    }
    Ok(result.into_terms())
}
