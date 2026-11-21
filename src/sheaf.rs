use std::cell::RefCell;
use std::collections::{btree_map::Entry, BTreeMap};
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Poset, TypeId, Value};
use crate::odd_cobar::Prime;
use crate::odd_ext;
use crate::prime_linear::{self, Basis, Coordinate, QuotientBasis};
use crate::syntax::Label;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) enum Form {
    Cochain,
    Class,
}

impl Form {
    pub(crate) fn total_name(self) -> &'static str {
        match self {
            Self::Cochain => "ShTot",
            Self::Class => "ShHyper",
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Cochain => "ShC",
            Self::Class => "ShH",
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Differential,
    Class,
    Representative,
    Scale(i64),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Differential => "shD",
            Self::Class => "shClass",
            Self::Representative => "shRep",
            Self::Scale(_) => "shScale",
        }
    }
}

pub(crate) struct Presentation {
    pub(crate) stalks: Vec<(Label, usize)>,
    pub(crate) maps: Vec<Cover>,
}

pub(crate) struct Cover {
    pub(crate) source: Label,
    pub(crate) target: Label,
    pub(crate) entries: Vec<i64>,
}

#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Generator {
    pub(crate) vertical: usize,
    pub(crate) chain: Rc<[usize]>,
    pub(crate) index: usize,
}

pub(crate) type Term = (Generator, u32);
pub(crate) type Support = Rc<[Term]>;

pub(crate) struct Sheaf {
    pub(crate) prime: Prime,
    base: TypeId,
    pub(crate) poset: Rc<Poset>,
    dims: Vec<usize>,
    maps: Vec<Option<Rc<Matrix>>>,
    lower: Vec<Vec<usize>>,
    covers: Vec<Vec<usize>>,
    longest: Vec<usize>,
    pub(crate) height: usize,
    pub(crate) cohomology: Cohomology,
}

#[derive(Default, PartialEq, Eq)]
pub(crate) struct Matrix(Vec<(usize, Box<[Coordinate]>)>);

impl Matrix {
    pub(crate) fn literal(
        model: &mut Model,
        prime: Prime,
        width: usize,
        values: Vec<i64>,
    ) -> Result<Self> {
        model.spend(values.len())?;
        let mut columns: BTreeMap<usize, Vec<Coordinate>> = BTreeMap::new();
        let mut cells = 1;
        // An empty input has no entries even when the matrix has width zero.
        for (index, value) in values.into_iter().enumerate() {
            let coefficient = residue(prime, value);
            if coefficient != 0 {
                let column = index % width;
                let row = index / width;
                let rows = match columns.entry(column) {
                    Entry::Occupied(slot) => slot.into_mut(),
                    Entry::Vacant(slot) => {
                        cells += 3;
                        slot.insert(Vec::new())
                    }
                };
                cells += 2;
                check_cells(cells)?;
                rows.push((row, coefficient));
            }
        }
        Ok(Self(
            columns
                .into_iter()
                .map(|(i, rows)| (i, rows.into()))
                .collect(),
        ))
    }

    pub(crate) fn column(&self, index: usize) -> &[Coordinate] {
        self.0
            .binary_search_by_key(&index, |&(column, _)| column)
            .map_or(&[], |column| &self.0[column].1)
    }

    pub(crate) fn cells(&self) -> usize {
        1 + self
            .0
            .iter()
            .map(|(_, rows)| 3 + 2 * rows.len())
            .sum::<usize>()
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    pub(crate) fn compose(
        model: &mut Model,
        prime: Prime,
        after: &Self,
        before: &Self,
    ) -> Result<Self> {
        if after.0.is_empty() || before.0.is_empty() {
            return Ok(Self::default());
        }
        let mut columns = Vec::new();
        let mut cells = 1;
        for (index, middle) in &before.0 {
            model.spend(1 + middle.len())?;
            let mut column = Vec::new();
            for &(j, coefficient) in middle {
                column =
                    prime_linear::add_scaled(model, prime, &column, after.column(j), coefficient)?;
            }
            if !column.is_empty() {
                cells += 3 + 2 * column.len();
                check_cells(cells)?;
                columns.push((*index, column.into()));
            }
        }
        Ok(Self(columns))
    }
}

impl Sheaf {
    pub(crate) fn build(
        model: &mut Model,
        base: TypeId,
        characteristic: usize,
        presentation: Presentation,
    ) -> Result<Self> {
        let prime = match characteristic {
            2 => Prime::validated(2),
            n if n < 2 || n % 2 == 0 => {
                return Err(Error::plain(format!(
                    "sheaf characteristic {n} is not prime"
                )));
            }
            n => Prime::new(model, n)?,
        };
        let poset = model.poset(base).expect("named base poset");
        let size = poset.names.len();
        let mut dims = vec![None; size];
        for (label, dimension) in presentation.stalks {
            let point = model.point(base, &label)?;
            if dims[point].replace(dimension).is_some() {
                return Err(Error::new(
                    label.span,
                    format!("duplicate stalk for '{}'", label.text),
                ));
            }
        }
        let dims: Vec<_> = dims
            .into_iter()
            .enumerate()
            .map(|(point, dimension)| {
                dimension.ok_or_else(|| {
                    Error::plain(format!(
                        "missing stalk dimension for '{}'",
                        poset.names[point]
                    ))
                })
            })
            .collect::<Result<_>>()?;

        let order = PosetData::build(model, &poset)?;
        let mut given = vec![None; size * size];
        let empty = Rc::new(Matrix::default());
        let mut map_cells = 1;
        for map in presentation.maps {
            let source = model.point(base, &map.source)?;
            let target = model.point(base, &map.target)?;
            if order.covers[source].binary_search(&target).is_err() {
                return Err(Error::new(
                    map.source.span,
                    format!(
                        "'{}<{}' is not a Hasse cover",
                        map.source.text, map.target.text
                    ),
                ));
            }
            let slot = &mut given[source * size + target];
            if slot.is_some() {
                return Err(Error::new(
                    map.source.span,
                    format!(
                        "duplicate restriction on cover '{}<{}'",
                        map.source.text, map.target.text
                    ),
                ));
            }
            let expected = dims[source].checked_mul(dims[target]).ok_or_else(|| {
                Error::new(
                    map.target.span,
                    "sheaf restriction matrix size is out of range",
                )
            })?;
            if map.entries.len() != expected {
                return Err(Error::new(map.target.span, format!(
                    "restriction '{}<{}' requires {} rows and {} columns ({expected} entries), found {} entries",
                    map.source.text, map.target.text, dims[target], dims[source], map.entries.len()
                )));
            }
            let matrix = Matrix::literal(model, prime, dims[source], map.entries)?;
            if matrix.0.is_empty() {
                *slot = Some(empty.clone());
            } else {
                map_cells += matrix.cells();
                check_cells(map_cells)?;
                *slot = Some(Rc::new(matrix));
            }
        }
        for (source, targets) in order.covers.iter().enumerate() {
            for &target in targets {
                if given[source * size + target].is_none() {
                    return Err(Error::plain(format!(
                        "missing restriction on cover '{}<{}'",
                        poset.names[source], poset.names[target]
                    )));
                }
            }
        }

        let maps = compose_paths(model, prime, &poset, &order, &given, &empty, &mut map_cells)?;
        let height = dims
            .iter()
            .zip(&order.longest)
            .filter_map(|(&dimension, &length)| (dimension != 0).then_some(length))
            .max()
            .unwrap_or(0);
        model.allocate(
            4 * size
                + maps.len()
                + map_cells
                + order
                    .lower
                    .iter()
                    .chain(&order.covers)
                    .map(Vec::len)
                    .sum::<usize>(),
        )?;
        Ok(Self {
            prime,
            base,
            poset,
            dims,
            maps,
            lower: order.lower,
            covers: order.covers,
            longest: order.longest,
            height,
            cohomology: Cohomology::default(),
        })
    }

    pub(crate) fn term(
        &self,
        model: &mut Model,
        chain: &[Label],
        index: usize,
        coefficient: i64,
    ) -> Result<Value> {
        if chain.is_empty() {
            return Err(Error::plain("shTerm requires a nonempty strict chain"));
        }
        model.spend(chain.len())?;
        let points: Vec<_> = chain
            .iter()
            .map(|p| model.point(self.base, p))
            .collect::<Result<_>>()?;
        for (adjacent, labels) in points.windows(2).zip(chain.windows(2)) {
            if !self.poset.strict_less(adjacent[0], adjacent[1]) {
                return Err(Error::new(
                    labels[1].span,
                    "shTerm requires a strict increasing chain",
                ));
            }
        }
        let last = *points.last().unwrap();
        if index >= self.dims[last] {
            return Err(Error::plain(format!(
                "terminal stalk coordinate {index} is out of range for '{}' of dimension {}",
                self.poset.names[last], self.dims[last]
            )));
        }
        let coefficient = residue(self.prime, coefficient);
        if coefficient == 0 {
            return Ok(Value::SheafC([].into()));
        }
        persist(
            model,
            vec![(
                Generator {
                    vertical: 0,
                    chain: points.into(),
                    index,
                },
                coefficient,
            )],
        )
    }

    pub(crate) fn base(&self) -> TypeId {
        self.base
    }

    pub(crate) fn dims(&self) -> &[usize] {
        &self.dims
    }

    pub(crate) fn covers(&self) -> &[Vec<usize>] {
        &self.covers
    }

    pub(crate) fn restriction(&self, source: usize, target: usize) -> &Matrix {
        self.maps[source * self.dims.len() + target]
            .as_deref()
            .expect("comparable map")
    }

    pub(crate) fn rank(&self, model: &mut Model, degree: usize) -> Result<usize> {
        if degree >= self.height {
            return Ok(0);
        }
        self.cohomology.rank(model, self, degree)
    }
    fn chains(
        &self,
        model: &mut Model,
        remaining: usize,
        dimension: usize,
        reversed: &mut Vec<usize>,
        output: &mut Layer,
        limit: usize,
    ) -> Result<()> {
        let point = *reversed.last().unwrap();
        if remaining == 0 {
            output.size = output
                .size
                .checked_add(dimension)
                .filter(|&size| size <= limit)
                .ok_or_else(|| layer_limit(limit))?;
            output.cells += 4 + reversed.len();
            check_cells(output.cells)?;
            model.spend(1 + reversed.len())?;
            output.blocks.push(Block {
                vertical: 0,
                chain: reversed.iter().rev().copied().collect(),
                start: 0,
            });
            return Ok(());
        }
        model.spend(1 + self.lower[point].len())?;
        for &lower in &self.lower[point] {
            if self.longest[lower] >= remaining {
                reversed.push(lower);
                self.chains(model, remaining - 1, dimension, reversed, output, limit)?;
                reversed.pop();
            }
        }
        Ok(())
    }
}

impl Cochains for Sheaf {
    fn prime(&self) -> Prime {
        self.prime
    }

    fn cofaces(
        &self,
        model: &mut Model,
        generator: &Generator,
        coefficient: u32,
        emit: &mut Sink<'_>,
    ) -> Result<()> {
        let chain = &generator.chain;
        if chain.len() >= self.height {
            return Ok(());
        }
        let size = self.dims.len();
        for slot in 0..=chain.len() {
            model.spend(size)?;
            let factor = if slot % 2 == 0 {
                coefficient
            } else {
                self.prime.negate(coefficient)
            };
            for point in 0..size {
                if (slot > 0 && !self.poset.strict_less(chain[slot - 1], point))
                    || (slot < chain.len() && !self.poset.strict_less(point, chain[slot]))
                {
                    continue;
                }
                if slot < chain.len() {
                    let chain = inserted(model, chain, slot, point)?;
                    emit(
                        model,
                        Generator {
                            vertical: generator.vertical,
                            chain,
                            index: generator.index,
                        },
                        factor,
                    )?;
                } else {
                    let source = chain[slot - 1];
                    let map = self.maps[source * size + point]
                        .as_ref()
                        .expect("comparable map");
                    let column = map.column(generator.index);
                    if column.is_empty() {
                        continue;
                    }
                    let chain = inserted(model, chain, slot, point)?;
                    for &(index, scalar) in column {
                        emit(
                            model,
                            Generator {
                                vertical: generator.vertical,
                                chain: chain.clone(),
                                index,
                            },
                            self.prime.mul(factor, scalar),
                        )?;
                    }
                }
            }
        }
        Ok(())
    }

    fn layer(&self, model: &mut Model, degree: usize, limit: usize) -> Result<Layer> {
        if degree >= self.height {
            return Ok(Layer::default());
        }
        let mut output = Layer::default();
        let mut reversed = Vec::new();
        for (end, &dimension) in self.dims.iter().enumerate() {
            if dimension != 0 && self.longest[end] > degree {
                reversed.push(end);
                self.chains(model, degree, dimension, &mut reversed, &mut output, limit)?;
                reversed.pop();
            }
        }
        model.spend(sort_work(output.blocks.len()).saturating_mul(degree + 1))?;
        output
            .blocks
            .sort_unstable_by(|left, right| left.chain.cmp(&right.chain));
        let mut next = 0;
        for block in &mut output.blocks {
            block.start = next;
            next += self.dims[*block.chain.last().unwrap()];
        }
        Ok(output)
    }
}

pub(crate) type Sink<'a> = dyn FnMut(&mut Model, Generator, u32) -> Result<()> + 'a;

pub(crate) trait Cochains {
    fn prime(&self) -> Prime;

    fn cofaces(
        &self,
        model: &mut Model,
        generator: &Generator,
        coefficient: u32,
        emit: &mut Sink<'_>,
    ) -> Result<()>;

    fn layer(&self, model: &mut Model, degree: usize, limit: usize) -> Result<Layer>;

    fn differential_support(&self, model: &mut Model, terms: &[Term]) -> Result<Vec<Term>> {
        let mut result = Cochain::default();
        for (generator, coefficient) in terms {
            self.cofaces(
                model,
                generator,
                *coefficient,
                &mut |model, target, scalar| result.add(model, self.prime(), target, scalar),
            )?;
        }
        Ok(result.into_terms())
    }

    fn differential(&self, model: &mut Model, terms: &[Term]) -> Result<Value> {
        let terms = self.differential_support(model, terms)?;
        persist(model, terms)
    }

    fn indexed_column(
        &self,
        model: &mut Model,
        source: &Generator,
        target: &Layer,
    ) -> Result<Vec<Coordinate>> {
        if target.size == 0 {
            return Ok(Vec::new());
        }
        let mut column = Vec::new();
        self.cofaces(model, source, 1, &mut |model, generator, coefficient| {
            check_cells(2 * (column.len() + 1))?;
            column.push((target.index(model, &generator)?, coefficient));
            Ok(())
        })?;
        model.spend(sort_work(column.len()))?;
        column.sort_unstable_by_key(|&(index, _)| index);
        Ok(column)
    }

    fn elements(&self, model: &mut Model, degree: usize, limit: usize) -> Result<Vec<Value>> {
        let order = self.prime().get() as usize;
        let mut max_rank = 0;
        let mut cardinality = 1;
        while cardinality <= limit / order {
            cardinality *= order;
            max_rank += 1;
        }
        let layer = self.layer(model, degree, max_rank)?;
        let size = order.pow(layer.size as u32);
        model.spend(size * (1 + layer.size))?;
        model.allocate(size)?;
        let mut values = Vec::with_capacity(size);
        for mut choices in 0..size {
            let mut support = Vec::new();
            for index in 0..layer.size {
                let coefficient = (choices % order) as u32;
                choices /= order;
                if coefficient != 0 {
                    support.push((layer.generator(index), coefficient));
                }
            }
            values.push(persist(model, support)?);
        }
        Ok(values)
    }
}

struct PosetData {
    order: Vec<usize>,
    covers: Vec<Vec<usize>>,
    lower: Vec<Vec<usize>>,
    longest: Vec<usize>,
}

impl PosetData {
    fn build(model: &mut Model, poset: &Poset) -> Result<Self> {
        let n = poset.names.len();
        let words = n.div_ceil(u64::BITS as usize);
        let mut upper_masks = vec![vec![0u64; words]; n];
        let mut lower_masks = vec![vec![0u64; words]; n];
        let mut lower = vec![Vec::new(); n];
        model.spend(n * n)?;
        for (a, upper) in upper_masks.iter_mut().enumerate() {
            for (b, below) in lower.iter_mut().enumerate() {
                if poset.strict_less(a, b) {
                    upper[b / 64] |= 1u64 << (b % 64);
                    lower_masks[b][a / 64] |= 1u64 << (a % 64);
                    below.push(a);
                }
            }
        }
        let mut covers = vec![Vec::new(); n];
        for (b, predecessors) in lower.iter().enumerate() {
            for &a in predecessors {
                model.spend(1 + words)?;
                if upper_masks[a]
                    .iter()
                    .zip(&lower_masks[b])
                    .all(|(&x, &y)| x & y == 0)
                {
                    covers[a].push(b);
                }
            }
        }
        let mut order: Vec<_> = (0..n).collect();
        order.sort_unstable_by_key(|&point| (lower[point].len(), point));
        let mut longest = vec![1; n];
        for &x in &order {
            for &y in &covers[x] {
                longest[y] = longest[y].max(1 + longest[x]);
            }
        }
        Ok(Self {
            order,
            covers,
            lower,
            longest,
        })
    }
}

fn compose_paths(
    model: &mut Model,
    prime: Prime,
    poset: &Poset,
    order: &PosetData,
    given: &[Option<Rc<Matrix>>],
    empty: &Rc<Matrix>,
    cells: &mut usize,
) -> Result<Vec<Option<Rc<Matrix>>>> {
    let size = order.order.len();
    let mut maps: Vec<Option<Rc<Matrix>>> = vec![None; size * size];
    for &source in &order.order {
        model.spend(size)?;
        let row = &mut maps[source * size..(source + 1) * size];
        for &mid in &order.order {
            let previous = row[mid].clone();
            if mid != source && previous.is_none() {
                continue;
            }
            for &target in &order.covers[mid] {
                model.spend(1)?;
                let edge = given[mid * size + target]
                    .as_ref()
                    .expect("validated cover map");
                if let Some(previous) = previous.as_ref() {
                    let inferred = Matrix::compose(model, prime, edge, previous)?;
                    if let Some(existing) = &row[target] {
                        model.spend(existing.cells() + inferred.cells())?;
                        if **existing != inferred {
                            return Err(Error::plain(format!(
                                "restriction paths from '{}' to '{}' do not commute",
                                poset.names[source], poset.names[target]
                            )));
                        }
                    } else if inferred.0.is_empty() {
                        row[target] = Some(empty.clone());
                    } else {
                        *cells += inferred.cells();
                        check_cells(*cells)?;
                        row[target] = Some(Rc::new(inferred));
                    }
                } else {
                    // The identity at the source stays implicit.
                    row[target] = Some(edge.clone());
                }
            }
        }
    }
    Ok(maps)
}

#[derive(Default)]
pub(crate) struct Layer {
    blocks: Vec<Block>,
    pub(crate) size: usize,
    cells: usize,
}

struct Block {
    vertical: usize,
    chain: Rc<[usize]>,
    start: usize,
}

struct Degree {
    layer: Layer,
    classes: QuotientBasis,
}

impl Layer {
    pub(crate) fn append_layer(
        &mut self,
        model: &mut Model,
        mut layer: Layer,
        q: usize,
        limit: usize,
    ) -> Result<()> {
        let size = self
            .size
            .checked_add(layer.size)
            .filter(|&size| size <= limit)
            .ok_or_else(|| layer_limit(limit))?;
        check_cells(self.cells + layer.cells)?;
        model.spend(layer.blocks.len())?;
        for block in &mut layer.blocks {
            block.vertical = q;
            block.start += self.size;
        }
        self.size = size;
        self.cells += layer.cells;
        self.blocks.extend(layer.blocks);
        Ok(())
    }

    pub(crate) fn cells(&self) -> usize {
        self.cells
    }

    pub(crate) fn vertical_start(&self, q: usize) -> usize {
        let index = self.blocks.partition_point(|block| block.vertical < q);
        self.blocks
            .get(index)
            .map_or(self.size, |block| block.start)
    }

    pub(crate) fn generator(&self, index: usize) -> Generator {
        let block = &self.blocks[self.blocks.partition_point(|block| block.start <= index) - 1];
        Generator {
            vertical: block.vertical,
            chain: block.chain.clone(),
            index: index - block.start,
        }
    }

    pub(crate) fn index(&self, model: &mut Model, generator: &Generator) -> Result<usize> {
        model.spend(1 + generator.chain.len() * (1 + self.blocks.len().max(1).ilog2() as usize))?;
        let block = self
            .blocks
            .binary_search_by(|block| {
                (block.vertical, &block.chain).cmp(&(generator.vertical, &generator.chain))
            })
            .map_err(|_| Error::plain("sheaf cochain chain is outside its degree"))?;
        Ok(self.blocks[block].start + generator.index)
    }
}

#[derive(Default)]
pub(crate) struct Cohomology {
    degrees: RefCell<BTreeMap<usize, Rc<Degree>>>,
}

impl Cohomology {
    pub(crate) fn rank(
        &self,
        model: &mut Model,
        complex: &(impl Cochains + ?Sized),
        degree: usize,
    ) -> Result<usize> {
        Ok(self.degree(model, complex, degree)?.classes.rank())
    }

    pub(crate) fn generator(
        &self,
        model: &mut Model,
        complex: &(impl Cochains + ?Sized),
        degree: usize,
        index: usize,
    ) -> Result<Value> {
        let rank = self.rank(model, complex, degree)?;
        if index >= rank {
            return Err(Error::plain(format!(
                "sheaf cohomology generator {index} is out of range for rank {rank}"
            )));
        }
        odd_ext::persist(model, vec![(index, 1)])
    }

    pub(crate) fn class(
        &self,
        model: &mut Model,
        complex: &(impl Cochains + ?Sized),
        degree: usize,
        raw: &[Term],
    ) -> Result<Value> {
        if raw.is_empty() {
            return Ok(Value::OddClass([].into()));
        }
        let basis = self.degree(model, complex, degree)?;
        let raw = raw
            .iter()
            .map(|(generator, coefficient)| {
                Ok((basis.layer.index(model, generator)?, *coefficient))
            })
            .collect::<Result<_>>()?;
        let coordinates = basis
            .classes
            .project(model, complex.prime(), raw)?
            .ok_or_else(|| Error::plain("shClass requires a closed sheaf cochain"))?;
        odd_ext::persist(model, coordinates)
    }

    pub(crate) fn representative_support(
        &self,
        model: &mut Model,
        complex: &(impl Cochains + ?Sized),
        degree: usize,
        class: &[Coordinate],
    ) -> Result<Vec<Term>> {
        if class.is_empty() {
            return Ok(Vec::new());
        }
        let basis = self.degree(model, complex, degree)?;
        let raw = basis
            .classes
            .representative(model, complex.prime(), class)?;
        Ok(raw
            .into_iter()
            .map(|(index, coefficient)| (basis.layer.generator(index), coefficient))
            .collect())
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        complex: &(impl Cochains + ?Sized),
        degree: usize,
        class: &[Coordinate],
    ) -> Result<Value> {
        let raw = self.representative_support(model, complex, degree, class)?;
        persist(model, raw)
    }

    fn degree(
        &self,
        model: &mut Model,
        complex: &(impl Cochains + ?Sized),
        degree: usize,
    ) -> Result<Rc<Degree>> {
        if let Some(layer) = self.degrees.borrow().get(&degree) {
            return Ok(layer.clone());
        }
        let layer = complex.layer(model, degree, MAX_PIVOT_CELLS)?;
        let prime = complex.prime();
        let classes = if layer.size == 0 {
            QuotientBasis::default()
        } else {
            let mut boundaries = Basis::default();
            if degree > 0 {
                let previous = complex.layer(model, degree - 1, MAX_PIVOT_CELLS)?;
                check_cells(layer.cells + previous.cells)?;
                for index in 0..previous.size {
                    let column =
                        complex.indexed_column(model, &previous.generator(index), &layer)?;
                    boundaries.insert(model, prime, column)?;
                }
            }
            let next = complex.layer(model, degree + 1, MAX_PIVOT_CELLS)?;
            check_cells(layer.cells + next.cells)?;
            QuotientBasis::build(model, prime, boundaries, layer.size, |model, index| {
                complex.indexed_column(model, &layer.generator(index), &next)
            })?
        };
        model.allocate(3 + layer.cells)?;
        let result = Rc::new(Degree { layer, classes });
        self.degrees.borrow_mut().insert(degree, result.clone());
        Ok(result)
    }
}

#[derive(Default)]
pub(crate) struct Cochain {
    terms: BTreeMap<Generator, u32>,
    cells: usize,
}

impl Cochain {
    pub(crate) fn into_terms(self) -> Vec<Term> {
        self.terms.into_iter().collect()
    }

    pub(crate) fn add(
        &mut self,
        model: &mut Model,
        prime: Prime,
        generator: Generator,
        coefficient: u32,
    ) -> Result<()> {
        if coefficient == 0 {
            return Ok(());
        }
        let term_cells = 5 + generator.chain.len();
        model.spend(term_cells * (1 + self.terms.len().max(1).ilog2() as usize))?;
        match self.terms.entry(generator) {
            Entry::Vacant(slot) => {
                self.cells += term_cells;
                check_cells(self.cells)?;
                slot.insert(coefficient);
            }
            Entry::Occupied(mut slot) => {
                let value = prime.add(*slot.get(), coefficient);
                if value == 0 {
                    self.cells -= term_cells;
                    slot.remove();
                } else {
                    *slot.get_mut() = value;
                }
            }
        }
        Ok(())
    }
}

pub(crate) fn add(model: &mut Model, prime: Prime, left: &[Term], right: &[Term]) -> Result<Value> {
    let mut result = Cochain::default();
    for (generator, coefficient) in left.iter().chain(right) {
        result.add(model, prime, generator.clone(), *coefficient)?;
    }
    persist(model, result.into_terms())
}

pub(crate) fn scale(model: &mut Model, prime: Prime, scalar: u32, terms: &[Term]) -> Result<Value> {
    if scalar == 0 {
        return Ok(Value::SheafC([].into()));
    }
    model.spend(terms.len())?;
    persist(
        model,
        terms
            .iter()
            .map(|(generator, coefficient)| (generator.clone(), prime.mul(scalar, *coefficient)))
            .collect(),
    )
}

pub(crate) fn compare(
    model: &mut Model,
    left: &[Term],
    right: &[Term],
) -> Result<std::cmp::Ordering> {
    // Field digits at larger chain coordinates are more significant.
    for (a, b) in left.iter().rev().zip(right.iter().rev()) {
        model.spend(1 + a.0.chain.len().min(b.0.chain.len()))?;
        let comparison = a.cmp(b);
        if !comparison.is_eq() {
            return Ok(comparison);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

pub(crate) fn residue(prime: Prime, value: i64) -> u32 {
    value.rem_euclid(i64::from(prime.get())) as u32
}

fn inserted(model: &mut Model, chain: &[usize], slot: usize, point: usize) -> Result<Rc<[usize]>> {
    model.spend(1 + chain.len())?;
    let mut output = Vec::with_capacity(chain.len() + 1);
    output.extend_from_slice(&chain[..slot]);
    output.push(point);
    output.extend_from_slice(&chain[slot..]);
    Ok(output.into())
}

pub(crate) fn persist(model: &mut Model, terms: Vec<Term>) -> Result<Value> {
    let cells = terms
        .iter()
        .map(|(generator, _)| 5 + generator.chain.len())
        .sum();
    check_cells(cells)?;
    model.allocate(cells)?;
    Ok(Value::SheafC(terms.into()))
}

fn layer_limit(limit: usize) -> Error {
    if limit == MAX_PIVOT_CELLS {
        Error::plain("sheaf cochain rank exceeds its storage limit")
    } else {
        Error::plain("sheaf cochains exceed the finite type enumeration limit")
    }
}

pub(crate) fn check_cells(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain("sheaf calculation exceeds the storage limit"));
    }
    Ok(())
}

fn sort_work(len: usize) -> usize {
    len * (1 + len.max(1).ilog2() as usize)
}
