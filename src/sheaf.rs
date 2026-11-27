use std::cell::RefCell;
use std::collections::{btree_map::Entry, BTreeMap, HashMap};
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
    pub(crate) place: Place,
    pub(crate) index: usize,
}

#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) enum Place {
    Chain(Rc<[usize]>),
    Cell { degree: usize, ordinal: usize },
}

impl Place {
    pub(crate) fn chain(&self) -> &[usize] {
        match self {
            Self::Chain(chain) => chain,
            Self::Cell { .. } => unreachable!("expected an ordinary Roos coordinate"),
        }
    }

    pub(crate) fn degree(&self) -> usize {
        match self {
            Self::Chain(chain) => chain.len() - 1,
            Self::Cell { degree, .. } => *degree,
        }
    }

    fn payload_cells(&self) -> usize {
        match self {
            Self::Chain(chain) => chain.len(),
            Self::Cell { .. } => 2,
        }
    }
}

impl Generator {
    pub(crate) fn cells(&self) -> usize {
        5 + self.place.payload_cells()
    }
}

pub(crate) type Term = (Generator, u32);
pub(crate) type Support = Rc<[Term]>;

pub(crate) struct Sheaf {
    pub(crate) prime: Prime,
    base: TypeId,
    pub(crate) poset: Rc<Poset>,
    dims: Vec<usize>,
    maps: Restrictions,
    order: Rc<PosetData>,
    pub(crate) height: usize,
    pub(crate) cohomology: Cohomology,
}

#[derive(Default, PartialEq, Eq, Hash)]
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

    pub(crate) fn from_columns(columns: Vec<(usize, Box<[Coordinate]>)>) -> Result<Self> {
        let matrix = Self(columns);
        check_cells(matrix.cells())?;
        Ok(matrix)
    }

    pub(crate) fn from_action(
        model: &mut Model,
        width: usize,
        mut image: impl FnMut(&mut Model, usize) -> Result<Vec<Coordinate>>,
    ) -> Result<Self> {
        model.spend(width)?;
        let mut columns = Vec::new();
        let mut cells = 1;
        for index in 0..width {
            let coordinates = image(model, index)?;
            if !coordinates.is_empty() {
                cells += 3 + 2 * coordinates.len();
                check_cells(cells)?;
                columns.push((index, coordinates.into()));
            }
        }
        Self::from_columns(columns)
    }

    pub(crate) fn identity(model: &mut Model, rank: usize) -> Result<Self> {
        check_cells(1usize.saturating_add(5usize.saturating_mul(rank)))?;
        model.spend(rank)?;
        Self::from_columns(
            (0..rank)
                .map(|index| (index, vec![(index, 1)].into()))
                .collect(),
        )
    }

    pub(crate) fn columns(&self) -> impl Iterator<Item = (usize, &[Coordinate])> {
        self.0.iter().map(|(i, rows)| (*i, rows.as_ref()))
    }

    pub(crate) fn row(&self, model: &mut Model, index: usize) -> Result<Vec<Coordinate>> {
        let mut result = Vec::new();
        for (column, rows) in self.columns() {
            model.spend(1 + rows.len().ilog2() as usize)?;
            if let Ok(position) = rows.binary_search_by_key(&index, |&(row, _)| row) {
                result.push((column, rows[position].1));
            }
        }
        Ok(result)
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

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        prime: Prime,
        input: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let mut column = Vec::new();
        for &(index, coefficient) in input {
            column =
                prime_linear::add_scaled(model, prime, &column, self.column(index), coefficient)?;
        }
        Ok(column)
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
            let column = after.apply(model, prime, middle)?;
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
        let prime = field(model, characteristic)?;
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

        let mut table = CoverTable::new(model, base)?;
        for map in presentation.maps {
            let source = model.point(base, &map.source)?;
            let target = model.point(base, &map.target)?;
            let Ok(edge) = table.order.covers[source].binary_search(&target) else {
                return Err(Error::new(
                    map.source.span,
                    format!(
                        "'{}<{}' is not a Hasse cover",
                        map.source.text, map.target.text
                    ),
                ));
            };
            if table.given[source][edge].is_some() {
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
            table.insert(model, source, edge, Rc::new(matrix))?;
        }
        for (source, targets) in table.order.covers.iter().enumerate() {
            for (edge, &target) in targets.iter().enumerate() {
                if table.given[source][edge].is_none() {
                    return Err(Error::plain(format!(
                        "missing restriction on cover '{}<{}'",
                        poset.names[source], poset.names[target]
                    )));
                }
            }
        }

        Self::assemble(model, base, prime, dims, table)
    }

    pub(crate) fn from_cover_maps(
        model: &mut Model,
        base: TypeId,
        prime: Prime,
        dims: Vec<usize>,
        mut make: impl FnMut(&mut Model, usize, usize) -> Result<Rc<Matrix>>,
    ) -> Result<Self> {
        let mut table = CoverTable::new(model, base)?;
        debug_assert_eq!(dims.len(), table.poset.names.len());
        for source in 0..table.order.covers.len() {
            for edge in 0..table.order.covers[source].len() {
                let target = table.order.covers[source][edge];
                let matrix = make(model, source, target)?;
                table.insert(model, source, edge, matrix)?;
            }
        }
        Self::assemble(model, base, prime, dims, table)
    }

    fn assemble(
        model: &mut Model,
        base: TypeId,
        prime: Prime,
        dims: Vec<usize>,
        table: CoverTable,
    ) -> Result<Self> {
        let CoverTable {
            poset,
            order,
            given,
            pool,
        } = table;
        let maps = Restrictions::build(model, prime, &poset, &order, &given, pool)?;
        let size = poset.names.len();
        let height = dims
            .iter()
            .zip(&order.longest)
            .filter_map(|(&dimension, &length)| (dimension != 0).then_some(length))
            .max()
            .unwrap_or(0);
        model.allocate(size + 4)?;
        Ok(Self {
            prime,
            base,
            poset,
            dims,
            maps,
            order,
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
                    place: Place::Chain(points.into()),
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

    pub(crate) fn lower(&self, point: usize) -> &[usize] {
        &self.order.lower[point]
    }

    pub(crate) fn topological(&self) -> &[usize] {
        &self.order.topological
    }

    pub(crate) fn covers(&self) -> &[Vec<usize>] {
        &self.order.covers
    }

    pub(crate) fn restriction(&self, source: usize, target: usize) -> &Matrix {
        debug_assert!(self.poset.strict_less(source, target));
        self.maps
            .at(self.order.position[source], self.order.position[target])
            .as_ref()
    }

    pub(crate) fn restriction_rc(&self, source: usize, target: usize) -> Rc<Matrix> {
        debug_assert!(self.poset.strict_less(source, target));
        self.maps
            .at(self.order.position[source], self.order.position[target])
            .clone()
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
        bounds: Bounds<'_>,
    ) -> Result<()> {
        let point = *reversed.last().unwrap();
        if remaining == 0 {
            if bounds.starts.is_some_and(|starts| !starts[point]) {
                return Ok(());
            }
            return output.push_block(
                model,
                reversed.iter().rev().copied().collect(),
                dimension,
                bounds.limit,
            );
        }
        model.spend(1 + self.order.lower[point].len())?;
        for &lower in &self.order.lower[point] {
            if bounds.longest[lower] >= remaining {
                reversed.push(lower);
                self.chains(model, remaining - 1, dimension, reversed, output, bounds)?;
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
        let chain = generator.place.chain();
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
                            place: Place::Chain(chain),
                            index: generator.index,
                        },
                        factor,
                    )?;
                } else {
                    let source = chain[slot - 1];
                    let map = self.restriction(source, point);
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
                                place: Place::Chain(chain.clone()),
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
        self.layer_on(
            model,
            degree,
            Bounds {
                longest: &self.order.longest,
                starts: None,
                limit,
            },
            self.height,
        )
    }
}

#[derive(Clone, Copy)]
struct Bounds<'a> {
    longest: &'a [usize],
    starts: Option<&'a [bool]>,
    limit: usize,
}

impl Sheaf {
    fn layer_on(
        &self,
        model: &mut Model,
        degree: usize,
        bounds: Bounds<'_>,
        height: usize,
    ) -> Result<Layer> {
        if degree >= height {
            return Ok(Layer::default());
        }
        let mut output = Layer::default();
        let mut reversed = Vec::new();
        for (end, &dimension) in self.dims.iter().enumerate() {
            if dimension != 0 && bounds.longest[end] > degree {
                reversed.push(end);
                self.chains(model, degree, dimension, &mut reversed, &mut output, bounds)?;
                reversed.pop();
            }
        }
        output.order_blocks(model, degree, |chain| self.dims[*chain.last().unwrap()])?;
        Ok(output)
    }
}

pub(crate) struct Open<'a> {
    sheaf: &'a Sheaf,
    longest: Vec<usize>,
    starts: Option<Rc<[bool]>>,
    pub(crate) height: usize,
}

impl<'a> Open<'a> {
    pub(crate) fn new(model: &mut Model, sheaf: &'a Sheaf, included: &[bool]) -> Result<Self> {
        Self::from_start_set(model, sheaf, included, None)
    }

    pub(crate) fn supported(
        model: &mut Model,
        sheaf: &'a Sheaf,
        starts: Rc<[bool]>,
    ) -> Result<Self> {
        Self::from_start_set(model, sheaf, &starts, Some(starts.clone()))
    }

    fn from_start_set(
        model: &mut Model,
        sheaf: &'a Sheaf,
        included: &[bool],
        starts: Option<Rc<[bool]>>,
    ) -> Result<Self> {
        let mut longest: Vec<_> = included
            .iter()
            .map(|&included| usize::from(included))
            .collect();
        model.spend(longest.len())?;
        for &point in &sheaf.order.topological {
            if longest[point] != 0 {
                model.spend(sheaf.covers()[point].len())?;
                for &next in &sheaf.covers()[point] {
                    // Ordinary opens are upward closed. For supports, propagate the
                    // maximum length from a permitted initial point up the whole open.
                    debug_assert!(starts.is_some() || included[next]);
                    longest[next] = longest[next].max(longest[point] + 1);
                }
            }
        }
        let height = longest
            .iter()
            .zip(&sheaf.dims)
            .filter_map(|(&length, &size)| (size != 0).then_some(length))
            .max()
            .unwrap_or(0);
        Ok(Self {
            sheaf,
            longest,
            starts,
            height,
        })
    }

    pub(crate) fn contains(&self, generator: &Generator) -> bool {
        let start = generator.place.chain()[0];
        self.starts
            .as_ref()
            .map_or(self.longest[start] != 0, |starts| starts[start])
    }
}

impl Cochains for Open<'_> {
    fn prime(&self) -> Prime {
        self.sheaf.prime
    }

    fn cofaces(
        &self,
        model: &mut Model,
        generator: &Generator,
        coefficient: u32,
        emit: &mut Sink<'_>,
    ) -> Result<()> {
        if generator.place.chain().len() >= self.height {
            return Ok(());
        }
        self.sheaf.cofaces(
            model,
            generator,
            coefficient,
            &mut |model, target, scalar| {
                if self.contains(&target) {
                    emit(model, target, scalar)?;
                }
                Ok(())
            },
        )
    }

    fn layer(&self, model: &mut Model, degree: usize, limit: usize) -> Result<Layer> {
        self.sheaf.layer_on(
            model,
            degree,
            Bounds {
                longest: &self.longest,
                starts: self.starts.as_deref(),
                limit,
            },
            self.height,
        )
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

struct CoverTable {
    poset: Rc<Poset>,
    order: Rc<PosetData>,
    given: Vec<Vec<Option<usize>>>,
    pool: MatrixPool,
}

impl CoverTable {
    fn new(model: &mut Model, base: TypeId) -> Result<Self> {
        let poset = model.poset(base).expect("named base poset");
        let order = model.sheaf_site_data(base)?;
        let given = order
            .covers
            .iter()
            .map(|row| vec![None; row.len()])
            .collect();
        Ok(Self {
            poset,
            order,
            given,
            pool: MatrixPool::new(),
        })
    }

    fn insert(
        &mut self,
        model: &mut Model,
        source: usize,
        edge: usize,
        matrix: Rc<Matrix>,
    ) -> Result<()> {
        self.given[source][edge] = Some(self.pool.intern(model, matrix)?);
        Ok(())
    }
}

struct MatrixPool {
    values: Vec<Rc<Matrix>>,
    lookup: HashMap<Rc<Matrix>, usize>,
    cells: usize,
}

impl MatrixPool {
    fn new() -> Self {
        let zero = Rc::new(Matrix::default());
        Self {
            lookup: HashMap::new(),
            values: vec![zero],
            cells: 4,
        }
    }

    fn intern(&mut self, model: &mut Model, matrix: Rc<Matrix>) -> Result<usize> {
        if matrix.is_zero() {
            return Ok(0);
        }
        let cells = matrix.cells();
        model.spend(cells)?;
        if let Some(&index) = self.lookup.get(matrix.as_ref()) {
            return Ok(index);
        }
        self.cells += 3 + cells;
        check_cells(self.cells)?;
        let index = self.values.len();
        self.lookup.insert(matrix.clone(), index);
        self.values.push(matrix);
        Ok(index)
    }
}

pub(crate) struct PosetData {
    topological: Vec<usize>,
    position: Vec<usize>,
    covers: Vec<Vec<usize>>,
    lower: Vec<Vec<usize>>,
    longest: Vec<usize>,
}

impl PosetData {
    pub(crate) fn covers(&self) -> &[Vec<usize>] {
        &self.covers
    }

    pub(crate) fn build(model: &mut Model, poset: &Poset) -> Result<Self> {
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
        model.spend(sort_work(n))?;
        let mut topological: Vec<_> = (0..n).collect();
        topological.sort_unstable_by_key(|&point| (lower[point].len(), point));
        let mut position = vec![0; n];
        let mut longest = vec![1; n];
        for (index, &x) in topological.iter().enumerate() {
            position[x] = index;
            for &y in &covers[x] {
                longest[y] = longest[y].max(1 + longest[x]);
            }
        }
        let cells = 9 * n + lower.iter().chain(&covers).map(Vec::len).sum::<usize>();
        check_cells(cells)?;
        model.allocate(cells)?;
        Ok(Self {
            topological,
            position,
            covers,
            lower,
            longest,
        })
    }
}

struct Restrictions {
    rows: Vec<Box<[TransportRun]>>,
    matrices: Vec<Rc<Matrix>>,
}

struct TransportRun {
    start: usize,
    end: usize,
    matrix: usize,
}

impl Restrictions {
    fn build(
        model: &mut Model,
        prime: Prime,
        poset: &Poset,
        order: &PosetData,
        given: &[Vec<Option<usize>>],
        mut pool: MatrixPool,
    ) -> Result<Self> {
        let size = order.topological.len();
        let mut rows = Vec::with_capacity(size);
        let mut storage = 3 * size;
        let mut row: Vec<Option<usize>> = vec![None; size];
        for (position, &source) in order.topological.iter().enumerate() {
            model.spend(3 * size - position)?;
            row.fill(None);
            for &mid in &order.topological[position..] {
                let previous = row[mid];
                if mid != source && previous.is_none() {
                    continue;
                }
                for (&target, &edge) in order.covers[mid].iter().zip(&given[mid]) {
                    model.spend(1)?;
                    let edge = edge.expect("validated cover map");
                    if let Some(previous) = previous {
                        let inferred = Matrix::compose(
                            model,
                            prime,
                            pool.values[edge].as_ref(),
                            pool.values[previous].as_ref(),
                        )?;
                        // Reached zero maps also certify paths; only unreached points have None.
                        if let Some(existing) = row[target] {
                            let existing = pool.values[existing].as_ref();
                            model.spend(existing.cells() + inferred.cells())?;
                            if *existing != inferred {
                                return Err(Error::plain(format!(
                                    "restriction paths from '{}' to '{}' do not commute",
                                    poset.names[source], poset.names[target]
                                )));
                            }
                        } else {
                            row[target] = Some(pool.intern(model, Rc::new(inferred))?);
                        }
                    } else {
                        // The identity is implicit at the source. Covers have no indirect path.
                        row[target] = Some(edge);
                    }
                }
            }
            let mut runs: Vec<TransportRun> = Vec::new();
            for (index, &target) in order.topological.iter().enumerate() {
                if let Some(matrix) = row[target].filter(|&matrix| matrix != 0) {
                    if let Some(last) = runs
                        .last_mut()
                        .filter(|last| last.matrix == matrix && last.end == index)
                    {
                        last.end += 1;
                    } else {
                        storage += 3;
                        check_cells(storage + pool.cells)?;
                        runs.push(TransportRun {
                            start: index,
                            end: index + 1,
                            matrix,
                        });
                    }
                }
            }
            rows.push(runs.into());
        }
        check_cells(storage + pool.cells)?;
        model.allocate(storage + pool.cells)?;
        Ok(Self {
            rows,
            matrices: pool.values,
        })
    }

    fn at(&self, source: usize, target: usize) -> &Rc<Matrix> {
        let row = &self.rows[source];
        let next = row.partition_point(|run| run.end <= target);
        let matrix = row
            .get(next)
            .filter(|run| run.start <= target)
            .map_or(0, |run| run.matrix);
        &self.matrices[matrix]
    }
}

#[derive(Default)]
pub(crate) struct Layer {
    blocks: Vec<Block>,
    pub(crate) size: usize,
    cells: usize,
}

struct Block {
    vertical: usize,
    place: Place,
    start: usize,
}

struct Degree {
    layer: Layer,
    classes: QuotientBasis,
}

impl Layer {
    pub(crate) fn push_block(
        &mut self,
        model: &mut Model,
        chain: Rc<[usize]>,
        dimension: usize,
        limit: usize,
    ) -> Result<()> {
        let size = self
            .size
            .checked_add(dimension)
            .filter(|&size| size <= limit)
            .ok_or_else(|| layer_limit(limit))?;
        self.cells += 4 + chain.len();
        check_cells(self.cells)?;
        model.spend(1 + chain.len())?;
        self.blocks.push(Block {
            vertical: 0,
            place: Place::Chain(chain),
            start: self.size,
        });
        self.size = size;
        Ok(())
    }

    pub(crate) fn order_blocks(
        &mut self,
        model: &mut Model,
        degree: usize,
        dimension: impl Fn(&[usize]) -> usize,
    ) -> Result<()> {
        model.spend(sort_work(self.blocks.len()).saturating_mul(degree + 1))?;
        self.blocks
            .sort_unstable_by(|left, right| left.place.cmp(&right.place));
        let mut next = 0;
        for block in &mut self.blocks {
            block.start = next;
            next += dimension(block.place.chain());
        }
        Ok(())
    }

    pub(crate) fn cellular(
        model: &mut Model,
        degree: usize,
        entries: impl Iterator<Item = (usize, usize)>,
        limit: usize,
    ) -> Result<Self> {
        let mut output = Self::default();
        for (ordinal, dimension) in entries {
            model.spend(1)?;
            if dimension == 0 {
                continue;
            }
            let start = output.size;
            output.size = start
                .checked_add(dimension)
                .filter(|&size| size <= limit)
                .ok_or_else(|| layer_limit(limit))?;
            output.cells += 6;
            check_cells(output.cells)?;
            output.blocks.push(Block {
                vertical: 0,
                place: Place::Cell { degree, ordinal },
                start,
            });
        }
        Ok(output)
    }

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
            place: block.place.clone(),
            index: index - block.start,
        }
    }

    pub(crate) fn index(&self, model: &mut Model, generator: &Generator) -> Result<usize> {
        model.spend(
            1 + generator.place.payload_cells() * (1 + self.blocks.len().max(1).ilog2() as usize),
        )?;
        let block = self
            .blocks
            .binary_search_by(|block| {
                (block.vertical, &block.place).cmp(&(generator.vertical, &generator.place))
            })
            .map_err(|_| Error::plain("sheaf cochain coordinate is outside its degree"))?;
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
        let coordinates = self.class_coordinates(model, complex, degree, raw)?;
        odd_ext::persist(model, coordinates)
    }

    pub(crate) fn class_coordinates(
        &self,
        model: &mut Model,
        complex: &(impl Cochains + ?Sized),
        degree: usize,
        raw: &[Term],
    ) -> Result<Vec<Coordinate>> {
        if raw.is_empty() {
            return Ok(Vec::new());
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
        Ok(coordinates)
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
        let term_cells = generator.cells();
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
        model.spend(1 + a.0.place.payload_cells().min(b.0.place.payload_cells()))?;
        let comparison = a.cmp(b);
        if !comparison.is_eq() {
            return Ok(comparison);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

pub(crate) fn field(model: &mut Model, characteristic: usize) -> Result<Prime> {
    match characteristic {
        2 => Ok(Prime::validated(2)),
        n if n < 2 || n % 2 == 0 => Err(Error::plain(format!(
            "sheaf characteristic {n} is not prime"
        ))),
        n => Prime::new(model, n),
    }
}

pub(crate) fn residue(prime: Prime, value: i64) -> u32 {
    value.rem_euclid(i64::from(prime.get())) as u32
}

pub(crate) fn inserted(
    model: &mut Model,
    chain: &[usize],
    slot: usize,
    point: usize,
) -> Result<Rc<[usize]>> {
    model.spend(1 + chain.len())?;
    let mut output = Vec::with_capacity(chain.len() + 1);
    output.extend_from_slice(&chain[..slot]);
    output.push(point);
    output.extend_from_slice(&chain[slot..]);
    Ok(output.into())
}

pub(crate) fn persist(model: &mut Model, terms: Vec<Term>) -> Result<Value> {
    let cells = terms.iter().map(|(generator, _)| generator.cells()).sum();
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
