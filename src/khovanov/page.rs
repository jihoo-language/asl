use std::collections::{btree_map::Entry, BTreeMap};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::{add_scaled, invert, Coordinate};

use super::{integral, persist, persist_class, Braid, Cell, Chain, Class, Context, Grade, Slice};

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct Spec {
    pub(crate) prime: Prime,
    pub(crate) r: usize,
    pub(crate) k: i64,
}

impl Spec {
    pub(crate) fn shifted(self, direction: i64) -> Result<Self> {
        let delta = direction as i128 * 2 * self.r as i128;
        let k = i64::try_from(i128::from(self.k) + delta)
            .map_err(|_| Error::plain("annular page filtration is out of range"))?;
        Ok(Self { k, ..self })
    }
}

pub(crate) fn field(model: &mut Model, p: usize) -> Result<Prime> {
    match p {
        2 => Ok(Prime::validated(2)),
        n if n < 2 || n.is_multiple_of(2) => Err(Error::plain(format!(
            "Khovanov characteristic {n} is not prime"
        ))),
        n => Prime::new(model, n),
    }
}

fn residue(prime: Prime, value: i128) -> u32 {
    value.rem_euclid(i128::from(prime.get())) as u32
}

pub(super) struct Cochains {
    ordered: Vec<(i64, Cell)>,
    by_cell: Vec<usize>,
}

impl Cochains {
    fn build(model: &mut Model, braid: &Braid, grade: Grade) -> Result<Self> {
        let mut context = Context::with_annular(braid, true);
        let basis = context.basis(model, grade, None)?;
        let width = basis.iter().map(|cell| cell.cells() + 2).sum::<usize>() + 1;
        integral::storage(width)?;
        let mut ordered = Vec::with_capacity(basis.len());
        for cell in basis {
            ordered.push((context.annular_k(model, &cell)?, cell));
        }
        model.spend(width * (1 + ordered.len().max(1).ilog2() as usize))?;
        ordered.sort_unstable();
        let mut by_cell: Vec<_> = (0..ordered.len()).collect();
        by_cell.sort_unstable_by(|&a, &b| ordered[a].1.cmp(&ordered[b].1));
        model.allocate(width)?;
        Ok(Self { ordered, by_cell })
    }

    fn index(&self, model: &mut Model, cell: &Cell) -> Result<usize> {
        model.spend(cell.cells() * (1 + self.by_cell.len().max(1).ilog2() as usize))?;
        let location = self
            .by_cell
            .binary_search_by(|&i| self.ordered[i].1.cmp(cell))
            .map_err(|_| Error::plain("annular Khovanov cell has the wrong bidegree"))?;
        Ok(self.by_cell[location])
    }

    fn coordinates(&self, model: &mut Model, chain: &Chain, k: i64) -> Result<Vec<Coordinate>> {
        let mut result = Vec::with_capacity(chain.len());
        for (cell, coefficient) in chain.iter() {
            let index = self.index(model, cell)?;
            if self.weight(index) > k {
                return Err(Error::plain(
                    "annular page lift has a term above its filtration",
                ));
            }
            result.push((index, *coefficient as u32));
        }
        model.spend(result.len() * (1 + result.len().max(1).ilog2() as usize))?;
        result.sort_unstable();
        Ok(result)
    }

    fn at(&self, k: i64) -> std::ops::Range<usize> {
        self.ordered.partition_point(|(weight, _)| *weight < k)
            ..self.ordered.partition_point(|(weight, _)| *weight <= k)
    }

    fn weight(&self, index: usize) -> i64 {
        self.ordered[index].0
    }

    fn chain(&self, model: &mut Model, vector: Vec<Coordinate>) -> Result<Chain> {
        let mut terms = Vec::with_capacity(vector.len());
        let mut width = 0;
        for (coordinate, coefficient) in vector {
            let cell = &self.ordered[coordinate].1;
            width += cell.cells();
            integral::storage(width)?;
            terms.push((cell.clone(), i128::from(coefficient)));
        }
        model.spend(width * (1 + terms.len().max(1).ilog2() as usize))?;
        terms.sort_unstable();
        persist(model, terms)
    }
}

struct Column {
    lift: Vec<Coordinate>,
    image: Vec<Coordinate>,
    leading_inverse: u32,
}

pub(super) struct Differential {
    source: Rc<Cochains>,
    target: Rc<Cochains>,
    columns: Vec<Column>,
    owners: Vec<Option<usize>>,
}

impl Differential {
    fn build(model: &mut Model, braid: &Braid, prime: Prime, grade: Grade) -> Result<Self> {
        let source = braid.page_cochains(model, grade)?;
        let target = braid.page_cochains(model, grade.shift(1)?)?;
        let mut context = Context::with_annular(braid, true);
        let mut cells = 1 + 8 * source.ordered.len() + 2 * target.ordered.len();
        integral::storage(cells)?;
        let mut columns: Vec<Column> = Vec::with_capacity(source.ordered.len());
        let mut owners: Vec<Option<usize>> = vec![None; target.ordered.len()];
        let slice = Some(Slice::AtMost(braid.ends.len() as i64));
        for (index, (_, cell)) in source.ordered.iter().enumerate() {
            let mut original = BTreeMap::new();
            if !target.ordered.is_empty() {
                context.edges(model, cell, false, slice, |model, term, sign| {
                    let row = target.index(model, &term)?;
                    insert(prime, &mut original, row, residue(prime, sign));
                    integral::storage(6 * original.len())
                })?;
            }
            let mut image: Vec<_> = original.into_iter().collect();
            let mut lift = vec![(index, 1)];
            let mut leading_inverse = 1;
            // Earlier columns have only cells of smaller filtration, or
            // earlier cells at equal filtration. Thus each tracked lift
            // retains a distinct leading source cell and each nonzero image
            // has a distinct leading target cell.
            model.spend(1 + image.len())?;
            while let Some(&(pivot, coefficient)) = image.last() {
                if let Some(owner) = owners[pivot] {
                    let known = &columns[owner];
                    let scalar = prime.negate(coefficient);
                    image = add_scaled(model, prime, &image, &known.image, scalar)?;
                    lift = add_scaled(model, prime, &lift, &known.lift, scalar)?;
                    continue;
                }
                owners[pivot] = Some(index);
                leading_inverse = coefficient;
                if coefficient != 1 {
                    let scalar = invert(model, prime, coefficient)?;
                    model.spend(image.len() + lift.len())?;
                    for (_, coefficient) in image.iter_mut().chain(&mut lift) {
                        *coefficient = prime.mul(*coefficient, scalar);
                    }
                }
                break;
            }
            cells += 2 * (image.len() + lift.len());
            integral::storage(cells)?;
            columns.push(Column {
                lift,
                image,
                leading_inverse,
            });
        }
        model.allocate(cells)?;
        Ok(Self {
            source,
            target,
            columns,
            owners,
        })
    }

    fn delay(&self, source: usize, target: usize) -> Result<usize> {
        let drop = self.source.weight(source) - self.target.weight(target);
        if drop < 0 || drop % 2 != 0 {
            return Err(Error::plain(
                "Khovanov differential violates annular filtration",
            ));
        }
        Ok((drop / 2) as usize)
    }
}

fn insert<K: Ord>(prime: Prime, output: &mut BTreeMap<K, u32>, key: K, value: u32) {
    if value == 0 {
        return;
    }
    match output.entry(key) {
        Entry::Vacant(entry) => {
            entry.insert(value);
        }
        Entry::Occupied(mut entry) => {
            let value = prime.add(*entry.get(), value);
            if value == 0 {
                entry.remove();
            } else {
                *entry.get_mut() = value;
            }
        }
    }
}

#[derive(Clone, Copy)]
enum Partner {
    Earlier { index: usize, delay: usize },
    Later { index: usize, delay: usize },
    Permanent,
}

struct Page {
    spec: Spec,
    outgoing: Rc<Differential>,
    incoming: Rc<Differential>,
    generators: Vec<usize>,
}

impl Page {
    fn build(model: &mut Model, braid: &Braid, grade: Grade, spec: Spec) -> Result<Self> {
        let outgoing = braid.page_differential(model, spec.prime, grade)?;
        let incoming = braid.page_differential(model, spec.prime, grade.shift(-1)?)?;
        let mut result = Self {
            spec,
            outgoing,
            incoming,
            generators: Vec::new(),
        };
        for index in result.outgoing.source.at(spec.k) {
            model.spend(1)?;
            let live = match result.partner(index)? {
                Partner::Earlier { delay, .. } | Partner::Later { delay, .. } => delay >= spec.r,
                Partner::Permanent => true,
            };
            if live {
                result.generators.push(index);
            }
        }
        integral::storage(result.generators.len())?;
        Ok(result)
    }

    fn partner(&self, index: usize) -> Result<Partner> {
        if let Some(&(target, _)) = self.outgoing.columns[index].image.last() {
            if self.incoming.owners[index].is_some() {
                return Err(Error::plain(
                    "the signed filtered differential does not square to zero",
                ));
            }
            return Ok(Partner::Later {
                index: target,
                delay: self.outgoing.delay(index, target)?,
            });
        }
        if let Some(source) = self.incoming.owners[index] {
            return Ok(Partner::Earlier {
                index: source,
                delay: self.incoming.delay(source, index)?,
            });
        }
        Ok(Partner::Permanent)
    }

    fn generator(&self, index: usize) -> Result<usize> {
        self.generators
            .get(index)
            .copied()
            .ok_or_else(|| Error::plain("annular page generator index is out of range"))
    }

    fn position(&self, model: &mut Model, index: usize) -> Result<usize> {
        model.spend(1 + self.generators.len().max(1).ilog2() as usize)?;
        self.generators
            .binary_search(&index)
            .map_err(|_| Error::plain("incompatible filtered Khovanov normal forms"))
    }

    fn adapted(&self, index: usize) -> (&[Coordinate], u32) {
        // A boundary pivot is led by a cell that led an outgoing kernel
        // vector. Replace that vector with the paired image. The resulting
        // triangular basis splits the filtered complex into actual pairs
        // db=z and permanent cycles.
        if let Some(source) = self.incoming.owners[index] {
            (&self.incoming.columns[source].image, 1)
        } else {
            let column = &self.outgoing.columns[index];
            (&column.lift, column.leading_inverse)
        }
    }

    pub(crate) fn rank(&self) -> usize {
        self.generators.len()
    }

    pub(crate) fn class(&self, model: &mut Model, chain: &Chain) -> Result<Class> {
        let mut vector = self
            .outgoing
            .source
            .coordinates(model, chain, self.spec.k)?;
        let mut class = Vec::new();
        let limit = i128::from(self.spec.k) - 2 * self.spec.r as i128;
        // Test the full lifted cochain, including all terms below k. The
        // nonzero images have distinct pivots, so the greatest filtration
        // among their pivots survives in their sum.
        while let Some(&(index, leading)) = vector.last() {
            let (adapted, inverse) = self.adapted(index);
            let coefficient = self.spec.prime.mul(leading, inverse);
            let partner = self.partner(index)?;
            if let Partner::Later { index: target, .. } = partner {
                if i128::from(self.outgoing.target.weight(target)) > limit {
                    return Err(Error::plain(
                        "annular page lift has a differential above its required cutoff",
                    ));
                }
            }
            if self.outgoing.source.weight(index) == self.spec.k {
                let killed =
                    matches!(partner, Partner::Earlier { delay, .. } if delay < self.spec.r);
                if !killed {
                    class.push((self.position(model, index)?, i128::from(coefficient)));
                }
            }
            vector = add_scaled(
                model,
                self.spec.prime,
                &vector,
                adapted,
                self.spec.prime.negate(coefficient),
            )?;
        }
        persist_class(model, class.into_iter().collect())
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &Class) -> Result<Chain> {
        let mut vector = Vec::new();
        for &(generator, coefficient) in class.iter() {
            let (adapted, _) = self.adapted(self.generator(generator)?);
            vector = add_scaled(model, self.spec.prime, &vector, adapted, coefficient as u32)?;
        }
        self.outgoing.source.chain(model, vector)
    }

    pub(crate) fn map(
        &self,
        model: &mut Model,
        class: &Class,
        operation: Map,
        target: impl FnOnce(&mut Model) -> Result<Self>,
    ) -> Result<Class> {
        let mut images = Vec::new();
        model.spend(class.len())?;
        for &(g, coefficient) in class.iter() {
            let source = self.generator(g)?;
            let partner = self.partner(source)?;
            let image = match (operation, partner) {
                (Map::Differential, Partner::Later { index, delay }) if delay == self.spec.r => {
                    Some(index)
                }
                (Map::Differential, _) => None,
                (Map::Boundary, Partner::Earlier { index, delay }) if delay == self.spec.r => {
                    Some(index)
                }
                (Map::Boundary, _) => {
                    return Err(Error::plain(format!(
                        "annular page class is not a d_{} boundary",
                        self.spec.r
                    )));
                }
                (Map::Next, Partner::Later { delay, .. }) if delay == self.spec.r => {
                    return Err(Error::plain(format!(
                        "annular page class has nonzero d_{}",
                        self.spec.r
                    )));
                }
                (Map::Next, Partner::Earlier { delay, .. }) if delay == self.spec.r => None,
                (Map::Next | Map::Previous, _) => Some(source),
            };
            if let Some(index) = image {
                images.push((index, coefficient));
            }
        }
        if images.is_empty() {
            return Ok(Rc::default());
        }
        let target = target(model)?;
        let mut output = BTreeMap::new();
        for (index, coefficient) in images {
            output.insert(target.position(model, index)?, coefficient);
        }
        persist_class(model, output)
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Map {
    Differential,
    Boundary,
    Next,
    Previous,
}

impl Braid {
    fn page_cochains(&self, model: &mut Model, grade: Grade) -> Result<Rc<Cochains>> {
        if let Some(cached) = self.page_cochains.borrow().get(&grade) {
            return Ok(cached.clone());
        }
        let cochains = Rc::new(Cochains::build(model, self, grade)?);
        self.page_cochains
            .borrow_mut()
            .insert(grade, cochains.clone());
        Ok(cochains)
    }

    fn page_differential(
        &self,
        model: &mut Model,
        prime: Prime,
        grade: Grade,
    ) -> Result<Rc<Differential>> {
        if let Some(cached) = self.page_differentials.borrow().get(&(prime, grade)) {
            return Ok(cached.clone());
        }
        let differential = Rc::new(Differential::build(model, self, prime, grade)?);
        self.page_differentials
            .borrow_mut()
            .insert((prime, grade), differential.clone());
        Ok(differential)
    }

    fn page(&self, model: &mut Model, grade: Grade, spec: Spec) -> Result<Page> {
        Page::build(model, self, grade, spec)
    }

    pub(crate) fn page_rank(&self, model: &mut Model, grade: Grade, spec: Spec) -> Result<usize> {
        if !self.possible(grade, super::Form::Page(spec)) {
            return Ok(0);
        }
        Ok(self.page(model, grade, spec)?.rank())
    }

    pub(crate) fn page_generator(
        &self,
        model: &mut Model,
        grade: Grade,
        spec: Spec,
        index: usize,
    ) -> Result<Class> {
        if !self.possible(grade, super::Form::Page(spec)) {
            return Err(Error::plain("annular page generator index is out of range"));
        }
        self.page(model, grade, spec)?.generator(index)?;
        persist_class(model, [(index, 1)].into())
    }

    pub(crate) fn modular_cutoff(
        &self,
        model: &mut Model,
        prime: Prime,
        input: &Chain,
        k: i64,
    ) -> Result<Chain> {
        let mut context = Context::with_annular(self, true);
        let mut output = Vec::new();
        for (cell, value) in input.iter() {
            let coefficient = residue(prime, *value);
            if coefficient != 0 && context.annular_k(model, cell)? <= k {
                output.push((cell.clone(), i128::from(coefficient)));
            }
        }
        persist(model, output)
    }

    pub(crate) fn modular_differential(
        &self,
        model: &mut Model,
        prime: Prime,
        input: &Chain,
    ) -> Result<Chain> {
        let mut context = Context::with_annular(self, true);
        let mut sum = ChainSum::default();
        let slice = Some(Slice::AtMost(self.ends.len() as i64));
        for (cell, coefficient) in input.iter() {
            context.edges(model, cell, false, slice, |model, term, sign| {
                let scalar = prime.mul(*coefficient as u32, residue(prime, sign));
                sum.add(model, prime, term, scalar)
            })?;
        }
        sum.finish(model)
    }

    pub(crate) fn modular_elements(
        &self,
        model: &mut Model,
        grade: Grade,
        prime: Prime,
        k: i64,
        limit: usize,
    ) -> Result<Vec<Value>> {
        if !self.possible(grade, super::Form::Modular(prime, k)) {
            return Ok(vec![Value::KhChain(Rc::default())]);
        }
        let cochains = self.page_cochains(model, grade)?;
        model.spend(cochains.by_cell.len())?;
        let selected: Vec<_> = cochains
            .by_cell
            .iter()
            .copied()
            .filter(|&i| cochains.weight(i) <= k)
            .collect();
        elements(model, prime, selected.len(), limit, |model, coordinates| {
            let vector = coordinates
                .into_iter()
                .map(|(i, c)| (selected[i], c))
                .collect();
            cochains.chain(model, vector).map(Value::KhChain)
        })
    }

    pub(crate) fn page_elements(
        &self,
        model: &mut Model,
        grade: Grade,
        spec: Spec,
        limit: usize,
    ) -> Result<Vec<Value>> {
        let rank = self.page_rank(model, grade, spec)?;
        elements(model, spec.prime, rank, limit, |model, coordinates| {
            persist_field_class(model, coordinates).map(Value::KhClass)
        })
    }
}

fn elements(
    model: &mut Model,
    prime: Prime,
    dimension: usize,
    limit: usize,
    mut make: impl FnMut(&mut Model, Vec<Coordinate>) -> Result<Value>,
) -> Result<Vec<Value>> {
    let field = prime.get() as usize;
    let size = field
        .checked_pow(dimension as u32)
        .filter(|&n| n <= limit)
        .ok_or_else(|| Error::plain("modular Khovanov space exceeds the finite element limit"))?;
    model.allocate(size)?;
    model.spend(size * (1 + dimension))?;
    let mut codes = (0..size)
        .map(|mut code| {
            let mut coordinates = Vec::new();
            for i in 0..dimension {
                let digit = code % field;
                code /= field;
                if digit != 0 {
                    coordinates.push((i, digit as u32));
                }
            }
            coordinates
        })
        .collect::<Vec<_>>();
    // The interpreter expects the same lexicographic order as value comparison.
    model.spend(size * (1 + dimension) * (1 + size.max(1).ilog2() as usize))?;
    codes.sort_unstable();
    codes.into_iter().map(|code| make(model, code)).collect()
}

fn persist_field_class(model: &mut Model, coordinates: Vec<Coordinate>) -> Result<Class> {
    persist_class(
        model,
        coordinates
            .into_iter()
            .map(|(i, value)| (i, i128::from(value)))
            .collect(),
    )
}

#[derive(Default)]
struct ChainSum {
    terms: BTreeMap<Cell, u32>,
    cells: usize,
}

impl ChainSum {
    fn add(&mut self, model: &mut Model, prime: Prime, cell: Cell, value: u32) -> Result<()> {
        let width = cell.cells();
        let before = self.terms.len();
        model.spend(width * (1 + before.max(1).ilog2() as usize))?;
        insert(prime, &mut self.terms, cell, value);
        match self.terms.len().cmp(&before) {
            std::cmp::Ordering::Greater => {
                self.cells += width;
                integral::storage(self.cells)?;
            }
            std::cmp::Ordering::Less => self.cells -= width,
            std::cmp::Ordering::Equal => {}
        }
        Ok(())
    }

    fn finish(self, model: &mut Model) -> Result<Chain> {
        persist(
            model,
            self.terms
                .into_iter()
                .map(|(cell, value)| (cell, i128::from(value)))
                .collect(),
        )
    }
}

pub(crate) fn add_chain(model: &mut Model, prime: Prime, a: &Chain, b: &Chain) -> Result<Chain> {
    let mut sum = ChainSum::default();
    for (cell, coefficient) in a.iter().chain(b.iter()) {
        sum.add(model, prime, cell.clone(), *coefficient as u32)?;
    }
    sum.finish(model)
}

pub(crate) fn scale_chain(
    model: &mut Model,
    prime: Prime,
    a: &Chain,
    factor: i128,
) -> Result<Chain> {
    let scalar = residue(prime, factor);
    if scalar == 0 {
        return Ok(Rc::default());
    }
    model.spend(a.len())?;
    persist(
        model,
        a.iter()
            .map(|(cell, value)| (cell.clone(), i128::from(prime.mul(*value as u32, scalar))))
            .collect(),
    )
}

pub(crate) fn add_class(model: &mut Model, prime: Prime, a: &Class, b: &Class) -> Result<Class> {
    let a: Vec<_> = a.iter().map(|&(g, value)| (g, value as u32)).collect();
    let b: Vec<_> = b.iter().map(|&(g, value)| (g, value as u32)).collect();
    let sum = add_scaled(model, prime, &a, &b, 1)?;
    persist_field_class(model, sum)
}

pub(crate) fn scale_class(
    model: &mut Model,
    prime: Prime,
    a: &Class,
    factor: i128,
) -> Result<Class> {
    let scalar = residue(prime, factor);
    if scalar == 0 {
        return Ok(Rc::default());
    }
    model.spend(a.len())?;
    persist_field_class(
        model,
        a.iter()
            .map(|&(g, value)| (g, prime.mul(value as u32, scalar)))
            .collect(),
    )
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Reduce(usize, i64),
    Class(usize, i64),
    Representative,
    Map(Map),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Reduce(..) => "khMod",
            Self::Class(..) => "khPageClass",
            Self::Representative => "khPageRep",
            Self::Map(Map::Differential) => "khPageD",
            Self::Map(Map::Boundary) => "khPageBoundary",
            Self::Map(Map::Next) => "khPageNext",
            Self::Map(Map::Previous) => "khPagePrev",
        }
    }

    pub(crate) fn target(
        self,
        model: &mut Model,
        grade: Grade,
        form: super::Form,
    ) -> Result<(Grade, super::Form)> {
        use super::Form;
        match (self, form) {
            (
                Self::Reduce(p, k),
                Form::Chain | Form::AnnularChain(_) | Form::Filtered(_) | Form::Modular(..),
            ) => {
                let prime = field(model, p)?;
                let source = form
                    .modular()
                    .map(|(prime, _)| prime)
                    .or_else(|| matches!(form, Form::Filtered(_)).then(|| Prime::validated(2)));
                if source.is_some_and(|source| source != prime) {
                    return Err(Error::plain(
                        "khMod requires the same characteristic on modular cochains",
                    ));
                }
                Ok((grade, Form::Modular(prime, k)))
            }
            (Self::Class(r, k), Form::Modular(prime, _)) => {
                Ok((grade, Form::Page(Spec { prime, r, k })))
            }
            (Self::Representative, Form::Page(spec)) => {
                Ok((grade, Form::Modular(spec.prime, spec.k)))
            }
            (Self::Map(map), Form::Page(spec)) => {
                let (grade, spec) = match map {
                    Map::Differential => (grade.shift(1)?, spec.shifted(-1)?),
                    Map::Boundary => (grade.shift(-1)?, spec.shifted(1)?),
                    Map::Next | Map::Previous => {
                        let r = (match map {
                            Map::Next => spec.r.checked_add(1),
                            _ => spec.r.checked_sub(1),
                        })
                        .ok_or_else(|| Error::plain("annular page index is out of range"))?;
                        (grade, Spec { r, ..spec })
                    }
                };
                Ok((grade, Form::Page(spec)))
            }
            _ => {
                let expected = match self {
                    Self::Reduce(..) => "an integral or modular Khovanov cochain",
                    Self::Class(..) => "a KhPC filtered cochain",
                    _ => "a KhPage class",
                };
                Err(Error::plain(format!("{} requires {expected}", self.name())))
            }
        }
    }

    pub(crate) fn apply(
        self,
        model: &mut Model,
        braid: &Braid,
        input: (Grade, super::Form, Value),
        output: (Grade, super::Form),
    ) -> Result<Value> {
        use super::Form;
        let (grade, form, value) = input;
        let (target_grade, target_form) = output;
        match (self, form, value, target_form) {
            (Self::Reduce(..), _, Value::KhChain(chain), Form::Modular(prime, k)) => Ok(
                Value::KhChain(braid.modular_cutoff(model, prime, &chain, k)?),
            ),
            (Self::Class(..), _, Value::KhChain(chain), Form::Page(spec)) => {
                if chain.is_empty() {
                    return Ok(Value::KhClass(Rc::default()));
                }
                Ok(Value::KhClass(
                    braid.page(model, grade, spec)?.class(model, &chain)?,
                ))
            }
            (Self::Representative, Form::Page(spec), Value::KhClass(class), _) => {
                if class.is_empty() {
                    return Ok(Value::KhChain(Rc::default()));
                }
                Ok(Value::KhChain(
                    braid
                        .page(model, grade, spec)?
                        .representative(model, &class)?,
                ))
            }
            (Self::Map(map), Form::Page(spec), Value::KhClass(class), Form::Page(target)) => {
                if class.is_empty() {
                    return Ok(Value::KhClass(Rc::default()));
                }
                let source = braid.page(model, grade, spec)?;
                let output = source.map(model, &class, map, |model| {
                    braid.page(model, target_grade, target)
                })?;
                Ok(Value::KhClass(output))
            }
            _ => unreachable!(),
        }
    }
}
