use std::cmp::Ordering;
use std::collections::{hash_map::Entry, BTreeSet, HashMap};
use std::ops::RangeInclusive;
use std::rc::Rc;

use crate::cobar::{self, Word};
use crate::cohomology;
use crate::comodule::Generator;
use crate::complex::xor_sorted;
use crate::error::{Error, Result};
use crate::linear::{xor, Basis, QuotientBasis};
use crate::map_ext;
use crate::model::{Kind, Model, TypeId, Value};
use crate::profile::Profile;
use crate::steenrod;
use crate::v_koszul;

const MAX_BASIS_TERMS: usize = 32_768;
const MAX_CELLS: usize = 300_000;

pub(crate) type Indices = Rc<[usize]>;
pub(crate) type Support = Rc<[Term]>;
pub(crate) type Powers = Rc<[i128]>;

#[derive(Clone, PartialEq, Eq, Hash)]
pub(crate) struct Grade<S = usize> {
    pub(crate) source: TypeId,
    pub(crate) target: TypeId,
    pub(crate) selected: Indices,
    pub(crate) s: S,
    pub(crate) t: i64,
}

#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Term {
    pub(crate) source: Generator,
    pub(crate) target: Generator,
    pub(crate) powers: Powers,
}

impl Term {
    pub(crate) fn cells(&self) -> usize {
        5 + power_cells(self.powers.len())
    }
}

pub(crate) fn power_cells(width: usize) -> usize {
    width * std::mem::size_of::<i128>().div_ceil(std::mem::size_of::<usize>())
}

pub(crate) fn indices(model: &mut Model, indices: Vec<usize>) -> Result<Indices> {
    model.spend(1 + indices.len())?;
    if indices.windows(2).any(|pair| pair[0] >= pair[1]) {
        return Err(Error::plain(
            "Milnor primitive indices must be strictly increasing",
        ));
    }
    if indices
        .last()
        .is_some_and(|&index| index >= usize::BITS as usize)
    {
        return Err(Error::plain("Milnor primitive degree is out of range"));
    }
    model.allocate(1 + indices.len())?;
    Ok(indices.into())
}

pub(crate) fn weight(index: usize) -> usize {
    usize::MAX >> (usize::BITS as usize - 1 - index)
}

pub(crate) fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain("Koszul calculation exceeds the storage limit"));
    }
    Ok(())
}

fn basis_limit() -> Error {
    Error::plain("Koszul basis exceeds the enumeration limit")
}

#[derive(Default)]
pub(crate) struct QExt {
    terms: Vec<Term>,
    quotient: QuotientBasis,
}

impl QExt {
    pub(crate) fn build(
        model: &mut Model,
        source: TypeId,
        target: TypeId,
        selected: &[usize],
        s: usize,
        t: i64,
    ) -> Result<Self> {
        Self::build_in(
            model,
            source,
            target,
            selected,
            Filtration::Polynomial(s),
            t,
        )
    }

    pub(crate) fn build_localized(
        model: &mut Model,
        source: TypeId,
        target: TypeId,
        selected: &[usize],
        s: i64,
        t: i64,
    ) -> Result<Self> {
        v_koszul::top(selected)?;
        Self::build_in(
            model,
            source,
            target,
            selected,
            Filtration::Localized(i128::from(s)),
            t,
        )
    }

    fn build_in(
        model: &mut Model,
        source: TypeId,
        target: TypeId,
        selected: &[usize],
        filtration: Filtration,
        t: i64,
    ) -> Result<Self> {
        let mut complex = Complex::new(model, source, target, selected)?;
        let terms = complex.basis(model, filtration, t, MAX_BASIS_TERMS)?;
        if terms.is_empty() {
            return Ok(Self::default());
        }

        let mut boundaries = Basis::default();
        if let Some(previous) = filtration.previous() {
            for term in complex.basis(model, previous, t, MAX_BASIS_TERMS)? {
                let image = complex.differential(model, std::slice::from_ref(&term))?;
                let column = term_indices(model, &terms, &image)?;
                boundaries.insert(model, column)?;
            }
        }

        let mut rows = HashMap::new();
        let mut row_cells = 0;
        let quotient = QuotientBasis::build(model, boundaries, terms.len(), |model, index| {
            let mut column = Vec::new();
            for term in complex.differential(model, std::slice::from_ref(&terms[index]))? {
                let next = rows.len();
                let cells = term.cells();
                model.spend(cells)?;
                let coordinate = match rows.entry(term) {
                    Entry::Occupied(entry) => *entry.get(),
                    Entry::Vacant(entry) => {
                        row_cells += cells;
                        check_storage(row_cells)?;
                        *entry.insert(next)
                    }
                };
                column.push(coordinate);
            }
            model.spend(column.len() * (1 + column.len().max(1).ilog2() as usize))?;
            column.sort_unstable();
            Ok(column)
        })?;
        model.allocate(terms.iter().map(Term::cells).sum())?;
        Ok(Self { terms, quotient })
    }

    pub(crate) fn rank(&self) -> usize {
        self.quotient.rank()
    }

    pub(crate) fn project(&self, model: &mut Model, terms: &[Term]) -> Result<Vec<usize>> {
        let coordinates = term_indices(model, &self.terms, terms)?;
        self.quotient
            .project(model, coordinates)?
            .ok_or_else(|| Error::plain("qClass requires a closed Koszul cochain"))
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<Term>> {
        let indices = self.quotient.representative(model, class)?;
        Ok(indices
            .into_iter()
            .map(|index| self.terms[index].clone())
            .collect())
    }
}

fn term_indices(model: &mut Model, basis: &[Term], terms: &[Term]) -> Result<Vec<usize>> {
    let search = 1 + basis.len().max(1).ilog2() as usize;
    terms
        .iter()
        .map(|term| {
            model.spend(search * term.cells())?;
            basis
                .binary_search(term)
                .map_err(|_| Error::plain("term is outside its Koszul bidegree"))
        })
        .collect()
}

pub(crate) struct Complex {
    source: TypeId,
    target: TypeId,
    source_bound: Option<usize>,
    target_bound: Option<usize>,
    primitives: Vec<usize>,
    weights: Vec<usize>,
    actions: HashMap<(TypeId, usize, usize), Action>,
    action_cells: usize,
}

struct Action {
    columns: Vec<Rc<[usize]>>,
    rows: Vec<Vec<usize>>,
}

#[derive(Clone, Copy)]
enum Filtration {
    Polynomial(usize),
    Localized(i128),
}

impl Filtration {
    fn previous(self) -> Option<Self> {
        match self {
            Self::Polynomial(s) => s.checked_sub(1).map(Self::Polynomial),
            Self::Localized(s) => s.checked_sub(1).map(Self::Localized),
        }
    }

    fn word_bounds(self, weights: &[usize]) -> (i128, i128) {
        let last = weights.last().copied().unwrap_or(0);
        match self {
            Self::Polynomial(s) => (
                wide_degree(s, weights.first().copied().unwrap_or(0)),
                wide_degree(s, last),
            ),
            Self::Localized(s) => {
                let upper = s.saturating_mul(last as i128);
                let lower = if weights.len() == 1 { upper } else { i128::MIN };
                (lower, upper)
            }
        }
    }

    fn monomials(
        self,
        model: &mut Model,
        weights: &[usize],
        degree: i128,
        limit: usize,
    ) -> Result<Vec<Powers>> {
        match self {
            Self::Polynomial(s) => polynomial_basis(model, weights, s, degree as u128, limit),
            Self::Localized(s) => v_koszul::laurent_basis(model, weights, s, degree, limit),
        }
    }
}

impl Complex {
    pub(crate) fn new(
        model: &mut Model,
        source: TypeId,
        target: TypeId,
        selected: &[usize],
    ) -> Result<Self> {
        Ok(Self {
            source,
            target,
            source_bound: endpoint_bound(model, source)?,
            target_bound: endpoint_bound(model, target)?,
            primitives: selected.to_vec(),
            weights: selected.iter().map(|&index| weight(index)).collect(),
            actions: HashMap::new(),
            action_cells: 0,
        })
    }

    pub(crate) fn coefficient_basis(&self, model: &mut Model) -> Result<Vec<Term>> {
        if self.source_bound.is_none() || self.target_bound.is_none() {
            return Ok(Vec::new());
        }
        let sources = margolis_groups(model, self.source)?;
        if sources.is_empty() {
            return Ok(Vec::new());
        }
        let targets = margolis_groups(model, self.target)?;
        let total = |ranks: &[(usize, usize)]| ranks.iter().map(|&(_, rank)| rank).sum::<usize>();
        let size = total(&sources).saturating_mul(total(&targets));
        if size > MAX_BASIS_TERMS {
            return Err(basis_limit());
        }
        let cells = size * (5 + power_cells(self.primitives.len()));
        check_storage(cells)?;
        model.spend(1 + cells)?;
        let zero: Powers = vec![0; self.primitives.len()].into();
        let mut basis = Vec::with_capacity(size);
        for (q, source_rank) in sources {
            for index in 0..source_rank {
                let source = Generator { degree: q, index };
                for &(p, target_rank) in &targets {
                    for index in 0..target_rank {
                        let target = Generator { degree: p, index };
                        basis.push(Term {
                            source,
                            target,
                            powers: zero.clone(),
                        });
                    }
                }
            }
        }
        Ok(basis)
    }

    fn basis(
        &self,
        model: &mut Model,
        filtration: Filtration,
        t: i64,
        limit: usize,
    ) -> Result<Vec<Term>> {
        let (Some(source_bound), Some(target_bound)) = (self.source_bound, self.target_bound)
        else {
            return Ok(Vec::new());
        };
        if self.weights.is_empty() && matches!(filtration, Filtration::Polynomial(s) if s != 0) {
            return Ok(Vec::new());
        }
        let (min_word, max_word) = filtration.word_bounds(&self.weights);
        let t = i128::from(t);
        let first_source = min_word.saturating_sub(t).max(0);
        let max_source = max_word
            .saturating_add(target_bound as i128)
            .saturating_sub(t);
        if first_source > (source_bound as i128).min(max_source) {
            return Ok(Vec::new());
        }
        let source_bound = exact_bound(model, self.source, source_bound)?;
        let target_bound = exact_bound(model, self.target, target_bound)?;
        let last_source = (source_bound as i128).min(
            max_word
                .saturating_add(target_bound as i128)
                .saturating_sub(t),
        );
        if first_source > last_source {
            return Ok(Vec::new());
        }
        let first_target = (first_source + t).saturating_sub(max_word).max(0);
        let last_target = (last_source + t)
            .saturating_sub(min_word)
            .min(target_bound as i128);
        if first_target > last_target {
            return Ok(Vec::new());
        }
        let sources = graded_ranks(
            model,
            self.source,
            first_source as usize..=last_source as usize,
        )?;
        if sources.is_empty() {
            return Ok(Vec::new());
        }
        let targets = graded_ranks(
            model,
            self.target,
            first_target as usize..=last_target as usize,
        )?;
        let mut catalogue = HashMap::<i128, Vec<Powers>>::new();
        let mut catalogue_cells = 0;
        let mut terms = Vec::new();
        let mut cells = 0;
        for (q, source_rank) in sources {
            let total = q as i128 + t;
            let min_p = total.saturating_sub(max_word);
            let max_p = total.saturating_sub(min_word);
            let first = targets.partition_point(|&(degree, _)| (degree as i128) < min_p);
            for &(p, target_rank) in &targets[first..] {
                model.spend(1)?;
                if p as i128 > max_p {
                    break;
                }
                let degree = total - p as i128;
                let copies = source_rank.saturating_mul(target_rank);
                let available = (limit - terms.len()) / copies;
                let powers = match catalogue.entry(degree) {
                    Entry::Occupied(entry) => entry.into_mut(),
                    Entry::Vacant(entry) => {
                        let powers =
                            filtration.monomials(model, &self.weights, degree, available)?;
                        catalogue_cells += 1 + powers.len() * (1 + power_cells(self.weights.len()));
                        check_storage(catalogue_cells)?;
                        entry.insert(powers)
                    }
                };
                if powers.len() > available {
                    return Err(basis_limit());
                }
                if powers.is_empty() {
                    continue;
                }
                let added = copies * powers.len() * (5 + power_cells(self.weights.len()));
                cells += added;
                check_storage(cells)?;
                model.spend(added)?;
                for source in 0..source_rank {
                    for target in 0..target_rank {
                        for powers in powers.iter() {
                            terms.push(Term {
                                source: Generator {
                                    degree: q,
                                    index: source,
                                },
                                target: Generator {
                                    degree: p,
                                    index: target,
                                },
                                powers: powers.clone(),
                            });
                        }
                    }
                }
            }
        }
        model.spend(cells * (1 + terms.len().max(1).ilog2() as usize))?;
        terms.sort_unstable();
        Ok(terms)
    }

    pub(crate) fn differential(&mut self, model: &mut Model, terms: &[Term]) -> Result<Vec<Term>> {
        let mut support = Accumulator::default();
        let (Some(source_bound), Some(target_bound)) = (self.source_bound, self.target_bound)
        else {
            return Ok(Vec::new());
        };
        if terms.is_empty() || self.primitives.is_empty() {
            return Ok(Vec::new());
        }
        let source_bound = exact_bound(model, self.source, source_bound)?;
        let target_bound = exact_bound(model, self.target, target_bound)?;
        let source_is_module = model.q_module(self.source).is_some();
        let target_is_module = model.q_module(self.target).is_some();
        for term in terms {
            model.spend(term.cells())?;
            for variable in 0..self.primitives.len() {
                model.spend(1)?;
                let index = self.primitives[variable];
                let weight = self.weights[variable];
                let source_degree = term.source.degree;
                if let Some(degree) = source_degree
                    .checked_add(weight)
                    .filter(|&d| d <= source_bound)
                {
                    if source_degree != 0 || source_is_module {
                        let action =
                            self.action(model, self.source, index, source_degree, degree)?;
                        let affected = &action.columns[term.source.index];
                        if !affected.is_empty() {
                            let powers = increment(model, &term.powers, variable)?;
                            for &index in affected.iter() {
                                support.toggle(
                                    model,
                                    Term {
                                        source: Generator { degree, index },
                                        target: term.target,
                                        powers: powers.clone(),
                                    },
                                )?;
                            }
                        }
                    }
                }
                let target_degree = term.target.degree;
                if let Some(degree) = target_degree
                    .checked_sub(weight)
                    .filter(|&d| (d != 0 || target_is_module) && d <= target_bound)
                {
                    let action = self.action(model, self.target, index, degree, target_degree)?;
                    let affected = &action.rows[term.target.index];
                    if !affected.is_empty() {
                        let powers = increment(model, &term.powers, variable)?;
                        for &index in affected {
                            support.toggle(
                                model,
                                Term {
                                    source: term.source,
                                    target: Generator { degree, index },
                                    powers: powers.clone(),
                                },
                            )?;
                        }
                    }
                }
            }
        }
        Ok(support.finish())
    }

    pub(crate) fn margolis_image(
        &mut self,
        model: &mut Model,
        base: TypeId,
        primitive: usize,
        lower: usize,
        generator: Generator,
    ) -> Result<&[usize]> {
        if lower == 0 && model.q_module(base).is_none() {
            return Ok(&[]);
        }
        let action = self.action(model, base, primitive, lower, generator.degree)?;
        Ok(&action.rows[generator.index])
    }

    fn action(
        &mut self,
        model: &mut Model,
        base: TypeId,
        primitive: usize,
        lower: usize,
        upper: usize,
    ) -> Result<&Action> {
        Ok(match self.actions.entry((base, primitive, lower)) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => {
                let columns = rank(model, base, lower)?;
                let rows = rank(model, base, upper)?;
                self.action_cells += 1 + rows + columns;
                check_storage(self.action_cells)?;
                model.spend(1 + rows + columns)?;
                let mut action = Action {
                    columns: Vec::with_capacity(columns),
                    rows: vec![Vec::new(); rows],
                };
                if let Some(module) = model.q_module(base) {
                    let mut transpose = vec![Vec::new(); columns];
                    for (index, row) in action.rows.iter_mut().enumerate() {
                        let image = module.image(
                            primitive,
                            Generator {
                                degree: upper,
                                index,
                            },
                        );
                        self.action_cells += 2 * image.len();
                        check_storage(self.action_cells)?;
                        model.spend(1 + image.len())?;
                        row.extend_from_slice(image);
                        for &output in image {
                            transpose[output].push(index);
                        }
                    }
                    action.columns = transpose.into_iter().map(Rc::from).collect();
                    return Ok(entry.insert(action));
                }
                for index in 0..columns {
                    let image = if rows == 0 {
                        [].into()
                    } else {
                        primitive_action(model, base, primitive, lower, [index].into())?
                    };
                    for &output in image.iter() {
                        action.rows[output].push(index);
                    }
                    self.action_cells += 2 * image.len();
                    check_storage(self.action_cells)?;
                    action.columns.push(image);
                }
                entry.insert(action)
            }
        })
    }
}

fn rank(model: &mut Model, base: TypeId, degree: usize) -> Result<usize> {
    if let Some(module) = model.q_module(base) {
        return Ok(module.rank(degree));
    }
    let ty = model.cohomology_type(base, degree)?;
    Ok(model.cohomology(ty)?.rank())
}

fn graded_ranks(
    model: &mut Model,
    base: TypeId,
    range: RangeInclusive<usize>,
) -> Result<Vec<(usize, usize)>> {
    if let Some(module) = model.q_module(base) {
        let ranks = module.graded_ranks(range);
        model.spend(1 + ranks.len())?;
        return Ok(ranks);
    }
    let mut result = Vec::new();
    for degree in range {
        model.spend(1)?;
        let rank = rank(model, base, degree)?;
        if rank != 0 {
            result.push((degree, rank));
        }
    }
    Ok(result)
}

pub(crate) fn margolis_groups(model: &mut Model, base: TypeId) -> Result<Vec<(usize, usize)>> {
    let Some(bound) = endpoint_bound(model, base)? else {
        return Ok(Vec::new());
    };
    let bound = exact_bound(model, base, bound)?;
    graded_ranks(model, base, 0..=bound)
}

fn endpoint_bound(model: &mut Model, base: TypeId) -> Result<Option<usize>> {
    if let Some(module) = model.q_module(base) {
        return Ok(module.bound());
    }
    Ok(model.height(base, false)?.checked_sub(1))
}

fn exact_bound(model: &mut Model, base: TypeId, bound: usize) -> Result<usize> {
    if bound == usize::MAX - 1 && model.q_module(base).is_none() {
        Ok(model.height(base, true)?.saturating_sub(1))
    } else {
        Ok(bound)
    }
}

fn wide_degree(s: usize, weight: usize) -> i128 {
    i128::try_from(s as u128 * weight as u128).unwrap_or(i128::MAX)
}

fn increment(model: &mut Model, powers: &[i128], variable: usize) -> Result<Powers> {
    model.spend(1 + powers.len())?;
    let mut powers = powers.to_vec();
    powers[variable] = powers[variable]
        .checked_add(1)
        .ok_or_else(|| Error::plain("Koszul exponent is out of range"))?;
    Ok(powers.into())
}

fn primitive_action(
    model: &mut Model,
    base: TypeId,
    index: usize,
    degree: usize,
    input: Rc<[usize]>,
) -> Result<Rc<[usize]>> {
    model.spend(1)?;
    if input.is_empty() || degree == 0 {
        return Ok([].into());
    }
    if index == 0 {
        return square(model, base, degree, 1, input);
    }
    // Q_i = Sq^(2^i) Q_(i-1) + Q_(i-1) Sq^(2^i).
    let shift = 1usize << index;
    let q_first = primitive_action(model, base, index - 1, degree, input.clone())?;
    let left = square(model, base, degree + shift - 1, shift, q_first)?;
    let sq_first = square(model, base, degree, shift, input)?;
    let right = primitive_action(model, base, index - 1, degree + shift, sq_first)?;
    let result = xor(model, &left, &right)?;
    model.allocate(result.len())?;
    Ok(result.into())
}

fn square(
    model: &mut Model,
    base: TypeId,
    degree: usize,
    index: usize,
    input: Rc<[usize]>,
) -> Result<Rc<[usize]>> {
    if input.is_empty() || index > degree {
        return Ok([].into());
    }
    let word: steenrod::Word = [index].into();
    let Value::Class(output) = cohomology::act(model, base, degree, &[word], input)? else {
        unreachable!()
    };
    Ok(output)
}

fn polynomial_basis(
    model: &mut Model,
    weights: &[usize],
    s: usize,
    degree: u128,
    limit: usize,
) -> Result<Vec<Powers>> {
    let first = weights.first().copied().unwrap_or(0);
    let mut divisor = 0;
    let divisors = weights
        .iter()
        .map(|&weight| {
            divisor = gcd(divisor, weight - first);
            divisor as u128
        })
        .collect();
    let mut basis = PolynomialBasis {
        powers: vec![0; weights.len()],
        terms: Vec::new(),
        divisors,
        limit,
    };
    basis.enumerate(model, weights, s as u128, degree)?;
    Ok(basis.terms)
}

struct PolynomialBasis {
    powers: Vec<usize>,
    terms: Vec<Powers>,
    divisors: Vec<u128>,
    limit: usize,
}

impl PolynomialBasis {
    fn enumerate(
        &mut self,
        model: &mut Model,
        weights: &[usize],
        s: u128,
        degree: u128,
    ) -> Result<()> {
        model.spend(1)?;
        let Some((&last, rest)) = weights.split_last() else {
            if s == 0 && degree == 0 {
                self.push(model)?;
            }
            return Ok(());
        };
        let min = weights[0] as u128;
        let max = last as u128;
        if degree < min * s || degree > max * s {
            return Ok(());
        }
        let excess = degree - min * s;
        let divisor = self.divisors[rest.len()];
        if divisor != 0 && !excess.is_multiple_of(divisor) {
            return Ok(());
        }
        match rest {
            [] => {
                self.powers[0] = s as usize;
                self.push(model)
            }
            [_] => {
                let step = max - min;
                let high = excess / step;
                self.powers[0] = (s - high) as usize;
                self.powers[1] = high as usize;
                self.push(model)
            }
            [.., next] => {
                let next = *next as u128;
                let start = degree.saturating_sub(next * s).div_ceil(max - next);
                let end = excess / (max - min);
                for high in start..=end {
                    self.powers[rest.len()] = high as usize;
                    self.enumerate(model, rest, s - high, degree - high * max)?;
                }
                Ok(())
            }
        }
    }

    fn push(&mut self, model: &mut Model) -> Result<()> {
        if self.terms.len() == self.limit {
            return Err(basis_limit());
        }
        let cells = 1 + power_cells(self.powers.len());
        check_storage(cells * (1 + self.terms.len()))?;
        model.spend(cells)?;
        self.terms
            .push(self.powers.iter().map(|&e| e as i128).collect());
        Ok(())
    }
}

pub(crate) fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

pub(crate) fn basis(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    selected: &[usize],
    s: usize,
    t: i64,
    limit: usize,
) -> Result<Vec<Term>> {
    Complex::new(model, source, target, selected)?.basis(model, Filtration::Polynomial(s), t, limit)
}

pub(crate) fn basis_localized(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    selected: &[usize],
    s: i128,
    t: i64,
    limit: usize,
) -> Result<Vec<Term>> {
    v_koszul::top(selected)?;
    Complex::new(model, source, target, selected)?.basis(model, Filtration::Localized(s), t, limit)
}

pub(crate) fn differential(
    model: &mut Model,
    source: TypeId,
    target: TypeId,
    selected: &[usize],
    terms: &[Term],
) -> Result<Support> {
    if terms.is_empty() || selected.is_empty() {
        return Ok([].into());
    }
    let terms = Complex::new(model, source, target, selected)?.differential(model, terms)?;
    persist(model, terms)
}

pub(crate) fn restrict_words(
    model: &mut Model,
    selected: &[usize],
    words: &[Word],
) -> Result<Vec<Term>> {
    let unit = Generator {
        degree: 0,
        index: 0,
    };
    restrict_cobar(model, selected, words.iter().map(|word| (unit, unit, word)))
}

pub(crate) fn restrict_maps(
    model: &mut Model,
    selected: &[usize],
    terms: &[map_ext::Term],
) -> Result<Vec<Term>> {
    restrict_cobar(
        model,
        selected,
        terms
            .iter()
            .map(|term| (term.source, term.target, &term.word)),
    )
}

fn restrict_cobar<'a>(
    model: &mut Model,
    selected: &[usize],
    terms: impl Iterator<Item = (Generator, Generator, &'a Word)>,
) -> Result<Vec<Term>> {
    let mut support = Accumulator::default();
    for (source, target, word) in terms {
        model.spend(cobar::word_cells(word) + selected.len())?;
        let mut powers = vec![0; selected.len()];
        let mut keep = true;
        for factor in word.iter() {
            let Some((&1, zeros)) = factor.split_last() else {
                keep = false;
                break;
            };
            if zeros.iter().any(|&power| power != 0) {
                keep = false;
                break;
            }
            let Ok(variable) = selected.binary_search(&zeros.len()) else {
                keep = false;
                break;
            };
            powers[variable] += 1;
        }
        if keep {
            support.toggle(
                model,
                Term {
                    source,
                    target,
                    powers: powers.into(),
                },
            )?;
        }
    }
    Ok(support.finish())
}

pub(crate) fn restrict_subset(
    model: &mut Model,
    source: &[usize],
    target: &[usize],
    terms: &[Term],
) -> Result<Vec<Term>> {
    let positions: Vec<_> = source
        .iter()
        .map(|index| target.binary_search(index).is_ok())
        .collect();
    let mut output = Vec::new();
    for term in terms {
        model.spend(term.cells())?;
        if term
            .powers
            .iter()
            .zip(&positions)
            .any(|(&power, &keep)| power != 0 && !keep)
        {
            continue;
        }
        let powers = term
            .powers
            .iter()
            .zip(&positions)
            .filter_map(|(&power, &keep)| keep.then_some(power))
            .collect();
        output.push(Term {
            powers,
            ..term.clone()
        });
    }
    Ok(output)
}

pub(crate) fn restriction_type(model: &mut Model, ty: TypeId, selected: Indices) -> Result<TypeId> {
    let kind = model.kind(ty);
    let raw = kind.cobar_grade();
    let class = kind.ext_grade();
    let (grade, closed) = if let Some((profile, s, t)) = raw.or(class) {
        if let Profile::An(n) = profile {
            if selected.last().is_some_and(|&index| index > n) {
                return Err(Error::plain(
                    "qRes requires Milnor primitives in the source Steenrod subalgebra",
                ));
            }
        }
        let one = model.named("One").unwrap();
        let t =
            i64::try_from(t).map_err(|_| Error::plain("qRes internal degree is out of range"))?;
        (
            Grade {
                source: one,
                target: one,
                selected,
                s,
                t,
            },
            class.is_some(),
        )
    } else {
        match kind {
            &Kind::MapCobar(source, target, s, t) | &Kind::MapExt(source, target, s, t) => (
                Grade {
                    source,
                    target,
                    selected,
                    s,
                    t,
                },
                matches!(kind, Kind::MapExt(..)),
            ),
            Kind::QKoszul(grade) | Kind::QExt(grade) => {
                if selected
                    .iter()
                    .any(|index| grade.selected.binary_search(index).is_err())
                {
                    return Err(Error::plain(
                        "qRes requires a subset of the source Milnor primitives",
                    ));
                }
                (
                    Grade {
                        selected,
                        ..grade.clone()
                    },
                    matches!(kind, Kind::QExt(..)),
                )
            }
            _ => {
                return Err(Error::plain(
                    "qRes requires a mod-two cobar or Koszul cochain or class",
                ))
            }
        }
    };
    if closed {
        model.q_ext_type(grade)
    } else {
        model.q_koszul_type(grade)
    }
}

pub(crate) fn restrict(
    model: &mut Model,
    source: TypeId,
    value: Value,
    target: TypeId,
) -> Result<Value> {
    if source == target {
        return Ok(value);
    }
    let (grade, closed) = match model.kind(target) {
        Kind::QExt(grade) => (grade, true),
        Kind::QKoszul(grade) => (grade, false),
        _ => unreachable!(),
    };
    let selected = grade.selected.clone();
    let terms = match value {
        Value::Cobar(words) => restrict_words(model, &selected, &words)?,
        Value::MapCobar(terms) => restrict_maps(model, &selected, &terms)?,
        Value::QKoszul(terms) => {
            let Kind::QKoszul(grade) = model.kind(source) else {
                unreachable!()
            };
            let original = grade.selected.clone();
            restrict_subset(model, &original, &selected, &terms)?
        }
        Value::Class(class) if class.is_empty() => Vec::new(),
        Value::Class(class) => match model.kind(source) {
            Kind::QExt(grade) => {
                let original = grade.selected.clone();
                let terms = model.q_ext(source)?.representative(model, &class)?;
                restrict_subset(model, &original, &selected, &terms)?
            }
            Kind::MapExt(..) => {
                let terms = model.map_ext(source)?.representative(model, &class)?;
                restrict_maps(model, &selected, &terms)?
            }
            _ => {
                let words = model.adams(source)?.representative(model, &class)?;
                restrict_words(model, &selected, &words)?
            }
        },
        _ => unreachable!(),
    };
    if closed {
        class(model, target, &terms)
    } else {
        Ok(Value::QKoszul(persist(model, terms)?))
    }
}

pub(crate) fn composition_type(model: &mut Model, outer: TypeId, inner: TypeId) -> Result<TypeId> {
    let (outer, inner, closed) = match (model.kind(outer), model.kind(inner)) {
        (Kind::QExt(g), Kind::QExt(f)) => (g, f, true),
        (Kind::QKoszul(g), Kind::QKoszul(f)) => (g, f, false),
        (Kind::VExt(g), Kind::VExt(f)) => {
            let grade = v_koszul::composed(g, f)?;
            return model.v_ext_type(grade);
        }
        (Kind::VKoszul(g), Kind::VKoszul(f)) => {
            let grade = v_koszul::composed(g, f)?;
            return model.v_koszul_type(grade);
        }
        _ => {
            return Err(Error::plain(
                "ecomp requires two Koszul cochains or two Koszul classes in the same family",
            ))
        }
    };
    matching_endpoints(outer, inner)?;
    let s = outer
        .s
        .checked_add(inner.s)
        .ok_or_else(|| Error::plain("ecomp degree is out of range"))?;
    let t = outer
        .t
        .checked_add(inner.t)
        .ok_or_else(|| Error::plain("ecomp signed internal degree is out of range"))?;
    let grade = Grade {
        source: inner.source,
        target: outer.target,
        selected: outer.selected.clone(),
        s,
        t,
    };
    if closed {
        model.q_ext_type(grade)
    } else {
        model.q_koszul_type(grade)
    }
}

pub(crate) fn matching_endpoints<S>(outer: &Grade<S>, inner: &Grade<S>) -> Result<()> {
    if outer.source != inner.target {
        return Err(Error::plain("ecomp requires composable Koszul endpoints"));
    }
    if outer.selected != inner.selected {
        return Err(Error::plain("ecomp requires the same Milnor primitives"));
    }
    Ok(())
}

pub(crate) fn identity(model: &mut Model, ty: TypeId) -> Result<Value> {
    let group = model.q_ext(ty)?;
    model.spend(group.terms.len())?;
    let diagonal = group
        .terms
        .iter()
        .enumerate()
        .filter_map(|(index, term)| (term.source == term.target).then_some(index))
        .collect();
    let coordinates = group
        .quotient
        .project(model, diagonal)?
        .ok_or_else(|| Error::plain("Koszul identity is not closed"))?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let coordinates = model.q_ext(ty)?.project(model, terms)?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

pub(crate) fn representative(model: &mut Model, ty: TypeId, class: &[usize]) -> Result<Support> {
    if class.is_empty() {
        return Ok([].into());
    }
    let terms = model.q_ext(ty)?.representative(model, class)?;
    persist(model, terms)
}

pub(crate) fn compose_classes(
    model: &mut Model,
    (outer_ty, outer): (TypeId, &[usize]),
    (inner_ty, inner): (TypeId, &[usize]),
    ty: TypeId,
) -> Result<Value> {
    if outer.is_empty() || inner.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let outer = model.q_ext(outer_ty)?.representative(model, outer)?;
    let inner = model.q_ext(inner_ty)?.representative(model, inner)?;
    let terms = compose_terms(model, &outer, &inner)?;
    class(model, ty, &terms)
}

pub(crate) fn compose(model: &mut Model, outer: &[Term], inner: &[Term]) -> Result<Support> {
    let terms = compose_terms(model, outer, inner)?;
    persist(model, terms)
}

pub(crate) fn compose_terms(
    model: &mut Model,
    outer: &[Term],
    inner: &[Term],
) -> Result<Vec<Term>> {
    let mut support = Accumulator::default();
    let search = 1 + outer.len().max(1).ilog2() as usize;
    for f in inner {
        model.spend(search)?;
        let start = outer.partition_point(|g| g.source < f.target);
        for g in &outer[start..] {
            model.spend(1)?;
            if g.source != f.target {
                break;
            }
            model.spend(f.powers.len())?;
            let powers = f
                .powers
                .iter()
                .zip(g.powers.iter())
                .map(|(&a, &b)| {
                    a.checked_add(b)
                        .ok_or_else(|| Error::plain("Koszul exponent is out of range"))
                })
                .collect::<Result<Powers>>()?;
            support.toggle(
                model,
                Term {
                    source: f.source,
                    target: g.target,
                    powers,
                },
            )?;
        }
    }
    Ok(support.finish())
}

pub(crate) fn add(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Support> {
    model.spend(left.iter().chain(right).map(Term::cells).sum())?;
    persist(model, xor_sorted(left, right))
}

pub(crate) fn compare(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Ordering> {
    for (a, b) in left.iter().rev().zip(right.iter().rev()) {
        model.spend(a.cells())?;
        let comparison = a.cmp(b);
        if comparison != Ordering::Equal {
            return Ok(comparison);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

pub(crate) fn persist(model: &mut Model, terms: Vec<Term>) -> Result<Support> {
    let cells = terms.iter().map(Term::cells).sum();
    check_storage(cells)?;
    model.spend(cells)?;
    model.allocate(cells)?;
    Ok(terms.into())
}

#[derive(Default)]
struct Accumulator {
    terms: BTreeSet<Term>,
    cells: usize,
}

impl Accumulator {
    fn toggle(&mut self, model: &mut Model, term: Term) -> Result<()> {
        let cells = term.cells();
        model.spend(cells * (1 + self.terms.len().max(1).ilog2() as usize))?;
        if self.terms.remove(&term) {
            self.cells -= cells;
        } else {
            self.cells += cells;
            check_storage(self.cells)?;
            self.terms.insert(term);
        }
        Ok(())
    }

    fn finish(self) -> Vec<Term> {
        self.terms.into_iter().collect()
    }
}
