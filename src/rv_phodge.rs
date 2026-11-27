use std::cell::OnceCell;
use std::cmp::{Ordering, Reverse};
use std::collections::{btree_map::Entry, BTreeMap, BTreeSet};
use std::fmt;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_ore::gcd;
use crate::rv_qp_quad::{Element, Field};

pub(crate) struct Specification {
    pub(crate) prime: usize,
    pub(crate) eigenvalues: Vec<(i128, i128)>,
    pub(crate) monodromy: Vec<(usize, usize, i128, i128)>,
    pub(crate) filtration: Vec<(i64, Vec<i128>)>,
}

pub(crate) struct QuadraticSpec<T> {
    pub(crate) eigenvalues: Vec<T>,
    pub(crate) monodromy: Vec<(usize, usize, T)>,
    pub(crate) filtration: Vec<(i64, Vec<T>)>,
}

pub(crate) enum Query {
    Data(Option<Vec<usize>>),
    WeakAdmissibility,
    HarderNarasimhan,
}

type Column<T = i128> = Vec<(usize, T)>;
type WeightedColumns<T> = Vec<(i64, Column<T>)>;

enum Filtration {
    Rational(WeightedColumns<i128>),
    Quadratic(Field, WeightedColumns<Element>),
}

struct Operators {
    newton: Vec<i64>,
    required: Vec<Option<usize>>,
    chains: Vec<Vec<usize>>,
}

pub(crate) struct Module {
    newton: Vec<i64>,
    required: Vec<Option<usize>>,
    chains: Vec<Vec<usize>>,
    filtration: Filtration,
    polygon: OnceCell<Vec<Step>>,
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Rational {
    numerator: i128,
    denominator: i128,
}

impl Rational {
    fn new(numerator: i128, denominator: i128) -> Result<Self> {
        if denominator <= 0 {
            return Err(Error::plain(
                "filtered module rational denominator must be positive",
            ));
        }
        let content = gcd(numerator.unsigned_abs(), denominator as u128) as i128;
        Ok(Self {
            numerator: numerator / content,
            denominator: denominator / content,
        })
    }

    fn is_prime_multiple_of(self, other: Self, prime: Prime) -> bool {
        let p = i128::from(prime.get());
        let common = gcd(p as u128, other.denominator as u128) as i128;
        self.denominator == other.denominator / common
            && other.numerator.checked_mul(p / common) == Some(self.numerator)
    }

    fn valuation(self, model: &mut Model, prime: Prime) -> Result<i64> {
        let p = u128::from(prime.get());
        let a = valuation(model, self.numerator.unsigned_abs(), p)?;
        let b = valuation(model, self.denominator as u128, p)?;
        Ok(a - b)
    }
}

fn valuation(model: &mut Model, mut value: u128, p: u128) -> Result<i64> {
    let mut order = 0;
    while value.is_multiple_of(p) {
        model.spend(1)?;
        value /= p;
        order += 1;
    }
    Ok(order)
}

struct Data {
    hodge: Vec<i64>,
    newton: Vec<i64>,
}

impl Data {
    fn totals(&self) -> (i128, i128) {
        // The dimension is bounded by the finite matrix storage limit.
        (
            self.hodge.iter().map(|&x| i128::from(x)).sum(),
            self.newton.iter().map(|&x| i128::from(x)).sum(),
        )
    }

    fn degree(&self) -> i128 {
        let (h, n) = self.totals();
        h - n
    }
}

#[derive(Clone, Copy)]
struct Slope {
    degree: i128,
    dimension: usize,
}

impl Slope {
    fn compare(self, other: Self) -> Result<Ordering> {
        let a = self
            .degree
            .checked_mul(other.dimension as i128)
            .ok_or_else(arithmetic)?;
        let b = other
            .degree
            .checked_mul(self.dimension as i128)
            .ok_or_else(arithmetic)?;
        Ok(a.cmp(&b))
    }
}

impl fmt::Display for Slope {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        let divisor = gcd(self.degree.unsigned_abs(), self.dimension as u128) as i128;
        let numerator = self.degree / divisor;
        let denominator = self.dimension as i128 / divisor;
        if denominator == 1 {
            write!(out, "{numerator}")
        } else {
            write!(out, "{numerator}/{denominator}")
        }
    }
}

struct State {
    prefixes: Vec<usize>,
    dimension: usize,
    degree: i128,
}

struct Step {
    coordinates: Vec<usize>,
    slope: Slope,
}

impl Module {
    pub(crate) fn build(model: &mut Model, spec: Specification) -> Result<Self> {
        let p = spec.prime;
        if p < 2 || (p != 2 && p.is_multiple_of(2)) {
            return Err(Error::plain(
                "filtered module residue characteristic must be prime",
            ));
        }
        let prime = if p == 2 {
            Prime::validated(2)
        } else {
            Prime::new(model, p)?
        };
        let n = spec.eigenvalues.len();
        validate_width(n, spec.filtration.len(), 2)?;
        let monodromy = spec
            .monodromy
            .into_iter()
            .map(|(i, j, a, b)| (i, j, (a, b)))
            .collect();
        let operators = validate_operators(
            model,
            spec.eigenvalues,
            monodromy,
            |(a, b)| Rational::new(a, b),
            |c| c.numerator == 0,
            |model, c| c.valuation(model, prime),
            |source, target| source.is_prime_multiple_of(target, prime),
        )?;
        let mut span = RationalBasis::default();
        let columns = build_flag(
            model,
            &operators,
            spec.filtration,
            3,
            |entries| primitive(sparse_column(entries, |c| c == 0)),
            |model, vector, cells| span.insert(model, vector, cells),
        )?;
        Ok(Self::with_filtration(
            operators,
            Filtration::Rational(columns),
        ))
    }

    pub(crate) fn build_quadratic(
        model: &mut Model,
        field: Field,
        spec: QuadraticSpec<Element>,
    ) -> Result<Self> {
        let n = spec.eigenvalues.len();
        validate_width(n, spec.filtration.len(), 6)?;
        let operators = validate_operators(
            model,
            spec.eigenvalues,
            spec.monodromy,
            Ok,
            Element::is_zero,
            |model, c| {
                field
                    .valuation(model, c)?
                    .ok_or_else(|| Error::plain("zero Frobenius eigenvalue"))
            },
            |source, target| source.is_prime_multiple_of(target, field.prime()),
        )?;
        let mut span = QuadraticBasis::default();
        let columns = build_flag(
            model,
            &operators,
            spec.filtration,
            7,
            |entries| sparse_column(entries, Element::is_zero),
            |model, vector, cells| span.insert(model, field, vector, cells),
        )?;
        Ok(Self::with_filtration(
            operators,
            Filtration::Quadratic(field, columns),
        ))
    }

    fn with_filtration(operators: Operators, filtration: Filtration) -> Self {
        Self {
            newton: operators.newton,
            required: operators.required,
            chains: operators.chains,
            filtration,
            polygon: OnceCell::new(),
        }
    }

    fn dimension(&self) -> usize {
        self.newton.len()
    }

    fn checked_subset(&self, model: &mut Model, mut coordinates: Vec<usize>) -> Result<Vec<usize>> {
        model.spend(coordinates.len() * (1 + coordinates.len().max(1).ilog2() as usize))?;
        coordinates.sort_unstable();
        if coordinates.last().is_some_and(|&i| i >= self.dimension()) {
            return Err(Error::plain("filtered subspace coordinate is out of range"));
        }
        if coordinates.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(Error::plain("filtered subspace coordinate is repeated"));
        }
        for &i in &coordinates {
            if let Some(j) = self.required[i] {
                if coordinates.binary_search(&j).is_err() {
                    return Err(Error::plain(format!(
                        "filtered subspace is not phi,N-stable: coordinate {i} requires {j}"
                    )));
                }
            }
        }
        Ok(coordinates)
    }

    fn data(&self, model: &mut Model, coordinates: &[usize], retained: usize) -> Result<Data> {
        let n = self.dimension();
        let width = coordinates.len();
        model.spend(1 + n + width * (1 + width.max(1).ilog2() as usize))?;
        let mut newton: Vec<_> = coordinates.iter().map(|&i| self.newton[i]).collect();
        newton.sort_unstable();
        if coordinates.is_empty() {
            return Ok(Data {
                hodge: Vec::new(),
                newton,
            });
        }
        if coordinates.len() == n {
            let hodge = match &self.filtration {
                Filtration::Rational(columns) => all_weights(columns),
                Filtration::Quadratic(_, columns) => all_weights(columns),
            };
            return Ok(Data { hodge, newton });
        }
        let mut selected = vec![false; n];
        for &i in coordinates {
            selected[i] = true;
        }
        let retained = retained + 5 * n;
        let hodge = match &self.filtration {
            Filtration::Rational(columns) => {
                let mut span = RationalBasis::default();
                intersection_weights(
                    model,
                    columns,
                    &selected,
                    retained,
                    |model, vector, cells| span.insert(model, vector, cells),
                )?
            }
            Filtration::Quadratic(field, columns) => {
                let mut span = QuadraticBasis::default();
                intersection_weights(
                    model,
                    columns,
                    &selected,
                    retained,
                    |model, vector, cells| span.insert(model, *field, vector, cells),
                )?
            }
        };
        Ok(Data { hodge, newton })
    }

    fn coordinates(&self, model: &mut Model, prefixes: &[usize]) -> Result<Vec<usize>> {
        let mut selected = Vec::new();
        for (chain, &len) in self.chains.iter().zip(prefixes) {
            selected.extend_from_slice(&chain[..len]);
        }
        model.spend(
            1 + prefixes.len() + selected.len() * (1 + selected.len().max(1).ilog2() as usize),
        )?;
        selected.sort_unstable();
        Ok(selected)
    }

    fn states(&self, model: &mut Model) -> Result<Vec<State>> {
        let count = self.chains.iter().try_fold(1usize, |count, chain| {
            count.checked_mul(chain.len() + 1).ok_or_else(storage_error)
        })?;
        let storage = count.saturating_mul(7 + self.chains.len());
        check_storage(storage)?;
        model.spend(count)?;
        let mut states = Vec::with_capacity(count);
        let mut current = vec![0; self.chains.len()];
        for _ in 0..count {
            let coordinates = self.coordinates(model, &current)?;
            let degree = self.data(model, &coordinates, storage)?.degree();
            states.push(State {
                prefixes: current.clone(),
                dimension: coordinates.len(),
                degree,
            });
            for (prefix, chain) in current.iter_mut().zip(&self.chains) {
                if *prefix < chain.len() {
                    *prefix += 1;
                    break;
                }
                *prefix = 0;
            }
        }
        Ok(states)
    }

    fn h_n(&self, model: &mut Model) -> Result<&[Step]> {
        if let Some(polygon) = self.polygon.get() {
            return Ok(polygon);
        }
        let states = self.states(model)?;
        let state_cells = states.len() * (7 + self.chains.len());
        let mut flag_cells = 0;
        let mut flag = Vec::new();
        let mut lower = &states[0];
        while lower.dimension < self.dimension() {
            let mut best: Option<(&State, Slope)> = None;
            for state in &states {
                model.spend(1 + state.prefixes.len())?;
                if state.dimension <= lower.dimension
                    || state
                        .prefixes
                        .iter()
                        .zip(&lower.prefixes)
                        .any(|(a, b)| a < b)
                {
                    continue;
                }
                let slope = Slope {
                    degree: state.degree - lower.degree,
                    dimension: state.dimension - lower.dimension,
                };
                let update = if let Some((old, rate)) = best {
                    match slope.compare(rate)? {
                        Ordering::Greater => true,
                        Ordering::Equal => state.dimension > old.dimension,
                        Ordering::Less => false,
                    }
                } else {
                    true
                };
                if update {
                    best = Some((state, slope));
                }
            }
            let Some((next, slope)) = best else {
                return Err(Error::plain(
                    "filtered module flag cannot reach the full module",
                ));
            };
            // Maximal-slope states are closed under union by supermodularity
            // of degree. The largest dimension therefore selects a unique
            // greatest state among all maximizers.
            let coordinates = self.coordinates(model, &next.prefixes)?;
            flag_cells += 7 + coordinates.len();
            check_storage(state_cells + flag_cells)?;
            flag.push(Step { coordinates, slope });
            lower = next;
        }
        model.allocate(flag_cells)?;
        Ok(self.polygon.get_or_init(|| flag))
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, query: Query) -> Result<String> {
        match query {
            Query::Data(subset) => {
                let full = subset.is_none();
                let subset = match subset {
                    Some(coordinates) => self.checked_subset(model, coordinates)?,
                    None => (0..self.dimension()).collect(),
                };
                let data = self.data(model, &subset, 0)?;
                let (h, n) = data.totals();
                let name = if full {
                    name.to_owned()
                } else {
                    format!("{name}, {subset:?}")
                };
                Ok(format!(
                    "rvPHodgeData[{name}] = tH {h}, tN {n}, degree {}, Hodge {:?}, Newton {:?}",
                    h - n,
                    data.hodge,
                    data.newton
                ))
            }
            Query::WeakAdmissibility => {
                let full: Vec<_> = (0..self.dimension()).collect();
                let (h, n) = self.data(model, &full, 0)?.totals();
                if h != n {
                    return Ok(format!(
                        "rvPHodgeWeak[{name}] = degree equality fails on {full:?} (tH {h}, tN {n})"
                    ));
                }
                if let Some(step) = self.h_n(model)?.first().filter(|s| s.slope.degree > 0) {
                    let (h, n) = self.data(model, &step.coordinates, 0)?.totals();
                    return Ok(format!(
                        "rvPHodgeWeak[{name}] = destabilized by {:?} (tH {h}, tN {n})",
                        step.coordinates
                    ));
                }
                Ok(format!("rvPHodgeWeak[{name}] = weakly admissible"))
            }
            Query::HarderNarasimhan => {
                let mut shown = Vec::new();
                for step in self.h_n(model)? {
                    model.spend(1 + step.coordinates.len())?;
                    shown.push(format!("({:?}, {})", step.coordinates, step.slope));
                }
                Ok(format!("rvPHodgeHN[{name}] = [{}]", shown.join(", ")))
            }
        }
    }
}

fn validate_width(n: usize, columns: usize, entry_cells: usize) -> Result<()> {
    check_storage(n.saturating_mul(n.saturating_mul(entry_cells).saturating_add(16)))?;
    if columns != n {
        return Err(Error::plain(format!(
            "filtered module requires {n} filtration columns"
        )));
    }
    Ok(())
}

fn validate_operators<A, T: Copy + Ord>(
    model: &mut Model,
    input: Vec<A>,
    monodromy: Vec<(usize, usize, A)>,
    normalize: impl Fn(A) -> Result<T>,
    zero: impl Fn(T) -> bool,
    mut valuation: impl FnMut(&mut Model, T) -> Result<i64>,
    ratio: impl Fn(T, T) -> bool,
) -> Result<Operators> {
    let n = input.len();
    let mut distinct = BTreeSet::new();
    let mut newton = Vec::with_capacity(n);
    let mut eigenvalues = Vec::with_capacity(n);
    for entry in input {
        let eigenvalue = normalize(entry)?;
        model.spend(1 + distinct.len().max(1).ilog2() as usize)?;
        if zero(eigenvalue) || !distinct.insert(eigenvalue) {
            return Err(Error::plain(
                "filtered module requires distinct nonzero Frobenius eigenvalues",
            ));
        }
        newton.push(valuation(model, eigenvalue)?);
        eigenvalues.push(eigenvalue);
    }
    let mut required = vec![None; n];
    let mut higher = vec![None; n];
    let mut supplied = BTreeSet::new();
    for (row, column, entry) in monodromy {
        if row >= n || column >= n {
            return Err(Error::plain(format!(
                "monodromy entry ({row},{column}) is outside dimension {n}"
            )));
        }
        model.spend(1 + supplied.len().max(1).ilog2() as usize)?;
        if !supplied.insert((row, column)) {
            return Err(Error::plain(format!(
                "duplicate monodromy entry ({row},{column})"
            )));
        }
        if zero(normalize(entry)?) {
            continue;
        }
        if !ratio(eigenvalues[column], eigenvalues[row]) {
            return Err(Error::plain(format!(
                "monodromy fails N phi = p phi N at entry ({row},{column})"
            )));
        }
        required[column] = Some(row);
        higher[row] = Some(column);
    }

    // Each eigenvalue has at most one p-multiple and one p-division
    // among the others. Following required coordinates strictly lowers
    // the valuation. The graph is therefore a disjoint union of chains.
    let mut chains = Vec::new();
    model.spend(2 * n)?;
    for (low, required) in required.iter().enumerate() {
        if required.is_some() {
            continue;
        }
        let mut chain = vec![low];
        let mut current = low;
        while let Some(index) = higher[current] {
            chain.push(index);
            current = index;
        }
        chains.push(chain);
    }

    debug_assert_eq!(chains.iter().map(Vec::len).sum::<usize>(), n);
    Ok(Operators {
        newton,
        required,
        chains,
    })
}

fn build_flag<T: Copy>(
    model: &mut Model,
    operators: &Operators,
    filtration: Vec<(i64, Vec<T>)>,
    entry_cells: usize,
    sparse: impl Fn(Vec<T>) -> Column<T>,
    mut independent: impl FnMut(&mut Model, Column<T>, usize) -> Result<bool>,
) -> Result<WeightedColumns<T>> {
    let n = operators.newton.len();
    let mut columns = Vec::with_capacity(n);
    let mut cells = 8 + 8 * n + 3 * operators.chains.len();
    for (weight, entries) in filtration {
        if entries.len() != n {
            return Err(Error::plain(format!(
                "each filtration column must have {n} entries"
            )));
        }
        model.spend(n)?;
        let column = sparse(entries);
        cells += 3 + entry_cells * column.len();
        if !independent(model, column.clone(), cells)? {
            return Err(Error::plain(
                "filtration columns must be a basis over characteristic zero",
            ));
        }
        columns.push((weight, column));
    }
    model.spend(n * (1 + n.max(1).ilog2() as usize))?;
    columns.sort_by_key(|(weight, _)| Reverse(*weight));
    check_storage(cells)?;
    model.allocate(cells)?;
    Ok(columns)
}

fn sparse_column<T: Copy>(entries: Vec<T>, zero: impl Fn(T) -> bool) -> Column<T> {
    entries
        .into_iter()
        .enumerate()
        .filter(|&(_, c)| !zero(c))
        .collect()
}

fn all_weights<T>(columns: &[(i64, Column<T>)]) -> Vec<i64> {
    columns.iter().rev().map(|(w, _)| *w).collect()
}

fn intersection_weights<T: Copy>(
    model: &mut Model,
    columns: &[(i64, Column<T>)],
    selected: &[bool],
    retained: usize,
    mut independent: impl FnMut(&mut Model, Column<T>, usize) -> Result<bool>,
) -> Result<Vec<i64>> {
    let mut hodge = Vec::new();
    for (weight, column) in columns {
        model.spend(1 + column.len())?;
        let projection = column
            .iter()
            .copied()
            .filter(|&(i, _)| !selected[i])
            .collect();
        // An added column increases the flag dimension by one. If its
        // projection is dependent, the intersection gains this weight.
        if !independent(model, projection, retained)? {
            hodge.push(*weight);
        }
    }
    hodge.reverse();
    Ok(hodge)
}

#[derive(Default)]
struct QuadraticBasis {
    pivots: BTreeMap<usize, Column<Element>>,
    cells: usize,
}

impl QuadraticBasis {
    fn insert(
        &mut self,
        model: &mut Model,
        field: Field,
        mut vector: Column<Element>,
        retained: usize,
    ) -> Result<bool> {
        while let Some(&(coordinate, factor)) = vector.last() {
            model.spend(1 + self.pivots.len().max(1).ilog2() as usize)?;
            check_storage(retained + self.cells + 7 * vector.len())?;
            match self.pivots.entry(coordinate) {
                Entry::Vacant(entry) => {
                    self.cells += 4 + 7 * vector.len();
                    check_storage(retained + self.cells)?;
                    let last = vector.len() - 1;
                    if last != 0 && factor != Element::ONE {
                        let inverse = field.inverse(model, factor)?;
                        for (_, value) in &mut vector[..last] {
                            *value = field.multiply(model, *value, inverse)?;
                        }
                    }
                    vector[last].1 = Element::ONE;
                    entry.insert(vector);
                    return Ok(true);
                }
                Entry::Occupied(entry) => {
                    let pivot = entry.get();
                    vector.pop();
                    if pivot.len() == 1 {
                        continue;
                    }
                    let width = vector.len() + pivot.len();
                    check_storage(retained + self.cells + 14 * width)?;
                    model.spend(width * (1 + width.max(1).ilog2() as usize))?;
                    let mut result: BTreeMap<_, _> = vector.into_iter().collect();
                    for &(i, c) in &pivot[..pivot.len() - 1] {
                        let term = field.multiply(model, c, factor)?;
                        let value = result.entry(i).or_insert(Element::ZERO);
                        *value = field.subtract(model, *value, term)?;
                    }
                    vector = result.into_iter().filter(|&(_, c)| !c.is_zero()).collect();
                }
            }
        }
        Ok(false)
    }
}

#[derive(Default)]
struct RationalBasis {
    pivots: BTreeMap<usize, Column>,
    cells: usize,
}

impl RationalBasis {
    fn insert(&mut self, model: &mut Model, vector: Column, retained: usize) -> Result<bool> {
        let mut vector = primitive(vector);
        while let Some(&(coordinate, factor)) = vector.last() {
            model.spend(1 + self.pivots.len().max(1).ilog2() as usize)?;
            check_storage(retained + self.cells + 3 * vector.len())?;
            match self.pivots.entry(coordinate) {
                Entry::Vacant(entry) => {
                    self.cells += 4 + 3 * vector.len();
                    check_storage(retained + self.cells)?;
                    entry.insert(vector);
                    return Ok(true);
                }
                Entry::Occupied(entry) => {
                    let pivot = entry.get();
                    vector.pop();
                    if pivot.len() == 1 {
                        continue;
                    }
                    let leading = pivot.last().unwrap().1;
                    let common = gcd(leading.unsigned_abs(), factor.unsigned_abs());
                    let scale = divide_content(leading, common);
                    let subtract = divide_content(factor, common);
                    let width = vector.len() + pivot.len();
                    check_storage(retained + self.cells + 8 * width)?;
                    model.spend(width * (1 + width.max(1).ilog2() as usize))?;
                    let mut result = BTreeMap::new();
                    for &(i, c) in &vector {
                        result.insert(i, c.checked_mul(scale).ok_or_else(arithmetic)?);
                    }
                    for &(i, c) in &pivot[..pivot.len() - 1] {
                        let term = c.checked_mul(subtract).ok_or_else(arithmetic)?;
                        let value = result.entry(i).or_insert(0);
                        *value = value.checked_sub(term).ok_or_else(arithmetic)?;
                    }
                    vector = primitive(result.into_iter().filter(|(_, c)| *c != 0).collect());
                }
            }
        }
        Ok(false)
    }
}

fn primitive(mut vector: Column) -> Column {
    let content = vector
        .iter()
        .fold(0, |common, (_, c)| gcd(common, c.unsigned_abs()));
    if content > 1 {
        for (_, c) in &mut vector {
            *c = divide_content(*c, content);
        }
    }
    vector
}

fn divide_content(value: i128, divisor: u128) -> i128 {
    // 2^127 divides only i128::MIN among the nonzero signed inputs.
    i128::try_from(divisor).map_or(-1, |divisor| value / divisor)
}

fn arithmetic() -> Error {
    Error::plain("filtered module exact integer arithmetic is out of range")
}

fn storage_error() -> Error {
    Error::plain("filtered module calculation exceeds the storage limit")
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(storage_error());
    }
    Ok(())
}
