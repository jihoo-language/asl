use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::{btree_map::Entry, BTreeMap, HashMap};
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;

pub(crate) type Literal = (Vec<i128>, i128);

#[derive(Clone, Copy, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Level {
    pub(crate) prime: u32,
    pub(crate) precision: usize,
    pub(crate) modulus: u128,
}

impl Level {
    pub(crate) fn new(model: &mut Model, prime: usize, precision: usize) -> Result<Self> {
        if precision == 0 {
            return Err(Error::plain("Cartier precision must be positive"));
        }
        let prime = match prime {
            2 => 2,
            n if n < 2 || n % 2 == 0 => return Err(Error::plain(format!("{n} is not prime"))),
            n => Prime::new(model, n)?.get(),
        };
        let mut modulus = 1u128;
        for _ in 0..precision {
            model.spend(1)?;
            modulus = modulus
                .checked_mul(u128::from(prime))
                .ok_or_else(|| Error::plain("the Cartier prime-power modulus exceeds u128"))?;
        }
        Ok(Self {
            prime,
            precision,
            modulus,
        })
    }

    fn denominator(self) -> u128 {
        self.modulus / u128::from(self.prime)
    }

    pub(crate) fn check_digit(self, digit: usize) -> Result<u32> {
        if digit >= self.prime as usize {
            return Err(Error::plain(format!(
                "Cartier digit {digit} is outside 0..{}",
                self.prime
            )));
        }
        Ok(digit as u32)
    }

    pub(crate) fn add(self, a: u128, b: u128) -> u128 {
        let gap = self.modulus - b;
        if a >= gap {
            a - gap
        } else {
            a + b
        }
    }

    fn reduce(self, integer: i128) -> u128 {
        let residue = integer.unsigned_abs() % self.modulus;
        if integer < 0 && residue != 0 {
            self.modulus - residue
        } else {
            residue
        }
    }
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) period: Rc<str>,
    pub(crate) level: Level,
}

#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Monomial {
    pub(crate) power: u128,
    pub(crate) exponents: Rc<[i128]>,
}

impl Monomial {
    fn unit(variables: usize) -> Self {
        Self {
            power: 0,
            exponents: vec![0; variables].into(),
        }
    }

    fn cells(&self) -> usize {
        6 + 2 * self.exponents.len()
    }

    fn product(&self, other: &Self) -> Result<Self> {
        let power = self
            .power
            .checked_add(other.power)
            .ok_or_else(exponent_overflow)?;
        let exponents = self
            .exponents
            .iter()
            .zip(other.exponents.iter())
            .map(|(&x, &y)| x.checked_add(y).ok_or_else(exponent_overflow))
            .collect::<Result<Rc<[_]>>>()?;
        Ok(Self { power, exponents })
    }

    fn sections(&self, other: &Self, prime: u32, digit: u32) -> Result<Self> {
        let p = i128::from(prime);
        let power = (self
            .power
            .checked_add(other.power)
            .ok_or_else(exponent_overflow)?
            - u128::from(digit))
            / u128::from(prime);
        let exponents = self
            .exponents
            .iter()
            .zip(other.exponents.iter())
            .map(|(&x, &y)| {
                let carry = (x.rem_euclid(p) + y.rem_euclid(p)) / p;
                x.div_euclid(p)
                    .checked_add(y.div_euclid(p))
                    .and_then(|q| q.checked_add(carry))
                    .ok_or_else(exponent_overflow)
            })
            .collect::<Result<Rc<[_]>>>()?;
        Ok(Self { power, exponents })
    }
}

pub(crate) type State = Rc<[(Monomial, u128)]>;

pub(crate) struct Period {
    variables: usize,
    laurent: Vec<(Monomial, i128)>,
    bounds: Vec<(i128, i128)>,
    initials: RefCell<BTreeMap<Level, State>>,
    multipliers: RefCell<BTreeMap<(Level, u32), Rc<Multiplier>>>,
}

impl Period {
    pub(crate) fn build(model: &mut Model, variables: usize, terms: Vec<Literal>) -> Result<Self> {
        check_cells(6usize.saturating_add(4usize.saturating_mul(variables)))?;
        let mut exact = Accumulator::default();
        for (exponents, coefficient) in terms {
            if exponents.len() != variables {
                return Err(Error::plain(format!(
                    "period term requires {variables} spatial exponents, found {}",
                    exponents.len()
                )));
            }
            exact.add(
                model,
                Monomial {
                    power: 0,
                    exponents: exponents.into(),
                },
                coefficient,
                Exact,
            )?;
        }
        let laurent = exact.finish();
        let mut bounds = vec![(0i128, 0i128); variables];
        for (term, _) in &laurent {
            model.spend(1 + variables)?;
            for ((lower, upper), &exponent) in bounds.iter_mut().zip(term.exponents.iter()) {
                *lower = (*lower).min(exponent);
                *upper = (*upper).max(exponent);
            }
        }
        model.allocate(6 + 4 * variables + cells(&laurent))?;
        Ok(Self {
            variables,
            laurent,
            bounds,
            initials: RefCell::default(),
            multipliers: RefCell::default(),
        })
    }

    pub(crate) fn exact(&self, model: &mut Model, exponent: u128) -> Result<i128> {
        if exponent == 0 {
            return Ok(1);
        }
        model.spend(cells(&self.laurent))?;
        let output = power(model, self.variables, exponent, self.laurent.clone(), Exact)?;
        constant(model, self.variables, &output)
    }

    fn reduced(
        &self,
        model: &mut Model,
        level: Level,
        denominator: bool,
    ) -> Result<Vec<(Monomial, u128)>> {
        model.spend(cells(&self.laurent))?;
        let mut result = Vec::new();
        if denominator {
            result.push((Monomial::unit(self.variables), 1));
        }
        for (monomial, coefficient) in &self.laurent {
            let mut value = level.reduce(*coefficient);
            if value == 0 {
                continue;
            }
            let mut term = monomial.clone();
            if denominator {
                term.power = 1;
                value = level.modulus - value;
            }
            result.push((term, value));
        }
        check_cells(cells(&result))?;
        Ok(result)
    }

    pub(crate) fn initial(&self, model: &mut Model, level: Level) -> Result<State> {
        if let Some(state) = self.initials.borrow().get(&level) {
            return Ok(state.clone());
        }
        let exponent = level.denominator() - 1;
        let base = if exponent == 0 {
            Vec::new()
        } else {
            self.reduced(model, level, true)?
        };
        let output = power(model, self.variables, exponent, base, level)?;
        let state = persist(model, output)?;
        model.allocate(6)?;
        self.initials.borrow_mut().insert(level, state.clone());
        Ok(state)
    }

    fn multiplier(&self, model: &mut Model, level: Level, digit: u32) -> Result<Rc<Multiplier>> {
        let over_field = level.precision == 1;
        let key = (level, if over_field { digit } else { 0 });
        if let Some(multiplier) = self.multipliers.borrow().get(&key) {
            return Ok(multiplier.clone());
        }
        // When M=1, the z^d coefficient of Q^(p-1) is L^d modulo p.
        let exponent = if over_field {
            u128::from(digit)
        } else {
            level.modulus - level.denominator()
        };
        let base = if exponent == 0 {
            Vec::new()
        } else {
            self.reduced(model, level, !over_field)?
        };
        let terms = power(model, self.variables, exponent, base, level)?;
        let multiplier = Rc::new(Multiplier::new(model, level.prime, terms)?);
        self.multipliers
            .borrow_mut()
            .insert(key, multiplier.clone());
        Ok(multiplier)
    }

    pub(crate) fn section(
        &self,
        model: &mut Model,
        level: Level,
        digit: u32,
        state: &[(Monomial, u128)],
    ) -> Result<Vec<(Monomial, u128)>> {
        if state.is_empty() {
            return Ok(Vec::new());
        }
        let multiplier = self.multiplier(model, level, digit)?;
        let digit = if level.precision == 1 { 0 } else { digit };
        let mut output = Accumulator::default();
        for (monomial, a) in state {
            model.spend(
                (1 + self.variables) * (1 + multiplier.bins.len().max(1).ilog2() as usize),
            )?;
            let key = residues(monomial, level.prime, Some(digit));
            if let Some(bin) = multiplier.bins.get(&key) {
                for (other, b) in bin {
                    model.spend(1)?;
                    let value = level.mul(model, *a, *b)?;
                    if value != 0 {
                        let exponent = monomial.sections(other, level.prime, digit)?;
                        output.add(model, exponent, value, level)?;
                    }
                }
            }
        }
        Ok(output.finish())
    }

    pub(crate) fn read(&self, model: &mut Model, state: &[(Monomial, u128)]) -> Result<u128> {
        constant(model, self.variables, state)
    }

    pub(crate) fn query(&self, model: &mut Model, level: Level, mut index: u128) -> Result<u128> {
        if index == 0 {
            return Ok(1);
        }
        let prime = u128::from(level.prime);
        let mut state = self.initial(model, level)?;
        // A finite-state path can visit the same numerator with the same
        // digit again. Cache verified exact transitions within this query.
        let mut seen: HashMap<(u32, State), State> = HashMap::new();
        let mut stored = 0;
        while index != 0 && !state.is_empty() {
            let digit = (index % prime) as u32;
            let weight = cells(&state);
            model.spend(1 + weight)?;
            let key = (digit, state.clone());
            state = if let Some(after) = seen.get(&key) {
                after.clone()
            } else {
                let after = self.section(model, level, digit, &state)?;
                let size = stored + 8 + weight + cells(&after);
                let after: State = after.into();
                if size <= MAX_PIVOT_CELLS {
                    model.spend(1 + weight)?;
                    seen.insert(key, after.clone());
                    stored = size;
                }
                after
            };
            index /= prime;
        }
        self.read(model, &state)
    }

    pub(crate) fn elements(
        &self,
        model: &mut Model,
        level: Level,
        limit: usize,
    ) -> Result<Vec<Value>> {
        let order = usize::try_from(level.modulus)
            .ok()
            .filter(|&m| m <= limit)
            .ok_or_else(enumeration_limit)?;
        let mut slots = 0;
        let mut cardinality = 1;
        while cardinality <= limit / order {
            cardinality *= order;
            slots += 1;
        }
        let denominator = level.denominator();
        let mut dimensions = usize::try_from(denominator)
            .ok()
            .filter(|&x| x <= slots)
            .ok_or_else(enumeration_limit)?;
        let multiplier = i128::try_from(denominator).map_err(|_| enumeration_limit())?;
        let mut ranges = Vec::new();
        for &(lower, upper) in &self.bounds {
            let lo = lower
                .checked_mul(multiplier)
                .ok_or_else(enumeration_limit)?;
            let hi = upper
                .checked_mul(multiplier)
                .ok_or_else(enumeration_limit)?;
            let size = hi
                .checked_sub(lo)
                .and_then(|d| d.checked_add(1))
                .and_then(|n| usize::try_from(n).ok())
                .ok_or_else(enumeration_limit)?;
            dimensions = dimensions
                .checked_mul(size)
                .filter(|&n| n <= slots)
                .ok_or_else(enumeration_limit)?;
            ranges.push((lo, hi));
        }
        check_cells(dimensions.saturating_mul(6 + 2 * self.variables))?;
        let mut monomials = Vec::with_capacity(dimensions);
        let mut location: Vec<_> = ranges.iter().map(|&(lo, _)| lo).collect();
        for z in 0..denominator {
            loop {
                model.spend(1 + self.variables)?;
                monomials.push(Monomial {
                    power: z,
                    exponents: location.as_slice().into(),
                });
                let Some(i) = (0..location.len())
                    .rev()
                    .find(|&i| location[i] < ranges[i].1)
                else {
                    for (entry, &(lo, _)) in location.iter_mut().zip(&ranges) {
                        *entry = lo;
                    }
                    break;
                };
                location[i] += 1;
                for (entry, &(lo, _)) in location[i + 1..].iter_mut().zip(&ranges[i + 1..]) {
                    *entry = lo;
                }
            }
        }
        let count = order.pow(dimensions as u32);
        model.spend(count * (1 + dimensions))?;
        model.allocate(count)?;
        (0..count)
            .map(|mut code| {
                let mut state = Vec::new();
                for monomial in &monomials {
                    let c = code % order;
                    code /= order;
                    if c != 0 {
                        state.push((monomial.clone(), c as u128));
                    }
                }
                Ok(Value::PeriodState(persist(model, state)?))
            })
            .collect()
    }
}

struct Multiplier {
    bins: BTreeMap<Vec<u32>, Vec<(Monomial, u128)>>,
}

impl Multiplier {
    fn new(model: &mut Model, prime: u32, terms: Vec<(Monomial, u128)>) -> Result<Self> {
        let mut stored = 6 + cells(&terms);
        let mut bins: BTreeMap<Vec<u32>, Vec<_>> = BTreeMap::new();
        for (monomial, coefficient) in terms {
            model
                .spend((1 + monomial.exponents.len()) * (1 + bins.len().max(1).ilog2() as usize))?;
            let key = residues(&monomial, prime, None);
            let bin = match bins.entry(key) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => {
                    stored += 7 + monomial.exponents.len();
                    check_cells(stored)?;
                    entry.insert(Vec::new())
                }
            };
            bin.push((monomial, coefficient));
        }
        model.allocate(stored)?;
        Ok(Self { bins })
    }
}

fn residues(monomial: &Monomial, prime: u32, needed: Option<u32>) -> Vec<u32> {
    let z = (monomial.power % u128::from(prime)) as u32;
    let spatial = monomial
        .exponents
        .iter()
        .map(|&a| a.rem_euclid(i128::from(prime)) as u32);
    let target = |wanted: u32, found: u32| {
        if wanted >= found {
            wanted - found
        } else {
            prime - (found - wanted)
        }
    };
    match needed {
        Some(digit) => std::iter::once(target(digit, z))
            .chain(spatial.map(|r| target(0, r)))
            .collect(),
        None => std::iter::once(z).chain(spatial).collect(),
    }
}

fn constant<C: Copy + Default>(
    model: &mut Model,
    variables: usize,
    terms: &[(Monomial, C)],
) -> Result<C> {
    model.spend((1 + variables) * (1 + terms.len().max(1).ilog2() as usize))?;
    let index = terms.binary_search_by(|(monomial, _)| {
        monomial.power.cmp(&0).then_with(|| {
            monomial
                .exponents
                .iter()
                .find(|&&x| x != 0)
                .map_or(Ordering::Equal, |x| x.cmp(&0))
        })
    });
    Ok(index.map_or_else(|_| C::default(), |i| terms[i].1))
}

pub(crate) fn persist(model: &mut Model, value: Vec<(Monomial, u128)>) -> Result<State> {
    let size = cells(&value);
    check_cells(size)?;
    model.allocate(size)?;
    Ok(value.into())
}

pub(crate) fn add(model: &mut Model, level: Level, left: &State, right: &State) -> Result<State> {
    let mut result = Accumulator::default();
    for (term, coefficient) in left.iter().chain(right.iter()) {
        result.add(model, term.clone(), *coefficient, level)?;
    }
    persist(model, result.finish())
}

pub(crate) fn compare(model: &mut Model, left: &State, right: &State) -> Result<Ordering> {
    // Enumeration increments the least Laurent monomial first.
    for (a, b) in left.iter().rev().zip(right.iter().rev()) {
        model.spend(a.0.cells())?;
        let cmp = a.cmp(b);
        if !cmp.is_eq() {
            return Ok(cmp);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

trait Ring<C: Copy>: Copy {
    fn one(self) -> C;
    fn add(self, a: C, b: C) -> Result<C>;
    fn mul(self, model: &mut Model, a: C, b: C) -> Result<C>;
}

#[derive(Clone, Copy)]
struct Exact;

impl Ring<i128> for Exact {
    fn one(self) -> i128 {
        1
    }
    fn add(self, a: i128, b: i128) -> Result<i128> {
        a.checked_add(b).ok_or_else(coefficient_overflow)
    }
    fn mul(self, _: &mut Model, a: i128, b: i128) -> Result<i128> {
        a.checked_mul(b).ok_or_else(coefficient_overflow)
    }
}

impl Ring<u128> for Level {
    fn one(self) -> u128 {
        1
    }
    fn add(self, a: u128, b: u128) -> Result<u128> {
        Ok(self.add(a, b))
    }
    fn mul(self, model: &mut Model, a: u128, b: u128) -> Result<u128> {
        if let Some(product) = a.checked_mul(b) {
            return Ok(product % self.modulus);
        }
        let mut value = a.max(b);
        let mut scalar = a.min(b);
        model.spend(1 + scalar.ilog2() as usize)?;
        let mut product = 0;
        while scalar != 0 {
            if scalar % 2 == 1 {
                product = self.add(product, value);
            }
            scalar /= 2;
            value = self.add(value, value);
        }
        Ok(product)
    }
}

#[derive(Default)]
struct Accumulator<C> {
    terms: BTreeMap<Monomial, C>,
    cells: usize,
}

impl<C: Copy + Default + Eq> Accumulator<C> {
    fn add(
        &mut self,
        model: &mut Model,
        term: Monomial,
        value: C,
        ring: impl Ring<C>,
    ) -> Result<()> {
        if value == C::default() {
            return Ok(());
        }
        let weight = term.cells();
        model.spend((1 + term.exponents.len()) * (1 + self.terms.len().max(1).ilog2() as usize))?;
        match self.terms.entry(term) {
            Entry::Vacant(entry) => {
                self.cells += weight;
                check_cells(self.cells)?;
                entry.insert(value);
            }
            Entry::Occupied(mut entry) => {
                let value = ring.add(*entry.get(), value)?;
                if value == C::default() {
                    self.cells -= weight;
                    entry.remove();
                } else {
                    *entry.get_mut() = value;
                }
            }
        }
        Ok(())
    }

    fn finish(self) -> Vec<(Monomial, C)> {
        self.terms.into_iter().collect()
    }
}

fn power<C: Copy + Default + Eq>(
    model: &mut Model,
    variables: usize,
    mut exponent: u128,
    mut base: Vec<(Monomial, C)>,
    ring: impl Ring<C>,
) -> Result<Vec<(Monomial, C)>> {
    let mut result = vec![(Monomial::unit(variables), ring.one())];
    // Small Laurent factors stay narrow. Multiplying by them directly avoids
    // a convolution between two wide intermediate powers.
    if exponent <= 32 && (2..=16).contains(&base.len()) {
        for _ in 0..exponent {
            model.spend(1)?;
            result = multiply(model, &result, &base, ring)?;
        }
        return Ok(result);
    }
    while exponent != 0 {
        model.spend(1)?;
        if exponent % 2 == 1 {
            result = multiply(model, &result, &base, ring)?;
        }
        exponent /= 2;
        if exponent != 0 {
            base = multiply(model, &base, &base, ring)?;
        }
    }
    Ok(result)
}

fn multiply<C: Copy + Default + Eq>(
    model: &mut Model,
    left: &[(Monomial, C)],
    right: &[(Monomial, C)],
    ring: impl Ring<C>,
) -> Result<Vec<(Monomial, C)>> {
    let mut result = Accumulator::default();
    for (a, x) in left {
        for (b, y) in right {
            model.spend(1)?;
            let coefficient = ring.mul(model, *x, *y)?;
            if coefficient != C::default() {
                result.add(model, a.product(b)?, coefficient, ring)?;
            }
        }
    }
    Ok(result.finish())
}

fn cells<C>(terms: &[(Monomial, C)]) -> usize {
    terms.iter().map(|(term, _)| term.cells()).sum()
}

fn check_cells(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        Err(Error::plain(
            "Laurent period exceeds its sparse storage limit",
        ))
    } else {
        Ok(())
    }
}

fn coefficient_overflow() -> Error {
    Error::plain("exact Laurent period coefficient exceeds i128")
}
fn exponent_overflow() -> Error {
    Error::plain("Laurent period exponent is out of range")
}
fn enumeration_limit() -> Error {
    Error::plain("Cartier state exceeds the finite type enumeration limit")
}
