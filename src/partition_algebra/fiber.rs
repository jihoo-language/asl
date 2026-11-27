use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::{Model, Value};
use crate::syntax::ParsedExpr;

use super::local::LocalAlgebra;
use super::{diagram, evaluate_polynomial, Algebra, Diagram, Integer};

mod elimination;
pub(crate) mod gabriel;
mod loewy;

pub(crate) const MAX_RANK: usize = 3;
const FULL_LOEWY_RANK: usize = 2;
type Coordinates = BTreeMap<usize, Rational>;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Parameter {
    pub(crate) rank: usize,
    pub(crate) n: u128,
}

impl Parameter {
    pub(crate) fn new(rank: usize, n: u128) -> Result<Self> {
        if rank > MAX_RANK {
            return Err(Error::plain(
                "specialized partition algebras are limited to rank 3",
            ));
        }
        Ok(Self { rank, n })
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) enum Mode {
    Full,
    Semisimple,
    TensorImage,
}

impl Mode {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Full => "RVPartitionFiber",
            Self::Semisimple => "RVPartitionSemisimple",
            Self::TensorImage => "RVPartitionTensorImage",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Grade {
    pub(crate) parameter: Parameter,
    pub(crate) mode: Mode,
}

impl Grade {
    pub(crate) fn is_zero_ring(self) -> bool {
        self.mode == Mode::TensorImage && self.parameter.n == 0 && self.parameter.rank != 0
    }

    pub(crate) fn render(self) -> String {
        format!(
            "{}({},{})",
            self.mode.name(),
            self.parameter.rank,
            self.parameter.n
        )
    }
}

pub(crate) enum Expression {
    Specialize(u128, Box<ParsedExpr>),
    Project(Mode, Box<ParsedExpr>),
    Scale(i128, i128, Box<ParsedExpr>),
    RadicalBasis(Parameter, u128, u128),
    NegligibleBasis(Parameter, u128),
    GabrielGenerator(Parameter, gabriel::Generator),
}

pub(crate) enum Report {
    Loewy(Parameter),
    Gabriel(Parameter, bool),
    Witness(Parameter, u128, u128),
    Coefficient(Vec<u128>, Box<ParsedExpr>),
    Contains(Option<u128>, Box<ParsedExpr>),
    Pair(bool, Box<[ParsedExpr; 2]>),
}

pub(crate) struct Element {
    pub(crate) grade: Grade,
    diagrams: Rc<[Diagram]>,
    coefficients: Coordinates,
}

impl Element {
    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(4 + 7 * self.coefficients.len())?;
        Ok(Value::RVPartitionFiber(Rc::new(self)))
    }

    pub(crate) fn render(&self) -> String {
        let entries: Vec<_> = self
            .coefficients
            .iter()
            .map(|(&i, coefficient)| format!("({:?},{coefficient})", self.diagrams[i].labels()))
            .collect();
        format!("{}{{{}}}", self.grade.render(), entries.join(", "))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(1 + 4 * self.coefficients.len().min(rhs.coefficients.len()))?;
        Ok(self
            .grade
            .cmp(&rhs.grade)
            .then_with(|| self.coefficients.cmp(&rhs.coefficients)))
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let mut sum = self.coefficients.clone();
        model.spend(self.coefficients.len() + rhs.coefficients.len())?;
        add_multiple(&mut sum, &rhs.coefficients, Rational::integer(1))?;
        Ok(Self {
            grade: self.grade,
            diagrams: self.diagrams.clone(),
            coefficients: sum,
        })
    }

    pub(crate) fn scale(&self, model: &mut Model, n: i128, d: i128) -> Result<Self> {
        let factor = ratio(n, d)?;
        let mut coefficients = Coordinates::new();
        model.spend(self.coefficients.len())?;
        add_multiple(&mut coefficients, &self.coefficients, factor)?;
        Ok(Self {
            grade: self.grade,
            diagrams: self.diagrams.clone(),
            coefficients,
        })
    }
}

#[derive(Clone, Copy)]
struct Product {
    index: u8,
    loops: u8,
}

pub(crate) struct DiagramTable {
    rank: usize,
    diagrams: Rc<[Diagram]>,
    index: BTreeMap<Diagram, usize>,
    products: Vec<Vec<Product>>,
}

impl DiagramTable {
    pub(crate) fn new(model: &mut Model, rank: usize) -> Result<Self> {
        let mut diagrams = Vec::new();
        diagram::all_patterns(model, 2 * rank, (2 * rank) as u128, &mut |_, pattern| {
            diagrams.push(Diagram::literal(rank, pattern.to_vec())?);
            Ok(true)
        })?;
        let size = diagrams.len();
        let index: BTreeMap<_, _> = diagrams.iter().cloned().zip(0..size).collect();
        let mut products = Vec::with_capacity(size);
        for top in &diagrams {
            model.spend(size * (1 + 8 * rank))?;
            products.push(
                diagrams
                    .iter()
                    .map(|bottom| {
                        let (diagram, loops) = top.compose(bottom, rank);
                        Product {
                            index: index[&diagram] as u8,
                            loops: loops as u8,
                        }
                    })
                    .collect(),
            );
        }
        model.allocate(2 * size * size + 2 * size * (2 * rank + 2))?;
        Ok(Self {
            rank,
            diagrams: diagrams.into(),
            index,
            products,
        })
    }

    fn orbit(&self, model: &mut Model, i: usize) -> Result<Coordinates> {
        Ok(self.diagrams[i]
            .orbit(model)?
            .into_iter()
            .map(|((diagram, _), mu)| (self.index[&diagram], Rational::integer(mu)))
            .collect())
    }

    fn diagram_index(&self, labels: Vec<u128>) -> Result<usize> {
        Ok(self.index[&Diagram::literal(self.rank, labels)?])
    }
}

pub(crate) struct Fiber {
    pub(crate) parameter: Parameter,
    table: Rc<DiagramTable>,
    traces: RefCell<Option<Rc<[i128]>>>,
    kernel: RefCell<Option<Rc<elimination::Kernel>>>,
    powers: RefCell<Option<Rc<loewy::Powers>>>,
    gabriel: RefCell<Option<Rc<gabriel::Gabriel>>>,
}

impl Fiber {
    pub(crate) fn new(parameter: Parameter, table: Rc<DiagramTable>) -> Self {
        Self {
            parameter,
            table,
            traces: RefCell::default(),
            kernel: RefCell::default(),
            powers: RefCell::default(),
            gabriel: RefCell::default(),
        }
    }

    fn element(&self, mode: Mode, coefficients: Coordinates) -> Element {
        Element {
            grade: Grade {
                parameter: self.parameter,
                mode,
            },
            diagrams: self.table.diagrams.clone(),
            coefficients,
        }
    }

    pub(crate) fn zero(&self, mode: Mode) -> Element {
        self.element(mode, Coordinates::new())
    }

    fn size(&self) -> usize {
        self.table.diagrams.len()
    }

    fn loop_factor(&self, loops: u8) -> Result<i128> {
        self.parameter
            .n
            .checked_pow(u32::from(loops))
            .and_then(|n| i128::try_from(n).ok())
            .ok_or_else(arithmetic)
    }

    fn factors(&self) -> Result<Vec<i128>> {
        (0..=self.parameter.rank)
            .map(|loops| self.loop_factor(loops as u8))
            .collect()
    }

    pub(crate) fn specialize_integral(
        &self,
        model: &mut Model,
        algebra: &Algebra,
    ) -> Result<Element> {
        self.specialize(model, algebra, Integer::from_i128(1))
    }

    pub(crate) fn specialize_local(
        &self,
        model: &mut Model,
        algebra: &LocalAlgebra,
    ) -> Result<Element> {
        let (numerator, denominator) = algebra.specialization(model, self.parameter.n)?;
        self.specialize(model, numerator, denominator)
    }

    fn specialize(
        &self,
        model: &mut Model,
        algebra: &Algebra,
        denominator: Integer,
    ) -> Result<Element> {
        let mut coefficients = Coordinates::new();
        for (diagram, row) in algebra.diagram_groups(model)? {
            let polynomial = row
                .into_iter()
                .map(|(degree, value)| (degree, Integer::from_i128(value)))
                .collect();
            let value = evaluate_polynomial(model, polynomial, self.parameter.n)?;
            if !value.is_zero() {
                coefficients.insert(
                    self.table.index[diagram],
                    from_wide(model, value, denominator)?,
                );
            }
        }
        Ok(self.element(Mode::Full, coefficients))
    }

    pub(crate) fn compose(&self, model: &mut Model, a: &Element, b: &Element) -> Result<Element> {
        let mode = a.grade.mode;
        if a.coefficients.is_empty() || b.coefficients.is_empty() {
            return Ok(self.zero(mode));
        }
        let coefficients = if mode == Mode::TensorImage {
            let a = self.lift_orbit(model, &a.coefficients)?;
            let b = self.lift_orbit(model, &b.coefficients)?;
            self.multiply(model, &a, &b)?
        } else {
            self.multiply(model, &a.coefficients, &b.coefficients)?
        };
        self.project_coordinates(model, mode, coefficients)
    }

    fn multiply(&self, model: &mut Model, a: &Coordinates, b: &Coordinates) -> Result<Coordinates> {
        let mut result = Coordinates::new();
        if a.is_empty() || b.is_empty() {
            return Ok(result);
        }
        let factors = self.factors()?;
        model.spend(a.len() * b.len())?;
        for (&i, &x) in a {
            for (&j, &y) in b {
                let Product { index, loops } = self.table.products[i][j];
                let scale = Rational::integer(factors[usize::from(loops)]);
                if scale.numerator != 0 {
                    accumulate(
                        &mut result,
                        usize::from(index),
                        multiply(multiply(x, y)?, scale)?,
                    )?;
                }
            }
        }
        Ok(result)
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        mode: Mode,
        value: &Element,
    ) -> Result<Element> {
        self.project_coordinates(model, mode, value.coefficients.clone())
    }

    fn project_coordinates(
        &self,
        model: &mut Model,
        mode: Mode,
        mut coefficients: Coordinates,
    ) -> Result<Element> {
        if coefficients.is_empty() {
            return Ok(self.zero(mode));
        }
        match mode {
            Mode::Full => (),
            Mode::Semisimple => {
                let kernel = self.radical(model)?;
                for (&i, row) in &kernel.rows {
                    if let Some(&value) = coefficients.get(&i) {
                        model.spend(row.len())?;
                        add_multiple(&mut coefficients, row, negate(value)?)?;
                    }
                }
            }
            Mode::TensorImage => coefficients = self.orbit_image(model, &coefficients)?,
        }
        Ok(self.element(mode, coefficients))
    }

    fn lift_orbit(&self, model: &mut Model, coefficients: &Coordinates) -> Result<Coordinates> {
        let mut lifted = Coordinates::new();
        for (&index, &coefficient) in coefficients {
            add_multiple(&mut lifted, &self.table.orbit(model, index)?, coefficient)?;
        }
        Ok(lifted)
    }

    fn orbit_image(&self, model: &mut Model, coefficients: &Coordinates) -> Result<Coordinates> {
        let mut image = Coordinates::new();
        for (i, orbit) in self.table.diagrams.iter().enumerate() {
            if orbit.blocks() as u128 > self.parameter.n {
                continue;
            }
            let labels: Vec<_> = orbit.labels().iter().copied().map(u128::from).collect();
            let mut coefficient = Rational::ZERO;
            model.spend(coefficients.len() * (1 + 2 * self.parameter.rank))?;
            for (&j, &value) in coefficients {
                if self.table.diagrams[j].supports(&labels) {
                    coefficient = add(coefficient, value)?;
                }
            }
            accumulate(&mut image, i, coefficient)?;
        }
        Ok(image)
    }

    pub(crate) fn coefficient(&self, value: &Element, labels: Vec<u128>) -> Result<String> {
        let index = self.table.diagram_index(labels)?;
        let coefficient = value
            .coefficients
            .get(&index)
            .copied()
            .unwrap_or(Rational::ZERO);
        Ok(format!(
            "rvPartitionFiberCoefficient[{:?}] = {coefficient}",
            self.table.diagrams[index].labels()
        ))
    }

    pub(crate) fn pair(
        &self,
        model: &mut Model,
        regular: bool,
        a: &Element,
        b: &Element,
    ) -> Result<String> {
        let value = if regular {
            self.regular_pair(model, &a.coefficients, &b.coefficients)?
        } else {
            let a = self.orbit_image(model, &a.coefficients)?;
            let b = self.orbit_image(model, &b.coefficients)?;
            let mut norm = Rational::ZERO;
            for (&j, &x) in &a {
                if let Some(&y) = b.get(&j) {
                    let blocks = self.table.diagrams[j].blocks() as u128;
                    let weight = (0..blocks).try_fold(1i128, |weight, k| {
                        let factor =
                            i128::try_from(self.parameter.n - k).map_err(|_| arithmetic())?;
                        checked_mul(weight, factor)
                    })?;
                    norm = add(norm, multiply(multiply(x, y)?, Rational::integer(weight))?)?;
                }
            }
            norm
        };
        let label = if regular {
            "rvPartitionRegularPair"
        } else {
            "rvPartitionCategoricalPair"
        };
        Ok(format!("{label} = {value}"))
    }

    fn regular_pair(
        &self,
        model: &mut Model,
        a: &Coordinates,
        b: &Coordinates,
    ) -> Result<Rational> {
        let traces = self.regular_traces(model)?;
        let factors = self.factors()?;
        let mut value = Rational::ZERO;
        model.spend(a.len() * b.len())?;
        for (&i, &x) in a {
            for (&j, &y) in b {
                let Product { index, loops } = self.table.products[i][j];
                let t = checked_mul(factors[usize::from(loops)], traces[usize::from(index)])?;
                value = add(value, multiply(multiply(x, y)?, Rational::integer(t))?)?;
            }
        }
        Ok(value)
    }

    fn regular_traces(&self, model: &mut Model) -> Result<Rc<[i128]>> {
        if let Some(cached) = self.traces.borrow().as_ref() {
            return Ok(cached.clone());
        }
        let factors = self.factors()?;
        let mut traces = Vec::with_capacity(self.size());
        for row in &self.table.products {
            model.spend(row.len())?;
            let mut trace = 0i128;
            for (j, &Product { index, loops }) in row.iter().enumerate() {
                if usize::from(index) == j {
                    trace = trace
                        .checked_add(factors[usize::from(loops)])
                        .ok_or_else(arithmetic)?;
                }
            }
            traces.push(trace);
        }
        model.allocate(traces.len() * 2)?;
        let traces: Rc<[i128]> = traces.into();
        *self.traces.borrow_mut() = Some(traces.clone());
        Ok(traces)
    }

    fn radical(&self, model: &mut Model) -> Result<Rc<elimination::Kernel>> {
        if let Some(kernel) = self.kernel.borrow().as_ref() {
            return Ok(kernel.clone());
        }
        let factors = self.factors()?;
        let traces = self.regular_traces(model)?;
        model.spend(self.size() * self.size())?;
        let gram = self
            .table
            .products
            .iter()
            .map(|row| {
                row.iter()
                    .map(|&Product { index, loops }| {
                        checked_mul(factors[usize::from(loops)], traces[usize::from(index)])
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let kernel = Rc::new(elimination::certify(model, &gram)?);
        model.allocate(
            2 * kernel.rows.len() + 6 * kernel.rows.values().map(BTreeMap::len).sum::<usize>(),
        )?;
        *self.kernel.borrow_mut() = Some(kernel.clone());
        Ok(kernel)
    }

    fn powers(&self, model: &mut Model) -> Result<Rc<loewy::Powers>> {
        if let Some(powers) = self.powers.borrow().as_ref() {
            return Ok(powers.clone());
        }
        let radical = self.radical(model)?;
        let powers = Rc::new(loewy::Powers::new(model, self, &radical)?);
        *self.powers.borrow_mut() = Some(powers.clone());
        Ok(powers)
    }

    fn gabriel(&self, model: &mut Model) -> Result<Rc<gabriel::Gabriel>> {
        if let Some(data) = self.gabriel.borrow().as_ref() {
            return Ok(data.clone());
        }
        let data = Rc::new(gabriel::Gabriel::new(model, self)?);
        *self.gabriel.borrow_mut() = Some(data.clone());
        Ok(data)
    }

    pub(crate) fn gabriel_generator(
        &self,
        model: &mut Model,
        generator: gabriel::Generator,
    ) -> Result<Element> {
        self.gabriel(model)?.generator(self, generator)
    }

    pub(crate) fn gabriel_report(&self, model: &mut Model, relations: bool) -> Result<String> {
        let data = self.gabriel(model)?;
        Ok(if relations {
            data.relations()
        } else {
            data.data()
        })
    }

    fn negligible_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.table
            .diagrams
            .iter()
            .enumerate()
            .filter_map(|(i, diagram)| (diagram.blocks() as u128 > self.parameter.n).then_some(i))
    }

    pub(crate) fn negligible_basis(&self, model: &mut Model, index: u128) -> Result<Element> {
        let i = usize::try_from(index)
            .ok()
            .and_then(|j| self.negligible_indices().nth(j));
        let i =
            i.ok_or_else(|| Error::plain("partition negligible basis index is out of range"))?;
        Ok(self.element(Mode::Full, self.table.orbit(model, i)?))
    }

    pub(crate) fn radical_basis(
        &self,
        model: &mut Model,
        power: u128,
        index: u128,
    ) -> Result<Element> {
        let range = || Error::plain("partition radical basis index is out of range");
        let i = usize::try_from(index).map_err(|_| range())?;
        let coefficients = match power {
            0 if i < self.size() => BTreeMap::from([(i, Rational::integer(1))]),
            0 => return Err(range()),
            1 => self
                .radical(model)?
                .rows
                .values()
                .nth(i)
                .cloned()
                .ok_or_else(range)?,
            _ => self
                .powers(model)?
                .basis(power)
                .get(i)
                .cloned()
                .ok_or_else(range)?,
        };
        Ok(self.element(Mode::Full, coefficients))
    }

    pub(crate) fn contains(
        &self,
        model: &mut Model,
        power: Option<u128>,
        value: &Element,
    ) -> Result<String> {
        if let Some(power) = power {
            let belongs = match power {
                0 => true,
                1 => self
                    .project(model, Mode::Semisimple, value)?
                    .coefficients
                    .is_empty(),
                _ if value.coefficients.is_empty() => true,
                _ => self
                    .powers(model)?
                    .contains(model, power, &value.coefficients)?,
            };
            Ok(format!("rvPartitionRadicalContains[{power}] = {belongs}"))
        } else {
            let belongs = self.orbit_image(model, &value.coefficients)?.is_empty();
            Ok(format!("rvPartitionNegligibleContains = {belongs}"))
        }
    }

    pub(crate) fn loewy_report(&self, model: &mut Model) -> Result<String> {
        let radical = self.radical(model)?.rows.len();
        let dimension = self.size();
        let negligible = self.negligible_indices().count();
        let powers = if self.parameter.rank <= FULL_LOEWY_RANK || radical == 0 {
            let dims = self.powers(model)?.dimensions(dimension, radical);
            format!("powers={dims:?}")
        } else {
            format!("powers=[{dimension}, {radical}] (through J^1)")
        };
        Ok(format!("rvPartitionLoewy[{},{}] = dimension={dimension}, negligible={negligible}, radical={radical}, semisimple={}, tensorImage={}, {powers}",
            self.parameter.rank, self.parameter.n, dimension - radical, dimension - negligible))
    }

    pub(crate) fn witness(&self, model: &mut Model, power: u128, index: u128) -> Result<String> {
        let (a, b) = self.powers(model)?.witness(power, index)?;
        Ok(format!(
            "rvPartitionLoewyWitness[{},{},{power},{index}] = [({},{}), (1,{b})]",
            self.parameter.rank,
            self.parameter.n,
            power - 1,
            a
        ))
    }
}

fn arithmetic() -> Error {
    Error::plain(
        "specialized partition rational arithmetic exceeds the signed 128-bit coefficient range",
    )
}
fn checked_mul(a: i128, b: i128) -> Result<i128> {
    a.checked_mul(b).ok_or_else(arithmetic)
}
fn ratio(n: i128, d: i128) -> Result<Rational> {
    Rational::ratio(n, d).map_err(|_| arithmetic())
}
fn from_wide(model: &mut Model, value: Integer, denominator: Integer) -> Result<Rational> {
    let divisor = value.magnitude.gcd(model, denominator.magnitude)?;
    let numerator = value
        .divide_exact(
            model,
            Integer {
                magnitude: divisor,
                negative: denominator.negative,
            },
        )?
        .and_then(Integer::signed)
        .ok_or_else(arithmetic)?;
    let denominator = denominator
        .magnitude
        .divide_exact(model, divisor)?
        .and_then(|value| value.small())
        .and_then(|value| i128::try_from(value).ok())
        .ok_or_else(arithmetic)?;
    ratio(numerator, denominator)
}
fn add(a: Rational, b: Rational) -> Result<Rational> {
    a.add(b).map_err(|_| arithmetic())
}
fn multiply(a: Rational, b: Rational) -> Result<Rational> {
    a.multiply(b).map_err(|_| arithmetic())
}
fn negate(value: Rational) -> Result<Rational> {
    ratio(
        value.numerator.checked_neg().ok_or_else(arithmetic)?,
        value.denominator,
    )
}
fn accumulate(target: &mut Coordinates, i: usize, coefficient: Rational) -> Result<()> {
    if coefficient.numerator == 0 {
        return Ok(());
    }
    let sum = add(
        target.get(&i).copied().unwrap_or(Rational::ZERO),
        coefficient,
    )?;
    if sum.numerator == 0 {
        target.remove(&i);
    } else {
        target.insert(i, sum);
    }
    Ok(())
}
fn add_multiple(target: &mut Coordinates, source: &Coordinates, factor: Rational) -> Result<()> {
    if factor.numerator != 0 {
        for (&i, &value) in source {
            accumulate(target, i, multiply(value, factor)?)?;
        }
    }
    Ok(())
}
