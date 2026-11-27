use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::super::Field;
use super::matrix::Matrix;
use super::{coordinate_vector, digits, Class, Groupoid};

mod gauge;
mod series;
mod solve;
use series::{homogeneous, Monomial, Series, AXES};
use solve::Solver;

const MAX_ORDER: u8 = 6;
const LOOKAHEAD_LIMIT: usize = 50_000;
type Coefficients = BTreeMap<Monomial, Vec<u32>>;
pub(crate) type Terms = Vec<(Vec<u128>, Vec<u128>)>;

pub(crate) enum Expression {
    Constant(Label, bool),
    Start(Label, Vec<ParsedExpr>),
    Transform(Operation, Box<ParsedExpr>),
}

pub(crate) enum Operation {
    Lift,
    Correct(Terms),
    Gauge(Terms),
    Slice,
    First(u128),
}

pub(crate) enum Report {
    Declare(Label, Label, [u128; 2]),
    Frame(Label),
    Obstruction(Box<ParsedExpr>),
    Lookahead(Box<ParsedExpr>),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Jet {
    pub(crate) frame: Rc<str>,
    axes: usize,
    order: u8,
    coefficients: Coefficients,
}

impl Jet {
    pub(crate) fn render(&self) -> String {
        format!(
            "RVPartitionFormal({}) order={}, corrections={}",
            self.frame,
            self.order,
            display(self.axes, &self.coefficients)
        )
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(
            5 + self
                .coefficients
                .values()
                .map(|c| c.len() + AXES)
                .sum::<usize>(),
        )?;
        Ok(Value::RVPartitionFormal(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<std::cmp::Ordering> {
        model.spend(1 + (AXES + 9) * self.coefficients.len().min(rhs.coefficients.len()))?;
        Ok(self.cmp(rhs))
    }
}

pub(crate) struct Frame {
    name: Rc<str>,
    group: Rc<Groupoid>,
    axes: usize,
    order: u8,
    equation_width: usize,
    derivative: Solver,
    boundary: Solver,
}

struct Analysis {
    residual: Coefficients,
    particular: Coefficients,
    blocked: Coefficients,
}

impl Frame {
    pub(crate) fn new(
        model: &mut Model,
        name: Rc<str>,
        group: Rc<Groupoid>,
        options: [u128; 2],
    ) -> Result<Self> {
        let [axes, order] = options;
        if axes == 0 || axes > AXES as u128 || order == 0 || order > u128::from(MAX_ORDER) {
            return Err(Error::plain(
                "partition formal frames require 1..=4 variables and order 1..=6",
            ));
        }
        if group.quotient != group.submodule {
            return Err(Error::plain(
                "partition formal frame requires a self-extension space",
            ));
        }
        let hall = &group.hall;
        let field = &hall.field;
        let width = group
            .setup
            .relations(model, hall, &vec![0; group.setup.variables])?
            .len();
        let columns = (0..group.setup.variables)
            .map(|i| {
                group
                    .setup
                    .relations(model, hall, &coordinate_vector(group.setup.variables, i))
            })
            .collect::<Result<_>>()?;
        let derivative = Solver::new(model, field, width, columns)?;
        let columns = (0..group.setup.gauges)
            .map(|i| {
                group
                    .setup
                    .boundary(model, hall, &coordinate_vector(group.setup.gauges, i))
            })
            .collect::<Result<_>>()?;
        let boundary = Solver::new(model, field, group.setup.variables, columns)?;
        if derivative.rank != group.equations.rank || boundary.rank != group.boundaries.rank {
            return Err(Error::plain(
                "partition formal derivative disagrees with the self-extension chart",
            ));
        }
        model.allocate(
            16 + 3 * width * group.setup.variables + 3 * group.setup.variables * group.setup.gauges,
        )?;
        let frame = Self {
            name,
            group,
            axes: axes as usize,
            order: order as u8,
            equation_width: width,
            derivative,
            boundary,
        };
        frame.verify(model, &frame.constant(false))?;
        Ok(frame)
    }

    pub(crate) fn extension_name(&self) -> &Rc<str> {
        &self.group.name
    }

    pub(crate) fn arity(&self) -> usize {
        self.axes
    }

    pub(crate) fn frame_report(&self) -> String {
        let q = self.group.hall.field.order();
        let tangent = self.group.setup.variables - self.derivative.rank;
        format!("rvPartitionFormalFrame[{}] = q={q}, axes={}, order={}, arrowEntries={}, gaugeEntries={}, relationEntries={}, derivativeRank={}, tangent={tangent}, selfExt1={}",self.name,self.axes,self.order,self.group.setup.variables,self.group.setup.gauges,self.equation_width,self.derivative.rank,self.group.dimension())
    }

    pub(crate) fn constant(&self, full: bool) -> Jet {
        self.jet(if full { self.order } else { 0 }, Coefficients::new())
    }

    fn jet(&self, order: u8, coefficients: Coefficients) -> Jet {
        Jet {
            frame: self.name.clone(),
            axes: self.axes,
            order,
            coefficients,
        }
    }

    pub(crate) fn start(&self, model: &mut Model, directions: &[Rc<Class>]) -> Result<Jet> {
        if directions.len() != self.axes {
            return Err(Error::plain(
                "partition formal start requires one direction per variable",
            ));
        }
        let mut terms = Coefficients::new();
        for (i, direction) in directions.iter().enumerate() {
            if direction.space != self.group.name {
                return Err(Error::plain(
                    "partition formal direction has the wrong nominal self-extension space",
                ));
            }
            let lift = self.group.lift(model, direction)?;
            if lift.iter().any(|&x| x != 0) {
                terms.insert(Monomial::linear(i), lift);
            }
        }
        let jet = self.jet(1, terms);
        self.verify(model, &jet)?;
        Ok(jet)
    }

    pub(crate) fn first(&self, model: &mut Model, jet: &Jet, i: u128) -> Result<Class> {
        let axis = usize::try_from(i)
            .ok()
            .filter(|&i| i < self.axes)
            .ok_or_else(|| Error::plain("partition formal variable index is out of range"))?;
        let lift = jet
            .coefficients
            .get(&Monomial::linear(axis))
            .cloned()
            .unwrap_or_else(|| vec![0; self.group.setup.variables]);
        self.group.project(model, lift)
    }

    fn next(&self, jet: &Jet) -> Result<u8> {
        if jet.order == self.order {
            return Err(Error::plain(
                "partition formal jet is already at the frame order",
            ));
        }
        Ok(jet.order + 1)
    }

    fn verify(&self, model: &mut Model, jet: &Jet) -> Result<()> {
        let arrows = self.arrows(model, jet)?;
        let values = self.relations(model, &arrows, jet.order)?;
        if let Some((power, values)) = values.first_key_value() {
            return Err(Error::plain(format!(
                "partition formal relation is nonzero at {:?}: {values:?}",
                &power.0[..self.axes]
            )));
        }
        Ok(())
    }

    fn arrows(&self, model: &mut Model, jet: &Jet) -> Result<Vec<Option<Series>>> {
        let field = &self.group.hall.field;
        let mut arrows = Vec::new();
        for (source, block) in self.group.setup.m.arrows.iter().enumerate() {
            let Some(block) = block else {
                arrows.push(None);
                continue;
            };
            let range = self.group.setup.arrows[source].as_ref().unwrap();
            let mut path = Series::constant(block.clone());
            for (&monomial, entries) in &jet.coefficients {
                let coefficients =
                    Matrix::new(block.rows, block.cols, entries[range.clone()].to_vec());
                path.add_term(model, field, monomial, &coefficients, 1)?;
            }
            arrows.push(Some(path));
        }
        Ok(arrows)
    }

    fn relations(
        &self,
        model: &mut Model,
        arrows: &[Option<Series>],
        order: u8,
    ) -> Result<Coefficients> {
        let hall = &self.group.hall;
        let field = &hall.field;
        let mut coefficients = Coefficients::new();
        let mut offset = 0;
        for (start, heads) in hall.path_heads.iter().enumerate() {
            let mut at = start;
            let mut path = Series::identity(self.group.setup.m.dims[start]);
            let mut exists = true;
            for _ in 0..heads.len() {
                let Some(next) = hall.successor[at] else {
                    exists = false;
                    break;
                };
                path = arrows[at]
                    .as_ref()
                    .unwrap()
                    .compose(model, field, &path, order)?;
                at = next;
            }
            if !exists {
                continue;
            }
            for (&monomial, matrix) in &path.terms {
                let entries = coefficients
                    .entry(monomial)
                    .or_insert_with(|| vec![0; self.equation_width]);
                for i in 0..matrix.rows {
                    for j in 0..matrix.cols {
                        entries[offset + i * matrix.cols + j] = matrix.get(i, j);
                    }
                }
            }
            offset += path.rows * path.cols;
        }
        if offset != self.equation_width {
            return Err(Error::plain(
                "partition formal relations do not match the cocycle derivative axes",
            ));
        }
        Ok(coefficients)
    }

    fn analysis(&self, model: &mut Model, jet: &Jet) -> Result<Analysis> {
        let next = self.next(jet)?;
        let arrows = self.arrows(model, jet)?;
        let residual = self.relations(model, &arrows, next)?;
        let field = &self.group.hall.field;
        let mut particular = Coefficients::new();
        let mut blocked = Coefficients::new();
        for (&monomial, value) in &residual {
            if monomial.degree() != next {
                return Err(Error::plain(
                    "partition formal jet has a residual below its certified order",
                ));
            }
            let rhs: Vec<_> = value.iter().map(|&x| field.sub(0, x)).collect();
            let (solution, reduced) = self.derivative.solve(model, field, &rhs)?;
            if reduced.iter().any(|&x| x != 0) {
                blocked.insert(monomial, reduced);
            } else if solution.iter().any(|&x| x != 0) {
                particular.insert(monomial, solution);
            }
        }
        Ok(Analysis {
            residual,
            particular,
            blocked,
        })
    }

    fn extend(&self, model: &mut Model, jet: &Jet, values: Coefficients) -> Result<Jet> {
        let mut terms = jet.coefficients.clone();
        terms.extend(values);
        let result = self.jet(self.next(jet)?, terms);
        self.verify(model, &result)?;
        Ok(result)
    }

    pub(crate) fn lift(&self, model: &mut Model, jet: &Jet) -> Result<Jet> {
        let analysis = self.analysis(model, jet)?;
        if let Some((power, _)) = analysis.blocked.first_key_value() {
            return Err(Error::plain(format!(
                "partition formal lift is obstructed at {:?}",
                &power.0[..self.axes]
            )));
        }
        self.extend(model, jet, analysis.particular)
    }

    pub(crate) fn correct(&self, model: &mut Model, jet: &Jet, terms: Terms) -> Result<Jet> {
        let order = self.next(jet)?;
        let terms = self.read_terms(terms, order..=order, self.group.setup.variables)?;
        self.extend(model, jet, terms)
    }

    pub(crate) fn obstruction_report(&self, model: &mut Model, jet: &Jet) -> Result<String> {
        let next = self.next(jet)?;
        let analysis = self.analysis(model, jet)?;
        let field = &self.group.hall.field;
        let mut witnesses = Vec::new();
        for (&power, reduced) in &analysis.blocked {
            let rhs: Vec<_> = analysis.residual[&power]
                .iter()
                .map(|&x| field.sub(0, x))
                .collect();
            let (row, value) = self.derivative.obstruction(model, field, &rhs, reduced)?;
            witnesses.push((power.0[..self.axes].to_vec(), row, field.sub(0, value)));
        }
        let blocks = homogeneous(self.axes, next).len();
        let q = field.order();
        let (raw, slice) = if witnesses.is_empty() {
            (
                power_text(
                    q,
                    blocks * (self.group.setup.variables - self.derivative.rank),
                ),
                power_text(q, blocks * self.group.dimension()),
            )
        } else {
            ("0".into(), "0".into())
        };
        Ok(format!("rvPartitionFormalObstruction[{}] = nextOrder={next}, residual={}, obstructions={witnesses:?}, lifts={raw}, liftsModuloNewGauge={slice}",self.name,display(self.axes,&analysis.residual)))
    }

    pub(crate) fn lookahead(&self, model: &mut Model, jet: &Jet) -> Result<String> {
        if jet.order + 2 > self.order {
            return Err(Error::plain(
                "partition formal lookahead needs two remaining orders",
            ));
        }
        let analysis = self.analysis(model, jet)?;
        let next = self.next(jet)?;
        let field = &self.group.hall.field;
        let q = field.order();
        if !analysis.blocked.is_empty() {
            return Ok(format!(
                "rvPartitionFormalLookahead[{}] = order={next}, lifts=0, extendable=0",
                self.name
            ));
        }
        let kernel = self.group.equations.kernel(model, field)?;
        let powers = homogeneous(self.axes, next);
        let dimension = kernel.len() * powers.len();
        let count = (q as usize)
            .checked_pow(dimension as u32)
            .filter(|&n| n <= LOOKAHEAD_LIMIT)
            .ok_or_else(|| {
                Error::plain("partition formal lookahead exceeds 50,000 next-order lifts")
            })?;
        let mut extendable = 0;
        for i in 0..count {
            let mut choices = analysis.particular.clone();
            let values = digits(i, q, dimension);
            model.spend(1 + dimension)?;
            if !kernel.is_empty() {
                for (&power, weights) in powers.iter().zip(values.chunks(kernel.len())) {
                    let mut vector = choices
                        .remove(&power)
                        .unwrap_or_else(|| vec![0; self.group.setup.variables]);
                    for (row, &scale) in kernel.iter().zip(weights) {
                        if scale == 0 {
                            continue;
                        }
                        for (x, &y) in vector.iter_mut().zip(row) {
                            *x = field.add(*x, field.mul(y, scale));
                        }
                    }
                    if vector.iter().any(|&x| x != 0) {
                        choices.insert(power, vector);
                    }
                }
            }
            let trial = self.extend(model, jet, choices)?;
            if self.analysis(model, &trial)?.blocked.is_empty() {
                extendable += 1;
            }
        }
        Ok(format!(
            "rvPartitionFormalLookahead[{}] = order={next}, lifts={count}, extendable={extendable}",
            self.name
        ))
    }

    fn read_terms(
        &self,
        terms: Terms,
        degrees: std::ops::RangeInclusive<u8>,
        width: usize,
    ) -> Result<Coefficients> {
        let mut output = Coefficients::new();
        for (powers, values) in terms {
            if powers.len() != self.axes || powers.iter().any(|&x| x > u128::from(MAX_ORDER)) {
                return Err(Error::plain("partition formal monomial has the wrong number of variables or exceeds the order bound"));
            }
            let mut exponents = [0; AXES];
            for (i, p) in powers.into_iter().enumerate() {
                exponents[i] = p as u8;
            }
            let monomial = Monomial(exponents);
            if !degrees.contains(&monomial.degree()) {
                return Err(Error::plain(
                    "partition formal term has an invalid total degree for this operation",
                ));
            }
            let value = self.group.valid_codes(values, width)?;
            if output.insert(monomial, value).is_some() {
                return Err(Error::plain("partition formal input repeats a monomial"));
            }
        }
        output.retain(|_, v| v.iter().any(|&x| x != 0));
        Ok(output)
    }
}

fn display(axes: usize, terms: &Coefficients) -> String {
    let list: Vec<_> = terms
        .iter()
        .map(|(power, value)| (&power.0[..axes], value))
        .collect();
    format!("{list:?}")
}

fn power_text(q: u32, degree: usize) -> String {
    u128::from(q)
        .checked_pow(degree as u32)
        .map_or_else(|| format!("{q}^{degree}"), |n| n.to_string())
}
