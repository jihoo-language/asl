use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::syntax::ParsedExpr;

use super::{Nakayama, Parameter, Path, State};

mod field;
pub(crate) mod stability;
mod subspaces;
use field::Field;
pub(crate) use field::{FieldName, FieldSpec};

const MAX_MODULE_LENGTH: usize = 64;
const MAX_HALL_LENGTH: usize = 6;
type Dims = Vec<usize>;
type Part = (u8, u8, u16);
type Distribution = BTreeMap<(Shape, Shape), u128>;

#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Key {
    pub(crate) parameter: Parameter,
    pub(crate) field: FieldName,
}

pub(crate) struct KeySpec {
    pub(crate) parameter: Parameter,
    pub(crate) field: FieldSpec,
}

impl KeySpec {
    pub(crate) fn resolve(&self, model: &mut Model) -> Result<Key> {
        Ok(Key {
            parameter: self.parameter,
            field: self.field.resolve(model)?,
        })
    }
}

#[derive(Clone, Copy, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) enum Mode {
    Module,
    Algebra,
}

#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Grade {
    pub(crate) key: Key,
    pub(crate) mode: Mode,
}

impl Grade {
    pub(crate) fn render(&self) -> String {
        let Parameter { rank, n } = self.key.parameter;
        let field = self.key.field.render();
        let name = if self.mode == Mode::Module {
            "RVPartitionHallModule"
        } else {
            "RVPartitionHall"
        };
        format!("{name}({rank},{n},{field})")
    }
}

#[derive(Clone, Default, Eq, Ord, PartialEq, PartialOrd)]
struct Shape(Vec<Part>);

impl Shape {
    fn length(&self) -> usize {
        self.0
            .iter()
            .map(|&(_, l, m)| usize::from(l) * usize::from(m))
            .sum()
    }
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
enum Data {
    Module(Shape),
    Algebra(BTreeMap<Shape, i128>),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Record {
    pub(crate) grade: Grade,
    data: Data,
}

impl Record {
    pub(crate) fn zero(grade: Grade) -> Self {
        let data = match grade.mode {
            Mode::Module => Data::Module(Shape::default()),
            Mode::Algebra => Data::Algebra(BTreeMap::new()),
        };
        Self { grade, data }
    }

    fn module(&self) -> &Shape {
        let Data::Module(shape) = &self.data else {
            unreachable!()
        };
        shape
    }

    fn terms(&self) -> &BTreeMap<Shape, i128> {
        let Data::Algebra(terms) = &self.data else {
            unreachable!()
        };
        terms
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        let cells = match &self.data {
            Data::Module(shape) => 4 + 3 * shape.0.len(),
            Data::Algebra(terms) => 4 + terms.keys().map(|s| 4 + 3 * s.0.len()).sum::<usize>(),
        };
        model.allocate(cells)?;
        Ok(Value::RVPartitionHall(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        let cells = match (&self.data, &rhs.data) {
            (Data::Module(a), Data::Module(b)) => a.0.len().min(b.0.len()),
            (Data::Algebra(a), Data::Algebra(b)) => a.len().min(b.len()) * 12,
            _ => 1,
        };
        model.spend(cells + 1)?;
        Ok(self.cmp(rhs))
    }

    pub(crate) fn render(&self) -> String {
        match &self.data {
            Data::Module(shape) => format!("{}{:?}", self.grade.render(), shape.0),
            Data::Algebra(terms) => {
                let entries = terms
                    .iter()
                    .map(|(x, c)| format!("({:?},{c})", x.0))
                    .collect::<Vec<_>>()
                    .join(", ");
                format!("{}{{{entries}}}", self.grade.render())
            }
        }
    }

    pub(crate) fn as_basis(&self) -> Self {
        Self {
            grade: Grade {
                key: self.grade.key.clone(),
                mode: Mode::Algebra,
            },
            data: Data::Algebra(BTreeMap::from([(self.module().clone(), 1)])),
        }
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let mut out = self.terms().clone();
        model.spend(self.terms().len() + rhs.terms().len())?;
        for (u, &c) in rhs.terms() {
            accumulate(&mut out, u.clone(), c)?;
        }
        Ok(Self {
            grade: self.grade.clone(),
            data: Data::Algebra(out),
        })
    }

    pub(crate) fn neg(&self, model: &mut Model) -> Result<Self> {
        model.spend(self.terms().len())?;
        let terms = self
            .terms()
            .iter()
            .map(|(u, &c)| Ok((u.clone(), c.checked_neg().ok_or_else(arithmetic)?)))
            .collect::<Result<_>>()?;
        Ok(Self {
            grade: self.grade.clone(),
            data: Data::Algebra(terms),
        })
    }

    pub(crate) fn read_term(&self, module: &Self) -> String {
        format!(
            "rvPartitionHallTerm = {}",
            self.terms().get(module.module()).copied().unwrap_or(0)
        )
    }
}

pub(crate) enum Expression {
    Module(KeySpec, Vec<(u128, u128, u128)>),
    Unit(KeySpec),
    Basis(Box<ParsedExpr>),
    Product(Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Copy)]
pub(crate) enum StabilityQuery {
    Stability,
    Hn,
    HnCheck,
}

pub(crate) enum Report {
    Coefficient(Box<[ParsedExpr; 3]>),
    Flag(Box<[ParsedExpr; 4]>),
    Term(Box<[ParsedExpr; 2]>),
    Stability(StabilityQuery, Vec<i128>, Box<ParsedExpr>),
    Chamber(Parameter, u128),
    Wall(Parameter, [u128; 2]),
}

type Cache = BTreeMap<(Shape, Dims), Rc<Distribution>>;

pub(crate) struct Hall {
    key: Key,
    field: Field,
    vertices: usize,
    shapes: Vec<State>,
    dimensions: Vec<Dims>,
    path_heads: Vec<Vec<usize>>,
    successor: Vec<Option<usize>>,
    predecessor: Vec<Option<usize>>,
    max_length: usize,
    coefficients: RefCell<Cache>,
}

impl Hall {
    pub(super) fn new(model: &mut Model, algebra: &Nakayama, key: &Key) -> Result<Self> {
        let vertices = algebra.paths.len();
        let field = Field::new(model, &key.field)?;
        let shapes: Vec<_> = algebra.resolutions.keys().copied().collect();
        let path_heads: Vec<Vec<_>> = (0..vertices)
            .map(|i| {
                (0..algebra.paths[i].len())
                    .map(|length| algebra.end(Path { start: i, length }))
                    .collect()
            })
            .collect();
        let dimensions = shapes
            .iter()
            .map(|state| {
                let mut dim = vec![0; vertices];
                for &vertex in &path_heads[state.head][..state.length] {
                    dim[vertex] += 1;
                }
                dim
            })
            .collect();
        let mut successor = vec![None; vertices];
        let mut predecessor = vec![None; vertices];
        for &a in &algebra.arrow_paths {
            let target = algebra.end(a);
            successor[a.start] = Some(target);
            predecessor[target] = Some(a.start);
        }
        let max_length = path_heads.iter().map(Vec::len).max().unwrap();
        Ok(Self {
            key: key.clone(),
            field,
            vertices,
            shapes,
            dimensions,
            path_heads,
            successor,
            predecessor,
            max_length,
            coefficients: RefCell::default(),
        })
    }

    pub(crate) fn literal(
        &self,
        model: &mut Model,
        parts: Vec<(u128, u128, u128)>,
    ) -> Result<Record> {
        let mut shape = BTreeMap::<(u8, u8), u16>::new();
        let mut length = 0;
        for (i, l, m) in parts {
            model.spend(1)?;
            let (i, l, m) = self.part(i, l, m)?;
            if m == 0 {
                continue;
            }
            length += usize::from(l) * usize::from(m);
            if length > MAX_MODULE_LENGTH {
                return Err(module_bound());
            }
            *shape.entry((i, l)).or_default() += m;
        }
        Ok(Record {
            grade: Grade {
                key: self.key.clone(),
                mode: Mode::Module,
            },
            data: Data::Module(Shape(
                shape.into_iter().map(|((i, l), m)| (i, l, m)).collect(),
            )),
        })
    }

    pub(crate) fn unit(&self) -> Record {
        Record {
            grade: Grade {
                key: self.key.clone(),
                mode: Mode::Algebra,
            },
            data: Data::Algebra(BTreeMap::from([(Shape::default(), 1)])),
        }
    }

    fn part(&self, i: u128, l: u128, m: u128) -> Result<Part> {
        let head = usize::try_from(i).unwrap_or(usize::MAX);
        let length = usize::try_from(l).unwrap_or(usize::MAX);
        if self.shapes.binary_search(&State { head, length }).is_err() {
            return Err(Error::plain(
                "partition Hall indecomposable top or length is out of range",
            ));
        }
        let multiplicity = u16::try_from(m)
            .ok()
            .filter(|&m| usize::from(m) <= MAX_MODULE_LENGTH)
            .ok_or_else(module_bound)?;
        Ok((head as u8, length as u8, multiplicity))
    }

    fn dim(&self, x: &Shape) -> Dims {
        let mut dim = vec![0; self.vertices];
        for &(i, l, m) in &x.0 {
            for &v in &self.path_heads[usize::from(i)][..usize::from(l)] {
                dim[v] += usize::from(m);
            }
        }
        dim
    }

    fn candidate_shapes(&self, model: &mut Model, dim: Dims) -> Result<Vec<Shape>> {
        let mut shapes = Vec::new();
        self.fill_shapes(model, &mut shapes, &mut Vec::new(), dim)?;
        Ok(shapes)
    }

    fn fill_shapes(
        &self,
        model: &mut Model,
        output: &mut Vec<Shape>,
        parts: &mut Vec<Part>,
        dim: Dims,
    ) -> Result<()> {
        model.spend(1 + self.vertices)?;
        if dim.iter().all(|&d| d == 0) {
            output.push(Shape(parts.clone()));
            return Ok(());
        }
        let start = parts.last().map_or(0, |&(i, l, _)| {
            self.shapes.partition_point(|state| {
                (state.head, state.length) <= (usize::from(i), usize::from(l))
            })
        });
        for s in start..self.shapes.len() {
            let factor = &self.dimensions[s];
            let stop = factor
                .iter()
                .zip(&dim)
                .filter(|(v, _)| **v != 0)
                .map(|(v, d)| d / v)
                .min()
                .unwrap();
            let state = self.shapes[s];
            for count in 1..=stop {
                let remainder = dim.iter().zip(factor).map(|(d, f)| d - count * f).collect();
                parts.push((state.head as u8, state.length as u8, count as u16));
                self.fill_shapes(model, output, parts, remainder)?;
                parts.pop();
            }
        }
        Ok(())
    }

    pub(crate) fn product(&self, model: &mut Model, a: &Record, b: &Record) -> Result<Record> {
        let mut terms = BTreeMap::new();
        for (m, &x) in a.terms() {
            for (n, &y) in b.terms() {
                let scalar = x.checked_mul(y).ok_or_else(arithmetic)?;
                if m.0.is_empty() || n.0.is_empty() {
                    let shape = if m.0.is_empty() { n } else { m };
                    accumulate(&mut terms, shape.clone(), scalar)?;
                    continue;
                }
                require_hall_length(m.length() + n.length())?;
                let dim = self
                    .dim(m)
                    .into_iter()
                    .zip(self.dim(n))
                    .map(|(a, b)| a + b)
                    .collect();
                for l in self.candidate_shapes(model, dim)? {
                    let coefficient = self.coefficient(model, &l, m, n)?;
                    let coefficient = i128::try_from(coefficient).map_err(|_| arithmetic())?;
                    accumulate(
                        &mut terms,
                        l,
                        scalar.checked_mul(coefficient).ok_or_else(arithmetic)?,
                    )?;
                }
            }
        }
        Ok(Record {
            grade: a.grade.clone(),
            data: Data::Algebra(terms),
        })
    }

    pub(crate) fn coefficient_report(
        &self,
        model: &mut Model,
        l: &Record,
        m: &Record,
        n: &Record,
    ) -> Result<String> {
        let value = self.coefficient(model, l.module(), m.module(), n.module())?;
        Ok(format!(
            "rvPartitionHallCoefficient[q={}] = {value}",
            self.field.order()
        ))
    }

    pub(crate) fn flag_report(
        &self,
        model: &mut Model,
        target: &Record,
        factors: [&Record; 3],
    ) -> Result<String> {
        let l = target.module();
        let factors: Vec<_> = factors
            .into_iter()
            .map(Record::module)
            .filter(|x| !x.0.is_empty())
            .collect();
        let number = match factors.as_slice() {
            [] => u128::from(l.0.is_empty()),
            [m] => u128::from(l == *m),
            [m, n] => self.coefficient(model, l, m, n)?,
            [m, n, k] => {
                let mn: Dims = self
                    .dim(m)
                    .into_iter()
                    .zip(self.dim(n))
                    .map(|(a, b)| a + b)
                    .collect();
                let nk: Dims = self
                    .dim(n)
                    .into_iter()
                    .zip(self.dim(k))
                    .map(|(a, b)| a + b)
                    .collect();
                if !self
                    .dim(l)
                    .iter()
                    .zip(&mn)
                    .zip(self.dim(k))
                    .all(|((&a, &b), c)| a == b + c)
                {
                    0
                } else {
                    require_hall_length(l.length())?;
                    let mut left = 0u128;
                    for x in self.candidate_shapes(model, mn)? {
                        let outer = self.coefficient(model, l, &x, k)?;
                        let inner = self.coefficient(model, &x, m, n)?;
                        left = outer
                            .checked_mul(inner)
                            .and_then(|x| x.checked_add(left))
                            .ok_or_else(arithmetic)?;
                    }
                    let mut right = 0u128;
                    for y in self.candidate_shapes(model, nk)? {
                        let outer = self.coefficient(model, l, m, &y)?;
                        let inner = self.coefficient(model, &y, n, k)?;
                        right = outer
                            .checked_mul(inner)
                            .and_then(|x| x.checked_add(right))
                            .ok_or_else(arithmetic)?;
                    }
                    if left != right {
                        return Err(Error::plain(
                            "partition Hall flag counts disagree between the two bracketings",
                        ));
                    }
                    left
                }
            }
            _ => unreachable!(),
        };
        Ok(format!(
            "rvPartitionHallFlag[q={}] = {number}",
            self.field.order()
        ))
    }

    fn coefficient(&self, model: &mut Model, l: &Shape, m: &Shape, n: &Shape) -> Result<u128> {
        if m.0.is_empty() {
            return Ok(u128::from(l == n));
        }
        if n.0.is_empty() {
            return Ok(u128::from(l == m));
        }
        let dims = self.dim(n);
        if !self
            .dim(l)
            .iter()
            .zip(self.dim(m))
            .zip(&dims)
            .all(|((&a, b), &c)| a == b + c)
        {
            return Ok(0);
        }
        require_hall_length(l.length())?;
        let table = self.distribution(model, l, &dims)?;
        Ok(table.get(&(m.clone(), n.clone())).copied().unwrap_or(0))
    }

    fn distribution(
        &self,
        model: &mut Model,
        module: &Shape,
        dim: &Dims,
    ) -> Result<Rc<Distribution>> {
        let key = (module.clone(), dim.clone());
        if let Some(value) = self.coefficients.borrow().get(&key) {
            return Ok(value.clone());
        }
        require_hall_length(module.length())?;
        let ans = Rc::new(subspaces::count(model, self, module, dim)?);
        model.allocate(
            8 + dim.len()
                + 3 * module.0.len()
                + ans
                    .keys()
                    .map(|(a, b)| 6 + 3 * (a.0.len() + b.0.len()))
                    .sum::<usize>(),
        )?;
        self.coefficients.borrow_mut().insert(key, ans.clone());
        Ok(ans)
    }
}

fn accumulate(entries: &mut BTreeMap<Shape, i128>, term: Shape, coefficient: i128) -> Result<()> {
    let sum = entries
        .get(&term)
        .copied()
        .unwrap_or(0)
        .checked_add(coefficient)
        .ok_or_else(arithmetic)?;
    if sum == 0 {
        entries.remove(&term);
    } else {
        entries.insert(term, sum);
    }
    Ok(())
}

fn arithmetic() -> Error {
    Error::plain("partition Hall exact coefficient arithmetic overflow")
}
fn module_bound() -> Error {
    Error::plain("partition Hall module length exceeds the limit of 64")
}
fn require_hall_length(length: usize) -> Result<()> {
    if length > MAX_HALL_LENGTH {
        return Err(Error::plain(
            "partition Hall submodule enumeration is limited to total length 6",
        ));
    }
    Ok(())
}
