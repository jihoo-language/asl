use std::cell::RefCell;
use std::ops::Range;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::linear::{Row, Span};
use super::{require_hall_length, Data, Grade, Hall, Key, Mode, Record, Shape};

pub(crate) mod formal;
pub(crate) mod ginzburg;
mod mass;
mod matrix;
use matrix::{Matrix, Representation};

const CLASS_LIMIT: usize = 50_000;

pub(crate) enum Expression {
    Coordinates(Label, Vec<u128>, bool),
    Middle(Box<ParsedExpr>),
}

pub(crate) enum Report {
    Declare(Label, Box<[ParsedExpr; 2]>),
    Groupoid(Label),
    Distribution(Label),
    Boundary(Label, Vec<u128>),
    Lift(Box<ParsedExpr>),
    Hom(Box<[ParsedExpr; 2]>),
    Automorphisms(Box<ParsedExpr>),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Class {
    pub(crate) space: Rc<str>,
    coordinates: Row,
}

impl Class {
    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(3 + self.coordinates.len())?;
        Ok(Value::RVPartitionExtension(Rc::new(self)))
    }

    pub(crate) fn render(&self) -> String {
        format!("RVPartitionExtension({}){:?}", self.space, self.coordinates)
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<std::cmp::Ordering> {
        model.spend(1 + self.coordinates.len().min(rhs.coordinates.len()))?;
        Ok(self.cmp(rhs))
    }
}

pub(crate) struct Groupoid {
    name: Rc<str>,
    hall: Rc<Hall>,
    quotient: Shape,
    submodule: Shape,
    setup: Setup,
    equations: Span,
    boundaries: Span,
    classes: Span,
    hom: usize,
    distribution: RefCell<Option<String>>,
}

struct Setup {
    m: Representation,
    n: Representation,
    arrows: Vec<Option<Range<usize>>>,
    splittings: Vec<Range<usize>>,
    variables: usize,
    gauges: usize,
}

impl Setup {
    fn new(hall: &Hall, quotient: &Shape, submodule: &Shape) -> Self {
        let m = Representation::new(hall, quotient);
        let n = Representation::new(hall, submodule);
        let mut variables = 0;
        let arrows = hall
            .successor
            .iter()
            .enumerate()
            .map(|(s, &target)| {
                let t = target?;
                let width = m.dims[s] * n.dims[t];
                let range = variables..variables + width;
                variables += width;
                Some(range)
            })
            .collect();
        let mut gauges = 0;
        let splittings = m
            .dims
            .iter()
            .zip(&n.dims)
            .map(|(a, b)| {
                let width = a * b;
                let range = gauges..gauges + width;
                gauges += width;
                range
            })
            .collect();
        Self {
            m,
            n,
            arrows,
            splittings,
            variables,
            gauges,
        }
    }

    fn middle(&self, xi: &[u32]) -> Representation {
        let dims = self
            .m
            .dims
            .iter()
            .zip(&self.n.dims)
            .map(|(m, n)| m + n)
            .collect();
        let arrows = self
            .arrows
            .iter()
            .zip(&self.m.arrows)
            .zip(&self.n.arrows)
            .map(|((range, m), n)| {
                Some(Matrix::upper_block(
                    n.as_ref()?,
                    m.as_ref()?,
                    &xi[range.as_ref()?.clone()],
                ))
            })
            .collect();
        Representation { dims, arrows }
    }

    fn relations(&self, model: &mut Model, hall: &Hall, xi: &[u32]) -> Result<Row> {
        let middle = self.middle(xi);
        let mut equations = Vec::new();
        for (start, heads) in hall.path_heads.iter().enumerate() {
            let Some((target, path)) = middle.path(model, hall, start, heads.len())? else {
                continue;
            };
            for i in 0..path.rows {
                for j in 0..path.cols {
                    let value = path.get(i, j);
                    if i < self.n.dims[target] && j >= self.n.dims[start] {
                        equations.push(value);
                    } else if value != 0 {
                        return Err(Error::plain("partition extension diagonal blocks failed a defining zero-path relation"));
                    }
                }
            }
        }
        Ok(equations)
    }

    fn boundary(&self, model: &mut Model, hall: &Hall, gauge: &[u32]) -> Result<Row> {
        let splittings: Vec<_> = self
            .splittings
            .iter()
            .enumerate()
            .map(|(i, range)| {
                Matrix::new(
                    self.n.dims[i],
                    self.m.dims[i],
                    gauge[range.clone()].to_vec(),
                )
            })
            .collect();
        let mut xi = Vec::with_capacity(self.variables);
        for (source, &target) in hall.successor.iter().enumerate() {
            let Some(target) = target else {
                continue;
            };
            let m = self.m.arrows[source].as_ref().unwrap();
            let n = self.n.arrows[source].as_ref().unwrap();
            let left = n.compose(model, &hall.field, &splittings[source])?;
            let right = splittings[target].compose(model, &hall.field, m)?;
            for i in 0..left.rows {
                for j in 0..left.cols {
                    xi.push(hall.field.sub(left.get(i, j), right.get(i, j)));
                }
            }
        }
        Ok(xi)
    }
}

impl Groupoid {
    pub(crate) fn new(
        model: &mut Model,
        name: Rc<str>,
        hall: Rc<Hall>,
        m: &Record,
        n: &Record,
    ) -> Result<Self> {
        let quotient = m.module().clone();
        let submodule = n.module().clone();
        require_hall_length(quotient.length() + submodule.length())?;
        let setup = Setup::new(&hall, &quotient, &submodule);
        let field = &hall.field;
        let mut equations = Span::new(setup.variables);
        let base = setup.relations(model, &hall, &vec![0; setup.variables])?;
        if base.iter().any(|&x| x != 0) {
            return Err(Error::plain(
                "the split partition extension does not satisfy its relations",
            ));
        }
        let mut rows = vec![vec![0; setup.variables]; base.len()];
        // A product of triangular arrow blocks can use an off-diagonal
        // block at most once. Evaluation at a coordinate vector thus gives
        // the corresponding linear column of the actual path equation.
        for j in 0..setup.variables {
            let xi = coordinate_vector(setup.variables, j);
            let images = setup.relations(model, &hall, &xi)?;
            for (row, entry) in rows.iter_mut().zip(images) {
                row[j] = entry;
            }
        }
        for row in rows {
            equations.insert(model, field, row)?;
        }
        let mut boundaries = Span::new(setup.variables);
        for j in 0..setup.gauges {
            let xi = setup.boundary(model, &hall, &coordinate_vector(setup.gauges, j))?;
            if !equations.vanishes_on(model, field, &xi)? {
                return Err(Error::plain(
                    "partition vertex splitting failed the cocycle relations",
                ));
            }
            boundaries.insert(model, field, xi)?;
        }
        let hom = hall.intertwiner_dimension(&quotient, &submodule);
        if hom != setup.gauges - boundaries.rank {
            return Err(Error::plain(
                "partition extension stabilizer disagrees with the Hom dimension",
            ));
        }
        let mut classes = Span::new(setup.variables);
        for xi in equations.kernel(model, field)? {
            let reduced = boundaries.reduce(model, field, xi)?;
            classes.insert(model, field, reduced)?;
        }
        if classes.rank + boundaries.rank + equations.rank != setup.variables {
            return Err(Error::plain(
                "partition extension quotient dimensions do not add up",
            ));
        }
        if classes.rank != hall.projective_ext_dimension(&quotient, &submodule) {
            return Err(Error::plain(
                "partition cocycle Ext dimension disagrees with the projective presentation",
            ));
        }
        model.allocate(16 + 5 * setup.variables.pow(2) + 4 * setup.gauges)?;
        Ok(Self {
            name,
            hall,
            quotient,
            submodule,
            setup,
            equations,
            boundaries,
            classes,
            hom,
            distribution: RefCell::default(),
        })
    }

    pub(crate) fn key(&self) -> &Key {
        &self.hall.key
    }

    fn dimension(&self) -> usize {
        self.classes.rank
    }

    pub(crate) fn cardinality(&self) -> usize {
        (self.hall.field.order() as usize).saturating_pow(self.dimension() as u32)
    }

    pub(crate) fn values(&self, model: &mut Model) -> Result<Vec<Value>> {
        let count = self.cardinality();
        if count > 4096 {
            return Err(Error::plain(
                "partition extension space exceeds the finite type enumeration limit",
            ));
        }
        (0..count)
            .map(|i| {
                self.class(digits(i, self.hall.field.order(), self.dimension()))
                    .persist(model)
            })
            .collect()
    }

    fn class(&self, coordinates: Row) -> Class {
        Class {
            space: self.name.clone(),
            coordinates,
        }
    }

    pub(crate) fn zero(&self) -> Class {
        self.class(vec![0; self.dimension()])
    }

    pub(crate) fn literal(&self, model: &mut Model, values: Vec<u128>, raw: bool) -> Result<Class> {
        let size = if raw {
            self.setup.variables
        } else {
            self.dimension()
        };
        let values = self.valid_codes(values, size)?;
        if raw {
            self.project(model, values)
        } else {
            Ok(self.class(values))
        }
    }

    fn valid_codes(&self, values: Vec<u128>, size: usize) -> Result<Row> {
        if values.len() != size {
            return Err(Error::plain(format!(
                "partition extension vector requires {size} coordinates"
            )));
        }
        let q = u128::from(self.hall.field.order());
        values
            .into_iter()
            .map(|i| {
                (i < q).then_some(i as u32).ok_or_else(|| {
                    Error::plain("partition extension coordinate is outside the finite field codes")
                })
            })
            .collect()
    }

    fn project(&self, model: &mut Model, xi: Row) -> Result<Class> {
        let field = &self.hall.field;
        if !self.equations.vanishes_on(model, field, &xi)? {
            return Err(Error::plain(
                "partition off-diagonal arrow vector is not a cocycle",
            ));
        }
        let reduced = self.boundaries.reduce(model, field, xi)?;
        let coordinates = self
            .classes
            .coordinates(model, field, reduced)?
            .ok_or_else(|| {
                Error::plain("partition cocycle did not lie in the certified extension quotient")
            })?;
        Ok(self.class(coordinates))
    }

    fn lift(&self, model: &mut Model, class: &Class) -> Result<Row> {
        let field = &self.hall.field;
        let mut lift = vec![0; self.setup.variables];
        for (row, &coefficient) in self.classes.rows.iter().flatten().zip(&class.coordinates) {
            if coefficient == 0 {
                continue;
            }
            model.spend(row.len())?;
            for (x, &y) in lift.iter_mut().zip(row) {
                *x = field.add(*x, field.mul(coefficient, y));
            }
        }
        Ok(lift)
    }

    pub(crate) fn lift_report(&self, model: &mut Model, class: &Class) -> Result<String> {
        Ok(format!(
            "rvPartitionExtensionLift[{}] = {:?}",
            self.name,
            self.lift(model, class)?
        ))
    }

    pub(crate) fn boundary_report(&self, model: &mut Model, values: Vec<u128>) -> Result<String> {
        let values = self.valid_codes(values, self.setup.gauges)?;
        let xi = self.setup.boundary(model, &self.hall, &values)?;
        let class = self.project(model, xi.clone())?;
        if class != self.zero() {
            return Err(Error::plain("partition coboundary did not project to zero"));
        }
        Ok(format!(
            "rvPartitionExtensionCoboundary[{}] = cocycle={xi:?}, class={:?}",
            self.name, class.coordinates
        ))
    }

    pub(crate) fn add(&self, model: &mut Model, a: &Class, b: &Class) -> Result<Class> {
        model.spend(self.dimension())?;
        let field = &self.hall.field;
        Ok(self.class(
            a.coordinates
                .iter()
                .zip(&b.coordinates)
                .map(|(&x, &y)| field.add(x, y))
                .collect(),
        ))
    }

    pub(crate) fn negate(&self, model: &mut Model, a: &Class) -> Result<Class> {
        model.spend(self.dimension())?;
        Ok(self.class(
            a.coordinates
                .iter()
                .map(|&x| self.hall.field.sub(0, x))
                .collect(),
        ))
    }

    pub(crate) fn middle(&self, model: &mut Model, class: &Class) -> Result<Record> {
        let shape = self.middle_shape(model, class)?;
        Ok(Record {
            grade: Grade {
                key: self.hall.key.clone(),
                mode: Mode::Module,
            },
            data: Data::Module(shape),
        })
    }

    fn middle_shape(&self, model: &mut Model, class: &Class) -> Result<Shape> {
        let field = &self.hall.field;
        let xi = self.lift(model, class)?;
        let middle = self.setup.middle(&xi);
        let mut ranks = vec![vec![0; self.hall.max_length + 2]; self.hall.vertices];
        for (v, row) in ranks.iter_mut().enumerate() {
            for (k, rank) in row.iter_mut().enumerate().take(self.hall.max_length) {
                if let Some((_, path)) = middle.path(model, &self.hall, v, k)? {
                    *rank = path.rank(model, field)?;
                }
            }
            if let Some((_, zero)) =
                middle.path(model, &self.hall, v, self.hall.path_heads[v].len())?
            {
                if !zero.is_zero() {
                    return Err(Error::plain(
                        "partition extension middle matrix failed a defining relation",
                    ));
                }
            }
        }
        super::subspaces::classify(&self.hall, &ranks)
    }
}

fn coordinate_vector(size: usize, index: usize) -> Row {
    let mut row = vec![0; size];
    row[index] = 1;
    row
}

fn digits(mut number: usize, field: u32, length: usize) -> Row {
    let q = field as usize;
    (0..length)
        .map(|_| {
            let cell = (number % q) as u32;
            number /= q;
            cell
        })
        .collect()
}
