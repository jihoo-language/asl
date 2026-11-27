use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_ideal::Natural;
use crate::syntax::{Label, ParsedExpr};

use super::cubic_division::{Element, Extension, Generator, Order};

const Q: usize = 27;

pub(crate) enum Expression {
    Root(Label),
    At(Label, usize, usize),
    Normalizer(Label),
    Line(usize, Box<[ParsedExpr; 2]>),
    Matrix(Box<[ParsedExpr; 4]>),
    VertexUnary(VertexUnary, Box<ParsedExpr>),
    VertexBinary(VertexBinary, Box<[ParsedExpr; 2]>),
    Act(bool, Box<[ParsedExpr; 2]>),
    Compose(Box<[ParsedExpr; 2]>),
    Inverse(Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum VertexUnary {
    Ancestor(usize),
    Child(usize),
    Translate,
    TranslateBack,
}

#[derive(Clone, Copy)]
pub(crate) enum VertexBinary {
    Lca,
    Adjacent,
}

pub(crate) enum Report {
    Sphere(Label, usize),
    Level(Box<ParsedExpr>),
    Path(Box<ParsedExpr>),
    Distance(Box<[ParsedExpr; 2]>),
    Hecke(usize, Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Vertex {
    radius: usize,
    second: bool,
    coordinate: Element,
}

impl Vertex {
    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(48)?;
        Ok(Value::CubicTreeVertex(Rc::new(self)))
    }

    pub(crate) fn description(&self) -> String {
        if self.radius == 0 {
            return "radius=0, root".to_owned();
        }
        let chart = if self.second { "(x,1)" } else { "(1,x)" };
        format!(
            "radius={}, chart={chart}, x={:?}",
            self.radius,
            self.coordinate.rows()
        )
    }

    pub(crate) fn level(&self) -> usize {
        self.radius
    }
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Operator {
    entries: [Element; 4],
    power: usize,
}

impl Operator {
    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(182)?;
        Ok(Value::CubicTreeOperator(Rc::new(self)))
    }

    pub(crate) fn description(&self) -> String {
        format!(
            "normalizer_power={}, entries={:?}",
            self.power,
            self.entries.each_ref().map(Element::rows)
        )
    }
}

pub(crate) struct Tree<'a> {
    order: Order<'a>,
    depth: usize,
}

impl<'a> Tree<'a> {
    pub(crate) fn new(model: &mut Model, extension: &'a Extension) -> Result<Self> {
        Ok(Self {
            order: extension.ring(model)?,
            depth: extension.precision,
        })
    }

    pub(crate) fn require_radius(depth: usize, radius: usize) -> Result<()> {
        if radius > depth {
            return Err(Error::plain(
                "cubic tree radius exceeds the named division precision",
            ));
        }
        Ok(())
    }

    pub(crate) fn root(&self) -> Vertex {
        Vertex {
            radius: 0,
            second: false,
            coordinate: self.order.zero(),
        }
    }

    fn vertex(&self, radius: usize, second: bool, coordinate: Element) -> Result<Vertex> {
        Self::require_radius(self.depth, radius)?;
        if radius == 0 {
            return Ok(self.root());
        }
        if second && self.order.is_unit(&coordinate) {
            return Err(inconsistent());
        }
        Ok(Vertex {
            radius,
            second,
            coordinate: self.order.truncate(coordinate, radius)?,
        })
    }

    pub(crate) fn line(
        &self,
        model: &mut Model,
        radius: usize,
        a: &Element,
        b: &Element,
    ) -> Result<Vertex> {
        Self::require_radius(self.depth, radius)?;
        let (second, pivot, other) = if self.order.is_unit(a) {
            (false, a, b)
        } else if self.order.is_unit(b) {
            (true, b, a)
        } else {
            return Err(Error::plain(
                "cubic tree projective line requires a primitive right column",
            ));
        };
        if radius == 0 {
            return Ok(self.root());
        }
        let inverse = self.order.inverse(model, pivot)?;
        let coordinate = self.order.multiply(model, other, &inverse)?;
        self.vertex(radius, second, coordinate)
    }

    pub(crate) fn at(&self, model: &mut Model, n: usize, mut index: usize) -> Result<Vertex> {
        Self::require_radius(self.depth, n)?;
        if n == 0 {
            return if index == 0 {
                Ok(self.root())
            } else {
                Err(sphere_index())
            };
        }
        let first = (0..n).try_fold(1usize, |a, _| a.checked_mul(Q));
        let second = first.is_some_and(|first| index >= first);
        if let Some(first) = first.filter(|_| second) {
            index -= first;
        }
        let start = usize::from(second);
        let mut coordinate = self.order.zero();
        for grade in start..n {
            let digit = self.order.graded_digit(model, grade, index % Q)?;
            coordinate = self.order.add(model, &coordinate, &digit)?;
            index /= Q;
            if index == 0 {
                break;
            }
        }
        if index != 0 {
            return Err(sphere_index());
        }
        self.vertex(n, second, coordinate)
    }

    pub(crate) fn sphere_report(
        model: &mut Model,
        name: &str,
        depth: usize,
        n: usize,
    ) -> Result<String> {
        Self::require_radius(depth, n)?;
        let mut size = Natural::unit();
        if n > 0 {
            let mut first = Natural::default();
            first.add_multiple(model, &size, (Q + 1) as u128)?;
            size = first;
            for _ in 1..n {
                let mut next = Natural::default();
                next.add_multiple(model, &size, Q as u128)?;
                size = next;
            }
        }
        Ok(format!(
            "rvCubicTreeSphere[{name},{n}] = size={}",
            size.decimal()
        ))
    }

    pub(crate) fn ancestor(&self, x: &Vertex, n: usize) -> Result<Vertex> {
        if n > x.radius {
            return Err(Error::plain(
                "cubic tree ancestor radius exceeds the vertex level",
            ));
        }
        self.vertex(n, x.second, x.coordinate.clone())
    }

    pub(crate) fn child(&self, model: &mut Model, x: &Vertex, index: usize) -> Result<Vertex> {
        if x.radius == 0 {
            return self.at(model, 1, index);
        }
        let radius = x.radius + 1;
        Self::require_radius(self.depth, radius)?;
        if index >= Q {
            return Err(Error::plain(
                "cubic tree nonroot child index must be in 0..27",
            ));
        }
        let tail = self.order.graded_digit(model, x.radius, index)?;
        let coordinate = self.order.add(model, &x.coordinate, &tail)?;
        self.vertex(radius, x.second, coordinate)
    }

    fn common_level(&self, a: &Vertex, b: &Vertex) -> Result<usize> {
        if a.second != b.second {
            return Ok(0);
        }
        let mut lower = 0;
        let mut upper = a.radius.min(b.radius) + 1;
        while upper - lower > 1 {
            let mid = (lower + upper) / 2;
            if self.ancestor(a, mid)? == self.ancestor(b, mid)? {
                lower = mid;
            } else {
                upper = mid;
            }
        }
        Ok(lower)
    }

    pub(crate) fn lca(&self, a: &Vertex, b: &Vertex) -> Result<Vertex> {
        self.ancestor(a, self.common_level(a, b)?)
    }

    pub(crate) fn distance(&self, a: &Vertex, b: &Vertex) -> Result<usize> {
        Ok(a.radius + b.radius - 2 * self.common_level(a, b)?)
    }

    pub(crate) fn path_report(&self, model: &mut Model, x: &Vertex) -> Result<String> {
        use std::fmt::Write as _;
        let mut result = format!("rvCubicTreePath = {}", self.root().description());
        for level in 1..=x.radius {
            model.spend(1)?;
            let parent = self.ancestor(x, level)?;
            let _ = write!(result, " -> {}", parent.description());
        }
        Ok(result)
    }

    pub(crate) fn hecke_report(
        &self,
        model: &mut Model,
        a: &Vertex,
        b: &Vertex,
        r: usize,
    ) -> Result<String> {
        let distance = self.distance(a, b)?;
        let mut paths = 0;
        if a.radius > 0 {
            let parent = self.ancestor(a, a.radius - 1)?;
            paths += usize::from(self.distance(&parent, b)? == r);
        }
        let children = if a.radius == 0 { Q + 1 } else { Q };
        if a.radius < self.depth {
            for j in 0..children {
                let child = self.child(model, a, j)?;
                paths += usize::from(self.distance(&child, b)? == r);
            }
        } else if distance + 1 == r {
            // Every outward neighbor of a boundary vertex is farther
            // from any vertex in the certified ball.
            paths += children;
        }
        let expected = if distance.checked_sub(1) == Some(r) {
            1
        } else if distance + 1 == r {
            if distance == 0 {
                Q + 1
            } else {
                Q
            }
        } else {
            0
        };
        if paths != expected {
            return Err(inconsistent());
        }
        Ok(format!(
            "rvCubicTreeHecke[{r}] = distance={distance}, paths={paths}"
        ))
    }

    fn identity_entries(&self) -> [Element; 4] {
        [
            self.order.unit(1),
            self.order.zero(),
            self.order.zero(),
            self.order.unit(1),
        ]
    }

    pub(crate) fn normalizer(&self) -> Operator {
        Operator {
            entries: self.identity_entries(),
            power: 1,
        }
    }

    pub(crate) fn matrix(&self, model: &mut Model, entries: [Element; 4]) -> Result<Operator> {
        self.matrix_inverse(model, &entries)?;
        Ok(Operator { entries, power: 0 })
    }

    fn matrix_inverse(&self, model: &mut Model, m: &[Element; 4]) -> Result<[Element; 4]> {
        let [a, b, c, d] = m.clone();
        let one = self.order.unit(1);
        let zero = self.order.zero();
        let mut rows = [[a, b, one.clone(), zero.clone()], [c, d, zero, one]];
        if !self.order.is_unit(&rows[0][0]) {
            rows.swap(0, 1);
        }
        for pivot in [0, 1] {
            let (upper, lower) = rows.split_at_mut(1);
            let (pivot_row, other_row) = if pivot == 0 {
                (&mut upper[0], &mut lower[0])
            } else {
                (&mut lower[0], &mut upper[0])
            };
            if !self.order.is_unit(&pivot_row[pivot]) {
                return Err(Error::plain(
                    "cubic tree matrix must have invertible reduction modulo Pi",
                ));
            }
            let inv = self.order.inverse(model, &pivot_row[pivot])?;
            for entry in pivot_row.iter_mut() {
                *entry = self.order.multiply(model, &inv, entry)?;
            }
            let factor = other_row[pivot].clone();
            for (out, entry) in other_row.iter_mut().zip(pivot_row) {
                let product = self.order.multiply(model, &factor, entry)?;
                let negated = self.order.negate(model, &product)?;
                *out = self.order.add(model, out, &negated)?;
            }
        }
        let [[_, _, a, b], [_, _, c, d]] = rows;
        let inverse = [a, b, c, d];
        let identity = self.identity_entries();
        if self.matrix_product(model, m, &inverse)? != identity
            || self.matrix_product(model, &inverse, m)? != identity
        {
            return Err(inconsistent());
        }
        Ok(inverse)
    }

    fn dot(
        &self,
        model: &mut Model,
        a: &Element,
        b: &Element,
        x: &Element,
        y: &Element,
    ) -> Result<Element> {
        let ax = self.order.multiply(model, a, x)?;
        let by = self.order.multiply(model, b, y)?;
        self.order.add(model, &ax, &by)
    }

    fn matrix_product(
        &self,
        model: &mut Model,
        a: &[Element; 4],
        b: &[Element; 4],
    ) -> Result<[Element; 4]> {
        Ok([
            self.dot(model, &a[0], &a[1], &b[0], &b[2])?,
            self.dot(model, &a[0], &a[1], &b[1], &b[3])?,
            self.dot(model, &a[2], &a[3], &b[0], &b[2])?,
            self.dot(model, &a[2], &a[3], &b[1], &b[3])?,
        ])
    }

    fn twist(
        &self,
        model: &mut Model,
        power: usize,
        mut entries: [Element; 4],
    ) -> Result<[Element; 4]> {
        if power != 0 {
            for entry in &mut entries {
                *entry = self.order.conjugate_pi_power(model, power, entry)?;
            }
        }
        Ok(entries)
    }

    pub(crate) fn compose(
        &self,
        model: &mut Model,
        a: &Operator,
        b: &Operator,
    ) -> Result<Operator> {
        let rhs = self.twist(model, a.power, b.entries.clone())?;
        let entries = self.matrix_product(model, &a.entries, &rhs)?;
        Ok(Operator {
            entries,
            power: (a.power + b.power) % 3,
        })
    }

    pub(crate) fn inverse(&self, model: &mut Model, a: &Operator) -> Result<Operator> {
        let power = (3 - a.power) % 3;
        let inverse = self.matrix_inverse(model, &a.entries)?;
        let entries = self.twist(model, power, inverse)?;
        Ok(Operator { entries, power })
    }

    fn image_column(&self, model: &mut Model, g: &Operator, x: &Vertex) -> Result<[Element; 2]> {
        let coordinate = if g.power == 0 {
            x.coordinate.clone()
        } else {
            self.order
                .conjugate_pi_power(model, g.power, &x.coordinate)?
        };
        let one = self.order.unit(1);
        let (first, second) = if x.second {
            (&coordinate, &one)
        } else {
            (&one, &coordinate)
        };
        let [a, b, c, d] = &g.entries;
        Ok([
            self.dot(model, a, b, first, second)?,
            self.dot(model, c, d, first, second)?,
        ])
    }

    pub(crate) fn apply(&self, model: &mut Model, g: &Operator, x: &Vertex) -> Result<Vertex> {
        if x.radius == 0 {
            return Ok(self.root());
        }
        let [top, bottom] = self.image_column(model, g, x)?;
        self.line(model, x.radius, &top, &bottom)
    }

    pub(crate) fn fixes(&self, model: &mut Model, g: &Operator, x: &Vertex) -> Result<bool> {
        if x.radius == 0 {
            return Ok(true);
        }
        let [top, bottom] = self.image_column(model, g, x)?;
        let (pivot, other) = if x.second {
            (bottom, top)
        } else {
            (top, bottom)
        };
        if !self.order.is_unit(&pivot) {
            return Ok(false);
        }
        let expected = self.order.multiply(model, &x.coordinate, &pivot)?;
        Ok(self.order.truncate(expected, x.radius)? == self.order.truncate(other, x.radius)?)
    }

    pub(crate) fn translate(&self, model: &mut Model, x: &Vertex) -> Result<Vertex> {
        if x.radius == 0 {
            return self.vertex(1, true, self.order.zero());
        }
        if !x.second && !self.order.is_unit(&x.coordinate) {
            let divided = self.order.divide_right_pi(model, &x.coordinate)?;
            return self.vertex(x.radius - 1, false, divided);
        }
        if x.radius == self.depth {
            return Err(Error::plain(
                "cubic tree translation needs one more outward digit of ambient precision",
            ));
        }
        let pi = self.order.generator(Generator::Pi);
        let coordinate = if x.second {
            x.coordinate.clone()
        } else {
            self.order.inverse(model, &x.coordinate)?
        };
        let shifted = self.order.multiply(model, &pi, &coordinate)?;
        self.vertex(x.radius + 1, true, shifted)
    }

    pub(crate) fn translate_back(&self, model: &mut Model, x: &Vertex) -> Result<Vertex> {
        let zero = self.order.zero();
        let one = self.order.unit(1);
        let swap = Operator {
            entries: [zero.clone(), one.clone(), one, zero],
            power: 0,
        };
        let reflected = self.apply(model, &swap, x)?;
        let moved = self.translate(model, &reflected)?;
        let reflected = self.apply(model, &swap, &moved)?;
        let back_twist = Operator {
            entries: self.identity_entries(),
            power: 2,
        };
        self.apply(model, &back_twist, &reflected)
    }
}

fn inconsistent() -> Error {
    Error::plain("inconsistent cubic tree lattice or operator identity")
}

fn sphere_index() -> Error {
    Error::plain("cubic tree sphere index is out of range")
}
