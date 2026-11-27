use std::cell::RefCell;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::Model;
use crate::syntax::ParsedExpr;

use super::super::{add, multiply, ratio, Coordinates, Element, Fiber, Mode, Parameter};
use super::Gabriel;

pub(crate) mod auslander;
mod bar;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Grade {
    pub(crate) parameter: Parameter,
    pub(crate) source: u128,
    pub(crate) target: u128,
    pub(crate) degree: u128,
}

impl Grade {
    pub(crate) fn render(self) -> String {
        let Parameter { rank, n } = self.parameter;
        format!(
            "RVPartitionExt({rank},{n},{},{},{})",
            self.source, self.target, self.degree
        )
    }

    pub(crate) fn compose(beta: Self, alpha: Self) -> Result<Self> {
        if beta.parameter != alpha.parameter || beta.source != alpha.target {
            return Err(Error::plain(
                "Yoneda products require the same partition fiber and matching intermediate simple",
            ));
        }
        let degree = alpha
            .degree
            .checked_add(beta.degree)
            .ok_or_else(|| Error::plain("partition Ext degree overflow"))?;
        Ok(Self {
            parameter: alpha.parameter,
            source: alpha.source,
            target: beta.target,
            degree,
        })
    }
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Class {
    pub(crate) grade: Grade,
    pub(crate) coefficient: Rational,
}

impl Class {
    pub(crate) fn zero(grade: Grade) -> Self {
        Self {
            grade,
            coefficient: Rational::ZERO,
        }
    }

    pub(crate) fn add(self, rhs: Self) -> Result<Self> {
        Ok(Self {
            grade: self.grade,
            coefficient: add(self.coefficient, rhs.coefficient)?,
        })
    }

    pub(crate) fn scale(self, numerator: i128, denominator: i128) -> Result<Self> {
        Ok(Self {
            grade: self.grade,
            coefficient: multiply(self.coefficient, ratio(numerator, denominator)?)?,
        })
    }

    pub(crate) fn render(self) -> String {
        format!("{}{{{}}}", self.grade.render(), self.coefficient)
    }
}

pub(crate) enum Expression {
    Generator(Grade),
    Product(Box<[ParsedExpr; 2]>),
    Scale(i128, i128, Box<ParsedExpr>),
    Differential(Parameter, [u128; 3]),
    Comparison(Grade, u128),
    Massey(Parameter, [u128; 3]),
}

pub(crate) enum Report {
    Dimension(Grade, bool),
    Resolution(Parameter, [u128; 2]),
    Global(Parameter),
    Coefficient(Box<ParsedExpr>),
    Massey(Parameter, [u128; 3]),
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct Path {
    start: usize,
    length: usize,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct State {
    head: usize,
    length: usize,
}

struct LiftedPath {
    end: usize,
    value: Coordinates,
}

pub(crate) struct Nakayama {
    parameter: Parameter,
    paths: Vec<Vec<LiftedPath>>,
    arrow_paths: Vec<Path>,
    simple_dimensions: Vec<usize>,
    auslander: RefCell<Option<Rc<auslander::Auslander>>>,
    resolutions: BTreeMap<State, Resolution>,
}

impl Nakayama {
    pub(in crate::partition_algebra::fiber) fn new(
        model: &mut Model,
        fiber: &Fiber,
        data: &Gabriel,
    ) -> Result<Self> {
        let vertices = data.groups.len();
        let mut outgoing = vec![None; vertices];
        let mut incoming = vec![None; vertices];
        let mut arrow_paths = Vec::new();
        for (i, arrow) in data.arrows.iter().enumerate() {
            if outgoing[arrow.source].replace(i).is_some()
                || incoming[arrow.target].replace(i).is_some()
            {
                return Err(Error::plain(
                    "partition basic quiver does not have Nakayama incidence",
                ));
            }
            arrow_paths.push(Path {
                start: arrow.source,
                length: 1,
            });
        }
        let mut paths = Vec::with_capacity(vertices);
        for (i, group) in data.groups.iter().enumerate() {
            let mut entries = vec![LiftedPath {
                end: i,
                value: data.primitive[group[0]].clone(),
            }];
            let length: usize = data.cartan.iter().map(|row| row[i]).sum();
            for _ in 1..length {
                let last = entries.last().unwrap();
                let arrow = &data.arrows[outgoing[last.end]
                    .ok_or_else(|| Error::plain("partition path basis terminates prematurely"))?];
                let value = fiber.multiply(model, &arrow.lift, &last.value)?;
                if value.is_empty() {
                    return Err(Error::plain(
                        "partition path basis has a premature zero product",
                    ));
                }
                entries.push(LiftedPath {
                    end: arrow.target,
                    value,
                });
            }
            let mut span = super::Span::default();
            let mut counts = vec![0; vertices];
            for term in &entries {
                if !span.insert(model, term.value.clone())? {
                    return Err(Error::plain(
                        "partition projective paths are linearly dependent",
                    ));
                }
                counts[term.end] += 1;
            }
            if counts
                .iter()
                .enumerate()
                .any(|(j, &count)| count != data.cartan[j][i])
            {
                return Err(Error::plain(
                    "partition projective paths do not exhaust their actual corners",
                ));
            }
            paths.push(entries);
        }
        let mut algebra = Self {
            parameter: fiber.parameter,
            paths,
            arrow_paths,
            simple_dimensions: data.groups.iter().map(Vec::len).collect(),
            auslander: RefCell::default(),
            resolutions: BTreeMap::new(),
        };
        algebra.certify_products(model, fiber)?;
        for i in 0..vertices {
            for length in 1..=algebra.paths[i].len() {
                let state = State { head: i, length };
                let resolution = Resolution::new(&algebra, state)?;
                model.allocate(3 + 4 * resolution.states.len())?;
                algebra.resolutions.insert(state, resolution);
            }
        }
        model.allocate(
            algebra
                .paths
                .iter()
                .flatten()
                .map(|p| 6 + 6 * p.value.len())
                .sum(),
        )?;
        Ok(algebra)
    }

    pub(crate) fn auslander(&self, model: &mut Model) -> Result<Rc<auslander::Auslander>> {
        if let Some(data) = self.auslander.borrow().as_ref() {
            return Ok(data.clone());
        }
        let data = Rc::new(auslander::Auslander::new(model, self)?);
        *self.auslander.borrow_mut() = Some(data.clone());
        Ok(data)
    }

    fn certify_products(&self, model: &mut Model, fiber: &Fiber) -> Result<()> {
        for a in self.all_paths() {
            for b in self.all_paths() {
                let actual = fiber.multiply(model, &self.lifted(b).value, &self.lifted(a).value)?;
                let valid = match self.then(a, b) {
                    Some(c) => actual == self.lifted(c).value,
                    None => actual.is_empty(),
                };
                if !valid {
                    return Err(Error::plain(
                        "partition Nakayama path multiplication failed in the diagram algebra",
                    ));
                }
            }
        }
        Ok(())
    }

    fn all_paths(&self) -> impl Iterator<Item = Path> + '_ {
        self.paths
            .iter()
            .enumerate()
            .flat_map(|(start, row)| (0..row.len()).map(move |length| Path { start, length }))
    }

    fn lifted(&self, path: Path) -> &LiftedPath {
        &self.paths[path.start][path.length]
    }

    fn end(&self, path: Path) -> usize {
        self.lifted(path).end
    }

    // Composition is in traversal order. These structure constants were
    // checked against the actual diagram multiplication, including zeros.
    fn then(&self, a: Path, b: Path) -> Option<Path> {
        let length = a.length + b.length;
        (self.end(a) == b.start && length < self.paths[a.start].len()).then_some(Path {
            start: a.start,
            length,
        })
    }

    fn vertex(&self, i: u128) -> Result<usize> {
        usize::try_from(i)
            .ok()
            .filter(|&i| i < self.paths.len())
            .ok_or_else(|| Error::plain("partition simple vertex is out of range"))
    }

    fn simple(&self, i: usize) -> &Resolution {
        &self.resolutions[&State { head: i, length: 1 }]
    }

    fn resolution(&self, head: u128, length: u128) -> Result<&Resolution> {
        let head = self.vertex(head)?;
        let length = usize::try_from(length)
            .ok()
            .filter(|&l| l > 0 && l <= self.paths[head].len())
            .ok_or_else(|| {
                Error::plain("partition module length must be in 1..=projective length")
            })?;
        Ok(&self.resolutions[&State { head, length }])
    }

    pub(crate) fn dimension(&self, grade: Grade) -> Result<usize> {
        let i = self.vertex(grade.source)?;
        let j = self.vertex(grade.target)?;
        Ok(usize::from(
            self.simple(i)
                .at(grade.degree)
                .is_some_and(|state| state.head == j),
        ))
    }

    pub(crate) fn generator(&self, grade: Grade) -> Result<Class> {
        if self.dimension(grade)? == 0 {
            return Err(Error::plain("partition Ext group has dimension zero"));
        }
        Ok(Class {
            grade,
            coefficient: Rational::integer(1),
        })
    }

    pub(crate) fn dimension_report(
        &self,
        model: &mut Model,
        grade: Grade,
        bar: bool,
    ) -> Result<String> {
        let dimension = self.dimension(grade)?;
        let name = if bar {
            let independent = self.bar_dimension(model, grade)?;
            if independent != dimension {
                return Err(Error::plain(
                    "partition relative bar and projective Ext dimensions disagree",
                ));
            }
            "rvPartitionBarExtDim"
        } else {
            "rvPartitionExtDim"
        };
        let Parameter { rank, n } = self.parameter;
        Ok(format!(
            "{name}[{rank},{n},{},{},{}] = {dimension}",
            grade.source, grade.target, grade.degree
        ))
    }

    pub(crate) fn resolution_report(&self, head: u128, length: u128) -> Result<String> {
        let resolution = self.resolution(head, length)?;
        let states: Vec<_> = resolution
            .states
            .iter()
            .map(|s| (s.head, s.length))
            .collect();
        let Parameter { rank, n } = self.parameter;
        let mut out =
            format!("rvPartitionResolution[{rank},{n},{head},{length}] = syzygies={states:?}");
        if let Some(start) = resolution.repeat {
            write!(
                out,
                ", repeatFrom={start}, period={}",
                resolution.states.len() - start
            )
            .unwrap();
        } else {
            write!(out, ", projectiveDimension={}", resolution.states.len() - 1).unwrap();
        }
        Ok(out)
    }

    pub(crate) fn global_report(&self) -> String {
        let max = (0..self.paths.len())
            .map(|i| self.simple(i))
            .try_fold(0, |largest, r| {
                r.repeat
                    .is_none()
                    .then_some(largest.max(r.states.len() - 1))
            });
        let Parameter { rank, n } = self.parameter;
        let value = max.map_or_else(|| "infinity".into(), |v| v.to_string());
        format!("rvPartitionGlobalDimension[{rank},{n}] = {value}")
    }

    pub(crate) fn differential(&self, fiber: &Fiber, key: [u128; 3]) -> Result<Element> {
        let [head, length, step] = key;
        if step == 0 {
            return Err(Error::plain(
                "partition differential index must be positive",
            ));
        }
        let resolution = self.resolution(head, length)?;
        let value = if resolution.at(step).is_some() {
            let path = resolution.differential(step).unwrap();
            self.lifted(path).value.clone()
        } else {
            Coordinates::new()
        };
        Ok(fiber.element(Mode::Full, value))
    }

    pub(crate) fn comparison_diagram(
        &self,
        model: &mut Model,
        fiber: &Fiber,
        grade: Grade,
        step: u128,
    ) -> Result<Element> {
        let value = self
            .comparison(model, grade, step)?
            .map_or_else(Coordinates::new, |path| self.lifted(path).value.clone());
        Ok(fiber.element(Mode::Full, value))
    }

    fn comparison(&self, model: &mut Model, grade: Grade, cap: u128) -> Result<Option<Path>> {
        if self.dimension(grade)? == 0 {
            return Err(Error::plain(
                "partition comparison lift requires a nonzero Ext group",
            ));
        }
        let a = self.simple(self.vertex(grade.source)?);
        let b = self.simple(self.vertex(grade.target)?);
        let offset = grade.degree;
        let mut lift = Path {
            start: self.vertex(grade.target)?,
            length: 0,
        };
        if a.at(offset).map(|s| s.head) != Some(lift.start) {
            return Err(Error::plain(
                "partition comparison lift has an incorrect initial head",
            ));
        }
        let final_degree = offset
            .checked_add(cap)
            .ok_or_else(|| Error::plain("partition Ext degree overflow"))?;
        if a.at(final_degree).is_none() || b.at(cap).is_none() {
            return Ok(None);
        }
        let mut visited = BTreeMap::new();
        let mut k = 0;
        while k < cap {
            model.spend(1)?;
            let key = (a.phase(offset + k), b.phase(k), lift);
            if let Some(&prior) = visited.get(&key) {
                let period = k - prior;
                let steps = (cap - k) / period * period;
                if steps != 0 {
                    k += steps;
                    if k == cap {
                        break;
                    }
                }
            } else {
                visited.insert(key, k);
            }
            k += 1;
            let an = a.differential(offset + k).ok_or_else(|| {
                Error::plain("partition source comparison resolution terminates early")
            })?;
            let bn = b.differential(k).ok_or_else(|| {
                Error::plain("partition target comparison resolution terminates early")
            })?;
            let Some(left) = self.then(lift, an) else {
                // Zero at this step extends by zero and has zero terminal evaluation.
                return Ok(None);
            };
            let length = left.length.checked_sub(bn.length).ok_or_else(|| {
                Error::plain(
                    "partition comparison path has no factor through the target differential",
                )
            })?;
            let next = Path {
                start: self.end(bn),
                length,
            };
            if next.length >= self.paths[next.start].len()
                || self.then(bn, next) != Some(left)
                || self.end(next) != self.end(an)
            {
                return Err(Error::plain(
                    "partition path comparison failed its chain-map identity",
                ));
            }
            lift = next;
        }
        Ok(Some(lift))
    }

    pub(crate) fn yoneda(&self, model: &mut Model, beta: Class, alpha: Class) -> Result<Class> {
        let grade = Grade::compose(beta.grade, alpha.grade)?;
        if beta.coefficient == Rational::ZERO
            || alpha.coefficient == Rational::ZERO
            || self.dimension(grade)? == 0
        {
            return Ok(Class::zero(grade));
        }
        let lift = self.comparison(model, alpha.grade, beta.grade.degree)?;
        let coefficient = if lift.is_some_and(|path| path.length == 0) {
            multiply(alpha.coefficient, beta.coefficient)?
        } else {
            Rational::ZERO
        };
        Ok(Class { grade, coefficient })
    }
}

struct Resolution {
    states: Vec<State>,
    repeat: Option<usize>,
}

impl Resolution {
    fn new(algebra: &Nakayama, initial: State) -> Result<Self> {
        let mut current = initial;
        let mut states = Vec::new();
        let mut previous = None;
        loop {
            if let Some(start) = states.iter().position(|s| *s == current) {
                let next = Path {
                    start: current.head,
                    length: current.length,
                };
                if previous.is_some_and(|old| algebra.then(old, next).is_some()) {
                    return Err(Error::plain(
                        "partition periodic differential failed square-zero at the cycle join",
                    ));
                }
                return Ok(Self {
                    states,
                    repeat: Some(start),
                });
            }
            states.push(current);
            let length = algebra.paths[current.head].len();
            if current.length == length {
                return Ok(Self {
                    states,
                    repeat: None,
                });
            }
            let path = Path {
                start: current.head,
                length: current.length,
            };
            let next = State {
                head: algebra.end(path),
                length: length - current.length,
            };
            let mut image_length = 0;
            for l in 0..algebra.paths[next.head].len() {
                let actual = algebra.then(
                    path,
                    Path {
                        start: next.head,
                        length: l,
                    },
                );
                if l < next.length {
                    if actual
                        != Some(Path {
                            start: current.head,
                            length: current.length + l,
                        })
                    {
                        return Err(Error::plain(
                            "partition projective cover missed a syzygy basis path",
                        ));
                    }
                    image_length += 1;
                } else if actual.is_some() {
                    return Err(Error::plain(
                        "partition projective differential kernel is incorrect",
                    ));
                }
            }
            if image_length != next.length
                || previous.is_some_and(|old| algebra.then(old, path).is_some())
            {
                return Err(Error::plain(
                    "partition resolution failed exactness or square-zero",
                ));
            }
            previous = Some(path);
            current = next;
        }
    }

    fn phase(&self, degree: u128) -> Option<usize> {
        if degree < self.states.len() as u128 {
            return Some(degree as usize);
        }
        let start = self.repeat?;
        Some(start + ((degree - start as u128) % (self.states.len() - start) as u128) as usize)
    }

    fn at(&self, degree: u128) -> Option<State> {
        self.phase(degree).map(|i| self.states[i])
    }

    fn differential(&self, degree: u128) -> Option<Path> {
        self.at(degree)?;
        let prior = self.at(degree.checked_sub(1)?)?;
        Some(Path {
            start: prior.head,
            length: prior.length,
        })
    }
}
