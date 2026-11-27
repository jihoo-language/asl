use std::collections::HashMap;
use std::rc::Rc;

use crate::dependent_sets::Future;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::set_family::{self, Diagram};
use crate::set_logic::Omega;
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Command {
    Subset(Label, Label, Vec<Label>),
    FromMap(Label, Box<ParsedExpr>),
    Lattice(Label, Label, Label, bool),
    Support(Label),
    JoinSteps(Label, Label),
    Reflect(Label, Label, Label),
    Check(Label, Label),
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Map,
    Unit,
    Inverse,
    Extend,
}

impl Operation {
    pub(crate) fn detect(name: &str) -> Option<Self> {
        Some(match name {
            "setTopMap" => Self::Map,
            "setReflectUnit" => Self::Unit,
            "setReflectInverse" => Self::Inverse,
            "setReflectExtend" => Self::Extend,
            _ => return None,
        })
    }
}

pub(crate) struct Topology {
    pub(crate) omega_name: Rc<str>,
    omega: Rc<Omega>,
    selected: Box<[bool]>,
    rows: Rc<[usize]>,
}

impl Topology {
    pub(crate) fn subset(
        model: &mut Model,
        omega_name: Rc<str>,
        omega: Rc<Omega>,
        labels: &[Label],
    ) -> Result<Self> {
        let site = model.poset(omega.base).expect("named site");
        let mut selected = vec![false; site.names.len()];
        model.spend(site.names.len() + labels.len())?;
        for label in labels {
            selected[model.point(omega.base, label)?] = true;
        }
        Self::build(model, omega_name, omega, selected)
    }

    fn build(
        model: &mut Model,
        omega_name: Rc<str>,
        omega: Rc<Omega>,
        selected: Vec<bool>,
    ) -> Result<Self> {
        let rows = omega.nucleus_rows(model, &selected)?;
        model.allocate(1 + selected.len() + rows.len())?;
        Ok(Self {
            omega_name,
            omega,
            selected: selected.into(),
            rows: rows.into(),
        })
    }

    pub(crate) fn from_map(
        model: &mut Model,
        omega_name: Rc<str>,
        omega: Rc<Omega>,
        diagram: &Diagram,
        rows: Rc<[usize]>,
    ) -> Result<Self> {
        let selected = omega.nucleus_support(model, diagram, &rows)?;
        model.allocate(1 + selected.len())?;
        Ok(Self {
            omega_name,
            omega,
            selected: selected.into(),
            rows,
        })
    }

    pub(crate) fn base(&self) -> TypeId {
        self.omega.base
    }

    pub(crate) fn map(&self) -> Value {
        Value::SetNat(self.rows.clone())
    }

    pub(crate) fn members(&self, model: &mut Model) -> Result<Vec<String>> {
        model.spend(self.selected.len())?;
        let site = model.poset(self.base()).expect("named site");
        Ok(site
            .names
            .iter()
            .zip(self.selected.iter())
            .filter(|(_, in_support)| **in_support)
            .map(|(name, _)| name.clone())
            .collect())
    }

    fn combined_support(&self, model: &mut Model, right: &Self, meet: bool) -> Result<Vec<bool>> {
        if self.base() != right.base() {
            return Err(Error::plain(
                "set topologies must be on the same named site",
            ));
        }
        model.spend(self.selected.len())?;
        Ok(self
            .selected
            .iter()
            .zip(right.selected.iter())
            .map(|(&a, &b)| if meet { a || b } else { a && b })
            .collect())
    }

    pub(crate) fn combine(&self, model: &mut Model, right: &Self, meet: bool) -> Result<Self> {
        let selected = self.combined_support(model, right, meet)?;
        Self::build(model, self.omega_name.clone(), self.omega.clone(), selected)
    }

    pub(crate) fn join_steps(
        &self,
        model: &mut Model,
        diagram: &Diagram,
        right: &Self,
    ) -> Result<usize> {
        let selected = self.combined_support(model, right, false)?;
        let join = self.omega.nucleus_rows(model, &selected)?;
        let mut rounds = 0;
        for x in 0..self.selected.len() {
            let (a, b) = (
                diagram.components(x, &self.rows),
                diagram.components(x, &right.rows),
            );
            let join = diagram.components(x, &join);
            model.spend(a.len())?;
            let mut distances = vec![0; a.len()];
            let mut fixed = vec![0; a.len()];
            for i in (0..a.len()).rev() {
                let next = a[b[i]];
                // Each pass is extensive. Strict inclusion increases the
                // sorted bitmask index, so larger indices are already solved.
                if next < i {
                    return Err(Error::plain("topology composition was not extensive"));
                }
                if next == i {
                    fixed[i] = i;
                } else {
                    distances[i] = 1 + distances[next];
                    fixed[i] = fixed[next];
                }
                if fixed[i] != join[i] {
                    return Err(Error::plain(
                        "alternating closure did not match the topology join",
                    ));
                }
                rounds = rounds.max(distances[i]);
            }
        }
        Ok(rounds)
    }
}

struct Shape {
    futures: Vec<Future>,
    cells: usize,
    count: usize,
}

impl Shape {
    fn build(model: &mut Model, topology: &Topology, family: &Diagram) -> Result<Self> {
        if topology.base() != family.base {
            return Err(Error::plain(
                "reflection requires a family on the topology's named site",
            ));
        }
        let n = family.site.names.len();
        let mut cells = n + 1;
        let mut count = 0;
        let mut futures = Vec::with_capacity(n);
        for x in 0..n {
            model.spend(n)?;
            let scope = (0..n)
                .filter(|&y| topology.selected[y] && (x == y || family.site.strict_less(x, y)))
                .collect();
            let future = Future::build(model, family, scope)?;
            count += future.sections.len();
            if count > set_family::MAX_OCCURRENCES {
                return Err(Error::plain(
                    "topology reflection exceeds the limit of 4096 section occurrences",
                ));
            }
            cells += future.cells();
            set_family::check_cells(cells)?;
            futures.push(future);
        }
        Ok(Self {
            futures,
            cells,
            count,
        })
    }

    fn unit(&self, model: &mut Model, family: &Diagram) -> Result<Vec<usize>> {
        let mut rows = Vec::with_capacity(family.count());
        for (x, future) in self.futures.iter().enumerate() {
            let mut choices = Vec::with_capacity(future.points.len());
            model.spend(1 + family.size(x) * (1 + future.points.len()))?;
            for value in 0..family.size(x) {
                choices.clear();
                choices.extend(future.points.iter().map(|&y| family.apply(x, y, value)));
                rows.push(future.index(model, &choices)?);
            }
        }
        Ok(rows)
    }

    fn matching_sizes(&self, model: &mut Model, family: &Diagram) -> Result<bool> {
        model.spend(self.futures.len())?;
        Ok(self
            .futures
            .iter()
            .enumerate()
            .all(|(x, future)| future.sections.len() == family.size(x)))
    }
}

pub(crate) struct Reflection {
    pub(crate) topology: Rc<str>,
    pub(crate) source: Rc<str>,
    shape: Shape,
    unit: Rc<[usize]>,
}

impl Reflection {
    pub(crate) fn build(
        model: &mut Model,
        topology_name: Rc<str>,
        topology: &Topology,
        source: Rc<str>,
        family: &Diagram,
    ) -> Result<(Diagram, Self)> {
        let shape = Shape::build(model, topology, family)?;
        let sizes = shape
            .futures
            .iter()
            .map(|future| future.sections.len())
            .collect();
        let reflected = Diagram::generated(model, family.base, sizes, |model, x, y| {
            shape.futures[x].restrict(model, &shape.futures[y])
        })?;
        let unit = shape.unit(model, family)?;
        model.allocate(shape.cells + unit.len() + 1)?;
        Ok((
            reflected,
            Self {
                topology: topology_name,
                source,
                shape,
                unit: unit.into(),
            },
        ))
    }

    pub(crate) fn unit(&self) -> Value {
        Value::SetNat(self.unit.clone())
    }

    pub(crate) fn inverse(&self, model: &mut Model, source: &Diagram) -> Result<Value> {
        let inverse = inverse(model, source, &self.shape, &self.unit)?
            .ok_or_else(|| Error::plain("source of reflection is not a sheaf for its topology"))?;
        model.allocate(inverse.len())?;
        Ok(Value::SetNat(inverse.into()))
    }

    pub(crate) fn extend(
        &self,
        model: &mut Model,
        source: &Diagram,
        target: &Diagram,
        certificate: &Certificate,
        map: &[usize],
    ) -> Result<Value> {
        let mut rows = Vec::with_capacity(self.shape.count);
        for (x, future) in self.shape.futures.iter().enumerate() {
            let target_future = &certificate.shape.futures[x];
            let mut choices = Vec::with_capacity(future.points.len());
            model.spend(1 + (future.points.len() + 1) * future.sections.len())?;
            for section in &future.sections {
                choices.clear();
                choices.extend(
                    future
                        .points
                        .iter()
                        .zip(section)
                        .map(|(&point, &value)| map[source.occurrence(point, value)]),
                );
                let index = target_future.index(model, &choices)?;
                rows.push(certificate.inverse[target.occurrence(x, index)]);
            }
        }
        model.allocate(rows.len())?;
        Ok(Value::SetNat(rows.into()))
    }
}

pub(crate) type Certificates = HashMap<(Rc<str>, Rc<str>), Option<Rc<Certificate>>>;

pub(crate) struct Certificate {
    shape: Shape,
    inverse: Vec<usize>,
}

impl Certificate {
    pub(crate) fn build(
        model: &mut Model,
        topology: &Topology,
        family: &Diagram,
    ) -> Result<Option<Self>> {
        let shape = Shape::build(model, topology, family)?;
        if !shape.matching_sizes(model, family)? {
            return Ok(None);
        }
        let unit = shape.unit(model, family)?;
        let Some(inverse) = inverse(model, family, &shape, &unit)? else {
            return Ok(None);
        };
        model.allocate(shape.cells + inverse.len() + 1)?;
        Ok(Some(Self { shape, inverse }))
    }
}

fn inverse(
    model: &mut Model,
    family: &Diagram,
    shape: &Shape,
    unit: &[usize],
) -> Result<Option<Vec<usize>>> {
    if !shape.matching_sizes(model, family)? {
        return Ok(None);
    }
    model.spend(family.count())?;
    let mut inverse = vec![usize::MAX; family.count()];
    for x in 0..shape.futures.len() {
        for (source, &target) in family.components(x, unit).iter().enumerate() {
            let position = family.occurrence(x, target);
            let replaced = std::mem::replace(&mut inverse[position], source);
            if replaced != usize::MAX {
                return Ok(None);
            }
        }
    }
    Ok(Some(inverse))
}
