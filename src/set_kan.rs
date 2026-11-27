use std::rc::Rc;

use crate::dependent_sets::{pullback, Future, Pullback};
use crate::error::{Error, Result};
use crate::model::{Model, Poset, TypeId, Value};
use crate::set_family::{self, Diagram};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    LanLift,
    LanRestrict,
    LanUnit,
    LanCounit,
    RanTranspose,
    RanUntranspose,
    RanUnit,
    RanCounit,
    LanComparison,
    RanComparison,
}

impl Operation {
    pub(crate) fn detect(name: &str) -> Option<Self> {
        Some(match name {
            "setLanLift" => Self::LanLift,
            "setLanRestrict" => Self::LanRestrict,
            "setLanUnit" => Self::LanUnit,
            "setLanCounit" => Self::LanCounit,
            "setRanTranspose" => Self::RanTranspose,
            "setRanUntranspose" => Self::RanUntranspose,
            "setRanUnit" => Self::RanUnit,
            "setRanCounit" => Self::RanCounit,
            "setLanBC" => Self::LanComparison,
            "setRanBC" => Self::RanComparison,
            _ => return None,
        })
    }

    pub(crate) fn second_label(self) -> bool {
        matches!(
            self,
            Self::LanUnit | Self::RanCounit | Self::LanComparison | Self::RanComparison
        )
    }

    pub(crate) fn takes_input(self) -> bool {
        matches!(
            self,
            Self::LanLift | Self::LanRestrict | Self::RanTranspose | Self::RanUntranspose
        )
    }
}

pub(crate) struct Along {
    pub(crate) body: Rc<str>,
    source: TypeId,
    target: TypeId,
    positions: Vec<usize>,
}

impl Along {
    pub(crate) fn new(
        model: &mut Model,
        body: Rc<str>,
        source: TypeId,
        target: TypeId,
        map: &[Value],
    ) -> Result<Self> {
        model.spend(map.len())?;
        let positions = map
            .iter()
            .map(|value| {
                let Value::Point(index) = *value else {
                    unreachable!("map of named sites")
                };
                index
            })
            .collect();
        Ok(Self {
            body,
            source,
            target,
            positions,
        })
    }

    pub(crate) fn check_pullback(&self, model: &mut Model, pullback: &Pullback) -> Result<()> {
        let family = model
            .set_family_named(&pullback.family)
            .expect("pullback's family");
        model.spend(self.positions.len() + 1)?;
        if self.source != pullback.base
            || self.target != family.base
            || self.positions != pullback.positions()
        {
            return Err(Error::plain(
                "expected a set family pulled back along the Kan extension's map",
            ));
        }
        Ok(())
    }
}

fn le(site: &Poset, from: usize, to: usize) -> bool {
    from == to || site.strict_less(from, to)
}

struct Components {
    parents: Vec<usize>,
    sizes: Vec<usize>,
}

impl Components {
    fn new(count: usize) -> Self {
        Self {
            parents: (0..count).collect(),
            sizes: vec![1; count],
        }
    }

    fn root(&mut self, model: &mut Model, mut node: usize) -> Result<usize> {
        while self.parents[node] != node {
            model.spend(1)?;
            self.parents[node] = self.parents[self.parents[node]];
            node = self.parents[node];
        }
        Ok(node)
    }

    fn merge(&mut self, model: &mut Model, left: usize, right: usize) -> Result<()> {
        let mut a = self.root(model, left)?;
        let mut b = self.root(model, right)?;
        if a == b {
            return Ok(());
        }
        if self.sizes[a] < self.sizes[b] {
            std::mem::swap(&mut a, &mut b);
        }
        self.parents[b] = a;
        self.sizes[a] += self.sizes[b];
        Ok(())
    }
}

struct Past {
    representatives: Vec<(usize, usize)>,
    classes: Vec<(usize, usize)>,
}

impl Past {
    fn build(model: &mut Model, body: &Diagram, points: &[usize]) -> Result<Self> {
        model.spend(3 * body.count() + points.len() * points.len())?;
        let mut sets = Components::new(body.count());
        for &x in points {
            if body.size(x) == 0 {
                continue;
            }
            for &z in points {
                if !body.site.strict_less(x, z) {
                    continue;
                }
                model.spend(body.size(x))?;
                for value in 0..body.size(x) {
                    sets.merge(
                        model,
                        body.occurrence(x, value),
                        body.occurrence(z, body.apply(x, z, value)),
                    )?;
                }
            }
        }
        let mut by_root = vec![None; body.count()];
        let mut representatives = Vec::new();
        let mut classes = Vec::new();
        for &x in points {
            for value in 0..body.size(x) {
                let slot = body.occurrence(x, value);
                let root = sets.root(model, slot)?;
                let class = *by_root[root].get_or_insert_with(|| {
                    let index = representatives.len();
                    representatives.push((x, value));
                    index
                });
                classes.push((slot, class));
            }
        }
        Ok(Self {
            representatives,
            classes,
        })
    }

    fn index(
        &self,
        model: &mut Model,
        body: &Diagram,
        point: usize,
        value: usize,
    ) -> Result<usize> {
        model.spend(1 + self.classes.len().max(1).ilog2() as usize)?;
        let slot = body.occurrence(point, value);
        let i = self
            .classes
            .binary_search_by_key(&slot, |&(key, _)| key)
            .expect("occurrence in past");
        Ok(self.classes[i].1)
    }

    fn cells(&self) -> usize {
        1 + 2 * (self.representatives.len() + self.classes.len())
    }
}

pub(crate) struct Lan {
    pub(crate) along: Along,
    pasts: Vec<Past>,
}

impl Lan {
    pub(crate) fn build(
        model: &mut Model,
        along: Along,
        body: &Diagram,
    ) -> Result<(Diagram, Self)> {
        let target = model.poset(along.target).expect("named target site");
        let mut pasts = Vec::with_capacity(target.names.len());
        let mut cells = along.positions.len() + 1;
        for y in 0..target.names.len() {
            model.spend(along.positions.len())?;
            let points: Vec<_> = along
                .positions
                .iter()
                .enumerate()
                .filter_map(|(x, &fx)| le(&target, fx, y).then_some(x))
                .collect();
            let past = Past::build(model, body, &points)?;
            cells += past.cells();
            set_family::check_cells(cells)?;
            pasts.push(past);
        }
        let sizes = pasts
            .iter()
            .map(|past| past.representatives.len())
            .collect();
        let diagram = Diagram::generated(model, along.target, sizes, |model, y, z| {
            pasts[y]
                .representatives
                .iter()
                .map(|&(x, a)| pasts[z].index(model, body, x, a))
                .collect()
        })?;
        model.allocate(cells)?;
        Ok((diagram, Self { along, pasts }))
    }

    fn lift(
        &self,
        model: &mut Model,
        body: &Diagram,
        target: &Diagram,
        rows: Option<&[usize]>,
    ) -> Result<Vec<usize>> {
        let mut images = Vec::new();
        for (y, past) in self.pasts.iter().enumerate() {
            model.spend(1 + past.representatives.len())?;
            images.extend(past.representatives.iter().map(|&(x, a)| {
                let value = component(rows, body, x, a);
                target.apply(self.along.positions[x], y, value)
            }));
        }
        Ok(images)
    }

    fn restrict(
        &self,
        model: &mut Model,
        extended: &Diagram,
        body: &Diagram,
        rows: Option<&[usize]>,
    ) -> Result<Vec<usize>> {
        let mut images = Vec::with_capacity(body.count());
        for (x, &y) in self.along.positions.iter().enumerate() {
            model.spend(1 + body.size(x))?;
            for value in 0..body.size(x) {
                let i = self.pasts[y].index(model, body, x, value)?;
                images.push(component(rows, extended, y, i));
            }
        }
        Ok(images)
    }

    fn compare(
        &self,
        model: &mut Model,
        old: &Self,
        square: &Square,
        body: &Diagram,
    ) -> Result<Vec<usize>> {
        let mut images = Vec::new();
        for (y, past) in self.pasts.iter().enumerate() {
            model.spend(1 + past.representatives.len())?;
            let gy = square.outside.positions()[y];
            for &(x, value) in &past.representatives {
                let gx = square.inside.positions()[x];
                images.push(old.pasts[gy].index(model, body, gx, value)?);
            }
        }
        Ok(images)
    }
}

pub(crate) struct Ran {
    pub(crate) along: Along,
    futures: Vec<Future>,
    evaluation: Vec<usize>,
}

impl Ran {
    pub(crate) fn build(
        model: &mut Model,
        along: Along,
        body: &Diagram,
    ) -> Result<(Diagram, Self)> {
        let target = model.poset(along.target).expect("named target site");
        let mut evaluation = vec![0; along.positions.len()];
        let mut futures = Vec::with_capacity(target.names.len());
        let mut cells = 2 * along.positions.len() + 1;
        let mut count = 0;
        for y in 0..target.names.len() {
            let mut points = Vec::new();
            model.spend(along.positions.len())?;
            for (x, &fx) in along.positions.iter().enumerate() {
                if y == fx {
                    evaluation[x] = points.len();
                }
                if le(&target, y, fx) {
                    points.push(x);
                }
            }
            let future = Future::build(model, body, points)?;
            count += future.sections.len();
            if count > set_family::MAX_OCCURRENCES {
                return Err(Error::plain(
                    "right Kan extension exceeds the limit of 4096 section occurrences",
                ));
            }
            cells += future.cells();
            set_family::check_cells(cells)?;
            futures.push(future);
        }
        let sizes = futures.iter().map(|future| future.sections.len()).collect();
        let diagram = Diagram::generated(model, along.target, sizes, |model, y, z| {
            futures[y].restrict(model, &futures[z])
        })?;
        model.allocate(cells)?;
        Ok((
            diagram,
            Self {
                along,
                futures,
                evaluation,
            },
        ))
    }

    fn transpose(
        &self,
        model: &mut Model,
        source: &Diagram,
        pulled: &Diagram,
        rows: Option<&[usize]>,
    ) -> Result<Vec<usize>> {
        let mut images = Vec::with_capacity(source.count());
        let mut choices = Vec::with_capacity(self.along.positions.len());
        for (y, future) in self.futures.iter().enumerate() {
            model.spend(1)?;
            for value in 0..source.size(y) {
                model.spend(1 + future.points.len())?;
                choices.clear();
                choices.extend(future.points.iter().map(|&x| {
                    let transported = source.apply(y, self.along.positions[x], value);
                    component(rows, pulled, x, transported)
                }));
                images.push(future.index(model, &choices)?);
            }
        }
        Ok(images)
    }

    fn untranspose(
        &self,
        model: &mut Model,
        source: &Diagram,
        rows: Option<&[usize]>,
    ) -> Result<Vec<usize>> {
        let mut images = Vec::new();
        for (x, &y) in self.along.positions.iter().enumerate() {
            model.spend(1 + source.size(y))?;
            for i in 0..source.size(y) {
                let j = component(rows, source, y, i);
                images.push(self.futures[y].sections[j][self.evaluation[x]]);
            }
        }
        Ok(images)
    }

    fn compare(&self, model: &mut Model, old: &Self, square: &Square) -> Result<Vec<usize>> {
        let mut images = Vec::new();
        for (y, future) in self.futures.iter().enumerate() {
            model.spend(1 + future.points.len())?;
            let old_future = &old.futures[square.outside.positions()[y]];
            let coordinates = future
                .points
                .iter()
                .map(|&x| old_future.position(model, square.inside.positions()[x]))
                .collect::<Result<Vec<_>>>()?;
            images.extend(old_future.select_into(model, future, &coordinates)?);
        }
        Ok(images)
    }
}

fn component(rows: Option<&[usize]>, domain: &Diagram, point: usize, value: usize) -> usize {
    rows.map_or(value, |rows| rows[domain.occurrence(point, value)])
}

pub(crate) struct Square {
    inside: Rc<Pullback>,
    outside: Rc<Pullback>,
}

impl Square {
    pub(crate) fn check(
        model: &mut Model,
        old: &Along,
        new: &Along,
        inside: Rc<Pullback>,
        outside: Rc<Pullback>,
    ) -> Result<Self> {
        if inside.family != old.body || inside.base != new.source || outside.base != new.target {
            return Err(Error::plain(
                "Kan comparison requires the two pullbacks of the given input and extension",
            ));
        }
        model.spend(1 + new.positions.len())?;
        if new
            .positions
            .iter()
            .zip(inside.positions())
            .any(|(&fx, &gx)| outside.positions()[fx] != old.positions[gx])
        {
            return Err(Error::plain(
                "Kan comparison requires a commuting square of monotone maps",
            ));
        }
        Ok(Self { inside, outside })
    }
}

pub(crate) enum Action {
    LanLift(Rc<Lan>, Rc<Diagram>, Rc<Diagram>),
    LanRestrict(Rc<Lan>, Rc<Diagram>, Rc<Diagram>),
    RanTranspose(Rc<Ran>, Rc<Diagram>, Rc<Diagram>),
    RanUntranspose(Rc<Ran>, Rc<Diagram>),
    LanComparison(Rc<Lan>, Rc<Lan>, Square, Rc<Diagram>),
    RanComparison(Rc<Ran>, Rc<Ran>, Square),
}

impl Action {
    pub(crate) fn evaluate(&self, model: &mut Model, rows: Option<&[usize]>) -> Result<Value> {
        let images = match self {
            Self::LanLift(lan, body, target) => lan.lift(model, body, target, rows)?,
            Self::LanRestrict(lan, extension, body) => {
                lan.restrict(model, extension, body, rows)?
            }
            Self::RanTranspose(ran, source, pulled) => {
                ran.transpose(model, source, pulled, rows)?
            }
            Self::RanUntranspose(ran, source) => ran.untranspose(model, source, rows)?,
            Self::LanComparison(new, old, square, body) => new.compare(model, old, square, body)?,
            Self::RanComparison(new, old, square) => new.compare(model, old, square)?,
        };
        set_family::check_cells(images.len())?;
        model.allocate(images.len())?;
        Ok(Value::SetNat(images.into()))
    }
}

fn lan(model: &Model, name: &str) -> Result<Rc<Lan>> {
    match model.set_origin(name) {
        Some(crate::dependent_sets::Origin::Lan(origin)) => Ok(origin),
        _ => Err(Error::plain(format!(
            "'{name}' is not a left Kan extension"
        ))),
    }
}

fn ran(model: &Model, name: &str) -> Result<Rc<Ran>> {
    match model.set_origin(name) {
        Some(crate::dependent_sets::Origin::Ran(origin)) => Ok(origin),
        _ => Err(Error::plain(format!(
            "'{name}' is not a right Kan extension"
        ))),
    }
}

fn family(model: &Model, name: &str) -> Rc<Diagram> {
    model.set_family_named(name).expect("registered set family")
}

impl Operation {
    pub(crate) fn plan(
        self,
        model: &mut Model,
        first: Rc<str>,
        second: Option<Rc<str>>,
        ends: Option<(Rc<str>, Rc<str>)>,
    ) -> Result<(TypeId, Action)> {
        let (source, target) = ends.unwrap_or_else(|| (Rc::from(""), Rc::from("")));
        let (from, to, action) = match self {
            Self::LanLift | Self::LanCounit => {
                let origin = lan(model, &first)?;
                let input = if matches!(self, Self::LanLift) {
                    if source != origin.along.body {
                        return Err(Error::plain(
                            "setLanLift requires a map from the left Kan extension's input",
                        ));
                    }
                    target
                } else {
                    origin.along.body.clone()
                };
                let pulled = pullback(model, &input)?;
                origin.along.check_pullback(model, &pulled)?;
                let body = family(model, &origin.along.body);
                let target = family(model, &pulled.family);
                (
                    first,
                    pulled.family.clone(),
                    Action::LanLift(origin, body, target),
                )
            }
            Self::LanRestrict | Self::LanUnit => {
                let (extension, pulled, output) = if matches!(self, Self::LanRestrict) {
                    (source, first, target)
                } else {
                    (first.clone(), second.expect("parsed unit pullback"), first)
                };
                let origin = lan(model, &extension)?;
                let lifted = pullback(model, &pulled)?;
                origin.along.check_pullback(model, &lifted)?;
                if lifted.family != output {
                    return Err(Error::plain(
                        "left Kan unit or restriction requires the pullback of the map's target",
                    ));
                }
                let body = family(model, &origin.along.body);
                let derived = family(model, &extension);
                (
                    origin.along.body.clone(),
                    pulled,
                    Action::LanRestrict(origin, derived, body),
                )
            }
            Self::RanTranspose | Self::RanUnit => {
                let origin = ran(model, &first)?;
                let input = if matches!(self, Self::RanTranspose) {
                    if target != origin.along.body {
                        return Err(Error::plain(
                            "setRanTranspose requires a map to the right Kan extension's input",
                        ));
                    }
                    source
                } else {
                    origin.along.body.clone()
                };
                let lifted = pullback(model, &input)?;
                origin.along.check_pullback(model, &lifted)?;
                let source = family(model, &lifted.family);
                let pulled = family(model, &input);
                (
                    lifted.family.clone(),
                    first,
                    Action::RanTranspose(origin, source, pulled),
                )
            }
            Self::RanUntranspose | Self::RanCounit => {
                let (extension, pulled, original) = if matches!(self, Self::RanUntranspose) {
                    (target, first, source)
                } else {
                    (
                        first.clone(),
                        second.expect("parsed counit pullback"),
                        first,
                    )
                };
                let origin = ran(model, &extension)?;
                let lifted = pullback(model, &pulled)?;
                origin.along.check_pullback(model, &lifted)?;
                if lifted.family != original {
                    return Err(Error::plain(
                        "right Kan counit or untranspose requires the pullback of the map's source",
                    ));
                }
                let source = family(model, &original);
                (
                    pulled,
                    origin.along.body.clone(),
                    Action::RanUntranspose(origin, source),
                )
            }
            Self::LanComparison => {
                let pulled = second.expect("parsed left comparison pullback");
                let new = lan(model, &first)?;
                let outside = pullback(model, &pulled)?;
                let old = lan(model, &outside.family)?;
                let inside = pullback(model, &new.along.body)?;
                let square = Square::check(model, &old.along, &new.along, inside, outside)?;
                let body = family(model, &old.along.body);
                (first, pulled, Action::LanComparison(new, old, square, body))
            }
            Self::RanComparison => {
                let pulled = second.expect("parsed right comparison pullback");
                let new = ran(model, &first)?;
                let outside = pullback(model, &pulled)?;
                let old = ran(model, &outside.family)?;
                let inside = pullback(model, &new.along.body)?;
                let square = Square::check(model, &old.along, &new.along, inside, outside)?;
                (pulled, first, Action::RanComparison(new, old, square))
            }
        };
        let ty = model.register_set_nat_type(from, to, !self.takes_input())?;
        Ok((ty, action))
    }
}
