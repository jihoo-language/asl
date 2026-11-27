use std::rc::Rc;

use crate::dependent_sets::{pullback, Future, Origin, Pullback};
use crate::error::{Error, Result};
use crate::model::{Model, Poset, TypeId, Value};
use crate::set_family::{self, Diagram};
use crate::syntax::Label;

#[derive(Clone, Copy)]
pub(crate) enum Side {
    Closed,
    Open,
}

#[derive(Clone, Copy)]
enum Point {
    Closed(usize),
    Open(usize),
}

#[derive(PartialEq, Eq)]
pub(crate) struct Embedding {
    pub(crate) domain: TypeId,
    pub(crate) images: Vec<usize>,
}

impl Embedding {
    pub(crate) fn from_map(domain: TypeId, values: &[Value]) -> Self {
        let images = values
            .iter()
            .map(|value| match value {
                Value::Point(index) => *index,
                _ => unreachable!("map of named sites"),
            })
            .collect();
        Self { domain, images }
    }

    fn matches(&self, pulled: &Pullback) -> bool {
        self.domain == pulled.base && self.images == pulled.positions()
    }
}

pub(crate) struct Recollement {
    base: TypeId,
    closed: Embedding,
    open: Embedding,
    points: Vec<Point>,
}

impl Recollement {
    pub(crate) fn build(
        model: &mut Model,
        base: TypeId,
        closed: Embedding,
        open: Embedding,
    ) -> Result<Self> {
        let site = model.poset(base).expect("named ambient site");
        let n = site.names.len();
        let mut points = vec![None; n];
        for (side, embedding) in [(Side::Closed, &closed), (Side::Open, &open)] {
            let stratum = model.poset(embedding.domain).expect("named stratum");
            model.spend(embedding.images.len() * (embedding.images.len() + 1))?;
            for (x, &image) in embedding.images.iter().enumerate() {
                let location = match side {
                    Side::Closed => Point::Closed(x),
                    Side::Open => Point::Open(x),
                };
                if points[image].replace(location).is_some() {
                    return Err(Error::plain(format!(
                        "recollement inclusions are not injective and disjoint at '{}'",
                        site.names[image]
                    )));
                }
                check_order(&stratum, &site, embedding, x, image)?;
            }
        }
        let points = points
            .into_iter()
            .zip(&site.names)
            .map(|(location, name)| {
                location.ok_or_else(|| {
                    Error::plain(format!("recollement omits ambient point '{name}'"))
                })
            })
            .collect::<Result<Vec<_>>>()?;
        model.spend(n * n)?;
        for &x in &open.images {
            for &y in &closed.images {
                if site.strict_less(x, y) {
                    return Err(Error::plain(format!(
                        "recollement open image is not upward closed: '{} < {}'",
                        site.names[x], site.names[y]
                    )));
                }
            }
        }
        model.allocate(n + closed.images.len() + open.images.len() + 1)?;
        Ok(Self {
            base,
            closed,
            open,
            points,
        })
    }

    fn same_inclusions(&self, model: &mut Model, other: &Self) -> Result<bool> {
        model.spend(self.points.len().min(other.points.len()) + 1)?;
        Ok(self.base == other.base && self.closed == other.closed && self.open == other.open)
    }

    pub(crate) fn restrict(
        &self,
        model: &mut Model,
        name: Rc<str>,
        body: &Diagram,
        side: Side,
    ) -> Result<(Diagram, Pullback)> {
        if body.base != self.base {
            return Err(Error::plain(
                "setrestrict requires a family on the ambient site",
            ));
        }
        let embedding = match side {
            Side::Closed => &self.closed,
            Side::Open => &self.open,
        };
        model.spend(embedding.images.len())?;
        let values: Vec<_> = embedding.images.iter().map(|&i| Value::Point(i)).collect();
        Pullback::build(model, embedding.domain, name, body, &values)
    }

    fn check_strata(&self, model: &mut Model, closed: &Pullback, open: &Pullback) -> Result<()> {
        model.spend(self.points.len() + 1)?;
        if !self.closed.matches(closed) || !self.open.matches(open) || closed.family != open.family
        {
            return Err(Error::plain(
                "attachment requires the closed and open restrictions of the very same ambient family",
            ));
        }
        Ok(())
    }
}

fn check_order(
    stratum: &Poset,
    site: &Poset,
    embedding: &Embedding,
    x: usize,
    image: usize,
) -> Result<()> {
    for (y, &other) in embedding.images.iter().enumerate() {
        if stratum.strict_less(x, y) != site.strict_less(image, other) {
            return Err(Error::plain(format!(
                "recollement inclusion must reflect the order between '{}' and '{}'",
                stratum.names[x], stratum.names[y]
            )));
        }
    }
    Ok(())
}

pub(crate) struct Matching {
    recollement: Rc<Recollement>,
    body: Rc<str>,
    futures: Vec<Future>,
}

impl Matching {
    pub(crate) fn build(
        model: &mut Model,
        recollement: Rc<Recollement>,
        name: Rc<str>,
        body: &Diagram,
    ) -> Result<(Diagram, Self)> {
        if body.base != recollement.open.domain {
            return Err(Error::plain(
                "setgluelimit requires a family on the open stratum",
            ));
        }
        let ambient = model.poset(recollement.base).expect("named ambient site");
        let mut futures = Vec::with_capacity(recollement.closed.images.len());
        let mut cells = 1;
        let mut count = 0;
        for &x in &recollement.closed.images {
            model.spend(1 + recollement.open.images.len())?;
            let points = recollement
                .open
                .images
                .iter()
                .enumerate()
                .filter_map(|(u, &y)| ambient.strict_less(x, y).then_some(u))
                .collect();
            let future = Future::build(model, body, points)?;
            cells += future.cells();
            count += future.sections.len();
            if count > set_family::MAX_OCCURRENCES {
                return Err(Error::plain(
                    "set gluing limit exceeds the limit of 4096 section occurrences",
                ));
            }
            set_family::check_cells(cells)?;
            futures.push(future);
        }
        let sizes = futures.iter().map(|future| future.sections.len()).collect();
        let diagram =
            Diagram::generated(model, recollement.closed.domain, sizes, |model, x, y| {
                futures[x].restrict(model, &futures[y])
            })?;
        model.allocate(cells)?;
        Ok((
            diagram,
            Self {
                recollement,
                body: name,
                futures,
            },
        ))
    }

    fn attachment(
        &self,
        model: &mut Model,
        closed_name: &str,
    ) -> Result<(Rc<Diagram>, Vec<usize>)> {
        let closed = pullback(model, closed_name)?;
        let open = pullback(model, &self.body)?;
        self.recollement.check_strata(model, &closed, &open)?;
        let original = family(model, &closed.family);
        if original.base != self.recollement.base {
            return Err(Error::plain(
                "attachment requires the original family on the named ambient site",
            ));
        }
        let closed = family(model, closed_name);
        let mut rows = Vec::with_capacity(closed.count());
        let mut values = Vec::new();
        for (c, &x) in self.recollement.closed.images.iter().enumerate() {
            let future = &self.futures[c];
            for i in 0..closed.size(c) {
                model.spend(future.points.len() + 1)?;
                values.clear();
                values.extend(
                    future
                        .points
                        .iter()
                        .map(|&u| original.apply(x, self.recollement.open.images[u], i)),
                );
                rows.push(future.index(model, &values)?);
            }
        }
        Ok((closed, rows))
    }

    fn check_same_site(&self, model: &mut Model, target: &Self) -> Result<()> {
        if !self
            .recollement
            .same_inclusions(model, &target.recollement)?
        {
            return Err(Error::plain(
                "set gluing requires the same closed and open inclusions for both diagrams",
            ));
        }
        Ok(())
    }
}

struct SectionMap<'a> {
    source: &'a Matching,
    target: &'a Matching,
    open: &'a Diagram,
    map: &'a [usize],
    choices: Vec<usize>,
}

impl<'a> SectionMap<'a> {
    fn new(
        source: &'a Matching,
        target: &'a Matching,
        open: &'a Diagram,
        map: &'a [usize],
    ) -> Self {
        Self {
            source,
            target,
            open,
            map,
            choices: Vec::new(),
        }
    }

    fn section(&mut self, model: &mut Model, c: usize, index: usize) -> Result<usize> {
        let future = &self.source.futures[c];
        let section = &future.sections[index];
        model.spend(1 + section.len())?;
        self.choices.clear();
        self.choices.extend(
            future
                .points
                .iter()
                .zip(section)
                .map(|(&u, &i)| self.map[self.open.occurrence(u, i)]),
        );
        self.target.futures[c].index(model, &self.choices)
    }

    fn rows(mut self, model: &mut Model) -> Result<Vec<usize>> {
        let mut rows = Vec::new();
        model.spend(self.source.futures.len())?;
        for (c, future) in self.source.futures.iter().enumerate() {
            for index in 0..future.sections.len() {
                rows.push(self.section(model, c, index)?);
            }
        }
        Ok(rows)
    }
}

pub(crate) struct Glued {
    closed: Rc<str>,
    matching: Rc<Matching>,
    alpha: Rc<[usize]>,
}

impl Glued {
    pub(crate) fn build(
        model: &mut Model,
        recollement: Rc<Recollement>,
        closed_name: Rc<str>,
        limit_name: &str,
        alpha: Rc<[usize]>,
    ) -> Result<(Diagram, Self)> {
        let matching = matching(model, limit_name)?;
        if !recollement.same_inclusions(model, &matching.recollement)? {
            return Err(Error::plain(
                "setglue attachment must target a limit of this recollement",
            ));
        }
        let closed = family(model, &closed_name);
        let open = family(model, &matching.body);
        let sizes = recollement
            .points
            .iter()
            .map(|point| match *point {
                Point::Closed(c) => closed.size(c),
                Point::Open(u) => open.size(u),
            })
            .collect();
        let diagram = Diagram::generated(model, recollement.base, sizes, |model, x, y| {
            match (recollement.points[x], recollement.points[y]) {
                (Point::Closed(c), Point::Closed(d)) => {
                    model.spend(closed.size(c))?;
                    Ok((0..closed.size(c)).map(|i| closed.apply(c, d, i)).collect())
                }
                (Point::Open(u), Point::Open(v)) => {
                    model.spend(open.size(u))?;
                    Ok((0..open.size(u)).map(|i| open.apply(u, v, i)).collect())
                }
                (Point::Closed(c), Point::Open(u)) => {
                    let future = &matching.futures[c];
                    let coordinate = future.position(model, u)?;
                    model.spend(closed.size(c))?;
                    Ok(closed
                        .components(c, &alpha)
                        .iter()
                        .map(|&index| future.sections[index][coordinate])
                        .collect())
                }
                (Point::Open(_), Point::Closed(_)) => unreachable!("open image is upward closed"),
            }
        })?;
        model.allocate(1)?;
        Ok((
            diagram,
            Self {
                closed: closed_name,
                matching,
                alpha,
            },
        ))
    }

    fn recover(&self, model: &mut Model, original_name: &str) -> Result<Rc<Diagram>> {
        let restriction = pullback(model, &self.closed)?;
        if restriction.family.as_ref() != original_name {
            return Err(Error::plain(
                "reconstruction requires the original ambient family's restrictions",
            ));
        }
        let (_, canonical) = self.matching.attachment(model, &self.closed)?;
        model.spend(canonical.len())?;
        if canonical.as_slice() != self.alpha.as_ref() {
            return Err(Error::plain(
                "reconstruction requires the original canonical attachment",
            ));
        }
        Ok(family(model, original_name))
    }
}

fn family(model: &Model, name: &str) -> Rc<Diagram> {
    model.set_family_named(name).expect("registered family")
}

fn matching(model: &Model, name: &str) -> Result<Rc<Matching>> {
    match model.set_origin(name) {
        Some(Origin::GlueLimit(origin)) => Ok(origin),
        _ => Err(Error::plain(format!("'{name}' is not a set gluing limit"))),
    }
}

fn glued(model: &Model, name: &str) -> Result<Rc<Glued>> {
    match model.set_origin(name) {
        Some(Origin::Glued(origin)) => Ok(origin),
        _ => Err(Error::plain(format!("'{name}' is not a glued set family"))),
    }
}

pub(crate) type Ends = (Rc<str>, Rc<str>);

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Attach,
    LimitMap,
    GlueMap,
    Closed,
    Open,
    Recover,
    Forget,
}

impl Operation {
    pub(crate) fn detect(text: &str) -> Option<Self> {
        Some(match text {
            "setAttach" => Self::Attach,
            "setGlueLimitMap" => Self::LimitMap,
            "setGlueMap" => Self::GlueMap,
            "setGlueClosed" => Self::Closed,
            "setGlueOpen" => Self::Open,
            "setGlueRecover" => Self::Recover,
            "setGlueForget" => Self::Forget,
            _ => return None,
        })
    }

    pub(crate) fn arity(self) -> usize {
        match self {
            Self::Attach | Self::Recover | Self::Forget => 0,
            Self::LimitMap | Self::Closed | Self::Open => 1,
            Self::GlueMap => 2,
        }
    }

    pub(crate) fn plan(
        self,
        model: &mut Model,
        [first, second]: [Label; 2],
        input: &[Ends],
        known: Option<&[Rc<[usize]>]>,
    ) -> Result<(TypeId, Action)> {
        let first: Rc<str> = first.text.into();
        let second: Rc<str> = second.text.into();
        let (from, to, action) = match self {
            Self::Attach => {
                let matching = matching(model, &second)?;
                let (source, rows) = matching.attachment(model, &first)?;
                let target = family(model, &second);
                let value = source.certify(model, &target, rows)?;
                (first, second, Action::Constant(value))
            }
            Self::Recover | Self::Forget => {
                let origin = glued(model, &first)?;
                let original = origin.recover(model, &second)?;
                let value = original.identity(model)?;
                if matches!(self, Self::Recover) {
                    (second, first, Action::Constant(value))
                } else {
                    (first, second, Action::Constant(value))
                }
            }
            Self::LimitMap => {
                let source = matching(model, &first)?;
                let target = matching(model, &second)?;
                source.check_same_site(model, &target)?;
                check_ends(&input[0], &source.body, &target.body)?;
                let open = family(model, &source.body);
                (first, second, Action::LimitMap(source, target, open))
            }
            Self::GlueMap => {
                let source = glued(model, &first)?;
                let target = glued(model, &second)?;
                source.matching.check_same_site(model, &target.matching)?;
                check_ends(&input[0], &source.closed, &target.closed)?;
                check_ends(&input[1], &source.matching.body, &target.matching.body)?;
                let ambient = family(model, &first);
                let ambient_target = family(model, &second);
                (
                    first,
                    second,
                    Action::GlueMap(Box::new(GlueMap {
                        source,
                        target,
                        ambient,
                        ambient_target,
                    })),
                )
            }
            Self::Closed | Self::Open => {
                let source = glued(model, &first)?;
                let target = glued(model, &second)?;
                source.matching.check_same_site(model, &target.matching)?;
                check_ends(&input[0], &first, &second)?;
                let (side, from, to) = if matches!(self, Self::Closed) {
                    (Side::Closed, source.closed.clone(), target.closed.clone())
                } else {
                    (
                        Side::Open,
                        source.matching.body.clone(),
                        target.matching.body.clone(),
                    )
                };
                let ambient = family(model, &first);
                (
                    from,
                    to,
                    Action::Part(source.matching.recollement.clone(), ambient, side),
                )
            }
        };
        // A supplied map proves inhabitance directly after checking its square.
        let action = if let Some(known) = known.filter(|maps| !maps.is_empty()) {
            Action::Constant(action.evaluate(model, known)?)
        } else {
            action
        };
        let ty = model.register_set_nat_type(from, to, matches!(action, Action::Constant(_)))?;
        Ok((ty, action))
    }
}

fn check_ends(input: &Ends, source: &str, target: &str) -> Result<()> {
    if input.0.as_ref() != source || input.1.as_ref() != target {
        return Err(Error::plain(format!(
            "set gluing expected SetNat({source}, {target})"
        )));
    }
    Ok(())
}

pub(crate) enum Action {
    Constant(Value),
    LimitMap(Rc<Matching>, Rc<Matching>, Rc<Diagram>),
    GlueMap(Box<GlueMap>),
    Part(Rc<Recollement>, Rc<Diagram>, Side),
}

pub(crate) struct GlueMap {
    source: Rc<Glued>,
    target: Rc<Glued>,
    ambient: Rc<Diagram>,
    ambient_target: Rc<Diagram>,
}

impl Action {
    pub(crate) fn evaluate(&self, model: &mut Model, inputs: &[Rc<[usize]>]) -> Result<Value> {
        let rows = match self {
            Self::Constant(value) => return Ok(value.clone()),
            Self::LimitMap(source, target, open) => {
                SectionMap::new(source, target, open, &inputs[0]).rows(model)?
            }
            Self::Part(recollement, ambient, side) => {
                let embedding = match side {
                    Side::Closed => &recollement.closed,
                    Side::Open => &recollement.open,
                };
                model.spend(embedding.images.len())?;
                embedding
                    .images
                    .iter()
                    .flat_map(|&x| ambient.components(x, &inputs[0]).iter().copied())
                    .collect()
            }
            Self::GlueMap(operation) => return operation.evaluate(model, &inputs[0], &inputs[1]),
        };
        model.spend(rows.len())?;
        model.allocate(rows.len())?;
        Ok(Value::SetNat(rows.into()))
    }
}

impl GlueMap {
    fn evaluate(&self, model: &mut Model, a: &[usize], b: &[usize]) -> Result<Value> {
        let closed = family(model, &self.source.closed);
        let closed_target = family(model, &self.target.closed);
        let open = family(model, &self.source.matching.body);
        let mut map = SectionMap::new(&self.source.matching, &self.target.matching, &open, b);
        for c in 0..closed.site.names.len() {
            for (i, &chosen) in closed.components(c, &self.source.alpha).iter().enumerate() {
                let mapped = map.section(model, c, chosen)?;
                let moved_closed = a[closed.occurrence(c, i)];
                let expected = self.target.alpha[closed_target.occurrence(c, moved_closed)];
                if mapped != expected {
                    return Err(Error::plain(format!(
                        "set gluing attachment square fails at closed point '{}' and element {i}",
                        closed.site.names[c]
                    )));
                }
            }
        }
        model.spend(self.ambient.count() + self.source.matching.recollement.points.len())?;
        let rows = self
            .source
            .matching
            .recollement
            .points
            .iter()
            .flat_map(|point| match *point {
                Point::Closed(c) => closed.components(c, a),
                Point::Open(u) => open.components(u, b),
            })
            .copied()
            .collect();
        self.ambient.certify(model, &self.ambient_target, rows)
    }
}
