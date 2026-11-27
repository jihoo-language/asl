use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::set_family::{self, Diagram, Elements};

#[derive(Clone)]
pub(crate) enum Origin {
    Pullback(Rc<Pullback>),
    Sum(Rc<Sum>),
    Product(Rc<Product>),
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    PullNatural,
    Evaluation,
    Transpose,
    Untranspose,
    SumLift,
    SumRestrict,
}

impl Operation {
    pub(crate) fn detect(name: &str) -> Option<Self> {
        Some(match name {
            "setPullNat" => Self::PullNatural,
            "setEval" => Self::Evaluation,
            "setTranspose" => Self::Transpose,
            "setUntranspose" => Self::Untranspose,
            "setSumLift" => Self::SumLift,
            "setSumRestrict" => Self::SumRestrict,
            _ => return None,
        })
    }
}

pub(crate) struct Pullback {
    pub(crate) family: Rc<str>,
    pub(crate) base: TypeId,
    positions: Vec<usize>,
}

impl Pullback {
    pub(crate) fn build(
        model: &mut Model,
        base: TypeId,
        family: Rc<str>,
        source: &Diagram,
        map: &[Value],
    ) -> Result<(Diagram, Self)> {
        let positions: Vec<_> = map
            .iter()
            .map(|value| {
                let Value::Point(index) = *value else {
                    unreachable!("map of named sites")
                };
                index
            })
            .collect();
        let sizes = positions.iter().map(|&x| source.size(x)).collect();
        let diagram = Diagram::generated(model, base, sizes, |model, x, y| {
            let x = positions[x];
            let y = positions[y];
            model.spend(source.size(x))?;
            Ok((0..source.size(x)).map(|i| source.apply(x, y, i)).collect())
        })?;
        model.allocate(positions.len() + 1)?;
        Ok((
            diagram,
            Self {
                family,
                base,
                positions,
            },
        ))
    }

    fn matches_projection(&self, base: TypeId, elements: &Elements) -> bool {
        self.base == base
            && self
                .positions
                .iter()
                .copied()
                .eq(elements.points.iter().map(|&(x, _)| x))
    }

    pub(crate) fn same_map(&self, other: &Self) -> bool {
        self.base == other.base && self.positions == other.positions
    }

    pub(crate) fn natural(
        &self,
        model: &mut Model,
        source: &Diagram,
        rows: &[usize],
    ) -> Result<Value> {
        let count = self
            .positions
            .iter()
            .map(|&x| source.size(x))
            .sum::<usize>();
        model.spend(count + self.positions.len())?;
        model.allocate(count)?;
        Ok(Value::SetNat(
            self.positions
                .iter()
                .flat_map(|&x| source.components(x, rows).iter().copied())
                .collect(),
        ))
    }
}

pub(crate) struct Projection {
    pub(crate) total: TypeId,
    pub(crate) elements: Rc<Elements>,
    pub(crate) body: Rc<str>,
}

impl Projection {
    pub(crate) fn check(&self, pullback: &Pullback) -> Result<()> {
        if !pullback.matches_projection(self.total, &self.elements) {
            return Err(Error::plain(
                "expected a set family pulled back along this elements projection",
            ));
        }
        Ok(())
    }
}

pub(crate) struct Sum {
    pub(crate) projection: Projection,
    pairs: Vec<Vec<(usize, usize)>>,
    starts: Vec<usize>,
}

impl Sum {
    pub(crate) fn build(
        model: &mut Model,
        projection: Projection,
        base: &Diagram,
        body: &Diagram,
    ) -> Result<(Diagram, Self)> {
        model.spend(body.count() + projection.elements.points.len())?;
        let mut pairs = vec![Vec::new(); base.site.names.len()];
        let mut starts = Vec::new();
        for (point, &(x, _)) in projection.elements.points.iter().enumerate() {
            starts.push(pairs[x].len());
            pairs[x].extend((0..body.size(point)).map(|value| (point, value)));
        }
        let sizes = pairs.iter().map(Vec::len).collect();
        let diagram = Diagram::generated(model, base.base, sizes, |model, x, y| {
            model.spend(pairs[x].len())?;
            Ok(pairs[x]
                .iter()
                .map(|&(point, value)| {
                    let (_, i) = projection.elements.points[point];
                    let target = base.occurrence(y, base.apply(x, y, i));
                    starts[target] + body.apply(point, target, value)
                })
                .collect())
        })?;
        model.allocate(2 * body.count() + starts.len() + pairs.len() + 1)?;
        Ok((
            diagram,
            Self {
                projection,
                pairs,
                starts,
            },
        ))
    }

    pub(crate) fn lift(&self, model: &mut Model, body: &Diagram, rows: &[usize]) -> Result<Value> {
        model.spend(body.count() + self.pairs.len())?;
        model.allocate(body.count())?;
        Ok(Value::SetNat(
            self.pairs
                .iter()
                .flatten()
                .map(|&(point, value)| rows[body.occurrence(point, value)])
                .collect(),
        ))
    }

    pub(crate) fn restrict(
        &self,
        model: &mut Model,
        sum: &Diagram,
        body: &Diagram,
        rows: &[usize],
    ) -> Result<Value> {
        model.spend(body.count() + self.starts.len())?;
        model.allocate(body.count())?;
        let images = self
            .projection
            .elements
            .points
            .iter()
            .enumerate()
            .flat_map(|(point, &(x, _))| {
                let start = sum.occurrence(x, self.starts[point]);
                rows[start..start + body.size(point)].iter().copied()
            })
            .collect();
        Ok(Value::SetNat(images))
    }
}

struct Future {
    points: Vec<usize>,
    sections: Vec<Vec<usize>>,
}

impl Future {
    fn index(&self, model: &mut Model, choices: &[usize]) -> Result<usize> {
        let search_steps = self.sections.len().max(1).ilog2() as usize + 1;
        model.spend((choices.len() + 1) * search_steps)?;
        self.sections
            .binary_search_by(|section| section.as_slice().cmp(choices))
            .map_err(|_| Error::plain("dependent product choices do not form a future section"))
    }
}

pub(crate) struct Product {
    pub(crate) projection: Projection,
    futures: Vec<Future>,
    // Coordinate of each point in the future of its own projection.
    evaluation: Vec<usize>,
}

impl Product {
    pub(crate) fn build(
        model: &mut Model,
        projection: Projection,
        base: &Diagram,
        body: &Diagram,
    ) -> Result<(Diagram, Self)> {
        let points = &projection.elements.points;
        let mut evaluation = vec![0; points.len()];
        let mut futures = Vec::new();
        let mut cells = points.len() + 1;
        let mut count = 0;
        for x in 0..base.site.names.len() {
            let mut future = Vec::new();
            model.spend(points.len())?;
            for (point, &(y, _)) in points.iter().enumerate() {
                if x == y {
                    evaluation[point] = future.len();
                }
                if x == y || base.site.strict_less(x, y) {
                    future.push(point);
                }
            }
            let sections = body.sections_on(model, &future)?;
            count += sections.len();
            if count > set_family::MAX_OCCURRENCES {
                return Err(Error::plain(
                    "dependent product exceeds the limit of 4096 section occurrences",
                ));
            }
            cells += future.len() + sections.len() * (future.len() + 1) + 1;
            set_family::check_cells(cells)?;
            futures.push(Future {
                points: future,
                sections,
            });
        }
        let sizes = futures.iter().map(|future| future.sections.len()).collect();
        let diagram = Diagram::generated(model, base.base, sizes, |model, x, y| {
            let source = &futures[x];
            let target = &futures[y];
            let mut retained = Vec::with_capacity(target.points.len());
            let mut scan = 0;
            model.spend(source.points.len() + target.points.len())?;
            for point in &target.points {
                while &source.points[scan] != point {
                    scan += 1;
                }
                retained.push(scan);
            }
            let mut choices = Vec::with_capacity(retained.len());
            source
                .sections
                .iter()
                .map(|section| {
                    model.spend(retained.len() + 1)?;
                    choices.clear();
                    choices.extend(retained.iter().map(|&position| section[position]));
                    target.index(model, &choices)
                })
                .collect()
        })?;
        model.allocate(cells)?;
        Ok((
            diagram,
            Self {
                projection,
                futures,
                evaluation,
            },
        ))
    }

    pub(crate) fn evaluate(&self, model: &mut Model) -> Result<Value> {
        let rows: Vec<_> = self
            .projection
            .elements
            .points
            .iter()
            .enumerate()
            .flat_map(|(point, &(x, _))| {
                self.futures[x]
                    .sections
                    .iter()
                    .map(move |section| section[self.evaluation[point]])
            })
            .collect();
        model.spend(rows.len())?;
        model.allocate(rows.len())?;
        Ok(Value::SetNat(rows.into()))
    }

    pub(crate) fn transpose(
        &self,
        model: &mut Model,
        source: &Diagram,
        pulled: &Diagram,
        rows: &[usize],
    ) -> Result<Value> {
        let points = &self.projection.elements.points;
        let mut images = Vec::with_capacity(source.count());
        let mut choices = Vec::with_capacity(points.len());
        for (x, future) in self.futures.iter().enumerate() {
            for value in 0..source.size(x) {
                model.spend(future.points.len() + 1)?;
                choices.clear();
                for &point in &future.points {
                    let (y, _) = points[point];
                    let image = source.apply(x, y, value);
                    choices.push(rows[pulled.occurrence(point, image)]);
                }
                images.push(future.index(model, &choices)?);
            }
        }
        model.allocate(images.len())?;
        Ok(Value::SetNat(images.into()))
    }

    pub(crate) fn untranspose(
        &self,
        model: &mut Model,
        source: &Diagram,
        pulled: &Diagram,
        rows: &[usize],
    ) -> Result<Value> {
        model.spend(pulled.count() + self.evaluation.len())?;
        model.allocate(pulled.count())?;
        Ok(Value::SetNat(
            self.projection
                .elements
                .points
                .iter()
                .enumerate()
                .flat_map(|(point, &(x, _))| {
                    source
                        .components(x, rows)
                        .iter()
                        .map(move |&image| self.futures[x].sections[image][self.evaluation[point]])
                })
                .collect(),
        ))
    }
}
