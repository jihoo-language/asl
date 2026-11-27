use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::set_family::{self, Diagram};

struct Axis {
    point: usize,
    radix: usize,
    stride: usize,
}

struct Future {
    axes: Vec<Axis>,
    size: usize,
}

impl Future {
    fn new(family: &Diagram, points: impl DoubleEndedIterator<Item = usize>) -> Self {
        let mut axes = Vec::new();
        let mut size = 1usize;
        for point in points.rev() {
            let radix = family.size(point);
            axes.push(Axis {
                point,
                radix,
                stride: size,
            });
            size = size.saturating_mul(radix);
        }
        axes.reverse();
        Self { axes, size }
    }

    fn tuple(&self, mut coordinate: impl FnMut(&Axis) -> usize) -> usize {
        self.axes
            .iter()
            .map(|axis| axis.stride * coordinate(axis))
            .sum()
    }

    fn project(&self, axis: &Axis, code: usize) -> usize {
        code / axis.stride % axis.radix
    }

    fn unit(&self, family: &Diagram, point: usize, value: usize) -> usize {
        self.tuple(|axis| family.apply(point, axis.point, value))
    }
}

struct Shape {
    futures: Vec<Future>,
}

impl Shape {
    fn build(model: &mut Model, family: &Diagram) -> Result<Self> {
        let n = family.site.names.len();
        model.spend(2 * n * n)?;
        let maxima: Vec<_> = (0..n)
            .filter(|&x| (0..n).all(|y| !family.site.strict_less(x, y)))
            .collect();
        let futures = (0..n)
            .map(|x| {
                Future::new(
                    family,
                    maxima
                        .iter()
                        .copied()
                        .filter(|&y| x == y || family.site.strict_less(x, y)),
                )
            })
            .collect();
        Ok(Self { futures })
    }

    fn cells(&self) -> usize {
        self.futures
            .iter()
            .map(|future| 3 * future.axes.len() + 2)
            .sum()
    }
}

pub(crate) struct Sheafification {
    pub(crate) family: Rc<str>,
    shape: Shape,
    unit: Rc<[usize]>,
}

impl Sheafification {
    pub(crate) fn build(
        model: &mut Model,
        name: Rc<str>,
        family: &Diagram,
    ) -> Result<(Diagram, Self)> {
        let shape = Shape::build(model, family)?;
        let sizes: Vec<_> = shape.futures.iter().map(|future| future.size).collect();
        if sizes.iter().copied().fold(0usize, usize::saturating_add) > set_family::MAX_OCCURRENCES {
            return Err(Error::plain(
                "double-negation sheafification exceeds the limit of 4096 occurrences",
            ));
        }
        let result = Diagram::generated(model, family.base, sizes, |model, x, y| {
            let source = &shape.futures[x];
            let target = &shape.futures[y];
            model.spend(
                source.axes.len() + target.axes.len() + source.size * (target.axes.len() + 1),
            )?;
            let mut scan = 0;
            let mut matching = Vec::new();
            for axis in &target.axes {
                while source.axes[scan].point != axis.point {
                    scan += 1;
                }
                matching.push(&source.axes[scan]);
            }
            Ok((0..source.size)
                .map(|i| {
                    target
                        .axes
                        .iter()
                        .zip(&matching)
                        .map(|(axis, matched)| axis.stride * source.project(matched, i))
                        .sum()
                })
                .collect())
        })?;
        let mut unit = Vec::with_capacity(family.count());
        for (x, future) in shape.futures.iter().enumerate() {
            model.spend(family.size(x) * (future.axes.len() + 1))?;
            unit.extend((0..family.size(x)).map(|i| future.unit(family, x, i)));
        }
        model.allocate(shape.cells() + unit.len() + 1)?;
        Ok((
            result,
            Self {
                family: name,
                shape,
                unit: unit.into(),
            },
        ))
    }

    pub(crate) fn unit(&self) -> Value {
        Value::SetNat(self.unit.clone())
    }

    pub(crate) fn extend(
        &self,
        model: &mut Model,
        source: &Diagram,
        target: &Diagram,
        certificate: &Certificate,
        map: &[usize],
    ) -> Result<Value> {
        let count = self.shape.futures.iter().map(|future| future.size).sum();
        model.allocate(count)?;
        let mut rows = Vec::with_capacity(count);
        for (x, future) in self.shape.futures.iter().enumerate() {
            let target_future = &certificate.shape.futures[x];
            model.spend((future.axes.len() + 1) * (future.size + 1))?;
            for i in 0..future.size {
                let code = target_future
                    .axes
                    .iter()
                    .zip(&future.axes)
                    .map(|(axis, from)| {
                        let value = future.project(from, i);
                        let output = map[source.occurrence(axis.point, value)];
                        axis.stride * output
                    })
                    .sum();
                rows.push(certificate.inverse[target.occurrence(x, code)]);
            }
        }
        Ok(Value::SetNat(rows.into()))
    }
}

pub(crate) struct Certificate {
    shape: Shape,
    inverse: Rc<[usize]>,
}

impl Certificate {
    pub(crate) fn build(model: &mut Model, family: &Diagram) -> Result<Option<Self>> {
        let shape = Shape::build(model, family)?;
        let mut inverse = vec![None; family.count()];
        for (x, future) in shape.futures.iter().enumerate() {
            if family.size(x) != future.size {
                return Ok(None);
            }
            model.spend(family.size(x) * (future.axes.len() + 1))?;
            for i in 0..family.size(x) {
                let code = future.unit(family, x, i);
                let position = family.occurrence(x, code);
                if inverse[position].replace(i).is_some() {
                    return Ok(None);
                }
            }
        }
        model.allocate(shape.cells() + inverse.len() + 1)?;
        let inverse = inverse
            .into_iter()
            .map(|i| i.expect("finite bijection"))
            .collect();
        Ok(Some(Self { shape, inverse }))
    }

    pub(crate) fn inverse(&self) -> Value {
        Value::SetNat(self.inverse.clone())
    }
}
