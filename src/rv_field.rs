use crate::binary_field::{BinaryField, Embedding};
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::rv_hom::{self, Grade, Support};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Frobenius(usize),
    Extend(usize, usize),
    Descend(usize, usize),
    Trace(usize, usize),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Frobenius(_) => "rvFrob",
            Self::Extend(..) => "rvExtend",
            Self::Descend(..) => "rvDescend",
            Self::Trace(..) => "rvTrace",
        }
    }
}

pub(crate) enum Map {
    Frobenius(BinaryField, usize),
    Extend(Embedding),
    Descend(Embedding),
    Trace(Embedding),
}

impl Map {
    fn new(model: &mut Model, field: BinaryField, operation: Operation) -> Result<Self> {
        let map = match operation {
            Operation::Frobenius(iterations) => Self::Frobenius(field, iterations),
            Operation::Extend(modulus, root) => {
                let large = BinaryField::new(model, modulus)?;
                Self::Extend(Embedding::new(model, field, large, root)?)
            }
            Operation::Descend(modulus, root) | Operation::Trace(modulus, root) => {
                let small = BinaryField::new(model, modulus)?;
                let inclusion = Embedding::new(model, small, field, root)?;
                match operation {
                    Operation::Trace(..) => Self::Trace(inclusion),
                    _ => Self::Descend(inclusion),
                }
            }
        };
        Ok(map)
    }

    fn output_field(&self) -> BinaryField {
        match self {
            Self::Frobenius(field, _) => *field,
            Self::Extend(inclusion) => inclusion.large(),
            Self::Descend(inclusion) | Self::Trace(inclusion) => inclusion.small(),
        }
    }

    fn point_coordinate(
        &self,
        model: &mut Model,
        value: usize,
        operation: Operation,
    ) -> Result<usize> {
        match self {
            Self::Descend(inclusion) | Self::Trace(inclusion) => {
                inclusion.preimage(model, value)?.ok_or_else(|| {
                    Error::plain(format!(
                        "{} requires the affine point to lie in the specified subfield",
                        operation.name()
                    ))
                })
            }
            _ => self.coefficient(model, value),
        }
    }

    fn coefficient(&self, model: &mut Model, value: usize) -> Result<usize> {
        match self {
            Self::Frobenius(field, iterations) => field.frobenius(model, value, *iterations),
            Self::Extend(inclusion) => inclusion.image(model, value),
            Self::Descend(inclusion) => inclusion.preimage(model, value)?.ok_or_else(|| {
                Error::plain(
                    "rvDescend requires every map coefficient to lie in the specified subfield",
                )
            }),
            Self::Trace(inclusion) => inclusion.trace(model, value),
        }
    }

    pub(crate) fn apply(&self, model: &mut Model, support: &Support) -> Result<Value> {
        // The Margolis reductions use ordered leading coordinates. A field
        // embedding preserves their supports and commutes with normalizing
        // each pivot, so the induced Hom basis uses the same coordinates.
        model.spend(support.len())?;
        let mut output = Vec::new();
        for &(index, value) in support.iter() {
            let coefficient = self.coefficient(model, value)?;
            if coefficient != 0 {
                output.push((index, coefficient));
            }
        }
        rv_hom::persist(model, output)
    }
}

pub(crate) fn plan(model: &mut Model, grade: Grade, operation: Operation) -> Result<(TypeId, Map)> {
    let map = Map::new(model, grade.field, operation)?;
    model.spend(grade.point.len())?;
    let point = grade
        .point
        .iter()
        .map(|&value| map.point_coordinate(model, value, operation))
        .collect::<Result<Vec<_>>>()?;
    model.allocate(1 + point.len())?;
    let ty = model.rv_hom_type(Grade {
        field: map.output_field(),
        point: point.into(),
        ..grade
    })?;
    Ok((ty, map))
}
