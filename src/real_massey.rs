use std::cmp::Ordering;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::linear::Basis;
use crate::model::{Model, Value};
use crate::motivic::Grade;
use crate::real_ext;
use crate::real_motivic::{self, Term};

pub(crate) type Input<'a> = (&'a [Term], bool);

pub(crate) struct Layout {
    products: [Grade; 2],
    pub(crate) target: Grade,
}

impl Layout {
    pub(crate) fn new([a, b, c]: [Grade; 3]) -> Result<Self> {
        let ab = a.product(b)?;
        let bc = b.product(c)?;
        let target = ab.product(c)?.previous().ok_or_else(|| {
            Error::plain("real-motivic Massey products require positive total cobar length")
        })?;
        Ok(Self {
            products: [ab, bc],
            target,
        })
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        inputs: [Input<'_>; 3],
    ) -> Result<Value> {
        let terms = self.defining_cycle(model, inputs)?;
        real_motivic::persist(model, terms)
    }

    fn defining_cycle(&self, model: &mut Model, [a, b, c]: [Input<'_>; 3]) -> Result<Vec<Term>> {
        for input in [a, b, c] {
            require_cycle(model, input)?;
        }
        let [ab, bc] = self.products;
        let u = product_null(model, ab, a.0, b.0, "A cup B")?;
        let v = product_null(model, bc, b.0, c.0, "B cup C")?;
        let left = real_motivic::concat_support(model, &u, c.0)?;
        let right = real_motivic::concat_support(model, a.0, &v)?;
        real_motivic::add(model, &left, &right)
    }

    pub(crate) fn coset(&self, model: &mut Model, inputs: [Input<'_>; 3]) -> Result<Coset> {
        let terms = self.defining_cycle(model, inputs)?;
        let [a, _, c] = inputs;
        let [ab, bc] = self.products;
        let ty = model.rmext_type(self.target)?;
        let target = model.real_ext(ty)?;
        let mut directions = Basis::default();

        for (product, fixed, on_left) in [(bc, a.0, true), (ab, c.0, false)] {
            let Some(grade) = product.previous() else {
                continue;
            };
            if fixed.is_empty() || target.rank() == 0 {
                continue;
            }
            let ty = model.rmext_type(grade)?;
            let group = model.real_ext(ty)?;
            for index in 0..group.rank() {
                let generator = group.representative(model, &[index])?;
                let (left, right) = if on_left {
                    (fixed, generator.as_slice())
                } else {
                    (generator.as_slice(), fixed)
                };
                let image = real_motivic::concat_support(model, left, right)?;
                let coordinates = target.project(model, &image)?;
                directions.insert(model, coordinates)?;
            }
        }

        // In pivot order, reduce each direction by all earlier pivots to
        // obtain the unique reduced basis of the indeterminacy subspace.
        let mut basis = Basis::default();
        for direction in directions.vectors() {
            let direction = basis.reduce(model, direction.to_vec())?;
            basis.insert(model, direction)?;
        }
        let origin = target.project(model, &terms)?;
        let origin: Rc<[usize]> = basis.reduce(model, origin)?.into();
        model.allocate(basis.cells() + origin.len())?;
        Ok(Coset { basis, origin })
    }
}

pub(crate) fn adjacent_null(
    model: &mut Model,
    grade: Grade,
    left: Input<'_>,
    right: Input<'_>,
) -> Result<Value> {
    require_cycle(model, left)?;
    require_cycle(model, right)?;
    let terms = product_null(model, grade, left.0, right.0, "the adjacent product")?;
    real_motivic::persist(model, terms)
}

fn require_cycle(model: &mut Model, (terms, known): Input<'_>) -> Result<()> {
    if !known && !real_motivic::differential_support(model, terms)?.is_empty() {
        return Err(Error::plain(
            "real-motivic Massey products require closed inputs",
        ));
    }
    Ok(())
}

fn product_null(
    model: &mut Model,
    grade: Grade,
    left: &[Term],
    right: &[Term],
    label: &str,
) -> Result<Vec<Term>> {
    let product = real_motivic::concat_support(model, left, right)?;
    real_ext::real_preimage(model, grade, &product, true)?.ok_or_else(|| {
        Error::plain(format!(
            "real-motivic Massey product is undefined because {label} is not a boundary"
        ))
    })
}

pub(crate) struct Coset {
    basis: Basis,
    origin: Rc<[usize]>,
}

impl Coset {
    pub(crate) fn base(&self) -> Value {
        Value::Class(self.origin.clone())
    }

    pub(crate) fn dimension(&self) -> usize {
        self.basis.rank()
    }

    pub(crate) fn direction(&self, model: &mut Model, index: usize) -> Result<Value> {
        let direction = self.basis.vectors().nth(index).ok_or_else(|| {
            Error::plain(format!(
                "rmTripleDir index {index} is outside the indeterminacy basis"
            ))
        })?;
        model.spend(index + direction.len())?;
        model.allocate(direction.len())?;
        Ok(Value::Class(direction.into()))
    }

    pub(crate) fn contains(&self, model: &mut Model, candidate: &[usize]) -> Result<bool> {
        let reduced = self.basis.reduce(model, candidate.to_vec())?;
        Ok(reduced.as_slice() == self.origin.as_ref())
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        model.spend(
            self.basis.cells() + other.basis.cells() + self.origin.len() + other.origin.len(),
        )?;
        Ok(self
            .basis
            .vectors()
            .cmp(other.basis.vectors())
            .then_with(|| self.origin.cmp(&other.origin)))
    }
}
