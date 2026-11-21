use std::fmt;

use crate::bp::{self, Operation, Shape};
use crate::bp_classes;
use crate::bp_ext;
use crate::bp_polynomial::{Monomial, Polynomial};
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::odd_cobar::Prime;
use crate::odd_ext;
use crate::prime_linear::{Basis, Coordinate, QuotientBasis};

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) shape: Shape,
    pub(crate) t: usize,
}

impl fmt::Display for Grade {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Shape {
            prime,
            height,
            arrows,
            depth,
            thickness,
        } = self.shape;
        if let Some(j) = thickness {
            write!(
                out,
                "BPThickExt({prime}, {height}, {arrows}, {}, {j}, {})",
                depth + 1,
                self.t
            )
        } else {
            write!(
                out,
                "BPIExt({prime}, {height}, {arrows}, {depth}, {})",
                self.t
            )
        }
    }
}

impl Grade {
    pub(crate) fn prime(self) -> Prime {
        Prime::validated(self.shape.prime)
    }

    pub(crate) fn cup(self, other: Self) -> Result<Self> {
        let shape = self.shape.cup(other.shape)?;
        let t = self
            .t
            .checked_add(other.t)
            .ok_or_else(|| Error::plain("bpCup internal degree is out of range"))?;
        Ok(Self { shape, t })
    }

    pub(crate) fn map(self, operation: Operation) -> Result<Self> {
        let shape = operation.target(self.shape)?;
        let power = match operation {
            Operation::ThickDelta => self.shape.thickness.map(|j| (self.shape.depth, j)),
            Operation::Delta if self.shape.depth > 1 => Some((shape.depth, 1)),
            Operation::Delta | Operation::Reduce(_) | Operation::Thick(_) => None,
            _ => {
                return Err(Error::plain(format!(
                    "{} requires a BP cochain",
                    operation.name()
                )))
            }
        };
        let mut t = self.t;
        if let Some((index, exponent)) = power {
            let drop = u128::from(self.shape.prime)
                .checked_pow(index as u32)
                .and_then(|degree| (degree - 1).checked_mul(2 * exponent as u128))
                .and_then(|degree| usize::try_from(degree).ok());
            t = drop.and_then(|drop| t.checked_sub(drop)).ok_or_else(|| {
                Error::plain(format!(
                    "{} has negative output internal degree",
                    operation.name()
                ))
            })?;
        }
        Ok(Self { shape, t })
    }
}

pub(crate) struct Group {
    prime: Prime,
    monomials: Vec<Monomial>,
    quotient: QuotientBasis,
}

impl Group {
    pub(crate) fn build(model: &mut Model, grade: Grade) -> Result<Self> {
        let prime = grade.prime();
        let shape = grade.shape;
        let monomials = bp_ext::basis(model, shape, grade.t)?;
        let quotient = if monomials.is_empty() {
            QuotientBasis::default()
        } else {
            let mut boundaries = Basis::default();
            if shape.arrows != 0 {
                let previous = bp_ext::previous_shape(shape)?;
                for monomial in bp_ext::basis(model, previous, grade.t)? {
                    let column = differential(model, previous, &monomial, &monomials, prime)?;
                    boundaries.insert(model, prime, column)?;
                }
            }
            let next = Shape {
                arrows: shape.arrows + 1,
                ..shape
            };
            let outgoing = bp_ext::basis(model, next, grade.t)?;
            QuotientBasis::build(model, prime, boundaries, monomials.len(), |model, index| {
                differential(model, shape, &monomials[index], &outgoing, prime)
            })?
        };
        model.allocate(
            monomials
                .iter()
                .map(|monomial| 2 + 2 * monomial.len())
                .sum(),
        )?;
        Ok(Self {
            prime,
            monomials,
            quotient,
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.quotient.rank()
    }

    pub(crate) fn project(&self, model: &mut Model, cycle: &Polynomial) -> Result<Vec<Coordinate>> {
        let coordinates = field_coordinates(model, cycle, &self.monomials, self.prime)?;
        self.quotient
            .project(model, self.prime, coordinates)?
            .ok_or_else(|| {
                Error::plain("mod-p BP Ext projection requires a closed normalized cochain")
            })
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        class: &[Coordinate],
    ) -> Result<Polynomial> {
        let coordinates = self.quotient.representative(model, self.prime, class)?;
        Polynomial::from_terms(
            model,
            coordinates
                .into_iter()
                .map(|(i, c)| (&self.monomials[i], i128::from(c))),
        )
    }
}

fn field_coordinates(
    model: &mut Model,
    cycle: &Polynomial,
    basis: &[Monomial],
    prime: Prime,
) -> Result<Vec<Coordinate>> {
    let mut result = Vec::new();
    for (monomial, coefficient) in cycle.terms() {
        let index = bp_ext::index_of(model, basis, monomial)?;
        let coefficient = coefficient.rem_euclid(i128::from(prime.get())) as u32;
        if coefficient != 0 {
            result.push((index, coefficient));
        }
    }
    Ok(result)
}

fn differential(
    model: &mut Model,
    shape: Shape,
    monomial: &Monomial,
    target: &[Monomial],
    prime: Prime,
) -> Result<Vec<Coordinate>> {
    if target.is_empty() {
        return Ok(Vec::new());
    }
    let input = Polynomial::from_terms(model, [(monomial, 1)])?;
    let output = bp::differential(model, shape, &input)?;
    field_coordinates(model, &output, target, prime)
}

pub(crate) fn class(model: &mut Model, ty: TypeId, cochain: &Polynomial) -> Result<Value> {
    if cochain.is_zero() {
        return Ok(odd_ext::zero());
    }
    let group = model.bp_mod_ext(ty)?;
    let coordinates = group.project(model, cochain)?;
    odd_ext::persist(model, coordinates)
}

pub(crate) fn cup(
    model: &mut Model,
    left: (TypeId, &[Coordinate]),
    right: (TypeId, &[Coordinate]),
    target: TypeId,
) -> Result<Value> {
    if left.1.is_empty() || right.1.is_empty() {
        return Ok(odd_ext::zero());
    }
    let (&Kind::BPModExt(a), &Kind::BPModExt(b)) = (model.kind(left.0), model.kind(right.0)) else {
        unreachable!()
    };
    let f = model.bp_mod_ext(left.0)?.representative(model, left.1)?;
    let g = model.bp_mod_ext(right.0)?.representative(model, right.1)?;
    let product = bp::cup(model, a.shape, b.shape, &f, &g)?;
    class(model, target, &product)
}

pub(crate) fn apply(
    model: &mut Model,
    operation: Operation,
    source: TypeId,
    target: TypeId,
    coordinates: &[Coordinate],
) -> Result<Value> {
    let integral = matches!(model.kind(target), Kind::BPExt(..));
    if coordinates.is_empty() {
        return Ok(if integral {
            bp_classes::zero()
        } else {
            odd_ext::zero()
        });
    }
    let Kind::BPModExt(grade) = *model.kind(source) else {
        unreachable!()
    };
    let cochain = model
        .bp_mod_ext(source)?
        .representative(model, coordinates)?;
    let image = operation.apply(model, grade.shape, &cochain)?;
    if integral {
        bp_classes::class(model, target, &image)
    } else {
        class(model, target, &image)
    }
}
