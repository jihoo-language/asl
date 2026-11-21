use std::cell::OnceCell;
use std::cmp::Ordering;
use std::rc::Rc;

use crate::bp::{self, Shape};
use crate::bp_ext::{self, Domain};
use crate::bp_polynomial::{Monomial, Polynomial};
use crate::bp_rel;
use crate::bp_smith::{PrimaryFrame, Projector};
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct Grade {
    pub(crate) shape: Shape,
    pub(crate) t: usize,
}

impl Grade {
    pub(crate) fn cup(self, other: Self) -> Result<Self> {
        let shape = self.shape.cup(other.shape)?;
        let t = self
            .t
            .checked_add(other.t)
            .ok_or_else(|| Error::plain("bpCup internal degree is out of range"))?;
        Ok(Self { shape, t })
    }
}

pub(crate) type Support = Rc<[(usize, i128)]>;

pub(crate) struct Cohomology {
    domain: Domain,
    t: usize,
    pub(crate) invariants: bp_ext::Group,
    frame: OnceCell<Frame>,
}

struct Frame {
    basis: Vec<Monomial>,
    generators: Vec<Polynomial>,
    projectors: Vec<Projector>,
}

impl Cohomology {
    pub(crate) fn build(model: &mut Model, grade: Grade) -> Result<Self> {
        Self::for_domain(model, Domain::absolute(grade.shape), grade.t)
    }

    pub(crate) fn build_relative(model: &mut Model, grade: bp_rel::Grade) -> Result<Self> {
        Self::for_domain(model, grade.shape.domain(), grade.t)
    }

    fn for_domain(model: &mut Model, domain: Domain, t: usize) -> Result<Self> {
        let invariants = domain.group(model, t)?;
        model.allocate(2 * invariants.torsion.len())?;
        Ok(Self {
            domain,
            t,
            invariants,
            frame: OnceCell::new(),
        })
    }

    fn is_free(&self) -> bool {
        self.invariants.free != 0
    }

    pub(crate) fn rank(&self) -> usize {
        self.invariants.free + self.invariants.torsion.len()
    }

    fn frame(&self, model: &mut Model) -> Result<&Frame> {
        if let Some(frame) = self.frame.get() {
            return Ok(frame);
        }
        let frame = self.build_frame(model)?;
        Ok(self.frame.get_or_init(|| frame))
    }

    fn build_frame(&self, model: &mut Model) -> Result<Frame> {
        let previous = self.domain.previous()?;
        let (source, rank) = previous.incoming_basis(model, self.t)?;
        let basis = self.domain.basis(model, self.t)?;
        let incoming = previous.matrix(model, &source, &basis)?;
        let orders = &self.invariants.torsion;
        let prime = self.domain.shape.prime;
        let smith = PrimaryFrame::build(model, &incoming, prime, rank, orders)?;
        let mut cells = basis
            .iter()
            .map(|monomial| 1 + 2 * monomial.len())
            .sum::<usize>();
        let mut projectors = Vec::with_capacity(orders.len());
        for i in 0..orders.len() {
            let projector = smith.projector(model, i)?;
            cells += 1 + projector.cells();
            check_frame(cells)?;
            projectors.push(projector);
        }
        let mut generators = Vec::with_capacity(orders.len());
        for (i, &modulus) in orders.iter().enumerate() {
            let preimage = smith.lift(model, i)?;
            let mut cycle = incoming.apply(model, &preimage)?;
            // The exact incoming image is divisible by the local Smith pivot.
            // Dividing gives an integral cycle with the corresponding unit coordinate.
            model.spend(cycle.len())?;
            for coefficient in &mut cycle {
                if *coefficient % modulus != 0 {
                    return Err(Error::plain("BP Ext generator is not exactly divisible"));
                }
                *coefficient /= modulus;
            }
            for (j, projector) in projectors.iter().enumerate() {
                if projector.evaluate(model, &cycle)? != i128::from(i == j) {
                    return Err(Error::plain(
                        "BP Ext generator does not have the expected coordinates",
                    ));
                }
            }
            let representative = Polynomial::from_terms(model, basis.iter().zip(cycle))?;
            cells += 1 + representative.cells();
            check_frame(cells)?;
            generators.push(representative);
        }
        model.allocate(cells)?;
        Ok(Frame {
            basis,
            generators,
            projectors,
        })
    }

    fn project(&self, model: &mut Model, cycle: &Polynomial) -> Result<Support> {
        if cycle.is_zero() {
            return Ok([].into());
        }
        if self.is_free() {
            let coefficient = cycle.terms().next().unwrap().1;
            return coordinates(model, vec![(0, coefficient)]);
        }
        if self.invariants.torsion.is_empty() {
            return Ok([].into());
        }
        let frame = self.frame(model)?;
        let vector = bp_ext::coordinates(model, cycle, &frame.basis)?;
        let mut support = Vec::new();
        for (index, projector) in frame.projectors.iter().enumerate() {
            let coefficient = projector.evaluate(model, &vector)?;
            if coefficient != 0 {
                support.push((index, coefficient));
            }
        }
        coordinates(model, support)
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        support: &[(usize, i128)],
    ) -> Result<Polynomial> {
        if support.is_empty() {
            return Ok(Polynomial::default());
        }
        if self.is_free() {
            return Ok(Polynomial::integer(support[0].1));
        }
        let frame = self.frame(model)?;
        let mut result = Polynomial::default();
        for &(coordinate, coefficient) in support {
            let mut term = frame.generators[coordinate].clone();
            term.scale(model, coefficient)?;
            result.add_assign(model, &term, false)?;
        }
        Ok(result)
    }

    pub(crate) fn elements(&self, model: &mut Model, limit: usize) -> Result<Vec<Value>> {
        if self.is_free() {
            return Err(Error::plain(
                "cannot enumerate the infinite free BP Ext group",
            ));
        }
        let orders = &self.invariants.torsion;
        let mut count = 1usize;
        for &order in orders {
            let order = usize::try_from(order).map_err(|_| too_many_elements(limit))?;
            count = count
                .checked_mul(order)
                .filter(|&count| count <= limit)
                .ok_or_else(|| too_many_elements(limit))?;
        }
        model.spend(count * (1 + orders.len()))?;
        model.allocate(count)?;
        let mut result = Vec::with_capacity(count);
        for mut digits in 0..count {
            let mut support = Vec::new();
            for (index, &order) in orders.iter().enumerate() {
                let coefficient = digits % order as usize;
                digits /= order as usize;
                if coefficient != 0 {
                    support.push((index, coefficient as i128));
                }
            }
            result.push(Value::BPExt(coordinates(model, support)?));
        }
        Ok(result)
    }
}

fn too_many_elements(limit: usize) -> Error {
    Error::plain(format!(
        "type exceeds the enumeration limit of {limit} elements"
    ))
}

fn check_frame(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "BP Ext coordinate frame exceeds the storage limit",
        ));
    }
    Ok(())
}

fn coordinates(model: &mut Model, support: Vec<(usize, i128)>) -> Result<Support> {
    model.allocate(3 * support.len())?;
    Ok(support.into())
}

pub(crate) fn zero() -> Value {
    Value::BPExt([].into())
}

pub(crate) fn generator(model: &mut Model, index: usize) -> Result<Value> {
    Ok(Value::BPExt(coordinates(model, vec![(index, 1)])?))
}

pub(crate) fn class(model: &mut Model, ty: TypeId, cycle: &Polynomial) -> Result<Value> {
    if cycle.is_zero() {
        return Ok(zero());
    }
    let (domain, t) = match *model.kind(ty) {
        Kind::BPExt(grade) => (Domain::absolute(grade.shape), grade.t),
        Kind::BPRelExt(grade) => (grade.shape.domain(), grade.t),
        _ => unreachable!(),
    };
    domain.validate_cycle(model, t, cycle)?;
    let group = model.bp_ext(ty)?;
    Ok(Value::BPExt(group.project(model, cycle)?))
}

pub(crate) fn representative(model: &mut Model, ty: TypeId, support: &Support) -> Result<Value> {
    let relative = matches!(model.kind(ty), Kind::BPRelExt(..));
    if support.is_empty() {
        return Ok(if relative { bp_rel::zero() } else { bp::zero() });
    }
    let group = model.bp_ext(ty)?;
    let polynomial = group.representative(model, support)?;
    if relative {
        bp_rel::persist(model, polynomial)
    } else {
        bp::persist(model, polynomial)
    }
}

pub(crate) fn add(model: &mut Model, ty: TypeId, left: &Support, right: &Support) -> Result<Value> {
    if left.is_empty() || right.is_empty() {
        return Ok(Value::BPExt(if left.is_empty() {
            right.clone()
        } else {
            left.clone()
        }));
    }
    let group = model.bp_ext(ty)?;
    let mut left = left.as_ref();
    let mut right = right.as_ref();
    model.spend(left.len() + right.len())?;
    let mut support = Vec::new();
    while let (Some((&(i, a), rest_left)), Some((&(j, b), rest_right))) =
        (left.split_first(), right.split_first())
    {
        match i.cmp(&j) {
            Ordering::Less => {
                support.push((i, a));
                left = rest_left;
            }
            Ordering::Greater => {
                support.push((j, b));
                right = rest_right;
            }
            Ordering::Equal => {
                let sum = a.checked_add(b).ok_or_else(coefficient_overflow)?;
                let sum = group
                    .invariants
                    .torsion
                    .get(i)
                    .map_or(sum, |&order| sum.rem_euclid(order));
                if sum != 0 {
                    support.push((i, sum));
                }
                left = rest_left;
                right = rest_right;
            }
        }
    }
    support.extend_from_slice(left);
    support.extend_from_slice(right);
    Ok(Value::BPExt(coordinates(model, support)?))
}

pub(crate) fn scale(
    model: &mut Model,
    ty: TypeId,
    support: &Support,
    scalar: i128,
) -> Result<Value> {
    if scalar == 0 || support.is_empty() {
        return Ok(zero());
    }
    if scalar == 1 {
        return Ok(Value::BPExt(support.clone()));
    }
    let group = model.bp_ext(ty)?;
    model.spend(support.len())?;
    let mut result = Vec::new();
    for &(index, coefficient) in support.iter() {
        let product = if let Some(&order) = group.invariants.torsion.get(index) {
            let scalar = scalar.rem_euclid(order) as u128;
            (scalar * coefficient as u128 % order as u128) as i128
        } else {
            coefficient
                .checked_mul(scalar)
                .ok_or_else(coefficient_overflow)?
        };
        if product != 0 {
            result.push((index, product));
        }
    }
    Ok(Value::BPExt(coordinates(model, result)?))
}

pub(crate) fn cup(
    model: &mut Model,
    (left_ty, left): (TypeId, &Support),
    (right_ty, right): (TypeId, &Support),
    target: TypeId,
) -> Result<Value> {
    if left.is_empty() || right.is_empty() {
        return Ok(zero());
    }
    let left_group = model.bp_ext(left_ty)?;
    let right_group = model.bp_ext(right_ty)?;
    // The empty arrow word with integer coefficient acts by its scalar.
    if left_group.is_free() {
        return scale(model, target, right, left[0].1);
    }
    if right_group.is_free() {
        return scale(model, target, left, right[0].1);
    }
    let f = left_group.representative(model, left)?;
    let g = right_group.representative(model, right)?;
    let product = bp::cup(
        model,
        left_group.domain.shape,
        right_group.domain.shape,
        &f,
        &g,
    )?;
    class(model, target, &product)
}

fn coefficient_overflow() -> Error {
    Error::plain("BP Ext coefficient exceeds the signed i128 range")
}
