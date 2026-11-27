use crate::error::{Error, Result};
use crate::model::Model;

use super::local::{Denominator, LocalAlgebra};
use super::{Algebra, Binary, Unary};

pub(crate) const PROJECTORS: usize = 7;
const ANCHORS: [usize; 7] = [0, 1, 1, 0, 1, 5, 6];

pub(crate) struct Square {
    projectors: [LocalAlgebra; 7],
    forward: [LocalAlgebra; 7],
    backward: [LocalAlgebra; 7],
}

type Diagram = LocalAlgebra;

impl Square {
    pub(crate) fn new(model: &mut Model) -> Result<Self> {
        let one = Diagram::identity(1);
        let identity = Diagram::identity(2);
        let j = basis(model, 1, &[0, 1])?;
        let u = divide(model, &j, [0, 1, 0, 0])?;
        let w = subtract(model, &one, &u)?;
        let a = u.binary(model, Binary::Tensor, &u)?;
        let b = u.binary(model, Binary::Tensor, &w)?;
        let c = w.binary(model, Binary::Tensor, &u)?;
        let p = w.binary(model, Binary::Tensor, &w)?;
        let swap = basis(model, 2, &[0, 1, 1, 0])?;
        let delta_delta = basis(model, 2, &[0, 0, 0, 0])?;
        let delta_j_delta = basis(model, 2, &[0, 0, 1, 1])?;
        let trivial = sandwich(model, &p, &delta_j_delta, &p)?;
        let g = divide(model, &trivial, [0, 0, 1, 0])?;
        let delta_u_delta = divide(model, &delta_j_delta, [0, 1, 0, 0])?;
        let delta_w_delta = subtract(model, &delta_delta, &delta_u_delta)?;
        let standard = sandwich(model, &p, &delta_w_delta, &p)?;
        let standard = standard.unary(model, Unary::Scale(1, 1))?;
        let h = divide(model, &standard, [0, 0, 0, 1])?;
        let plus = identity.add(model, &swap)?;
        let minus = subtract(model, &identity, &swap)?;
        let sym = p.binary(model, Binary::Compose, &plus)?;
        let alt = p.binary(model, Binary::Compose, &minus)?;
        let sym = divide(model, &sym, [1, 0, 0, 0])?;
        let l = divide(model, &alt, [1, 0, 0, 0])?;
        let remaining = subtract(model, &sym, &g)?;
        let k = subtract(model, &remaining, &h)?;
        let projectors = [a, b, c, g, h, k, l];
        certify_projectors(model, &projectors, &identity)?;

        // The first bridge is |P δu><u⊗u|, of squared length t²(t−1).
        let trivial_bridge = basis(model, 2, &[0, 0, 1, 2])?;
        let trivial_bridge = sandwich(model, &p, &trivial_bridge, &projectors[0])?;
        let trivial_adjoint = trivial_bridge.unary(model, Unary::Reflect)?;
        let trivial_inverse = divide(model, &trivial_adjoint, [0, 2, 1, 0])?;
        // This bridge sends u⊗w to t P δw, of squared length t−2
        // relative to the norm of u⊗w.
        let standard_bridge = basis(model, 2, &[0, 0, 1, 0])?;
        let standard_bridge = sandwich(model, &p, &standard_bridge, &projectors[1])?;
        let standard_adjoint = standard_bridge.unary(model, Unary::Reflect)?;
        let standard_inverse = divide(model, &standard_adjoint, [0, 0, 0, 1])?;
        let swap_bridge = swap.binary(model, Binary::Compose, &projectors[1])?;
        let swap_inverse = swap_bridge.unary(model, Unary::Reflect)?;
        let forward = [
            projectors[0].clone(),
            projectors[1].clone(),
            swap_bridge,
            trivial_bridge,
            standard_bridge,
            projectors[5].clone(),
            projectors[6].clone(),
        ];
        let backward = [
            projectors[0].clone(),
            projectors[1].clone(),
            swap_inverse,
            trivial_inverse,
            standard_inverse,
            projectors[5].clone(),
            projectors[6].clone(),
        ];
        for (i, &anchor) in ANCHORS.iter().enumerate() {
            let source = backward[i].binary(model, Binary::Compose, &forward[i])?;
            let target = forward[i].binary(model, Binary::Compose, &backward[i])?;
            let f = forward[i].binary(model, Binary::Compose, &projectors[anchor])?;
            let b = projectors[anchor].binary(model, Binary::Compose, &backward[i])?;
            if source != projectors[anchor]
                || target != projectors[i]
                || f != forward[i]
                || b != backward[i]
            {
                return Err(Error::plain(
                    "Deligne square matrix unit failed its corner identities",
                ));
            }
        }
        let units = ANCHORS
            .iter()
            .flat_map(|a| ANCHORS.iter().map(move |b| usize::from(a == b)))
            .sum::<usize>();
        let (dimension, _) = super::dimension::dimensions(model, 4, 0)?;
        if dimension.magnitude.small() != Some(units as u128) {
            return Err(Error::plain(
                "Deligne square matrix units do not exhaust the diagram dimension",
            ));
        }
        for element in projectors.iter().chain(&forward).chain(&backward) {
            element.charge(model)?;
        }
        Ok(Self {
            projectors,
            forward,
            backward,
        })
    }

    pub(crate) fn unit(&self, model: &mut Model, i: usize, j: usize) -> Result<Diagram> {
        if i == j {
            return Ok(self.projectors[i].clone());
        }
        if ANCHORS[i] != ANCHORS[j] {
            return Ok(Diagram::zero(2));
        }
        self.forward[i].binary(model, Binary::Compose, &self.backward[j])
    }
}

fn certify_projectors(
    model: &mut Model,
    projectors: &[Diagram; 7],
    identity: &Diagram,
) -> Result<()> {
    let mut sum = Diagram::zero(2);
    for (i, left) in projectors.iter().enumerate() {
        if left.is_zero() {
            return Err(Error::plain("Deligne square projector is zero"));
        }
        if left.unary(model, Unary::Reflect)? != *left {
            return Err(Error::plain("Deligne square projector is not self-adjoint"));
        }
        sum = sum.add(model, left)?;
        for (j, right) in projectors.iter().enumerate() {
            let product = left.binary(model, Binary::Compose, right)?;
            if (i == j && product != *left) || (i != j && !product.is_zero()) {
                return Err(Error::plain(
                    "Deligne square projectors failed orthogonal idempotence",
                ));
            }
        }
    }
    if sum != *identity {
        return Err(Error::plain("Deligne square projectors are incomplete"));
    }
    Ok(())
}

fn basis(model: &mut Model, rank: usize, labels: &[u128]) -> Result<Diagram> {
    Diagram::embed(Algebra::basis(model, rank, labels.to_vec(), false)?)
}

fn subtract(model: &mut Model, a: &Diagram, b: &Diagram) -> Result<Diagram> {
    let negative = b.unary(model, Unary::Scale(-1, 0))?;
    a.add(model, &negative)
}

fn divide(model: &mut Model, a: &Diagram, powers: [u128; 4]) -> Result<Diagram> {
    a.divide(model, Denominator::new(powers)?)
}

fn sandwich(
    model: &mut Model,
    left: &Diagram,
    middle: &Diagram,
    right: &Diagram,
) -> Result<Diagram> {
    let mr = middle.binary(model, Binary::Compose, right)?;
    left.binary(model, Binary::Compose, &mr)
}
