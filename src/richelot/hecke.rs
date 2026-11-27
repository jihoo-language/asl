use std::rc::Rc;

use crate::ec_poly::Field;
use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::jacobian_level::{binary_mask_image, Basis};
use crate::model::{Model, Value};
use crate::syntax::Label;

use super::{register, transform, validate_splitting, Splitting};

pub(crate) fn report(model: &mut Model, basis: &Label, base: &Label) -> Result<String> {
    let neighbors = Neighbors::new(model, basis, base)?;
    let field = Field {
        ring: &neighbors.curve.ring,
    };
    let mut stable = Vec::with_capacity(neighbors.stable.len());
    for plane in &neighbors.stable {
        let factors = neighbors.splitting(model, plane)?;
        validate_splitting(model, &neighbors.curve, &factors)?;
        let determinant_zero = transform(model, field, &factors)?.is_none();
        stable.push(format!(
            "{{masks={:?}, rational={}, determinant_zero={determinant_zero}}}",
            plane.masks, plane.rational
        ));
    }
    Ok(format!(
        "rvJacRichelotKernels[{}, {}] = orbits={:?}, stable=[{}]",
        basis.text,
        base.text,
        neighbors.orbits,
        stable.join(", ")
    ))
}

pub(crate) fn kernel(
    model: &mut Model,
    basis: &Label,
    base: &Label,
    plane: usize,
    index: usize,
) -> Result<(Rc<str>, Value)> {
    let neighbors = Neighbors::new(model, basis, base)?;
    let plane = neighbors.plane(basis, plane)?;
    let mask = plane.masks.get(index).ok_or_else(|| {
        Error::new(
            basis.span,
            "Richelot kernel generator index is out of range",
        )
    })?;
    let divisor = neighbors.basis.binary_combination(model, *mask)?;
    Ok((
        neighbors.basis.curve_name.clone(),
        divisor.persist(model, &neighbors.curve)?,
    ))
}

pub(crate) fn define(
    model: &mut Model,
    target: Label,
    basis: &Label,
    base: &Label,
    index: usize,
) -> Result<()> {
    let neighbors = Neighbors::new(model, basis, base)?;
    let plane = neighbors.plane(basis, index)?;
    let factors = neighbors.splitting(model, plane)?;
    let source = Label {
        text: neighbors.basis.curve_name.to_string(),
        span: basis.span,
    };
    register(model, target, source, &neighbors.curve, factors)
}

type Masks = [u8; 3];

struct Plane {
    masks: Masks,
    rational: usize,
}

struct Neighbors {
    basis: Rc<Basis>,
    curve: Rc<Curve>,
    orbits: Vec<usize>,
    stable: Vec<Plane>,
}

impl Neighbors {
    fn new(model: &mut Model, basis_name: &Label, base: &Label) -> Result<Self> {
        let basis = model.jacobian_level(basis_name)?;
        if basis.prime.get() != 2 || basis.len() != 4 {
            return Err(Error::new(
                basis_name.span,
                "Richelot kernels require a certified genus-two level-2 basis",
            ));
        }
        let curve = model.jacobian_curve(&basis.curve_name)?;
        let base_curve = model
            .jacobian_curve(&base.text)
            .map_err(|error| error.at(base.span))?;
        let frobenius = basis
            .frobenius(model, &base_curve, 1)
            .map_err(|error| error.at(base.span))?;
        model.spend(16 * 16)?;
        let masks = isotropic_planes();
        let images: [u8; 16] =
            std::array::from_fn(|x| binary_mask_image(&frobenius.matrix, x as u8));
        let permutation = masks
            .iter()
            .map(|triple| {
                let mut image = triple.map(|x| images[usize::from(x)]);
                image.sort_unstable();
                masks
                    .binary_search(&image)
                    .expect("certified Frobenius preserves Lagrangians")
            })
            .collect::<Vec<_>>();
        let orbits = orbit_lengths(&permutation);
        let stable = masks
            .into_iter()
            .enumerate()
            .filter(|(index, _)| *index == permutation[*index])
            .map(|(_, masks)| Plane {
                masks,
                rational: 1 + masks
                    .iter()
                    .filter(|&&x| x == images[usize::from(x)])
                    .count(),
            })
            .collect();
        Ok(Self {
            basis,
            curve,
            orbits,
            stable,
        })
    }

    fn plane(&self, source: &Label, index: usize) -> Result<&Plane> {
        self.stable
            .get(index)
            .ok_or_else(|| Error::new(source.span, "stable Richelot plane index is out of range"))
    }

    fn splitting(&self, model: &mut Model, plane: &Plane) -> Result<Splitting> {
        let [a, b, c] = plane.masks;
        Ok([
            self.basis.binary_combination(model, a)?.u,
            self.basis.binary_combination(model, b)?.u,
            self.basis.binary_combination(model, c)?.u,
        ])
    }
}

fn isotropic_planes() -> Vec<Masks> {
    let mut planes = Vec::with_capacity(15);
    for a in 1..16u8 {
        for b in a + 1..16 {
            let c = a ^ b;
            let pairing = ((a & (b >> 2)) ^ (b & (a >> 2))).count_ones() % 2;
            // a < b < c keeps one representative per plane, in lexicographic order.
            if b < c && pairing == 0 {
                planes.push([a, b, c]);
            }
        }
    }
    planes
}

fn orbit_lengths(permutation: &[usize]) -> Vec<usize> {
    let mut lengths = Vec::new();
    let mut seen = vec![false; permutation.len()];
    for index in 0..permutation.len() {
        let mut length = 0;
        let mut current = index;
        while !seen[current] {
            seen[current] = true;
            length += 1;
            current = permutation[current];
        }
        if length != 0 {
            lengths.push(length);
        }
    }
    lengths.sort_unstable();
    lengths
}
