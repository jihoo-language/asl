use std::collections::HashMap;
use std::rc::Rc;

use crate::ec_poly::Scalar;
use crate::error::{Error, Result};
use crate::hw_algebra::{self, Matrix};
use crate::hw_closed::Curve;
use crate::hyperelliptic::Residues;
use crate::hyperelliptic_jacobian::{miller, Divisor, Pairing};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_unram::{mul_mod, Ring};
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Command {
    Define(Label, Label, usize, Box<ParsedExpr>, Vec<ParsedExpr>),
    Frobenius(Label, Label, usize),
}

pub(crate) struct Basis {
    pub(crate) curve_name: Rc<str>,
    pub(crate) prime: Prime,
    curve: Rc<Curve>,
    generators: Box<[Rc<Divisor>]>,
    logarithm: Logarithm,
}

pub(crate) struct Frobenius {
    pub(crate) matrix: Matrix<u64>,
    pub(crate) multiplier: u64,
    pub(crate) characteristic: Vec<u64>,
}

impl Basis {
    pub(crate) fn new(
        model: &mut Model,
        curve_name: Rc<str>,
        curve: Rc<Curve>,
        prime: Prime,
        root: &[u64],
        generators: Vec<Rc<Divisor>>,
    ) -> Result<Self> {
        let order = u64::from(prime.get());
        let genus = curve.grade.genus;
        if generators.len() != 2 * genus {
            return Err(Error::plain(format!(
                "Jacobian torsion basis requires exactly {} generators",
                2 * genus
            )));
        }
        if order == curve.ring.prime() {
            return Err(Error::plain(
                "Jacobian torsion level must differ from the residue characteristic",
            ));
        }
        let one = curve.ring.scalar(1);
        if root == one || curve.ring.power(model, root, order)? != one {
            return Err(Error::plain(format!(
                "Jacobian torsion basis requires a primitive root of unity of order {order}"
            )));
        }
        for (index, divisor) in generators.iter().enumerate() {
            if !is_torsion(model, &curve, order, divisor)? {
                return Err(Error::plain(format!(
                    "Jacobian torsion basis generator {index} must be killed by {order}"
                )));
            }
        }
        // The Weil pairing is alternating, so the entries above the diagonal
        // determine the full Gram matrix.
        for (i, left) in generators.iter().enumerate() {
            for (j, right) in generators.iter().enumerate().skip(i + 1) {
                let expected = if j == i + genus { root } else { &one };
                let actual = pairing(model, &curve, order, left, right)?;
                if actual != expected {
                    return Err(Error::plain(format!(
                        "Jacobian torsion basis has incorrect symplectic pairing at ({i}, {j})"
                    )));
                }
            }
        }
        model.allocate(8 + curve_name.len() + generators.len())?;
        let logarithm = Logarithm::new(model, &curve.ring, order, root)?;
        Ok(Self {
            curve_name,
            prime,
            curve,
            generators: generators.into_boxed_slice(),
            logarithm,
        })
    }

    pub(crate) fn len(&self) -> usize {
        self.generators.len()
    }

    pub(crate) fn coordinate(&self, model: &mut Model, t: &Divisor, index: usize) -> Result<u32> {
        let order = u64::from(self.prime.get());
        if !is_torsion(model, &self.curve, order, t)? {
            return Err(Error::plain(format!(
                "rvJacTorsionCoord requires a class killed by {order}"
            )));
        }
        self.coordinate_of_torsion(model, t, index)
    }

    fn coordinate_of_torsion(&self, model: &mut Model, t: &Divisor, index: usize) -> Result<u32> {
        let genus = self.curve.grade.genus;
        let (left, right) = if index < genus {
            (t, self.generators[genus + index].as_ref())
        } else {
            (self.generators[index - genus].as_ref(), t)
        };
        let order = u64::from(self.prime.get());
        let value = pairing(model, &self.curve, order, left, right)?;
        self.logarithm
            .evaluate(model, &self.curve.ring, order, value)
    }

    pub(crate) fn frobenius(
        &self,
        model: &mut Model,
        base: &Curve,
        power: usize,
    ) -> Result<Frobenius> {
        let curve = &self.curve;
        let genus = curve.grade.genus;
        let width = self.len();
        let order = u64::from(self.prime.get());
        let ring = Residues::new(self.prime, 1)?;
        let cells = ring.matrix_size(width)?;
        model.spend(cells)?;
        let mut entries = vec![0; cells];
        if width == 0 {
            // Even an empty basis requires a valid base extension.
            Divisor::zero(curve).frobenius(model, base, curve, power)?;
        }
        for (column, generator) in self.generators.iter().enumerate() {
            let image = generator.frobenius(model, base, curve, power)?;
            let mut reconstruction = Divisor::zero(curve);
            for (row, basis_element) in self.generators.iter().enumerate() {
                let coefficient = self.coordinate_of_torsion(model, &image, row)?;
                entries[row * width + column] = u64::from(coefficient);
                if coefficient != 0 {
                    let term = basis_element.scale(model, curve, i128::from(coefficient))?;
                    reconstruction = reconstruction.add(model, curve, &term)?;
                }
            }
            if !reconstruction.compare(model, &image)?.is_eq() {
                return Err(Error::plain(
                    "Jacobian torsion Frobenius column failed to reconstruct its divisor class",
                ));
            }
        }
        let matrix = Matrix::from_entries(width, entries)?;
        let q = power_mod(model, base.ring.prime(), base.ring.degree(), order)?;
        let multiplier = power_mod(model, q, power, order)?;
        verify_similitude(model, self.prime, genus, multiplier as u32, &matrix)?;
        let characteristic = hw_algebra::characteristic(model, &ring, &matrix)?;
        Ok(Frobenius {
            matrix,
            multiplier,
            characteristic,
        })
    }
}

fn pairing(
    model: &mut Model,
    curve: &Curve,
    order: u64,
    left: &Divisor,
    right: &Divisor,
) -> Result<Scalar> {
    if order == 2 {
        Pairing::Weil.evaluate(model, curve, left, right)
    } else {
        miller::weil(model, curve, order, left, right)
    }
}

fn is_torsion(model: &mut Model, curve: &Curve, order: u64, divisor: &Divisor) -> Result<bool> {
    if order == 2 {
        // Cantor negation sends the canonical v to −v in odd characteristic.
        return Ok(divisor.v.is_empty());
    }
    Ok(divisor.scale(model, curve, i128::from(order))?.u.len() == 1)
}

struct Logarithm {
    babies: HashMap<Scalar, u32>,
    step: u64,
    inverse_step: Scalar,
}

impl Logarithm {
    fn new(model: &mut Model, ring: &Ring, order: u64, root: &[u64]) -> Result<Self> {
        let step = (order - 1).isqrt() + 1;
        let width = ring.degree() + 7;
        model.allocate((step as usize + 1) * width)?;
        let mut babies = HashMap::with_capacity(step as usize);
        let mut value = ring.scalar(1);
        for j in 0..step {
            let next = ring.multiply(model, &value, root)?;
            model.spend(width)?;
            babies.insert(value, j as u32);
            value = next;
        }
        let inverse_step = ring.power(model, root, order - step)?;
        Ok(Self {
            babies,
            step,
            inverse_step,
        })
    }

    fn evaluate(
        &self,
        model: &mut Model,
        ring: &Ring,
        order: u64,
        mut value: Scalar,
    ) -> Result<u32> {
        for block in 0..order.div_ceil(self.step) {
            model.spend(ring.degree() + 1)?;
            if let Some(&baby) = self.babies.get(&value) {
                return Ok(((block * self.step + u64::from(baby)) % order) as u32);
            }
            value = ring.multiply(model, &value, &self.inverse_step)?;
        }
        Err(Error::plain(
            "Jacobian Weil pairing lies outside the certified roots of unity",
        ))
    }
}

fn power_mod(model: &mut Model, mut value: u64, mut power: usize, modulus: u64) -> Result<u64> {
    let mut result = 1;
    value %= modulus;
    while power != 0 {
        model.spend(1)?;
        if power & 1 == 1 {
            result = mul_mod(result, value, modulus);
        }
        value = mul_mod(value, value, modulus);
        power >>= 1;
    }
    Ok(result)
}

fn verify_similitude(
    model: &mut Model,
    prime: Prime,
    genus: usize,
    multiplier: u32,
    matrix: &Matrix<u64>,
) -> Result<()> {
    let width = matrix.size();
    for i in 0..width {
        for j in i + 1..width {
            model.spend(1 + genus)?;
            let mut pair = 0;
            for k in 0..genus {
                let forward = prime.mul(
                    *matrix.entry(k, i) as u32,
                    *matrix.entry(k + genus, j) as u32,
                );
                let reverse = prime.mul(
                    *matrix.entry(k + genus, i) as u32,
                    *matrix.entry(k, j) as u32,
                );
                pair = prime.add(pair, prime.add(forward, prime.negate(reverse)));
            }
            let expected = if j == i + genus { multiplier } else { 0 };
            if pair != expected {
                return Err(Error::plain(
                    "Jacobian torsion Frobenius failed its symplectic similitude check",
                ));
            }
        }
    }
    Ok(())
}
