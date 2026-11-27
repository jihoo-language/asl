use std::cell::OnceCell;
use std::collections::HashMap;
use std::rc::Rc;

use crate::ec_poly::Scalar;
use crate::error::{Error, Result};
use crate::hw_algebra::{self, Matrix};
use crate::hw_closed::Curve;
use crate::hyperelliptic::Residues;
use crate::hyperelliptic_jacobian::{miller, Divisor, Pairing};
use crate::hyperelliptic_kummer_surface::theta;
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_unram::{mul_mod, Ring};
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod brauer;
pub(crate) mod cm_cycles;
pub(crate) mod cm_hom;
pub(crate) mod linking;
pub(crate) mod richelot;
pub(crate) mod transfers;
pub(crate) mod twisted_mukai;

pub(crate) enum Command {
    Define(Label, Label, usize, Box<ParsedExpr>, Vec<ParsedExpr>),
    PowerDefine(Label, Label, usize, usize, Box<ParsedExpr>, Vec<ParsedExpr>),
    PowerFrobenius(Label, Label, usize),
    PowerTransition(Label, Label),
    Brauer(Label, Label, usize, bool),
    Frobenius(Label, Label, usize),
    RichelotMatrix(Label, Label),
}

pub(crate) struct Basis {
    pub(crate) curve_name: Rc<str>,
    pub(crate) prime: Prime,
    pub(crate) height: usize,
    pub(crate) order: u64,
    curve: Rc<Curve>,
    generators: Box<[Rc<Divisor>]>,
    root: Scalar,
    logarithm: Logarithm,
    theta: OnceCell<Rc<theta::Group>>,
}

pub(crate) struct Frobenius {
    pub(crate) matrix: Matrix<u64>,
    pub(crate) multiplier: u64,
    pub(crate) characteristic: Vec<u64>,
}

// The four columns of a genus-two level-2 matrix are the generator images.
pub(crate) fn binary_mask_image(matrix: &Matrix<u64>, mask: u8) -> u8 {
    let mut output = 0;
    for row in 0..4 {
        let mut bit = 0;
        for column in 0..4 {
            bit ^= ((mask >> column) & 1) & (*matrix.entry(row, column) as u8);
        }
        output |= bit << row;
    }
    output
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
        Self::new_power(
            model,
            curve_name,
            curve,
            Residues::new(prime, 1)?,
            root,
            generators,
        )
    }

    pub(crate) fn new_power(
        model: &mut Model,
        curve_name: Rc<str>,
        curve: Rc<Curve>,
        level: Residues,
        root: &[u64],
        generators: Vec<Rc<Divisor>>,
    ) -> Result<Self> {
        let Residues {
            prime,
            precision: height,
            modulus: order,
        } = level;
        let genus = curve.grade.genus;
        if generators.len() != 2 * genus {
            return Err(Error::plain(format!(
                "Jacobian torsion basis requires exactly {} generators",
                2 * genus
            )));
        }
        if u64::from(prime.get()) == curve.ring.prime() {
            return Err(Error::plain(
                "Jacobian torsion level must differ from the residue characteristic",
            ));
        }
        let one = curve.ring.scalar(1);
        if curve
            .ring
            .power(model, root, order / u64::from(prime.get()))?
            == one
            || curve.ring.power(model, root, order)? != one
        {
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
        model.allocate(9 + curve_name.len() + generators.len() + root.len())?;
        let logarithm = Logarithm::new(model, &curve.ring, order, root)?;
        Ok(Self {
            curve_name,
            prime,
            height,
            order,
            curve,
            generators: generators.into_boxed_slice(),
            root: root.to_vec(),
            logarithm,
            theta: OnceCell::new(),
        })
    }

    pub(crate) fn len(&self) -> usize {
        self.generators.len()
    }

    pub(crate) fn theta(&self, model: &mut Model) -> Result<Rc<theta::Group>> {
        if let Some(group) = self.theta.get() {
            return Ok(group.clone());
        }
        let group = theta::Group::build(model, self, self.curve.clone())?;
        Ok(self.theta.get_or_init(|| Rc::new(group)).clone())
    }

    pub(crate) fn require_base_extension(&self, model: &mut Model, base: &Curve) -> Result<()> {
        crate::hyperelliptic_jacobian::require_extension(model, base, &self.curve)
    }

    pub(crate) fn binary_combination(&self, model: &mut Model, mut mask: u8) -> Result<Divisor> {
        let mut result = Divisor::zero(&self.curve);
        for generator in &self.generators {
            if mask == 0 {
                break;
            }
            if mask & 1 != 0 {
                result = result.add(model, &self.curve, generator)?;
            }
            mask >>= 1;
        }
        Ok(result)
    }

    pub(crate) fn coordinate(&self, model: &mut Model, t: &Divisor, index: usize) -> Result<u32> {
        self.require_torsion(model, t, "rvJacTorsionCoord")?;
        Ok(self.coordinate_of_torsion(model, t, index)? as u32)
    }

    pub(crate) fn power_coordinate(
        &self,
        model: &mut Model,
        t: &Divisor,
        index: usize,
    ) -> Result<u64> {
        self.require_torsion(model, t, "Jacobian torsion coordinate")?;
        self.coordinate_of_torsion(model, t, index)
    }

    pub(crate) fn translation_coordinates(
        &self,
        model: &mut Model,
        t: &Divisor,
    ) -> Result<Vec<u32>> {
        self.require_torsion(model, t, "generalized Kummer translation")?;
        (0..self.len())
            .map(|index| {
                self.coordinate_of_torsion(model, t, index)
                    .map(|value| value as u32)
            })
            .collect()
    }

    fn require_torsion(&self, model: &mut Model, t: &Divisor, context: &str) -> Result<()> {
        let order = self.order;
        if !is_torsion(model, &self.curve, order, t)? {
            return Err(Error::plain(format!(
                "{context} requires a class killed by {order}"
            )));
        }
        Ok(())
    }

    fn coordinate_of_torsion(&self, model: &mut Model, t: &Divisor, index: usize) -> Result<u64> {
        let genus = self.curve.grade.genus;
        let (left, right) = if index < genus {
            (t, self.generators[genus + index].as_ref())
        } else {
            (self.generators[index - genus].as_ref(), t)
        };
        let order = self.order;
        let value = pairing(model, &self.curve, order, left, right)?;
        self.logarithm
            .evaluate(model, &self.curve.ring, order, value)
    }

    pub(super) fn generators(&self) -> &[Rc<Divisor>] {
        &self.generators
    }

    pub(super) fn polarization_scale(
        &self,
        model: &mut Model,
        source: &Self,
        degree: u64,
    ) -> Result<u64> {
        let value = self.curve.ring.power(model, &source.root, degree)?;
        self.logarithm
            .evaluate(model, &self.curve.ring, self.order, value)
    }

    pub(super) fn checked_coordinates(
        &self,
        model: &mut Model,
        value: &Divisor,
    ) -> Result<Vec<u64>> {
        let curve = &self.curve;
        let mut reconstruction = Divisor::zero(curve);
        let mut coefficients = Vec::with_capacity(self.len());
        for (row, generator) in self.generators.iter().enumerate() {
            let coefficient = self.coordinate_of_torsion(model, value, row)?;
            coefficients.push(coefficient);
            if coefficient != 0 {
                let term = generator.scale(model, curve, i128::from(coefficient))?;
                reconstruction = reconstruction.add(model, curve, &term)?;
            }
        }
        if !reconstruction.compare(model, value)?.is_eq() {
            return Err(Error::plain(
                "Jacobian torsion matrix column failed to reconstruct its divisor class",
            ));
        }
        Ok(coefficients)
    }

    /// Express images of `source` generators in this basis, checking every
    /// column by recombining the actual divisor classes on the target curve.
    pub(crate) fn action_matrix(
        &self,
        model: &mut Model,
        source: &Self,
        mut image_of: impl FnMut(&mut Model, &Divisor) -> Result<Divisor>,
    ) -> Result<Matrix<u64>> {
        let width = self.len();
        if width != source.len() || self.order != source.order || self.prime != source.prime {
            return Err(Error::plain(
                "Jacobian torsion action requires bases of equal rank and level",
            ));
        }
        let ring = Residues::new(self.prime, self.height)?;
        let cells = ring.matrix_size(width)?;
        model.spend(cells)?;
        let mut entries = vec![0; cells];
        for (column, generator) in source.generators.iter().enumerate() {
            let image = image_of(model, generator)?;
            for (row, coefficient) in self
                .checked_coordinates(model, &image)?
                .into_iter()
                .enumerate()
            {
                entries[row * width + column] = coefficient;
            }
        }
        Matrix::from_entries(width, entries)
    }

    pub(crate) fn transition(
        &self,
        model: &mut Model,
        source: &Self,
    ) -> Result<(Matrix<u64>, Matrix<u64>, u64)> {
        if self.curve_name != source.curve_name || self.order != source.order {
            return Err(Error::plain(
                "Jacobian level transition requires the same named curve and level",
            ));
        }
        let matrix = self.action_matrix(model, source, |_, divisor| Ok(divisor.clone()))?;
        let reverse = source.action_matrix(model, self, |_, divisor| Ok(divisor.clone()))?;
        let ring = Residues::new(self.prime, self.height)?;
        let inverse = hw_algebra::inverse(model, &ring, &matrix)?
            .ok_or_else(|| Error::plain("Jacobian level transition is singular"))?;
        if inverse.entries() != reverse.entries() {
            return Err(Error::plain(
                "Jacobian level inverse disagrees with the reverse divisor coordinates",
            ));
        }
        let multiplier =
            self.logarithm
                .evaluate(model, &self.curve.ring, self.order, source.root.clone())?;
        verify_similitude_mod(
            model,
            self.order,
            self.curve.grade.genus,
            multiplier,
            &matrix,
        )?;
        Ok((matrix, reverse, multiplier))
    }

    pub(crate) fn power_frobenius(
        &self,
        model: &mut Model,
        base: &Curve,
        power: usize,
    ) -> Result<(Frobenius, Matrix<u64>)> {
        let direct = self.frobenius(model, base, power)?;
        let ring = Residues::new(self.prime, self.height)?;
        if power != 1 {
            let mut expected = hw_algebra::identity(model, &ring, self.len())?;
            let mut exponent = power;
            if exponent > 0 {
                let mut current = self.frobenius(model, base, 1)?.matrix;
                while exponent > 0 {
                    if exponent & 1 != 0 {
                        expected = hw_algebra::multiply(model, &ring, &expected, &current)?;
                    }
                    exponent >>= 1;
                    if exponent > 0 {
                        current = hw_algebra::multiply(model, &ring, &current, &current)?;
                    }
                }
            }
            if direct.matrix.entries() != expected.entries() {
                return Err(Error::plain(
                    "Jacobian Frobenius power disagrees with matrix composition",
                ));
            }
        }
        let inverse = hw_algebra::inverse(model, &ring, &direct.matrix)?
            .ok_or_else(|| Error::plain("Jacobian torsion Frobenius is singular"))?;
        Ok((direct, inverse))
    }

    pub(crate) fn frobenius(
        &self,
        model: &mut Model,
        base: &Curve,
        power: usize,
    ) -> Result<Frobenius> {
        let curve = &self.curve;
        let p_power =
            crate::hyperelliptic_jacobian::relative_frobenius_power(model, base, curve, power)?;
        let matrix = self.action_matrix(model, self, |model, generator| {
            generator.frobenius_coefficients(model, curve, p_power)
        })?;
        let order = self.order;
        let q = power_mod(model, base.ring.prime(), base.ring.degree(), order)?;
        let multiplier = power_mod(model, q, power, order)?;
        verify_similitude_mod(model, order, curve.grade.genus, multiplier, &matrix)?;
        let ring = Residues::new(self.prime, self.height)?;
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
    babies: HashMap<Scalar, u64>,
    step: u64,
    inverse_step: Scalar,
}

impl Logarithm {
    fn new(model: &mut Model, ring: &Ring, order: u64, root: &[u64]) -> Result<Self> {
        let step = (order - 1).isqrt() + 1;
        let width = ring.degree() + 7;
        let capacity = usize::try_from(step)
            .map_err(|_| Error::plain("Jacobian torsion logarithm exceeds the supported size"))?;
        let entries = capacity
            .checked_add(1)
            .and_then(|entries| entries.checked_mul(width))
            .ok_or_else(|| Error::plain("Jacobian torsion logarithm exceeds the supported size"))?;
        model.allocate(entries)?;
        let mut babies = HashMap::with_capacity(capacity);
        let mut value = ring.scalar(1);
        for j in 0..step {
            let next = ring.multiply(model, &value, root)?;
            model.spend(width)?;
            babies.insert(value, j);
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
    ) -> Result<u64> {
        for block in 0..order.div_ceil(self.step) {
            model.spend(ring.degree() + 1)?;
            if let Some(&baby) = self.babies.get(&value) {
                let index = u128::from(block) * u128::from(self.step) + u128::from(baby);
                return Ok((index % u128::from(order)) as u64);
            }
            value = ring.multiply(model, &value, &self.inverse_step)?;
        }
        Err(Error::plain(
            "Jacobian Weil pairing lies outside the certified roots of unity",
        ))
    }
}

pub(crate) fn power_mod(
    model: &mut Model,
    mut value: u64,
    mut power: usize,
    modulus: u64,
) -> Result<u64> {
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

pub(crate) fn verify_similitude(
    model: &mut Model,
    prime: Prime,
    genus: usize,
    multiplier: u32,
    matrix: &Matrix<u64>,
) -> Result<()> {
    verify_similitude_mod(
        model,
        u64::from(prime.get()),
        genus,
        u64::from(multiplier),
        matrix,
    )
}

fn verify_similitude_mod(
    model: &mut Model,
    modulus: u64,
    genus: usize,
    multiplier: u64,
    matrix: &Matrix<u64>,
) -> Result<()> {
    use crate::rv_unram::{add_mod, sub_mod};
    let width = matrix.size();
    for i in 0..width {
        for j in i + 1..width {
            model.spend(1 + genus)?;
            let mut pair = 0;
            for k in 0..genus {
                let forward = mul_mod(*matrix.entry(k, i), *matrix.entry(k + genus, j), modulus);
                let reverse = mul_mod(*matrix.entry(k + genus, i), *matrix.entry(k, j), modulus);
                pair = add_mod(pair, sub_mod(forward, reverse, modulus), modulus);
            }
            let expected = if j == i + genus { multiplier } else { 0 };
            if pair != expected {
                return Err(Error::plain(
                    "Jacobian torsion matrix failed its symplectic similitude check",
                ));
            }
        }
    }
    Ok(())
}
