use std::rc::Rc;

use super::brauer::{identity, smith_frame, Exterior};
use super::linking::{column, image, product, Fixed, Pairing};
use crate::error::{Error, Result};
use crate::hyperelliptic::Residues;
use crate::hyperelliptic_kummer_surface::k3;
use crate::model::Model;
use crate::rv_unram::{add_mod, mul_mod, sub_mod};
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Statement {
    Define(Label, Label, usize, Option<Label>),
    Data(Label),
}

pub(crate) enum Expression {
    Class(Label, Box<ParsedExpr>),
    Algebraic(Label, Box<ParsedExpr>),
    Generator(Label, usize),
    Lift(Box<ParsedExpr>),
    Transfer(Label, Direction, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Direction {
    Restrict,
    Corestrict,
}

pub(crate) type Class = [u64; 6];

pub(crate) struct Quotient {
    pub(crate) pairing: Rc<Pairing>,
    pub(crate) pairing_name: Label,
    pub(crate) level: usize,
    modulus: u64,
    free: usize,
    divisorial_degree: usize,
    basis: Exterior,
    inverse: Exterior,
    steps: [u64; 6],
    orders: [u64; 6],
}

impl Quotient {
    pub(crate) fn new(
        model: &mut Model,
        pairing_name: Label,
        pairing: Rc<Pairing>,
        level: usize,
        reference: Option<&Pairing>,
    ) -> Result<Self> {
        let modulus = pairing.modulus(level)?;
        let reference = reference.unwrap_or(&pairing);
        extension_ratio(reference, &pairing)?;
        let divisorial_degree = reference.degree;
        let curve = model.kummer_curve(&pairing.base.text)?;
        let q = curve.ring.size_of_residue_field()?;
        let weil = curve.weil(model)?;
        let free = k3::rational_abelian_rank(model, q, weil, pairing.degree as u128)?;
        let reference_free = k3::rational_abelian_rank(model, q, weil, reference.degree as u128)?;
        if free != reference_free {
            return Err(Error::plain(
                "Kummer Brauer reference Frobenius has a different abelian fixed rank",
            ));
        }
        let finite = 6usize.checked_sub(free).ok_or_else(inconsistent)?;
        let stable = smith_frame(model, reference.ring, reference.difference, 6)?;
        let height = reference.ring.precision;
        if stable.exponents[finite..].iter().any(|&e| e != height) {
            return Err(inconsistent());
        }
        let max = stable.exponents[..finite].last().copied().unwrap_or(0);
        if max >= height || level + max > height {
            return Err(Error::plain("Kummer Brauer quotient needs more marking height to recover the integral fixed lattice"));
        }
        // From the earlier power with the same fixed rank, its saturated
        // integral fixed lattice is also the present integral fixed lattice.
        // The last columns of this Smith basis project to its full reduction.
        model.spend(6usize.pow(3))?;
        let defect = product(&pairing.difference, &stable.basis, modulus);
        if defect
            .iter()
            .any(|row| row[finite..].iter().any(|&entry| entry != 0))
        {
            return Err(inconsistent());
        }

        // Diagonalize the current defect on the complementary columns.
        // Zero modulo this level is enough here: free columns are already
        // separated from finite factors by the earlier integral certificate.
        let ring = Residues::new(pairing.ring.prime, level)?;
        let local = smith_frame(model, ring, defect, finite)?;
        model.spend(2 * 6usize.pow(3))?;
        let basis = product(&stable.basis, &local.basis, modulus);
        let inverse = product(&local.inverse, &stable.inverse, modulus);
        let prime = u64::from(ring.prime.get());
        let mut steps = [1; 6];
        let mut orders = [1; 6];
        for (i, &exponent) in local.exponents[..finite].iter().enumerate() {
            steps[i] = prime.pow((level - exponent) as u32);
            orders[i] = prime.pow(exponent as u32);
        }
        model.allocate(110 + pairing_name.text.len())?;
        Ok(Self {
            pairing,
            pairing_name,
            level,
            modulus,
            free,
            divisorial_degree,
            basis,
            inverse,
            steps,
            orders,
        })
    }

    pub(crate) fn project(&self, model: &mut Model, fixed: Fixed) -> Result<Class> {
        self.pairing.require_fixed(model, self.modulus, fixed)?;
        model.spend(36)?;
        let frame = image(&self.inverse, fixed, self.modulus);
        let mut coords = [0; 6];
        for (i, c) in coords.iter_mut().enumerate().take(6 - self.free) {
            if !frame[i].is_multiple_of(self.steps[i]) {
                return Err(inconsistent());
            }
            *c = frame[i] / self.steps[i];
        }
        Ok(coords)
    }

    pub(crate) fn lift(&self, model: &mut Model, class: Class) -> Result<Fixed> {
        model.spend(36)?;
        let frame = std::array::from_fn(|i| self.steps[i] * class[i]);
        let lift = image(&self.basis, frame, self.modulus);
        self.pairing.require_fixed(model, self.modulus, lift)?;
        Ok(lift)
    }

    pub(crate) fn add(&self, left: Class, right: Class) -> Class {
        std::array::from_fn(|i| add_mod(left[i], right[i], self.orders[i]))
    }

    pub(crate) fn scale(&self, class: Class, coefficient: i128) -> Class {
        std::array::from_fn(|i| {
            let residue = coefficient.rem_euclid(i128::from(self.orders[i])) as u64;
            mul_mod(class[i], residue, self.orders[i])
        })
    }

    pub(super) fn order(&self, class: Class) -> u64 {
        // Every invariant factor is a power of the same prime.
        self.orders
            .into_iter()
            .zip(class)
            .map(|(modulus, coefficient)| {
                modulus / crate::rv_ore::gcd(u128::from(modulus), u128::from(coefficient)) as u64
            })
            .max()
            .unwrap_or(1)
    }

    pub(crate) fn generator(&self, index: usize) -> Result<Class> {
        let slot = (0..6)
            .filter(|&i| self.orders[i] > 1)
            .nth(index)
            .ok_or_else(|| {
                Error::plain("Kummer Brauer quotient generator index is out of range")
            })?;
        let mut result = [0; 6];
        result[slot] = 1;
        Ok(result)
    }

    pub(crate) fn cardinality_bound(&self) -> usize {
        self.orders
            .iter()
            .try_fold(1usize, |count, &order| {
                count.checked_mul(usize::try_from(order).ok()?)
            })
            .unwrap_or(usize::MAX)
    }

    pub(crate) fn elements(&self, model: &mut Model, count: usize) -> Result<Vec<Class>> {
        model.spend(count.saturating_mul(6))?;
        model.allocate(count.saturating_mul(6))?;
        Ok((0..count)
            .map(|mut code| {
                self.orders.map(|order| {
                    let value = (code as u64) % order;
                    code /= order as usize;
                    value
                })
            })
            .collect())
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str) -> Result<String> {
        let divisors: Vec<_> = (6 - self.free..6).map(|i| column(&self.basis, i)).collect();
        let orders: Vec<_> = self
            .orders
            .iter()
            .copied()
            .filter(|&order| order > 1)
            .collect();
        let lifts = (0..orders.len())
            .map(|j| self.lift(model, self.generator(j)?))
            .collect::<Result<Vec<_>>>()?;
        Ok(format!(
            "rvKummerBrauerQuotientData[{name}] = level={}, divisor_rank={}, reference_degree={}, invariant_factors={orders:?}, generator_lifts={lifts:?}, divisor_basis={divisors:?}",
            self.modulus, self.free, self.divisorial_degree,
        ))
    }
}

pub(crate) struct Transfer {
    matrix: Exterior,
}

impl Transfer {
    pub(crate) fn new(
        model: &mut Model,
        from: &Pairing,
        to: &Pairing,
        direction: Direction,
    ) -> Result<Self> {
        let (base, extension) = match direction {
            Direction::Restrict => (from, to),
            Direction::Corestrict => (to, from),
        };
        let ratio = extension_ratio(base, extension)?;
        let m = base.ring.modulus;
        let action = action(base, m);
        let (power, norm) = power_and_sum(model, action, ratio, m)?;
        if power != self::action(extension, m) {
            return Err(Error::plain(
                "Kummer Brauer extension Frobenius does not match the matrix power",
            ));
        }
        let matrix = match direction {
            Direction::Restrict => identity(),
            Direction::Corestrict => norm,
        };
        model.allocate(36)?;
        Ok(Self { matrix })
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        to: &Pairing,
        level: usize,
        fixed: Fixed,
    ) -> Result<Fixed> {
        let m = to.modulus(level)?;
        model.spend(36)?;
        let image = image(&self.matrix, fixed, m);
        to.require_fixed(model, m, image)?;
        Ok(image)
    }
}

pub(super) fn action(pairing: &Pairing, modulus: u64) -> Exterior {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            sub_mod(
                u64::from(i == j),
                pairing.difference[i][j] % modulus,
                modulus,
            )
        })
    })
}

pub(super) fn extension_ratio(base: &Pairing, extension: &Pairing) -> Result<usize> {
    if base.base.text != extension.base.text
        || !Rc::ptr_eq(&base.basis, &extension.basis)
        || !extension.degree.is_multiple_of(base.degree)
    {
        return Err(Error::plain(
            "Kummer Brauer transfer requires the same base and marking with nested degrees",
        ));
    }
    Ok(extension.degree / base.degree)
}

// (U^d, 1+U+...+U^(d-1)) computed in the matrix monoid.
fn power_and_sum(
    model: &mut Model,
    action: Exterior,
    mut exponent: usize,
    modulus: u64,
) -> Result<(Exterior, Exterior)> {
    let mut result = (identity(), [[0; 6]; 6]);
    let mut current = (action, identity());
    while exponent != 0 {
        if exponent & 1 == 1 {
            result = compose(model, &result, &current, modulus)?;
        }
        exponent >>= 1;
        if exponent > 0 {
            current = compose(model, &current, &current, modulus)?;
        }
    }
    Ok(result)
}

fn compose(
    model: &mut Model,
    a: &(Exterior, Exterior),
    b: &(Exterior, Exterior),
    m: u64,
) -> Result<(Exterior, Exterior)> {
    model.spend(2 * 6usize.pow(3) + 36)?;
    let power = product(&a.0, &b.0, m);
    let mut sum = product(&a.0, &b.1, m);
    for (row, previous) in sum.iter_mut().zip(a.1) {
        for (entry, previous) in row.iter_mut().zip(previous) {
            *entry = add_mod(*entry, previous, m);
        }
    }
    Ok((power, sum))
}

fn inconsistent() -> Error {
    Error::plain("Kummer Brauer quotient disagrees with the certified integral fixed lattice")
}
