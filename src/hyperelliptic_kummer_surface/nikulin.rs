use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::jacobian_level::{binary_mask_image, Basis};
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};

const PAIRS: [(u8, u8); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
pub(crate) const ORDER: u8 = 1 << PAIRS.len();
const NODES: usize = 16;

// The truth table of coordinate x_i has bit v set exactly when v_i=1.
const AFFINE: [u16; 5] = [u16::MAX, 0xaaaa, 0xcccc, 0xf0f0, 0xff00];
const GEOMETRIC_DETERMINANT: i128 = (-2i128).pow(NODES as u32) / (1i128 << (2 * AFFINE.len()));
const QUADRATIC: [u16; 6] = [
    AFFINE[1] & AFFINE[2],
    AFFINE[1] & AFFINE[3],
    AFFINE[1] & AFFINE[4],
    AFFINE[2] & AFFINE[3],
    AFFINE[2] & AFFINE[4],
    AFFINE[3] & AFFINE[4],
];

pub(crate) enum Expression {
    Element(Label, u8),
    Quadratic(Box<ParsedExpr>),
    Pair(Box<[ParsedExpr; 2]>),
    Frobenius(Label, usize, Box<ParsedExpr>),
}

pub(crate) struct Command(pub(crate) Label, pub(crate) Label, pub(crate) usize);

pub(crate) fn basis(model: &Model, name: &Label) -> Result<Rc<Basis>> {
    basis_for(
        model,
        name,
        "Kummer discriminants require a certified genus-two level-2 basis",
    )
}

pub(crate) fn basis_for(
    model: &Model,
    name: &Label,
    diagnostic: &'static str,
) -> Result<Rc<Basis>> {
    let basis = model.jacobian_level(name)?;
    if basis.prime.get() != 2 || basis.len() != 4 {
        return Err(Error::new(name.span, diagnostic));
    }
    Ok(basis)
}

fn combine<T>(mask: u8, basis: &[T]) -> T
where
    T: Copy + Default + std::ops::BitXor<Output = T>,
{
    basis
        .iter()
        .enumerate()
        .filter(|(i, _)| mask & (1 << i) != 0)
        .fold(T::default(), |sum, (_, &term)| sum ^ term)
}

pub(crate) fn quadratic(mask: u8) -> u32 {
    // The lift (Σ f(v)E_v)/2 has square −wt(f)/2 in Z/2Z.
    (combine(mask, &QUADRATIC).count_ones() / 2) & 1
}

pub(crate) fn pair(left: u8, right: u8) -> u32 {
    quadratic(left ^ right) ^ quadratic(left) ^ quadratic(right)
}

pub(crate) fn affine_code() -> impl Iterator<Item = u16> {
    (0..1u8 << AFFINE.len()).map(|mask| combine(mask, &AFFINE))
}

pub(crate) struct NodePermutation([u8; NODES]);

impl NodePermutation {
    pub(crate) fn new(
        model: &mut Model,
        basis: &Basis,
        base: &Curve,
        power: usize,
    ) -> Result<Self> {
        let frobenius = basis.frobenius(model, base, power)?;
        model.spend(NODES * NODES)?;
        let images = std::array::from_fn(|node| binary_mask_image(&frobenius.matrix, node as u8));
        let seen = images.iter().fold(0u16, |seen, &image| seen | (1 << image));
        if seen != u16::MAX {
            return Err(Error::plain(
                "certified Frobenius did not permute the Kummer nodes",
            ));
        }
        Ok(Self(images))
    }

    pub(crate) fn image(&self, node: usize) -> u8 {
        self.0[node]
    }

    pub(crate) fn on_truth_table(&self, table: u16) -> u16 {
        // Carry the coefficient of E_v to E_{Mv}. At a given node this is
        // the pullback of the function by M^−1.
        self.0.iter().enumerate().fold(0, |image, (node, &target)| {
            image | (((table >> node) & 1) << target)
        })
    }

    fn orbit_lengths(&self) -> Vec<usize> {
        let mut orbits = Vec::new();
        let mut seen = 0u16;
        for start in 0..NODES {
            let mut node = start;
            let mut length = 0;
            while seen & (1 << node) == 0 {
                seen |= 1 << node;
                length += 1;
                node = usize::from(self.image(node));
            }
            if length != 0 {
                orbits.push(length);
            }
        }
        orbits.sort_unstable();
        orbits
    }
}

pub(crate) struct Action([u8; 6]);

impl Action {
    pub(crate) fn new(
        model: &mut Model,
        basis: &Basis,
        base: &Curve,
        power: usize,
    ) -> Result<Self> {
        let permutation = NodePermutation::new(model, basis, base, power)?;
        model.spend(QUADRATIC.len() * (NODES + PAIRS.len()))?;
        let images =
            QUADRATIC.map(|table| quadratic_coefficients(permutation.on_truth_table(table)));
        Ok(Self(images))
    }

    pub(crate) fn apply(&self, mask: u8) -> u8 {
        combine(mask, &self.0)
    }
}

fn quadratic_coefficients(table: u16) -> u8 {
    let value = |node| ((table >> node) & 1u16) as u8;
    PAIRS.iter().enumerate().fold(0, |mask, (index, &(i, j))| {
        let first = 1u8 << i;
        let second = 1u8 << j;
        let coefficient = value(0) ^ value(first) ^ value(second) ^ value(first | second);
        mask | (coefficient << index)
    })
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<String> {
        let Self(name, base, power) = self;
        let basis = basis(model, &name)?;
        let curve = model
            .jacobian_curve(&base.text)
            .map_err(|error| error.at(base.span))?;
        let permutation = NodePermutation::new(model, &basis, &curve, power)
            .map_err(|error| error.at(base.span))?;
        let orbits = permutation.orbit_lengths();
        let rank = orbits.len();
        model.spend(NODES * (1 << AFFINE.len()))?;
        let invariants = affine_code()
            .filter(|&table| permutation.on_truth_table(table) == table)
            .count();
        if !invariants.is_power_of_two() {
            return Err(Error::plain(
                "Kummer affine-code invariants do not form a binary subspace",
            ));
        }
        let dimension = invariants.ilog2();
        let orbit_determinant =
            (-2i128).pow(rank as u32) * orbits.iter().map(|&n| n as i128).product::<i128>();
        // Each invariant affine codeword lifts to a fixed half-sum of orbits.
        // It enlarges the orbit-sum lattice by precisely the code index.
        let index_squared = (invariants as i128).pow(2);
        if orbit_determinant % index_squared != 0 {
            return Err(Error::plain(
                "Kummer fixed lattice determinant division is not exact",
            ));
        }
        let determinant = orbit_determinant / index_squared;
        Ok(format!(
            "rvKummerNodeLattice[{}, {}, {power}] = orbits={orbits:?}, rank={rank}, affine invariants={dimension}, index={invariants}, determinant={determinant}, geometric rank={NODES}, geometric determinant={GEOMETRIC_DETERMINANT}",
            name.text, base.text,
        ))
    }
}
