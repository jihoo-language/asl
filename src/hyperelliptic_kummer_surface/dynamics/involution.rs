use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_qp_quad::is_coefficient_overflow;

use super::super::{
    nikulin::binary_rref,
    unimodular::{self, Integer, Vector, RANK},
};
use super::{half, incident_tropes, signed, Generator, Word};

impl Word {
    pub(crate) fn apply_lattice(&self, model: &mut Model, input: &Vector) -> Result<Vector> {
        let (w, b) = input.coordinates();
        let mut w = w.map(Integer::from_i128);
        let mut b = b.map(Integer::from_i128);
        // A word can leave the signed 128-bit coordinate range along the
        // way and return to it. Narrow only when forming the final value.
        for &generator in &self.generators {
            generator
                .act_on_lattice(model, self.origin, &mut w, &mut b)
                .map_err(|error| {
                    if is_coefficient_overflow(&error) {
                        Error::plain("K3 lattice automorphism exceeds the exact intermediate range")
                    } else {
                        error
                    }
                })?;
        }
        let mut wedge = [0; 6];
        let mut exceptional = [0; 16];
        for (output, entry) in wedge.iter_mut().zip(w) {
            *output = signed(entry)?;
        }
        for (output, entry) in exceptional.iter_mut().zip(b) {
            *output = signed(entry)?;
        }
        Vector::new(wedge, exceptional)
    }

    pub(super) fn nikulin(&self, model: &mut Model, name: &str) -> Result<String> {
        let basis = unimodular::integral_basis(model)?;
        let images = basis
            .iter()
            .map(|vector| self.apply_lattice(model, vector))
            .collect::<Result<Vec<_>>>()?;
        let mut columns = [[0i128; RANK]; RANK];
        for (column, image) in columns.iter_mut().zip(&images) {
            *column = image.basis_coordinates(model)?;
        }
        verify_involution(model, &columns)?;

        let mut diagonal = Integer::from_i128(0);
        let mut mod_two = [0u32; RANK];
        let mut delta = 0;
        for (i, column) in columns.iter().enumerate() {
            diagonal = diagonal.add(Integer::from_i128(column[i]))?;
            mod_two[i] = column.iter().enumerate().fold(1 << i, |bits, (j, entry)| {
                bits ^ (((entry & 1) as u32) << j)
            });
            // For a projected dual vector (x+Gx)/2, its square is
            // (x²+(x,Gx))/2. Modulo two the cross terms cancel because
            // the involutive isometry G is self-adjoint.
            let parity = basis[i].pair(model, &images[i])? & 1;
            delta |= parity;
            for (j, image) in images[..=i].iter().enumerate() {
                if images[i].wide_pair(model, image)?.signed()
                    != Some(basis[i].pair(model, &basis[j])?)
                {
                    return Err(Error::plain(
                        "K3 automorphism did not preserve the lattice pairing",
                    ));
                }
            }
        }

        // In a unimodular lattice the discriminant of the fixed lattice
        // records the quotient of its projected dual by itself. Its binary
        // dimension is the rank of I+G on the integral lattice modulo two.
        let a = binary_rref(&mut mod_two, RANK).len() as i128;
        let twice_rank = RANK as i128 + signed(diagonal)?;
        if !(0..=2 * RANK as i128).contains(&twice_rank) || twice_rank % 2 != 0 {
            return Err(inconsistent_invariants());
        }
        let r = twice_rank / 2;
        let (kind, locus) = fixed_locus(self.reverses_primitive(), r, a, delta)?;
        Ok(format!(
            "rvK3AutoNikulin[{name}] = {kind}, r={r}, a={a}, delta={delta}, discriminant=(Z/2)^{a}, fixed locus={locus}"
        ))
    }
}

impl Generator {
    fn act_on_lattice(
        self,
        model: &mut Model,
        origin: u16,
        w: &mut [Integer; 6],
        b: &mut [Integer; 16],
    ) -> Result<()> {
        // The hyperplane coefficient of the marked projection is
        // h=(w[1]+w[4])/2. Switches and projections negate the primitive
        // part, so add h'+h in the two hyperplane wedge positions.
        let shift = match self {
            Self::Translate(mask) => {
                model.spend(RANK)?;
                *b = std::array::from_fn(|node| b[node ^ usize::from(mask)]);
                return Ok(());
            }
            Self::Project(mask) => {
                model.spend(RANK * 13)?;
                let hyperplane_twice = w[1].add(w[4])?;
                let node = usize::from(mask);
                let old = b[node];
                let shift = hyperplane_twice.add(old)?.multiply(2)?;
                b[node] = hyperplane_twice.multiply(-2)?.add(old.multiply(-3)?)?;
                shift
            }
            Self::Switch(mask) => {
                model.spend(RANK * 16 * 13)?;
                let hyperplane_twice = w[1].add(w[4])?;
                let total = b
                    .iter()
                    .try_fold(Integer::from_i128(0), |sum, &entry| sum.add(entry))?;
                let old = *b;
                for (node, output) in b.iter_mut().enumerate() {
                    let sum = incident_tropes(origin, mask, node)
                        .try_fold(hyperplane_twice, |sum, trope| sum.add(old[trope]))?;
                    *output = half(sum)?.negate();
                }
                hyperplane_twice.multiply(2)?.add(half(total)?)?
            }
        };
        *w = w.map(Integer::negate);
        w[1] = w[1].add(shift)?;
        w[4] = w[4].add(shift)?;
        Ok(())
    }
}

fn verify_involution(model: &mut Model, columns: &[[i128; RANK]; RANK]) -> Result<()> {
    model.spend(RANK * RANK * RANK * 13)?;
    for (i, column) in columns.iter().enumerate() {
        for j in 0..RANK {
            let square = columns
                .iter()
                .zip(column)
                .try_fold(Integer::from_i128(0), |sum, (other, &coefficient)| {
                    sum.add(Integer::product(other[j], coefficient))
                })?;
            if square.signed() != Some(i128::from(i == j)) {
                return Err(Error::plain(
                    "rvK3AutoNikulin requires a K3 automorphism whose square is the identity",
                ));
            }
        }
    }
    Ok(())
}

fn fixed_locus(
    antisymplectic: bool,
    r: i128,
    a: i128,
    delta: i128,
) -> Result<(&'static str, String)> {
    if !antisymplectic {
        return match (r, a, delta) {
            (22, 0, 0) => Ok(("identity", "whole K3 surface".into())),
            (14, 8, 0) => Ok(("symplectic", "8 isolated points".into())),
            _ => Err(inconsistent_invariants()),
        };
    }
    let locus = match (r, a, delta) {
        (10, 10, 0) => "empty".into(),
        (10, 8, 0) => "2 elliptic curves".into(),
        _ => {
            let twice_genus = 22 - r - a;
            let twice_extra = r - a;
            if twice_genus < 0 || twice_extra < 0 || (twice_genus | twice_extra) % 2 != 0 {
                return Err(inconsistent_invariants());
            }
            let genus = twice_genus / 2;
            let extra = twice_extra / 2;
            match (genus, extra) {
                (0, 0) => "1 rational curve".into(),
                (0, _) => format!("{} rational curves", extra + 1),
                (_, 0) => format!("1 curve of genus {genus}"),
                (_, 1) => format!("1 curve of genus {genus} and 1 rational curve"),
                _ => format!("1 curve of genus {genus} and {extra} rational curves"),
            }
        }
    };
    Ok(("antisymplectic", locus))
}

fn inconsistent_invariants() -> Error {
    Error::plain("K3 involution has inconsistent invariant lattice data")
}
