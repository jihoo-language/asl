use super::Basis;
use crate::error::{Error, Result};
use crate::hw_algebra::{CoefficientRing, Matrix};
use crate::hw_closed::Curve;
use crate::hyperelliptic::Residues;
use crate::hyperelliptic_kummer_surface::k3;
use crate::model::Model;
use crate::rv_unram::{add_mod, mul_mod, sub_mod};

pub(super) const PAIRS: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
pub(super) type Exterior = [[u64; 6]; 6];

pub(crate) fn report(
    model: &mut Model,
    basis: &Basis,
    base: &Curve,
    name: &str,
    base_name: &str,
    degree: usize,
    with_artin_tate: bool,
) -> Result<String> {
    if basis.prime.get() == 2 {
        return Err(Error::plain(
            "Kummer K3 Brauer calculation requires an odd torsion prime",
        ));
    }
    if degree == 0 {
        return Err(Error::plain(
            "Kummer K3 Brauer extension degree must be positive",
        ));
    }
    let frobenius = basis.frobenius(model, base, degree)?;
    if frobenius.matrix.size() != 4 {
        return Err(Error::plain(
            "Kummer K3 Brauer calculation requires a genus-two marking",
        ));
    }

    let q = base.ring.size_of_residue_field()?;
    let orbits = k3::node_orbits(model, base, q)?;
    let degree_wide = degree as u128;
    let weil = base.weil(model)?;
    let (rank, _) = k3::picard_ranks(model, q, weil, &orbits, degree_wide)?;
    let exceptional_rank = k3::exceptional_picard_rank(&orbits, degree_wide);
    let free_rank = (rank - exceptional_rank) as usize;
    let ring = Residues::new(basis.prime, basis.height)?;
    let matrix = twisted_exterior_square(&frobenius.matrix, ring.modulus, frobenius.multiplier);
    let smith = smith_exponents(model, ring, matrix)?;
    let exponents = finite_factors(&smith, free_rank, ring.precision)?;
    let prime = u64::from(ring.prime.get());
    let factors: Vec<_> = exponents.iter().map(|&e| prime.pow(e as u32)).collect();
    let complete = !exponents.contains(&ring.precision);
    let query = if with_artin_tate {
        "rvKummerK3BrauerArtinTate"
    } else {
        "rvKummerK3Brauer"
    };
    let mut output = format!(
        "{query}[{name}, {base_name}, {degree}] = level={}, prime={prime}, height={}, exceptional_rank={exceptional_rank}, abelian_free_rank={free_rank}, smith_exponents={smith:?}, brauer_invariant_factors={factors:?}, complete_primary={complete}",
        ring.modulus, ring.precision,
    );
    if with_artin_tate {
        let magnitude = k3::artin_tate_magnitude(model, q, weil, &orbits, degree_wide)?;
        let mut remaining = magnitude;
        let mut upper = 0;
        while remaining.is_multiple_of(u128::from(prime)) {
            model.spend(1)?;
            remaining /= u128::from(prime);
            upper += 1;
        }
        let detected: usize = exponents.iter().sum();
        if detected > upper {
            return Err(Error::plain(
                "Kummer K3 Brauer level order exceeds the Artin–Tate bound",
            ));
        }
        let certified = complete || detected == upper;
        output.push_str(&format!(
            ", artin_tate_product_magnitude={magnitude}, artin_tate_valuation={upper}, detected_order_valuation={detected}, complete_with_artin_tate={certified}"
        ));
    }
    Ok(output)
}

pub(super) fn exterior_square(matrix: &Matrix<u64>, modulus: u64) -> Exterior {
    std::array::from_fn(|i| {
        let (r, s) = PAIRS[i];
        std::array::from_fn(|j| {
            let (u, v) = PAIRS[j];
            let forward = mul_mod(*matrix.entry(r, u), *matrix.entry(s, v), modulus);
            let reverse = mul_mod(*matrix.entry(r, v), *matrix.entry(s, u), modulus);
            sub_mod(forward, reverse, modulus)
        })
    })
}

fn twisted_exterior_square(matrix: &Matrix<u64>, modulus: u64, multiplier: u64) -> Exterior {
    let mut exterior = exterior_square(matrix, modulus);
    for (i, row) in exterior.iter_mut().enumerate() {
        row[i] = sub_mod(row[i], multiplier, modulus);
    }
    exterior
}

pub(super) fn identity() -> Exterior {
    std::array::from_fn(|i| std::array::from_fn(|j| u64::from(i == j)))
}

pub(super) struct Smith {
    pub(super) exponents: [usize; 6],
    pub(super) basis: Exterior,
    pub(super) inverse: Exterior,
}

fn smith_exponents(model: &mut Model, ring: Residues, rows: Exterior) -> Result<[usize; 6]> {
    Ok(smith_frame(model, ring, rows, 6)?.exponents)
}

pub(super) fn smith_frame(
    model: &mut Model,
    ring: Residues,
    mut rows: Exterior,
    width: usize,
) -> Result<Smith> {
    let prime = u64::from(ring.prime.get());
    let mut result = [ring.precision; 6];
    let mut basis = identity();
    let mut inverse_basis = identity();
    for (i, slot) in result.iter_mut().enumerate().take(width) {
        let mut pivot = None;
        for (row, items) in rows.iter().enumerate().skip(i) {
            for (col, &entry) in items.iter().enumerate().take(width).skip(i) {
                model.spend(1)?;
                if entry != 0 {
                    let exponent = valuation(entry, prime);
                    if pivot.is_none_or(|(best, _, _)| exponent < best) {
                        pivot = Some((exponent, row, col));
                    }
                }
            }
        }
        let Some((exponent, row, col)) = pivot else {
            break;
        };
        *slot = exponent;
        rows.swap(i, row);
        for (row, transform) in rows.iter_mut().zip(&mut basis) {
            row.swap(i, col);
            transform.swap(i, col);
        }
        inverse_basis.swap(i, col);
        let divisor = prime.pow(exponent as u32);
        let pivot_row = rows[i];
        let inverse = ring.inverse_unit(model, &(pivot_row[i] / divisor))?;
        // Every entry in the trailing block is divisible by this pivot's
        // prime power. Clear its column by integral row operations.
        for row in rows.iter_mut().skip(i + 1) {
            let quotient = mul_mod(row[i] / divisor, inverse, ring.modulus);
            for (entry, &leading) in row.iter_mut().zip(&pivot_row).skip(i + 1) {
                model.spend(1)?;
                let product = mul_mod(quotient, leading, ring.modulus);
                *entry = sub_mod(*entry, product, ring.modulus);
            }
            row[i] = 0;
        }
        // The cleared column is zero off its pivot. The column operations
        // only change this row, and we track their domain basis and inverse.
        for (j, &entry) in pivot_row.iter().enumerate().take(width).skip(i + 1) {
            let quotient = mul_mod(entry / divisor, inverse, ring.modulus);
            rows[i][j] = 0;
            if quotient == 0 {
                continue;
            }
            model.spend(12)?;
            for row in &mut basis {
                row[j] = sub_mod(
                    row[j],
                    mul_mod(quotient, row[i], ring.modulus),
                    ring.modulus,
                );
            }
            let source = inverse_basis[j];
            for (entry, source) in inverse_basis[i].iter_mut().zip(source) {
                *entry = add_mod(
                    *entry,
                    mul_mod(quotient, source, ring.modulus),
                    ring.modulus,
                );
            }
        }
    }
    Ok(Smith {
        exponents: result,
        basis,
        inverse: inverse_basis,
    })
}

fn valuation(mut nonzero: u64, prime: u64) -> usize {
    let mut exponent = 0;
    while nonzero.is_multiple_of(prime) {
        nonzero /= prime;
        exponent += 1;
    }
    exponent
}

fn finite_factors(smith: &[usize; 6], free_rank: usize, height: usize) -> Result<Vec<usize>> {
    let end = smith.len().checked_sub(free_rank).ok_or_else(rank_error)?;
    if smith[end..].iter().any(|&e| e != height) {
        return Err(rank_error());
    }
    Ok(smith[..end].iter().copied().filter(|&e| e != 0).collect())
}

fn rank_error() -> Error {
    Error::plain("Kummer K3 torsion matrix disagrees with the rational Picard rank")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::odd_cobar::Prime;

    #[test]
    fn nonunit_pivots_and_free_factors_over_a_prime_power() {
        let ring = Residues::new(Prime::validated(3), 3)
            .ok()
            .expect("modulus 27");
        let n = ring.modulus;
        let mut diagonally_equivalent = [[0; 6]; 6];
        diagonally_equivalent[5][1] = 2;
        diagonally_equivalent[4][2] = 3;
        diagonally_equivalent[1][0] = 18;
        let add = crate::rv_unram::add_mod;
        let source = diagonally_equivalent[5];
        for (entry, addend) in diagonally_equivalent[4].iter_mut().zip(source) {
            *entry = add(*entry, addend, n);
        }
        for row in &mut diagonally_equivalent {
            row[1] = add(row[1], row[2], n);
            row[2] = add(row[2], row[0], n);
        }
        let mut model = Model::new();
        let smith = smith_exponents(&mut model, ring, diagonally_equivalent)
            .ok()
            .expect("valid Smith reduction");
        assert_eq!(smith, [0, 1, 2, 3, 3, 3]);
        assert_eq!(
            finite_factors(&smith, 2, 3)
                .ok()
                .expect("remove two free factors"),
            vec![1, 2, 3]
        );
        assert_eq!(
            finite_factors(&smith, 3, 3)
                .ok()
                .expect("remove three free factors"),
            vec![1, 2]
        );
        assert!(finite_factors(&smith, 4, 3).is_err());
    }

    #[test]
    fn the_domain_frame_recovers_mixed_torsion_and_noncoordinate_free_columns() {
        use crate::jacobian_level::linking::{column, image, product};

        let ring = Residues::new(Prime::validated(3), 4)
            .ok()
            .expect("modulus 81");
        let mut defect = [[0; 6]; 6];
        defect[2] = [9, 0, 9, 0, 0, 0];
        defect[4] = [0, 6, 6, 0, 0, 0];
        let mut model = Model::new();
        let smith = smith_frame(&mut model, ring, defect, 6)
            .ok()
            .expect("mixed finite factors");
        assert_eq!(smith.exponents, [1, 2, 4, 4, 4, 4]);
        assert_eq!(
            product(&smith.basis, &smith.inverse, ring.modulus),
            identity()
        );
        assert_eq!(
            image(&defect, column(&smith.basis, 2), ring.modulus),
            [0; 6]
        );
        assert_eq!(column(&smith.basis, 2), [80, 80, 1, 0, 0, 0]);

        let level = 9;
        let coordinates = [6, 8, 7, 0, 0, 0];
        let representative = image(&smith.basis, coordinates, level);
        assert_eq!(image(&defect, representative, level), [0; 6]);
        assert_eq!(image(&smith.inverse, representative, level), coordinates);
        assert_eq!(coordinates[0] / 3, 2);
        assert_eq!(coordinates[1], 8);
    }
}
