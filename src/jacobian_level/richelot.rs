use crate::error::{Error, Result, Span};
use crate::hw_algebra::{self, Matrix};
use crate::hyperelliptic::Residues;
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear;
use crate::richelot::correspondence::Correspondence;
use crate::syntax::Label;

use super::{verify_similitude, Basis};

pub(crate) struct Report {
    pub(crate) prime: Prime,
    pub(crate) root_exponent: u32,
    pub(crate) multiplier: u32,
    pub(crate) rank: usize,
    pub(crate) matrix: Matrix<u64>,
    pub(crate) adjoint: Matrix<u64>,
}

impl Basis {
    pub(crate) fn richelot_matrix(
        &self,
        model: &mut Model,
        target: &Self,
        target_span: Span,
    ) -> Result<Report> {
        if self.curve.grade.genus != 2 || target.curve.grade.genus != 2 {
            return Err(Error::plain(
                "Richelot torsion matrices require genus-two bases",
            ));
        }
        if self.prime != target.prime {
            return Err(Error::plain(
                "Richelot torsion bases must have the same prime level",
            ));
        }
        let correspondence = Correspondence::resolve(
            model,
            &Label {
                text: target.curve_name.to_string(),
                span: target_span,
            },
        )?;
        if self.curve_name.as_ref() != correspondence.source_name() {
            return Err(Error::plain(format!(
                "Richelot torsion matrix requires a basis of its named source curve '{}'",
                correspondence.source_name()
            )));
        }
        if self.curve.grade.ring != target.curve.grade.ring {
            return Err(Error::plain(
                "Richelot torsion bases must use the same named field",
            ));
        }

        let prime = self.prime;
        let order = u64::from(prime.get());
        model.spend(self.root.len())?;
        let root_exponent =
            target
                .logarithm
                .evaluate(model, &target.curve.ring, order, self.root.clone())?;
        let inverse_exponent = prime_linear::invert(model, prime, root_exponent)?;
        let two = prime.add(1, 1);
        let multiplier = prime.mul(two, root_exponent);
        let matrix = target.action_matrix(model, self, |model, divisor| {
            correspondence.push(model, divisor)
        })?;
        verify_similitude(model, prime, 2, multiplier, &matrix)?;
        let adjoint = polarized_adjoint(model, prime, inverse_exponent, &matrix)?;
        check_compositions(model, prime, two, &matrix, &adjoint)?;
        let rank = matrix_rank(model, prime, &matrix)?;
        let expected_rank = if order == 2 { 2 } else { 4 };
        if rank != expected_rank {
            return Err(Error::plain(format!(
                "Richelot torsion matrix has rank {rank}, expected {expected_rank}"
            )));
        }
        Ok(Report {
            prime,
            root_exponent,
            multiplier,
            rank,
            matrix,
            adjoint,
        })
    }
}

fn polarized_adjoint(
    model: &mut Model,
    prime: Prime,
    inverse_exponent: u32,
    matrix: &Matrix<u64>,
) -> Result<Matrix<u64>> {
    let width = matrix.size();
    let genus = width / 2;
    let mut entries = Vec::with_capacity(matrix.entries().len());
    model.spend(matrix.entries().len())?;
    // J^{-1} M^t J swaps both halves, transposes, and negates the off-diagonal blocks.
    for i in 0..width {
        for j in 0..width {
            let value = *matrix.entry((j + genus) % width, (i + genus) % width) as u32;
            let value = if (i < genus) == (j < genus) {
                value
            } else {
                prime.negate(value)
            };
            entries.push(u64::from(prime.mul(inverse_exponent, value)));
        }
    }
    Matrix::from_entries(width, entries)
}

fn check_compositions(
    model: &mut Model,
    prime: Prime,
    two: u32,
    matrix: &Matrix<u64>,
    adjoint: &Matrix<u64>,
) -> Result<()> {
    let ring = Residues::new(prime, 1)?;
    let width = matrix.size();
    let mut scalar = vec![0; matrix.entries().len()];
    model.spend(scalar.len())?;
    for i in 0..width {
        scalar[i * (width + 1)] = u64::from(two);
    }
    let expected = Matrix::from_entries(width, scalar)?;
    let left = hw_algebra::multiply(model, &ring, adjoint, matrix)?;
    let right = hw_algebra::multiply(model, &ring, matrix, adjoint)?;
    if left != expected || right != expected {
        return Err(Error::plain(
            "Richelot matrix and adjoint must compose to multiplication by 2",
        ));
    }
    Ok(())
}

fn matrix_rank(model: &mut Model, prime: Prime, matrix: &Matrix<u64>) -> Result<usize> {
    let mut basis = prime_linear::Basis::default();
    for row in matrix.entries().chunks_exact(matrix.size()) {
        let coordinates = row
            .iter()
            .enumerate()
            .filter_map(|(index, &value)| (value != 0).then_some((index, value as u32)))
            .collect();
        basis.insert(model, prime, coordinates)?;
    }
    Ok(basis.rank())
}
