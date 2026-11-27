use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear;

use super::linear::{add, mul, negate};
use super::{certificate, polynomial, Block, FactorKind, Order};

impl Order {
    pub(super) fn congruence(
        &self,
        model: &mut Model,
        scalar: usize,
        algebraic: usize,
        ell: usize,
    ) -> Result<String> {
        let prime = if ell == 2 {
            Prime::validated(2)
        } else {
            Prime::new(model, ell)?
        };
        let line = self.block(scalar)?;
        let field = self.block(algebraic)?;
        if line.factor.len() != 2
            || line.kernel.len() != 1
            || field.factor.len() <= 2
            || field.kind != FactorKind::Irreducible
        {
            return Err(Error::plain(
                "Brandt congruence requires a rational eigenline and a certified higher-degree eigenfield",
            ));
        }
        let Some(images) = &field.images else {
            return Err(Error::plain(
                "Brandt congruence requires a simple eigenfield block",
            ));
        };
        let root = negate(line.factor[0])?;
        let scalar_values: Vec<_> = line.restrictions.iter().map(|a| a[0][0]).collect();
        let p = i128::from(prime.get());
        let mut valuations = vec![(
            "minimal polynomial",
            valuation(polynomial::evaluate(&field.factor, root)?, p),
        )];
        for ((&name, image), &value) in self.operator_names.iter().zip(images).zip(&scalar_values) {
            if image.denominator % p == 0 {
                return Err(Error::plain(format!(
                    "{name} has a nonunit denominator at {ell}; congruence in this polynomial presentation is unsupported"
                )));
            }
            let residual = add(
                polynomial::evaluate(&image.numerator, root)?,
                negate(mul(image.denominator, value)?)?,
            )?;
            valuations.push((name, valuation(residual, p)));
        }
        let depth = valuations
            .iter()
            .filter_map(|(_, v)| *v)
            .min()
            .ok_or_else(|| certificate("distinct characteristic-zero eigencharacters"))?;
        let valuations = valuations
            .iter()
            .map(|(name, v)| {
                let order = v.map_or_else(|| "∞".to_string(), |value| value.to_string());
                format!("{name}={order}")
            })
            .collect::<Vec<_>>()
            .join(", ");
        let summary = format!(
            "prime {ell}; T eigenvalue = {root}; residue root = {}; exact joint congruence depth = {depth}; residual valuations: {valuations}",
            root.rem_euclid(p),
        );
        if depth == 0 {
            return Ok(summary);
        }
        let scalar_reduction = eigenvector(model, line, prime, &scalar_values)?;
        let algebraic_reduction = eigenvector(model, field, prime, &scalar_values)?;
        Ok(format!(
            "{summary}\nscalar block ambient eigenvector = {scalar_reduction:?}\n\
             algebraic block ambient eigenvector = {algebraic_reduction:?}"
        ))
    }
}

fn valuation(mut value: i128, p: i128) -> Option<usize> {
    if value == 0 {
        return None;
    }
    let mut order = 0;
    while value % p == 0 {
        order += 1;
        value /= p;
    }
    Some(order)
}

fn eigenvector(
    model: &mut Model,
    block: &Block,
    prime: Prime,
    eigenvalues: &[i128],
) -> Result<Vec<u32>> {
    let p = i128::from(prime.get());
    let d = block.kernel.len();
    let residue = |value: i128| value.rem_euclid(p) as u32;
    let mut equations = Vec::with_capacity(d * eigenvalues.len());
    for (matrix, &scalar) in block.restrictions.iter().zip(eigenvalues) {
        for (i, row) in matrix.iter().enumerate() {
            let mut coefficients: Vec<_> = row.iter().copied().map(residue).collect();
            coefficients[i] = prime.add(coefficients[i], prime.negate(residue(scalar)));
            equations.push(coefficients);
        }
    }
    let original = equations.clone();
    model.spend(equations.len() * d * d)?;
    let mut pivots = Vec::new();
    let mut free = None;
    for j in 0..d {
        let i = pivots.len();
        let Some(pivot) = (i..equations.len()).find(|&row| equations[row][j] != 0) else {
            free.get_or_insert(j);
            continue;
        };
        equations.swap(i, pivot);
        let inverse = prime_linear::invert(model, prime, equations[i][j])?;
        for a in &mut equations[i][j..] {
            *a = prime.mul(*a, inverse);
        }
        let pivot_row = equations[i].clone();
        for (k, row) in equations.iter_mut().enumerate() {
            if k != i {
                let coefficient = prime.negate(row[j]);
                for (a, &b) in row[j..].iter_mut().zip(&pivot_row[j..]) {
                    *a = prime.add(*a, prime.mul(coefficient, b));
                }
            }
        }
        pivots.push(j);
    }
    let free = free.ok_or_else(|| certificate("nonzero common residue eigenvector"))?;
    let mut solution = vec![0; d];
    solution[free] = 1;
    for (row, j) in pivots.into_iter().enumerate() {
        solution[j] = prime.negate(equations[row][free]);
    }
    if original.iter().any(|row| {
        row.iter()
            .zip(&solution)
            .fold(0, |sum, (&a, &b)| prime.add(sum, prime.mul(a, b)))
            != 0
    }) {
        return Err(certificate("modular simultaneous eigenvector equations"));
    }
    let mut ambient = vec![0; block.kernel[0].len()];
    for (vector, &scale) in block.kernel.iter().zip(&solution) {
        for (a, &b) in ambient.iter_mut().zip(vector) {
            *a = prime.add(*a, prime.mul(scale, residue(b)));
        }
    }
    let first = ambient
        .iter()
        .find(|&&a| a != 0)
        .ok_or_else(|| certificate("saturated eigenlattice injects modulo ell"))?;
    let scale = prime_linear::invert(model, prime, *first)?;
    for a in &mut ambient {
        *a = prime.mul(*a, scale);
    }
    Ok(ambient)
}
