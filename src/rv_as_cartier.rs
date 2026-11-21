use crate::binary_field::BinaryField;
use crate::binary_linear::{self, Entry};
use crate::error::{Error, Result};
use crate::linear::{self, Basis};
use crate::model::{Model, Value};
use crate::rv_as_cover::Cover;
use crate::rv_ideal::check_storage;
use crate::rv_poly::Polynomial;

pub(crate) fn omega(model: &mut Model, cover: &Cover, index: usize) -> Result<Value> {
    if index >= cover.genus() {
        return Err(Error::plain(format!(
            "rvASOmega index {index} is out of range for genus {}",
            cover.genus()
        )));
    }
    Ok(Value::ASDiff(binary_linear::persist(
        model,
        vec![(index, 1)],
    )?))
}

pub(crate) fn apply(
    model: &mut Model,
    cover: &Cover,
    field: BinaryField,
    input: &[Entry],
) -> Result<Value> {
    if input.is_empty() {
        return Ok(Value::ASDiff([].into()));
    }
    let denominator = cover.cartier_denominator(model)?;
    let mut result = Vec::new();
    for &(index, coefficient) in input {
        let vector: Vec<_> = column(model, &denominator, index, cover.genus())?
            .into_iter()
            .map(|coordinate| (coordinate, 1))
            .collect();
        let root = field.frobenius(model, coefficient, field.degree() - 1)?;
        result = binary_linear::add_multiple(model, field, &result, &vector, root)?;
    }
    Ok(Value::ASDiff(binary_linear::persist(model, result)?))
}

pub(crate) fn rank(model: &mut Model, cover: &Cover, power: usize) -> Result<usize> {
    let genus = cover.genus();
    let power = power.min(genus);
    if power == 0 {
        return Ok(genus);
    }
    let denominator = cover.cartier_denominator(model)?;
    if let [monomial] = denominator.as_ref() {
        return Ok(monomial_rank(genus, monomial[0], power));
    }

    model.spend(genus)?;
    let mut matrix_cells = genus;
    let mut previous_cells = genus.saturating_mul(2);
    check_storage(matrix_cells.saturating_add(previous_cells))?;
    let mut matrix = Vec::with_capacity(genus);
    for index in 0..genus {
        let coefficients = column(model, &denominator, index, genus)?;
        matrix_cells += coefficients.len();
        check_storage(matrix_cells + previous_cells)?;
        matrix.push(coefficients);
    }
    let mut previous: Vec<_> = (0..genus).map(|index| vec![index]).collect();
    for _ in 0..power {
        let old_rank = previous.len();
        let mut span = Basis::default();
        let mut next = Vec::new();
        let mut next_cells = 0;
        for vector in &previous {
            let mut image = Vec::new();
            for &index in vector {
                image = linear::xor(model, &image, &matrix[index])?;
                check_storage(
                    matrix_cells + previous_cells + next_cells + span.cells() + image.len(),
                )?;
            }
            if span.insert(model, image.clone())? {
                next_cells += 1 + image.len();
                check_storage(matrix_cells + previous_cells + next_cells + span.cells())?;
                next.push(image);
            }
        }
        // The image subspaces are nested. Equality of dimensions makes all
        // following ranks equal, including for an arbitrarily large power.
        if next.len() == old_rank || next.is_empty() {
            return Ok(next.len());
        }
        previous = next;
        previous_cells = next_cells;
    }
    Ok(previous.len())
}

fn column(
    model: &mut Model,
    denominator: &Polynomial,
    index: usize,
    genus: usize,
) -> Result<Vec<usize>> {
    model.spend(1 + denominator.len())?;
    let mut result = Vec::new();
    // C(x^j dx/Q) = C(x^j Q dx)/Q. Only odd powers survive Cartier.
    for term in denominator.iter() {
        let exponent = term[0] as u128 + index as u128;
        if exponent & 1 == 0 {
            continue;
        }
        let target = (exponent / 2) as usize;
        if target >= genus {
            return Err(Error::plain(
                "Cartier image exceeds the holomorphic differential rank",
            ));
        }
        result.push(target);
    }
    Ok(result)
}

fn monomial_rank(genus: usize, exponent: usize, power: usize) -> usize {
    // For Q=x^e, a surviving column is j -> (j+e-1)/2. After k steps
    // its original index must be e-1 modulo 2^k; those images are distinct.
    if power >= usize::BITS as usize {
        return usize::from(exponent > 0 && exponent - 1 < genus);
    }
    let modulus = 1usize << power;
    let first = exponent.wrapping_sub(1) & (modulus - 1);
    if first >= genus {
        0
    } else {
        1 + (genus - 1 - first) / modulus
    }
}
