use super::{
    apply_matrix, arithmetic, divide_matrix, gcd, inconsistent, multiply, transpose, Matrix,
    Quotient, Vector, Wide,
};
use crate::error::{Error, Result};
use crate::integer_smith::{self, Track};
use crate::model::{Model, Value};
use std::rc::Rc;

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Endomorphism(pub(crate) Matrix);

impl Quotient {
    fn centralizer(&self, model: &mut Model) -> Result<&[Endomorphism; 8]> {
        if let Some(known) = self.centralizer.get() {
            return Ok(known);
        }
        let a = self.omega_numerator;
        let mut equations = integer_smith::Matrix::new(model, 16, 16)?;
        for (i, row) in a.iter().enumerate() {
            for j in 0..4 {
                for (k, &entry) in row.iter().enumerate() {
                    for (l, coefficients) in a.iter().enumerate() {
                        let first = if i == k { coefficients[j] } else { 0 };
                        let second = if l == j { entry } else { 0 };
                        let value = first.checked_sub(second).ok_or_else(arithmetic)?;
                        equations.write(model, 4 * i + j, 4 * k + l, value)?;
                    }
                }
            }
        }
        let reduction = integer_smith::reduce(model, equations, Track::None, Track::Basis)?;
        if reduction.factors.len() != 8 {
            return Err(inconsistent());
        }
        let mut integral: [[i128; 16]; 8] = std::array::from_fn(|i| {
            let slot = reduction.columns.basis(8 + i);
            std::array::from_fn(|j| slot.get(&j).copied().unwrap_or(0))
        });
        hermite_basis(model, &mut integral)?;
        let mut basis = [Endomorphism([[0; 4]; 4]); 8];
        for (vector, entries) in basis.iter_mut().zip(integral) {
            let matrix = std::array::from_fn(|i| std::array::from_fn(|j| entries[4 * i + j]));
            *vector = self.endomorphism(matrix)?;
        }
        model.allocate(128)?;
        Ok(self.centralizer.get_or_init(|| basis))
    }

    pub(crate) fn order_report(&self, model: &mut Model) -> Result<String> {
        let basis = self.centralizer(model)?.map(|endomorphism| endomorphism.0);
        Ok(format!(
            "rvCMQuotientEndoData = rank=8, basis={basis:?}, center_conductor={}, geometric_prime_to_characteristic_model={}",
            self.conductor,
            !(self.subgroup.members().len() as u64).is_multiple_of(self.kernel.characteristic),
        ))
    }

    pub(crate) fn endomorphism_basis(&self, model: &mut Model, i: usize) -> Result<Endomorphism> {
        self.centralizer(model)?.get(i).copied().ok_or_else(|| {
            Error::plain("CM quotient endomorphism basis index must be between zero and seven")
        })
    }

    pub(crate) fn endomorphism(&self, matrix: Matrix) -> Result<Endomorphism> {
        if multiply(matrix, self.omega_numerator)? != multiply(self.omega_numerator, matrix)? {
            return Err(Error::plain(
                "CM quotient endomorphism must commute with the diagonal CM field action",
            ));
        }
        Ok(Endomorphism(matrix))
    }

    pub(crate) fn compose_endomorphism(
        &self,
        a: Endomorphism,
        b: Endomorphism,
    ) -> Result<Endomorphism> {
        self.endomorphism(multiply(a.0, b.0)?)
    }

    fn rosati_pair(&self, map: Endomorphism) -> Result<(Matrix, i128)> {
        let alternating = self.form;
        let (adjugate, pfaffian) = alternating_inverse(alternating)?;
        let right = multiply(transpose(map.0), alternating)?;
        let numerator = multiply(adjugate, right)?;
        let sign = pfaffian.signum();
        let denominator = pfaffian
            .checked_abs()
            .filter(|&x| x != 0)
            .ok_or_else(arithmetic)?;
        let numerator = multiply(super::diagonal(sign), numerator)?;
        let content = numerator
            .iter()
            .flatten()
            .fold(denominator as u128, |c, a| gcd(c, a.unsigned_abs()))
            as i128;
        let numerator = divide_matrix(numerator, content).ok_or_else(inconsistent)?;
        let denominator = denominator / content;
        if multiply(alternating, numerator)? != multiply(super::diagonal(denominator), right)?
            || multiply(self.omega_numerator, numerator)?
                != multiply(numerator, self.omega_numerator)?
        {
            return Err(inconsistent());
        }
        Ok((numerator, denominator))
    }

    pub(crate) fn rosati_report(&self, map: Endomorphism) -> Result<String> {
        let (numerator, denominator) = self.rosati_pair(map)?;
        Ok(format!(
            "rvCMQuotientRosati = numerator={numerator:?}, denominator={denominator}, integral={}",
            denominator == 1
        ))
    }

    pub(crate) fn rosati_integral(&self, map: Endomorphism) -> Result<Endomorphism> {
        let (numerator, denominator) = self.rosati_pair(map)?;
        if denominator != 1 {
            return Err(Error::plain(
                "CM quotient Rosati involution can be rational at nonprincipal polarization",
            ));
        }
        self.endomorphism(numerator)
    }
}

// Row Hermite reduction keeps the saturated Smith kernel and reduces its
// coefficients before exposing the integral endomorphism basis.
pub(super) fn hermite_basis<const N: usize, const R: usize>(
    model: &mut Model,
    rows: &mut [[i128; N]; R],
) -> Result<()> {
    let mut rank = 0;
    for column in 0..N {
        let Some(pivot) = (rank..rows.len()).find(|&i| rows[i][column] != 0) else {
            continue;
        };
        rows.swap(rank, pivot);
        for lower in rank + 1..rows.len() {
            while rows[lower][column] != 0 {
                positive_pivot(&mut rows[rank], column)?;
                let pivot = rows[rank];
                reduce_row(model, &mut rows[lower], &pivot, column)?;
                if rows[lower][column] != 0 {
                    rows.swap(rank, lower);
                }
            }
        }
        positive_pivot(&mut rows[rank], column)?;
        let (above, current) = rows.split_at_mut(rank);
        for row in above {
            reduce_row(model, row, &current[0], column)?;
        }
        rank += 1;
    }
    if rank != rows.len() {
        return Err(inconsistent());
    }
    Ok(())
}

fn positive_pivot(row: &mut [i128], column: usize) -> Result<()> {
    if row[column] < 0 {
        for x in &mut row[column..] {
            *x = x.checked_neg().ok_or_else(arithmetic)?;
        }
    }
    Ok(())
}

fn reduce_row(model: &mut Model, row: &mut [i128], pivot: &[i128], column: usize) -> Result<()> {
    model.spend(row.len() - column)?;
    let q = row[column].div_euclid(pivot[column]);
    if q != 0 {
        for (x, &y) in row[column..].iter_mut().zip(&pivot[column..]) {
            *x = Wide::<8>::from_i128(*x)
                .add(Wide::product(q, y).negate())?
                .signed()
                .ok_or_else(arithmetic)?;
        }
    }
    Ok(())
}

pub(super) fn alternating_inverse(alternating: Matrix) -> Result<(Matrix, i128)> {
    let [_, a, b, c] = alternating[0];
    let [_, _, d, e] = alternating[1];
    let f = alternating[2][3];
    let pfaffian = Wide::<8>::product(a, f)
        .add(Wide::product(b, e).negate())?
        .add(Wide::product(c, d))?
        .signed()
        .ok_or_else(arithmetic)?;
    let opposite = |x: i128| x.checked_neg().ok_or_else(arithmetic);
    let adjugate = [
        [0, opposite(f)?, e, opposite(d)?],
        [f, 0, opposite(c)?, b],
        [opposite(e)?, c, 0, opposite(a)?],
        [d, opposite(b)?, a, 0],
    ];
    Ok((adjugate, pfaffian))
}

impl Endomorphism {
    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(32)?;
        Ok(Value::CMQuotientEnd(Rc::new(self)))
    }

    pub(crate) const fn identity() -> Self {
        Self([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    }

    pub(crate) fn apply(self, x: Vector) -> Result<Vector> {
        apply_matrix(self.0, x)
    }
}
