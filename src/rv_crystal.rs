use std::cell::OnceCell;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::Coordinate;
use crate::{eo, rv_flag, rv_newton};

pub(crate) type Entries = Vec<(usize, usize, i128)>;
type Column = Rc<[Coordinate]>;

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Characteristic,
    Newton,
    Matrices(usize),
}

pub(crate) struct Crystal {
    prime: Prime,
    genus: usize,
    frobenius: Box<[i128]>,
    verschiebung: Box<[i128]>,
    characteristic: OnceCell<Rc<[i128]>>,
}

impl Crystal {
    pub(crate) fn build(
        model: &mut Model,
        p: usize,
        genus: usize,
        entries: Entries,
    ) -> Result<Self> {
        if p < 2 || (p != 2 && p.is_multiple_of(2)) {
            return Err(Error::plain("F-crystal characteristic must be prime"));
        }
        let prime = if p == 2 {
            Prime::validated(2)
        } else {
            Prime::new(model, p)?
        };
        let n = genus
            .checked_mul(2)
            .ok_or_else(|| Error::plain("F-crystal dimension is out of range"))?;
        let slots = n.saturating_mul(n);
        check_storage(slots.saturating_mul(5))?;
        model.spend(slots)?;
        let mut f = vec![0; slots];
        let mut seen = vec![false; slots];
        for (row, column, value) in entries {
            model.spend(1)?;
            if row >= n || column >= n {
                return Err(Error::plain(format!(
                    "F-crystal matrix entry ({row},{column}) is outside dimension {n}"
                )));
            }
            let index = row * n + column;
            if seen[index] {
                return Err(Error::plain(format!(
                    "duplicate F-crystal matrix entry ({row},{column})"
                )));
            }
            seen[index] = true;
            f[index] = value;
        }
        drop(seen);

        // In the basis e_1,...,e_g,f_1,...,f_g, V is the J-adjoint F.
        model.spend(slots)?;
        let mut v = Vec::with_capacity(slots);
        for row in 0..n {
            for column in 0..n {
                let value = f[((column + genus) % n) * n + (row + genus) % n];
                v.push(if (row < genus) == (column < genus) {
                    value
                } else {
                    value.checked_neg().ok_or_else(arithmetic_error)?
                });
            }
        }
        // Testing VF=pI is exactly F^t J F=pJ. Over Q_p the resulting
        // inverse F^-1=p^-1 V is two-sided, so FV=pI as well.
        for row in 0..n {
            for column in 0..n {
                let row_entry = v[row * n..(row + 1) * n].iter().copied();
                let col_entry = f[column..].iter().step_by(n).copied();
                let value = pairing(model, row_entry, col_entry)?;
                let expected = if row == column {
                    i128::from(prime.get())
                } else {
                    0
                };
                if value != expected {
                    return Err(Error::plain(format!(
                        "F-crystal Frobenius fails F^t J F = pJ at adjoint product entry ({row},{column})"
                    )));
                }
            }
        }
        model.allocate(8 + 4 * slots)?;
        Ok(Self {
            prime,
            genus,
            frobenius: f.into_boxed_slice(),
            verschiebung: v.into_boxed_slice(),
            characteristic: OnceCell::new(),
        })
    }

    pub(crate) fn genus(&self) -> usize {
        self.genus
    }

    fn size(&self) -> usize {
        self.genus * 2
    }

    fn characteristic(&self, model: &mut Model) -> Result<Rc<[i128]>> {
        if let Some(result) = self.characteristic.get() {
            return Ok(result.clone());
        }
        let n = self.size();
        let slots = self.frobenius.len();
        check_storage(8 * slots + 2 * (n + 1))?;
        // Faddeev–LeVerrier applies Newton's identities to the whole
        // Frobenius matrix. Each division of a trace by k is exact over Z.
        // Working on the whole matrix also avoids the growth of unrelated
        // principal-block eigenvalues after a change of symplectic basis.
        let mut coefficients = vec![1];
        let mut current = vec![0; slots];
        for i in 0..n {
            current[i * n + i] = 1;
        }
        for k in 1..=n {
            let mut product = Vec::with_capacity(slots);
            let mut trace = 0i128;
            for row in 0..n {
                let source = &self.frobenius[row * n..(row + 1) * n];
                for column in 0..n {
                    let value = pairing(
                        model,
                        source.iter().copied(),
                        current[column..].iter().step_by(n).copied(),
                    )?;
                    if row == column {
                        trace = trace.checked_add(value).ok_or_else(arithmetic_error)?;
                    }
                    product.push(value);
                }
            }
            let k = k as i128;
            if trace % k != 0 {
                return Err(Error::plain(
                    "F-crystal characteristic trace division is not integral",
                ));
            }
            let coefficient = (trace / k).checked_neg().ok_or_else(arithmetic_error)?;
            coefficients.push(coefficient);
            for i in 0..n {
                let diagonal = &mut product[i * n + i];
                *diagonal = diagonal
                    .checked_add(coefficient)
                    .ok_or_else(arithmetic_error)?;
            }
            current = product;
        }
        model.allocate(1 + 2 * coefficients.len())?;
        Ok(self
            .characteristic
            .get_or_init(|| coefficients.into())
            .clone())
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, query: Query) -> Result<String> {
        match query {
            Query::Characteristic => {
                let characteristic = self.characteristic(model)?;
                model.spend(characteristic.len())?;
                Ok(format!("rvCrystalChar[{name}] = {characteristic:?}"))
            }
            Query::Newton => {
                let characteristic = self.characteristic(model)?;
                let mut valuations = Vec::with_capacity(characteristic.len());
                for &value in characteristic.iter() {
                    valuations.push(if value == 0 {
                        None
                    } else {
                        Some(valuation(model, self.prime, value.unsigned_abs())?)
                    });
                }
                let slopes = rv_newton::scaled_slopes(model, valuations, 1)?;
                Ok(format!("rvCrystalNewton[{name}] = [{}]", slopes.join(", ")))
            }
            Query::Matrices(level) => {
                validate_level(level)?;
                let modulus = u32::try_from(level)
                    .ok()
                    .and_then(|level| i128::from(self.prime.get()).checked_pow(level))
                    .ok_or_else(|| Error::plain("F-crystal matrix modulus is out of range"))?;
                let n = self.size();
                model.spend(2 * self.frobenius.len())?;
                check_storage(4 * self.frobenius.len() + 6 * n)?;
                let f = reduced_rows(&self.frobenius, n, modulus);
                let v = reduced_rows(&self.verschiebung, n, modulus);
                Ok(format!(
                    "rvCrystalMatrices[{name}, {level}] = F {f:?}, V {v:?}"
                ))
            }
        }
    }

    pub(crate) fn congruent(&self, model: &mut Model, other: &Self, level: usize) -> Result<bool> {
        validate_level(level)?;
        if (self.prime, self.genus) != (other.prime, other.genus) {
            return Err(Error::plain(
                "F-crystal congruence requires the same prime and genus",
            ));
        }
        let a = self.frobenius.iter().chain(self.verschiebung.iter());
        let b = other.frobenius.iter().chain(other.verschiebung.iter());
        for (&a, &b) in a.zip(b) {
            model.spend(1)?;
            if a != b && (valuation(model, self.prime, a.abs_diff(b))? as usize) < level {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub(crate) fn eo_type(&self, model: &mut Model) -> Result<Value> {
        let (mut f, mut v) = (Vec::new(), Vec::new());
        let mut cells = 0;
        for i in 0..self.size() {
            f.push(mod_column(
                model,
                self.prime,
                self.size(),
                &self.frobenius,
                i,
                &mut cells,
            )?);
            v.push(mod_column(
                model,
                self.prime,
                self.size(),
                &self.verschiebung,
                i,
                &mut cells,
            )?);
        }
        let sequence = rv_flag::final_type(
            model,
            self.prime,
            self.size(),
            || cells,
            |_, i| Ok(f[i].clone()),
            |_, i| Ok(v[i].clone()),
        )?;
        eo::persist(model, sequence)
    }
}

fn pairing(
    model: &mut Model,
    a: impl ExactSizeIterator<Item = i128>,
    b: impl Iterator<Item = i128>,
) -> Result<i128> {
    model.spend(a.len())?;
    a.zip(b).try_fold(0i128, |sum, (c, d)| {
        let term = c.checked_mul(d).ok_or_else(arithmetic_error)?;
        sum.checked_add(term).ok_or_else(arithmetic_error)
    })
}

fn reduced_rows(matrix: &[i128], n: usize, modulus: i128) -> Vec<Vec<i128>> {
    matrix
        .chunks(n.max(1))
        .map(|row| row.iter().map(|a| a.rem_euclid(modulus)).collect())
        .collect()
}

fn mod_column(
    model: &mut Model,
    prime: Prime,
    n: usize,
    matrix: &[i128],
    index: usize,
    cells: &mut usize,
) -> Result<Column> {
    model.spend(n)?;
    let p = i128::from(prime.get());
    let terms = matrix[index..]
        .iter()
        .step_by(n)
        .enumerate()
        .filter_map(|(row, &entry)| {
            let value = entry.rem_euclid(p) as u32;
            (value != 0).then_some((row, value))
        })
        .collect::<Vec<_>>();
    *cells += 3 + 2 * terms.len();
    check_storage(*cells)?;
    Ok(terms.into())
}

fn valuation(model: &mut Model, prime: Prime, mut nonzero: u128) -> Result<i64> {
    let p = u128::from(prime.get());
    let mut v = 0;
    while nonzero.is_multiple_of(p) {
        model.spend(1)?;
        nonzero /= p;
        v += 1;
    }
    Ok(v)
}

fn validate_level(level: usize) -> Result<()> {
    if level == 0 {
        return Err(Error::plain("F-crystal lattice level must be positive"));
    }
    Ok(())
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain("F-crystal matrices exceed the storage limit"));
    }
    Ok(())
}

fn arithmetic_error() -> Error {
    Error::plain("F-crystal integer arithmetic is out of range")
}
