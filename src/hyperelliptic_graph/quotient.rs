use crate::error::{Error, Result};
use crate::integer_smith::{self, Vector};
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::wide_integer::Natural;

pub(super) fn factors(model: &mut Model, periods: &[Vec<i128>]) -> Result<Vec<i128>> {
    let order = determinant::<8>(model, periods, 127)?
        .ok_or_else(|| Error::plain("hyperelliptic component-group order exceeds signed i128"))?;
    Ok(presentation(model, periods, None, order)?.factors)
}

pub(super) struct Quotient {
    pub(super) factors: Vec<i128>,
    pub(super) characters: Vec<Vec<i128>>,
}

pub(super) fn characters(
    model: &mut Model,
    periods: &[Vec<i128>],
    action: &[Vector],
) -> Result<Quotient> {
    if let Some(order) = determinant::<8>(model, periods, 127)? {
        return presentation(model, periods, Some(action), order);
    }
    // The fixed subgroup can fit even when the geometric order does not.
    if let Some(order) = determinant::<16>(model, periods, 255)? {
        return presentation(model, periods, Some(action), order);
    }
    let order = determinant::<64>(model, periods, 1023)?
        .ok_or_else(|| Error::plain("hyperelliptic period determinant exceeds the work range"))?;
    presentation(model, periods, Some(action), order)
}

fn determinant<const N: usize>(
    model: &mut Model,
    periods: &[Vec<i128>],
    bits: usize,
) -> Result<Option<Natural<N>>> {
    let one = Wide::from_i128(1);
    let mut previous = one;
    let rank = periods.len();
    integer_smith::storage(rank * (3 + (N + 1) * rank))?;
    model.spend(rank * rank)?;
    let mut matrix: Vec<Vec<_>> = periods
        .iter()
        .map(|row| row.iter().copied().map(Wide::from_i128).collect())
        .collect();
    for i in 0..rank {
        let pivot = matrix[i][i];
        if pivot.negative || pivot.is_zero() {
            return Err(Error::plain(
                "integral hyperelliptic monodromy pairing is degenerate",
            ));
        }
        if pivot.magnitude.bits() > bits {
            return Ok(None);
        }
        for j in i + 1..rank {
            for k in j..rank {
                model.spend(1 + limb_work(pivot.magnitude, matrix[j][k].magnitude, N))?;
                let diagonal = pivot.multiply_wide(matrix[j][k])?;
                model.spend(limb_work(matrix[j][i].magnitude, matrix[i][k].magnitude, N))?;
                let corner = matrix[j][i].multiply_wide(matrix[i][k])?;
                let value = diagonal
                    .add(corner.negate())?
                    .divide_exact(model, previous)?
                    .ok_or_else(|| Error::plain("nonintegral hyperelliptic Bareiss minor"))?;
                // Q >= I in the fundamental cycle basis. Every principal
                // minor is at most det Q, and the corresponding mixed
                // minors are bounded by Cauchy-Schwarz on exterior powers.
                if value.magnitude.bits() > bits {
                    return Ok(None);
                }
                matrix[j][k] = value;
                matrix[k][j] = value;
            }
        }
        previous = pivot;
    }
    Ok(Some(previous.magnitude))
}

fn limb_work<const N: usize>(a: Natural<N>, b: Natural<N>, width: usize) -> usize {
    width + a.bits().div_ceil(32) * b.bits().div_ceil(32)
}

// Both quotients are killed by det Q. Do the ordinary integer pivot
// reduction with det Q times the active coordinate vectors implicit.
// Reduction modulo det Q therefore preserves the presented lattice.
fn presentation<const N: usize>(
    model: &mut Model,
    periods: &[Vec<i128>],
    action: Option<&[Vector]>,
    modulus: Natural<N>,
) -> Result<Quotient> {
    let rank = periods.len();
    let one = Natural::from_small(1);
    if modulus == one {
        return Ok(Quotient {
            factors: Vec::new(),
            characters: Vec::new(),
        });
    }
    let input_width = rank * (1 + usize::from(action.is_some()));
    let frame_width = if action.is_some() { rank } else { 0 };
    let width = input_width + rank;
    integer_smith::storage((width + frame_width) * (3 + N * rank) + 2 * rank * frame_width)?;
    model.spend(rank * width)?;
    let mut smith = FiniteSmith {
        modulus,
        columns: vec![vec![Natural::default(); rank]; width],
        frame: action.map(|_| {
            (0..rank)
                .map(|i| {
                    let mut row = vec![Natural::default(); rank];
                    row[i] = one;
                    row
                })
                .collect()
        }),
    };
    for (i, row) in periods.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            smith.columns[j][i] = residue(model, value, modulus)?;
        }
    }
    if let Some(action) = action {
        for (i, column) in action.iter().enumerate() {
            for j in 0..rank {
                let value = i128::from(i == j) - column.get(&j).copied().unwrap_or(0);
                smith.columns[rank + j][i] = residue(model, value, modulus)?;
            }
        }
    }
    let mut output = Quotient {
        factors: Vec::new(),
        characters: Vec::new(),
    };
    let mut remaining = modulus;
    let mut previous = one;
    for i in 0..rank {
        // Add this one exact column before reducing row i. The remaining
        // modulus columns are still implicit in the trailing sublattice.
        smith.columns[input_width + i][i] = modulus;
        let pivot = smith.pivot(model, i)?;
        if pivot != one {
            let (quotient, remainder) = remaining.div_rem(model, pivot)?;
            if !remainder.is_zero() || !pivot.div_rem(model, previous)?.1.is_zero() {
                return Err(Error::plain("invalid hyperelliptic Smith pivot"));
            }
            remaining = quotient;
            previous = pivot;
            let d = pivot.signed(false).ok_or_else(|| {
                Error::plain("hyperelliptic component-group factor exceeds signed i128")
            })?;
            output.factors.push(d);
            if let Some(frame) = &smith.frame {
                let character = frame[i]
                    .iter()
                    .map(|&value| {
                        let (_, residue) = value.divide(model, d as u128)?;
                        Ok(residue as i128)
                    })
                    .collect::<Result<_>>()?;
                output.characters.push(character);
            }
        }
    }
    if action.is_none() && remaining != one {
        return Err(Error::plain(
            "hyperelliptic Smith order differs from the period determinant",
        ));
    }
    Ok(output)
}

struct FiniteSmith<const N: usize> {
    modulus: Natural<N>,
    columns: Vec<Vec<Natural<N>>>,
    frame: Option<Vec<Vec<Natural<N>>>>,
}

impl<const N: usize> FiniteSmith<N> {
    fn pivot(&mut self, model: &mut Model, i: usize) -> Result<Natural<N>> {
        let rank = self.columns[i].len();
        loop {
            self.clear_row(model, i)?;
            let pivot = self.columns[i][i];
            let mut remainder = false;
            for j in i + 1..rank {
                let entry = self.columns[i][j];
                if entry.is_zero() {
                    continue;
                }
                let (q, r) = entry.div_rem(model, pivot)?;
                self.columns[i][j] = r;
                // All other active entries of row i are already zero.
                if let Some(frame) = &mut self.frame {
                    let (head, tail) = frame.split_at_mut(j);
                    subtract_multiple(model, &mut tail[0], &head[i], q, self.modulus)?;
                }
                if !r.is_zero() {
                    model.spend(self.columns.len())?;
                    for column in &mut self.columns {
                        column.swap(i, j);
                    }
                    if let Some(frame) = &mut self.frame {
                        frame.swap(i, j);
                    }
                    remainder = true;
                    break;
                }
            }
            if remainder {
                continue;
            }
            let Some(row) = self.indivisible_row(model, i, pivot)? else {
                return Ok(pivot);
            };
            model.spend(self.columns.len() - i - 1)?;
            for column in &mut self.columns[i + 1..] {
                column[i] = add_mod(column[i], column[row], self.modulus);
            }
            if let Some(frame) = &mut self.frame {
                let (head, tail) = frame.split_at_mut(row);
                model.spend(rank)?;
                for (entry, &other) in head[i].iter_mut().zip(&tail[0]) {
                    *entry = add_mod(*entry, other, self.modulus);
                }
            }
        }
    }

    fn clear_row(&mut self, model: &mut Model, i: usize) -> Result<()> {
        if self.columns[i][i] != Natural::from_small(1) {
            let mut best = i;
            for j in i + 1..self.columns.len() {
                model.spend(1)?;
                let candidate = self.columns[j][i];
                let pivot = self.columns[best][i];
                if !candidate.is_zero() && (pivot.is_zero() || candidate.cmp(pivot).is_lt()) {
                    best = j;
                    if candidate == Natural::from_small(1) {
                        break;
                    }
                }
            }
            self.columns.swap(i, best);
        }
        for j in i + 1..self.columns.len() {
            model.spend(1)?;
            while !self.columns[j][i].is_zero() {
                let a = self.columns[i][i];
                if a.is_zero() {
                    self.columns.swap(i, j);
                    break;
                }
                let (q, r) = self.columns[j][i].div_rem(model, a)?;
                self.columns[j][i] = r;
                let (head, tail) = self.columns.split_at_mut(j);
                subtract_multiple(
                    model,
                    &mut tail[0][i + 1..],
                    &head[i][i + 1..],
                    q,
                    self.modulus,
                )?;
                if !r.is_zero() {
                    self.columns.swap(i, j);
                }
            }
        }
        Ok(())
    }

    fn indivisible_row(
        &self,
        model: &mut Model,
        i: usize,
        pivot: Natural<N>,
    ) -> Result<Option<usize>> {
        if pivot == Natural::from_small(1) {
            return Ok(None);
        }
        for column in &self.columns[i + 1..] {
            for (row, &value) in column.iter().enumerate().skip(i + 1) {
                model.spend(1)?;
                if !value.is_zero() && !value.div_rem(model, pivot)?.1.is_zero() {
                    return Ok(Some(row));
                }
            }
        }
        Ok(None)
    }
}

fn residue<const N: usize>(
    model: &mut Model,
    value: i128,
    modulus: Natural<N>,
) -> Result<Natural<N>> {
    let abs = Natural::from_small(value.unsigned_abs());
    let r = abs.div_rem(model, modulus)?.1;
    Ok(if value < 0 && !r.is_zero() {
        modulus.difference(r)
    } else {
        r
    })
}

fn add_mod<const N: usize>(mut a: Natural<N>, b: Natural<N>, modulus: Natural<N>) -> Natural<N> {
    a.add_multiple(b, 1)
        .expect("twice the modulus fits the scratch width");
    if a.cmp(modulus).is_lt() {
        a
    } else {
        a.difference(modulus)
    }
}

fn subtract_multiple<const N: usize>(
    model: &mut Model,
    target: &mut [Natural<N>],
    source: &[Natural<N>],
    factor: Natural<N>,
    modulus: Natural<N>,
) -> Result<()> {
    if factor.is_zero() {
        return Ok(());
    }
    model.spend(source.len())?;
    for (entry, &value) in target.iter_mut().zip(source) {
        if value.is_zero() {
            continue;
        }
        let change = if factor == Natural::from_small(1) {
            value
        } else {
            model.spend(limb_work(value, factor, N))?;
            let product = value
                .checked_product(factor)
                .expect("the product of two residues fits the scratch width");
            product.div_rem(model, modulus)?.1
        };
        *entry = if entry.cmp(change).is_lt() {
            modulus.difference(change.difference(*entry))
        } else {
            entry.difference(change)
        };
    }
    Ok(())
}
