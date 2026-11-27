use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::rv_lens::{self, Rational as Scalar};

use super::super::{Morphism, Write};
use super::field::Field;
use super::{check_report, Defect, Error, Model, Object, Rc, Result, Value};

type Vector = Vec<Scalar>;

struct Rect {
    columns: usize,
    rows: Vec<Vector>,
}

impl Rect {
    fn zero(rows: usize, columns: usize) -> Self {
        Self {
            columns,
            rows: vec![vec![Scalar::zero(); columns]; rows],
        }
    }

    fn identity(n: usize) -> Self {
        let mut matrix = Self::zero(n, n);
        for (i, row) in matrix.rows.iter_mut().enumerate() {
            row[i] = Scalar::integer(1);
        }
        matrix
    }

    fn compose(&self, model: &mut Model, field: &Field, before: &Self) -> Result<Self> {
        let mut result = Self::zero(self.rows.len(), before.columns);
        for (row, out) in self.rows.iter().zip(&mut result.rows) {
            for (a, b) in row.iter().zip(&before.rows) {
                let minus = rv_lens::negative(model, a)?;
                subtract_multiple(model, field, out, &minus, b, 0)?;
            }
        }
        Ok(result)
    }

    fn apply(&self, model: &mut Model, field: &Field, v: &[Scalar]) -> Result<Vector> {
        let mut result = Vec::with_capacity(self.rows.len());
        for row in &self.rows {
            model.spend(self.columns)?;
            let mut coefficient = Scalar::zero();
            for (a, b) in row.iter().zip(v) {
                let term = field.multiply(model, a, b)?;
                coefficient = field.add(model, &coefficient, &term)?;
            }
            result.push(coefficient);
        }
        Ok(result)
    }

    fn echelon(&mut self, model: &mut Model, field: &Field) -> Result<Vec<usize>> {
        let mut columns = Vec::new();
        for col in 0..self.columns {
            let rank = columns.len();
            model.spend(1 + self.rows.len() - rank)?;
            let Some(pivot) = (rank..self.rows.len()).find(|&r| !self.rows[r][col].is_zero())
            else {
                continue;
            };
            self.rows.swap(rank, pivot);
            let inverse = field.inverse_unit(model, &self.rows[rank][col])?;
            for entry in &mut self.rows[rank][col..] {
                *entry = field.multiply(model, entry, &inverse)?;
            }
            let pivot = self.rows[rank].clone();
            for (i, row) in self.rows.iter_mut().enumerate() {
                if i != rank {
                    let factor = row[col].clone();
                    subtract_multiple(model, field, row, &factor, &pivot, col)?;
                }
            }
            columns.push(col);
        }
        Ok(columns)
    }

    fn kernel(&self, model: &mut Model, pivots: &[usize]) -> Result<Vec<Vector>> {
        let mut basis = Vec::with_capacity(self.columns - pivots.len());
        for j in 0..self.columns {
            if pivots.binary_search(&j).is_ok() {
                continue;
            }
            let mut vector = vec![Scalar::zero(); self.columns];
            vector[j] = Scalar::integer(1);
            for (i, &col) in pivots.iter().enumerate() {
                vector[col] = rv_lens::negative(model, &self.rows[i][j])?;
            }
            basis.push(vector);
        }
        Ok(basis)
    }
}

fn subtract_multiple(
    model: &mut Model,
    field: &Field,
    destination: &mut [Scalar],
    factor: &Scalar,
    source: &[Scalar],
    start: usize,
) -> Result<()> {
    if factor.is_zero() {
        return Ok(());
    }
    model.spend(destination.len() - start)?;
    for (a, b) in destination[start..].iter_mut().zip(&source[start..]) {
        let product = field.multiply(model, factor, b)?;
        *a = field.subtract(model, a, &product)?;
    }
    Ok(())
}

fn extend(
    model: &mut Model,
    field: &Field,
    basis: &mut [Option<Vector>],
    mut v: Vector,
) -> Result<bool> {
    for (col, pivot) in basis.iter_mut().enumerate() {
        if v[col].is_zero() {
            continue;
        }
        let factor = v[col].clone();
        if let Some(row) = pivot {
            subtract_multiple(model, field, &mut v, &factor, row, col)?;
        } else {
            let inverse = field.inverse_unit(model, &factor)?;
            for coefficient in &mut v[col..] {
                *coefficient = field.multiply(model, coefficient, &inverse)?;
            }
            *pivot = Some(v);
            return Ok(true);
        }
    }
    Ok(false)
}

#[derive(Clone, Copy, Default)]
struct Location {
    grade: usize,
    index: usize,
}

struct Chain {
    grade: usize,
    vectors: Vec<Vector>,
}

#[derive(Clone, Copy)]
struct Block {
    start: usize,
    length: usize,
}

pub(super) struct Summand {
    j: usize,
    monomials: Vec<Location>,
    jordan: Vec<Location>,
    blocks: Vec<Block>,
    ranks: Vec<usize>,
    x: Vec<Rect>,
    basis: Vec<Matrix<Scalar>>,
    inverse: Vec<Matrix<Scalar>>,
}

fn dimensions(model: &mut Model, field: &Field, r: usize, j: usize) -> Result<Vec<usize>> {
    let dim = r * j;
    model.spend(dim)?;
    let dimensions: Vec<_> = (0..r + j - 1)
        .map(|h| (r - 1).min(h) - h.saturating_sub(j - 1) + 1)
        .collect();
    let matrices: usize = dimensions.iter().map(|&n| n * n).sum();
    let powers: usize = dimensions
        .iter()
        .enumerate()
        .map(|(g, &n)| n * dimensions[g..].iter().take(field.d + 1).sum::<usize>())
        .sum();
    if (powers + 6 * matrices + 10 * dim).saturating_mul(field.cells())
        > crate::complex::MAX_PIVOT_CELLS
    {
        return Err(Error::plain(
            "LG defect Jordan calculation exceeds the coefficient storage limit",
        ));
    }
    Ok(dimensions)
}

impl Summand {
    pub(super) fn new(model: &mut Model, field: &Field, defect: &Defect, j: usize) -> Result<Self> {
        let r = defect.indices.len();
        let dims = dimensions(model, field, r, j)?;
        let mut monomials = Vec::with_capacity(r * j);
        let mut counts = vec![0; dims.len()];
        for a in 0..r {
            for b in 0..j {
                let grade = a + b;
                monomials.push(Location {
                    grade,
                    index: counts[grade],
                });
                counts[grade] += 1;
            }
        }
        let mut x: Vec<_> = dims
            .iter()
            .enumerate()
            .map(|(g, &n)| Rect::zero(dims.get(g + 1).copied().unwrap_or(0), n))
            .collect();
        let negative_p: Vec<_> = defect.p[..r]
            .iter()
            .map(|c| {
                let c = rv_lens::fraction(model, c, 1)?;
                rv_lens::negative(model, &c)
            })
            .collect::<Result<_>>()?;
        for a in 0..r {
            for b in 0..j {
                let col = monomials[a * j + b];
                if a + 1 < r {
                    let row = monomials[(a + 1) * j + b];
                    x[col.grade].rows[row.index][col.index] = Scalar::integer(1);
                } else {
                    for (c, a) in negative_p.iter().enumerate() {
                        let y = b + r - c;
                        if y < j {
                            let row = monomials[c * j + y];
                            x[col.grade].rows[row.index][col.index] = a.clone();
                        }
                    }
                }
            }
        }
        let (chains, ranks) = chains(model, field, &dims, &x)?;
        let mut jordan = Vec::with_capacity(r * j);
        let mut blocks = Vec::with_capacity(chains.len());
        let mut positions = vec![0; dims.len()];
        let mut matrices: Vec<_> = dims.iter().map(|&n| Rect::zero(n, n)).collect();
        for chain in chains {
            blocks.push(Block {
                start: jordan.len(),
                length: chain.vectors.len(),
            });
            for (k, vector) in chain.vectors.into_iter().enumerate() {
                let grade = chain.grade + k;
                let index = positions[grade];
                if index >= dims[grade] {
                    return Err(Error::plain(
                        "LG defect produced too many Jordan vectors in one degree",
                    ));
                }
                for (row, value) in matrices[grade].rows.iter_mut().zip(vector) {
                    row[index] = value;
                }
                jordan.push(Location { grade, index });
                positions[grade] += 1;
            }
        }
        if positions != dims {
            return Err(Error::plain(
                "LG defect Jordan vectors do not span the quotient",
            ));
        }
        let mut basis = Vec::with_capacity(dims.len());
        let mut inverse = Vec::with_capacity(dims.len());
        for (n, matrix) in dims.into_iter().zip(matrices) {
            let matrix = Matrix::from_entries(n, matrix.rows.into_iter().flatten().collect())?;
            let Some(back) = hw_algebra::inverse(model, field, &matrix)? else {
                return Err(Error::plain("LG defect Jordan basis is singular"));
            };
            let identity = hw_algebra::identity(model, field, n)?;
            if hw_algebra::multiply(model, field, &matrix, &back)? != identity
                || hw_algebra::multiply(model, field, &back, &matrix)? != identity
            {
                return Err(Error::plain(
                    "LG defect Jordan basis inverse failed verification",
                ));
            }
            basis.push(matrix);
            inverse.push(back);
        }
        Ok(Self {
            j,
            monomials,
            jordan,
            blocks,
            ranks,
            x,
            basis,
            inverse,
        })
    }

    pub(super) fn lengths(&self) -> impl Iterator<Item = usize> + '_ {
        self.blocks.iter().map(|block| block.length)
    }

    fn cells(&self, field: &Field) -> usize {
        let cells = self
            .x
            .iter()
            .map(|m| m.columns * m.rows.len())
            .sum::<usize>()
            + 2 * self.basis.iter().map(|m| m.entries().len()).sum::<usize>();
        cells * field.cells() + 10 * self.monomials.len() + self.ranks.len()
    }

    fn entry(&self, kind: usize, row: usize, col: usize) -> Result<Scalar> {
        if row >= self.monomials.len() || col >= self.monomials.len() {
            return Err(Error::plain(
                "LG defect action matrix index is out of range",
            ));
        }
        let (r, c) = match kind {
            0 => (self.monomials[row], self.monomials[col]),
            1 => (self.monomials[row], self.jordan[col]),
            2 => (self.jordan[row], self.monomials[col]),
            _ => {
                return Err(Error::plain(
                    "LG defect matrix kind must be 0 (x), 1 (basis), or 2 (inverse)",
                ))
            }
        };
        match kind {
            0 if r.grade == c.grade + 1 => Ok(self.x[c.grade].rows[r.index][c.index].clone()),
            1 if r.grade == c.grade => Ok(self.basis[r.grade].entry(r.index, c.index).clone()),
            2 if r.grade == c.grade => Ok(self.inverse[r.grade].entry(r.index, c.index).clone()),
            _ => Ok(Scalar::zero()),
        }
    }
}

fn chains(
    model: &mut Model,
    field: &Field,
    dims: &[usize],
    x: &[Rect],
) -> Result<(Vec<Chain>, Vec<usize>)> {
    let mut powers = Vec::with_capacity(dims.len());
    for (g, &n) in dims.iter().enumerate() {
        let mut row = vec![Rect::identity(n)];
        for map in x[g..].iter().take(field.d) {
            row.push(map.compose(model, field, row.last().unwrap())?);
        }
        powers.push(row);
    }
    let mut endpoints: Vec<_> = dims.iter().map(|&n| vec![None; n]).collect();
    let mut chains = Vec::new();
    let mut ranks = vec![0; field.d + 1];
    ranks[0] = dims.iter().sum();
    // X^(r-1)(ker X^r) is the endpoint filtration in ker X. Independent
    // endpoints of longer chains are kept when choosing each shorter chain.
    for (r, rank) in ranks.iter_mut().enumerate().skip(1).rev() {
        for (g, row) in powers.iter_mut().enumerate() {
            if r >= row.len() {
                continue;
            }
            let pivots = row[r].echelon(model, field)?;
            *rank += pivots.len();
            for head in row[r].kernel(model, &pivots)? {
                let endpoint = row[r - 1].apply(model, field, &head)?;
                if !extend(model, field, &mut endpoints[g + r - 1], endpoint)? {
                    continue;
                }
                let vectors: Vec<_> = row[..r]
                    .iter()
                    .map(|power| power.apply(model, field, &head))
                    .collect::<Result<_>>()?;
                for (k, v) in vectors.iter().enumerate() {
                    let expected = x[g + k].apply(model, field, v)?;
                    if let Some(next) = vectors.get(k + 1) {
                        if expected != *next {
                            return Err(Error::plain(
                                "LG defect Jordan chain has an incorrect successor",
                            ));
                        }
                    } else if expected.iter().any(|c| !c.is_zero()) {
                        return Err(Error::plain("LG defect Jordan chain does not terminate"));
                    }
                }
                chains.push(Chain { grade: g, vectors });
            }
        }
    }
    if ranks[field.d] != 0
        || ranks.iter().enumerate().any(|(k, &rank)| {
            rank != chains
                .iter()
                .map(|chain| chain.vectors.len().saturating_sub(k))
                .sum::<usize>()
        })
    {
        return Err(Error::plain(
            "LG defect Jordan lengths failed the exact rank tests",
        ));
    }
    chains.sort_by_key(|chain| chain.vectors.len());
    Ok((chains, ranks))
}

pub(crate) struct Push {
    pub(crate) defect: Rc<str>,
    pub(crate) source: Rc<str>,
    mark: Rc<str>,
    pub(crate) d: usize,
    summands: Vec<Summand>,
}

impl Push {
    pub(crate) fn new(
        model: &mut Model,
        defect_name: Rc<str>,
        source_name: Rc<str>,
        defect: &Defect,
        source: &Object,
    ) -> Result<Self> {
        defect.check_source(source)?;
        let field = Field::new(model, defect.d)?;
        let mut summands = Vec::with_capacity(source.factors.len());
        model.allocate(12 + source.factors.len())?;
        for &j in &source.factors {
            let summand = Summand::new(model, &field, defect, j)?;
            model.allocate(summand.cells(&field))?;
            summands.push(summand);
        }
        Ok(Self {
            defect: defect_name,
            source: source_name,
            mark: source.mark.clone(),
            d: defect.d,
            summands,
        })
    }

    pub(crate) fn report(&self) -> Result<String> {
        let mut report = String::new();
        for (i, summand) in self.summands.iter().enumerate() {
            let stable: Vec<_> = summand.lengths().filter(|&j| j < self.d).collect();
            let projective = summand.lengths().filter(|&j| j == self.d).count();
            write!(
                &mut report,
                "{}P{} ranks={:?} stable={stable:?} projective={projective}",
                if i == 0 { "[" } else { ", " },
                summand.j,
                summand.ranks
            )
            .unwrap();
            check_report(&report)?;
        }
        if report.is_empty() {
            report.push('[');
        }
        report.push(']');
        Ok(report)
    }

    pub(crate) fn result(&self) -> Object {
        let factors = self
            .summands
            .iter()
            .flat_map(Summand::lengths)
            .filter(|&j| j < self.d)
            .collect();
        Object {
            mark: self.mark.clone(),
            d: self.d,
            factors,
        }
    }

    pub(crate) fn entry(
        &self,
        model: &mut Model,
        summand: usize,
        kind: usize,
        row: usize,
        col: usize,
    ) -> Result<Value> {
        let summand = self
            .summands
            .get(summand)
            .ok_or_else(|| Error::plain("LG defect summand index is out of range"))?;
        rv_lens::rational_value(model, summand.entry(kind, row, col)?)
    }

    fn block(&self, mut i: usize) -> Result<(usize, &Summand, Block)> {
        for (s, summand) in self.summands.iter().enumerate() {
            for &block in &summand.blocks {
                if block.length == self.d {
                    continue;
                }
                if i == 0 {
                    return Ok((s, summand, block));
                }
                i -= 1;
            }
        }
        Err(Error::plain("LG defect stable block index is out of range"))
    }

    pub(crate) fn map(
        &self,
        model: &mut Model,
        target: &Self,
        indices: [usize; 3],
        map: &Morphism,
    ) -> Result<Value> {
        let [row, col, k] = indices;
        let (t, to, output) = target.block(row)?;
        let (s, from, input) = self.block(col)?;
        let field = Field::new(model, self.d)?;
        if k >= output.length {
            return Err(Error::plain(
                "LG defect stable map power is out of range for the target block",
            ));
        }
        if k >= self.d - input.length {
            return rv_lens::rational_value(model, Scalar::zero());
        }
        let mut result = Scalar::zero();
        for (key, coefficient) in &map.terms {
            model.spend(1)?;
            if key.row != t || key.column != s || !key.slot || key.exponent >= to.j as u128 {
                continue;
            }
            let y = key.exponent as usize;
            let (n, d) = coefficient.small_parts().ok_or_else(|| {
                Error::plain("LG defect map coefficient exceeds the rational cyclotomic range")
            })?;
            let coefficient = rv_lens::fraction(model, &[n], d)?;
            for monomial in 0..from.monomials.len() {
                model.spend(1)?;
                let b = monomial % from.j + y;
                if b >= to.j {
                    continue;
                }
                let image = (monomial / from.j) * to.j + b;
                let start = from.entry(1, monomial, input.start)?;
                let end = to.entry(2, output.start + k, image)?;
                let value = field.multiply(model, &start, &end)?;
                let value = field.multiply(model, &value, &coefficient)?;
                result = field.add(model, &result, &value)?;
            }
        }
        if k < output.length.saturating_sub(input.length) && !result.is_zero() {
            return Err(Error::plain(
                "LG defect action produced a map that is not x-linear",
            ));
        }
        rv_lens::rational_value(model, result)
    }
}
