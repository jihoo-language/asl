use std::cell::OnceCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use super::{widen, Algebra, Class, Cusp, Element as Spin, Terms, DIMENSION, EVEN_RANK, RANK};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::rational_sparse::{
    decimal, Coefficients, Domain, Integer, Magnitude, WORK_LIMBS,
};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Expression {
    Unit(Label, usize, usize),
    Volume(Label),
    Embed(Label, Box<ParsedExpr>),
    Left(Box<ParsedExpr>),
    Compose(Box<[ParsedExpr; 2]>),
    Scale(i128, u128, Box<ParsedExpr>),
    Project(Box<ParsedExpr>),
}

pub(crate) enum Command {
    VolumeSquare(Label),
    Coefficient(usize, usize, Box<ParsedExpr>),
    Trace(Box<ParsedExpr>),
    TracePair(Box<[ParsedExpr; 2]>),
}

pub(crate) enum Endomorphism {}

impl Domain for Endomorphism {
    const NAME: &'static str = "Kuga–Satake endomorphism";
}

pub(crate) type End = Coefficients<Endomorphism>;
type Gram = [[Integer; RANK]; RANK];

impl End {
    pub(crate) fn unit(i: usize, j: usize) -> Result<Self> {
        Ok(Self::singleton(checked_key(i, j)?, 1))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        self.charge_storage(model)?;
        Ok(Value::KSEnd(Rc::new(self)))
    }

    pub(crate) fn entry(&self, model: &mut Model, i: usize, j: usize) -> Result<String> {
        self.coefficient(model, checked_key(i, j)?)
    }

    pub(crate) fn compose(&self, model: &mut Model, right: &Self) -> Result<Self> {
        let right_terms = right.terms();
        let offsets: [usize; EVEN_RANK + 1] =
            std::array::from_fn(|i| right_terms.partition_point(|&(index, _)| index < key(i, 0)));
        let mut output = Terms::new();
        for &(index, coefficient) in self.terms() {
            let (i, j) = indices(index);
            let row = &right_terms[offsets[j]..offsets[j + 1]];
            model.spend((1 + row.len()) * WORK_LIMBS)?;
            let coefficient: Integer = coefficient.resize()?;
            for &(index, entry) in row {
                let (_, k) = indices(index);
                Self::accumulate(
                    &mut output,
                    key(i, k),
                    coefficient.multiply_wide(entry.resize()?)?,
                )?;
            }
        }
        Self::from_terms(model, &output, self.product_denominator(right)?)
    }

    pub(crate) fn trace(&self, model: &mut Model) -> Result<String> {
        let mut numerator = Integer::from_i128(0);
        model.spend(self.terms().len() * WORK_LIMBS)?;
        for &(index, value) in self.terms() {
            if index == transpose(index) {
                numerator = numerator.add(value.resize()?)?;
            }
        }
        rational_text(model, numerator, self.denominator_wide()?)
    }

    pub(crate) fn trace_pair(&self, model: &mut Model, right: &Self) -> Result<String> {
        let terms = right.terms();
        let numerator = self.trace_numerator(model, terms.len(), |index| {
            terms
                .binary_search_by_key(&index, |&(key, _)| key)
                .ok()
                .map(|i| terms[i].1.resize())
                .transpose()
        })?;
        rational_text(model, numerator, self.product_denominator(right)?)
    }

    fn trace_numerator(
        &self,
        model: &mut Model,
        entries: usize,
        lookup: impl Fn(u32) -> Result<Option<Integer>>,
    ) -> Result<Integer> {
        let mut sum = Integer::from_i128(0);
        model.spend(self.terms().len() * (WORK_LIMBS + entries.max(1).ilog2() as usize))?;
        for &(index, value) in self.terms() {
            if let Some(other) = lookup(transpose(index))? {
                sum = sum.add(other.multiply_wide(value.resize()?)?)?;
            }
        }
        Ok(sum)
    }

    fn from_terms(model: &mut Model, terms: &Terms, denominator: Magnitude) -> Result<Self> {
        Self::normalized(
            model,
            terms.iter().map(|(&key, entry)| (key, entry)),
            denominator,
        )
    }
}

fn checked_key(i: usize, j: usize) -> Result<u32> {
    if i >= EVEN_RANK || j >= EVEN_RANK {
        return Err(Error::plain(
            "Kuga–Satake matrix indices must be between 0 and 63",
        ));
    }
    Ok(key(i, j))
}

fn key(i: usize, j: usize) -> u32 {
    (i * EVEN_RANK + j) as u32
}

fn indices(index: u32) -> (usize, usize) {
    let index = index as usize;
    (index / EVEN_RANK, index % EVEN_RANK)
}

fn transpose(index: u32) -> u32 {
    let (i, j) = indices(index);
    key(j, i)
}

pub(super) struct Marking {
    volume: Terms,
    gram: Gram,
    determinant: Integer,
    square: Integer,
    kernel: OnceCell<Box<Kernel>>,
}

struct Kernel {
    adjugate: Gram,
    images: [Terms; RANK],
}

impl Algebra {
    fn projector_marking(&self, model: &mut Model) -> Result<&Marking> {
        if let Some(marking) = self.projector.get() {
            return Ok(marking);
        }
        let mut gram = self.gram;
        for (i, row) in gram.iter_mut().enumerate() {
            row[i] = row[i].multiply(2)?;
        }
        let all = DIMENSION - 1;
        let determinant = determinant(model, &gram, all, all)?;
        if determinant.is_zero() {
            return Err(Error::plain(
                "polarized Kuga–Satake Gram matrix is singular",
            ));
        }
        let square = determinant
            .negate()
            .divide_exact(model, Integer::from_i128(2))?
            .ok_or_else(|| Error::plain("polarized rank-seven Gram determinant must be even"))?;

        let mut exterior = vec![Terms::new(); DIMENSION];
        exterior[0].insert(0, Integer::from_i128(1));
        for mask in 1..DIMENSION {
            let i = mask.trailing_zeros() as usize;
            let rest = mask ^ (1 << i);
            let mut word = self.generator_left(model, i, &exterior[rest])?;
            // Encode a k-blade as 2^floor(k/2) times its alternating product.
            // The wedge recursion subtracts contractions by G/2.
            if mask.count_ones().is_multiple_of(2) {
                for value in word.values_mut() {
                    *value = value.multiply(2)?;
                }
            }
            let mut remaining = rest;
            let mut negative = true;
            while remaining != 0 {
                let j = remaining.trailing_zeros() as usize;
                let pairing = gram[i][j];
                let contraction = &exterior[rest ^ (1 << j)];
                model.spend((1 + contraction.len()) * WORK_LIMBS)?;
                if !pairing.is_zero() {
                    for (&blade, &entry) in contraction {
                        let value = pairing.multiply_wide(entry)?;
                        Spin::accumulate(
                            &mut word,
                            blade,
                            if negative { value.negate() } else { value },
                        )?;
                    }
                }
                negative = !negative;
                remaining &= remaining - 1;
            }
            exterior[mask] = word;
        }
        let volume = exterior.pop().expect("the top Clifford volume is present");
        model.allocate((RANK * RANK + volume.len() + 2) * (WORK_LIMBS + 3))?;
        Ok(self.projector.get_or_init(|| {
            Box::new(Marking {
                volume,
                gram,
                determinant,
                square,
                kernel: OnceCell::new(),
            })
        }))
    }

    pub(crate) fn volume(&self, model: &mut Model) -> Result<Spin> {
        let marking = self.projector_marking(model)?;
        Spin::integral(model, &marking.volume)
    }

    pub(crate) fn volume_square(&self, model: &mut Model) -> Result<String> {
        let marking = self.projector_marking(model)?;
        rational_text(model, marking.square, Magnitude::from_small(1))
    }

    pub(crate) fn left_end(&self, model: &mut Model, spin: &Spin) -> Result<End> {
        let terms = self.left_terms(model, &widen(spin)?)?;
        End::from_terms(model, &terms, Magnitude::from_small(1))
    }

    pub(crate) fn embed(&self, model: &mut Model, cusp: &Cusp, vector: &Class) -> Result<End> {
        if !vector.wide_pair(model, &cusp.polarization)?.is_zero() {
            return Err(Error::plain(
                "Kuga–Satake embedding requires a polarized vector",
            ));
        }
        let vector = self.vector_terms(model, vector)?;
        let marking = self.projector_marking(model)?;
        let terms = self.product(model, &marking.volume, &vector)?;
        let output = self.left_terms(model, &terms)?;
        End::from_terms(model, &output, Magnitude::from_small(1))
    }

    pub(crate) fn project(&self, model: &mut Model, end: &End) -> Result<End> {
        if end.terms().is_empty() {
            return Ok(End::zero());
        }
        let marking = self.projector_marking(model)?;
        let kernel = self.projector_kernel(model, marking)?;
        let mut traces = [Integer::from_i128(0); RANK];
        for (trace, image) in traces.iter_mut().zip(&kernel.images) {
            *trace = end.trace_numerator(model, image.len(), |key| Ok(image.get(&key).copied()))?;
        }
        let divisor = marking
            .determinant
            .multiply_wide(marking.square)?
            .multiply(32)?;
        let denominator = divisor
            .magnitude
            .checked_product(end.denominator_wide()?)
            .ok_or_else(End::arithmetic)?;
        let mut output = Terms::new();
        model.spend(RANK * RANK * WORK_LIMBS)?;
        for (row, image) in kernel.adjugate.iter().zip(&kernel.images) {
            let mut coordinate = Integer::from_i128(0);
            for (&entry, &trace) in row.iter().zip(&traces) {
                coordinate = coordinate.add(entry.multiply_wide(trace)?)?;
            }
            if divisor.negative {
                coordinate = coordinate.negate();
            }
            if !coordinate.is_zero() {
                model.spend(image.len() * WORK_LIMBS)?;
                for (&key, &value) in image {
                    End::accumulate(&mut output, key, coordinate.multiply_wide(value)?)?;
                }
            }
        }
        End::from_terms(model, &output, denominator)
    }

    fn projector_kernel<'a>(&self, model: &mut Model, marking: &'a Marking) -> Result<&'a Kernel> {
        if let Some(kernel) = marking.kernel.get() {
            return Ok(kernel);
        }
        let mut adjugate = [[Integer::from_i128(0); RANK]; RANK];
        let all = DIMENSION - 1;
        for i in 0..RANK {
            let (previous, remaining) = adjugate.split_at_mut(i);
            let row = &mut remaining[0];
            let columns = all ^ (1 << i);
            row[i] = determinant(model, &marking.gram, columns, columns)?;
            for (j, previous_row) in previous.iter_mut().enumerate() {
                let mut entry = determinant(model, &marking.gram, all ^ (1 << j), columns)?;
                if !(i + j).is_multiple_of(2) {
                    entry = entry.negate();
                }
                row[j] = entry;
                previous_row[i] = entry;
            }
        }
        let mut images: [Terms; RANK] = std::array::from_fn(|_| BTreeMap::new());
        for (i, image) in images.iter_mut().enumerate() {
            let vector = Terms::from([(1 << i, Integer::from_i128(1))]);
            let even = self.product(model, &marking.volume, &vector)?;
            *image = self.left_terms(model, &even)?;
        }
        let entries = images.iter().map(BTreeMap::len).sum::<usize>() + RANK * RANK;
        model.allocate(entries * (WORK_LIMBS + 3))?;
        Ok(marking
            .kernel
            .get_or_init(|| Box::new(Kernel { adjugate, images })))
    }

    fn left_terms(&self, model: &mut Model, terms: &Terms) -> Result<Terms> {
        let columns = self.regular_columns(model, terms)?;
        let count = columns.iter().map(Vec::len).sum::<usize>();
        model.spend(count * WORK_LIMBS)?;
        Ok(columns
            .into_iter()
            .enumerate()
            .flat_map(|(j, column)| column.into_iter().map(move |(i, value)| (key(i, j), value)))
            .collect())
    }
}

// Laplace expansion on column subsets avoids division and retains wide cancellation.
fn determinant(model: &mut Model, gram: &Gram, rows: usize, columns: usize) -> Result<Integer> {
    let mut minors = [Integer::from_i128(0); DIMENSION];
    minors[0] = Integer::from_i128(1);
    for (length, (_, row)) in gram
        .iter()
        .enumerate()
        .rev()
        .filter(|(i, _)| rows & (1 << i) != 0)
        .enumerate()
    {
        let mut subset = columns;
        while subset != 0 {
            if subset.count_ones() as usize == length + 1 {
                model.spend((length + 1) * WORK_LIMBS)?;
                let mut remaining = subset;
                let mut sum = Integer::from_i128(0);
                let mut negative = false;
                while remaining != 0 {
                    let j = remaining.trailing_zeros() as usize;
                    let product = row[j].multiply_wide(minors[subset ^ (1 << j)])?;
                    sum = sum.add(if negative { product.negate() } else { product })?;
                    negative = !negative;
                    remaining &= remaining - 1;
                }
                minors[subset] = sum;
            }
            subset = (subset - 1) & columns;
        }
    }
    Ok(minors[columns])
}

fn rational_text(model: &mut Model, numerator: Integer, denominator: Magnitude) -> Result<String> {
    if numerator.is_zero() {
        return Ok("0".into());
    }
    let gcd = numerator.magnitude.gcd(model, denominator)?;
    let n = numerator
        .magnitude
        .divide_exact(model, gcd)?
        .expect("gcd divides numerator");
    let d = denominator
        .divide_exact(model, gcd)?
        .expect("gcd divides denominator");
    let sign = if numerator.negative { "-" } else { "" };
    if d.small() == Some(1) {
        Ok(format!("{sign}{}", decimal(n)))
    } else {
        Ok(format!("{sign}{}/{}", decimal(n), decimal(d)))
    }
}
