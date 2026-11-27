use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::mem::size_of;
use std::ops::Range;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;

pub(crate) type Coordinate = (usize, u32);

// Vectors contain only nonzero coefficients and are sorted by coordinate.
// The last coordinate of every stored row is a pivot with coefficient one.
#[derive(Default)]
pub(crate) struct Basis {
    rows: BTreeMap<usize, Row>,
    entries: Entries,
    cells: usize,
}

struct Row {
    start: u32,
    middle: u32,
    end: u32,
}

impl Row {
    fn image(&self) -> Range<usize> {
        self.start as usize..self.middle as usize
    }

    fn witness(&self) -> Range<usize> {
        self.middle as usize..self.end as usize
    }
}

// Store small-field coordinates in one u32. Large primes retain all 32
// coefficient bits in a pair. The arena stores row ranges instead of
// allocating two slice headers for every pivot.
#[derive(Default)]
struct Entries {
    narrow: Vec<u32>,
    wide: Vec<(u32, u32)>,
}

impl Entries {
    const INDEX_BITS: u32 = MAX_PIVOT_CELLS.ilog2() + 1;

    fn narrow(prime: Prime) -> bool {
        prime.get() <= 1 << (u32::BITS - Self::INDEX_BITS)
    }

    fn len(&self) -> usize {
        self.narrow.len() + self.wide.len()
    }

    fn cells(len: usize, prime: Prime) -> usize {
        let width = if Self::narrow(prime) {
            size_of::<u32>()
        } else {
            size_of::<(u32, u32)>()
        };
        (len * width).div_ceil(size_of::<usize>())
    }

    fn append(
        &mut self,
        prime: Prime,
        image: Vec<Coordinate>,
        witness: Vec<Coordinate>,
    ) -> Result<Row> {
        let narrow = Self::narrow(prime);
        let limit = if narrow {
            (1 << Self::INDEX_BITS) - 1
        } else {
            u32::MAX
        } as usize;
        if [image.last(), witness.last()]
            .into_iter()
            .flatten()
            .any(|&(index, _)| index > limit)
        {
            return Err(Error::plain(
                "odd-prime reduction coordinate exceeds its limit",
            ));
        }
        let start = self.len();
        let middle = start + image.len();
        let end = middle + witness.len();
        let entries = image.into_iter().chain(witness);
        if narrow {
            self.narrow
                .extend(entries.map(|(index, value)| index as u32 | value << Self::INDEX_BITS));
        } else {
            self.wide
                .extend(entries.map(|(index, value)| (index as u32, value)));
        }
        Ok(Row {
            start: start as u32,
            middle: middle as u32,
            end: end as u32,
        })
    }

    fn read(
        &self,
        prime: Prime,
        range: Range<usize>,
    ) -> impl ExactSizeIterator<Item = Coordinate> + '_ {
        let narrow = Self::narrow(prime);
        range.map(move |i| {
            if narrow {
                let packed = self.narrow[i];
                let mask = (1 << Self::INDEX_BITS) - 1;
                ((packed & mask) as usize, packed >> Self::INDEX_BITS)
            } else {
                let (index, value) = self.wide[i];
                (index as usize, value)
            }
        })
    }
}

impl Basis {
    pub(crate) fn rank(&self) -> usize {
        self.rows.len()
    }

    pub(crate) fn vectors(&self, prime: Prime) -> impl Iterator<Item = Vec<Coordinate>> + '_ {
        self.rows
            .values()
            .map(move |row| self.entries.read(prime, row.image()).collect())
    }

    pub(crate) fn has_pivot(&self, coordinate: usize) -> bool {
        self.rows.contains_key(&coordinate)
    }

    pub(crate) fn cells(&self) -> usize {
        self.cells
    }

    pub(crate) fn insert(
        &mut self,
        model: &mut Model,
        prime: Prime,
        image: Vec<Coordinate>,
    ) -> Result<()> {
        self.insert_with_witness(model, prime, image, Vec::new())?;
        Ok(())
    }

    // An independent column is normalized before storage. A dependent
    // column returns the accumulated witness whose image is zero.
    pub(crate) fn insert_with_witness(
        &mut self,
        model: &mut Model,
        prime: Prime,
        mut image: Vec<Coordinate>,
        mut witness: Vec<Coordinate>,
    ) -> Result<Option<Vec<Coordinate>>> {
        model.spend(1 + image.len() + witness.len())?;
        while let Some(&(leading, coefficient)) = image.last() {
            if let Some(row) = self.rows.get(&leading) {
                let factor = prime.negate(coefficient);
                image = merge(
                    model,
                    prime,
                    &image,
                    self.entries.read(prime, row.image()),
                    factor,
                )?;
                witness = merge(
                    model,
                    prime,
                    &witness,
                    self.entries.read(prime, row.witness()),
                    factor,
                )?;
            } else {
                // Count the pivot key, row ranges, and retained arena entries.
                let entries = image.len() + witness.len();
                let start = self.entries.len();
                let end = start + entries;
                let cells = self.cells
                    + 1
                    + size_of::<Row>().div_ceil(size_of::<usize>())
                    + Entries::cells(end, prime)
                    - Entries::cells(start, prime);
                check_storage(cells)?;
                if coefficient != 1 {
                    let inverse = invert(model, prime, coefficient)?;
                    model.spend(image.len() + witness.len())?;
                    for (_, value) in image.iter_mut().chain(&mut witness) {
                        *value = prime.mul(*value, inverse);
                    }
                }
                model.spend(entries)?;
                let row = self.entries.append(prime, image, witness)?;
                self.cells = cells;
                self.rows.insert(leading, row);
                return Ok(None);
            }
        }
        Ok(Some(witness))
    }

    pub(crate) fn reduce(
        &self,
        model: &mut Model,
        prime: Prime,
        image: Vec<Coordinate>,
    ) -> Result<Vec<Coordinate>> {
        self.remainder(model, prime, image, None)
    }

    pub(crate) fn reduce_with_witness(
        &self,
        model: &mut Model,
        prime: Prime,
        image: Vec<Coordinate>,
    ) -> Result<(Vec<Coordinate>, Vec<Coordinate>)> {
        let mut witness = Vec::new();
        let remainder = self.remainder(model, prime, image, Some(&mut witness))?;
        Ok((remainder, witness))
    }

    fn remainder(
        &self,
        model: &mut Model,
        prime: Prime,
        mut image: Vec<Coordinate>,
        mut witness: Option<&mut Vec<Coordinate>>,
    ) -> Result<Vec<Coordinate>> {
        let mut remainder = Vec::new();
        while let Some((leading, coefficient)) = image.pop() {
            if let Some(row) = self.rows.get(&leading) {
                image.push((leading, coefficient));
                image = merge(
                    model,
                    prime,
                    &image,
                    self.entries.read(prime, row.image()),
                    prime.negate(coefficient),
                )?;
                if let Some(witness) = witness.as_mut() {
                    // The original vector equals the remainder plus the
                    // images of these positively accumulated coefficients.
                    **witness = merge(
                        model,
                        prime,
                        witness,
                        self.entries.read(prime, row.witness()),
                        coefficient,
                    )?;
                }
            } else {
                model.spend(1)?;
                remainder.push((leading, coefficient));
            }
        }
        remainder.reverse();
        Ok(remainder)
    }
}

#[derive(Default)]
pub(crate) struct QuotientBasis {
    boundaries: Basis,
    cycles: Basis,
    generators: Vec<usize>,
}

impl QuotientBasis {
    pub(crate) fn rank_only(
        model: &mut Model,
        prime: Prime,
        boundaries: &Basis,
        size: usize,
        mut differential: impl FnMut(&mut Model, usize) -> Result<Vec<Coordinate>>,
    ) -> Result<usize> {
        let mut images = Basis::default();
        let mut rank = 0;
        for coordinate in 0..size {
            if boundaries.rows.contains_key(&coordinate) {
                continue;
            }
            // A dependence among the nonpivot coordinates adds one
            // independent class; its witness is unnecessary for the rank.
            let image = differential(model, coordinate)?;
            if images
                .insert_with_witness(model, prime, image, Vec::new())?
                .is_some()
            {
                rank += 1;
            }
        }
        Ok(rank)
    }

    pub(crate) fn build(
        model: &mut Model,
        prime: Prime,
        boundaries: Basis,
        size: usize,
        mut differential: impl FnMut(&mut Model, usize) -> Result<Vec<Coordinate>>,
    ) -> Result<Self> {
        let mut images = Basis::default();
        let mut cycles = Basis::default();
        let mut generators = Vec::new();
        for coordinate in 0..size {
            if boundaries.rows.contains_key(&coordinate) {
                continue;
            }
            // Nonpivot coordinates form a complement to the boundaries.
            // A dependence gives a cycle led by the new coordinate.
            let image = differential(model, coordinate)?;
            if let Some(cycle) =
                images.insert_with_witness(model, prime, image, vec![(coordinate, 1)])?
            {
                let result =
                    cycles.insert_with_witness(model, prime, cycle, vec![(generators.len(), 1)])?;
                debug_assert!(result.is_none());
                generators.push(coordinate);
            }
        }
        model.allocate(boundaries.cells() + cycles.cells() + generators.len())?;
        Ok(Self {
            boundaries,
            cycles,
            generators,
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.generators.len()
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        prime: Prime,
        support: Vec<Coordinate>,
    ) -> Result<Option<Vec<Coordinate>>> {
        let support = self.boundaries.reduce(model, prime, support)?;
        let (remainder, coordinates) = self.cycles.reduce_with_witness(model, prime, support)?;
        Ok(remainder.is_empty().then_some(coordinates))
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        prime: Prime,
        class: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let mut support = Vec::new();
        for &(index, coefficient) in class {
            let pivot = self.generators[index];
            let row = &self.cycles.rows[&pivot];
            support = merge(
                model,
                prime,
                &support,
                self.cycles.entries.read(prime, row.image()),
                coefficient,
            )?;
        }
        model.spend(support.len())?;
        Ok(support)
    }

    pub(crate) fn is_boundary(
        &self,
        model: &mut Model,
        prime: Prime,
        support: Vec<Coordinate>,
    ) -> Result<bool> {
        Ok(self.boundaries.reduce(model, prime, support)?.is_empty())
    }
}

pub(crate) fn add_scaled(
    model: &mut Model,
    prime: Prime,
    left: &[Coordinate],
    right: &[Coordinate],
    factor: u32,
) -> Result<Vec<Coordinate>> {
    merge(model, prime, left, right.iter().copied(), factor)
}

fn merge(
    model: &mut Model,
    prime: Prime,
    left: &[Coordinate],
    right: impl ExactSizeIterator<Item = Coordinate>,
    factor: u32,
) -> Result<Vec<Coordinate>> {
    if factor == 0 {
        model.spend(left.len())?;
        return Ok(left.to_vec());
    }
    let size = left.len() + right.len();
    model.spend(size)?;
    let mut output = Vec::with_capacity(size.min(MAX_PIVOT_CELLS / 2));
    let mut i = 0;
    let mut right = right.peekable();
    while let (Some(&(a, x)), Some(&(b, y))) = (left.get(i), right.peek()) {
        match a.cmp(&b) {
            Ordering::Less => {
                output.push((a, x));
                i += 1;
            }
            Ordering::Equal => {
                let coefficient = prime.add(x, prime.mul(factor, y));
                if coefficient != 0 {
                    output.push((a, coefficient));
                }
                i += 1;
                right.next();
            }
            Ordering::Greater => {
                output.push((b, prime.mul(factor, y)));
                right.next();
            }
        }
    }
    output.extend_from_slice(&left[i..]);
    output.extend(right.map(|(index, value)| (index, prime.mul(value, factor))));
    check_storage(2 * output.len())?;
    Ok(output)
}

// Fermat inversion uses at most two field multiplications per exponent bit.
pub(crate) fn invert(model: &mut Model, prime: Prime, mut value: u32) -> Result<u32> {
    let mut power = prime.get() - 2;
    let mut inverse = 1;
    while power != 0 {
        model.spend(2)?;
        if power & 1 != 0 {
            inverse = prime.mul(inverse, value);
        }
        value = prime.mul(value, value);
        power >>= 1;
    }
    Ok(inverse)
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "odd-prime linear reduction exceeds the storage limit",
        ));
    }
    Ok(())
}
