use std::cell::OnceCell;
use std::cmp::Ordering;
use std::collections::{hash_map::Entry, HashMap};
use std::mem::{size_of, size_of_val};

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::odd_cobar::{self, word_cells, Differential, Factor, Prime, Term, Word};
use crate::odd_profile;
use crate::prime_linear::{add_scaled, Basis, Coordinate, QuotientBasis};
use crate::profile::Profile;

const MAX_BASIS_WORDS: usize = 32_768;

pub(crate) struct OddExt {
    prime: Prime,
    profile: Profile,
    s: usize,
    t: usize,
    words: Vec<Word>,
    quotient: QuotientBasis,
    preimages: OnceCell<BoundaryPreimages>,
}

struct BoundaryPreimages {
    words: Vec<Word>,
    images: Basis,
}

impl OddExt {
    pub(crate) fn build(
        model: &mut Model,
        prime: Prime,
        profile: Profile,
        s: usize,
        t: usize,
    ) -> Result<Self> {
        let words = odd_cobar::basis_profile(model, prime, profile, s, t, MAX_BASIS_WORDS)?;
        let quotient = if words.is_empty() {
            QuotientBasis::default()
        } else {
            let boundaries = boundaries(model, prime, profile, s, t, &words)?;
            let mut rows = WordRows::default();
            let mut differential = Differential::new(prime, profile);
            QuotientBasis::build(model, prime, boundaries, words.len(), |model, index| {
                rows.differential(model, &mut differential, &words[index], true)
            })?
        };
        Ok(Self {
            prime,
            profile,
            s,
            t,
            words,
            quotient,
            preimages: OnceCell::new(),
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.quotient.rank()
    }

    pub(crate) fn project(&self, model: &mut Model, terms: &[Term]) -> Result<Vec<Coordinate>> {
        let column = word_indices(model, &self.words, terms)?;
        self.quotient
            .project(model, self.prime, column)?
            .ok_or_else(|| Error::plain("extClass requires a closed odd-prime cobar cochain"))
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        class: &[Coordinate],
    ) -> Result<Vec<Term>> {
        let support = self.quotient.representative(model, self.prime, class)?;
        Ok(decode(&self.words, support))
    }

    pub(crate) fn primitive(&self, model: &mut Model, terms: &[Term]) -> Result<Option<Vec<Term>>> {
        if terms.is_empty() {
            return Ok(Some(Vec::new()));
        }
        let column = word_indices(model, &self.words, terms)?;
        if !self
            .quotient
            .is_boundary(model, self.prime, column.clone())?
        {
            return Ok(None);
        }
        let preimages = if let Some(preimages) = self.preimages.get() {
            preimages
        } else {
            let mut images = Basis::default();
            let words = boundary_columns(
                model,
                self.prime,
                self.profile,
                self.s,
                self.t,
                &self.words,
                |model, index, differential| {
                    images.insert_with_witness(
                        model,
                        self.prime,
                        differential,
                        vec![(index, 1)],
                    )?;
                    Ok(())
                },
            )?;
            model.allocate(images.cells() + words.iter().map(word_cells).sum::<usize>())?;
            self.preimages
                .get_or_init(|| BoundaryPreimages { words, images })
        };
        let (remainder, support) = preimages
            .images
            .reduce_with_witness(model, self.prime, column)?;
        debug_assert!(remainder.is_empty());
        Ok(Some(decode(&preimages.words, support)))
    }
}

pub(crate) fn rank(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    s: usize,
    t: usize,
) -> Result<usize> {
    let words = odd_cobar::temporary_basis(model, prime, profile, s, t, MAX_BASIS_WORDS)?;
    if words.is_empty() {
        return Ok(0);
    }
    let boundaries = boundaries(model, prime, profile, s, t, &words)?;
    let mut rows = WordRows::default();
    let mut differential = Differential::new(prime, profile);
    QuotientBasis::rank_only(model, prime, &boundaries, words.len(), |model, index| {
        rows.differential(model, &mut differential, &words[index], true)
    })
}

pub(crate) fn canonical_primitive(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    (s, t): (usize, usize),
    terms: &[Term],
) -> Result<Option<Vec<Term>>> {
    if terms.is_empty() {
        return Ok(Some(Vec::new()));
    }
    let ty = model.profile_odd_ext_type(prime, profile, s, t)?;
    if let Some(group) = model.cached_odd_ext(ty) {
        return group.primitive(model, terms);
    }
    let Some(previous) = s.checked_sub(1) else {
        return Ok(None);
    };
    let source = odd_cobar::temporary_basis(model, prime, profile, previous, t, MAX_BASIS_WORDS)?;
    let mut rows = WordRows::default();
    let target = rows.encode(model, terms.to_vec())?;
    let mut images = Basis::default();
    let mut differential = Differential::new(prime, profile);
    for (index, word) in source.iter().enumerate() {
        let column = rows.differential(model, &mut differential, word, true)?;
        images.insert_with_witness(model, prime, column, vec![(index, 1)])?;
    }
    // Independent source columns are selected in the same order as in
    // OddExt::primitive. Their expression of a boundary is unique, even
    // though the target rows here are interned instead of enumerated.
    let (remainder, support) = images.reduce_with_witness(model, prime, target)?;
    Ok(remainder.is_empty().then(|| decode(&source, support)))
}

fn boundaries(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    s: usize,
    t: usize,
    target: &[Word],
) -> Result<Basis> {
    let mut basis = Basis::default();
    boundary_columns(model, prime, profile, s, t, target, |model, _, column| {
        basis.insert(model, prime, column)
    })?;
    Ok(basis)
}

fn boundary_columns(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    s: usize,
    t: usize,
    target: &[Word],
    mut process: impl FnMut(&mut Model, usize, Vec<Coordinate>) -> Result<()>,
) -> Result<Vec<Word>> {
    let Some(previous) = s.checked_sub(1) else {
        return Ok(Vec::new());
    };
    let source = odd_cobar::temporary_basis(model, prime, profile, previous, t, MAX_BASIS_WORDS)?;
    if source.is_empty() {
        return Ok(source);
    }
    let mut rows = WordRows::default();
    for (index, word) in target.iter().enumerate() {
        let factors = rows.factor_indices(model, word)?;
        let actual = rows.row(model, &factors, true)?;
        debug_assert_eq!(actual, index);
    }
    let mut differential = Differential::new(prime, profile);
    for (index, word) in source.iter().enumerate() {
        let column = rows.differential(model, &mut differential, word, false)?;
        process(model, index, column)?;
    }
    Ok(source)
}

pub(crate) fn word_indices(
    model: &mut Model,
    basis: &[Word],
    terms: &[Term],
) -> Result<Vec<Coordinate>> {
    let search_steps = 1 + basis.len().max(1).ilog2() as usize;
    terms
        .iter()
        .map(|(word, coefficient)| {
            model.spend(search_steps * word_cells(word))?;
            let index = basis
                .binary_search(word)
                .map_err(|_| Error::plain("odd-prime word is outside its cobar bidegree"))?;
            Ok((index, *coefficient))
        })
        .collect()
}

pub(crate) fn decode(words: &[Word], coordinates: Vec<Coordinate>) -> Vec<Term> {
    coordinates
        .into_iter()
        .map(|(index, coefficient)| (words[index].clone(), coefficient))
        .collect()
}

// Intern each monomial once. Outgoing reductions only index words that
// occur in a differential.
#[derive(Default)]
pub(crate) struct WordRows {
    factors: HashMap<Factor, u32>,
    short: HashMap<u64, usize>,
    long: HashMap<Box<[u32]>, usize>,
    cells: usize,
}

impl WordRows {
    fn store(cells: &mut usize, size: usize) -> Result<()> {
        if size > MAX_PIVOT_CELLS - *cells {
            return Err(Error::plain(
                "odd-prime cobar reduction exceeds the row storage limit",
            ));
        }
        *cells += size;
        Ok(())
    }

    fn factor_index(&mut self, model: &mut Model, factor: &Factor) -> Result<u32> {
        model.spend(1 + factor.cells())?;
        let next = (1 + self.factors.len()) as u32; // Bounded by the local storage cap.
        match self.factors.entry(factor.clone()) {
            Entry::Occupied(entry) => Ok(*entry.get()),
            Entry::Vacant(entry) => {
                Self::store(&mut self.cells, 3 + factor.cells())?;
                Ok(*entry.insert(next))
            }
        }
    }

    fn factor_indices(&mut self, model: &mut Model, word: &Word) -> Result<Vec<u32>> {
        word.iter()
            .map(|factor| self.factor_index(model, factor))
            .collect()
    }

    fn row(&mut self, model: &mut Model, indices: &[u32], insert: bool) -> Result<usize> {
        model.spend(1 + indices.len())?;
        let next = self.short.len() + self.long.len();
        let outside = || Error::plain("odd-prime word is outside its cobar bidegree");
        if indices.len() <= size_of::<u64>() && indices.iter().all(|&index| index <= u8::MAX as u32)
        {
            // Factor indices start at one, so padding with zero is unambiguous.
            let mut bytes = [0; size_of::<u64>()];
            for (byte, &index) in bytes.iter_mut().zip(indices) {
                *byte = index as u8;
            }
            let key = u64::from_le_bytes(bytes);
            return match self.short.entry(key) {
                Entry::Occupied(entry) => Ok(*entry.get()),
                Entry::Vacant(_) if !insert => Err(outside()),
                Entry::Vacant(entry) => {
                    Self::store(
                        &mut self.cells,
                        2 + size_of::<u64>().div_ceil(size_of::<usize>()),
                    )?;
                    Ok(*entry.insert(next))
                }
            };
        }
        if let Some(&index) = self.long.get(indices) {
            return Ok(index);
        }
        if !insert {
            return Err(outside());
        }
        let payload = size_of_val(indices).div_ceil(size_of::<usize>());
        Self::store(&mut self.cells, 4 + payload)?;
        self.long.insert(indices.into(), next);
        Ok(next)
    }

    fn differential(
        &mut self,
        model: &mut Model,
        differential: &mut Differential,
        word: &Word,
        insert: bool,
    ) -> Result<Vec<Coordinate>> {
        let factors = self.factor_indices(model, word)?;
        let mut split = Vec::with_capacity(word.len() + 1);
        let mut column = Vec::new();
        differential.for_each_split(model, word, |model, position, left, right, value| {
            if column.len() >= MAX_PIVOT_CELLS / 2 {
                return Err(Error::plain(
                    "odd-prime linear reduction exceeds the storage limit",
                ));
            }
            split.clear();
            split.extend_from_slice(&factors[..position]);
            split.push(self.factor_index(model, left)?);
            split.push(self.factor_index(model, right)?);
            split.extend_from_slice(&factors[position + 1..]);
            column.push((self.row(model, &split, insert)?, value));
            Ok(())
        })?;
        // Different split positions cannot give the same word: the first
        // changed factor always has strictly lower degree.
        Self::sort_column(model, &mut column)?;
        debug_assert!(column.windows(2).all(|pair| pair[0].0 != pair[1].0));
        Ok(column)
    }

    pub(crate) fn encode(
        &mut self,
        model: &mut Model,
        terms: Vec<Term>,
    ) -> Result<Vec<Coordinate>> {
        let mut column = Vec::with_capacity(terms.len());
        for (word, coefficient) in terms {
            let factors = self.factor_indices(model, &word)?;
            column.push((self.row(model, &factors, true)?, coefficient));
        }
        Self::sort_column(model, &mut column)?;
        Ok(column)
    }

    fn sort_column(model: &mut Model, column: &mut [Coordinate]) -> Result<()> {
        model.spend(column.len() * (1 + column.len().max(1).ilog2() as usize))?;
        column.sort_unstable();
        Ok(())
    }
}

pub(crate) fn zero() -> Value {
    Value::OddClass([].into())
}

pub(crate) fn persist(model: &mut Model, coordinates: Vec<Coordinate>) -> Result<Value> {
    model.allocate(2 * coordinates.len())?;
    Ok(Value::OddClass(coordinates.into()))
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(zero());
    }
    let coordinates = model.odd_ext(ty)?.project(model, terms)?;
    persist(model, coordinates)
}

pub(crate) fn restrict(
    model: &mut Model,
    source_ty: TypeId,
    coordinates: &[Coordinate],
    target_ty: TypeId,
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(zero());
    }
    let (p, Profile::An(n), _, _) = model.kind(target_ty).pext_grade().unwrap() else {
        unreachable!()
    };
    let representative = model
        .odd_ext(source_ty)?
        .representative(model, coordinates)?;
    let projected = odd_profile::restrict(model, p, n, &representative)?;
    class(model, target_ty, &projected)
}

pub(crate) fn add(
    model: &mut Model,
    prime: Prime,
    left: &[Coordinate],
    right: &[Coordinate],
) -> Result<Value> {
    let coordinates = add_scaled(model, prime, left, right, 1)?;
    persist(model, coordinates)
}

pub(crate) fn scale(
    model: &mut Model,
    prime: Prime,
    scalar: usize,
    class: &[Coordinate],
) -> Result<Value> {
    let coordinates = add_scaled(model, prime, &[], class, prime.reduce(scalar))?;
    persist(model, coordinates)
}

pub(crate) fn elements(
    model: &mut Model,
    prime: Prime,
    rank: usize,
    limit: usize,
) -> Result<Vec<Value>> {
    let p = prime.get() as usize;
    let mut size = 1;
    for _ in 0..rank {
        if size > limit / p {
            return Err(Error::plain(format!(
                "type exceeds the enumeration limit of {limit} elements"
            )));
        }
        size *= p;
    }
    model.allocate(size)?;
    model.spend(size * (rank + 1))?;
    let mut values = Vec::with_capacity(size);
    for mut index in 0..size {
        let mut coordinates = Vec::new();
        for generator in 0..rank {
            let coefficient = (index % p) as u32;
            if coefficient != 0 {
                coordinates.push((generator, coefficient));
            }
            index /= p;
        }
        values.push(persist(model, coordinates)?);
    }
    Ok(values)
}

pub(crate) fn compare(
    model: &mut Model,
    left: &[Coordinate],
    right: &[Coordinate],
) -> Result<Ordering> {
    model.spend(1 + left.len().min(right.len()))?;
    Ok(left.iter().rev().cmp(right.iter().rev()))
}

pub(crate) fn yoneda(
    model: &mut Model,
    left_ty: TypeId,
    left: &[Coordinate],
    right_ty: TypeId,
    right: &[Coordinate],
    target: TypeId,
) -> Result<Value> {
    if left.is_empty() || right.is_empty() {
        return Ok(zero());
    }
    let left_group = model.odd_ext(left_ty)?;
    let left = left_group.representative(model, left)?;
    let right = model.odd_ext(right_ty)?.representative(model, right)?;
    let terms = odd_cobar::concat_support(model, left_group.prime, &left, &right)?;
    class(model, target, &terms)
}
