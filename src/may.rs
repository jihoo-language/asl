use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::hash::Hash;

use crate::cobar::{self, word_cells, Word};
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::linear::{xor, Basis};
use crate::model::{Model, TypeId, Value};
use crate::profile::Profile;

const MAX_BASIS_WORDS: usize = 32_768;

// A basis and its differential columns are sorted, duplicate-free supports.
// The basis must include every term with filtration at least f, and the
// differential preserves the decreasing filtration.
pub(crate) trait FilteredComplex {
    type Term: Clone + Hash + Ord;
    const NAME: &'static str = "May";
    const ELEMENT_NAME: &'static str = "cobar";

    fn basis(&mut self, model: &mut Model, s: usize, f: usize) -> Result<Vec<Self::Term>>;
    fn differential(&mut self, model: &mut Model, term: &Self::Term) -> Result<Vec<Self::Term>>;
    fn differential_support(
        &mut self,
        model: &mut Model,
        terms: &[Self::Term],
    ) -> Result<Vec<Self::Term>>;
    fn filtration(&self, term: &Self::Term) -> Result<usize>;
    fn cells(term: &Self::Term) -> usize;
}

pub(crate) struct FilteredPage<C: FilteredComplex> {
    f: usize,
    cutoff: u128,
    terms: Vec<C::Term>,
    leading_start: usize,
    quotient: Option<Quotient>,
}

#[derive(Default)]
struct Quotient {
    span: Basis,
    representatives: Vec<Box<[usize]>>,
}

impl<C: FilteredComplex> FilteredPage<C> {
    pub(crate) fn empty(r: usize, f: usize) -> Self {
        Self {
            f,
            cutoff: f as u128 + r as u128,
            terms: Vec::new(),
            leading_start: 0,
            quotient: (r != 0).then(Quotient::default),
        }
    }

    pub(crate) fn build(
        model: &mut Model,
        complex: &mut C,
        r: usize,
        s: usize,
        f: usize,
    ) -> Result<Self> {
        let mut page = Self::empty(r, f);
        let mut leading = Vec::new();
        for term in complex.basis(model, s, f)? {
            model.spend(C::cells(&term))?;
            let term_filtration = complex.filtration(&term)?;
            match term_filtration.cmp(&f) {
                Ordering::Equal => leading.push(term),
                // Words at or above the cutoff cannot correct a differential
                // below it, so the canonical lifts need only the intervening words.
                Ordering::Greater if (term_filtration as u128) < page.cutoff => {
                    page.terms.push(term)
                }
                _ => {}
            }
        }
        if leading.is_empty() {
            page.terms = Vec::new();
            return Ok(page);
        }
        page.leading_start = page.terms.len();
        page.terms.extend(leading);

        if r != 0 {
            page.quotient = Some(page.build_quotient(model, complex, r, s)?);
        }
        Ok(page)
    }

    fn build_quotient(
        &self,
        model: &mut Model,
        complex: &mut C,
        r: usize,
        s: usize,
    ) -> Result<Quotient> {
        // In the leading slice of Z_r^f, the kernel of projection is
        // Z_{r-1}^{f+1}. Only the leading parts of boundaries remain to quotient.
        let mut span = self.boundaries(model, complex, r, s)?;
        let mut rows = Rows::new(C::NAME);
        let mut cycles = Basis::default();
        let mut representatives = Vec::new();
        let mut lift_cells = 0;
        for (coordinate, source) in self.terms.iter().enumerate() {
            let mut lower = Vec::new();
            for term in complex.differential(model, source)? {
                model.spend(C::cells(&term))?;
                if (complex.filtration(&term)? as u128) < self.cutoff {
                    lower.push(rows.coordinate(model, C::cells(&term), term)?);
                }
            }
            sort_support(model, &mut lower)?;
            // Processing F^{f+1} first makes kernel witnesses in that subspace
            // invisible to the quotient, but available to correct later lifts.
            let Some(lift) = cycles.insert_with_witness(model, lower, vec![coordinate])? else {
                continue;
            };
            if coordinate < self.leading_start {
                continue;
            }
            let image = lift
                .iter()
                .filter(|&&index| index >= self.leading_start)
                .map(|&index| index - self.leading_start)
                .collect();
            if span
                .insert_with_witness(model, image, vec![representatives.len()])?
                .is_none()
            {
                lift_cells += 1 + lift.len();
                check_storage(lift_cells + span.cells(), C::NAME)?;
                representatives.push(lift.into_boxed_slice());
            }
        }
        model.allocate(lift_cells + span.cells())?;
        Ok(Quotient {
            span,
            representatives,
        })
    }

    fn boundaries(&self, model: &mut Model, complex: &mut C, r: usize, s: usize) -> Result<Basis> {
        let mut boundaries = Basis::default();
        if s == 0 {
            return Ok(boundaries);
        }
        let source_filtration = self.f.saturating_sub(r - 1);
        let mut rows = Rows::new(C::NAME);
        let mut lower_images = Basis::default();
        for source in complex.basis(model, s - 1, source_filtration)? {
            model.spend(C::cells(&source))?;
            if !(source_filtration..=self.f).contains(&complex.filtration(&source)?) {
                continue;
            }
            let (mut lower, mut leading) = (Vec::new(), Vec::new());
            for term in complex.differential(model, &source)? {
                model.spend(C::cells(&term))?;
                match complex.filtration(&term)?.cmp(&self.f) {
                    Ordering::Less => lower.push(rows.coordinate(model, C::cells(&term), term)?),
                    Ordering::Equal => {
                        leading.push(self.leading_coordinate(model, &term)?);
                    }
                    Ordering::Greater => {}
                }
            }
            sort_support(model, &mut lower)?;
            // When the lower image vanishes, the accumulated leading image
            // is the leading part of a boundary lying in F^f.
            if let Some(image) = lower_images.insert_with_witness(model, lower, leading)? {
                boundaries.insert(model, image)?;
            }
        }
        Ok(boundaries)
    }

    fn leading_coordinate(&self, model: &mut Model, term: &C::Term) -> Result<usize> {
        let leading = &self.terms[self.leading_start..];
        model.spend(C::cells(term) * (1 + leading.len().max(1).ilog2() as usize))?;
        leading.binary_search(term).map_err(|_| {
            Error::plain(format!(
                "{} term is missing from the {} filtration slice",
                C::ELEMENT_NAME,
                C::NAME
            ))
        })
    }

    pub(crate) fn rank(&self) -> usize {
        match &self.quotient {
            Some(quotient) => quotient.representatives.len(),
            None => self.terms.len(),
        }
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        complex: &mut C,
        terms: &[C::Term],
    ) -> Result<Vec<usize>> {
        let mut leading = Vec::new();
        let mut below_cutoff = Vec::new();
        for term in terms {
            model.spend(C::cells(term))?;
            let term_filtration = complex.filtration(term)?;
            match term_filtration.cmp(&self.f) {
                Ordering::Less => {
                    return Err(Error::plain(format!(
                        "{} chain lies below the requested {} filtration",
                        C::ELEMENT_NAME,
                        C::NAME,
                    )));
                }
                Ordering::Equal => leading.push(self.leading_coordinate(model, term)?),
                Ordering::Greater => {}
            }
            if (term_filtration as u128) < self.cutoff {
                below_cutoff.push(term.clone());
            }
        }
        if let Some(quotient) = &self.quotient {
            for term in complex.differential_support(model, &below_cutoff)? {
                model.spend(C::cells(&term))?;
                if (complex.filtration(&term)? as u128) < self.cutoff {
                    return Err(Error::plain(format!(
                        "{} chain does not satisfy the {} page cycle condition",
                        C::ELEMENT_NAME,
                        C::NAME,
                    )));
                }
            }
            let (remainder, coordinates) =
                quotient
                    .span
                    .reduce_with_witness(model, leading, Vec::new())?;
            if !remainder.is_empty() {
                return Err(Error::plain(format!(
                    "{} cycle is missing from the {} page basis",
                    C::ELEMENT_NAME,
                    C::NAME,
                )));
            }
            leading = coordinates;
        }
        Ok(leading)
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        class: &[usize],
    ) -> Result<Vec<C::Term>> {
        let support = if let Some(quotient) = &self.quotient {
            let mut support = Vec::new();
            for &index in class {
                let lift = quotient
                    .representatives
                    .get(index)
                    .ok_or_else(generator_range::<C>)?;
                support = xor(model, &support, lift)?;
            }
            support
        } else {
            model.spend(class.len())?;
            class.to_vec()
        };
        let mut terms = Vec::with_capacity(support.len());
        for index in support {
            terms.push(
                self.terms
                    .get(index)
                    .ok_or_else(generator_range::<C>)?
                    .clone(),
            );
        }
        if self.quotient.is_some() {
            let cells: usize = terms.iter().map(C::cells).sum();
            model.spend(cells * (1 + terms.len().max(1).ilog2() as usize))?;
            terms.sort_unstable();
        }
        Ok(terms)
    }
}

// Number only terms that occur in the differential's projected rows.
pub(crate) struct Rows<T> {
    indices: HashMap<T, usize>,
    cells: usize,
    name: &'static str,
}

impl<T: Hash + Eq> Rows<T> {
    pub(crate) fn new(name: &'static str) -> Self {
        Self {
            indices: HashMap::new(),
            cells: 0,
            name,
        }
    }

    pub(crate) fn coordinate(
        &mut self,
        model: &mut Model,
        term_cells: usize,
        term: T,
    ) -> Result<usize> {
        let size = 1 + term_cells;
        model.spend(size)?;
        let next = self.indices.len();
        match self.indices.entry(term) {
            Entry::Occupied(entry) => Ok(*entry.get()),
            Entry::Vacant(entry) => {
                check_storage(self.cells + size, self.name)?;
                self.cells += size;
                Ok(*entry.insert(next))
            }
        }
    }
}

pub(crate) struct Sphere {
    profile: Profile,
    t: usize,
}

impl Sphere {
    pub(crate) fn new(profile: Profile, t: usize) -> Self {
        Self { profile, t }
    }
}

impl FilteredComplex for Sphere {
    type Term = Word;

    fn basis(&mut self, model: &mut Model, s: usize, f: usize) -> Result<Vec<Word>> {
        // Each positive-degree factor has May weight at least one.
        if s > self.t || f > self.t - s || (s == 0 && self.t != 0) {
            return Ok(Vec::new());
        }
        self.profile.basis(model, s, self.t, MAX_BASIS_WORDS)
    }

    fn differential(&mut self, model: &mut Model, word: &Word) -> Result<Vec<Word>> {
        self.profile.differential_word(model, word)
    }

    fn differential_support(&mut self, model: &mut Model, words: &[Word]) -> Result<Vec<Word>> {
        self.profile.differential_support(model, words)
    }

    fn filtration(&self, word: &Word) -> Result<usize> {
        Ok(filtration(self.t, word))
    }

    fn cells(word: &Word) -> usize {
        word_cells(word)
    }
}

pub(crate) struct MayPage {
    profile: Profile,
    s: usize,
    t: usize,
    page: FilteredPage<Sphere>,
}

impl MayPage {
    pub(crate) fn build(
        model: &mut Model,
        profile: Profile,
        r: usize,
        s: usize,
        t: usize,
        f: usize,
    ) -> Result<Self> {
        let page = FilteredPage::build(model, &mut Sphere { profile, t }, r, s, f)?;
        Ok(Self {
            profile,
            s,
            t,
            page,
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.page.rank()
    }

    pub(crate) fn project(&self, model: &mut Model, words: &[Word]) -> Result<Vec<usize>> {
        for word in words {
            model.spend(word_cells(word))?;
            if cobar::word_metadata(word)? != (self.s, self.t) {
                return Err(Error::plain("cobar chain has the wrong May bidegree"));
            }
        }
        let mut sphere = Sphere {
            profile: self.profile,
            t: self.t,
        };
        self.page.project(model, &mut sphere, words)
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<Word>> {
        self.page.representative(model, class)
    }
}

pub(crate) fn class(model: &mut Model, ty: TypeId, words: &[Word]) -> Result<Value> {
    if words.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let coordinates = model.may_page(ty)?.project(model, words)?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

pub(crate) fn differential(
    model: &mut Model,
    source: TypeId,
    coordinates: &[usize],
    target: TypeId,
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let page = model.may_page(source)?;
    let lift = page.representative(model, coordinates)?;
    let image = page.profile.differential_support(model, &lift)?;
    class(model, target, &image)
}

pub(crate) fn multiply(
    model: &mut Model,
    left_ty: TypeId,
    left: &[usize],
    right_ty: TypeId,
    right: &[usize],
    target: TypeId,
) -> Result<Value> {
    if left.is_empty() || right.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let left = model.may_page(left_ty)?.representative(model, left)?;
    let right = model.may_page(right_ty)?.representative(model, right)?;
    let words = cobar::concat_support(model, &left, &right)?;
    class(model, target, &words)
}

pub(crate) fn sort_support(model: &mut Model, support: &mut [usize]) -> Result<()> {
    model.spend(support.len() * (1 + support.len().max(1).ilog2() as usize))?;
    support.sort_unstable();
    Ok(())
}

// Coproduct splits of a generator power preserve or lower its May weight.
// Multiplying the resulting factors can only lower the weight by binary carries.
pub(crate) fn filtration(t: usize, word: &Word) -> usize {
    let weight: usize = word
        .iter()
        .flat_map(|factor| factor.iter().enumerate())
        .map(|(index, &exponent)| (2 * index + 1) * exponent.count_ones() as usize)
        .sum();
    t - weight
}

pub(crate) fn check_storage(cells: usize, name: &str) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(format!(
            "{name} calculation exceeds the storage limit"
        )));
    }
    Ok(())
}

fn generator_range<C: FilteredComplex>() -> Error {
    Error::plain(format!("{} generator index is out of range", C::NAME))
}
