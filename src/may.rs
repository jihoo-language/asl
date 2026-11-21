use std::cmp::Ordering;

use crate::cobar::{self, word_cells, Word, WordRows};
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::linear::{xor, Basis};
use crate::model::{Model, TypeId, Value};

const MAX_BASIS_WORDS: usize = 32_768;

pub(crate) struct MayPage {
    s: usize,
    t: usize,
    f: usize,
    cutoff: usize,
    words: Vec<Word>,
    leading_start: usize,
    quotient: Option<Quotient>,
}

#[derive(Default)]
struct Quotient {
    span: Basis,
    representatives: Vec<Box<[usize]>>,
}

impl MayPage {
    pub(crate) fn build(model: &mut Model, r: usize, s: usize, t: usize, f: usize) -> Result<Self> {
        let mut page = Self {
            s,
            t,
            f,
            cutoff: f.saturating_add(r),
            words: Vec::new(),
            leading_start: 0,
            quotient: (r != 0).then(Quotient::default),
        };
        // Each positive-degree factor has May weight at least one.
        if s > t || f > t - s || (s == 0 && t != 0) {
            return Ok(page);
        }
        let mut leading = Vec::new();
        for word in cobar::basis(model, s, t, MAX_BASIS_WORDS)? {
            model.spend(word_cells(&word))?;
            let word_filtration = filtration(t, &word);
            match word_filtration.cmp(&f) {
                Ordering::Equal => leading.push(word),
                // Words at or above the cutoff cannot correct a differential
                // below it, so the canonical lifts need only the intervening words.
                Ordering::Greater if word_filtration < page.cutoff => page.words.push(word),
                _ => {}
            }
        }
        if leading.is_empty() {
            page.words = Vec::new();
            return Ok(page);
        }
        page.leading_start = page.words.len();
        page.words.extend(leading);

        if r != 0 {
            page.quotient = Some(page.build_quotient(model, r)?);
        }
        Ok(page)
    }

    fn build_quotient(&self, model: &mut Model, r: usize) -> Result<Quotient> {
        // In the leading slice of Z_r^f, the kernel of projection is
        // Z_{r-1}^{f+1}. Only the leading parts of boundaries remain to quotient.
        let mut span = self.boundaries(model, r)?;
        let mut rows = WordRows::default();
        let mut cycles = Basis::default();
        let mut representatives = Vec::new();
        let mut lift_cells = 0;
        for (coordinate, word) in self.words.iter().enumerate() {
            let mut lower = Vec::new();
            for term in cobar::differential_word(model, word)? {
                model.spend(word_cells(&term))?;
                if filtration(self.t, &term) < self.cutoff {
                    lower.push(rows.coordinate(model, term)?);
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
                check_storage(lift_cells + span.cells())?;
                representatives.push(lift.into_boxed_slice());
            }
        }
        model.allocate(lift_cells + span.cells())?;
        Ok(Quotient {
            span,
            representatives,
        })
    }

    fn boundaries(&self, model: &mut Model, r: usize) -> Result<Basis> {
        let mut boundaries = Basis::default();
        if self.s == 0 {
            return Ok(boundaries);
        }
        let source_filtration = self.f.saturating_sub(r - 1);
        let mut rows = WordRows::default();
        let mut lower_images = Basis::default();
        for word in cobar::basis(model, self.s - 1, self.t, MAX_BASIS_WORDS)? {
            model.spend(word_cells(&word))?;
            if !(source_filtration..=self.f).contains(&filtration(self.t, &word)) {
                continue;
            }
            let (mut lower, mut leading) = (Vec::new(), Vec::new());
            for term in cobar::differential_word(model, &word)? {
                model.spend(word_cells(&term))?;
                match filtration(self.t, &term).cmp(&self.f) {
                    Ordering::Less => lower.push(rows.coordinate(model, term)?),
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

    fn leading_coordinate(&self, model: &mut Model, word: &Word) -> Result<usize> {
        let leading = &self.words[self.leading_start..];
        model.spend(word_cells(word) * (1 + leading.len().max(1).ilog2() as usize))?;
        leading
            .binary_search(word)
            .map_err(|_| Error::plain("cobar word is missing from the May filtration slice"))
    }

    pub(crate) fn rank(&self) -> usize {
        match &self.quotient {
            Some(quotient) => quotient.representatives.len(),
            None => self.words.len(),
        }
    }

    pub(crate) fn project(&self, model: &mut Model, words: &[Word]) -> Result<Vec<usize>> {
        let mut leading = Vec::new();
        let mut below_cutoff = Vec::new();
        for word in words {
            model.spend(word_cells(word))?;
            if cobar::word_metadata(word)? != (self.s, self.t) {
                return Err(Error::plain("cobar chain has the wrong May bidegree"));
            }
            let word_filtration = filtration(self.t, word);
            match word_filtration.cmp(&self.f) {
                Ordering::Less => {
                    return Err(Error::plain(
                        "cobar chain lies below the requested May filtration",
                    ));
                }
                Ordering::Equal => leading.push(self.leading_coordinate(model, word)?),
                Ordering::Greater => {}
            }
            if word_filtration < self.cutoff {
                below_cutoff.push(word.clone());
            }
        }
        if let Some(quotient) = &self.quotient {
            for term in cobar::differential_support(model, &below_cutoff)? {
                model.spend(word_cells(&term))?;
                if filtration(self.t, &term) < self.cutoff {
                    return Err(Error::plain(
                        "cobar chain does not satisfy the May page cycle condition",
                    ));
                }
            }
            let (remainder, coordinates) =
                quotient
                    .span
                    .reduce_with_witness(model, leading, Vec::new())?;
            if !remainder.is_empty() {
                return Err(Error::plain(
                    "cobar cycle is missing from the May page basis",
                ));
            }
            leading = coordinates;
        }
        Ok(leading)
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<Word>> {
        let support = if let Some(quotient) = &self.quotient {
            let mut support = Vec::new();
            for &index in class {
                let lift = quotient
                    .representatives
                    .get(index)
                    .ok_or_else(generator_range)?;
                support = xor(model, &support, lift)?;
            }
            support
        } else {
            model.spend(class.len())?;
            class.to_vec()
        };
        let mut words = Vec::with_capacity(support.len());
        for index in support {
            words.push(self.words.get(index).ok_or_else(generator_range)?.clone());
        }
        if self.quotient.is_some() {
            let cells: usize = words.iter().map(word_cells).sum();
            model.spend(cells * (1 + words.len().max(1).ilog2() as usize))?;
            words.sort_unstable();
        }
        Ok(words)
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
    let lift = model.may_page(source)?.representative(model, coordinates)?;
    let image = cobar::differential_support(model, &lift)?;
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

fn sort_support(model: &mut Model, support: &mut [usize]) -> Result<()> {
    model.spend(support.len() * (1 + support.len().max(1).ilog2() as usize))?;
    support.sort_unstable();
    Ok(())
}

// Coproduct splits of a generator power preserve or lower its May weight.
// Multiplying the resulting factors can only lower the weight by binary carries.
fn filtration(t: usize, word: &Word) -> usize {
    let weight: usize = word
        .iter()
        .flat_map(|factor| factor.iter().enumerate())
        .map(|(index, &exponent)| (2 * index + 1) * exponent.count_ones() as usize)
        .sum();
    t - weight
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain("May calculation exceeds the storage limit"));
    }
    Ok(())
}

fn generator_range() -> Error {
    Error::plain("May generator index is out of range")
}
