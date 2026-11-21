use std::cell::OnceCell;
use std::cmp::Ordering;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::odd_cobar::{self, word_cells, Prime, Term, Word};
use crate::odd_ext::{self, decode, word_indices, WordRows};
use crate::odd_profile;
use crate::prime_linear::{add_scaled, Basis, Coordinate};
use crate::profile::Profile;

const MAX_BASIS_WORDS: usize = 32_768;

pub(crate) struct OddMay {
    prime: Prime,
    profile: Profile,
    r: usize,
    s: usize,
    t: usize,
    f: usize,
    cutoff: usize,
    words: Vec<Word>,
    leading_start: usize,
    span: Basis,
    representatives: Vec<Box<[Coordinate]>>,
    transition: OnceCell<Basis>,
}

impl OddMay {
    pub(crate) fn build(
        model: &mut Model,
        prime: Prime,
        profile: Profile,
        r: usize,
        s: usize,
        t: usize,
        f: usize,
    ) -> Result<Self> {
        let mut page = Self {
            prime,
            profile,
            r,
            s,
            t,
            f,
            cutoff: f.saturating_add(r),
            words: Vec::new(),
            leading_start: 0,
            span: Basis::default(),
            representatives: Vec::new(),
            transition: OnceCell::new(),
        };
        let mut leading = Vec::new();
        for word in page.basis(model, s, f)? {
            let filtration = page.filtration(model, &word)?;
            if filtration == f {
                leading.push(word);
            } else if (f..page.cutoff).contains(&filtration) {
                page.words.push(word);
            }
        }
        if leading.is_empty() {
            page.words = Vec::new();
            return Ok(page);
        }
        page.leading_start = page.words.len();
        page.words.extend(leading);
        if r > 0 {
            page.build_quotient(model)?;
        }
        Ok(page)
    }

    fn basis(&self, model: &mut Model, s: usize, f: usize) -> Result<Vec<Word>> {
        // Each positive-degree factor has at least one unit of May weight.
        if s > self.t || f > self.t - s {
            return Ok(Vec::new());
        }
        odd_cobar::basis_profile(model, self.prime, self.profile, s, self.t, MAX_BASIS_WORDS)
    }

    fn build_quotient(&mut self, model: &mut Model) -> Result<()> {
        let prime = self.prime;
        let mut span = self.boundaries(model)?;
        let mut rows = WordRows::default();
        let mut cycles = Basis::default();
        let mut cells = 0;
        for (index, word) in self.words.iter().enumerate() {
            let mut lower = Vec::new();
            for term in odd_cobar::differential_profile_word(model, prime, self.profile, word)? {
                if self.filtration(model, &term.0)? < self.cutoff {
                    lower.push(term);
                }
            }
            let column = rows.encode(model, lower)?;
            // Process the higher-filtration correction words first. Their
            // dependencies contribute nothing to the leading quotient.
            let Some(lift) = cycles.insert_with_witness(model, prime, column, vec![(index, 1)])?
            else {
                continue;
            };
            if index < self.leading_start {
                continue;
            }
            let leading = lift
                .iter()
                .filter(|&&(index, _)| index >= self.leading_start)
                .map(|&(index, coefficient)| (index - self.leading_start, coefficient))
                .collect();
            let witness = vec![(self.representatives.len(), 1)];
            if span
                .insert_with_witness(model, prime, leading, witness)?
                .is_none()
            {
                cells += 1 + 2 * lift.len();
                check_storage(cells + span.cells())?;
                self.representatives.push(lift.into_boxed_slice());
            }
        }
        model.allocate(cells + span.cells())?;
        self.span = span;
        Ok(())
    }

    fn boundaries(&self, model: &mut Model) -> Result<Basis> {
        let mut boundaries = Basis::default();
        let Some(previous) = self.s.checked_sub(1) else {
            return Ok(boundaries);
        };
        let prime = self.prime;
        let input_filtration = self.f.saturating_sub(self.r - 1);
        let mut rows = WordRows::default();
        let mut lower_images = Basis::default();
        for word in self.basis(model, previous, input_filtration)? {
            if !(input_filtration..=self.f).contains(&self.filtration(model, &word)?) {
                continue;
            }
            let mut lower = Vec::new();
            let mut leading = Vec::new();
            for term in odd_cobar::differential_profile_word(model, prime, self.profile, &word)? {
                match self.filtration(model, &term.0)?.cmp(&self.f) {
                    Ordering::Less => lower.push(term),
                    Ordering::Equal => leading.push(term),
                    Ordering::Greater => {}
                }
            }
            let lower = rows.encode(model, lower)?;
            let leading = word_indices(model, self.leading_words(), &leading)?;
            // When the lower terms cancel, the witness gives the leading
            // coordinates of a boundary in the requested filtration.
            if let Some(image) = lower_images.insert_with_witness(model, prime, lower, leading)? {
                boundaries.insert(model, prime, image)?;
            }
        }
        Ok(boundaries)
    }

    fn leading_words(&self) -> &[Word] {
        &self.words[self.leading_start..]
    }

    pub(crate) fn rank(&self) -> usize {
        if self.r == 0 {
            self.words.len()
        } else {
            self.representatives.len()
        }
    }

    pub(crate) fn project(&self, model: &mut Model, terms: &[Term]) -> Result<Vec<Coordinate>> {
        let mut leading = Vec::new();
        let mut below_cutoff = Vec::new();
        for (word, coefficient) in terms {
            if odd_cobar::word_metadata(word)? != (self.s, self.t) {
                return Err(Error::plain(
                    "odd-prime cobar chain has the wrong May bidegree",
                ));
            }
            let filtration = self.filtration(model, word)?;
            if filtration < self.f {
                return Err(Error::plain(
                    "odd-prime cobar chain lies below the requested May filtration",
                ));
            }
            if filtration == self.f {
                leading.push((word.clone(), *coefficient));
            }
            if filtration < self.cutoff {
                below_cutoff.push((word.clone(), *coefficient));
            }
        }
        let leading = word_indices(model, self.leading_words(), &leading)?;
        if self.r == 0 {
            return Ok(leading);
        }
        for (word, _) in self.differential(model, &below_cutoff)? {
            if self.filtration(model, &word)? < self.cutoff {
                return Err(Error::plain(
                    "odd-prime cobar chain does not satisfy the May page cycle condition",
                ));
            }
        }
        let (remainder, coordinates) = self.span.reduce_with_witness(model, self.prime, leading)?;
        if !remainder.is_empty() {
            return Err(Error::plain(
                "cobar cycle is missing from the odd-prime May page basis",
            ));
        }
        Ok(coordinates)
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        class: &[Coordinate],
    ) -> Result<Vec<Term>> {
        if self.r == 0 {
            model.spend(class.len())?;
            return Ok(decode(&self.words, class.to_vec()));
        }
        let mut support = Vec::new();
        for &(index, coefficient) in class {
            support = add_scaled(
                model,
                self.prime,
                &support,
                &self.representatives[index],
                coefficient,
            )?;
        }
        let mut terms = decode(&self.words, support);
        let cells = terms.iter().map(odd_cobar::term_cells).sum::<usize>();
        model.spend(cells * (1 + terms.len().max(1).ilog2() as usize))?;
        terms.sort_unstable();
        Ok(terms)
    }

    fn filtration(&self, model: &mut Model, word: &Word) -> Result<usize> {
        // Each generator's weight is bounded by its internal degree, and
        // the base-p digit sum of an exponent is bounded by the exponent.
        model.spend(word_cells(word))?;
        let p = self.prime.get() as usize;
        let mut weight = 0;
        for monomial in word.iter() {
            for (index, &exponent) in monomial.exponents.iter().enumerate() {
                let mut quotient = exponent;
                while quotient > 0 {
                    model.spend(1)?;
                    weight += (2 * index + 1) * (quotient % p);
                    quotient /= p;
                }
            }
            for &index in monomial.taus.iter() {
                weight += 2 * index + 1;
            }
        }
        self.t
            .checked_sub(weight)
            .ok_or_else(|| Error::plain("odd-prime May weight exceeds the internal degree"))
    }

    fn differential(&self, model: &mut Model, terms: &[Term]) -> Result<Vec<Term>> {
        odd_cobar::differential_profile_support(model, self.prime, self.profile, terms)
    }

    fn build_transition(&self, model: &mut Model, target: TypeId) -> Result<Basis> {
        let mut cycles = Basis::default();
        if let (Some(s), Some(f)) = (self.s.checked_sub(1), self.f.checked_sub(self.r)) {
            let ty = model.profile_pmay_type(self.prime, self.profile, self.r, s, self.t, f)?;
            let incoming = model.pmay(ty)?;
            for index in 0..incoming.rank() {
                let lift = incoming.representative(model, &[(index, 1)])?;
                let image = self.differential(model, &lift)?;
                let coordinates = self.project(model, &image)?;
                cycles.insert(model, self.prime, coordinates)?;
            }
        }

        let next = model.pmay(target)?;
        for index in 0..next.rank() {
            let lift = next.representative(model, &[(index, 1)])?;
            let earlier = self.project(model, &lift)?;
            if cycles
                .insert_with_witness(model, self.prime, earlier, vec![(index, 1)])?
                .is_some()
            {
                return Err(Error::plain(
                    "odd-prime May page transition has dependent cycle lifts",
                ));
            }
        }
        model.allocate(cycles.cells())?;
        Ok(cycles)
    }
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(odd_ext::zero());
    }
    let coordinates = model.pmay(ty)?.project(model, terms)?;
    odd_ext::persist(model, coordinates)
}

pub(crate) fn differential(
    model: &mut Model,
    source: TypeId,
    coordinates: &[Coordinate],
    target: TypeId,
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(odd_ext::zero());
    }
    let page = model.pmay(source)?;
    let lift = page.representative(model, coordinates)?;
    let image = page.differential(model, &lift)?;
    class(model, target, &image)
}

pub(crate) fn multiply(
    model: &mut Model,
    left_ty: TypeId,
    left: &[Coordinate],
    right_ty: TypeId,
    right: &[Coordinate],
    target: TypeId,
) -> Result<Value> {
    if left.is_empty() || right.is_empty() {
        return Ok(odd_ext::zero());
    }
    let page = model.pmay(left_ty)?;
    let left = page.representative(model, left)?;
    let right = model.pmay(right_ty)?.representative(model, right)?;
    let product = odd_cobar::concat_support(model, page.prime, &left, &right)?;
    class(model, target, &product)
}

pub(crate) fn next(
    model: &mut Model,
    source: TypeId,
    coordinates: &[Coordinate],
    target: TypeId,
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(odd_ext::zero());
    }
    let page = model.pmay(source)?;
    let transition = if let Some(transition) = page.transition.get() {
        transition
    } else {
        let transition = page.build_transition(model, target)?;
        page.transition.get_or_init(|| transition)
    };
    let (remainder, class) =
        transition.reduce_with_witness(model, page.prime, coordinates.to_vec())?;
    if !remainder.is_empty() {
        return Err(Error::plain(format!(
            "mayNext requires a cycle for d_{}",
            page.r
        )));
    }
    odd_ext::persist(model, class)
}

pub(crate) fn previous(
    model: &mut Model,
    source: TypeId,
    coordinates: &[Coordinate],
    target: TypeId,
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(odd_ext::zero());
    }
    let lift = model.pmay(source)?.representative(model, coordinates)?;
    class(model, target, &lift)
}

pub(crate) fn restrict(
    model: &mut Model,
    source: TypeId,
    coordinates: &[Coordinate],
    target: TypeId,
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(odd_ext::zero());
    }
    let (_, Profile::An(n), ..) = model.kind(target).pmay_grade().unwrap() else {
        unreachable!()
    };
    let page = model.pmay(source)?;
    let lift = page.representative(model, coordinates)?;
    let terms = odd_profile::restrict(model, page.prime, n, &lift)?;
    class(model, target, &terms)
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "odd-prime May calculation exceeds the storage limit",
        ));
    }
    Ok(())
}
