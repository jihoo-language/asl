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

#[derive(Clone, Copy)]
pub(crate) enum Filtration {
    May(Profile),
    Ce(usize),
}

impl Filtration {
    fn profile(self) -> Profile {
        match self {
            Self::May(profile) => profile,
            Self::Ce(n) => Profile::An(n),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::May(_) => "May",
            Self::Ce(_) => "Cartan–Eilenberg",
        }
    }
}

pub(crate) struct OddPage {
    prime: Prime,
    filter: Filtration,
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

impl OddPage {
    pub(crate) fn build(
        model: &mut Model,
        prime: Prime,
        filter: Filtration,
        r: usize,
        s: usize,
        t: usize,
        f: usize,
    ) -> Result<Self> {
        let mut page = Self {
            prime,
            filter,
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
        if !page.can_lead() {
            return Ok(page);
        }
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

    fn can_lead(&self) -> bool {
        let Filtration::Ce(n) = self.filter else {
            return true;
        };
        let (Some(q), Some(excess)) = (self.s.checked_sub(self.f), self.t.checked_sub(self.s))
        else {
            return false;
        };
        if n == 0 {
            return self.f == 0 && excess == 0;
        }
        if q == 0 {
            return self.s > 0 || self.t == 0;
        }
        // Each of the q factors containing tau_n contributes at least
        // |tau_n|-1 more internal degree than the minimum per factor.
        odd_cobar::xi_degree(self.prime, n).is_ok_and(|extra| q <= excess / extra)
    }

    fn basis(&self, model: &mut Model, s: usize, f: usize) -> Result<Vec<Word>> {
        let Some(excess) = self.t.checked_sub(s) else {
            return Ok(Vec::new());
        };
        let max_filtration = match self.filter {
            Filtration::May(_) => excess,
            Filtration::Ce(_) => s,
        };
        if f > max_filtration {
            return Ok(Vec::new());
        }
        odd_cobar::basis_profile(
            model,
            self.prime,
            self.filter.profile(),
            s,
            self.t,
            MAX_BASIS_WORDS,
        )
    }

    fn build_quotient(&mut self, model: &mut Model) -> Result<()> {
        let prime = self.prime;
        let mut span = self.boundaries(model)?;
        let mut rows = WordRows::default();
        let mut cycles = Basis::default();
        let mut cells = 0;
        for (index, word) in self.words.iter().enumerate() {
            let mut lower = Vec::new();
            for term in
                odd_cobar::differential_profile_word(model, prime, self.filter.profile(), word)?
            {
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
            for term in
                odd_cobar::differential_profile_word(model, prime, self.filter.profile(), &word)?
            {
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
                return Err(Error::plain(format!(
                    "odd-prime cobar chain has the wrong {} bidegree",
                    self.filter.name()
                )));
            }
            let filtration = self.filtration(model, word)?;
            if filtration < self.f {
                return Err(Error::plain(format!(
                    "odd-prime cobar chain lies below the requested {} filtration",
                    self.filter.name()
                )));
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
                return Err(Error::plain(format!(
                    "odd-prime cobar chain does not satisfy the {} page cycle condition",
                    self.filter.name()
                )));
            }
        }
        let (remainder, coordinates) = self.span.reduce_with_witness(model, self.prime, leading)?;
        if !remainder.is_empty() {
            return Err(Error::plain(format!(
                "cobar cycle is missing from the odd-prime {} page basis",
                self.filter.name()
            )));
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
        model.spend(word_cells(word))?;
        if let Filtration::Ce(n) = self.filter {
            // The top exterior generator occurs at most once in each factor.
            return Ok(word
                .iter()
                .filter(|factor| factor.taus.last().copied() != Some(n))
                .count());
        }
        // Each generator's weight is bounded by its internal degree, and
        // the base-p digit sum of an exponent is bounded by the exponent.
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
        odd_cobar::differential_profile_support(model, self.prime, self.filter.profile(), terms)
    }

    fn build_transition(&self, model: &mut Model, target: TypeId) -> Result<Basis> {
        let mut cycles = Basis::default();
        if let (Some(s), Some(f)) = (self.s.checked_sub(1), self.f.checked_sub(self.r)) {
            let ty = match self.filter {
                Filtration::May(profile) => {
                    model.profile_pmay_type(self.prime, profile, self.r, s, self.t, f)?
                }
                Filtration::Ce(n) => model.pce_type(self.prime, n, self.r, s, self.t, f)?,
            };
            let incoming = model.odd_page(ty)?;
            for index in 0..incoming.rank() {
                let lift = incoming.representative(model, &[(index, 1)])?;
                let image = self.differential(model, &lift)?;
                let coordinates = self.project(model, &image)?;
                cycles.insert(model, self.prime, coordinates)?;
            }
        }

        let next = model.odd_page(target)?;
        for index in 0..next.rank() {
            let lift = next.representative(model, &[(index, 1)])?;
            let earlier = self.project(model, &lift)?;
            if cycles
                .insert_with_witness(model, self.prime, earlier, vec![(index, 1)])?
                .is_some()
            {
                return Err(Error::plain(format!(
                    "odd-prime {} page transition has dependent cycle lifts",
                    self.filter.name()
                )));
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
    let coordinates = model.odd_page(ty)?.project(model, terms)?;
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
    let page = model.odd_page(source)?;
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
    let page = model.odd_page(left_ty)?;
    let left = page.representative(model, left)?;
    let right = model.odd_page(right_ty)?.representative(model, right)?;
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
    let page = model.odd_page(source)?;
    let transition = if let Some(transition) = page.transition.get() {
        transition
    } else {
        let transition = page.build_transition(model, target)?;
        page.transition.get_or_init(|| transition)
    };
    let (remainder, class) =
        transition.reduce_with_witness(model, page.prime, coordinates.to_vec())?;
    if !remainder.is_empty() {
        let name = match page.filter {
            Filtration::May(_) => "mayNext",
            Filtration::Ce(_) => "ceNext",
        };
        return Err(Error::plain(format!(
            "{name} requires a cycle for d_{}",
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
    let lift = model.odd_page(source)?.representative(model, coordinates)?;
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
    let page = model.odd_page(source)?;
    let lift = page.representative(model, coordinates)?;
    let terms = odd_profile::restrict(model, page.prime, n, &lift)?;
    class(model, target, &terms)
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "odd-prime filtered page calculation exceeds the storage limit",
        ));
    }
    Ok(())
}
