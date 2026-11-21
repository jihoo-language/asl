use crate::cobar::{self, Word};
use crate::error::{Error, Result};
use crate::may::{FilteredComplex, FilteredPage};
use crate::model::{Model, TypeId, Value};
use crate::profile::Profile;

const MAX_BASIS_WORDS: usize = 32_768;

struct CentralExtension {
    n: usize,
    t: usize,
}

impl CentralExtension {
    fn profile(&self) -> Profile {
        Profile::An(self.n)
    }

    fn can_lead(&self, s: usize, f: usize) -> bool {
        let Some(q) = s.checked_sub(f) else {
            return false;
        };
        let Some(excess) = self.t.checked_sub(s) else {
            return false;
        };
        if self.n == 0 {
            // A(0) has only the top generator among positive-degree monomials.
            return f == 0 && excess == 0;
        }
        if q == 0 {
            return s != 0 || self.t == 0;
        }
        if self.n >= usize::BITS as usize {
            return false;
        }
        let top_degree = usize::MAX >> (usize::BITS as usize - 1 - self.n);
        q <= excess / (top_degree - 1)
    }
}

impl FilteredComplex for CentralExtension {
    type Term = Word;
    const NAME: &'static str = "Cartan–Eilenberg";

    fn basis(&mut self, model: &mut Model, s: usize, f: usize) -> Result<Vec<Word>> {
        if f > s || (self.n == 0 && f != 0) {
            return Ok(Vec::new());
        }
        self.profile().basis(model, s, self.t, MAX_BASIS_WORDS)
    }

    fn differential(&mut self, model: &mut Model, word: &Word) -> Result<Vec<Word>> {
        self.profile().differential_word(model, word)
    }

    fn differential_support(&mut self, model: &mut Model, words: &[Word]) -> Result<Vec<Word>> {
        self.profile().differential_support(model, words)
    }

    fn filtration(&self, word: &Word) -> Result<usize> {
        // The top generator is exterior. Exactly the factors missing it
        // contribute to the base filtration s-q.
        Ok(word
            .iter()
            .filter(|factor| factor.get(self.n).is_none_or(|&power| power == 0))
            .count())
    }

    fn cells(word: &Word) -> usize {
        cobar::word_cells(word)
    }
}

pub(crate) struct CePage {
    n: usize,
    s: usize,
    t: usize,
    page: FilteredPage<CentralExtension>,
}

impl CePage {
    pub(crate) fn build(
        model: &mut Model,
        n: usize,
        r: usize,
        s: usize,
        t: usize,
        f: usize,
    ) -> Result<Self> {
        let mut complex = CentralExtension { n, t };
        let page = if complex.can_lead(s, f) {
            FilteredPage::build(model, &mut complex, r, s, f)?
        } else {
            FilteredPage::empty(r, f)
        };
        Ok(Self { n, s, t, page })
    }

    pub(crate) fn rank(&self) -> usize {
        self.page.rank()
    }

    pub(crate) fn project(&self, model: &mut Model, words: &[Word]) -> Result<Vec<usize>> {
        for word in words {
            model.spend(cobar::word_cells(word))?;
            if cobar::word_metadata(word)? != (self.s, self.t) {
                return Err(Error::plain(
                    "cobar chain has the wrong Cartan–Eilenberg bidegree",
                ));
            }
        }
        let mut complex = CentralExtension {
            n: self.n,
            t: self.t,
        };
        self.page.project(model, &mut complex, words)
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<Word>> {
        self.page.representative(model, class)
    }
}

pub(crate) fn class(model: &mut Model, ty: TypeId, words: &[Word]) -> Result<Value> {
    if words.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let coordinates = model.ce_page(ty)?.project(model, words)?;
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
    let page = model.ce_page(source)?;
    let lift = page.representative(model, coordinates)?;
    let image = Profile::An(page.n).differential_support(model, &lift)?;
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
    let left = model.ce_page(left_ty)?.representative(model, left)?;
    let right = model.ce_page(right_ty)?.representative(model, right)?;
    let words = cobar::concat_temporary(model, &left, &right)?;
    class(model, target, &words)
}
