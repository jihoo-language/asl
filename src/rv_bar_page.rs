use std::rc::Rc;

use crate::cohomology;
use crate::error::{Error, Result};
use crate::may::{FilteredComplex, FilteredPage};
use crate::model::{Kind, Model, TypeId, Value};
use crate::page_transition::Transition;
use crate::rv_bar::{self, Support, Word};

const MAX_BASIS: usize = 32_768;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) chain: rv_bar::Grade,
    pub(crate) r: usize,
    pub(crate) length: i64,
}

impl Grade {
    fn span(self) -> usize {
        self.chain.end - self.chain.start
    }

    pub(crate) fn differential(self) -> Result<Self> {
        let length = i64::try_from(i128::from(self.length) - self.r as i128)
            .map_err(|_| Error::plain("rvBarD filtration is out of range"))?;
        Ok(Self {
            chain: self.chain.differential()?,
            length,
            ..self
        })
    }

    pub(crate) fn next(self) -> Result<Self> {
        let r = self
            .r
            .checked_add(1)
            .ok_or_else(|| Error::plain("rvBarNext page is out of range"))?;
        Ok(Self { r, ..self })
    }

    pub(crate) fn previous(self) -> Result<Self> {
        let r = self
            .r
            .checked_sub(1)
            .ok_or_else(|| Error::plain("rvBarPrev requires a positive page number"))?;
        Ok(Self { r, ..self })
    }

    fn incoming(self) -> Option<Self> {
        let length = i128::from(self.length) + self.r as i128;
        if !(0..=self.span() as i128).contains(&length) {
            return None;
        }
        let chain = rv_bar::Grade {
            b: self.chain.b.checked_add(1)?,
            ..self.chain
        };
        Some(Self {
            chain,
            length: length as i64,
            ..self
        })
    }
}

struct Complex {
    chain: rv_bar::Grade,
    span: usize,
}

impl Complex {
    fn new(grade: Grade) -> Self {
        Self {
            chain: grade.chain,
            span: grade.span(),
        }
    }
}

impl FilteredComplex for Complex {
    type Term = Word;
    const NAME: &'static str = "bar";
    const ELEMENT_NAME: &'static str = "bar";

    fn basis(&mut self, model: &mut Model, s: usize, f: usize) -> Result<Vec<Word>> {
        let Some(length) = self.span.checked_sub(f) else {
            return Ok(Vec::new());
        };
        // The generic reducer asks for its current degree (1) and its
        // incoming cohomological degree (0). The latter is b+1 here.
        let b = match s {
            0 => match self.chain.b.checked_add(1) {
                Some(b) => b,
                None => return Ok(Vec::new()),
            },
            1 => self.chain.b,
            _ => unreachable!(),
        };
        let basis = rv_bar::basis(model, rv_bar::Grade { b, ..self.chain }, length, MAX_BASIS)?;
        model.allocate(basis.iter().map(|word| rv_bar::word_cells(word)).sum())?;
        Ok(basis)
    }

    fn differential(&mut self, model: &mut Model, term: &Word) -> Result<Vec<Word>> {
        self.differential_support(model, std::slice::from_ref(term))
    }

    fn differential_support(&mut self, model: &mut Model, terms: &[Word]) -> Result<Vec<Word>> {
        rv_bar::differential(model, self.chain, terms)
    }

    fn filtration(&self, term: &Word) -> Result<usize> {
        Ok(self.span - term.len())
    }

    fn cells(term: &Word) -> usize {
        rv_bar::word_cells(term)
    }
}

pub(crate) struct Page {
    grade: Grade,
    filtered: Option<FilteredPage<Complex>>,
}

impl Page {
    pub(crate) fn build(model: &mut Model, grade: Grade) -> Result<Self> {
        let mut complex = Complex::new(grade);
        let filtered = usize::try_from(grade.length)
            .ok()
            .and_then(|l| grade.span().checked_sub(l))
            .map(|f| FilteredPage::build(model, &mut complex, grade.r, 1, f))
            .transpose()?;
        Ok(Self { grade, filtered })
    }

    pub(crate) fn rank(&self) -> usize {
        self.filtered.as_ref().map_or(0, FilteredPage::rank)
    }

    fn project(&self, model: &mut Model, terms: &[Word]) -> Result<Vec<usize>> {
        model.spend(terms.len())?;
        if terms
            .iter()
            .any(|word| word.len() as i128 > i128::from(self.grade.length))
        {
            return Err(Error::plain(
                "bar chain lies above the requested bar page length filtration",
            ));
        }
        if let Some(page) = &self.filtered {
            return page.project(model, &mut Complex::new(self.grade), terms);
        }
        if self.grade.r != 0 && !terms.is_empty() {
            let cutoff = i128::from(self.grade.length) - self.grade.r as i128;
            let image = rv_bar::differential(model, self.grade.chain, terms)?;
            if image.iter().any(|word| word.len() as i128 > cutoff) {
                return Err(Error::plain(
                    "bar chain does not satisfy the bar page cycle condition",
                ));
            }
        }
        Ok(Vec::new())
    }

    fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<Word>> {
        if class.is_empty() {
            return Ok(Vec::new());
        }
        self.filtered.as_ref().unwrap().representative(model, class)
    }
}

pub(crate) fn representative(model: &mut Model, ty: TypeId, class: &[usize]) -> Result<Support> {
    if class.is_empty() {
        return Ok([].into());
    }
    let terms = model.rv_bar_page(ty)?.representative(model, class)?;
    rv_bar::persist(model, terms)
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Word]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(cohomology::zero());
    }
    let coordinates = model.rv_bar_page(ty)?.project(model, terms)?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

pub(crate) fn differential(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    coordinates: &[usize],
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(cohomology::zero());
    }
    let page = model.rv_bar_page(input)?;
    let lift = page.representative(model, coordinates)?;
    let image = rv_bar::differential(model, page.grade.chain, &lift)?;
    class(model, output, &image)
}

pub(crate) fn previous(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    coordinates: &[usize],
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(cohomology::zero());
    }
    let lift = model
        .rv_bar_page(input)?
        .representative(model, coordinates)?;
    class(model, output, &lift)
}

pub(crate) fn transition(model: &mut Model, ty: TypeId) -> Result<Transition> {
    let Kind::RVBarE(grade) = *model.kind(ty) else {
        unreachable!();
    };
    let source = model.rv_bar_page(ty)?;
    let next = model.rv_bar_page_type(grade.next()?)?;
    let next = model.rv_bar_page(next)?;
    let incoming: Option<Rc<Page>> = grade
        .incoming()
        .map(|incoming| {
            let ty = model.rv_bar_page_type(incoming)?;
            model.rv_bar_page(ty)
        })
        .transpose()?;
    Transition::from_columns(
        model,
        grade.r,
        "rvBarNext",
        incoming.as_ref().map_or(0, |page| page.rank()),
        next.rank(),
        |model, index| {
            let incoming = incoming.as_ref().unwrap();
            let lift = incoming.representative(model, &[index])?;
            let image = rv_bar::differential(model, incoming.grade.chain, &lift)?;
            source.project(model, &image)
        },
        |model, index| {
            let lift = next.representative(model, &[index])?;
            source.project(model, &lift)
        },
    )
}
