use crate::cobar;
use crate::error::{Error, Result};
use crate::map_cochain::{self, Support};
use crate::map_ext::{self, MapCochains, Term};
use crate::may::{self, FilteredComplex, FilteredPage};
use crate::model::{Model, TypeId, Value};

const MAX_BASIS_TERMS: usize = 32_768;

pub(crate) struct MapMayPage {
    source: TypeId,
    target: TypeId,
    s: usize,
    t: i64,
    page: FilteredPage<MappingComplex>,
}

struct MappingComplex {
    source: TypeId,
    target: TypeId,
    t: i64,
    cochains: Option<MapCochains>,
}

impl MappingComplex {
    fn new(model: &mut Model, source: TypeId, target: TypeId, t: i64) -> Result<Self> {
        Ok(Self {
            source,
            target,
            t,
            cochains: MapCochains::new(model, source, target, t)?,
        })
    }
}

impl FilteredComplex for MappingComplex {
    type Term = Term;

    fn basis(&mut self, model: &mut Model, s: usize, f: usize) -> Result<Vec<Term>> {
        // A term with s factors has May weight at least s. A filtration
        // bound therefore gives a lower bound on q + t - p before enumeration.
        let Some(min_degree) = f.checked_add(s) else {
            return Ok(Vec::new());
        };
        let Some(cochains) = &self.cochains else {
            return Ok(Vec::new());
        };
        cochains.basis(model, s, MAX_BASIS_TERMS, min_degree)
    }

    fn differential(&mut self, model: &mut Model, term: &Term) -> Result<Vec<Term>> {
        let Some(cochains) = &mut self.cochains else {
            return Ok(Vec::new());
        };
        cochains.differential(model, term)
    }

    fn differential_support(&mut self, model: &mut Model, terms: &[Term]) -> Result<Vec<Term>> {
        map_ext::differential_support(model, self.source, self.target, terms)
    }

    fn filtration(&self, term: &Term) -> Result<usize> {
        let degree = term.source.degree as i128 + i128::from(self.t) - term.target.degree as i128;
        let degree = usize::try_from(degree)
            .map_err(|_| Error::plain("map cobar word degree is out of range"))?;
        Ok(may::filtration(degree, &term.word))
    }

    fn cells(term: &Term) -> usize {
        term.cells()
    }
}

impl MapMayPage {
    pub(crate) fn build(
        model: &mut Model,
        source: TypeId,
        target: TypeId,
        r: usize,
        s: usize,
        t: i64,
        f: usize,
    ) -> Result<Self> {
        let mut complex = MappingComplex::new(model, source, target, t)?;
        let page = FilteredPage::build(model, &mut complex, r, s, f)?;
        Ok(Self {
            source,
            target,
            s,
            t,
            page,
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.page.rank()
    }

    pub(crate) fn project(&self, model: &mut Model, terms: &[Term]) -> Result<Vec<usize>> {
        for term in terms {
            model.spend(term.cells())?;
            let (s, degree) = cobar::word_metadata(&term.word)?;
            if s != self.s
                || map_cochain::internal_degree(term.source.degree, term.target.degree, degree)?
                    != self.t
            {
                return Err(Error::plain("map cobar cochain has the wrong May bidegree"));
            }
        }
        let mut complex = MappingComplex::new(model, self.source, self.target, self.t)?;
        self.page.project(model, &mut complex, terms)
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Vec<Term>> {
        self.page.representative(model, class)
    }
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let coordinates = model.map_may(ty)?.project(model, terms)?;
    model.allocate(coordinates.len())?;
    Ok(Value::Class(coordinates.into()))
}

pub(crate) fn representative(model: &mut Model, ty: TypeId, class: &[usize]) -> Result<Support> {
    if class.is_empty() {
        return Ok([].into());
    }
    let terms = model.map_may(ty)?.representative(model, class)?;
    map_cochain::persist(model, terms)
}

pub(crate) fn differential(
    model: &mut Model,
    ty: TypeId,
    coordinates: &[usize],
    output: TypeId,
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let page = model.map_may(ty)?;
    let lift = page.representative(model, coordinates)?;
    let image = map_ext::differential_support(model, page.source, page.target, &lift)?;
    class(model, output, &image)
}

pub(crate) fn compose(
    model: &mut Model,
    (outer_ty, outer): (TypeId, &[usize]),
    (inner_ty, inner): (TypeId, &[usize]),
    output: TypeId,
) -> Result<Value> {
    if outer.is_empty() || inner.is_empty() {
        return Ok(Value::Class([].into()));
    }
    let outer = model.map_may(outer_ty)?.representative(model, outer)?;
    let inner = model.map_may(inner_ty)?.representative(model, inner)?;
    let terms = map_ext::compose_terms(model, &outer, &inner)?;
    class(model, output, &terms)
}
