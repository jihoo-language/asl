use std::cell::OnceCell;
use std::rc::Rc;

use crate::binary_field::BinaryField;
use crate::binary_linear::{add_multiple, check_storage, Basis, Entry, Vector};
use crate::error::{Error, Result};
use crate::jet_hom::{self, Complex, Term, WeightedTerm};
use crate::may::Rows;
use crate::model::{Model, TypeId, Value};
use crate::rank_jets::Path;

const MAX_BASIS: usize = 32_768;

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct Grade {
    pub(crate) chain: jet_hom::Grade,
    pub(crate) r: usize,
    pub(crate) f: usize,
}

impl Grade {
    pub(crate) fn differential(self) -> Result<Self> {
        let f = self
            .f
            .checked_add(self.r)
            .ok_or_else(|| Error::plain("jetD filtration is out of range"))?;
        Ok(Self {
            chain: self.chain.differential(),
            f,
            ..self
        })
    }

    pub(crate) fn next(self) -> Result<Self> {
        let r = self
            .r
            .checked_add(1)
            .ok_or_else(|| Error::plain("jetNext page is out of range"))?;
        Ok(Self { r, ..self })
    }

    pub(crate) fn prev(self) -> Result<Self> {
        let r = self
            .r
            .checked_sub(1)
            .ok_or_else(|| Error::plain("jetPrev requires a positive page number"))?;
        Ok(Self { r, ..self })
    }

    pub(crate) fn compose(outer: Self, inner: Self) -> Result<Self> {
        if outer.r != inner.r {
            return Err(Error::plain(
                "ecomp requires JetPage classes on the same page",
            ));
        }
        let f = outer
            .f
            .checked_add(inner.f)
            .ok_or_else(|| Error::plain("ecomp filtration is out of range"))?;
        Ok(Self {
            chain: jet_hom::Grade::compose(outer.chain, inner.chain)?,
            r: outer.r,
            f,
        })
    }
}

pub(crate) struct Page {
    grade: Grade,
    cutoff: usize,
    leading_only: bool,
    basis: Vec<Term>,
    leading_start: usize,
    span: Basis,
    lifts: Vec<Box<[Entry]>>,
    transition: OnceCell<Basis>,
}

impl Page {
    pub(crate) fn build(model: &mut Model, grade: Grade) -> Result<Self> {
        let path = model.jet_path(grade.chain.diagram.path);
        let mut page = Self {
            grade,
            cutoff: grade.f.saturating_add(grade.r).min(path.order),
            basis: Vec::new(),
            leading_only: true,
            leading_start: 0,
            span: Basis::default(),
            lifts: Vec::new(),
            transition: OnceCell::new(),
        };
        if grade.f >= path.order {
            return Ok(page);
        }
        let complex = model.jet_complex(grade.chain.diagram)?;
        page.leading_only = grade.r == 0 || complex.is_zero();
        let end = if page.leading_only {
            grade.f + 1
        } else {
            page.cutoff
        };
        let mut basis = complex.basis(model, grade.chain.parity, grade.f..end, MAX_BASIS)?;
        let leading_len = basis.partition_point(|term| term.power == grade.f);
        if leading_len == 0 {
            return Ok(page);
        }
        model.spend(basis.len())?;
        model.allocate(5 * basis.len())?;
        basis.rotate_left(leading_len);
        page.leading_start = basis.len() - leading_len;
        page.basis = basis;
        if !page.leading_only {
            page.build_quotient(model, &path, &complex)?;
        }
        Ok(page)
    }

    fn build_quotient(&mut self, model: &mut Model, path: &Path, complex: &Complex) -> Result<()> {
        let field = path.field;
        let mut span = self.boundaries(model, path, complex)?;
        let mut cycles = Basis::default();
        let mut rows = Rows::new("jet page");
        let mut lift_cells = 0;
        for (index, &term) in self.basis.iter().enumerate() {
            let image = complex.column(model, path, term)?;
            let lower = encode(
                model,
                &mut rows,
                image
                    .into_iter()
                    .filter(|(term, _)| term.power < self.cutoff),
            )?;
            // Process correction terms first so their dependencies have zero leading part.
            let Some(lift) = cycles.insert_with_witness(model, field, lower, vec![(index, 1)])?
            else {
                continue;
            };
            if index < self.leading_start {
                continue;
            }
            let leading = lift
                .iter()
                .filter_map(|&(index, coefficient)| {
                    index
                        .checked_sub(self.leading_start)
                        .map(|i| (i, coefficient))
                })
                .collect();
            if span
                .insert_with_witness(model, field, leading, vec![(self.lifts.len(), 1)])?
                .is_none()
            {
                lift_cells += 1 + 2 * lift.len();
                check_storage(lift_cells + span.cells())?;
                self.lifts.push(lift.into_boxed_slice());
            }
        }
        model.allocate(lift_cells + span.cells())?;
        self.span = span;
        Ok(())
    }

    fn boundaries(&self, model: &mut Model, path: &Path, complex: &Complex) -> Result<Basis> {
        let r = self.grade.r;
        let f = self.grade.f;
        let lo = f.saturating_sub(r - 1);
        let mut rows = Rows::new("jet page");
        let mut below = Basis::default();
        let mut boundaries = Basis::default();
        let previous = complex.basis(model, self.grade.chain.parity ^ 1, lo..f + 1, MAX_BASIS)?;
        for term in previous {
            let image = complex.column(model, path, term)?;
            let lower = encode(
                model,
                &mut rows,
                image.iter().copied().filter(|(term, _)| term.power < f),
            )?;
            let leading =
                self.leading_column(model, image.into_iter().filter(|(term, _)| term.power == f))?;
            if let Some(boundary) = below.insert_with_witness(model, path.field, lower, leading)? {
                boundaries.insert(model, path.field, boundary)?;
            }
        }
        Ok(boundaries)
    }

    fn leading_column(
        &self,
        model: &mut Model,
        entries: impl Iterator<Item = WeightedTerm>,
    ) -> Result<Vector> {
        let basis = &self.basis[self.leading_start..];
        let mut result = Vec::new();
        for (term, coefficient) in entries {
            model.spend(5 * (1 + basis.len().max(1).ilog2() as usize))?;
            let index = basis
                .binary_search(&term)
                .map_err(|_| Error::plain("term is outside the jet page filtration slice"))?;
            result.push((index, coefficient));
        }
        sort_coordinates(model, &mut result)?;
        Ok(result)
    }

    pub(crate) fn rank(&self) -> usize {
        if self.leading_only {
            self.basis.len()
        } else {
            self.lifts.len()
        }
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        path: &Path,
        complex: &Complex,
        terms: &[WeightedTerm],
    ) -> Result<Vector> {
        model.spend(terms.len())?;
        if terms
            .first()
            .is_some_and(|(term, _)| term.power < self.grade.f)
        {
            return Err(Error::plain(
                "cochain lies below the requested jet filtration",
            ));
        }
        let leading = self.leading_column(
            model,
            terms
                .iter()
                .copied()
                .take_while(|(term, _)| term.power == self.grade.f),
        )?;
        if self.leading_only {
            return Ok(leading);
        }
        let end = terms.partition_point(|(term, _)| term.power < self.cutoff);
        let image = complex.differential(model, path, &terms[..end])?;
        if image
            .first()
            .is_some_and(|(term, _)| term.power < self.cutoff)
        {
            return Err(Error::plain(
                "cochain does not satisfy the jet page cycle condition",
            ));
        }
        let mut coordinates = Vec::new();
        let remainder = self
            .span
            .reduce(model, path.field, leading, Some(&mut coordinates))?;
        if !remainder.is_empty() {
            return Err(Error::plain(
                "cochain cycle is missing from the jet page basis",
            ));
        }
        Ok(coordinates)
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        field: BinaryField,
        class: &[Entry],
    ) -> Result<Vec<WeightedTerm>> {
        let mut terms = if self.leading_only {
            model.spend(class.len())?;
            self.decode(class)
        } else {
            let mut support = Vec::new();
            for &(index, coefficient) in class {
                support = add_multiple(model, field, &support, &self.lifts[index], coefficient)?;
            }
            self.decode(&support)
        };
        model.spend(5 * terms.len() * (1 + terms.len().max(1).ilog2() as usize))?;
        terms.sort_unstable();
        Ok(terms)
    }

    fn decode(&self, column: &[Entry]) -> Vec<WeightedTerm> {
        column
            .iter()
            .map(|&(index, coefficient)| (self.basis[index], coefficient))
            .collect()
    }

    fn build_transition(
        &self,
        model: &mut Model,
        path: &Path,
        complex: &Complex,
        later: TypeId,
    ) -> Result<Basis> {
        let mut cycles = Basis::default();
        if let Some(f) = self.grade.f.checked_sub(self.grade.r) {
            let incoming = Grade {
                chain: self.grade.chain.differential(),
                f,
                ..self.grade
            };
            let ty = model.jet_page_type(incoming)?;
            let page = model.jet_page(ty)?;
            for index in 0..page.rank() {
                let lift = page.representative(model, path.field, &[(index, 1)])?;
                let image = complex.differential(model, path, &lift)?;
                let coordinates = self.project(model, path, complex, &image)?;
                cycles.insert(model, path.field, coordinates)?;
            }
        }

        let next = model.jet_page(later)?;
        for index in 0..next.rank() {
            let lift = next.representative(model, path.field, &[(index, 1)])?;
            let old_class = self.project(model, path, complex, &lift)?;
            if cycles
                .insert_with_witness(model, path.field, old_class, vec![(index, 1)])?
                .is_some()
            {
                return Err(Error::plain(
                    "jet page transition has dependent next-page lifts",
                ));
            }
        }
        model.allocate(cycles.cells())?;
        Ok(cycles)
    }
}

fn encode(
    model: &mut Model,
    rows: &mut Rows<Term>,
    terms: impl Iterator<Item = WeightedTerm>,
) -> Result<Vector> {
    let mut coordinates = Vec::new();
    for (term, coefficient) in terms {
        coordinates.push((rows.coordinate(model, 5, term)?, coefficient));
    }
    sort_coordinates(model, &mut coordinates)?;
    Ok(coordinates)
}

fn sort_coordinates(model: &mut Model, coordinates: &mut [Entry]) -> Result<()> {
    model.spend(coordinates.len() * (1 + coordinates.len().max(1).ilog2() as usize))?;
    coordinates.sort_unstable();
    Ok(())
}

pub(crate) fn zero() -> Value {
    Value::JetClass([].into())
}

pub(crate) fn persist(model: &mut Model, values: Vector) -> Result<Value> {
    model.allocate(2 * values.len())?;
    Ok(Value::JetClass(values.into()))
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[WeightedTerm]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(zero());
    }
    let page = model.jet_page(ty)?;
    let path = model.jet_path(page.grade.chain.diagram.path);
    let complex = model.jet_complex(page.grade.chain.diagram)?;
    let class = page.project(model, &path, &complex, terms)?;
    persist(model, class)
}

pub(crate) fn representative(
    model: &mut Model,
    ty: TypeId,
    class: &[Entry],
) -> Result<Vec<WeightedTerm>> {
    if class.is_empty() {
        return Ok(Vec::new());
    }
    let page = model.jet_page(ty)?;
    let path = model.jet_path(page.grade.chain.diagram.path);
    page.representative(model, path.field, class)
}

pub(crate) fn differential(
    model: &mut Model,
    ty: TypeId,
    target: TypeId,
    class: &[Entry],
) -> Result<Value> {
    if class.is_empty() {
        return Ok(zero());
    }
    let page = model.jet_page(ty)?;
    let path = model.jet_path(page.grade.chain.diagram.path);
    let complex = model.jet_complex(page.grade.chain.diagram)?;
    let lift = page.representative(model, path.field, class)?;
    let image = complex.differential(model, &path, &lift)?;
    self::class(model, target, &image)
}

pub(crate) fn compose(
    model: &mut Model,
    outer: (TypeId, &[Entry]),
    inner: (TypeId, &[Entry]),
    target: TypeId,
) -> Result<Value> {
    if outer.1.is_empty() || inner.1.is_empty() {
        return Ok(zero());
    }
    let outer_page = model.jet_page(outer.0)?;
    let inner_page = model.jet_page(inner.0)?;
    let path = model.jet_path(outer_page.grade.chain.diagram.path);
    let outer = outer_page.representative(model, path.field, outer.1)?;
    let inner = inner_page.representative(model, path.field, inner.1)?;
    let terms = jet_hom::compose(model, &path, &outer, &inner)?;
    class(model, target, &terms)
}

pub(crate) fn next(
    model: &mut Model,
    ty: TypeId,
    target: TypeId,
    class: &[Entry],
) -> Result<Value> {
    if class.is_empty() {
        return Ok(zero());
    }
    let page = model.jet_page(ty)?;
    let path = model.jet_path(page.grade.chain.diagram.path);
    let complex = model.jet_complex(page.grade.chain.diagram)?;
    let transition = match page.transition.get() {
        Some(transition) => transition,
        None => {
            let transition = page.build_transition(model, &path, &complex, target)?;
            page.transition.get_or_init(|| transition)
        }
    };
    let mut coordinates = Vec::new();
    let remainder = transition.reduce(model, path.field, class.to_vec(), Some(&mut coordinates))?;
    if !remainder.is_empty() {
        return Err(Error::plain("jetNext requires a cycle on the current page"));
    }
    persist(model, coordinates)
}

pub(crate) fn prev(
    model: &mut Model,
    ty: TypeId,
    target: TypeId,
    values: &[Entry],
) -> Result<Value> {
    let lift = representative(model, ty, values)?;
    class(model, target, &lift)
}

pub(crate) fn elements(
    model: &mut Model,
    field: BinaryField,
    rank: usize,
    limit: usize,
) -> Result<Vec<Value>> {
    let max_rank = limit.ilog2() as usize / field.degree();
    if rank > max_rank {
        return Err(Error::plain(format!(
            "type exceeds the enumeration limit of {limit} elements"
        )));
    }
    let order = field.order();
    let count = order.pow(rank as u32);
    model.allocate(count + 2 * rank * (count - count / order))?;
    model.spend(count * (1 + rank))?;
    Ok((0..count)
        .map(|mut digits| {
            Value::JetClass(
                (0..rank)
                    .filter_map(|coordinate| {
                        let coefficient = digits % order;
                        digits /= order;
                        (coefficient != 0).then_some((coordinate, coefficient))
                    })
                    .collect(),
            )
        })
        .collect())
}

pub(crate) type Coordinates = Rc<[Entry]>;
