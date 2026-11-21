use std::cell::OnceCell;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt;

use crate::bp::{self, Shape};
use crate::bp_ext::{self, Domain};
use crate::bp_polynomial::{Monomial, Polynomial};
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::odd_cobar::Prime;
use crate::odd_ext;
use crate::prime_linear::{add_scaled, Basis, Coordinate};

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) shape: Shape,
    pub(crate) page: usize,
    pub(crate) t: usize,
    pub(crate) filtration: usize,
}

impl fmt::Display for Grade {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        let shape = self.shape;
        write!(
            out,
            "BPVPage({}, {}, {}, {}, {}, {}, {})",
            shape.prime,
            shape.height,
            self.page,
            shape.arrows,
            shape.depth,
            self.t,
            self.filtration
        )
    }
}

impl Grade {
    pub(crate) fn new(shape: Shape, page: usize, t: usize, filtration: usize) -> Result<Self> {
        if shape.depth == 0 || shape.depth > shape.height || shape.thickness.is_some() {
            return Err(Error::plain(
                "BPVPage requires a BPI quotient with an invariant generator v_h (1 <= h <= height)",
            ));
        }
        Ok(Self {
            shape,
            page,
            t,
            filtration,
        })
    }

    pub(crate) fn prime(self) -> Prime {
        Prime::validated(self.shape.prime)
    }

    fn exponent(self, monomial: &Monomial) -> usize {
        monomial
            .first()
            .filter(|&&(index, _)| index == self.shape.depth - 1)
            .map_or(0, |&(_, power)| power)
    }

    fn shifted_t(self, exponent: usize) -> Option<usize> {
        if exponent == 0 {
            return Some(self.t);
        }
        let power = u128::from(self.shape.prime).checked_pow(self.shape.depth as u32)?;
        let degree = (power - 1).checked_mul(2)?;
        let weight = degree.checked_mul(exponent as u128)?;
        self.t.checked_sub(usize::try_from(weight).ok()?)
    }

    fn cutoff(self) -> usize {
        self.filtration.saturating_add(self.page.max(1))
    }

    // Exponents in [first, first + width), keeping source variables divided
    // by v_h^first during enumeration. All arrow slots are still normalized.
    fn basis(self, model: &mut Model, first: usize, width: usize) -> Result<Vec<Monomial>> {
        let Some(t) = self.shifted_t(first) else {
            return Ok(Vec::new());
        };
        let mut monomials = bp_ext::basis(model, self.shape.thick(width)?, t)?;
        if first != 0 {
            for monomial in &mut monomials {
                model.spend(1 + monomial.len())?;
                if self.exponent(monomial) == 0 {
                    *monomial = std::iter::once((self.shape.depth - 1, first))
                        .chain(monomial.iter().copied())
                        .collect();
                } else {
                    monomial[0].1 += first;
                }
            }
            model.spend(monomials.len() * (1 + monomials.len().max(1).ilog2() as usize))?;
            monomials.sort_unstable();
        }
        Ok(monomials)
    }

    fn truncated_d(
        self,
        model: &mut Model,
        input: &Polynomial,
        first: usize,
        cutoff: usize,
    ) -> Result<Polynomial> {
        debug_assert!(cutoff > 0);
        let Some(t) = self.shifted_t(first) else {
            return Ok(Polynomial::default());
        };
        let shape = self.shape.thick(cutoff)?;
        if !Domain::absolute(shape.next()?).has_monomials(t) {
            return Ok(Polynomial::default());
        }
        bp::differential(model, shape, input)
    }

    pub(crate) fn cup(self, other: Self) -> Result<Self> {
        if self.page != other.page {
            return Err(Error::plain("bpCup requires the same BPVPage index"));
        }
        Ok(Self {
            shape: self.shape.cup(other.shape)?,
            t: self
                .t
                .checked_add(other.t)
                .ok_or_else(|| Error::plain("BPVPage cup internal degree is out of range"))?,
            filtration: self
                .filtration
                .checked_add(other.filtration)
                .ok_or_else(filtration_overflow)?,
            ..self
        })
    }
}

fn filtration_overflow() -> Error {
    Error::plain("BPVPage filtration is out of range")
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Rep,
    Differential,
    Next,
    Previous,
}

impl Operation {
    pub(crate) fn from_name(name: &str) -> Option<Self> {
        match name {
            "bpVRep" => Some(Self::Rep),
            "bpVD" => Some(Self::Differential),
            "bpVNext" => Some(Self::Next),
            "bpVPrev" => Some(Self::Previous),
            _ => None,
        }
    }

    pub(crate) fn target(self, grade: Grade) -> Result<Grade> {
        Ok(match self {
            Self::Rep => grade,
            Self::Differential => Grade {
                shape: grade.shape.next()?,
                filtration: grade
                    .filtration
                    .checked_add(grade.page)
                    .ok_or_else(filtration_overflow)?,
                ..grade
            },
            Self::Next => Grade {
                page: grade
                    .page
                    .checked_add(1)
                    .ok_or_else(|| Error::plain("BPVPage index is out of range"))?,
                ..grade
            },
            Self::Previous => Grade {
                page: grade
                    .page
                    .checked_sub(1)
                    .ok_or_else(|| Error::plain("BPVPage zero has no previous page"))?,
                ..grade
            },
        })
    }
}

#[derive(Default)]
struct Rows {
    index: BTreeMap<Monomial, usize>,
    cells: usize,
}

impl Rows {
    fn add(
        &mut self,
        model: &mut Model,
        monomial: &Monomial,
        coefficient: i128,
        column: &mut Vec<Coordinate>,
    ) -> Result<()> {
        model.spend((1 + monomial.len()) * (1 + self.index.len().max(1).ilog2() as usize))?;
        let index = if let Some(&index) = self.index.get(monomial) {
            index
        } else {
            self.cells += 4 + 2 * monomial.len();
            check_cells(self.cells)?;
            let index = self.index.len();
            self.index.insert(monomial.clone(), index);
            index
        };
        column.push((index, coefficient as u32));
        Ok(())
    }
}

fn sort_column(model: &mut Model, column: &mut [Coordinate]) -> Result<()> {
    model.spend(column.len() * (1 + column.len().max(1).ilog2() as usize))?;
    column.sort_unstable();
    Ok(())
}

fn check_cells(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "BPVPage reduction exceeds the row storage limit",
        ));
    }
    Ok(())
}

pub(crate) struct Page {
    grade: Grade,
    words: Vec<Monomial>,
    leading_start: usize,
    span: Basis,
    lifts: Vec<Box<[Coordinate]>>,
    transition: OnceCell<Basis>,
}

impl Page {
    pub(crate) fn build(model: &mut Model, grade: Grade) -> Result<Self> {
        let f = grade.filtration;
        let leading = grade.basis(model, f, 1)?;
        let mut page = Self {
            grade,
            words: Vec::new(),
            leading_start: 0,
            span: Basis::default(),
            lifts: Vec::new(),
            transition: OnceCell::new(),
        };
        if leading.is_empty() {
            return Ok(page);
        }
        if grade.page > 1 {
            page.words = grade.basis(model, f + 1, grade.page - 1)?;
        }
        page.leading_start = page.words.len();
        page.words.extend(leading);
        let word_cells = page.words.iter().map(|m| 2 + 2 * m.len()).sum();
        check_cells(word_cells)?;
        model.allocate(word_cells)?;
        if grade.page > 0 {
            page.build_quotient(model)?;
        }
        Ok(page)
    }

    fn leading(&self) -> &[Monomial] {
        &self.words[self.leading_start..]
    }

    fn build_quotient(&mut self, model: &mut Model) -> Result<()> {
        let grade = self.grade;
        let prime = grade.prime();
        let f = grade.filtration;
        let cutoff = grade.cutoff();
        let mut span = self.boundaries(model)?;
        let mut rows = Rows::default();
        let mut images = Basis::default();
        let mut lift_cells = 0;
        for (index, monomial) in self.words.iter().enumerate() {
            let input = Polynomial::from_terms(model, [(monomial, 1)])?;
            let output = grade.truncated_d(model, &input, f, cutoff)?;
            let mut column = Vec::new();
            for (term, coefficient) in output.terms() {
                rows.add(model, term, coefficient, &mut column)?;
            }
            sort_column(model, &mut column)?;
            let Some(lift) = images.insert_with_witness(model, prime, column, vec![(index, 1)])?
            else {
                check_cells(rows.cells + images.cells() + lift_cells + span.cells())?;
                continue;
            };
            if index < self.leading_start {
                continue;
            }
            // Higher-filtration corrections were processed first. Their
            // dependencies have zero component in the leading quotient.
            let leading = lift
                .iter()
                .filter(|&&(i, _)| i >= self.leading_start)
                .map(|&(i, coefficient)| (i - self.leading_start, coefficient))
                .collect();
            if span
                .insert_with_witness(model, prime, leading, vec![(self.lifts.len(), 1)])?
                .is_none()
            {
                lift_cells += 2 + 2 * lift.len();
                check_cells(rows.cells + images.cells() + lift_cells + span.cells())?;
                self.lifts.push(lift.into_boxed_slice());
            }
        }
        model.allocate(lift_cells + span.cells())?;
        self.span = span;
        Ok(())
    }

    fn boundaries(&self, model: &mut Model) -> Result<Basis> {
        let mut span = Basis::default();
        let grade = self.grade;
        if grade.shape.arrows == 0 {
            return Ok(span);
        }
        let grade = Grade {
            shape: bp_ext::previous_shape(grade.shape)?,
            ..grade
        };
        let prime = grade.prime();
        let f = grade.filtration;
        let first = f.saturating_sub(grade.page - 1);
        let mut rows = Rows::default();
        let mut images = Basis::default();
        for monomial in grade.basis(model, first, f - first + 1)? {
            let input = Polynomial::from_terms(model, [(&monomial, 1)])?;
            let output = grade.truncated_d(model, &input, first, f + 1)?;
            let mut lower = Vec::new();
            let mut leading = Vec::new();
            for (term, coefficient) in output.terms() {
                if grade.exponent(term) < f {
                    rows.add(model, term, coefficient, &mut lower)?;
                } else {
                    leading.push((
                        bp_ext::index_of(model, self.leading(), term)?,
                        coefficient as u32,
                    ));
                }
            }
            sort_column(model, &mut lower)?;
            // A dependence with image zero below f is a boundary in F^f,
            // and its witness records that boundary's component at f.
            if let Some(image) = images.insert_with_witness(model, prime, lower, leading)? {
                span.insert(model, prime, image)?;
            }
            check_cells(rows.cells + images.cells() + span.cells())?;
        }
        Ok(span)
    }

    pub(crate) fn rank(&self) -> usize {
        if self.grade.page == 0 {
            self.leading().len()
        } else {
            self.lifts.len()
        }
    }

    fn project(&self, model: &mut Model, cochain: &Polynomial) -> Result<Vec<Coordinate>> {
        if cochain.is_zero() {
            return Ok(Vec::new());
        }
        let grade = self.grade;
        bp_ext::validate_homogeneous(model, grade.shape, grade.t, cochain)?;
        let mut leading = Vec::new();
        for (monomial, coefficient) in cochain.terms() {
            match grade.exponent(monomial).cmp(&grade.filtration) {
                Ordering::Less => {
                    return Err(Error::plain(
                        "BP cochain lies below the requested v_h-adic filtration",
                    ))
                }
                Ordering::Equal => leading.push((
                    bp_ext::index_of(model, self.leading(), monomial)?,
                    coefficient as u32,
                )),
                Ordering::Greater => {}
            }
        }
        if grade.page == 0 {
            return Ok(leading);
        }
        let cutoff = grade.cutoff();
        if !grade
            .truncated_d(model, cochain, grade.filtration, cutoff)?
            .is_zero()
        {
            return Err(Error::plain(
                "BP cochain does not satisfy the v_h-adic page cycle condition",
            ));
        }
        let (remainder, coordinates) =
            self.span
                .reduce_with_witness(model, grade.prime(), leading)?;
        if !remainder.is_empty() {
            return Err(Error::plain(
                "BP cochain is missing from the v_h-adic page basis",
            ));
        }
        Ok(coordinates)
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        coordinates: &[Coordinate],
    ) -> Result<Polynomial> {
        if self.grade.page == 0 {
            return self.decode(model, coordinates);
        }
        let prime = self.grade.prime();
        let mut terms = Vec::new();
        for &(index, coefficient) in coordinates {
            terms = add_scaled(model, prime, &terms, &self.lifts[index], coefficient)?;
        }
        self.decode(model, &terms)
    }

    fn decode(&self, model: &mut Model, coordinates: &[Coordinate]) -> Result<Polynomial> {
        Polynomial::from_terms(
            model,
            coordinates
                .iter()
                .map(|&(i, c)| (&self.words[i], i128::from(c))),
        )
    }

    fn build_transition(&self, model: &mut Model, next: TypeId) -> Result<Basis> {
        let grade = self.grade;
        let prime = grade.prime();
        let mut span = Basis::default();
        if let (Some(a), Some(f)) = (
            grade.shape.arrows.checked_sub(1),
            grade.filtration.checked_sub(grade.page),
        ) {
            let incoming_ty = model.bp_vpage_type(Grade {
                shape: grade.shape.with_arrows(a)?,
                filtration: f,
                ..grade
            })?;
            let incoming = model.bp_vpage(incoming_ty)?;
            for index in 0..incoming.rank() {
                let lift = incoming.representative(model, &[(index, 1)])?;
                let image = incoming
                    .grade
                    .truncated_d(model, &lift, f, grade.cutoff())?;
                let image = self.project(model, &image)?;
                span.insert(model, prime, image)?;
            }
        }
        let next = model.bp_vpage(next)?;
        for index in 0..next.rank() {
            let lift = next.representative(model, &[(index, 1)])?;
            let image = self.project(model, &lift)?;
            if span
                .insert_with_witness(model, prime, image, vec![(index, 1)])?
                .is_some()
            {
                return Err(Error::plain("BPVPage transition has dependent cycle lifts"));
            }
        }
        model.allocate(span.cells())?;
        Ok(span)
    }
}

pub(crate) fn class(model: &mut Model, ty: TypeId, cochain: &Polynomial) -> Result<Value> {
    if cochain.is_zero() {
        return Ok(odd_ext::zero());
    }
    let page = model.bp_vpage(ty)?;
    let coordinates = page.project(model, cochain)?;
    odd_ext::persist(model, coordinates)
}

pub(crate) fn apply(
    model: &mut Model,
    operation: Operation,
    source: TypeId,
    target: TypeId,
    coordinates: &[Coordinate],
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(if matches!(operation, Operation::Rep) {
            bp::zero()
        } else {
            odd_ext::zero()
        });
    }
    let page = model.bp_vpage(source)?;
    if matches!(operation, Operation::Next) {
        let span = if let Some(span) = page.transition.get() {
            span
        } else {
            let span = page.build_transition(model, target)?;
            page.transition.get_or_init(|| span)
        };
        let (remainder, image) =
            span.reduce_with_witness(model, page.grade.prime(), coordinates.to_vec())?;
        if !remainder.is_empty() {
            return Err(Error::plain(format!(
                "bpVNext requires a cycle for d_{}",
                page.grade.page
            )));
        }
        return odd_ext::persist(model, image);
    }
    let lift = page.representative(model, coordinates)?;
    match operation {
        Operation::Rep => bp::persist(model, lift),
        Operation::Differential => {
            let Kind::BPVPage(grade) = *model.kind(target) else {
                unreachable!()
            };
            let image =
                page.grade
                    .truncated_d(model, &lift, page.grade.filtration, grade.cutoff())?;
            class(model, target, &image)
        }
        Operation::Previous => class(model, target, &lift),
        Operation::Next => unreachable!(),
    }
}

pub(crate) fn cup(
    model: &mut Model,
    left: (TypeId, &[Coordinate]),
    right: (TypeId, &[Coordinate]),
    target: TypeId,
) -> Result<Value> {
    if left.1.is_empty() || right.1.is_empty() {
        return Ok(odd_ext::zero());
    }
    let (&Kind::BPVPage(a), &Kind::BPVPage(b)) = (model.kind(left.0), model.kind(right.0)) else {
        unreachable!()
    };
    let f = model.bp_vpage(left.0)?.representative(model, left.1)?;
    let g = model.bp_vpage(right.0)?.representative(model, right.1)?;
    let product = bp::cup(model, a.shape, b.shape, &f, &g)?;
    class(model, target, &product)
}
