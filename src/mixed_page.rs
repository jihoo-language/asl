use std::cell::OnceCell;
use std::cmp::Ordering;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::mixed::{Complex, Differential, Grade as Chain};
use crate::model::{Kind, Model, TypeId, Value};
use crate::odd_ext;
use crate::prime_linear::{add_scaled, Basis, Coordinate};

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) chain: Chain,
    pub(crate) r: usize,
    pub(crate) f: usize,
}

impl Grade {
    pub(crate) fn differential(self) -> Result<Self> {
        let f = self
            .f
            .checked_add(self.r)
            .ok_or_else(|| Error::plain("mixed page filtration is out of range"))?;
        Ok(Self {
            chain: self.chain.shift(-1)?,
            f,
            ..self
        })
    }

    pub(crate) fn next(self) -> Result<Self> {
        let r = self
            .r
            .checked_add(1)
            .ok_or_else(|| Error::plain("mixed page number is out of range"))?;
        Ok(Self { r, ..self })
    }

    pub(crate) fn previous(self) -> Result<Self> {
        let r = self
            .r
            .checked_sub(1)
            .ok_or_else(|| Error::plain("mixPrev requires a positive page number"))?;
        Ok(Self { r, ..self })
    }

    fn cutoff(self) -> i128 {
        self.f as i128 + self.r as i128
    }

    fn incoming(self) -> Option<Self> {
        let f = self.f.checked_sub(self.r)?;
        let chain = self.chain.shift(1).ok()?;
        Some(Self { chain, f, ..self })
    }
}

pub(crate) struct Page {
    pub(crate) grade: Grade,
    indices: Vec<usize>,
    leading_start: usize,
    span: Basis,
    lifts: Vec<Box<[Coordinate]>>,
    next_basis: OnceCell<Basis>,
}

impl Page {
    pub(crate) fn build(model: &mut Model, grade: Grade) -> Result<Self> {
        let complex = model.mixed_complex(grade.chain.complex);
        let mut page = Self {
            grade,
            indices: Vec::new(),
            leading_start: 0,
            span: Basis::default(),
            lifts: Vec::new(),
            next_basis: OnceCell::new(),
        };
        let f = grade.f as i128;
        let b = i128::from(grade.chain.b);
        let mut leading = Vec::new();
        for i in complex.basis(model, b, true)? {
            let j = complex.filtration(i, b).unwrap();
            if j == f {
                leading.push(i);
            } else if f < j && j < grade.cutoff() {
                page.indices.push(i);
            }
        }
        if leading.is_empty() {
            page.indices.clear();
            return Ok(page);
        }
        page.leading_start = page.indices.len();
        page.indices.extend(leading);
        model.allocate(page.indices.len())?;
        if grade.r != 0 {
            page.build_quotient(model, &complex)?;
        }
        Ok(page)
    }

    fn leading(&self) -> &[usize] {
        &self.indices[self.leading_start..]
    }

    pub(crate) fn rank(&self) -> usize {
        if self.grade.r == 0 {
            self.leading().len()
        } else {
            self.lifts.len()
        }
    }

    fn build_quotient(&mut self, model: &mut Model, complex: &Complex) -> Result<()> {
        let mut span = self.boundaries(model, complex)?;
        let mut cycles = Basis::default();
        let prime = complex.prime;
        let b = i128::from(self.grade.chain.b);
        let f = self.grade.cutoff();
        let mut cells = 0;
        for (position, &index) in self.indices.iter().enumerate() {
            let mut column = complex.apply(model, Differential::Total, &[(index, 1)])?;
            model.spend(column.len())?;
            column.retain(|&(j, _)| complex.filtration(j, b - 1).unwrap() < f);
            let Some(lift) =
                cycles.insert_with_witness(model, prime, column, vec![(position, 1)])?
            else {
                continue;
            };
            if position < self.leading_start {
                continue;
            }
            let leading = lift
                .iter()
                .filter(|&&(i, _)| i >= self.leading_start)
                .map(|&(i, coefficient)| (i - self.leading_start, coefficient))
                .collect();
            if span
                .insert_with_witness(model, prime, leading, vec![(self.lifts.len(), 1)])?
                .is_none()
            {
                cells += 1 + 2 * lift.len();
                if cells + span.cells() > MAX_PIVOT_CELLS {
                    return Err(Error::plain(
                        "mixed page reduction exceeds the storage limit",
                    ));
                }
                self.lifts.push(lift.into());
            }
        }
        model.allocate(cells + span.cells())?;
        self.span = span;
        Ok(())
    }

    fn boundaries(&self, model: &mut Model, complex: &Complex) -> Result<Basis> {
        let b = i128::from(self.grade.chain.b);
        let f = self.grade.f as i128;
        let floor = self.grade.f.saturating_sub(self.grade.r - 1) as i128;
        let mut boundaries = Basis::default();
        let mut lower = Basis::default();
        for index in complex.basis(model, b + 1, true)? {
            let j = complex.filtration(index, b + 1).unwrap();
            if j < floor || j > f {
                continue;
            }
            let image = complex.apply(model, Differential::Total, &[(index, 1)])?;
            model.spend(image.len())?;
            let mut earlier = Vec::new();
            let mut leading = Vec::new();
            for (i, coefficient) in image {
                match complex.filtration(i, b).unwrap().cmp(&f) {
                    Ordering::Less => earlier.push((i, coefficient)),
                    Ordering::Equal => leading.push((self.index(model, i)?, coefficient)),
                    Ordering::Greater => {}
                }
            }
            if let Some(image) =
                lower.insert_with_witness(model, complex.prime, earlier, leading)?
            {
                boundaries.insert(model, complex.prime, image)?;
            }
        }
        Ok(boundaries)
    }

    fn index(&self, model: &mut Model, index: usize) -> Result<usize> {
        model.spend(1 + self.leading().len().max(1).ilog2() as usize)?;
        self.leading()
            .binary_search(&index)
            .map_err(|_| Error::plain("mixed page leading term is outside the basis"))
    }

    fn project(&self, model: &mut Model, input: &[Coordinate]) -> Result<Vec<Coordinate>> {
        let complex = model.mixed_complex(self.grade.chain.complex);
        let b = i128::from(self.grade.chain.b);
        let cutoff = self.grade.cutoff();
        let f = self.grade.f as i128;
        let mut leading = Vec::new();
        let mut terms = Vec::new();
        model.spend(input.len())?;
        for &(index, coefficient) in input {
            let j = complex
                .filtration(index, b)
                .ok_or_else(|| Error::plain("mixed chain term has the wrong total degree"))?;
            if j < f {
                return Err(Error::plain(
                    "mixed chain lies below the requested filtration",
                ));
            }
            if j == f {
                // A nonzero leading coefficient implies this layer exists.
                leading.push((self.index(model, index)?, coefficient));
            }
            if j < cutoff {
                terms.push((index, coefficient));
            }
        }
        if self.grade.r == 0 {
            return Ok(leading);
        }
        let image = complex.apply(model, Differential::Total, &terms)?;
        if image
            .iter()
            .any(|&(i, _)| complex.filtration(i, b - 1).unwrap() < cutoff)
        {
            return Err(Error::plain("mixed chain fails the page cycle condition"));
        }
        let (remainder, coordinates) =
            self.span
                .reduce_with_witness(model, complex.prime, leading)?;
        if !remainder.is_empty() {
            return Err(Error::plain(
                "mixed page cycle is missing from the computed basis",
            ));
        }
        Ok(coordinates)
    }

    fn representative(&self, model: &mut Model, class: &[Coordinate]) -> Result<Vec<Coordinate>> {
        let complex = model.mixed_complex(self.grade.chain.complex);
        if self.grade.r == 0 {
            model.spend(class.len())?;
            return Ok(class.iter().map(|&(i, c)| (self.leading()[i], c)).collect());
        }
        let mut support = Vec::new();
        for &(i, coefficient) in class {
            support = add_scaled(model, complex.prime, &support, &self.lifts[i], coefficient)?;
        }
        for (position, _) in &mut support {
            *position = self.indices[*position];
        }
        model.spend(support.len() * (1 + support.len().max(1).ilog2() as usize))?;
        support.sort_unstable();
        Ok(support)
    }

    fn transition(&self, model: &mut Model, target: TypeId) -> Result<Basis> {
        let complex = model.mixed_complex(self.grade.chain.complex);
        let mut cycles = Basis::default();
        if let Some(incoming) = self.grade.incoming() {
            let ty = model.mixed_page_type(incoming)?;
            let page = model.mixed_page(ty)?;
            for i in 0..page.rank() {
                let lift = page.representative(model, &[(i, 1)])?;
                let image = complex.apply(model, Differential::Total, &lift)?;
                let coords = self.project(model, &image)?;
                cycles.insert(model, complex.prime, coords)?;
            }
        }
        let next = model.mixed_page(target)?;
        for i in 0..next.rank() {
            let lift = next.representative(model, &[(i, 1)])?;
            let coords = self.project(model, &lift)?;
            if cycles
                .insert_with_witness(model, complex.prime, coords, vec![(i, 1)])?
                .is_some()
            {
                return Err(Error::plain("dependent mixed page transition lifts"));
            }
        }
        model.allocate(cycles.cells())?;
        Ok(cycles)
    }
}

pub(crate) fn class(model: &mut Model, ty: TypeId, input: &[Coordinate]) -> Result<Value> {
    if input.is_empty() {
        return Ok(odd_ext::zero());
    }
    let class = model.mixed_page(ty)?.project(model, input)?;
    odd_ext::persist(model, class)
}

pub(crate) fn representative(model: &mut Model, ty: TypeId, input: &[Coordinate]) -> Result<Value> {
    if input.is_empty() {
        return Ok(odd_ext::zero());
    }
    let terms = model.mixed_page(ty)?.representative(model, input)?;
    odd_ext::persist(model, terms)
}

pub(crate) fn differential(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    coordinates: &[Coordinate],
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(odd_ext::zero());
    }
    let Kind::MixedPage(grade) = *model.kind(input) else {
        unreachable!()
    };
    let complex = model.mixed_complex(grade.chain.complex);
    let terms = model
        .mixed_page(input)?
        .representative(model, coordinates)?;
    let image = complex.apply(model, Differential::Total, &terms)?;
    class(model, output, &image)
}

pub(crate) fn previous(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    coordinates: &[Coordinate],
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(odd_ext::zero());
    }
    let terms = model
        .mixed_page(input)?
        .representative(model, coordinates)?;
    class(model, output, &terms)
}

pub(crate) fn next(
    model: &mut Model,
    input: TypeId,
    output: TypeId,
    coordinates: &[Coordinate],
) -> Result<Value> {
    if coordinates.is_empty() {
        return Ok(odd_ext::zero());
    }
    let page = model.mixed_page(input)?;
    let complex = model.mixed_complex(page.grade.chain.complex);
    let transition = match page.next_basis.get() {
        Some(basis) => basis,
        None => {
            let basis = page.transition(model, output)?;
            page.next_basis.get_or_init(|| basis)
        }
    };
    let (remainder, class) =
        transition.reduce_with_witness(model, complex.prime, coordinates.to_vec())?;
    if !remainder.is_empty() {
        return Err(Error::plain(format!(
            "mixNext requires a cycle for d_{}",
            page.grade.r
        )));
    }
    odd_ext::persist(model, class)
}
