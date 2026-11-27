use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::mixed::{self, Complex, Differential};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::odd_ext;
use crate::prime_linear::{add_scaled, Basis, Coordinate, QuotientBasis};

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Id(pub(crate) usize);

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) minimal: Id,
    pub(crate) b: i64,
}

impl Grade {
    pub(crate) fn shift(self, offset: i128) -> Result<Self> {
        let b = i64::try_from(i128::from(self.b) + offset)
            .map_err(|_| Error::plain("minimal mixed model degree is out of range"))?;
        Ok(Self { b, ..self })
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Projection,
    Homotopy,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Projection => "mixMinP",
            Self::Homotopy => "mixMinH",
        }
    }
}

pub(crate) struct Minimal {
    pub(crate) source: mixed::Id,
    complex: Rc<Complex>,
    boundary_lifts: Basis,
    quotient: QuotientBasis,
    degrees: Vec<i64>,
}

impl Minimal {
    pub(crate) fn build(model: &mut Model, source: mixed::Id) -> Result<Self> {
        let complex = model.mixed_complex(source);
        if !complex.is_strict() {
            return Err(Error::plain("mixminimal requires a strict mixed complex"));
        }

        let prime = complex.prime;
        let mut boundary_lifts = Basis::default();
        let mut boundaries = Basis::default();
        for index in 0..complex.len() {
            let column = complex.column(Differential::Down, index);
            if boundary_lifts
                .insert_with_witness(model, prime, column.to_vec(), vec![(index, 1)])?
                .is_none()
            {
                boundaries.insert(model, prime, column.to_vec())?;
                if boundaries.cells() + boundary_lifts.cells() > MAX_PIVOT_CELLS {
                    return Err(Error::plain("mixed contraction exceeds the storage limit"));
                }
            }
        }
        let quotient = QuotientBasis::build(model, prime, boundaries, complex.len(), |_, i| {
            Ok(complex.column(Differential::Down, i).to_vec())
        })?;
        let mut degrees = Vec::with_capacity(quotient.rank());
        for i in 0..quotient.rank() {
            let representative = quotient.representative(model, prime, &[(i, 1)])?;
            degrees.push(complex.raw_degree(representative[0].0));
        }
        model.allocate(4 + boundary_lifts.cells() + degrees.len())?;
        Ok(Self {
            source,
            complex,
            boundary_lifts,
            quotient,
            degrees,
        })
    }

    pub(crate) fn prime(&self) -> Prime {
        self.complex.prime
    }

    pub(crate) fn degree(&self, generator: usize) -> i64 {
        self.degrees[generator]
    }

    pub(crate) fn filtration(&self, generator: usize, degree: i64) -> Option<i128> {
        let difference = i128::from(self.degree(generator)) - i128::from(degree);
        (difference >= 0 && difference % 2 == 0).then_some(difference / 2)
    }

    pub(crate) fn basis(&self, model: &mut Model, b: i64, total: bool) -> Result<Vec<usize>> {
        model.spend(self.degrees.len())?;
        Ok((0..self.degrees.len())
            .filter(|&index| {
                if total {
                    self.filtration(index, b).is_some()
                } else {
                    self.degree(index) == b
                }
            })
            .collect())
    }

    pub(crate) fn elements(
        &self,
        model: &mut Model,
        degree: i64,
        total: bool,
        limit: usize,
    ) -> Result<Vec<Value>> {
        let basis = self.basis(model, degree, total)?;
        let mut values = odd_ext::elements(model, self.prime(), basis.len(), limit)?;
        for value in &mut values {
            let Value::OddClass(support) = value else {
                unreachable!()
            };
            for (index, _) in Rc::make_mut(support) {
                *index = basis[*index];
            }
        }
        Ok(values)
    }

    fn primitive(&self, model: &mut Model, boundary: Vec<Coordinate>) -> Result<Vec<Coordinate>> {
        let (remainder, lift) =
            self.boundary_lifts
                .reduce_with_witness(model, self.prime(), boundary)?;
        if !remainder.is_empty() {
            return Err(Error::plain("mixed contraction expected a boundary"));
        }
        Ok(lift)
    }

    fn closed(&self, model: &mut Model, chain: &[Coordinate]) -> Result<Vec<Coordinate>> {
        let image = self.complex.apply(model, Differential::Down, chain)?;
        let lifted = self.primitive(model, image)?;
        add_scaled(model, self.prime(), chain, &lifted, self.prime().negate(1))
    }

    fn project_closed(
        &self,
        model: &mut Model,
        closed: Vec<Coordinate>,
    ) -> Result<Vec<Coordinate>> {
        self.quotient
            .project(model, self.prime(), closed)?
            .ok_or_else(|| Error::plain("mixed contraction expected a cycle"))
    }

    fn project(&self, model: &mut Model, chain: &[Coordinate]) -> Result<Vec<Coordinate>> {
        let closed = self.closed(model, chain)?;
        self.project_closed(model, closed)
    }

    fn include(&self, model: &mut Model, class: &[Coordinate]) -> Result<Vec<Coordinate>> {
        self.quotient.representative(model, self.prime(), class)
    }

    fn split(
        &self,
        model: &mut Model,
        chain: &[Coordinate],
    ) -> Result<(Vec<Coordinate>, Vec<Coordinate>)> {
        let closed = self.closed(model, chain)?;
        let class = self.project_closed(model, closed.clone())?;
        let inclusion = self.include(model, &class)?;
        let remainder = add_scaled(
            model,
            self.prime(),
            &closed,
            &inclusion,
            self.prime().negate(1),
        )?;
        Ok((class, self.primitive(model, remainder)?))
    }

    fn homotopy(&self, model: &mut Model, chain: &[Coordinate]) -> Result<Vec<Coordinate>> {
        Ok(self.split(model, chain)?.1)
    }

    fn negate(&self, model: &mut Model, chain: &mut [Coordinate]) -> Result<()> {
        model.spend(chain.len())?;
        for (_, coefficient) in chain {
            *coefficient = self.prime().negate(*coefficient);
        }
        Ok(())
    }

    pub(crate) fn delta(
        &self,
        model: &mut Model,
        order: usize,
        class: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        if order == 0 || class.is_empty() {
            return Ok(Vec::new());
        }
        let mut active = self.include(model, class)?;
        for _ in 1..order {
            let image = self.complex.apply(model, Differential::Up, &active)?;
            active = self.homotopy(model, &image)?;
            self.negate(model, &mut active)?;
            if active.is_empty() {
                return Ok(active);
            }
        }
        let image = self.complex.apply(model, Differential::Up, &active)?;
        self.project(model, &image)
    }

    pub(crate) fn inclusion(
        &self,
        model: &mut Model,
        class: &[Coordinate],
        total: bool,
    ) -> Result<Vec<Coordinate>> {
        if total {
            self.perturb_inclusion(model, class, false)
        } else {
            self.include(model, class)
        }
    }

    pub(crate) fn differential(
        &self,
        model: &mut Model,
        class: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        self.perturb_inclusion(model, class, true)
    }

    fn perturb_inclusion(
        &self,
        model: &mut Model,
        class: &[Coordinate],
        differential: bool,
    ) -> Result<Vec<Coordinate>> {
        let mut active = self.include(model, class)?;
        let mut result = Vec::new();
        while !active.is_empty() {
            let image = self.complex.apply(model, Differential::Up, &active)?;
            let (row, next) = if differential {
                self.split(model, &image)?
            } else {
                (active, self.homotopy(model, &image)?)
            };
            result = add_scaled(model, self.prime(), &result, &row, 1)?;
            active = next;
            self.negate(model, &mut active)?;
        }
        Ok(result)
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        chain: &[Coordinate],
        operation: Operation,
        total: bool,
    ) -> Result<Vec<Coordinate>> {
        if !total {
            return match operation {
                Operation::Projection => self.project(model, chain),
                Operation::Homotopy => self.homotopy(model, chain),
            };
        }
        let mut active = chain.to_vec();
        let mut result = Vec::new();
        while !active.is_empty() {
            let (row, lift) = match operation {
                Operation::Projection => self.split(model, &active)?,
                Operation::Homotopy => {
                    let lift = self.homotopy(model, &active)?;
                    (lift.clone(), lift)
                }
            };
            result = add_scaled(model, self.prime(), &result, &row, 1)?;
            active = self.complex.apply(model, Differential::Up, &lift)?;
            self.negate(model, &mut active)?;
        }
        Ok(result)
    }
}
