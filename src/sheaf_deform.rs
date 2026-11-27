use std::cell::OnceCell;
use std::cmp::Ordering;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::prime_linear::Basis;
use crate::sheaf::{self, Cochain, Cochains, Layer, Support, Term};
use crate::sheaf_ext::{self, Hom};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Term(usize),
    Curvature,
    Obstruction,
    Lift,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Term(_) => "shDefTerm",
            Self::Curvature => "shDefCurvature",
            Self::Obstruction => "shDefObstruction",
            Self::Lift => "shDefLift",
        }
    }
}

pub(crate) struct Jet {
    coefficients: Box<[Support]>,
    curvature: OnceCell<Support>,
}

impl Jet {
    pub(crate) fn new(model: &mut Model, coefficients: Vec<Support>) -> Result<Rc<Self>> {
        sheaf::check_cells(1 + coefficients.len())?;
        model.allocate(1 + coefficients.len())?;
        Ok(Rc::new(Self {
            coefficients: coefficients.into(),
            curvature: OnceCell::new(),
        }))
    }

    pub(crate) fn coefficients(&self) -> &[Support] {
        &self.coefficients
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        for (a, b) in self.coefficients.iter().zip(rhs.coefficients.iter()) {
            model.spend(1)?;
            if !Rc::ptr_eq(a, b) {
                let order = sheaf::compare(model, a, b)?;
                if !order.is_eq() {
                    return Ok(order);
                }
            }
        }
        Ok(self.coefficients.len().cmp(&rhs.coefficients.len()))
    }

    pub(crate) fn curvature(&self, model: &mut Model, hom: &Hom) -> Result<Support> {
        if let Some(value) = self.curvature.get() {
            return Ok(value.clone());
        }
        let mut result = Cochain::default();
        let prime = hom.prime();
        model.spend(self.coefficients.len())?;
        for (a, b) in self.coefficients.iter().zip(self.coefficients.iter().rev()) {
            // The sign of (1,1) Hom composition is negative.
            for (generator, coefficient) in sheaf_ext::compose(model, (hom, 1, a), (hom, 1, b))? {
                result.add(model, prime, generator, prime.negate(coefficient))?;
            }
        }
        let result = retained(model, result.into_terms())?;
        let _ = self.curvature.set(result.clone());
        Ok(result)
    }

    pub(crate) fn append(
        &self,
        model: &mut Model,
        hom: &Hom,
        coefficient: Support,
    ) -> Result<Rc<Self>> {
        let curvature = self.curvature(model, hom)?;
        let differential = hom.differential_support(model, &coefficient)?;
        if !cancel(model, hom, &curvature, &differential)? {
            return Err(Error::plain(
                "shDefAppend coefficient does not satisfy dA + curvature = 0",
            ));
        }
        self.extend(model, coefficient)
    }

    pub(crate) fn lift(&self, model: &mut Model, hom: &Hom) -> Result<Rc<Self>> {
        let curvature = self.curvature(model, hom)?;
        let coefficient = primitive(model, hom, &curvature)?;
        let coefficient = retained(model, coefficient)?;
        self.extend(model, coefficient)
    }

    fn extend(&self, model: &mut Model, coefficient: Support) -> Result<Rc<Self>> {
        sheaf::check_cells(2 + self.coefficients.len())?;
        model.spend(1 + self.coefficients.len())?;
        let mut coefficients = self.coefficients.to_vec();
        coefficients.push(coefficient);
        Self::new(model, coefficients)
    }

    pub(crate) fn gauge(
        &self,
        model: &mut Model,
        hom: &Hom,
        exponent: usize,
        gauge: &[Term],
        witness: bool,
    ) -> Result<Value> {
        let count = self.coefficients.len() + usize::from(witness);
        sheaf::check_cells(count)?;
        model.spend(count)?;
        let mut output = Vec::with_capacity(count);
        let mut temporary = count;
        let prime = hom.prime();
        for i in 1..=count {
            let old = self.coefficients.get(i - 1).map_or(&[][..], |entry| entry);
            if i < exponent || gauge.is_empty() {
                output.push(self.coefficients.get(i - 1).cloned().unwrap_or_default());
                continue;
            }
            let mut total = Cochain::default();
            for (generator, coefficient) in old {
                total.add(model, prime, generator.clone(), *coefficient)?;
            }
            if i == exponent {
                // Q R_0 - R_0 Q is the Roos differential of Q.
                for (generator, coefficient) in hom.differential_support(model, gauge)? {
                    total.add(model, prime, generator, coefficient)?;
                }
            } else {
                // R'_i = R_i + Q R_(i-e) - R'_(i-e) Q, from R'G = G R.
                let lower = i - exponent - 1;
                for (generator, coefficient) in
                    sheaf_ext::compose(model, (hom, 0, gauge), (hom, 1, &self.coefficients[lower]))?
                {
                    total.add(model, prime, generator, coefficient)?;
                }
                for (generator, coefficient) in
                    sheaf_ext::compose(model, (hom, 1, &output[lower]), (hom, 0, gauge))?
                {
                    total.add(model, prime, generator, prime.negate(coefficient))?;
                }
            }
            let terms = total.into_terms();
            temporary += terms
                .iter()
                .map(|(generator, _)| generator.cells())
                .sum::<usize>();
            sheaf::check_cells(temporary)?;
            output.push(terms.into());
        }
        if witness {
            let last = output.pop().expect("a next coefficient");
            model.allocate(last.iter().map(|(generator, _)| generator.cells()).sum())?;
            Ok(Value::SheafC(last))
        } else {
            model.allocate(temporary - count)?;
            Ok(Value::SheafDeformation(Self::new(model, output)?))
        }
    }
}

pub(crate) fn constant(model: &mut Model, hom: &Hom) -> Result<Value> {
    let data = &hom.source;
    let mut result = Cochain::default();
    for (target, &rows) in data.dims().iter().enumerate() {
        model.spend(1 + data.lower(target).len())?;
        for &source in data.lower(target) {
            let chain = Rc::<[usize]>::from([source, target]);
            for (column, entries) in data.restriction(source, target).columns() {
                model.spend(1 + entries.len())?;
                for &(row, coefficient) in entries {
                    let generator = sheaf_ext::coordinate(chain.clone(), rows, column, row)?;
                    result.add(model, hom.prime(), generator, coefficient)?;
                }
            }
        }
    }
    sheaf::persist(model, result.into_terms())
}

fn retained(model: &mut Model, terms: Vec<Term>) -> Result<Support> {
    let Value::SheafC(terms) = sheaf::persist(model, terms)? else {
        unreachable!()
    };
    Ok(terms)
}

fn cancel(model: &mut Model, hom: &Hom, a: &[Term], b: &[Term]) -> Result<bool> {
    if a.len() != b.len() {
        return Ok(false);
    }
    for ((x, c), (y, d)) in a.iter().zip(b) {
        model.spend(x.cells())?;
        if x != y || hom.prime().negate(*c) != *d {
            return Ok(false);
        }
    }
    Ok(true)
}

pub(crate) struct Boundaries {
    source: Layer,
    target: Layer,
    basis: Basis,
}

impl Boundaries {
    pub(crate) fn build(model: &mut Model, hom: &Hom, degree: usize) -> Result<Self> {
        let previous = degree.checked_sub(1).ok_or_else(|| {
            Error::plain("degree-zero sheaf cochains have no incoming differential")
        })?;
        let source = hom.layer(model, previous, MAX_PIVOT_CELLS)?;
        let target = hom.layer(model, degree, MAX_PIVOT_CELLS)?;
        let cells = source.cells() + target.cells();
        sheaf::check_cells(cells)?;
        let mut basis = Basis::default();
        for i in 0..source.size {
            let image = hom.indexed_column(model, &source.generator(i), &target)?;
            basis.insert_with_witness(model, hom.prime(), image, vec![(i, 1)])?;
            sheaf::check_cells(cells + basis.cells())?;
        }
        Ok(Self {
            source,
            target,
            basis,
        })
    }

    pub(crate) fn cells(&self) -> usize {
        self.source.cells() + self.target.cells() + self.basis.cells()
    }

    pub(crate) fn preimage(
        &self,
        model: &mut Model,
        hom: &Hom,
        raw: &[Term],
        factor: u32,
    ) -> Result<Option<Vec<Term>>> {
        if raw.is_empty() || factor == 0 {
            return Ok(Some(Vec::new()));
        }
        let prime = hom.prime();
        let target = raw
            .iter()
            .map(|(generator, coefficient)| {
                Ok((
                    self.target.index(model, generator)?,
                    prime.mul(*coefficient, factor),
                ))
            })
            .collect::<Result<_>>()?;
        let (remainder, witness) = self.basis.reduce_with_witness(model, prime, target)?;
        if !remainder.is_empty() {
            return Ok(None);
        }
        Ok(Some(
            witness
                .into_iter()
                .map(|(i, coefficient)| (self.source.generator(i), coefficient))
                .collect(),
        ))
    }
}

fn primitive(model: &mut Model, hom: &Hom, curvature: &[Term]) -> Result<Vec<Term>> {
    if curvature.is_empty() {
        return Ok(Vec::new());
    }
    Boundaries::build(model, hom, 2)?
        .preimage(model, hom, curvature, hom.prime().negate(1))?
        .ok_or_else(|| Error::plain("shDefLift has a nonzero obstruction in sheaf Ext degree two"))
}

pub(crate) fn elements(
    model: &mut Model,
    hom: &Hom,
    order: usize,
    limit: usize,
) -> Result<Vec<Value>> {
    model.spend(order)?;
    if order == 0 {
        return Ok(vec![Value::SheafDeformation(Jet::new(model, Vec::new())?)]);
    }
    let values = hom.elements(model, 1, limit)?;
    let mut cells = 3 * values.len();
    let mut candidates = Vec::with_capacity(values.len());
    for value in values {
        let Value::SheafC(coefficient) = value else {
            unreachable!()
        };
        let differential = hom.differential_support(model, &coefficient)?;
        cells += differential
            .iter()
            .map(|(generator, _)| generator.cells())
            .sum::<usize>();
        sheaf::check_cells(cells)?;
        candidates.push((differential, coefficient));
    }
    if candidates.len() == 1 {
        sheaf::check_cells(1 + order)?;
        let jet = Jet::new(model, vec![candidates[0].1.clone(); order])?;
        return Ok(vec![Value::SheafDeformation(jet)]);
    }
    let mut jets = vec![Jet::new(model, Vec::new())?];
    for _ in 0..order {
        let mut next = Vec::new();
        for jet in jets {
            let curvature = jet.curvature(model, hom)?;
            for (differential, coefficient) in &candidates {
                model.spend(1)?;
                if cancel(model, hom, &curvature, differential)? {
                    if next.len() == limit {
                        return Err(Error::plain(
                            "sheaf deformations exceed the finite type enumeration limit",
                        ));
                    }
                    next.push(jet.extend(model, coefficient.clone())?);
                }
            }
        }
        jets = next;
    }
    Ok(jets.into_iter().map(Value::SheafDeformation).collect())
}
