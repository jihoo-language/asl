use std::cell::OnceCell;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::mixed::{Complex, Differential, Grade, Id};
use crate::mixed_bars::Bars;
use crate::model::Model;
use crate::prime_linear::{add_scaled, Basis, Coordinate};

pub(crate) struct Slice {
    grade: Grade,
    bars: Rc<Bars>,
    living: Vec<usize>,
    decomposition: OnceCell<Basis>,
}

impl Slice {
    pub(crate) fn new(model: &mut Model, grade: Grade) -> Result<Self> {
        let bars = model.mixed_bars(grade.complex)?;
        model.spend(bars.bars().len())?;
        let living: Vec<_> = bars
            .bars()
            .iter()
            .enumerate()
            .filter_map(|(i, bar)| bar.lives_at(grade.b).then_some(i))
            .collect();
        model.allocate(living.len())?;
        Ok(Self {
            grade,
            bars,
            living,
            decomposition: OnceCell::new(),
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.living.len()
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        class: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let prime = model.mixed_complex(self.grade.complex).prime;
        let mut chain = Vec::new();
        model.spend(class.len())?;
        for &(i, coefficient) in class {
            let bar = &self.bars.bars()[self.living[i]];
            chain = add_scaled(model, prime, &chain, &bar.cycle, coefficient)?;
        }
        Ok(chain)
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        chain: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        if chain.is_empty() {
            return Ok(Vec::new());
        }
        let complex = model.mixed_complex(self.grade.complex);
        if !complex.apply(model, Differential::Total, chain)?.is_empty() {
            return Err(Error::plain("mixHClass requires a total cycle"));
        }
        let (remaining, coordinates) = self.projector(model, &complex)?.reduce_with_witness(
            model,
            complex.prime,
            chain.to_vec(),
        )?;
        if !remaining.is_empty() {
            return Err(Error::plain(
                "mixed bars did not span the cycle modulo boundaries",
            ));
        }
        Ok(coordinates)
    }

    fn projector(&self, model: &mut Model, complex: &Complex) -> Result<&Basis> {
        if self.decomposition.get().is_none() {
            let mut basis = incoming(model, self.grade, complex, false)?;
            // Boundary rows have empty labels. A pivot from a living bar
            // carries its coordinate in the chosen homology basis.
            for (index, &i) in self.living.iter().enumerate() {
                let image = self.bars.bars()[i].cycle.to_vec();
                if basis
                    .insert_with_witness(model, complex.prime, image, vec![(index, 1)])?
                    .is_some()
                {
                    return Err(Error::plain(
                        "a living mixed bar was dependent modulo boundaries",
                    ));
                }
            }
            model.allocate(basis.cells())?;
            let _ = self.decomposition.set(basis);
        }

        Ok(self
            .decomposition
            .get()
            .expect("mixed homology decomposition initialized"))
    }

    pub(crate) fn shift(
        &self,
        model: &mut Model,
        target: &Self,
        class: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        model.spend(class.len() * (1 + target.living.len().max(1).ilog2() as usize))?;
        Ok(class
            .iter()
            .filter_map(|&(i, c)| {
                target
                    .living
                    .binary_search(&self.living[i])
                    .ok()
                    .map(|index| (index, c))
            })
            .collect())
    }
}

fn incoming(model: &mut Model, grade: Grade, complex: &Complex, track: bool) -> Result<Basis> {
    let mut basis = Basis::default();
    for i in complex.basis(model, i128::from(grade.b) + 1, true)? {
        let image = complex.total_column(model, i)?;
        if !image.is_empty() {
            let witness = if track { vec![(i, 1)] } else { Vec::new() };
            basis.insert_with_witness(model, complex.prime, image, witness)?;
        }
    }
    Ok(basis)
}

pub(crate) fn primitive(
    model: &mut Model,
    grade: Grade,
    chain: &[Coordinate],
) -> Result<Vec<Coordinate>> {
    if chain.is_empty() {
        return Ok(Vec::new());
    }
    let complex = model.mixed_complex(grade.complex);
    let (remaining, preimage) = incoming(model, grade, &complex, true)?.reduce_with_witness(
        model,
        complex.prime,
        chain.to_vec(),
    )?;
    if !remaining.is_empty() {
        return Err(Error::plain("mixHPrimitive requires an exact total chain"));
    }
    Ok(preimage)
}

pub(crate) fn map(
    model: &mut Model,
    source: Grade,
    target: Grade,
    class: &[Coordinate],
    evaluate: impl FnOnce(&mut Model, &[Coordinate]) -> Result<Vec<Coordinate>>,
) -> Result<Vec<Coordinate>> {
    if class.is_empty() {
        return Ok(Vec::new());
    }
    let from = model.mixed_homology(source)?;
    let chain = from.representative(model, class)?;
    let image = evaluate(model, &chain)?;
    if image.is_empty() {
        return Ok(image);
    }
    model.mixed_homology(target)?.project(model, &image)
}

pub(crate) fn induce(
    model: &mut Model,
    map_id: Id,
    support: &[Coordinate],
    source: Grade,
    target: Grade,
    class: &[Coordinate],
) -> Result<Vec<Coordinate>> {
    let hom = model.mixed_hom(map_id)?;
    let complex = model.mixed_complex(map_id);
    if !complex
        .apply(model, Differential::Total, support)?
        .is_empty()
    {
        return Err(Error::plain("mixInduce requires a closed total Hom map"));
    }
    map(model, source, target, class, |model, chain| {
        hom.evaluate(model, support, chain)
    })
}
