use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::lee_linear::Echelon;
use crate::model::Model;
use crate::rv_cyclotomic;
use crate::rv_lens::{self, Rational};

use super::super::defects::relative::tube::{Element, Tube};
use super::{Double, Field, Move, State};

impl Double {
    fn torus_coordinates(&self, model: &mut Model, state: &State) -> Result<Vec<Rational>> {
        let catalog = self.surface_catalog(model)?;
        let mut field = Field::new(model, self.conductor);
        let mut result = vec![Rational::zero(); self.count()];
        let inv = rv_lens::fraction(field.model, &[1], catalog.order as u128)?;
        for (m, x) in &state.0 {
            let [g, b] = *m.as_slice() else {
                return Err(Error::plain(
                    "surface Fourier transform requires a genus-one marking",
                ));
            };
            if catalog.charge(m)? != 0 {
                return Err(Error::plain(
                    "surface Fourier transform requires the neutral torus",
                ));
            }
            let lo = self.flux_offsets[g];
            let hi = self.flux_offsets[g + 1];
            for (slot, s) in result[lo..hi].iter_mut().zip(&self.simples[lo..hi]) {
                if let Some(p) = s.phases[b] {
                    let v = field.shift(x, p)?;
                    let q = field.multiply(&v, &Rational::integer(s.dimension as i128))?;
                    *slot = field.add(slot, &q)?;
                }
            }
        }
        for x in &mut result {
            *x = field.multiply(x, &inv)?;
        }
        Ok(result)
    }

    fn torus_markings(&self, model: &mut Model, coords: &[Rational]) -> Result<State> {
        let catalog = self.surface_catalog(model)?;
        if coords.len() != self.count() {
            return Err(Error::plain(
                "torus Fourier coordinate count disagrees with the bulk",
            ));
        }
        let mut field = Field::new(model, self.conductor);
        let mut result = BTreeMap::new();
        for (s, x) in self
            .simples
            .iter()
            .zip(coords)
            .filter(|(_, x)| !x.is_zero())
        {
            let dx = field.multiply(x, &Rational::integer(s.dimension as i128))?;
            for (b, p) in s.phases.iter().copied().enumerate() {
                if let Some(phase) = p {
                    let m = vec![s.flux, b];
                    if catalog.charge(&m)? != 0 {
                        return Err(Error::plain(
                            "bulk projective character is supported on a charged torus",
                        ));
                    }
                    let c = field.shift(&dx, (self.conductor - phase) % self.conductor)?;
                    let cell = result.entry(m).or_insert_with(Rational::zero);
                    *cell = field.add(cell, &c)?;
                }
            }
        }
        Ok(State(
            result.into_iter().filter(|(_, x)| !x.is_zero()).collect(),
        ))
    }

    pub(crate) fn surface_to_tube(
        &self,
        model: &mut Model,
        tube: &Tube,
        state: &State,
    ) -> Result<Element> {
        let native = self.torus_coordinates(model, state)?;
        let embedded = native
            .iter()
            .map(|x| embed(model, self.conductor, tube.conductor, x))
            .collect::<Result<Vec<_>>>()?;
        tube.lift_center(model, &embedded)
    }

    pub(crate) fn surface_from_tube(
        &self,
        model: &mut Model,
        tube: &Tube,
        state: &Element,
    ) -> Result<State> {
        let full = tube.extract_center(model, state)?;
        let native = descend(model, self.conductor, tube.conductor, full)?;
        self.torus_markings(model, &native)
    }

    pub(crate) fn surface_verify_modular(&self, model: &mut Model) -> Result<String> {
        let catalog = self.surface_catalog(model)?;
        let dim = catalog.dimensions(model, self, 1)?[0];
        if dim != self.count() as u128 {
            return Err(Error::plain(
                "neutral torus size disagrees with the rank of its Fourier basis",
            ));
        }
        for a in 0..self.count() {
            let mut basis = vec![Rational::zero(); self.count()];
            basis[a] = Rational::integer(1);
            let marking = self.torus_markings(model, &basis)?;
            if self.torus_coordinates(model, &marking)? != basis {
                return Err(Error::plain(
                    "marked torus Fourier inverse does not recover a bulk color",
                ));
            }
            let t = marking.moved(model, self, Move::T(0, false), 0)?;
            let t = self.torus_coordinates(model, &t)?;
            let s = marking.moved(model, self, Move::S(0, false), 0)?;
            let s = self.torus_coordinates(model, &s)?;
            let mut field = Field::new(model, self.conductor);
            let mut expected_t = vec![Rational::zero(); self.count()];
            expected_t[a] = field.phase(self.twists[a])?;
            if t != expected_t {
                return Err(Error::plain(
                    "marked torus T chain disagrees with the bulk positive twist",
                ));
            }
            for (r, x) in s.iter().enumerate() {
                let y = self.s[r][a];
                let root = field.phase(y.phase)?;
                let weight = rv_lens::fraction(field.model, &[y.weight], catalog.order as u128)?;
                if *x != field.multiply(&root, &weight)? {
                    return Err(Error::plain(
                        "marked torus S chain disagrees with the bulk orientation",
                    ));
                }
            }
        }
        let s = Move::S(0, false);
        let t = Move::T(0, false);
        for a in 0..catalog.order {
            for b in 0..catalog.order {
                let m = [a, b];
                if catalog.charge(&m)? != 0 {
                    continue;
                }
                let ss = word(self, &m, &[s, s])?;
                let cube = word(self, &m, &[t, s, t, s, t, s])?;
                let fourth = word(self, &m, &[s, s, s, s])?;
                if ss != cube || fourth.0 != m || fourth.1 != 0 {
                    return Err(Error::plain(
                        "marked torus chain lifts fail the modular identities on a neutral state",
                    ));
                }
            }
        }
        Ok(format!("rvDoubleSurfaceModular = neutral_states={dim}, conductor={}, fourier=true, S4=true, ST3_equals_S2=true",self.conductor))
    }
}

fn word(double: &Double, start: &[usize], steps: &[Move]) -> Result<(Vec<usize>, usize)> {
    let mut value = start.to_vec();
    let mut phase = 0;
    for &step in steps {
        let (m, _, p) = step.apply(double, &value)?;
        phase = (phase + p) % double.conductor;
        value = m;
    }
    Ok((value, phase))
}

fn embed(model: &mut Model, source: usize, target: usize, value: &Rational) -> Result<Rational> {
    if !target.is_multiple_of(source) {
        return Err(Error::plain(
            "tube Fourier field does not contain the marked surface conductor",
        ));
    }
    if source == target || value.numerator.len() <= 1 {
        return Ok(value.clone());
    }
    let step = target / source;
    let mut raw = vec![0; target];
    for (i, &c) in value.numerator.iter().enumerate() {
        raw[step * i] = c;
    }
    let p = rv_cyclotomic::reduced(model, target, raw)?;
    rv_lens::fraction(model, &p, value.denominator)
}

fn descend(
    model: &mut Model,
    source: usize,
    target: usize,
    values: Vec<Rational>,
) -> Result<Vec<Rational>> {
    if source == target {
        return Ok(values);
    }
    if !target.is_multiple_of(source) {
        return Err(Error::plain(
            "tube Fourier field does not contain the marked surface conductor",
        ));
    }
    let d = model.cyclotomic(source)?.len() - 1;
    let mut span = Echelon::default();
    for j in 0..d {
        let mut raw = vec![0; target / source * j + 1];
        *raw.last_mut().unwrap() = 1;
        let vector = rv_cyclotomic::reduced(model, target, raw)?
            .into_iter()
            .enumerate()
            .filter(|&(_, v)| v != 0)
            .collect();
        if span.insert(model, vector, vec![(j, 1)])?.is_some() {
            return Err(Error::plain("cyclotomic surface embedding is singular"));
        }
    }
    values
        .into_iter()
        .map(|v| {
            if v.numerator.len() <= 1 {
                return Ok(v);
            }
            let image = v
                .numerator
                .iter()
                .copied()
                .enumerate()
                .filter(|&(_, x)| x != 0)
                .collect();
            let relation = span.relation(model, image, vec![(d, 1)])?.ok_or_else(|| {
                Error::plain("tube center amplitude does not descend to the marked surface field")
            })?;
            let &(_, n) = relation.last().filter(|(j, _)| *j == d).ok_or_else(|| {
                Error::plain("cyclotomic surface descent lost its marked coefficient")
            })?;
            let mut raw = vec![0; d];
            for &(j, c) in &relation[..relation.len() - 1] {
                raw[j] = if n < 0 {
                    c
                } else {
                    c.checked_neg().ok_or_else(|| {
                        Error::plain("cyclotomic surface descent exceeds its coefficient bound")
                    })?
                };
            }
            let x = rv_lens::fraction(model, &raw, n.unsigned_abs())?;
            let d = rv_lens::fraction(model, &[1], v.denominator)?;
            rv_lens::multiply(model, source, &x, &d)
        })
        .collect()
}
