use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hw_closed::Curve;
use crate::model::{Model, Value};

use super::super::resolution::{quartic, QuarticTerm, ResolvedPoint};
use super::super::{theta, Point};
use super::{Generator, Word};

const NODE_COUNT: usize = 16;
const DIRECT_STEPS: u128 = 16;

pub(super) struct Geometry {
    curve: Rc<Curve>,
    quartic: Vec<QuarticTerm>,
    translations: Option<theta::Group>,
    centers: [Option<Rc<Point>>; NODE_COUNT],
}

impl Geometry {
    fn compile(model: &mut Model, word: &Word) -> Result<Self> {
        let curve = model.kummer_curve(&word.basis.curve_name)?;
        let quartic = quartic(model, &curve)?;
        let needs_translation = word.generators.iter().any(|&generator| match generator {
            Generator::Translate(mask) | Generator::Switch(mask) => mask != 0,
            Generator::Project(_) => false,
        });
        let translations = needs_translation
            .then(|| theta::Group::build(model, &word.basis, curve.clone()))
            .transpose()?;
        let mut centers = std::array::from_fn(|_| None);
        model.spend(word.generators.len())?;
        for &generator in &word.generators {
            if let Generator::Project(mask) = generator {
                let center = &mut centers[usize::from(mask)];
                if center.is_none() {
                    let torsion = word.basis.binary_combination(model, mask)?;
                    let point = Point::project(model, &curve, &torsion)?;
                    if !point.is_node(model, &curve)? {
                        return Err(Error::plain(
                            "marked K3 projection center must be a Kummer node",
                        ));
                    }
                    *center = Some(point.keep(model)?);
                }
            }
        }
        model.allocate(NODE_COUNT + 4 + quartic.len() * (7 + curve.ring.degree()))?;
        Ok(Self {
            curve,
            quartic,
            translations,
            centers,
        })
    }

    fn translate(
        &self,
        model: &mut Model,
        mask: u8,
        point: ResolvedPoint,
    ) -> Result<ResolvedPoint> {
        if mask == 0 {
            return Ok(point);
        }
        let group = self
            .translations
            .as_ref()
            .expect("compiled translations for the word");
        point.projective_image_transient(
            model,
            &self.curve,
            Some(&self.quartic),
            |model, coordinates| group.apply(model, usize::from(mask), coordinates),
        )
    }

    fn apply(
        &self,
        model: &mut Model,
        generator: Generator,
        point: ResolvedPoint,
    ) -> Result<ResolvedPoint> {
        match generator {
            Generator::Translate(mask) => self.translate(model, mask, point),
            Generator::Switch(mask) => {
                let switched = point.gauss_switch_transient(model, &self.curve, &self.quartic)?;
                self.translate(model, mask, switched)
            }
            Generator::Project(mask) => {
                let center = self.centers[usize::from(mask)]
                    .as_ref()
                    .expect("compiled projection center for the word");
                point.node_involution_transient(model, &self.curve, center, &self.quartic)
            }
        }
    }
}

impl Word {
    pub(crate) fn curve_name(&self) -> Rc<str> {
        self.basis.curve_name.clone()
    }

    fn geometry(&self, model: &mut Model) -> Result<&Geometry> {
        if let Some(geometry) = self.point_geometry.get() {
            return Ok(geometry);
        }
        let geometry = Geometry::compile(model, self)?;
        Ok(self.point_geometry.get_or_init(|| geometry))
    }

    fn point_action(
        &self,
        model: &mut Model,
        geometry: &Geometry,
        mut point: ResolvedPoint,
        reverse: bool,
    ) -> Result<ResolvedPoint> {
        // Each listed generator is an involution, so reverse the word to invert it.
        let mut generators = self.generators.iter();
        while let Some(&generator) = if reverse {
            generators.next_back()
        } else {
            generators.next()
        } {
            model.spend(1)?;
            point = geometry.apply(model, generator, point)?;
        }
        Ok(point)
    }

    pub(crate) fn iterate(
        &self,
        model: &mut Model,
        point: Rc<ResolvedPoint>,
        power: i128,
    ) -> Result<Value> {
        if power == 0 || self.generators.is_empty() {
            return Ok(Value::KummerK3(point));
        }
        let geometry = self.geometry(model)?;
        let reverse = power < 0;
        let mut steps = power.unsigned_abs();
        let mut image = point.as_ref().clone();
        if steps > DIRECT_STEPS {
            let (returned, first_return) =
                self.walk_to_return(model, geometry, &point, reverse, steps)?;
            let Some(period) = first_return else {
                return returned.persist_image(model);
            };
            steps %= period;
            if steps == 0 {
                return Ok(Value::KummerK3(point));
            }
        }
        for _ in 0..steps {
            image = self.point_action(model, geometry, image, reverse)?;
        }
        image.persist_image(model)
    }

    pub(crate) fn period(&self, model: &mut Model, point: &ResolvedPoint) -> Result<u128> {
        if self.generators.is_empty() {
            return Ok(1);
        }
        let geometry = self.geometry(model)?;
        let (_, period) = self.walk_to_return(model, geometry, point, false, u128::MAX)?;
        period
            .ok_or_else(|| Error::plain("K3 automorphism orbit period exceeds the supported range"))
    }

    fn walk_to_return(
        &self,
        model: &mut Model,
        geometry: &Geometry,
        point: &ResolvedPoint,
        reverse: bool,
        steps: u128,
    ) -> Result<(ResolvedPoint, Option<u128>)> {
        // A permutation orbit has no tail. Equality with the start detects
        // its exact period without keeping the intermediate points.
        let mut current = point.clone();
        for step in 1..=steps {
            current = self.point_action(model, geometry, current, reverse)?;
            if current.compare(model, point)?.is_eq() {
                return Ok((current, Some(step)));
            }
        }
        Ok((current, None))
    }

    pub(super) fn cycles(&self, model: &mut Model, name: &str) -> Result<String> {
        let ty = model.kummer_k3_type(self.curve_name())?;
        let points = model.elements(ty)?;
        let total = points.len();
        model.spend(total)?;
        if self.generators.is_empty() {
            return Ok(format!(
                "rvK3AutoCycles[{name}] = total={total}, fixed={total}, cycles=[(1, {total})]"
            ));
        }
        let geometry = self.geometry(model)?;
        let mut seen = vec![false; total];
        let mut counts = BTreeMap::<usize, usize>::new();
        for start in 0..total {
            if seen[start] {
                continue;
            }
            let mut current = start;
            let mut length = 0;
            while !seen[current] {
                seen[current] = true;
                let Value::KummerK3(point) = &points[current] else {
                    unreachable!()
                };
                let image = self.point_action(model, geometry, point.as_ref().clone(), false)?;
                current = model.index(ty, &Value::KummerK3(Rc::new(image)))?;
                length += 1;
            }
            if current != start {
                return Err(Error::plain(
                    "K3 automorphism word did not permute the resolved points",
                ));
            }
            *counts.entry(length).or_default() += 1;
        }
        let fixed = counts.get(&1).copied().unwrap_or(0);
        let cycles: Vec<_> = counts.into_iter().collect();
        Ok(format!(
            "rvK3AutoCycles[{name}] = total={total}, fixed={fixed}, cycles={cycles:?}"
        ))
    }
}
