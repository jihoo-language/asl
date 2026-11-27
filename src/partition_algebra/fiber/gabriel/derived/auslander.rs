use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::Model;
use crate::syntax::ParsedExpr;

use super::super::super::{
    add, add_multiple, multiply, ratio, Coordinates, Element, Fiber, Mode, Parameter,
};
use super::{Nakayama, Path, State};

mod mesh;
pub(crate) mod tau;
use mesh::Mesh;

type Terms = BTreeMap<usize, Rational>;
const ONE: Rational = Rational::integer(1);

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Grade {
    pub(crate) parameter: Parameter,
    pub(crate) source: [u128; 2],
    pub(crate) target: [u128; 2],
}

impl Grade {
    pub(crate) fn render(self) -> String {
        let Parameter { rank, n } = self.parameter;
        let [i, l] = self.source;
        let [j, m] = self.target;
        format!("RVPartitionHom({rank},{n},{i},{l},{j},{m})")
    }

    pub(crate) fn compose(g: Self, f: Self) -> Result<Self> {
        if g.parameter != f.parameter || g.source != f.target {
            return Err(Error::plain("partition Hom composition requires the same fiber and matching intermediate module"));
        }
        Ok(Self {
            parameter: f.parameter,
            source: f.source,
            target: g.target,
        })
    }
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Map {
    pub(crate) grade: Grade,
    terms: Terms,
}

impl Map {
    pub(crate) fn zero(grade: Grade) -> Self {
        Self {
            grade,
            terms: Terms::new(),
        }
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<crate::model::Value> {
        model.allocate(4 + 6 * self.terms.len())?;
        Ok(crate::model::Value::RVPartitionHom(std::rc::Rc::new(self)))
    }

    pub(crate) fn render(&self) -> String {
        format!("{}{{{}}}", self.grade.render(), terms_display(&self.terms))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<std::cmp::Ordering> {
        model.spend(1 + self.terms.len().min(rhs.terms.len()))?;
        Ok(self.cmp(rhs))
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let mut terms = self.terms.clone();
        model.spend(self.terms.len() + rhs.terms.len())?;
        for (&layer, &coefficient) in &rhs.terms {
            accumulate(&mut terms, layer, coefficient)?;
        }
        Ok(Self {
            grade: self.grade,
            terms,
        })
    }

    pub(crate) fn scale(&self, model: &mut Model, n: i128, d: i128) -> Result<Self> {
        let coefficient = ratio(n, d)?;
        let mut value = Self::zero(self.grade);
        model.spend(self.terms.len())?;
        if coefficient != Rational::ZERO {
            for (&layer, &term) in &self.terms {
                value.terms.insert(layer, multiply(coefficient, term)?);
            }
        }
        Ok(value)
    }
}

pub(crate) enum Expression {
    Basis(Grade, u128),
    Compose(Box<[ParsedExpr; 2]>),
    Scale(i128, i128, Box<ParsedExpr>),
    Apply(u128, Box<ParsedExpr>),
    Arrow(Parameter, u128),
    Radical(Parameter, [u128; 2]),
    ModulePath(Parameter, [u128; 3]),
    MeshMap(Parameter, [u128; 4]),
}

pub(crate) enum Report {
    Indecomposables(Parameter),
    Algebra(Parameter),
    Mesh(Parameter, [u128; 2]),
    HomDim(Grade),
    Apply(u128, Box<ParsedExpr>),
    RadicalContains(u128, Box<ParsedExpr>),
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct Basis {
    source: usize,
    target: usize,
    shift: usize,
}

impl Basis {
    fn tuple(self) -> (usize, usize, usize) {
        (self.source, self.target, self.shift)
    }
}

pub(crate) struct Auslander {
    parameter: Parameter,
    modules: Vec<State>,
    dimensions: Vec<usize>,
    homs: Vec<Vec<Vec<usize>>>,
    powers: Vec<BTreeSet<Basis>>,
    arrows: Vec<Basis>,
    meshes: BTreeMap<usize, Mesh>,
    quiver_layers: Vec<usize>,
    tau_fan: RefCell<Option<Rc<tau::Fan>>>,
    projectives: BTreeSet<usize>,
}

impl Auslander {
    pub(super) fn new(model: &mut Model, algebra: &Nakayama) -> Result<Self> {
        let modules: Vec<_> = algebra.resolutions.keys().copied().collect();
        let dimensions = modules
            .iter()
            .map(|state| {
                (0..state.length)
                    .map(|length| {
                        algebra.simple_dimensions[algebra.end(Path {
                            start: state.head,
                            length,
                        })]
                    })
                    .sum()
            })
            .collect();
        let mut homs = vec![vec![Vec::new(); modules.len()]; modules.len()];
        let mut projectives = BTreeSet::new();
        for (source, &a) in modules.iter().enumerate() {
            if a.length == algebra.paths[a.head].len() {
                projectives.insert(source);
            }
            for (target, &b) in modules.iter().enumerate() {
                model.spend(b.length)?;
                homs[target][source] = (0..b.length)
                    .filter(|&shift| {
                        algebra.end(Path {
                            start: b.head,
                            length: shift,
                        }) == a.head
                            && shift + a.length >= b.length
                    })
                    .collect();
            }
        }
        let mut data = Self {
            parameter: algebra.parameter,
            modules,
            dimensions,
            homs,
            powers: Vec::new(),
            arrows: Vec::new(),
            meshes: BTreeMap::new(),
            quiver_layers: Vec::new(),
            projectives,
            tau_fan: RefCell::default(),
        };
        let mut basis = BTreeSet::new();
        for i in 0..data.modules.len() {
            for j in 0..data.modules.len() {
                basis.extend(data.basis(i, j));
            }
        }
        data.certify_homs(model, algebra, &basis)?;
        let rad: BTreeSet<_> = basis
            .iter()
            .copied()
            .filter(|&b| b.source != b.target || b.shift != 0)
            .collect();
        data.powers.push(basis);
        data.powers.push(rad.clone());
        while !data.powers.last().unwrap().is_empty() {
            let last = data.powers.last().unwrap();
            let mut next = BTreeSet::new();
            model.spend(last.len() * rad.len())?;
            for &a in last {
                for &b in &rad {
                    if let Some(c) = data.then(a, b) {
                        if !last.contains(&c) {
                            return Err(Error::plain("partition Auslander radical powers did not form a decreasing chain"));
                        }
                        next.insert(c);
                    }
                }
            }
            if next.len() == last.len() {
                return Err(Error::plain(
                    "partition Auslander radical did not become nilpotent",
                ));
            }
            data.powers.push(next);
        }
        let squared = data.powers.get(2);
        data.arrows = rad
            .into_iter()
            .filter(|a| squared.is_none_or(|row| !row.contains(a)))
            .collect();
        data.meshes = mesh::all_meshes(model, algebra, &data)?;
        data.quiver_layers = mesh::certify_presentation(model, &data)?;
        let work =
            data.powers.iter().map(BTreeSet::len).sum::<usize>() + 12 * data.modules.len().pow(2);
        model.allocate(4 * work)?;
        Ok(data)
    }

    pub(crate) fn tau_fan(&self, model: &mut Model, algebra: &Nakayama) -> Result<Rc<tau::Fan>> {
        if let Some(fan) = self.tau_fan.borrow().as_ref() {
            return Ok(fan.clone());
        }
        let fan = Rc::new(tau::Fan::new(model, algebra, self)?);
        *self.tau_fan.borrow_mut() = Some(fan.clone());
        Ok(fan)
    }

    fn certify_homs(
        &self,
        model: &mut Model,
        algebra: &Nakayama,
        basis: &BTreeSet<Basis>,
    ) -> Result<()> {
        for &a in basis {
            let source = self.modules[a.source];
            for l in source.length..algebra.paths[source.head].len() {
                if self.image_layer(algebra, a, l).is_some() {
                    return Err(Error::plain(
                        "partition Hom basis map does not factor through its source truncation",
                    ));
                }
            }
            for &b in basis.iter().filter(|b| a.target == b.source) {
                model.spend(1 + source.length)?;
                let product = self.then(a, b);
                if product.is_some_and(|c| !basis.contains(&c)) {
                    return Err(Error::plain(
                        "partition Hom composition escaped the canonical basis",
                    ));
                }
                if let Some(c) = product {
                    if c.source == c.target
                        && c.shift == 0
                        && (a.source != a.target || b.source != b.target)
                    {
                        return Err(Error::plain(
                            "partition Hom off-diagonal product has a split identity component",
                        ));
                    }
                }
                for l in 0..source.length {
                    let successive = self
                        .image_layer(algebra, a, l)
                        .and_then(|t| self.image_layer(algebra, b, t));
                    let together = product.and_then(|c| self.image_layer(algebra, c, l));
                    if successive != together {
                        return Err(Error::plain(
                            "partition Hom composition failed on an actual module path",
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    fn image_layer(&self, algebra: &Nakayama, map: Basis, layer: usize) -> Option<usize> {
        let i = self.modules[map.source];
        let j = self.modules[map.target];
        let product = algebra.then(
            Path {
                start: j.head,
                length: map.shift,
            },
            Path {
                start: i.head,
                length: layer,
            },
        )?;
        (product.length < j.length).then_some(product.length)
    }

    fn then(&self, a: Basis, b: Basis) -> Option<Basis> {
        let shift = a.shift + b.shift;
        (a.target == b.source && shift < self.modules[b.target].length).then_some(Basis {
            source: a.source,
            target: b.target,
            shift,
        })
    }

    fn basis(&self, source: usize, target: usize) -> impl Iterator<Item = Basis> + '_ {
        self.homs[target][source].iter().map(move |&shift| Basis {
            source,
            target,
            shift,
        })
    }

    fn module(&self, key: [u128; 2]) -> Result<usize> {
        let [head, length] = key.map(|x| usize::try_from(x).unwrap_or(usize::MAX));
        self.modules
            .binary_search(&State { head, length })
            .map_err(|_| {
                Error::plain("partition indecomposable module top or length is out of range")
            })
    }

    fn key(&self, module: usize) -> [u128; 2] {
        let state = self.modules[module];
        [state.head as u128, state.length as u128]
    }

    fn endpoints(&self, grade: Grade) -> Result<[usize; 2]> {
        Ok([self.module(grade.source)?, self.module(grade.target)?])
    }

    fn basis_map(&self, basis: Basis) -> Map {
        let grade = Grade {
            parameter: self.parameter,
            source: self.key(basis.source),
            target: self.key(basis.target),
        };
        Map {
            grade,
            terms: BTreeMap::from([(basis.shift, ONE)]),
        }
    }

    pub(crate) fn dimension(&self, grade: Grade) -> Result<usize> {
        let [source, target] = self.endpoints(grade)?;
        Ok(self.homs[target][source].len())
    }

    pub(crate) fn hom(&self, grade: Grade, layer: u128) -> Result<Map> {
        let [source, target] = self.endpoints(grade)?;
        let shift = usize::try_from(layer).unwrap_or(usize::MAX);
        if self.homs[target][source].binary_search(&shift).is_err() {
            return Err(Error::plain(
                "partition Hom head-to-layer map is not in this Hom basis",
            ));
        }
        Ok(self.basis_map(Basis {
            source,
            target,
            shift,
        }))
    }

    pub(crate) fn compose(&self, model: &mut Model, g: &Map, f: &Map) -> Result<Map> {
        let grade = Grade::compose(g.grade, f.grade)?;
        let [source, middle] = self.endpoints(f.grade)?;
        let target = self.module(g.grade.target)?;
        model.spend(f.terms.len() * g.terms.len())?;
        let mut result = Map::zero(grade);
        for (&k, &a) in &f.terms {
            for (&l, &b) in &g.terms {
                if let Some(c) = self.then(
                    Basis {
                        source,
                        target: middle,
                        shift: k,
                    },
                    Basis {
                        source: middle,
                        target,
                        shift: l,
                    },
                ) {
                    accumulate(&mut result.terms, c.shift, multiply(a, b)?)?;
                }
            }
        }
        Ok(result)
    }

    fn apply(&self, model: &mut Model, algebra: &Nakayama, f: &Map, l: u128) -> Result<Terms> {
        let [source, target] = self.endpoints(f.grade)?;
        let l = usize::try_from(l)
            .ok()
            .filter(|&l| l < self.modules[source].length)
            .ok_or_else(|| Error::plain("partition module basis layer is out of range"))?;
        let mut terms = Terms::new();
        model.spend(1 + f.terms.len())?;
        for (&shift, &coefficient) in &f.terms {
            if let Some(layer) = self.image_layer(
                algebra,
                Basis {
                    source,
                    target,
                    shift,
                },
                l,
            ) {
                terms.insert(layer, coefficient);
            }
        }
        Ok(terms)
    }

    pub(crate) fn apply_report(
        &self,
        model: &mut Model,
        algebra: &Nakayama,
        f: &Map,
        l: u128,
    ) -> Result<String> {
        let coords = self.apply(model, algebra, f, l)?;
        Ok(format!(
            "rvPartitionHomApply[{l}] = [{}]",
            terms_display(&coords)
        ))
    }

    pub(crate) fn apply_path(
        &self,
        model: &mut Model,
        algebra: &Nakayama,
        fiber: &Fiber,
        f: &Map,
        l: u128,
    ) -> Result<Element> {
        let coordinates = self.apply(model, algebra, f, l)?;
        let module = self.modules[self.module(f.grade.target)?];
        let mut result = Coordinates::new();
        for (length, coefficient) in coordinates {
            add_multiple(
                &mut result,
                &algebra
                    .lifted(Path {
                        start: module.head,
                        length,
                    })
                    .value,
                coefficient,
            )?;
        }
        Ok(fiber.element(Mode::Full, result))
    }

    pub(crate) fn module_path(
        &self,
        algebra: &Nakayama,
        fiber: &Fiber,
        key: [u128; 3],
    ) -> Result<Element> {
        let [i, l, t] = key;
        let module = self.modules[self.module([i, l])?];
        let length = usize::try_from(t)
            .ok()
            .filter(|&x| x < module.length)
            .ok_or_else(|| Error::plain("partition module basis layer is out of range"))?;
        let path = algebra.lifted(Path {
            start: module.head,
            length,
        });
        Ok(fiber.element(Mode::Full, path.value.clone()))
    }

    pub(crate) fn arrow(&self, index: u128) -> Result<Map> {
        let index = usize::try_from(index).ok();
        let basis = index
            .and_then(|i| self.arrows.get(i))
            .ok_or_else(|| Error::plain("partition Auslander arrow index is out of range"))?;
        Ok(self.basis_map(*basis))
    }

    pub(crate) fn radical_basis(&self, key: [u128; 2]) -> Result<Map> {
        let [power, index] = key.map(|x| usize::try_from(x).unwrap_or(usize::MAX));
        let basis = self
            .powers
            .get(power)
            .and_then(|row| row.iter().nth(index))
            .ok_or_else(|| {
                Error::plain("partition Auslander radical basis index is out of range")
            })?;
        Ok(self.basis_map(*basis))
    }

    pub(crate) fn contains(&self, f: &Map, power: u128) -> Result<String> {
        let [source, target] = self.endpoints(f.grade)?;
        let row = usize::try_from(power).ok().and_then(|i| self.powers.get(i));
        let yes = f.terms.keys().all(|&shift| {
            row.is_some_and(|row| {
                row.contains(&Basis {
                    source,
                    target,
                    shift,
                })
            })
        });
        Ok(format!(
            "rvPartitionAuslanderRadicalContains[{power}] = {yes}"
        ))
    }

    pub(crate) fn indecomposables(&self) -> String {
        let modules: Vec<_> = self
            .modules
            .iter()
            .enumerate()
            .map(|(index, module)| {
                (
                    module.head,
                    module.length,
                    self.dimensions[index],
                    self.projectives.contains(&index),
                )
            })
            .collect();
        let Parameter { rank, n } = self.parameter;
        format!("rvPartitionIndecomposables[{rank},{n}] = modules={modules:?}")
    }

    pub(crate) fn hom_dimension(&self, grade: Grade) -> Result<String> {
        let [source, target] = self.endpoints(grade)?;
        let Parameter { rank, n } = self.parameter;
        let [i, l] = grade.source;
        let [j, m] = grade.target;
        Ok(format!(
            "rvPartitionHomDim[{rank},{n},{i},{l},{j},{m}] = dimension={}, layers={:?}",
            self.homs[target][source].len(),
            self.homs[target][source]
        ))
    }

    pub(crate) fn algebra_report(&self) -> String {
        let dims: Vec<_> = self.powers.iter().map(BTreeSet::len).collect();
        let cartan: Vec<Vec<_>> = self
            .homs
            .iter()
            .map(|row| row.iter().map(Vec::len).collect())
            .collect();
        let arrows: Vec<_> = self.arrows.iter().map(|a| a.tuple()).collect();
        let Parameter { rank, n } = self.parameter;
        format!("rvPartitionAuslander[{rank},{n}] = indecomposables={}, dimension={}, powers={dims:?}, cartan={cartan:?}, arrows={arrows:?}, meshes={}, pathQuotient={:?}",self.modules.len(),dims[0],self.meshes.len(),self.quiver_layers)
    }

    pub(crate) fn mesh_map(&self, model: &mut Model, key: [u128; 4]) -> Result<Map> {
        let [i, l, side, term] = key;
        let target = self.module([i, l])?;
        let mesh = self.meshes.get(&target).ok_or_else(|| {
            Error::plain("partition projective module has no almost split sequence ending in it")
        })?;
        let (basis, sign) = mesh.map(side, term)?;
        self.basis_map(basis).scale(model, sign, 1)
    }

    pub(crate) fn mesh_report(&self, key: [u128; 2]) -> Result<String> {
        let module = self.module(key)?;
        let Parameter { rank, n } = self.parameter;
        let [i, l] = key;
        let label = format!("rvPartitionAR[{rank},{n},{i},{l}]");
        match self.meshes.get(&module) {
            Some(mesh) => Ok(format!("{label} = {}", mesh.render(self))),
            None => Ok(format!("{label} = projective")),
        }
    }
}

fn terms_display(terms: &Terms) -> String {
    terms
        .iter()
        .map(|(j, c)| format!("({j},{c})"))
        .collect::<Vec<_>>()
        .join(", ")
}

fn accumulate(target: &mut Terms, key: usize, value: Rational) -> Result<()> {
    let total = add(target.get(&key).copied().unwrap_or(Rational::ZERO), value)?;
    if total == Rational::ZERO {
        target.remove(&key);
    } else {
        target.insert(key, total);
    }
    Ok(())
}
