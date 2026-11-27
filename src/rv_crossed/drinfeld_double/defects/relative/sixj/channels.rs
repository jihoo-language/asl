use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::rv_lens::Rational;

use super::super::linear::{self, GradedMap};
use super::super::tensor::Layout;
use super::super::tree::{self, Branch, Node, Tree};
use super::super::{decompose, Ambient, Balanced, Hom, Image, Module};
use super::matrix::{self, Field, Matrix};
use super::{Category, Interval};

const MAX_RECOUPLING_CELLS: usize = 1_048_576;

pub(super) struct Work<'a> {
    pub(super) ambient: Ambient<'a>,
    pub(super) field: Field<'a>,
    categories: BTreeMap<Interval, &'a Category>,
    used_cells: usize,
    simples: BTreeMap<(Interval, usize), Rc<Module>>,
    spaces: BTreeMap<([usize; 3], [usize; 2]), Rc<Space>>,
    vertices: BTreeMap<VertexKey, Rc<[GradedMap]>>,
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
struct VertexKey {
    endpoints: [usize; 3],
    simples: [usize; 3],
}

struct Space {
    factors: [Rc<Module>; 2],
    product: Balanced,
    homs: Vec<Hom>,
}

pub(super) struct Channel {
    simple: usize,
    labels: Vec<(usize, usize)>,
    map: GradedMap,
}

impl Channel {
    pub(super) fn triple_path(&self) -> [usize; 3] {
        let [(internal, first), (_, second)] = *self.labels.as_slice() else {
            unreachable!("triple fusion has two trivalent vertices");
        };
        [internal, first, second]
    }
}

pub(super) struct Transition {
    pub(super) matrices: [Matrix; 2],
    pub(super) norms: [Vec<Rational>; 2],
    pub(super) grades: Vec<Vec<Image>>,
}

impl<'a> Work<'a> {
    pub(super) fn new(
        field: Field<'a>,
        ambient: Ambient<'a>,
        categories: BTreeMap<Interval, &'a Category>,
    ) -> Self {
        Self {
            ambient,
            field,
            categories,
            used_cells: 0,
            simples: BTreeMap::new(),
            spaces: BTreeMap::new(),
            vertices: BTreeMap::new(),
        }
    }

    fn reserve(&mut self, source: &Module, target: &Module) -> Result<()> {
        let cells = source
            .fibers
            .iter()
            .zip(&target.fibers)
            .map(|(s, t)| 1 + s * t)
            .sum::<usize>();
        self.used_cells = self
            .used_cells
            .checked_add(cells)
            .filter(|&n| n <= MAX_RECOUPLING_CELLS)
            .ok_or_else(|| {
                Error::plain("boundary recoupling exceeds its graded matrix storage limit")
            })?;
        Ok(())
    }

    fn simple(&mut self, interval: Interval, index: usize) -> Result<Rc<Module>> {
        if let Some(module) = self.simples.get(&(interval, index)) {
            return Ok(module.clone());
        }
        let module = Rc::new(self.ambient.induce(
            self.field.model,
            self.categories[&interval],
            index,
        )?);
        self.simples.insert((interval, index), module.clone());
        Ok(module)
    }

    fn space(
        &mut self,
        bounds @ [l, s, r]: [usize; 3],
        pair @ [a, b]: [usize; 2],
    ) -> Result<Rc<Space>> {
        if let Some(space) = self.spaces.get(&(bounds, pair)) {
            return Ok(space.clone());
        }
        let factors = [self.simple([l, s], a)?, self.simple([s, r], b)?];
        let product = Balanced::new(self.field.model, &self.ambient, &factors[0], &factors[1])?;
        let homs = decompose(
            self.field.model,
            &self.ambient,
            self.categories[&[l, r]],
            &product.module,
        )?;
        let space = Rc::new(Space {
            factors,
            product,
            homs,
        });
        self.spaces.insert((bounds, pair), space.clone());
        Ok(space)
    }

    fn vertices(
        &mut self,
        bounds @ [l, _, r]: [usize; 3],
        colors @ [a, b, c]: [usize; 3],
    ) -> Result<Rc<[GradedMap]>> {
        let key = VertexKey {
            endpoints: bounds,
            simples: colors,
        };
        if let Some(result) = self.vertices.get(&key) {
            return Ok(result.clone());
        }
        let space = self.space(bounds, [a, b])?;
        let h = &space.homs[c];
        if h.orbits.is_empty() {
            return Ok(Rc::from([]));
        }
        let simple = self.simple([l, r], c)?;
        let g = self.categories[&[l, r]].simple(c)?.simple.flux;
        let mut maps = Vec::new();
        for orbit in &h.orbits {
            let mut base = Matrix::zero(h.rows, h.cols)?;
            self.field.model.spend(orbit.len())?;
            for &(p, phase) in orbit {
                base.set(p / h.cols, p % h.cols, self.field.phase(phase)?);
            }
            self.reserve(&simple, &space.product.module)?;
            maps.push(GradedMap::extend(
                &mut self.field,
                &self.ambient,
                &simple,
                &space.product.module,
                g,
                &base,
            )?);
        }
        let maps: Rc<[GradedMap]> = maps.into();
        self.vertices.insert(key, maps.clone());
        Ok(maps)
    }

    pub(super) fn basis(
        &mut self,
        tree: &Tree<'_>,
        inputs: &[usize],
        wanted: Option<usize>,
    ) -> Result<Vec<Channel>> {
        let Node::Join(branch) = &tree.node else {
            let Node::Leaf(i, module) = &tree.node else {
                unreachable!()
            };
            self.reserve(module, module)?;
            let grades = module
                .fibers
                .iter()
                .map(|&dim| Matrix::identity(dim))
                .collect::<Result<Vec<_>>>()?;
            return Ok(vec![Channel {
                simple: inputs[*i],
                labels: Vec::new(),
                map: GradedMap { grades },
            }]);
        };
        let expected = decompose(
            self.field.model,
            &self.ambient,
            self.categories[&tree.interval],
            tree.module(),
        )?;
        let mut sectors = BTreeMap::new();
        for (c, h) in expected.iter().enumerate() {
            if !h.orbits.is_empty() && wanted.is_none_or(|d| c == d) {
                sectors.insert(c, h.orbits.len());
            }
        }
        if sectors.is_empty() {
            return Ok(Vec::new());
        }
        let [left, right] = &branch.children;
        let a = self.basis(left, inputs, None)?;
        let b = self.basis(right, inputs, None)?;
        let [l, r] = tree.interval;
        let bounds = [l, left.interval[1], r];
        let mut result = Vec::new();
        let mut count = BTreeMap::new();
        for x in &a {
            for y in &b {
                let space = self.space(bounds, [x.simple, y.simple])?;
                for &c in sectors.keys() {
                    let vertices = self.vertices(bounds, [x.simple, y.simple, c])?;
                    if vertices.is_empty() {
                        continue;
                    }
                    let g = self.categories[&[l, r]].simple(c)?.simple.flux;
                    let source = self.simple([l, r], c)?;
                    let tensor = self.tensor(branch, &space, [x, y], g)?;
                    for (b, elementary) in vertices.iter().enumerate() {
                        let base = tensor.multiply(&mut self.field, &elementary.grades[g])?;
                        self.reserve(&source, tree.module())?;
                        let map = GradedMap::extend(
                            &mut self.field,
                            &self.ambient,
                            &source,
                            tree.module(),
                            g,
                            &base,
                        )?;
                        let mut labels = x.labels.clone();
                        labels.extend_from_slice(&y.labels);
                        labels.push((c, b));
                        result.push(Channel {
                            simple: c,
                            labels,
                            map,
                        });
                        *count.entry(c).or_insert(0) += 1;
                    }
                }
            }
        }
        if count != sectors {
            return Err(Error::plain(
                "boundary fusion paths disagree with the independent balanced Hom dimensions",
            ));
        }
        result.sort_by(|a, b| (a.simple, &a.labels).cmp(&(b.simple, &b.labels)));
        Ok(result)
    }

    fn tensor(
        &mut self,
        branch: &Branch<'_>,
        space: &Space,
        inclusions: [&Channel; 2],
        grade: usize,
    ) -> Result<Matrix> {
        let from = Layout {
            product: &space.product,
            factors: space.factors.each_ref().map(Rc::as_ref),
        };
        let to = Layout {
            product: &branch.product,
            factors: branch.children.each_ref().map(Tree::module),
        };
        linear::tensor_at(
            &mut self.field,
            &self.ambient,
            from,
            to,
            inclusions.map(|x| &x.map),
            grade,
        )
    }

    pub(super) fn transition(
        &mut self,
        inputs: &[Module],
        ends: [(&Tree<'_>, &[Channel]); 2],
        rotation: [usize; 4],
        grade: usize,
    ) -> Result<Transition> {
        let [(from, source), (to, target)] = ends;
        let images = tree::rotation(self.field.model, &self.ambient, inputs, from, to, rotation)?;
        let source = source
            .iter()
            .map(|b| &b.map.grades[grade])
            .collect::<Vec<_>>();
        let target = target
            .iter()
            .map(|b| &b.map.grades[grade])
            .collect::<Vec<_>>();
        let moved = source
            .iter()
            .map(|b| {
                b.transformed(&mut self.field, |i| {
                    let x = images[grade][i];
                    (x.target, x.phase)
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let moved = moved.iter().collect::<Vec<_>>();
        let forward = matrix::coordinates(&mut self.field, &target, &moved)?;
        let backward = matrix::coordinates(&mut self.field, &moved, &target)?;
        let identity = Matrix::identity(source.len())?;
        if forward.multiply(&mut self.field, &backward)? != identity
            || backward.multiply(&mut self.field, &forward)? != identity
        {
            return Err(Error::plain(
                "boundary recoupling matrices fail the inverse test",
            ));
        }
        let left = norms(&mut self.field, &source)?;
        let right = norms(&mut self.field, &target)?;
        for (i, v) in right.iter().enumerate() {
            for (j, u) in left.iter().enumerate() {
                let adjoint = self.field.conjugate(forward.get(i, j))?;
                if self.field.multiply(&adjoint, v)?
                    != self.field.multiply(backward.get(j, i), u)?
                {
                    return Err(Error::plain(
                        "boundary recoupling violates the path Gram metrics",
                    ));
                }
            }
        }
        Ok(Transition {
            matrices: [forward, backward],
            norms: [left, right],
            grades: images,
        })
    }
}

fn norms(field: &mut Field<'_>, paths: &[&Matrix]) -> Result<Vec<Rational>> {
    paths.iter().map(|m| m.inner(field, m)).collect()
}
