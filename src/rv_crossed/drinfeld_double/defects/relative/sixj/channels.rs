use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::rv_lens::Rational;

use super::super::tree::{self, Branch, Node, Tree};
use super::super::{decompose, Ambient, Balanced, Hom, Image, Module, Pair};
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

struct GradedMap {
    grades: Vec<Matrix>,
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
        let mut result = Matrix::zero(
            branch.product.bases[grade].len(),
            space.product.bases[grade].len(),
        )?;
        for (col, &pair) in space.product.bases[grade].iter().enumerate() {
            self.tensor_terms(branch, inclusions, pair, |field, row, term| {
                result.set(row, col, field.add(result.get(row, col), &term)?);
                Ok(())
            })?;
        }
        // Independently map the noncanonical pairs in this marked grade.
        // Grade transport supplies all other grades.
        for u in 0..self.ambient.order() {
            let v = self.ambient.mul(self.ambient.inv(u), grade);
            for i in 0..space.factors[0].fibers[u] {
                for j in 0..space.factors[1].fibers[v] {
                    let pair = Pair {
                        grades: [u, v],
                        columns: [i, j],
                    };
                    let image = space.product.canonical(
                        &self.ambient,
                        &space.factors[0],
                        &space.factors[1],
                        pair,
                        0,
                    )?;
                    if pair == space.product.bases[grade][image.target] && image.phase == 0 {
                        continue;
                    }
                    self.field.model.spend(result.rows)?;
                    let mut column = vec![Rational::zero(); result.rows];
                    self.tensor_terms(branch, inclusions, pair, |field, row, term| {
                        column[row] = field.add(&column[row], &term)?;
                        Ok(())
                    })?;
                    for (row, x) in column.iter().enumerate() {
                        if *x
                            != self
                                .field
                                .shift(result.get(row, image.target), image.phase)?
                        {
                            return Err(Error::plain(
                                "boundary fusion inclusion fails to descend through balancing",
                            ));
                        }
                    }
                }
            }
        }
        Ok(result)
    }

    fn tensor_terms(
        &mut self,
        branch: &Branch<'_>,
        inclusions: [&Channel; 2],
        pair: Pair,
        mut term: impl FnMut(&mut Field<'_>, usize, Rational) -> Result<()>,
    ) -> Result<()> {
        let [left, right] = branch.children.each_ref().map(Tree::module);
        let [u, v] = pair.grades;
        let [i, j] = pair.columns;
        let a = &inclusions[0].map.grades[u];
        let b = &inclusions[1].map.grades[v];
        self.field.model.spend(a.rows * b.rows)?;
        for s in 0..a.rows {
            let p = a.get(s, i);
            if p.is_zero() {
                continue;
            }
            for t in 0..b.rows {
                let q = b.get(t, j);
                if q.is_zero() {
                    continue;
                }
                let scale = self.field.multiply(p, q)?;
                let image = branch.product.canonical(
                    &self.ambient,
                    left,
                    right,
                    Pair {
                        grades: [u, v],
                        columns: [s, t],
                    },
                    0,
                )?;
                let coefficient = self.field.shift(&scale, image.phase)?;
                term(&mut self.field, image.target, coefficient)?;
            }
        }
        Ok(())
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

impl GradedMap {
    fn extend(
        field: &mut Field<'_>,
        ambient: &Ambient<'_>,
        source: &Module,
        target: &Module,
        representative: usize,
        base: &Matrix,
    ) -> Result<Self> {
        let mut sections = vec![None; ambient.order()];
        let [left, right] = source.ends.map(|i| &ambient.boundaries[i].members);
        field.model.spend(left.len() * right.len())?;
        for &h in left {
            for &k in right {
                let g = ambient.mul(ambient.mul(h, representative), k);
                sections[g].get_or_insert([h, k]);
            }
        }
        let mut grades = Vec::with_capacity(ambient.order());
        for (g, section) in sections.into_iter().enumerate() {
            let (rows, cols) = (target.fibers[g], source.fibers[g]);
            field.model.spend(1 + rows * cols)?;
            let mut map = Matrix::zero(rows, cols)?;
            if let Some(section) = section {
                for i in 0..base.rows {
                    let x = transport(ambient, target, representative, i, section);
                    for j in 0..base.cols {
                        let value = base.get(i, j);
                        if value.is_zero() {
                            continue;
                        }
                        let y = transport(ambient, source, representative, j, section);
                        let exponent = ambient.sum([x.phase, ambient.neg(y.phase)]);
                        map.set(x.target, y.target, field.shift(value, exponent)?);
                    }
                }
            }
            grades.push(map);
        }
        for (side, &boundary) in source.ends.iter().enumerate() {
            for &x in &ambient.boundaries[boundary].members {
                for (g, m) in grades.iter().enumerate() {
                    let next = &grades[ambient.mul(g, x)];
                    field.model.spend(m.rows * m.cols)?;
                    for i in 0..m.rows {
                        for j in 0..m.cols {
                            let u = target.grade_action(ambient, side, x, g, i);
                            let v = source.grade_action(ambient, side, x, g, j);
                            let phase = ambient.sum([u.phase, ambient.neg(v.phase)]);
                            if field.shift(m.get(i, j), phase)? != *next.get(u.target, v.target) {
                                return Err(Error::plain("boundary fusion inclusion does not intertwine the induced grade actions"));
                            }
                        }
                    }
                }
            }
        }
        Ok(Self { grades })
    }
}

fn transport(
    ambient: &Ambient<'_>,
    module: &Module,
    g: usize,
    i: usize,
    [h, k]: [usize; 2],
) -> Image {
    let left = module.grade_action(ambient, 0, h, g, i);
    let right = module.grade_action(ambient, 1, k, ambient.mul(g, h), left.target);
    Image {
        target: right.target,
        phase: ambient.sum([left.phase, right.phase]),
    }
}
