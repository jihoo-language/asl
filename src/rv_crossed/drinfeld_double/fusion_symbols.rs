use std::collections::{BTreeMap, HashMap};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_lens;
use crate::syntax::Label;

use super::matrix::{self, coordinates, Field, Matrix};
use super::{negate_phase, representations, Double};

mod coherence;

const MAX_TENSOR: usize = 256;

pub(crate) enum Command {
    Space(Label, [usize; 3]),
    FData(Label, [usize; 4]),
    Pentagon(Label, [usize; 5]),
    Hexagon(Label, [usize; 4]),
}

pub(crate) enum Expression {
    Vertex(Label, [usize; 6]),
    F(Label, [usize; 6], bool),
    R(Label, [usize; 5], bool),
}

struct Work<'a> {
    double: &'a Double,
    field: Field<'a>,
    spaces: HashMap<[usize; 3], Rc<Space>>,
}

struct Space {
    rows: usize,
    cols: usize,
    maps: Vec<Matrix>,
}

#[derive(Clone, Copy, Eq, PartialEq)]
enum Tree {
    Leaf(usize),
    Join(&'static Tree, &'static Tree),
}

const A: Tree = Tree::Leaf(0);
const B: Tree = Tree::Leaf(1);
const C: Tree = Tree::Leaf(2);
const AB: Tree = Tree::Join(&A, &B);
const BC: Tree = Tree::Join(&B, &C);
const LEFT: Tree = Tree::Join(&AB, &C);
const RIGHT: Tree = Tree::Join(&A, &BC);

struct Channel {
    labels: Vec<(usize, usize)>,
    map: Matrix,
}

impl Channel {
    fn pair(&self) -> (usize, usize, usize) {
        let [(internal, i), (_, j)] = *self.labels.as_slice() else {
            unreachable!("a triple has two vertices")
        };
        (internal, i, j)
    }
}

impl Double {
    fn treeshift(&self, fluxes: [usize; 3], inverse: bool) -> usize {
        let phase = self.cocycle_phase(fluxes);
        if inverse {
            phase
        } else {
            negate_phase(self.conductor, phase)
        }
    }

    pub(crate) fn fusion_space(
        &self,
        model: &mut Model,
        name: &str,
        args @ [a, b, c]: [usize; 3],
    ) -> Result<String> {
        let space = Work::new(model, self).space(args)?;
        Ok(format!(
            "rvDoubleFusionSpace[{name},{a},{b},{c}] = multiplicity={}, rows={}, columns={}",
            space.maps.len(),
            space.rows,
            space.cols,
        ))
    }

    pub(crate) fn fusion_vertex(
        &self,
        model: &mut Model,
        [a, b, c, basis, row, col]: [usize; 6],
    ) -> Result<Value> {
        let mut work = Work::new(model, self);
        let space = work.space([a, b, c])?;
        let vertex = space
            .maps
            .get(basis)
            .ok_or_else(|| Error::plain("double trivalent channel index is out of range"))?;
        let value = vertex.entry([row, col])?;
        rv_lens::rational_value(work.field.model, value)
    }

    pub(crate) fn f_data(
        &self,
        model: &mut Model,
        name: &str,
        [a, b, c, d]: [usize; 4],
    ) -> Result<String> {
        let mut work = Work::new(model, self);
        let (l, r) = work.triple([a, b, c], d)?;
        work.f_move([a, b, c], &l, &r, false)?;
        let left = l.iter().map(Channel::pair).collect::<Vec<_>>();
        let right = r.iter().map(Channel::pair).collect::<Vec<_>>();
        Ok(format!(
            "rvDoubleFData[{name},{a},{b},{c},{d}] = dimension={}, left={left:?}, right={right:?}",
            left.len()
        ))
    }

    pub(crate) fn f_entry(
        &self,
        model: &mut Model,
        [a, b, c, d, row, col]: [usize; 6],
        inverse: bool,
    ) -> Result<Value> {
        let mut work = Work::new(model, self);
        let (left, right) = work.triple([a, b, c], d)?;
        let matrix = work.f_move([a, b, c], &left, &right, inverse)?;
        rv_lens::rational_value(work.field.model, matrix.entry([row, col])?)
    }

    pub(crate) fn r_entry(
        &self,
        model: &mut Model,
        [a, b, c, row, col]: [usize; 5],
        inverse: bool,
    ) -> Result<Value> {
        let mut work = Work::new(model, self);
        let matrix = work.r_move([a, b, c], inverse)?;
        rv_lens::rational_value(work.field.model, matrix.entry([row, col])?)
    }
}

impl<'a> Work<'a> {
    fn new(model: &'a mut Model, double: &'a Double) -> Self {
        Self {
            double,
            field: Field::new(model, double.conductor),
            spaces: HashMap::new(),
        }
    }

    fn check(&self, colors: &[usize]) -> Result<()> {
        for &x in colors {
            self.double.check_index(x)?;
        }
        Ok(())
    }

    fn space(&mut self, key: [usize; 3]) -> Result<Rc<Space>> {
        if let Some(space) = self.spaces.get(&key) {
            return Ok(space.clone());
        }
        self.check(&key)?;
        let space = Rc::new(self.solve(key)?);
        self.spaces.insert(key, space.clone());
        Ok(space)
    }

    fn solve(&mut self, [a, b, c]: [usize; 3]) -> Result<Space> {
        let d = self.double;
        let [a0, b0, c0] = [a, b, c].map(|i| &d.simples[i]);
        let (rows, cols) = (a0.dimension * b0.dimension, c0.dimension);
        let expected = d.verlinde(self.field.model, a, b, c)?;
        if d.group.mul(a0.flux, b0.flux) != c0.flux {
            if expected != 0 {
                return Err(Error::plain(
                    "double Verlinde coefficient violates the flux grading",
                ));
            }
            return Ok(Space {
                rows,
                cols,
                maps: Vec::new(),
            });
        }
        let phases = (0..d.group.order())
            .map(|x| d.tensor_phase(x, a0.flux, b0.flux))
            .collect::<Vec<_>>();
        let tensor = representations::tensor_action(
            self.field.model,
            &d.representations[a].action,
            &d.representations[b].action,
            &phases,
            d.conductor,
        )?;
        let orbits = representations::intertwiner_orbits(
            self.field.model,
            &tensor,
            &d.representations[c].action,
            d.conductor,
        )?;
        let mut maps = Vec::new();
        for orbit in orbits {
            let mut matrix = Matrix::zero(rows, cols)?;
            for (i, phase) in orbit {
                matrix.set(i / cols, i % cols, self.field.phase(phase)?);
            }
            maps.push(matrix);
        }
        if maps.len() as i128 != expected {
            return Err(Error::plain(
                "double trivalent intertwiner dimension disagrees with Verlinde",
            ));
        }
        Ok(Space { rows, cols, maps })
    }

    fn triple(&mut self, colors: [usize; 3], total: usize) -> Result<(Vec<Channel>, Vec<Channel>)> {
        self.check(&[colors[0], colors[1], colors[2], total])?;
        Ok((
            self.basis(&LEFT, &colors, total)?,
            self.basis(&RIGHT, &colors, total)?,
        ))
    }

    fn basis(&mut self, shape: &Tree, colors: &[usize], total: usize) -> Result<Vec<Channel>> {
        Ok(self
            .expand(shape, colors, Some(total))?
            .remove(&total)
            .unwrap_or_default())
    }

    fn expand(
        &mut self,
        shape: &Tree,
        colors: &[usize],
        total: Option<usize>,
    ) -> Result<BTreeMap<usize, Vec<Channel>>> {
        let mut result = BTreeMap::new();
        match shape {
            Tree::Leaf(i) => {
                let color = colors[*i];
                if total.is_none_or(|out| out == color) {
                    result.insert(
                        color,
                        vec![Channel {
                            labels: Vec::new(),
                            map: Matrix::identity(self.double.simples[color].dimension)?,
                        }],
                    );
                }
            }
            Tree::Join(left, right) => {
                let left = self.expand(left, colors, None)?;
                let right = self.expand(right, colors, None)?;
                let mut size = 0;
                for (&a, ls) in &left {
                    for (&b, rs) in &right {
                        let g = self
                            .double
                            .group
                            .mul(self.double.simples[a].flux, self.double.simples[b].flux);
                        let range = self.double.flux_offsets[g]..self.double.flux_offsets[g + 1];
                        for c in range.filter(|&c| total.is_none_or(|out| out == c)) {
                            let space = self.space([a, b, c])?;
                            if space.maps.is_empty() {
                                continue;
                            }
                            let channels = result.entry(c).or_insert_with(Vec::new);
                            for l in ls {
                                for r in rs {
                                    if l.map.rows * r.map.rows > MAX_TENSOR {
                                        return Err(Error::plain(
                                            "double coherence tensor exceeds dimension 256",
                                        ));
                                    }
                                    let lr = l.map.tensor(&mut self.field, &r.map)?;
                                    for (i, vertex) in space.maps.iter().enumerate() {
                                        size += 1;
                                        if size > MAX_TENSOR {
                                            return Err(Error::plain(
                                                "double coherence exceeds its channel limit of 256",
                                            ));
                                        }
                                        let map = lr.multiply(&mut self.field, vertex)?;
                                        let mut labels =
                                            Vec::with_capacity(l.labels.len() + r.labels.len() + 1);
                                        labels.extend(&l.labels);
                                        labels.extend(&r.labels);
                                        labels.push((c, i));
                                        channels.push(Channel { labels, map });
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(result)
    }

    fn move_by(&mut self, source: &[Channel], target: &[Channel], phase: usize) -> Result<Matrix> {
        let target = target.iter().map(|x| &x.map).collect::<Vec<_>>();
        let shifted = source
            .iter()
            .map(|x| x.map.transformed(&mut self.field, |i| (i, phase)))
            .collect::<Result<Vec<_>>>()?;
        let source = shifted.iter().collect::<Vec<_>>();
        coordinates(&mut self.field, &target, &source)
    }

    fn f_move(
        &mut self,
        colors: [usize; 3],
        left: &[Channel],
        right: &[Channel],
        inverse: bool,
    ) -> Result<Matrix> {
        let fluxes = colors.map(|i| self.double.simples[i].flux);
        let phase = self.double.treeshift(fluxes, inverse);
        if inverse {
            self.move_by(right, left, phase)
        } else {
            self.move_by(left, right, phase)
        }
    }

    fn r_move(&mut self, [a, b, c]: [usize; 3], inverse: bool) -> Result<Matrix> {
        let source = self.space([a, b, c])?;
        let target = self.space([b, a, c])?;
        let (x, y) = (&self.double.simples[a], &self.double.simples[b]);
        let (g, h, da, db) = (x.flux, y.flux, x.dimension, y.dimension);
        let transform = |row: usize| {
            let (i, j) = (row / db, row % db);
            if inverse {
                let image =
                    self.double.representations[a].inverse_image(h, i, self.double.conductor);
                (j * da + image.target, image.phase)
            } else {
                let image = self.double.representations[b].action[g][j];
                (image.target * da + i, image.phase)
            }
        };
        let maps = source
            .maps
            .iter()
            .map(|m| m.transformed(&mut self.field, transform))
            .collect::<Result<Vec<_>>>()?;
        coordinates(
            &mut self.field,
            &target.maps.iter().collect::<Vec<_>>(),
            &maps.iter().collect::<Vec<_>>(),
        )
    }
}
