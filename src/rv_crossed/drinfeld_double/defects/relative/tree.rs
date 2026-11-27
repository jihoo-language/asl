use std::borrow::Borrow;

use crate::error::{Error, Result};
use crate::model::Model;

use super::{Ambient, Balanced, Image, Module, Pair};

pub(super) type Interval = [usize; 2];

pub(super) enum Shape {
    Leaf(usize),
    Join(&'static Shape, &'static Shape),
}

const A: Shape = Shape::Leaf(0);
const B: Shape = Shape::Leaf(1);
const C: Shape = Shape::Leaf(2);
const D: Shape = Shape::Leaf(3);
const AB: Shape = Shape::Join(&A, &B);
const BC: Shape = Shape::Join(&B, &C);
const CD: Shape = Shape::Join(&C, &D);
const ABC_L: Shape = Shape::Join(&AB, &C);
const ABC_R: Shape = Shape::Join(&A, &BC);
const BCD_L: Shape = Shape::Join(&BC, &D);
const BCD_R: Shape = Shape::Join(&B, &CD);

pub(super) const TRIPLE: [Shape; 2] = [ABC_L, ABC_R];
// (((ab)c)d), ((a(bc))d), (a((bc)d)), (a(b(cd))), ((ab)(cd))
pub(super) const FOUR: [Shape; 5] = [
    Shape::Join(&ABC_L, &D),
    Shape::Join(&ABC_R, &D),
    Shape::Join(&A, &BCD_L),
    Shape::Join(&A, &BCD_R),
    Shape::Join(&AB, &CD),
];

pub(super) struct Raw {
    pub(super) grades: [usize; 4],
    columns: [usize; 4],
}

pub(super) struct Tree<'a> {
    pub(super) interval: Interval,
    pub(super) node: Node<'a>,
}

pub(super) enum Node<'a> {
    Leaf(usize, &'a Module),
    Join(Box<Branch<'a>>),
}

pub(super) struct Branch<'a> {
    pub(super) children: [Tree<'a>; 2],
    pub(super) product: Balanced,
}

impl<'a> Tree<'a> {
    pub(super) fn build(
        model: &mut Model,
        ambient: &Ambient<'_>,
        inputs: &'a [impl Borrow<Module>],
        shape: &Shape,
    ) -> Result<Self> {
        let (interval, node) = match *shape {
            Shape::Leaf(i) => ([i, i + 1], Node::Leaf(i, inputs[i].borrow())),
            Shape::Join(l, r) => {
                let left = Self::build(model, ambient, inputs, l)?;
                let right = Self::build(model, ambient, inputs, r)?;
                let interval = [left.interval[0], right.interval[1]];
                let product = Balanced::new(model, ambient, left.module(), right.module())?;
                (
                    interval,
                    Node::Join(Box::new(Branch {
                        children: [left, right],
                        product,
                    })),
                )
            }
        };
        Ok(Self { interval, node })
    }

    pub(super) fn module(&self) -> &Module {
        match &self.node {
            Node::Leaf(_, module) => module,
            Node::Join(branch) => &branch.product.module,
        }
    }

    fn expand(&self, grade: usize, column: usize, raw: &mut Raw) {
        match &self.node {
            Node::Leaf(i, _) => {
                raw.grades[*i] = grade;
                raw.columns[*i] = column;
            }
            Node::Join(branch) => {
                let Pair { grades, columns } = branch.product.bases[grade][column];
                for (i, child) in branch.children.iter().enumerate() {
                    child.expand(grades[i], columns[i], raw);
                }
            }
        }
    }

    fn reduce(&self, ambient: &Ambient<'_>, raw: &Raw) -> Result<(usize, Image)> {
        match &self.node {
            Node::Leaf(i, _) => Ok((
                raw.grades[*i],
                Image {
                    target: raw.columns[*i],
                    phase: 0,
                },
            )),
            Node::Join(branch) => {
                let [l, r] = &branch.children;
                let (u, a) = l.reduce(ambient, raw)?;
                let (v, b) = r.reduce(ambient, raw)?;
                let image = branch.product.canonical(
                    ambient,
                    l.module(),
                    r.module(),
                    Pair {
                        grades: [u, v],
                        columns: [a.target, b.target],
                    },
                    ambient.sum([a.phase, b.phase]),
                )?;
                Ok((ambient.mul(u, v), image))
            }
        }
    }
}

// Three blocks identified by their four consecutive endpoints. The local
// associator multiplies by alpha on these block grades, with identity on
// spectator leaves. Comparing balanced matrices also checks naturality.
pub(super) fn rotation(
    model: &mut Model,
    ambient: &Ambient<'_>,
    inputs: &[impl Borrow<Module>],
    from: &Tree<'_>,
    to: &Tree<'_>,
    [a, b, c, d]: [usize; 4],
) -> Result<Vec<Vec<Image>>> {
    let local_phase = |raw: &Raw| {
        let word = |i, j| {
            raw.grades[i..j]
                .iter()
                .copied()
                .fold(0, |u, v| ambient.mul(u, v))
        };
        ambient.alpha(word(a, b), word(b, c), word(c, d))
    };
    if from.module().fibers != to.module().fibers {
        return Err(Error::plain(
            "boundary recoupling rotation has unequal graded fibers",
        ));
    }
    let mut rows = Vec::with_capacity(ambient.order());
    for (g, &dim) in from.module().fibers.iter().enumerate() {
        model.spend(1 + dim)?;
        let mut row = Vec::with_capacity(dim);
        let mut seen = vec![false; dim];
        let mut raw = Raw {
            grades: [0; 4],
            columns: [0; 4],
        };
        for i in 0..dim {
            from.expand(g, i, &mut raw);
            let (grade, mut image) = to.reduce(ambient, &raw)?;
            image.phase = ambient.sum([image.phase, local_phase(&raw)]);
            if grade != g || std::mem::replace(&mut seen[image.target], true) {
                return Err(Error::plain(
                    "boundary recoupling rotation is not a graded bijection",
                ));
            }
            row.push(image);
        }
        rows.push(row);
    }
    let count = inputs
        .iter()
        .try_fold(1usize, |n, m| n.checked_mul(m.borrow().dimension()))
        .ok_or_else(|| {
            Error::plain("boundary recoupling raw tensor size exceeds the computation limit")
        })?;
    model.spend(count.checked_mul(4).ok_or_else(|| {
        Error::plain("boundary recoupling raw tensor size exceeds the computation limit")
    })?)?;
    all_raw(
        inputs,
        &mut Raw {
            grades: [0; 4],
            columns: [0; 4],
        },
        &mut |raw| {
            let (g, x) = from.reduce(ambient, raw)?;
            let (h, y) = to.reduce(ambient, raw)?;
            let mapped = rows[g][x.target];
            if g != h
                || mapped.target != y.target
                || ambient.sum([x.phase, mapped.phase]) != ambient.sum([y.phase, local_phase(raw)])
            {
                return Err(Error::plain(
                    "boundary recoupling rotation fails on a balanced relation",
                ));
            }
            Ok(())
        },
    )?;
    check_actions(model, ambient, from.module(), to.module(), &rows)?;
    Ok(rows)
}

fn all_raw(
    rest: &[impl Borrow<Module>],
    raw: &mut Raw,
    action: &mut impl FnMut(&Raw) -> Result<()>,
) -> Result<()> {
    let Some((module, rest)) = rest.split_last() else {
        return action(raw);
    };
    let i = rest.len();
    for (g, &dim) in module.borrow().fibers.iter().enumerate() {
        raw.grades[i] = g;
        for k in 0..dim {
            raw.columns[i] = k;
            all_raw(rest, raw, action)?;
        }
    }
    Ok(())
}

pub(super) fn check_actions(
    model: &mut Model,
    ambient: &Ambient<'_>,
    from: &Module,
    to: &Module,
    rows: &[Vec<Image>],
) -> Result<()> {
    for (side, &boundary) in from.ends.iter().enumerate() {
        let group = &ambient.boundaries[boundary].members;
        model.spend(group.len() * from.dimension())?;
        for &x in group {
            for (g, row) in rows.iter().enumerate() {
                for (i, iso) in row.iter().enumerate() {
                    let u = from.grade_action(ambient, side, x, g, i);
                    let v = to.grade_action(ambient, side, x, g, iso.target);
                    let w = rows[ambient.mul(g, x)][u.target];
                    if w.target != v.target
                        || ambient.sum([u.phase, w.phase]) != ambient.sum([iso.phase, v.phase])
                    {
                        return Err(Error::plain(
                            "boundary recoupling rotation changes an outer action",
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}
