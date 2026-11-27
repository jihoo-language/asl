use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_lens::{self, Rational};
use crate::syntax::{Label, ParsedExpr};

use super::super::super::matrix::{self, Field, Matrix};
use super::{check_categories, Ambient, Category, Double, Image, Module};

mod channels;
use super::tree::{self, Interval, Tree};
use channels::{Channel, Transition, Work};
pub(crate) const TRIPLE_INTERVALS: [Interval; 2] = [[0, 2], [1, 3]];
pub(crate) const FOUR_INTERVALS: [Interval; 5] = [[0, 2], [1, 3], [2, 4], [0, 3], [1, 4]];

pub(crate) enum Command {
    Data(Box<ParsedExpr>),
    Pentagon([Label; 5], Box<[ParsedExpr; 5]>),
}

pub(crate) enum Expression {
    Assemble([Label; 2], Box<[ParsedExpr; 4]>),
    Entry([usize; 2], bool, Box<ParsedExpr>),
    Norm([usize; 2], Box<ParsedExpr>),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct SixJ {
    conductor: usize,
    intermediate: [Rc<str>; 2],
    objects: [(Rc<str>, usize); 4],
    paths: [Vec<[usize; 3]>; 2],
    norms: [Vec<Rational>; 2],
    matrices: [Matrix; 2],
}

pub(crate) struct Problem<'a> {
    pub(crate) double: &'a Double,
    pub(crate) sources: &'a [&'a Category],
    pub(crate) intermediates: &'a [(Interval, &'a Category)],
    pub(crate) target: &'a Category,
}

impl<'a> Problem<'a> {
    pub(crate) fn check(&self) -> Result<()> {
        check_categories(self.sources, self.target)?;
        for &([l, r], category) in self.intermediates {
            check_categories(&self.sources[l..r], category)?;
        }
        Ok(())
    }

    fn categories(&self) -> BTreeMap<Interval, &'a Category> {
        let mut result = self
            .intermediates
            .iter()
            .copied()
            .collect::<BTreeMap<_, _>>();
        result.insert([0, self.sources.len()], self.target);
        for (i, &cat) in self.sources.iter().enumerate() {
            result.insert([i, i + 1], cat);
        }
        result
    }

    pub(crate) fn conductor(&self, model: &mut Model) -> Result<usize> {
        self.check()?;
        let cats = self.categories().into_values().collect::<Vec<_>>();
        Ok(Ambient::new(model, self.double, &cats, [])?.conductor)
    }

    fn work<'b>(&self, model: &'b mut Model) -> Result<Work<'b>>
    where
        'a: 'b,
    {
        self.check()?;
        let table = self.categories();
        let cats = table.values().copied().collect::<Vec<_>>();
        let ambient = Ambient::new(model, self.double, &cats, [])?;
        let field = Field::new(model, ambient.conductor);
        Ok(Work::new(field, ambient, table))
    }

    fn inputs(&self, work: &mut Work<'_>, simple_indices: &[usize]) -> Result<Vec<Module>> {
        self.sources
            .iter()
            .zip(simple_indices)
            .map(|(cat, &index)| work.ambient.induce(work.field.model, cat, index))
            .collect()
    }

    pub(crate) fn six_j(
        &self,
        model: &mut Model,
        intermediate: [Rc<str>; 2],
        objects: [(Rc<str>, usize); 4],
    ) -> Result<Value> {
        let mut work = self.work(model)?;
        let indices = objects.each_ref().map(|(_, i)| *i);
        let [a, b, c, output] = indices;
        let inputs = self.inputs(&mut work, &[a, b, c])?;
        let grade = self.target.simple(output)?.simple.flux;
        let [l, r] = tree::TRIPLE
            .each_ref()
            .map(|shape| Tree::build(work.field.model, &work.ambient, &inputs, shape));
        let [l, r] = [l?, r?];
        let left = work.basis(&l, &indices[..3], Some(output))?;
        let right = work.basis(&r, &indices[..3], Some(output))?;
        let transition =
            work.transition(&inputs, [(&l, &left), (&r, &right)], [0, 1, 2, 3], grade)?;
        let value = SixJ {
            conductor: work.ambient.conductor,
            intermediate,
            objects,
            paths: [
                left.iter().map(Channel::triple_path).collect(),
                right.iter().map(Channel::triple_path).collect(),
            ],
            norms: transition.norms,
            matrices: transition.matrices,
        };
        value.persist(work.field.model)
    }

    pub(crate) fn pentagon(&self, model: &mut Model, indices: [usize; 5]) -> Result<String> {
        let mut work = self.work(model)?;
        let output = indices[4];
        let grade = self.target.simple(output)?.simple.flux;
        let inputs = self.inputs(&mut work, &indices[..4])?;
        let mut trees = Vec::new();
        let mut bases = Vec::new();
        for shape in &tree::FOUR {
            let tree = Tree::build(work.field.model, &work.ambient, &inputs, shape)?;
            bases.push(work.basis(&tree, &indices[..4], Some(output))?);
            trees.push(tree);
        }
        // (((ab)c)d) -> ((a(bc))d) -> (a((bc)d)) -> (a(b(cd)))
        // (((ab)c)d) -> ((ab)(cd)) -----------------> (a(b(cd)))
        let mut edges = Vec::new();
        for ([from, to], rotation) in [
            ([0, 1], [0, 1, 2, 3]),
            ([1, 2], [0, 1, 3, 4]),
            ([2, 3], [1, 2, 3, 4]),
            ([0, 4], [0, 2, 3, 4]),
            ([4, 3], [0, 1, 2, 4]),
        ] {
            edges.push(work.transition(
                &inputs,
                [(&trees[from], &bases[from]), (&trees[to], &bases[to])],
                rotation,
                grade,
            )?);
        }
        let [ab, bc, cd, ae, ed] =
            <&[Transition; 5]>::try_from(edges.as_slice()).expect("five pentagon edges");
        let first_two = bc.matrices[0].multiply(&mut work.field, &ab.matrices[0])?;
        let three = cd.matrices[0].multiply(&mut work.field, &first_two)?;
        let two = ed.matrices[0].multiply(&mut work.field, &ae.matrices[0])?;
        if three != two {
            return Err(Error::plain(
                "boundary recoupling pentagon fails on the multiplicity path bases",
            ));
        }
        for (g, row) in ab.grades.iter().enumerate() {
            work.field.model.spend(row.len())?;
            for i in 0..row.len() {
                let a = follow(&work.ambient, g, i, [ab, bc, cd]);
                let b = follow(&work.ambient, g, i, [ae, ed]);
                if a != b {
                    return Err(Error::plain(
                        "boundary associator pentagon fails on a balanced homogeneous column",
                    ));
                }
            }
        }
        let dimensions = bases.iter().map(Vec::len).collect::<Vec<_>>();
        let mixed_edges = edges
            .iter()
            .filter(|edge| {
                let f = &edge.matrices[0];
                (0..f.rows).any(|row| {
                    (0..f.cols)
                        .filter(|&col| !f.get(row, col).is_zero())
                        .count()
                        > 1
                })
            })
            .count();
        Ok(format!(
            "rvDoubleDefectPentagon = dimensions={dimensions:?}, conductor={}, mixing_edges={mixed_edges}, verified=true",
            work.ambient.conductor,
        ))
    }
}

fn follow<const N: usize>(
    ambient: &Ambient<'_>,
    g: usize,
    i: usize,
    edges: [&Transition; N],
) -> Image {
    edges.into_iter().fold(
        Image {
            target: i,
            phase: 0,
        },
        |a, edge| {
            let b = edge.grades[g][a.target];
            Image {
                target: b.target,
                phase: ambient.sum([a.phase, b.phase]),
            }
        },
    )
}

impl SixJ {
    fn persist(self, model: &mut Model) -> Result<Value> {
        let cells = 16
            + self.matrices.iter().map(Matrix::cells).sum::<usize>()
            + self.paths.iter().map(|row| row.len() * 4).sum::<usize>()
            + self
                .norms
                .iter()
                .flatten()
                .map(|x| 3 + 2 * x.numerator.len())
                .sum::<usize>();
        model.allocate(cells)?;
        Ok(Value::DoubleDefectSixJ(Rc::new(self)))
    }

    pub(crate) fn report(&self) -> String {
        let [ab, bc] = &self.intermediate;
        let [l, r] = &self.paths;
        let [(a, ia), (b, ib), (c, ic), (d, id)] = &self.objects;
        format!(
            "rvDoubleDefectSixJData = dimension={}, conductor={}, objects=[{a}:{ia},{b}:{ib},{c}:{ic},{d}:{id}], left={ab}{l:?}, right={bc}{r:?}",
            l.len(), self.conductor,
        )
    }

    pub(crate) fn entry(&self, model: &mut Model, inverse: bool, key: [usize; 2]) -> Result<Value> {
        rv_lens::rational_value(model, self.matrices[usize::from(inverse)].entry(key)?)
    }

    pub(crate) fn norm(&self, model: &mut Model, [side, index]: [usize; 2]) -> Result<Value> {
        let value = self
            .norms
            .get(side)
            .and_then(|v| v.get(index))
            .ok_or_else(|| {
                Error::plain("boundary recoupling norm side or path index is out of range")
            })?;
        rv_lens::rational_value(model, value.clone())
    }
}
