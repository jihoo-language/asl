use crate::error::{Error, Result};
use crate::model::Model;

use super::graded::{basis_limit, Ambient, Module, MAX_BASIS};
use super::Image;

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(super) struct Pair {
    pub(super) grades: [usize; 2],
    pub(super) columns: [usize; 2],
}

pub(super) struct Balanced {
    pub(super) module: Module,
    pub(super) bases: Vec<Vec<Pair>>,
    minima: Vec<usize>,
}

impl Balanced {
    pub(super) fn new(
        model: &mut Model,
        ambient: &Ambient<'_>,
        m: &Module,
        n: &Module,
    ) -> Result<Self> {
        let middle = m.ends[1];
        if middle != n.ends[0] {
            return Err(Error::plain(
                "relative tensor requires the identical middle boundary",
            ));
        }
        let group = &ambient.boundaries[middle].members;
        let raw_size = m
            .dimension()
            .checked_mul(n.dimension())
            .ok_or_else(basis_limit)?;
        let size = raw_size / group.len();
        if size > MAX_BASIS {
            return Err(basis_limit());
        }
        model.spend(ambient.order() * group.len() + size)?;
        let minima = (0..ambient.order())
            .map(|z| group.iter().map(|&x| ambient.mul(z, x)).min().unwrap())
            .collect::<Vec<_>>();
        let mut bases = vec![Vec::new(); ambient.order()];
        for (z, &md) in m.fibers.iter().enumerate() {
            if minima[z] != z {
                continue;
            }
            for (v, &nd) in n.fibers.iter().enumerate() {
                let basis = &mut bases[ambient.mul(z, v)];
                for i in 0..md {
                    for j in 0..nd {
                        basis.push(Pair {
                            grades: [z, v],
                            columns: [i, j],
                        });
                    }
                }
            }
        }
        let fibers = bases.iter().map(Vec::len).collect::<Vec<_>>();
        if fibers.iter().sum::<usize>() * group.len() != raw_size {
            return Err(Error::plain(
                "balanced tensor fails the fiber-dimension identity",
            ));
        }
        let ends = [m.ends[0], n.ends[1]];
        let reducer = Reducer {
            bases: &bases,
            minima: &minima,
        };
        let module = Module::build(model, ambient, ends, fibers, |side, x, grade, column| {
            let Pair {
                grades: [z, v],
                columns: [i, j],
            } = bases[grade][column];
            let (pair, phase) = if side == 0 {
                let moved = m.grade_action(ambient, 0, x, z, i);
                (
                    Pair {
                        grades: [ambient.mul(x, z), v],
                        columns: [moved.target, j],
                    },
                    ambient.sum([moved.phase, ambient.neg(ambient.alpha(x, z, v))]),
                )
            } else {
                let moved = n.grade_action(ambient, 1, x, v, j);
                (
                    Pair {
                        grades: [z, ambient.mul(v, x)],
                        columns: [i, moved.target],
                    },
                    ambient.sum([moved.phase, ambient.alpha(z, v, x)]),
                )
            };
            reducer.canonical(ambient, m, n, pair, phase)
        })?;
        let result = Self {
            module,
            bases,
            minima,
        };
        result.verify_balancing(model, ambient, m, n, raw_size)?;
        Ok(result)
    }

    pub(super) fn canonical(
        &self,
        ambient: &Ambient<'_>,
        m: &Module,
        n: &Module,
        pair: Pair,
        phase: usize,
    ) -> Result<Image> {
        Reducer {
            bases: &self.bases,
            minima: &self.minima,
        }
        .canonical(ambient, m, n, pair, phase)
    }

    fn verify_balancing(
        &self,
        model: &mut Model,
        ambient: &Ambient<'_>,
        m: &Module,
        n: &Module,
        raw_size: usize,
    ) -> Result<()> {
        let middle = &ambient.boundaries[m.ends[1]].members;
        model.spend(raw_size * middle.len() * (1 + m.fibers.iter().copied().max().unwrap_or(0)))?;
        for (z, &md) in m.fibers.iter().enumerate() {
            for (v, &nd) in n.fibers.iter().enumerate() {
                for i in 0..md {
                    for j in 0..nd {
                        for &x in middle {
                            let p = m.grade_action(ambient, 1, x, z, i);
                            let q = n.grade_action(ambient, 0, x, v, j);
                            let left = self.canonical(
                                ambient,
                                m,
                                n,
                                Pair {
                                    grades: [ambient.mul(z, x), v],
                                    columns: [p.target, j],
                                },
                                p.phase,
                            )?;
                            let right = self.canonical(
                                ambient,
                                m,
                                n,
                                Pair {
                                    grades: [z, ambient.mul(x, v)],
                                    columns: [i, q.target],
                                },
                                ambient.sum([ambient.alpha(z, x, v), q.phase]),
                            )?;
                            if left != right {
                                return Err(Error::plain(
                                    "balanced defect relation fails on a homogeneous basis pair",
                                ));
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

struct Reducer<'a> {
    bases: &'a [Vec<Pair>],
    minima: &'a [usize],
}

impl Reducer<'_> {
    fn canonical(
        &self,
        ambient: &Ambient<'_>,
        m: &Module,
        n: &Module,
        pair: Pair,
        phase: usize,
    ) -> Result<Image> {
        let Pair {
            grades: [z, v],
            columns: [i, j],
        } = pair;
        let q = self.minima[z];
        let x = ambient.mul(z, ambient.inv(q));
        let i = m.inverse_right(ambient, x, q, i)?;
        let j = n.grade_action(ambient, 0, x, v, j);
        let corrected = Pair {
            grades: [q, ambient.mul(v, x)],
            columns: [i.target, j.target],
        };
        let grade = ambient.mul(z, v);
        let target = self.bases[grade].binary_search(&corrected).map_err(|_| {
            Error::plain("balanced tensor has no canonical representative of a basis pair")
        })?;
        Ok(Image {
            target,
            phase: ambient.sum([phase, i.phase, j.phase, ambient.alpha(q, x, v)]),
        })
    }
}
