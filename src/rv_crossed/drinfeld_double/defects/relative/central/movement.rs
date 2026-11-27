use crate::error::{Error, Result};
use crate::model::Model;

use super::super::tensor::Layout;
use super::super::{tree, Ambient, Balanced, Image, Module, Pair};
use super::branch::{bulk_image, Bulk};

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(in crate::rv_crossed::drinfeld_double::defects::relative) struct Map(
    pub(in crate::rv_crossed::drinfeld_double::defects::relative) Vec<Vec<Image>>,
);

impl Map {
    fn build(
        model: &mut Model,
        ambient: &Ambient<'_>,
        from: Layout<'_>,
        to: Layout<'_>,
        raw_map: impl Fn(Pair) -> Result<Image>,
    ) -> Result<Self> {
        if from.product.module.fibers != to.product.module.fibers {
            return Err(Error::plain("central defect map has unequal graded fibers"));
        }
        let mut map = Vec::with_capacity(ambient.order());
        for basis in &from.product.bases {
            model.spend(1 + basis.len())?;
            let mut used = vec![false; basis.len()];
            let row = basis
                .iter()
                .map(|&pair| {
                    let image = raw_map(pair)?;
                    if image.target >= used.len()
                        || std::mem::replace(&mut used[image.target], true)
                    {
                        return Err(Error::plain("central defect map is not a graded bijection"));
                    }
                    Ok(image)
                })
                .collect::<Result<Vec<_>>>()?;
            map.push(row);
        }
        let [l, r] = from.factors;
        let raw_size = l
            .dimension()
            .checked_mul(r.dimension())
            .and_then(|n| n.checked_mul(4))
            .ok_or_else(|| Error::plain("central defect map exceeds the raw tensor bound"))?;
        model.spend(raw_size)?;
        for (u, &a) in l.fibers.iter().enumerate() {
            for (v, &b) in r.fibers.iter().enumerate() {
                let g = ambient.mul(u, v);
                for i in 0..a {
                    for j in 0..b {
                        let pair = Pair {
                            grades: [u, v],
                            columns: [i, j],
                        };
                        let x = from.reduce(ambient, pair, 0)?;
                        let y = raw_map(pair)?;
                        let check = map[g][x.target];
                        if check.target != y.target
                            || ambient.sum([check.phase, x.phase]) != y.phase
                        {
                            return Err(Error::plain(
                                "central defect map fails to descend through balancing",
                            ));
                        }
                    }
                }
            }
        }
        tree::check_actions(
            model,
            ambient,
            &from.product.module,
            &to.product.module,
            &map,
        )?;
        Ok(Self(map))
    }

    pub(super) fn statistics(&self) -> (usize, usize) {
        self.0.iter().fold((0, 0), |(p, m), row| {
            (
                p + row.iter().filter(|x| x.phase != 0).count(),
                m + row
                    .iter()
                    .enumerate()
                    .filter(|(i, x)| *i != x.target)
                    .count(),
            )
        })
    }

    pub(super) fn followed(&self, ambient: &Ambient<'_>, next: &Self) -> Self {
        Self(
            self.0
                .iter()
                .zip(&next.0)
                .map(|(from, to)| {
                    from.iter()
                        .map(|x| {
                            let y = to[x.target];
                            Image {
                                target: y.target,
                                phase: ambient.sum([x.phase, y.phase]),
                            }
                        })
                        .collect()
                })
                .collect(),
        )
    }

    pub(super) fn inverse(&self, ambient: &Ambient<'_>) -> Self {
        Self(
            self.0
                .iter()
                .map(|row| {
                    let mut result = vec![
                        Image {
                            target: 0,
                            phase: 0
                        };
                        row.len()
                    ];
                    for (i, x) in row.iter().enumerate() {
                        result[x.target] = Image {
                            target: i,
                            phase: ambient.neg(x.phase),
                        };
                    }
                    result
                })
                .collect(),
        )
    }
}

pub(in crate::rv_crossed::drinfeld_double::defects::relative) struct Central {
    pub(in crate::rv_crossed::drinfeld_double::defects::relative) products: [Balanced; 2],
    pub(in crate::rv_crossed::drinfeld_double::defects::relative) maps: [Map; 2],
}

impl Central {
    pub(in crate::rv_crossed::drinfeld_double::defects::relative) fn new(
        model: &mut Model,
        ambient: &Ambient<'_>,
        bulk: [&Bulk; 2],
        defect: &Module,
    ) -> Result<Self> {
        let [bh, bk] = bulk;
        let g = bh.flux;
        let a = bh.simple;
        let source = Balanced::new(model, ambient, &bh.module, defect)?;
        let target = Balanced::new(model, ambient, defect, &bk.module)?;
        let from = Layout {
            product: &source,
            factors: [&bh.module, defect],
        };
        let to = Layout {
            product: &target,
            factors: [defect, &bk.module],
        };
        // c_{x,Z} is the projective action pi_Z(x). Moving Z past the
        // defect to its right therefore uses the inverse half-braiding.
        let forward = Map::build(model, ambient, from, to, |pair| {
            let [u, z] = pair.grades;
            let [i, j] = pair.columns;
            let h = ambient.mul(ambient.inv(g), u);
            let hz = ambient.mul(h, z);
            let moved = defect.grade_action(ambient, 0, h, z, j);
            let pi = bulk_image(ambient, a, hz, i, true);
            to.reduce(
                ambient,
                Pair {
                    grades: [hz, g],
                    columns: [moved.target, pi.target],
                },
                ambient.sum([ambient.alpha(g, h, z), moved.phase, pi.phase]),
            )
        })?;
        let inverse = Map::build(model, ambient, to, from, |pair| {
            let [z, v] = pair.grades;
            let [j, i] = pair.columns;
            let k = ambient.mul(ambient.inv(g), v);
            let zk = ambient.mul(z, k);
            let moved = defect.grade_action(ambient, 1, k, z, j);
            let pi = bulk_image(ambient, a, z, i, false);
            from.reduce(
                ambient,
                Pair {
                    grades: [g, zk],
                    columns: [pi.target, moved.target],
                },
                ambient.sum([
                    ambient.neg(ambient.alpha(z, g, k)),
                    ambient.alpha(g, z, k),
                    moved.phase,
                    pi.phase,
                ]),
            )
        })?;
        if forward.inverse(ambient) != inverse {
            return Err(Error::plain("central defect maps fail the inverse test"));
        }
        Ok(Self {
            products: [source, target],
            maps: [forward, inverse],
        })
    }

    pub(super) fn whisker(
        &self,
        model: &mut Model,
        ambient: &Ambient<'_>,
        from: Layout<'_>,
        to: Layout<'_>,
        side: usize,
    ) -> Result<Map> {
        if *from.factors[side] != self.products[0].module
            || *to.factors[side] != self.products[1].module
            || from.factors[1 - side] != to.factors[1 - side]
        {
            return Err(Error::plain(
                "central whiskering changed an unaffected tensor factor",
            ));
        }
        Map::build(model, ambient, from, to, |mut pair| {
            let moved = self.maps[0].0[pair.grades[side]][pair.columns[side]];
            pair.columns[side] = moved.target;
            to.reduce(ambient, pair, moved.phase)
        })
    }
}
