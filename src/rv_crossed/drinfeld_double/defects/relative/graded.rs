use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_cyclotomic;
use crate::rv_lens;

use super::super::super::boundaries::Boundary;
use super::super::{conductor_limit, negate_phase, Category, Double, Image};

pub(super) const MAX_BASIS: usize = 65_536;

pub(super) struct Ambient<'a> {
    pub(super) double: &'a Double,
    pub(super) conductor: usize,
    pub(super) boundaries: Vec<BoundaryTable>,
}

pub(super) struct BoundaryTable {
    source: Rc<Boundary>,
    pub(super) members: Vec<usize>,
    phases: Vec<Vec<usize>>,
}

impl Ambient<'_> {
    pub(super) fn new<'a>(
        model: &mut Model,
        double: &'a Double,
        categories: &[&Category],
        conductors: impl IntoIterator<Item = usize>,
    ) -> Result<Ambient<'a>> {
        let mut conductor = double.conductor as u128;
        for order in conductors {
            conductor = lcm(conductor, order as u128)?;
        }
        let mut sources = Vec::new();
        for &cat in categories {
            conductor = lcm(conductor, cat.conductor as u128)?;
            for boundary in &cat.boundaries {
                if sources.iter().any(|old| Rc::ptr_eq(old, boundary)) {
                    continue;
                }
                let members = boundary.members()?;
                model.spend(members.len().pow(2))?;
                for &x in members {
                    for &y in members {
                        let phase = boundary.cochain([x, y])?;
                        conductor =
                            lcm(conductor, phase.root_order().ok_or_else(conductor_limit)?)?;
                    }
                }
                sources.push(boundary.clone());
            }
        }
        let conductor = usize::try_from(conductor).map_err(|_| conductor_limit())?;
        rv_cyclotomic::check_order(conductor).map_err(|_| conductor_limit())?;
        let boundaries = sources
            .into_iter()
            .map(|source| {
                let members = source.members()?.to_vec();
                let phases = members
                    .iter()
                    .map(|&x| {
                        members
                            .iter()
                            .map(|&y| {
                                Ok(source
                                    .cochain([x, y])?
                                    .root_index(conductor)
                                    .expect("common conductor"))
                            })
                            .collect::<Result<Vec<_>>>()
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(BoundaryTable {
                    source,
                    members,
                    phases,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Ambient {
            double,
            conductor,
            boundaries,
        })
    }

    pub(super) fn ends(&self, cat: &Category) -> [usize; 2] {
        cat.boundaries.each_ref().map(|boundary| {
            self.boundaries
                .iter()
                .position(|x| Rc::ptr_eq(boundary, &x.source))
                .expect("known ambient boundary")
        })
    }
    pub(super) fn order(&self) -> usize {
        self.double.group.order()
    }
    pub(super) fn mul(&self, x: usize, y: usize) -> usize {
        self.double.group.mul(x, y)
    }
    pub(super) fn inv(&self, x: usize) -> usize {
        self.double.group.inv(x)
    }
    pub(super) fn neg(&self, p: usize) -> usize {
        negate_phase(self.conductor, p)
    }
    pub(super) fn sum(&self, terms: impl IntoIterator<Item = usize>) -> usize {
        terms.into_iter().sum::<usize>() % self.conductor
    }
    pub(super) fn alpha(&self, x: usize, y: usize, z: usize) -> usize {
        self.neg(self.double.cocycle_phase([x, y, z]) * (self.conductor / self.double.conductor))
    }

    fn mu(&self, boundary: usize, x: usize, y: usize) -> usize {
        let t = &self.boundaries[boundary];
        t.phases[t.local(x)][t.local(y)]
    }

    // Multiplier for (u,v) acting after (h,k) on grade z, where
    // P(u,v)=R_v L_u. All arguments are actual bulk group elements.
    fn couple(&self, ends: [usize; 2], [u, v]: [usize; 2], [h, k]: [usize; 2], z: usize) -> usize {
        let zh = self.mul(z, h);
        self.sum([
            self.neg(self.alpha(u, zh, k)),
            self.alpha(self.mul(zh, u), k, v),
            self.mu(ends[1], k, v),
            self.mu(ends[0], u, h),
            self.neg(self.alpha(u, h, z)),
        ])
    }

    pub(super) fn induce(&self, model: &mut Model, cat: &Category, index: usize) -> Result<Module> {
        let object = cat.simple(index)?;
        let g = object.simple.flux;
        let dim = object.simple.dimension;
        let ends = self.ends(cat);
        let mut sections = vec![None; self.order()];
        let [h, k] = ends.map(|i| self.boundaries[i].members.as_slice());
        for &u in h {
            for &v in k {
                let z = self.mul(self.mul(u, g), v);
                sections[z].get_or_insert([u, v]);
            }
        }
        if sections[g] != Some([0, 0]) {
            return Err(Error::plain(
                "induced defect section at its representative is not the identity",
            ));
        }
        let fibers = sections
            .iter()
            .map(|s| if s.is_some() { dim } else { 0 })
            .collect::<Vec<_>>();
        let scale = self.conductor / cat.conductor;
        let module = Module::build(model, self, ends, fibers, |side, a, z, i| {
            let [u, v] = if side == 0 { [a, 0] } else { [0, a] };
            let [h, k] = sections[z].expect("supported grade");
            let zz = self.mul(z, a);
            let [s, t] = sections[zz].expect("supported translated grade");
            let x = self.mul(self.mul(u, h), self.inv(s));
            let xx = self.inv(x);
            if self.mul(self.mul(v, k), self.inv(t)) != xx {
                return Err(Error::plain(
                    "defect induction section did not differ by the stabilizer",
                ));
            }
            let rep = object.representation.action[cat.local(x)?][i];
            let phase = self.sum([
                rep.phase * scale,
                self.alpha(g, x, xx),
                self.mu(ends[1], x, xx),
                self.couple(ends, [u, v], [h, k], g),
                self.neg(self.couple(ends, [s, t], [x, xx], g)),
            ]);
            Ok(Image {
                target: rep.target,
                phase,
            })
        })?;
        let witness = module.transport(model, self, cat, g)?;
        for (row, expected) in witness.iter().zip(&object.representation.action) {
            if row
                .iter()
                .zip(expected)
                .any(|(a, b)| a.target != b.target || a.phase != b.phase * scale)
            {
                return Err(Error::plain(
                    "induced bimodule changed its stabilizer on the marked fiber",
                ));
            }
        }
        Ok(module)
    }
}

impl BoundaryTable {
    pub(super) fn local(&self, x: usize) -> usize {
        self.members
            .binary_search(&x)
            .expect("known subgroup action")
    }
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(super) struct Module {
    pub(super) ends: [usize; 2],
    pub(super) fibers: Vec<usize>,
    pub(super) actions: [Vec<Vec<Vec<Image>>>; 2],
}

impl Module {
    pub(super) fn build(
        model: &mut Model,
        ambient: &Ambient<'_>,
        ends: [usize; 2],
        fibers: Vec<usize>,
        mut action: impl FnMut(usize, usize, usize, usize) -> Result<Image>,
    ) -> Result<Self> {
        let size = fibers.iter().sum::<usize>();
        if size > MAX_BASIS {
            return Err(basis_limit());
        }
        let mut actions = [Vec::new(), Vec::new()];
        for (side, rows) in actions.iter_mut().enumerate() {
            let group = &ambient.boundaries[ends[side]].members;
            model.spend(size * group.len())?;
            for &x in group {
                let mut grades = Vec::new();
                for (z, &dim) in fibers.iter().enumerate() {
                    let mut images = Vec::with_capacity(dim);
                    for i in 0..dim {
                        images.push(action(side, x, z, i)?);
                    }
                    grades.push(images);
                }
                rows.push(grades);
            }
        }
        let module = Self {
            ends,
            fibers,
            actions,
        };
        module.verify(model, ambient)?;
        Ok(module)
    }

    pub(super) fn dimension(&self) -> usize {
        self.fibers.iter().sum()
    }

    pub(super) fn grade_action(
        &self,
        ambient: &Ambient<'_>,
        side: usize,
        x: usize,
        z: usize,
        i: usize,
    ) -> Image {
        self.actions[side][ambient.boundaries[self.ends[side]].local(x)][z][i]
    }

    pub(super) fn inverse_right(
        &self,
        ambient: &Ambient<'_>,
        x: usize,
        z: usize,
        target: usize,
    ) -> Result<Image> {
        let row = &self.actions[1][ambient.boundaries[self.ends[1]].local(x)][z];
        let (column, image) = row
            .iter()
            .enumerate()
            .find(|(_, image)| image.target == target)
            .ok_or_else(|| {
                Error::plain("graded defect right action has no inverse on a basis column")
            })?;
        Ok(Image {
            target: column,
            phase: ambient.neg(image.phase),
        })
    }

    pub(super) fn transport(
        &self,
        model: &mut Model,
        ambient: &Ambient<'_>,
        cat: &Category,
        g: usize,
    ) -> Result<Vec<Vec<Image>>> {
        let dim = self.fibers[g];
        model.spend(cat.intersection.len() * dim * (1 + dim))?;
        cat.intersection
            .iter()
            .map(|&x| {
                (0..dim)
                    .map(|i| {
                        let left = self.grade_action(ambient, 0, x, g, i);
                        let back = self.inverse_right(ambient, x, g, left.target)?;
                        Ok(Image {
                            target: back.target,
                            phase: ambient.sum([left.phase, back.phase]),
                        })
                    })
                    .collect()
            })
            .collect()
    }

    pub(super) fn cells(&self) -> usize {
        8 + self.fibers.len()
            + self
                .actions
                .iter()
                .map(|side| {
                    3 * side.len()
                        + side
                            .iter()
                            .map(|row| 3 * row.len() + 2 * row.iter().map(Vec::len).sum::<usize>())
                            .sum::<usize>()
                })
                .sum::<usize>()
    }

    fn verify(&self, model: &mut Model, ambient: &Ambient<'_>) -> Result<()> {
        let [h, k] = self.ends.map(|i| ambient.boundaries[i].members.as_slice());
        model
            .spend(self.dimension() * (h.len().pow(2) + k.len().pow(2) + h.len() * k.len() + 2))?;
        for (z, &dim) in self.fibers.iter().enumerate() {
            for i in 0..dim {
                for side in 0..2 {
                    let unit = self.grade_action(ambient, side, 0, z, i);
                    if unit.target != i || unit.phase != 0 {
                        return Err(Error::plain("graded defect unit action fails"));
                    }
                }
                for &x in h {
                    let left = self.grade_action(ambient, 0, x, z, i);
                    let xz = ambient.mul(x, z);
                    for &y in h {
                        let after = self.grade_action(ambient, 0, y, xz, left.target);
                        let merged = self.grade_action(ambient, 0, ambient.mul(y, x), z, i);
                        let phase = ambient.sum([
                            merged.phase,
                            ambient.mu(self.ends[0], y, x),
                            ambient.neg(ambient.alpha(y, x, z)),
                        ]);
                        certify(
                            left.phase,
                            after,
                            merged.target,
                            phase,
                            ambient.conductor,
                            "left",
                        )?;
                    }
                    for &y in k {
                        let after = self.grade_action(ambient, 1, y, xz, left.target);
                        let right = self.grade_action(ambient, 1, y, z, i);
                        let remapped =
                            self.grade_action(ambient, 0, x, ambient.mul(z, y), right.target);
                        let phase =
                            ambient.sum([remapped.phase, right.phase, ambient.alpha(x, z, y)]);
                        certify(
                            left.phase,
                            after,
                            remapped.target,
                            phase,
                            ambient.conductor,
                            "commuting",
                        )?;
                    }
                }
                for &x in k {
                    let right = self.grade_action(ambient, 1, x, z, i);
                    let zx = ambient.mul(z, x);
                    for &y in k {
                        let after = self.grade_action(ambient, 1, y, zx, right.target);
                        let merged = self.grade_action(ambient, 1, ambient.mul(x, y), z, i);
                        let phase = ambient.sum([
                            merged.phase,
                            ambient.mu(self.ends[1], x, y),
                            ambient.alpha(z, x, y),
                        ]);
                        certify(
                            right.phase,
                            after,
                            merged.target,
                            phase,
                            ambient.conductor,
                            "right",
                        )?;
                    }
                }
            }
        }
        Ok(())
    }
}

fn certify(
    start: usize,
    after: Image,
    target: usize,
    phase: usize,
    conductor: usize,
    side: &str,
) -> Result<()> {
    if after.target != target || (start + after.phase) % conductor != phase {
        return Err(Error::plain(format!(
            "graded defect {side} action identity failed"
        )));
    }
    Ok(())
}

fn lcm(a: u128, b: u128) -> Result<u128> {
    (a / rv_lens::gcd(a, b))
        .checked_mul(b)
        .ok_or_else(conductor_limit)
}

pub(super) fn basis_limit() -> Error {
    Error::plain("relative defect tensor exceeds its balanced basis size limit")
}
