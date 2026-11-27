use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};

use super::super::super::super::Sector;
use super::super::{
    decompose, sector_counts, Ambient, Balanced, Category, Double, Image, Module, Pair, Relative,
};

pub(super) struct Bulk {
    pub(super) simple: usize,
    pub(super) flux: usize,
    pub(super) module: Module,
}

pub(crate) fn self_boundary(category: &Category) -> Result<()> {
    if !Rc::ptr_eq(&category.boundaries[0], &category.boundaries[1]) {
        return Err(Error::plain(
            "bulk boundary induction requires a nominal self-boundary category",
        ));
    }
    Ok(())
}

impl Bulk {
    pub(super) fn new(
        model: &mut Model,
        ambient: &Ambient<'_>,
        category: &Category,
        simple: usize,
    ) -> Result<Self> {
        self_boundary(category)?;
        let color = ambient.double.simple(simple)?;
        let g = color.flux;
        let ends = ambient.ends(category);
        let h = &ambient.boundaries[ends[0]].members;
        let mut fibers = vec![0; ambient.order()];
        for &x in h {
            fibers[ambient.mul(g, x)] = color.dimension;
        }
        let module = Module::build(model, ambient, ends, fibers, |side, x, z, i| {
            let a = ambient.mul(ambient.inv(g), z);
            if side == 0 {
                let pi = bulk_image(ambient, simple, x, i, false);
                Ok(Image {
                    target: pi.target,
                    phase: ambient.sum([
                        pi.phase,
                        ambient.neg(ambient.alpha(x, g, a)),
                        ambient.alpha(g, x, a),
                        ambient.mu(ends[0], x, a),
                    ]),
                })
            } else {
                Ok(Image {
                    target: i,
                    phase: ambient.sum([ambient.alpha(g, a, x), ambient.mu(ends[1], a, x)]),
                })
            }
        })?;
        let result = Self {
            simple,
            flux: g,
            module,
        };
        result.check_transport(model, ambient, category)?;
        Ok(result)
    }

    fn check_transport(
        &self,
        model: &mut Model,
        ambient: &Ambient<'_>,
        category: &Category,
    ) -> Result<()> {
        let g = self.flux;
        let representative = category.cosets[category.coset_index[g]].support[0];
        for r in [g, representative] {
            let h = ambient.mul(g, ambient.inv(r));
            let action = self.module.transport(model, ambient, category, r)?;
            for (&x, row) in category.intersection.iter().zip(&action) {
                let b = self.module.ends[0];
                let mixed = ambient.double.tensor_phase(x, r, h)
                    * (ambient.conductor / ambient.double.conductor);
                let gauge = ambient.sum([
                    ambient.neg(ambient.mu(b, x, h)),
                    ambient.mu(b, h, x),
                    ambient.neg(mixed),
                ]);
                for (i, got) in row.iter().enumerate() {
                    let pi = bulk_image(ambient, self.simple, x, i, false);
                    if got.target != pi.target || got.phase != ambient.sum([pi.phase, gauge]) {
                        return Err(Error::plain(
                            "induced bulk branch fails its representative stabilizer gauge",
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    pub(super) fn branches(
        &self,
        model: &mut Model,
        ambient: &Ambient<'_>,
        category: &Category,
    ) -> Result<Vec<(usize, usize)>> {
        Ok(sector_counts(&decompose(
            model,
            ambient,
            category,
            &self.module,
        )?))
    }

    pub(super) fn value(
        self,
        model: &mut Model,
        ambient: &Ambient<'_>,
        category: &Category,
    ) -> Result<Value> {
        let g = self.flux;
        let bases = self
            .module
            .fibers
            .iter()
            .enumerate()
            .map(|(z, &dimension)| {
                let h = ambient.mul(ambient.inv(g), z);
                (0..dimension)
                    .map(|i| Pair {
                        grades: [g, h],
                        columns: [i, 0],
                    })
                    .collect()
            })
            .collect();
        let homs = decompose(model, ambient, category, &self.module)?;
        Relative::from_graded(ambient.conductor, self.module, bases, homs).persist(model)
    }
}

pub(super) fn bulk_image(
    ambient: &Ambient<'_>,
    simple: usize,
    x: usize,
    i: usize,
    inverse: bool,
) -> Image {
    let double = ambient.double;
    let rep = &double.representations[simple];
    let image = if inverse {
        rep.inverse_image(x, i, double.conductor)
    } else {
        rep.action[x][i]
    };
    Image {
        target: image.target,
        phase: image.phase * (ambient.conductor / double.conductor),
    }
}

struct BranchTable {
    branches: Vec<Vec<(usize, usize)>>,
    vacuum: usize,
}

impl BranchTable {
    fn new(model: &mut Model, ambient: &Ambient<'_>, category: &Category) -> Result<Self> {
        let branches = (0..ambient.double.count())
            .map(|a| Bulk::new(model, ambient, category, a)?.branches(model, ambient, category))
            .collect::<Result<Vec<_>>>()?;
        let [(vacuum, 1)] = *branches[0].as_slice() else {
            return Err(Error::plain("bulk vacuum did not induce one boundary unit"));
        };
        let unit = ambient.induce(model, category, vacuum)?;
        let trivial = Bulk::new(model, ambient, category, 0)?;
        let group = &ambient.boundaries[unit.ends[0]].members;
        for &x in group {
            for &g in group {
                for (side, p) in [
                    ambient.mu(unit.ends[0], x, g),
                    ambient.mu(unit.ends[1], g, x),
                ]
                .into_iter()
                .enumerate()
                {
                    let expected = Image {
                        target: 0,
                        phase: p,
                    };
                    if trivial.module.grade_action(ambient, side, x, g, 0) != expected
                        || unit.grade_action(ambient, side, x, g, 0) != expected
                    {
                        return Err(Error::plain(
                            "bulk vacuum branch fails the regular boundary module",
                        ));
                    }
                }
            }
        }
        // The exact algebra left/right action gives an independent unitor
        // check elsewhere; the annular vacuum comparison uses every color.
        let induced = branches
            .iter()
            .enumerate()
            .filter_map(|(a, part)| {
                let n = part.iter().find_map(|&(b, n)| (b == vacuum).then_some(n))?;
                Some((a, n as i128))
            })
            .collect::<Vec<_>>();
        if induced.as_slice() != category.boundaries[0].sector()?.0 {
            return Err(Error::plain(
                "bulk-to-boundary vacuum branching disagrees with the condensate",
            ));
        }
        Ok(Self { branches, vacuum })
    }
}

type Multiplicities = Vec<(usize, usize)>;

pub(crate) fn fusion(
    model: &mut Model,
    double: &Double,
    category: &Category,
    one: Option<[usize; 2]>,
) -> Result<String> {
    self_boundary(category)?;
    let ambient = Ambient::new(model, double, &[category], [])?;
    let table = BranchTable::new(model, &ambient, category)?;
    let mut self_fusion: BTreeMap<[usize; 2], Multiplicities> = BTreeMap::new();
    let checks: Box<dyn Iterator<Item = [usize; 2]>> = match one {
        Some([a, b]) => {
            double.check_index(a)?;
            double.check_index(b)?;
            Box::new(std::iter::once([a, b]))
        }
        None => {
            Box::new((0..double.count()).flat_map(|a| (0..double.count()).map(move |b| [a, b])))
        }
    };
    let mut count = 0;
    let mut last = Vec::new();
    for [a, b] in checks {
        let direct = double.fuse(model, &Sector(vec![(a, 1)]), &Sector(vec![(b, 1)]))?;
        let mut lhs = vec![0; category.count()];
        for (c, n) in direct.0 {
            let n = usize::try_from(n).map_err(|_| arithmetic())?;
            accumulate(&mut lhs, &table.branches[c], n)?;
        }
        let mut rhs = vec![0; category.count()];
        for &(i, u) in &table.branches[a] {
            for &(j, v) in &table.branches[b] {
                let channels = match self_fusion.entry([i, j]) {
                    std::collections::btree_map::Entry::Occupied(cached) => cached.into_mut(),
                    std::collections::btree_map::Entry::Vacant(new) => {
                        let m = ambient.induce(model, category, i)?;
                        let n = ambient.induce(model, category, j)?;
                        let product = Balanced::new(model, &ambient, &m, &n)?;
                        new.insert(sector_counts(&decompose(
                            model,
                            &ambient,
                            category,
                            &product.module,
                        )?))
                    }
                };
                accumulate(&mut rhs, channels, u.checked_mul(v).ok_or_else(arithmetic)?)?;
            }
        }
        if lhs != rhs {
            return Err(Error::plain(format!(
                "bulk branching disagrees with boundary fusion for [{a},{b}]"
            )));
        }
        count += 1;
        if one.is_some() {
            last = rhs
                .into_iter()
                .enumerate()
                .filter(|&(_, n)| n != 0)
                .collect();
        }
    }
    if let Some([a, b]) = one {
        Ok(format!(
            "rvDoubleBulkBoundaryFusion[{a},{b}] = sectors={last:?}, verified=true"
        ))
    } else {
        Ok(format!("rvDoubleBulkBoundaryValidate = bulk_simples={}, boundary_simples={}, vacuum={}, fusion_pairs={count}, verified=true",
            double.count(),category.count(),table.vacuum))
    }
}

fn arithmetic() -> Error {
    Error::plain("bulk boundary fusion multiplicity exceeds the integer bound")
}

fn accumulate(target: &mut [usize], coeff: &[(usize, usize)], weight: usize) -> Result<()> {
    for &(i, n) in coeff {
        target[i] = n
            .checked_mul(weight)
            .and_then(|term| target[i].checked_add(term))
            .ok_or_else(arithmetic)?;
    }
    Ok(())
}
