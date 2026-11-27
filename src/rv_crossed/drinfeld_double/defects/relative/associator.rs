use crate::error::{Error, Result};
use crate::model::Model;

use super::graded::{Ambient, Module};
use super::tensor::{Balanced, Pair};
use super::{check_categories, decompose, sector_counts, Category, Double, Image, Operand};

pub(crate) fn associate(
    model: &mut Model,
    double: &Double,
    target: &Category,
    inputs: [(&Category, Operand); 3],
    entry: Option<[usize; 2]>,
) -> Result<String> {
    let [(a, x), (b, y), (c, z)] = inputs;
    check_categories(&[a, b, c], target)?;
    let ambient = Ambient::new(
        model,
        double,
        &[a, b, c, target],
        [x.conductor(), y.conductor(), z.conductor()]
            .into_iter()
            .flatten(),
    )?;
    let m = x.module(model, &ambient, a)?;
    let n = y.module(model, &ambient, b)?;
    let p = z.module(model, &ambient, c)?;
    let mn = Balanced::new(model, &ambient, &m, &n)?;
    let np = Balanced::new(model, &ambient, &n, &p)?;
    let left = Balanced::new(model, &ambient, &mn.module, &p)?;
    let right = Balanced::new(model, &ambient, &m, &np.module)?;
    let context = Parenthesizations {
        m: &m,
        n: &n,
        p: &p,
        mn: &mn,
        np: &np,
        left: &left,
        right: &right,
    };
    let images = context.associator(model, &ambient)?;
    let first = sector_counts(&decompose(model, &ambient, target, &left.module)?);
    let second = sector_counts(&decompose(model, &ambient, target, &right.module)?);
    if first != second {
        return Err(Error::plain(
            "relative associator changed the simple-sector decomposition",
        ));
    }
    if let Some([grade, column]) = entry {
        let image = images
            .get(grade)
            .and_then(|row| row.get(column))
            .ok_or_else(|| {
                Error::plain("relative associator grade or basis column is out of range")
            })?;
        Ok(format!(
            "rvDoubleRelativeAssocImage[{grade},{column}] = target={}, phase={}, conductor={}",
            image.target, image.phase, ambient.conductor
        ))
    } else {
        let nonzero = images.iter().flatten().filter(|x| x.phase != 0).count();
        Ok(format!(
            "rvDoubleRelativeAssociator = dimension={}, conductor={}, nonzero_phases={nonzero}, sectors={first:?}, verified=true",
            left.module.dimension(),ambient.conductor,
        ))
    }
}

struct Parenthesizations<'a> {
    m: &'a Module,
    n: &'a Module,
    p: &'a Module,
    mn: &'a Balanced,
    np: &'a Balanced,
    left: &'a Balanced,
    right: &'a Balanced,
}

impl Parenthesizations<'_> {
    fn associator(&self, model: &mut Model, ambient: &Ambient<'_>) -> Result<Vec<Vec<Image>>> {
        if self.left.module.fibers != self.right.module.fibers {
            return Err(Error::plain(
                "relative associator has unequal graded dimensions",
            ));
        }
        let mut images = Vec::new();
        for (g, basis) in self.left.bases.iter().enumerate() {
            model.spend(basis.len())?;
            let mut row = Vec::new();
            let mut seen = vec![false; basis.len()];
            for pair in basis {
                let Pair {
                    grades: [mn, w],
                    columns: [u, k],
                } = *pair;
                let Pair {
                    grades: [z, v],
                    columns: [i, j],
                } = self.mn.bases[mn][u];
                let image = self.on_right(ambient, [z, v, w], [i, j, k])?;
                if std::mem::replace(&mut seen[image.target], true) {
                    return Err(Error::plain(
                        "relative associator is not invertible on a graded basis",
                    ));
                }
                row.push(image);
            }
            if self.right.bases[g].len() != row.len() {
                return Err(Error::plain(
                    "relative associator has the wrong graded basis",
                ));
            }
            images.push(row);
        }
        self.check_descent(model, ambient, &images)?;
        self.check_actions(model, ambient, &images)?;
        Ok(images)
    }

    fn on_right(
        &self,
        ambient: &Ambient<'_>,
        [z, v, w]: [usize; 3],
        [i, j, k]: [usize; 3],
    ) -> Result<Image> {
        let image = self.np.canonical(
            ambient,
            self.n,
            self.p,
            Pair {
                grades: [v, w],
                columns: [j, k],
            },
            ambient.alpha(z, v, w),
        )?;
        self.right.canonical(
            ambient,
            self.m,
            &self.np.module,
            Pair {
                grades: [z, ambient.mul(v, w)],
                columns: [i, image.target],
            },
            image.phase,
        )
    }

    fn on_left(
        &self,
        ambient: &Ambient<'_>,
        [z, v, w]: [usize; 3],
        [i, j, k]: [usize; 3],
    ) -> Result<Image> {
        let image = self.mn.canonical(
            ambient,
            self.m,
            self.n,
            Pair {
                grades: [z, v],
                columns: [i, j],
            },
            0,
        )?;
        self.left.canonical(
            ambient,
            &self.mn.module,
            self.p,
            Pair {
                grades: [ambient.mul(z, v), w],
                columns: [image.target, k],
            },
            image.phase,
        )
    }

    fn check_descent(
        &self,
        model: &mut Model,
        ambient: &Ambient<'_>,
        images: &[Vec<Image>],
    ) -> Result<()> {
        model.spend(4 * self.m.dimension() * self.n.dimension() * self.p.dimension())?;
        for (z, &a) in self.m.fibers.iter().enumerate() {
            for (v, &b) in self.n.fibers.iter().enumerate() {
                for (w, &c) in self.p.fibers.iter().enumerate() {
                    let g = ambient.mul(ambient.mul(z, v), w);
                    for i in 0..a {
                        for j in 0..b {
                            for k in 0..c {
                                let first = self.on_left(ambient, [z, v, w], [i, j, k])?;
                                let last = self.on_right(ambient, [z, v, w], [i, j, k])?;
                                let assoc = images[g][first.target];
                                if assoc.target != last.target
                                    || ambient.sum([first.phase, assoc.phase]) != last.phase
                                {
                                    return Err(Error::plain("relative associator depends on the balanced representative"));
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn check_actions(
        &self,
        model: &mut Model,
        ambient: &Ambient<'_>,
        images: &[Vec<Image>],
    ) -> Result<()> {
        for (side, &boundary) in self.left.module.ends.iter().enumerate() {
            let group = &ambient.boundaries[boundary].members;
            model.spend(group.len() * self.left.module.dimension())?;
            for &x in group {
                for (g, row) in images.iter().enumerate() {
                    for (i, iso) in row.iter().enumerate() {
                        let left = self.left.module.grade_action(ambient, side, x, g, i);
                        let right = self
                            .right
                            .module
                            .grade_action(ambient, side, x, g, iso.target);
                        let after = images[ambient.mul(x, g)][left.target];
                        if after.target != right.target
                            || ambient.sum([left.phase, after.phase])
                                != ambient.sum([iso.phase, right.phase])
                        {
                            return Err(Error::plain("relative associator does not intertwine the outer boundary actions"));
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
