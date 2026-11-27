use crate::error::{Error, Result};
use crate::rv_lens::Rational;

use super::super::super::matrix::{Field, Matrix};
use super::tensor::Layout;
use super::{Ambient, Image, Module, Pair};

#[derive(Clone)]
pub(super) struct GradedMap {
    pub(super) grades: Vec<Matrix>,
}

impl GradedMap {
    pub(super) fn extend(
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

pub(super) fn tensor_at(
    field: &mut Field<'_>,
    ambient: &Ambient<'_>,
    from: Layout<'_>,
    to: Layout<'_>,
    inclusions: [&GradedMap; 2],
    grade: usize,
) -> Result<Matrix> {
    let source = &from.product.bases[grade];
    let mut result = Matrix::zero(to.product.bases[grade].len(), source.len())?;
    for (col, &pair) in source.iter().enumerate() {
        tensor_terms(field, ambient, to, inclusions, pair, |field, row, term| {
            result.set(row, col, field.add(result.get(row, col), &term)?);
            Ok(())
        })?;
    }
    // Independently compare images of noncanonical pairs with their quotient
    // representatives. This is also used for inclusions of the fused tube loop.
    for u in 0..ambient.order() {
        let v = ambient.mul(ambient.inv(u), grade);
        for i in 0..from.factors[0].fibers[u] {
            for j in 0..from.factors[1].fibers[v] {
                let pair = Pair {
                    grades: [u, v],
                    columns: [i, j],
                };
                let image = from.reduce(ambient, pair, 0)?;
                if pair == source[image.target] && image.phase == 0 {
                    continue;
                }
                field.model.spend(result.rows)?;
                let mut column = vec![Rational::zero(); result.rows];
                tensor_terms(field, ambient, to, inclusions, pair, |field, row, term| {
                    column[row] = field.add(&column[row], &term)?;
                    Ok(())
                })?;
                for (row, given) in column.iter().enumerate() {
                    if *given != field.shift(result.get(row, image.target), image.phase)? {
                        return Err(Error::plain(
                            "graded tensor map fails to descend through balancing",
                        ));
                    }
                }
            }
        }
    }
    Ok(result)
}

fn tensor_terms(
    field: &mut Field<'_>,
    ambient: &Ambient<'_>,
    to: Layout<'_>,
    inclusions: [&GradedMap; 2],
    pair: Pair,
    mut term: impl FnMut(&mut Field<'_>, usize, Rational) -> Result<()>,
) -> Result<()> {
    let [u, v] = pair.grades;
    let [i, j] = pair.columns;
    let a = &inclusions[0].grades[u];
    let b = &inclusions[1].grades[v];
    field.model.spend(a.rows * b.rows)?;
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
            let scale = field.multiply(p, q)?;
            let image = to.reduce(
                ambient,
                Pair {
                    grades: [u, v],
                    columns: [s, t],
                },
                0,
            )?;
            let coefficient = field.shift(&scale, image.phase)?;
            term(field, image.target, coefficient)?;
        }
    }
    Ok(())
}

impl GradedMap {
    pub(super) fn identity(module: &Module) -> Result<Self> {
        Ok(Self {
            grades: module
                .fibers
                .iter()
                .map(|&n| Matrix::identity(n))
                .collect::<Result<Vec<_>>>()?,
        })
    }

    pub(super) fn zero(source: &Module, target: &Module) -> Result<Self> {
        Ok(Self {
            grades: source
                .fibers
                .iter()
                .zip(&target.fibers)
                .map(|(&s, &t)| Matrix::zero(t, s))
                .collect::<Result<Vec<_>>>()?,
        })
    }

    pub(super) fn compose(&self, field: &mut Field<'_>, before: &Self) -> Result<Self> {
        Ok(Self {
            grades: self
                .grades
                .iter()
                .zip(&before.grades)
                .map(|(a, b)| a.multiply(field, b))
                .collect::<Result<Vec<_>>>()?,
        })
    }

    pub(super) fn adjoint(&self, field: &mut Field<'_>) -> Result<Self> {
        Ok(Self {
            grades: self
                .grades
                .iter()
                .map(|a| a.adjoint(field))
                .collect::<Result<Vec<_>>>()?,
        })
    }

    pub(super) fn add_assign(&mut self, field: &mut Field<'_>, right: &Self) -> Result<()> {
        for (a, b) in self.grades.iter_mut().zip(&right.grades) {
            a.add_assign(field, b)?;
        }
        Ok(())
    }

    pub(super) fn scale(&mut self, field: &mut Field<'_>, scalar: &Rational) -> Result<()> {
        for a in &mut self.grades {
            a.scale(field, scalar)?;
        }
        Ok(())
    }

    pub(super) fn inner(&self, field: &mut Field<'_>, other: &Self) -> Result<Rational> {
        let mut sum = Rational::zero();
        for (a, b) in self.grades.iter().zip(&other.grades) {
            let term = a.inner(field, b)?;
            sum = field.add(&sum, &term)?;
        }
        Ok(sum)
    }

    pub(super) fn after_monomial(
        &self,
        field: &mut Field<'_>,
        images: &[Vec<Image>],
    ) -> Result<Self> {
        Ok(Self {
            grades: self
                .grades
                .iter()
                .zip(images)
                .map(|(mat, row)| mat.transformed(field, |i| (row[i].target, row[i].phase)))
                .collect::<Result<Vec<_>>>()?,
        })
    }

    pub(super) fn monomial(field: &mut Field<'_>, images: &[Vec<Image>]) -> Result<Self> {
        let mut grades = Vec::new();
        for row in images {
            let mut mat = Matrix::zero(row.len(), row.len())?;
            for (i, image) in row.iter().enumerate() {
                mat.set(image.target, i, field.phase(image.phase)?);
            }
            grades.push(mat);
        }
        Ok(Self { grades })
    }
}

pub(super) fn tensor(
    field: &mut Field<'_>,
    ambient: &Ambient<'_>,
    from: Layout<'_>,
    to: Layout<'_>,
    maps: [&GradedMap; 2],
) -> Result<GradedMap> {
    Ok(GradedMap {
        grades: (0..ambient.order())
            .map(|g| tensor_at(field, ambient, from, to, maps, g))
            .collect::<Result<Vec<_>>>()?,
    })
}

pub(super) struct ChannelMap {
    pub(super) inclusion: GradedMap,
    pub(super) projection: GradedMap,
    pub(super) inverse_norm: Rational,
}

pub(super) fn hom_basis(
    field: &mut Field<'_>,
    ambient: &Ambient<'_>,
    source: &Module,
    target: &Module,
    representative: usize,
    hom: &super::Hom,
) -> Result<Vec<ChannelMap>> {
    hom.orbits
        .iter()
        .map(|orbit| {
            let mut base = Matrix::zero(hom.rows, hom.cols)?;
            for &(p, phase) in orbit {
                base.set(p / hom.cols, p % hom.cols, field.phase(phase)?);
            }
            let inclusion =
                GradedMap::extend(field, ambient, source, target, representative, &base)?;
            let norm = inclusion.inner(field, &inclusion)?;
            let inverse_norm = field.inverse(&norm)?;
            let scale = field.multiply(
                &inverse_norm,
                &Rational::integer(source.dimension() as i128),
            )?;
            let mut projection = inclusion.adjoint(field)?;
            projection.scale(field, &scale)?;
            let product = projection.compose(field, &inclusion)?;
            if product.grades != GradedMap::identity(source)?.grades {
                return Err(Error::plain(
                    "annular trivalent projection does not split its inclusion",
                ));
            }
            Ok(ChannelMap {
                inclusion,
                projection,
                inverse_norm,
            })
        })
        .collect()
}
