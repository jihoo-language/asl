use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::{self, Coordinate};
use crate::sheaf::{self, Matrix};

use super::{linear, Disk, Hom, Origin, Pair};

pub(crate) struct Sequence {
    pub(crate) prime: Prime,
    pub(crate) nearby: usize,
    pub(crate) terms: Vec<Rc<Disk>>,
    pub(crate) total: Matrix,
    pub(crate) braid: Option<Braid>,
}

pub(crate) struct Braid {
    pub(crate) source: Rc<str>,
    word: Vec<i64>,
    permutation: Vec<usize>,
    steps: Vec<Step>,
}

struct Step {
    index: usize,
    positive: bool,
    original: Rc<Disk>,
    around: Rc<Disk>,
}

impl Sequence {
    pub(crate) fn build(model: &mut Model, terms: Vec<Rc<Disk>>) -> Result<Self> {
        let Some(first) = terms.first() else {
            return Err(Error::plain("disk sequence requires at least one quiver"));
        };
        let (prime, nearby) = (first.prime, first.nearby);
        let mut total = Matrix::identity(model, nearby)?;
        for disk in &terms {
            model.spend(1)?;
            if (disk.prime, disk.nearby) != (prime, nearby) {
                return Err(Error::plain(
                    "disk sequence requires one based nearby space and prime field",
                ));
            }
            total = Matrix::compose(model, prime, &total, &disk.t)?;
        }
        model.allocate(terms.len() + total.cells() + 2)?;
        Ok(Self {
            prime,
            nearby,
            terms,
            total,
            braid: None,
        })
    }

    pub(crate) fn hurwitz(
        model: &mut Model,
        name: Rc<str>,
        source: &Self,
        word: Vec<i64>,
    ) -> Result<Self> {
        let n = source.terms.len();
        sheaf::check_cells(n.saturating_add(5usize.saturating_mul(word.len())))?;
        model.spend(n + word.len())?;
        let mut terms = source.terms.clone();
        let mut permutation: Vec<_> = (0..n).collect();
        let mut steps = Vec::with_capacity(word.len());
        for &generator in &word {
            let i = usize::try_from(generator.unsigned_abs())
                .ok()
                .and_then(|index| index.checked_sub(1))
                .filter(|&index| index < n - 1)
                .ok_or_else(|| {
                    Error::plain(format!(
                        "Hurwitz generator {generator} is out of range for {n} strands"
                    ))
                })?;
            let positive = generator > 0;
            let (first, next) = (terms[i].clone(), terms[i + 1].clone());
            let (original, around) = if positive {
                let conjugated = transform(model, &first, &next.t, &next.t_inverse)?;
                terms[i] = next.clone();
                terms[i + 1] = Rc::new(conjugated);
                (first, next)
            } else {
                let conjugated = Rc::new(transform(model, &next, &first.t_inverse, &first.t)?);
                terms[i] = conjugated.clone();
                terms[i + 1] = first.clone();
                (conjugated, first)
            };
            permutation.swap(i, i + 1);
            steps.push(Step {
                index: i,
                positive,
                original,
                around,
            });
        }
        let mut result = Self::build(model, terms)?;
        if !same_matrix(model, &source.total, &result.total)? {
            return Err(Error::plain(
                "Hurwitz word did not preserve total monodromy",
            ));
        }
        model.allocate(n + 5 * steps.len() + 1)?;
        result.braid = Some(Braid {
            source: name,
            word,
            permutation,
            steps,
        });
        Ok(result)
    }

    pub(crate) fn point(&self, i: usize) -> Result<Rc<Disk>> {
        self.terms
            .get(i)
            .cloned()
            .ok_or_else(|| Error::plain("disk strand index out of range"))
    }

    pub(crate) fn product_report(&self, model: &mut Model) -> Result<String> {
        linear::report(model, &self.total, self.nearby, self.nearby)
    }

    pub(crate) fn same_data(&self, model: &mut Model, other: &Self, quivers: bool) -> Result<bool> {
        if (self.prime, self.nearby, self.terms.len())
            != (other.prime, other.nearby, other.terms.len())
        {
            return Ok(false);
        }
        for (a, b) in self.terms.iter().zip(&other.terms) {
            let equal = if quivers {
                same_disk(model, a, b)?
            } else {
                same_matrix(model, &a.t, &b.t)?
            };
            if !equal {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

fn transform(
    model: &mut Model,
    disk: &Disk,
    right: &Matrix,
    right_inverse: &Matrix,
) -> Result<Disk> {
    let p = disk.prime;
    let u = Matrix::compose(model, p, &disk.u, right)?;
    let v = Matrix::compose(model, p, right_inverse, &disk.v)?;
    let output = Disk::new(
        model,
        p,
        disk.nearby,
        disk.vanishing,
        u,
        v,
        Origin::Declared,
    )?;
    if !same_matrix(model, &disk.s, &output.s)? {
        return Err(Error::plain(
            "Hurwitz move did not preserve vanishing monodromy",
        ));
    }
    Ok(output)
}

impl Braid {
    pub(crate) fn action(&self, model: &mut Model, rows: &mut [Vec<Coordinate>]) -> Result<()> {
        model.spend(rows.len())?;
        let mut entries: usize = rows.iter().map(Vec::len).sum();
        for step in &self.steps {
            let p = step.original.prime;
            let (before, after) = rows.split_at_mut(step.index + 1);
            let (a, b) = (&mut before[step.index], &mut after[0]);
            entries -= a.len() + b.len();
            if step.positive {
                let correction = step.original.t.apply(model, p, b)?;
                let correction = prime_linear::add_scaled(model, p, &correction, b, p.negate(1))?;
                let sum = prime_linear::add_scaled(model, p, a, &correction, 1)?;
                let image = step.around.t_inverse.apply(model, p, &sum)?;
                *a = std::mem::replace(b, image);
            } else {
                let image = step.around.t.apply(model, p, b)?;
                let correction = step.original.t.apply(model, p, a)?;
                let difference = prime_linear::add_scaled(model, p, a, &correction, p.negate(1))?;
                let image = prime_linear::add_scaled(model, p, &image, &difference, 1)?;
                *b = std::mem::replace(a, image);
            }
            entries += a.len() + b.len();
            sheaf::check_cells(2 * entries + rows.len())?;
        }
        Ok(())
    }
}

pub(crate) struct SequenceMap {
    pub(crate) source: Rc<str>,
    pub(crate) target: Rc<str>,
    pub(crate) a: Rc<Matrix>,
    vanishing: Vec<Rc<Matrix>>,
}

impl SequenceMap {
    pub(crate) fn build(
        model: &mut Model,
        names: (Rc<str>, Rc<str>),
        source: &Sequence,
        target: &Sequence,
        maps: Vec<(Rc<Hom>, Pair)>,
    ) -> Result<Self> {
        let (source_name, target_name) = names;
        let n = source.terms.len();
        if source.prime != target.prime || target.terms.len() != n || maps.len() != n {
            return Err(Error::plain(
                "disk sequence map needs a morphism per corresponding strand over the same field",
            ));
        }
        let a = maps[0].1.a.clone();
        let mut vanishing = Vec::with_capacity(n);
        for (i, ((left, right), (hom, map))) in
            source.terms.iter().zip(&target.terms).zip(maps).enumerate()
        {
            if !same_disk(model, left, &hom.source)? || !same_disk(model, right, &hom.target)? {
                return Err(Error::plain(format!(
                    "disk sequence morphism endpoints disagree at strand {i}"
                )));
            }
            if !same_matrix(model, &a, &map.a)? {
                return Err(Error::plain(format!(
                    "disk sequence morphisms differ on the nearby map at strand {i}"
                )));
            }
            vanishing.push(map.b);
        }
        let cells = 3 + a.cells() + vanishing.iter().map(|b| b.cells()).sum::<usize>();
        model.allocate(cells)?;
        Ok(Self {
            source: source_name,
            target: target_name,
            a,
            vanishing,
        })
    }

    pub(crate) fn braid(
        &self,
        model: &mut Model,
        names: (Rc<str>, Rc<str>),
        left: &Sequence,
        right: &Sequence,
    ) -> Result<Self> {
        let (Some(a), Some(b)) = (&left.braid, &right.braid) else {
            return Err(Error::plain(
                "braided disk map needs declared Hurwitz images at both endpoints",
            ));
        };
        model.spend(a.word.len() + 1)?;
        if a.source != self.source || b.source != self.target || a.word != b.word {
            return Err(Error::plain(
                "braided disk map needs the same signed word on the original endpoints",
            ));
        }
        let vanishing: Vec<_> = a
            .permutation
            .iter()
            .map(|&j| self.vanishing[j].clone())
            .collect();
        for ((source, target), matrix) in left.terms.iter().zip(&right.terms).zip(&vanishing) {
            let p = source.prime;
            let u1 = Matrix::compose(model, p, &target.u, &self.a)?;
            let u2 = Matrix::compose(model, p, matrix, &source.u)?;
            let v1 = Matrix::compose(model, p, &target.v, matrix)?;
            let v2 = Matrix::compose(model, p, &self.a, &source.v)?;
            if !same_matrix(model, &u1, &u2)? || !same_matrix(model, &v1, &v2)? {
                return Err(Error::plain(
                    "transported disk sequence map does not intertwine both arrows",
                ));
            }
        }
        model.allocate(2 + vanishing.len())?;
        Ok(Self {
            source: names.0,
            target: names.1,
            a: self.a.clone(),
            vanishing,
        })
    }

    pub(crate) fn report(
        &self,
        model: &mut Model,
        from: &Sequence,
        to: &Sequence,
        i: usize,
    ) -> Result<String> {
        let source = from.point(i)?;
        let target = to.point(i)?;
        let a = linear::report(model, &self.a, to.nearby, from.nearby)?;
        let b = linear::report(
            model,
            &self.vanishing[i],
            target.vanishing,
            source.vanishing,
        )?;
        Ok(format!("{a} / {b}"))
    }
}

fn same_disk(model: &mut Model, a: &Disk, b: &Disk) -> Result<bool> {
    Ok(
        (a.prime, a.nearby, a.vanishing) == (b.prime, b.nearby, b.vanishing)
            && same_matrix(model, &a.u, &b.u)?
            && same_matrix(model, &a.v, &b.v)?,
    )
}

fn same_matrix(model: &mut Model, a: &Matrix, b: &Matrix) -> Result<bool> {
    model.spend(a.cells().min(b.cells()))?;
    Ok(a == b)
}
