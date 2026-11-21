use crate::bp::Shape;
use crate::bp_ext::{self, Domain, Group};
use crate::bp_polynomial::Monomial;
use crate::bp_smith::{self, IntegralKernel, Matrix};
use crate::error::{Error, Result};
use crate::model::Model;

#[derive(Clone, Copy)]
pub(crate) struct Grade {
    pub(crate) top: Shape,
    pub(crate) r: usize,
    pub(crate) h: usize,
    pub(crate) t: usize,
}

impl Grade {
    pub(crate) fn validate(self) -> Result<Self> {
        if self.h > self.top.height {
            return Err(Error::plain(
                "BP height filtration exceeds the ambient height",
            ));
        }
        Ok(self)
    }

    fn shape(self) -> Shape {
        Shape {
            height: self.h,
            ..self.top
        }
    }

    fn leading(self) -> Domain {
        debug_assert!(self.h > 0);
        Domain::relative(self.shape(), self.h - 1)
    }

    fn has_positive_slice(self) -> bool {
        self.t != 0 && self.h != 0 && self.leading().has_monomials(self.t)
    }

    fn next(self) -> Option<Self> {
        let h = self.h.checked_sub(self.r)?;
        let top = Shape {
            arrows: self.top.arrows + 1,
            ..self.top
        };
        Some(Self { top, h, ..self })
    }
}

pub(crate) enum Image {
    Free(usize),
    Finite(i128),
}

pub(crate) fn page(model: &mut Model, grade: Grade) -> Result<Group> {
    if !grade.has_positive_slice() {
        let free = usize::from(grade.t == 0 && grade.h == 0 && grade.top.arrows == 0);
        return Ok(Group {
            free,
            torsion: Vec::new(),
        });
    }
    if grade.r == 0 {
        let free = grade.leading().basis(model, grade.t)?.len();
        return Ok(Group {
            free,
            torsion: Vec::new(),
        });
    }
    if grade.r == 1 {
        return grade.leading().group(model, grade.t);
    }
    let prepared = PageLattices::build(model, grade)?;
    let (_, factors) = prepared.reduce(model)?;
    Ok(Group {
        free: 0,
        torsion: factors.into_iter().filter(|&value| value != 1).collect(),
    })
}

pub(crate) fn image(model: &mut Model, grade: Grade) -> Result<Image> {
    let Some(target) = grade.next() else {
        return Ok(Image::Finite(1));
    };
    if !grade.has_positive_slice() || !target.has_positive_slice() {
        return Ok(Image::Finite(1));
    }
    if grade.r == 0 {
        let (_, rank) = grade.leading().incoming_basis(model, grade.t)?;
        if rank == 0 {
            return Ok(Image::Finite(1));
        }
        return Ok(Image::Free(rank));
    }
    let prepared = PageLattices::build(model, target)?;
    let (modulus, original) = prepared.reduce(model)?;
    if original.iter().all(|&value| value == 1) {
        return Ok(Image::Finite(1));
    }
    let Some(extra) = MappedCycles::boundary(
        model,
        grade.shape(),
        target.h,
        prepared.ambient_height,
        &prepared.ambient,
        grade.t,
    )?
    else {
        return Ok(Image::Finite(1));
    };
    let expanded = prepared.quotient_matrix(model, modulus, Some(&extra))?;
    let expanded =
        bp_smith::factors_at_precision(model, &expanded, grade.top.prime, modulus, prepared.rank)?
            .ok_or_else(|| Error::plain("BP height image reduction lost a certified pivot"))?;
    let prime = i128::from(grade.top.prime);
    model.spend(original.len() + expanded.len())?;
    let index =
        |factors: &[i128]| -> u32 { factors.iter().map(|&factor| factor.ilog(prime)).sum() };
    let exponent = index(&original)
        .checked_sub(index(&expanded))
        .ok_or_else(|| {
            Error::plain("BP height image indices disagree with the denominator inclusion")
        })?;
    let order = prime
        .checked_pow(exponent)
        .ok_or_else(|| Error::plain("BP height image order exceeds the signed i128 range"))?;
    Ok(Image::Finite(order))
}

struct PageLattices {
    prime: u32,
    rank: usize,
    ambient_height: usize,
    ambient: Vec<Monomial>,
    shadow: Option<MappedCycles>,
    boundaries: Option<MappedCycles>,
}

impl PageLattices {
    fn build(model: &mut Model, grade: Grade) -> Result<Self> {
        let shape = grade.shape();
        let cutoff = grade.h.saturating_sub(grade.r);
        let ambient = bp_ext::basis(model, shape, grade.t)?;
        let rank = cycle_rank(
            ambient.len(),
            constraint_rank(model, shape, cutoff, grade.t)?,
        )?;
        let prime = shape.prime;
        if rank == 0 {
            return Ok(Self {
                prime,
                rank,
                ambient_height: grade.h,
                ambient,
                shadow: None,
                boundaries: None,
            });
        }
        let lower_shape = Shape {
            height: grade.h - 1,
            ..shape
        };
        let shadow =
            MappedCycles::inclusion(model, lower_shape, cutoff, grade.h, &ambient, grade.t)?;
        let boundaries = if let Some(arrows) = shape.arrows.checked_sub(1) {
            // The source height is h + r - 1 before clamping to N.
            // Its cycle condition is d(x) in F_h, even after clamping.
            let height = grade.h.saturating_add(grade.r - 1).min(grade.top.height);
            let source = Shape {
                height,
                arrows,
                ..shape
            };
            MappedCycles::boundary(model, source, grade.h, grade.h, &ambient, grade.t)?
        } else {
            None
        };
        Ok(Self {
            prime,
            rank,
            ambient_height: grade.h,
            ambient,
            shadow,
            boundaries,
        })
    }

    fn reduce(&self, model: &mut Model) -> Result<(u64, Vec<i128>)> {
        let mut modulus = u64::from(self.prime) * u64::from(self.prime);
        if self.rank == 0 {
            return Ok((modulus, Vec::new()));
        }
        loop {
            let matrix = self.quotient_matrix(model, modulus, None)?;
            // Positive-weight E_1 is finite p-primary. Every later page
            // is a subquotient, so the denominator has the numerator's
            // rational rank. Its ambient invariant factors are valid
            // because the numerator is saturated in the ambient group.
            if let Some(factors) =
                bp_smith::factors_at_precision(model, &matrix, self.prime, modulus, self.rank)?
            {
                return Ok((modulus, factors));
            }
            modulus = bp_smith::next_precision(modulus, self.prime).ok_or_else(|| {
                Error::plain("BP height Smith precision exceeds the u64 modulus limit")
            })?;
        }
    }

    fn quotient_matrix(
        &self,
        model: &mut Model,
        modulus: u64,
        extra: Option<&MappedCycles>,
    ) -> Result<Matrix> {
        let maps = self
            .shadow
            .iter()
            .chain(self.boundaries.iter())
            .chain(extra);
        let columns = maps.clone().map(|map| map.cycles.dimension).sum();
        let mut matrix = Matrix::new(self.ambient.len(), columns)?;
        let mut first = 0;
        for map in maps {
            map.record(model, &mut matrix, first, modulus)?;
            first += map.cycles.dimension;
        }
        Ok(matrix)
    }
}

struct CycleLattice {
    source: Vec<Monomial>,
    condition: Matrix,
    prime: u32,
    rank: usize,
    exponent: u64,
    dimension: usize,
}

impl CycleLattice {
    fn build(model: &mut Model, shape: Shape, cutoff: usize, t: usize) -> Result<Self> {
        let source = bp_ext::basis(model, shape, t)?;
        let prime = shape.prime;
        let rank = if source.is_empty() {
            0
        } else {
            constraint_rank(model, shape, cutoff, t)?
        };
        let dimension = cycle_rank(source.len(), rank)?;
        let (condition, exponent) = if rank == 0 || dimension == 0 {
            (Matrix::new(0, source.len())?, 1)
        } else {
            let next = shape.next()?;
            let high = Domain::relative(next, cutoff).basis(model, t)?;
            let condition = Domain::relative(shape, cutoff).matrix(model, &source, &high)?;
            let exponent = bp_smith::primary_factors(model, &condition, prime, rank)?
                .last()
                .copied()
                .unwrap_or(1) as u64;
            (condition, exponent)
        };
        Ok(Self {
            source,
            condition,
            prime,
            rank,
            exponent,
            dimension,
        })
    }

    fn residues(&self, model: &mut Model, modulus: u64) -> Result<IntegralKernel> {
        IntegralKernel::build(
            model,
            &self.condition,
            self.prime,
            self.rank,
            self.exponent,
            modulus,
        )
    }
}

fn constraint_rank(model: &mut Model, shape: Shape, cutoff: usize, t: usize) -> Result<usize> {
    if cutoff >= shape.height {
        return Ok(0);
    }
    let next = Shape {
        arrows: shape.arrows + 1,
        ..shape
    };
    if !Domain::relative(next, cutoff).has_monomials(t) {
        return Ok(0);
    }
    // The filtered constraint is the relative differential. Heights
    // below zero and height zero agree in positive internal degree.
    let (_, rank) = Domain::relative(shape, cutoff).incoming_basis(model, t)?;
    Ok(rank)
}

fn cycle_rank(dimension: usize, constraint_rank: usize) -> Result<usize> {
    dimension
        .checked_sub(constraint_rank)
        .ok_or_else(|| Error::plain("BP height constraint rank exceeds the cochain dimension"))
}

struct MappedCycles {
    cycles: CycleLattice,
    map: CycleMap,
}

enum CycleMap {
    Inclusion(Vec<usize>),
    Differential(Matrix),
}

impl MappedCycles {
    fn inclusion(
        model: &mut Model,
        shape: Shape,
        cutoff: usize,
        ambient_height: usize,
        ambient: &[Monomial],
        t: usize,
    ) -> Result<Option<Self>> {
        let cycles = CycleLattice::build(model, shape, cutoff, t)?;
        if cycles.dimension == 0 {
            return Ok(None);
        }
        let mut rows = Vec::with_capacity(cycles.source.len());
        for monomial in &cycles.source {
            let embedded = widen(model, monomial, shape.height, ambient_height)?;
            rows.push(bp_ext::index_of(model, ambient, &embedded)?);
        }
        Ok(Some(Self {
            cycles,
            map: CycleMap::Inclusion(rows),
        }))
    }

    fn boundary(
        model: &mut Model,
        source: Shape,
        cutoff: usize,
        ambient_height: usize,
        ambient: &[Monomial],
        t: usize,
    ) -> Result<Option<Self>> {
        let cycles = CycleLattice::build(model, source, cutoff, t)?;
        if cycles.dimension == 0 {
            return Ok(None);
        }
        let embedded = ambient
            .iter()
            .map(|monomial| widen(model, monomial, ambient_height, source.height))
            .collect::<Result<Vec<_>>>()?;
        let differential =
            bp_ext::low_matrix(model, source, &cycles.source, &embedded, ambient_height)?;
        Ok(Some(Self {
            cycles,
            map: CycleMap::Differential(differential),
        }))
    }

    fn record(
        &self,
        model: &mut Model,
        matrix: &mut Matrix,
        first: usize,
        modulus: u64,
    ) -> Result<()> {
        if self.cycles.dimension == 0 {
            return Ok(());
        }
        let kernel = self.cycles.residues(model, modulus)?;
        for i in 0..self.cycles.dimension {
            let vector = kernel.generator(model, i)?;
            let column = first + i;
            match &self.map {
                CycleMap::Inclusion(rows) => {
                    model.spend(rows.len())?;
                    for (&row, coefficient) in rows.iter().zip(vector) {
                        matrix.set(row, column, i128::from(coefficient));
                    }
                }
                CycleMap::Differential(differential) => {
                    let image = differential.apply_residue(model, &vector, modulus)?;
                    for (row, coefficient) in image.into_iter().enumerate() {
                        matrix.set(row, column, i128::from(coefficient));
                    }
                }
            }
        }
        Ok(())
    }
}

fn widen(model: &mut Model, monomial: &Monomial, old: usize, new: usize) -> Result<Monomial> {
    model.spend(1 + monomial.len())?;
    Ok(monomial
        .iter()
        .map(|&(variable, power)| {
            let level = variable % old;
            let slot = variable / old;
            (slot * new + level, power)
        })
        .collect())
}
