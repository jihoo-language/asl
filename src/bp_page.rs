use std::fmt;

use crate::bp::{self, Shape};
use crate::bp_ext;
use crate::bp_lattice::{self, Lattice};
use crate::bp_polynomial::{Monomial, Polynomial};
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::{Basis, Coordinate};

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum Mode {
    Bockstein,
    Novikov,
}

impl Mode {
    pub(crate) fn dimension_name(self) -> &'static str {
        match self {
            Self::Bockstein => "bpBockdim",
            Self::Novikov => "bpNovdim",
        }
    }

    pub(crate) fn type_name(self) -> &'static str {
        match self {
            Self::Bockstein => "BPBock",
            Self::Novikov => "BPNov",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) mode: Mode,
    pub(crate) shape: Shape,
    pub(crate) page: usize,
    pub(crate) t: usize,
    pub(crate) filtration: usize,
}

impl fmt::Display for Grade {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self {
            mode,
            shape,
            page,
            t,
            filtration,
        } = *self;
        let Shape {
            prime,
            height,
            arrows,
            ..
        } = shape;
        write!(
            out,
            "{}({prime}, {height}, {page}, {arrows}, {t}, {filtration})",
            mode.type_name()
        )
    }
}

impl Grade {
    pub(crate) fn cup(self, other: Self) -> Result<Self> {
        if self.mode != other.mode || self.page != other.page {
            return Err(Error::plain(
                "bpCup requires the same BP page family and page index",
            ));
        }
        Ok(Self {
            shape: self.shape.cup(other.shape)?,
            t: self
                .t
                .checked_add(other.t)
                .ok_or_else(|| Error::plain("BP cup product internal degree is out of range"))?,
            filtration: self
                .filtration
                .checked_add(other.filtration)
                .ok_or_else(|| Error::plain("BP cup product filtration is out of range"))?,
            ..self
        })
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Representative,
    Differential,
    Next,
    Previous,
}

impl Operation {
    pub(crate) fn from_name(name: &str) -> Option<(Mode, Self)> {
        let (mode, tail) = if let Some(tail) = name.strip_prefix("bpBock") {
            (Mode::Bockstein, tail)
        } else {
            (Mode::Novikov, name.strip_prefix("bpNov")?)
        };
        let operation = match tail {
            "Rep" => Self::Representative,
            "D" => Self::Differential,
            "Next" => Self::Next,
            "Prev" => Self::Previous,
            _ => return None,
        };
        Some((mode, operation))
    }

    pub(crate) fn target(self, grade: Grade) -> Result<Grade> {
        Ok(match self {
            Self::Representative => grade,
            Self::Differential => Grade {
                shape: grade.shape.next()?,
                filtration: grade
                    .filtration
                    .checked_add(grade.page)
                    .ok_or_else(|| Error::plain("BP differential filtration is out of range"))?,
                ..grade
            },
            Self::Next => Grade {
                page: grade
                    .page
                    .checked_add(1)
                    .ok_or_else(|| Error::plain("BP page index is out of range"))?,
                ..grade
            },
            Self::Previous => Grade {
                page: grade
                    .page
                    .checked_sub(1)
                    .ok_or_else(|| Error::plain("BP page zero has no previous page"))?,
                ..grade
            },
        })
    }
}

pub(crate) struct Page {
    grade: Grade,
    prime: Prime,
    monomials: Vec<Monomial>,
    cycles: Lattice,
    denominator: Basis,
    generators: Vec<usize>,
    common: i128,
}

impl Page {
    pub(crate) fn build(model: &mut Model, grade: Grade) -> Result<Self> {
        let prime = Prime::validated(grade.shape.prime);
        let monomials = bp_ext::basis(model, grade.shape, grade.t)?;
        let (cycles, denominator, common) = if grade.page == 0 || monomials.is_empty() {
            let wt = weights(model, grade.mode, grade.shape.height, &monomials)?;
            let f = grade.filtration as i128;
            let common = wt.iter().max().map_or(0, |&w| (f - w).max(0));
            let cycles = Lattice::build(model, prime, &wt, f, common, None)?;
            let mut denominator = Basis::default();
            for (i, &w) in wt.iter().enumerate() {
                if w > f {
                    denominator.insert(model, prime, vec![(i, 1)])?;
                }
            }
            (cycles, denominator, common)
        } else {
            lattice_quotient(model, prime, grade, &monomials)?
        };
        model.spend(monomials.len())?;
        let generators: Vec<_> = (0..monomials.len())
            .filter(|&i| !denominator.has_pivot(i))
            .collect();
        model.allocate(
            cycles.cells()
                + denominator.cells()
                + generators.len()
                + monomials
                    .iter()
                    .map(|monomial| 2 + 2 * monomial.len())
                    .sum::<usize>(),
        )?;
        Ok(Self {
            grade,
            prime,
            monomials,
            cycles,
            denominator,
            generators,
            common,
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.generators.len()
    }

    pub(crate) fn class(&self, model: &mut Model, input: &Polynomial) -> Result<Vec<Coordinate>> {
        let coordinates = bp_ext::coordinates(model, input, &self.monomials)?;
        self.project(model, coordinates, 0)
    }

    fn reduce(&self, model: &mut Model, coordinates: Vec<Coordinate>) -> Result<Vec<Coordinate>> {
        let reduced = self.denominator.reduce(model, self.prime, coordinates)?;
        model.spend(reduced.len() * (1 + self.generators.len().max(1).ilog2() as usize))?;
        Ok(reduced
            .into_iter()
            .map(|(i, value)| (self.generators.binary_search(&i).unwrap(), value))
            .collect())
    }

    fn project(
        &self,
        model: &mut Model,
        mut coordinates: Vec<i128>,
        common: i128,
    ) -> Result<Vec<Coordinate>> {
        model.spend(coordinates.len())?;
        if coordinates.iter().all(|&value| value == 0) {
            return Ok(Vec::new());
        }
        let wrong_filtration =
            || Error::plain("BP page class requires a filtered cycle on the requested page");
        let difference = common - self.common;
        if difference > 0 {
            scale(model, self.prime, &mut coordinates, difference)?;
        } else if difference < 0 {
            let divisor =
                bp_lattice::power(self.prime, -difference).map_err(|_| wrong_filtration())?;
            for value in &mut coordinates {
                if *value % divisor != 0 {
                    return Err(wrong_filtration());
                }
                *value /= divisor;
            }
        }
        let cycles = self
            .cycles
            .encode(model, coordinates)?
            .ok_or_else(wrong_filtration)?;
        self.reduce(model, residues(self.prime, cycles))
    }

    fn lift_coordinates(&self, model: &mut Model, class: &[Coordinate]) -> Result<Vec<i128>> {
        model.spend(self.monomials.len() + class.len())?;
        let mut coordinates = vec![0; self.monomials.len()];
        for &(index, value) in class {
            coordinates[self.generators[index]] = i128::from(value);
        }
        self.cycles.decode(model, coordinates)
    }

    pub(crate) fn representative(
        &self,
        model: &mut Model,
        class: &[Coordinate],
        absolute: bool,
    ) -> Result<Polynomial> {
        if class.is_empty() {
            return Ok(Polynomial::default());
        }
        let mut coordinates = self.lift_coordinates(model, class)?;
        if absolute && self.common != 0 {
            scale(model, self.prime, &mut coordinates, self.common)?;
        }
        Polynomial::from_terms(model, self.monomials.iter().zip(coordinates))
    }

    pub(crate) fn differential(
        &self,
        model: &mut Model,
        target: &Self,
        class: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let representative = self.representative(model, class, false)?;
        let image = bp::differential(model, self.grade.shape, &representative)?;
        let coordinates = bp_ext::coordinates(model, &image, &target.monomials)?;
        target.project(model, coordinates, self.common)
    }

    pub(crate) fn cup(
        &self,
        model: &mut Model,
        other: &Self,
        target: &Self,
        left: &[Coordinate],
        right: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let f = self.representative(model, left, false)?;
        let g = other.representative(model, right, false)?;
        let product = bp::cup(model, self.grade.shape, other.grade.shape, &f, &g)?;
        let coordinates = bp_ext::coordinates(model, &product, &target.monomials)?;
        target.project(model, coordinates, self.common + other.common)
    }

    pub(crate) fn previous(
        &self,
        model: &mut Model,
        target: &Self,
        class: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let coordinates = self.lift_coordinates(model, class)?;
        target.project(model, coordinates, self.common)
    }

    pub(crate) fn next(
        &self,
        model: &mut Model,
        target: &Self,
        class: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        let mut lifts = Basis::default();
        for i in 0..target.monomials.len() {
            let coordinates = target.cycles.generator(model, i)?;
            let column = self.project(model, coordinates, target.common)?;
            lifts.insert_with_witness(model, self.prime, column, vec![(i, 1)])?;
        }
        let (remainder, witness) = lifts.reduce_with_witness(model, self.prime, class.to_vec())?;
        if !remainder.is_empty() {
            return Err(Error::plain(
                "BP next page requires a cycle for the page differential",
            ));
        }
        target.reduce(model, witness)
    }
}

fn scale(model: &mut Model, prime: Prime, coordinates: &mut [i128], exponent: i128) -> Result<()> {
    let factor = bp_lattice::power(prime, exponent)?;
    model.spend(coordinates.len())?;
    for value in coordinates {
        *value = value
            .checked_mul(factor)
            .ok_or_else(|| Error::plain("BP page representative exceeds the signed i128 range"))?;
    }
    Ok(())
}

fn weights(
    model: &mut Model,
    mode: Mode,
    height: usize,
    monomials: &[Monomial],
) -> Result<Vec<i128>> {
    model.spend(monomials.len())?;
    if mode == Mode::Bockstein {
        return Ok(vec![0; monomials.len()]);
    }
    monomials
        .iter()
        .map(|monomial| {
            model.spend(monomial.len())?;
            Ok(monomial
                .iter()
                .take_while(|&&(i, _)| i < height)
                .map(|&(_, e)| e as i128)
                .sum())
        })
        .collect()
}

pub(crate) fn rank(model: &mut Model, grade: Grade) -> Result<usize> {
    let monomials = bp_ext::basis(model, grade.shape, grade.t)?;
    if monomials.is_empty() {
        return Ok(0);
    }
    if grade.page == 0 {
        let weights = weights(model, grade.mode, grade.shape.height, &monomials)?;
        return Ok(weights
            .into_iter()
            .filter(|&weight| weight <= grade.filtration as i128)
            .count());
    }
    let prime = Prime::validated(grade.shape.prime);
    let (_, denominator, _) = lattice_quotient(model, prime, grade, &monomials)?;
    Ok(monomials.len() - denominator.rank())
}

fn lattice_quotient(
    model: &mut Model,
    prime: Prime,
    grade: Grade,
    monomials: &[Monomial],
) -> Result<(Lattice, Basis, i128)> {
    let Grade {
        mode,
        shape,
        page,
        t,
        filtration,
    } = grade;
    let f = filtration as i128;
    let r = page as i128;
    let current = weights(model, mode, shape.height, monomials)?;
    // The lower filtration index may be negative. It still participates
    // in the sum of filtration and page when computing its cycle condition.
    let lower_filtration = f - r + 1;
    let previous_shape = (shape.arrows != 0).then(|| Shape {
        arrows: shape.arrows - 1,
        ..shape
    });
    let previous = if let Some(previous_shape) = previous_shape {
        bp_ext::basis(model, previous_shape, t)?
    } else {
        Vec::new()
    };
    let lower_weights = weights(model, mode, shape.height, &previous)?;
    let min_scale = |f: i128, wt: &[i128]| wt.iter().max().map(|&w| (f - w).max(0));
    let common = min_scale(f, &current)
        .into_iter()
        .chain(min_scale(lower_filtration, &lower_weights))
        .min()
        .unwrap_or(0);

    // The homogeneous basis is empty in sufficiently large arity. If it
    // is nonempty, building its differential enforces the raw BP arity cap.
    let next_shape = Shape {
        arrows: shape.arrows + 1,
        ..shape
    };
    let next = bp_ext::basis(model, next_shape, t)?;
    let upper_weights = weights(model, mode, shape.height, &next)?;
    let outgoing = bp_ext::matrix(model, shape, monomials, &next)?;
    let upper_condition = (&outgoing, &upper_weights[..], f + r);
    let cycles = Lattice::build(model, prime, &current, f, common, Some(upper_condition))?;
    let upper = Lattice::build(
        model,
        prime,
        &current,
        f + 1,
        common,
        (page > 1).then_some(upper_condition),
    )?;
    let mut denominator = Basis::default();
    for i in 0..monomials.len() {
        let column = upper.generator(model, i)?;
        add_column(model, prime, &cycles, &mut denominator, column)?;
    }
    if let Some(previous_shape) = previous_shape {
        if !previous.is_empty() {
            let incoming = bp_ext::matrix(model, previous_shape, &previous, monomials)?;
            let lower = Lattice::build(
                model,
                prime,
                &lower_weights,
                lower_filtration,
                common,
                (page > 1).then_some((&incoming, &current[..], f)),
            )?;
            for i in 0..previous.len() {
                let column = lower.generator(model, i)?;
                let image = incoming.apply(model, &column)?;
                add_column(model, prime, &cycles, &mut denominator, image)?;
            }
        }
    }
    Ok((cycles, denominator, common))
}

fn add_column(
    model: &mut Model,
    prime: Prime,
    cycles: &Lattice,
    denominator: &mut Basis,
    column: Vec<i128>,
) -> Result<()> {
    let coordinates = cycles
        .encode(model, column)?
        .ok_or_else(|| Error::plain("BP filtered denominator is outside its cycle lattice"))?;
    denominator.insert(model, prime, residues(prime, coordinates))
}

fn residues(prime: Prime, coordinates: Vec<i128>) -> Vec<Coordinate> {
    let modulus = i128::from(prime.get());
    coordinates
        .into_iter()
        .enumerate()
        .filter_map(|(index, value)| {
            let coefficient = value.rem_euclid(modulus) as u32;
            (coefficient != 0).then_some((index, coefficient))
        })
        .collect()
}
