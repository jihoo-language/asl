use crate::bp::Shape;
use crate::bp_ext;
use crate::bp_lattice::Lattice;
use crate::bp_polynomial::Monomial;
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
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) mode: Mode,
    pub(crate) shape: Shape,
    pub(crate) page: usize,
    pub(crate) t: usize,
    pub(crate) filtration: usize,
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
