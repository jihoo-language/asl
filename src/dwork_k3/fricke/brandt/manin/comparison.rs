use std::rc::Rc;

use crate::error::{Error, Result};
use crate::lee_linear::Rational as Q;
use crate::model::{Model, Value};
use crate::syntax::Label;

use super::super::Basis;
use super::linear::{self as lin, Frame, Matrix, Vector, ZERO};
use super::{certificate, fraction, persist, rationals, Fraction, Operator, Space};

pub(crate) struct Comparison {
    pub(crate) manin: Label,
    pub(crate) brandt: Label,
    space: Rc<Space>,
    augmentation: Frame,
    forward: Matrix,
    backward: Matrix,
    determinant: Q,
    spectra: Vec<(&'static str, Vector)>,
}

impl Comparison {
    pub(super) fn new(
        model: &mut Model,
        manin: Label,
        brandt: Label,
        space: Rc<Space>,
        basis: &Basis,
    ) -> Result<Self> {
        let p = space.geometry.prime;
        if basis.p as usize != p {
            return Err(Error::plain(
                "Manin level must be twice the Brandt characteristic",
            ));
        }
        let n = basis.rank();
        let d = space.joint.vectors.len();
        if n != d + 1 {
            return Err(certificate(
                "cuspidal J+ W2+ and Brandt augmentation have equal dimension",
            ));
        }
        let augmentation = Frame::nullspace(model, n, vec![vec![lin::ONE; n]])?;
        let mut comparisons = Vec::new();
        for (ell, name) in [(3, "T3"), (5, "T5")] {
            if ell == p {
                continue;
            }
            let matrix = basis.operator(model, ell)?;
            let divisor_matrix: Matrix = (0..n)
                .map(|j| matrix.0.iter().map(|row| Q::integer(row[j])).collect())
                .collect();
            let hecke = space.operator(model, Operator::Hecke(ell))?;
            comparisons.push((
                name,
                space.joint.restrict(model, &hecke)?,
                augmentation.restrict(model, &divisor_matrix)?,
            ));
        }
        let mut frobenius = vec![vec![ZERO; n]; n];
        for (j, &i) in basis.frobenius.iter().enumerate() {
            frobenius[i][j] = lin::ONE;
        }
        let wp = space.operator(model, Operator::Atkin(p))?;
        let mut minus_wp = space.joint.restrict(model, &wp)?;
        for entry in minus_wp.iter_mut().flatten() {
            *entry = lin::neg(*entry)?;
        }
        comparisons.push((
            "-Wp / Frobenius",
            minus_wp,
            augmentation.restrict(model, &frobenius)?,
        ));
        let forward = find_intertwiner(model, &comparisons, d)?;
        let (backward, determinant) = lin::inverse(model, &forward)?
            .ok_or_else(|| certificate("invertible simultaneous Jacquet–Langlands intertwiner"))?;
        if lin::multiply(model, &forward, &backward)? != lin::identity(d)
            || lin::multiply(model, &backward, &forward)? != lin::identity(d)
        {
            return Err(certificate("the two comparison maps are inverse"));
        }
        let mut spectra = Vec::new();
        for (name, source, target) in &comparisons {
            if !intertwines(model, &forward, source, target)?
                || !intertwines(model, &backward, target, source)?
            {
                return Err(certificate(
                    "both directions intertwine the actual operators",
                ));
            }
            spectra.push((*name, lin::characteristic(model, source)?));
        }
        model.allocate(32 + 8 * d * d + 4 * n * d + 12 * d)?;
        Ok(Self {
            manin,
            brandt,
            space,
            augmentation,
            forward,
            backward,
            determinant,
            spectra,
        })
    }

    pub(crate) fn transfer(
        &self,
        model: &mut Model,
        to_brandt: bool,
        value: &Fraction,
    ) -> Result<Value> {
        let input = rationals(value)?;
        if to_brandt {
            let coordinate = self
                .space
                .joint
                .coordinates(model, &input)?
                .ok_or_else(|| {
                    Error::plain(
                        "Jacquet–Langlands transfer requires a cuspidal J+ W2+ Manin class",
                    )
                })?;
            let coordinate = lin::apply(model, &self.forward, &coordinate)?;
            let output = self.augmentation.expand(model, &coordinate)?;
            fraction(&output)?.persist(model)
        } else {
            let coordinate = self
                .augmentation
                .coordinates(model, &input)?
                .ok_or_else(|| {
                    Error::plain(
                        "Jacquet–Langlands transfer requires a Brandt divisor of degree zero",
                    )
                })?;
            let coordinate = lin::apply(model, &self.backward, &coordinate)?;
            let output = self.space.joint.expand(model, &coordinate)?;
            persist(model, &output)
        }
    }

    pub(super) fn describe(&self) -> String {
        let mut lines = vec![format!(
            "Jacquet–Langlands dimension {}; determinant = {}",
            self.forward.len(),
            lin::describe(&[self.determinant]).trim_matches(['[', ']'])
        )];
        for (name, characteristic) in &self.spectra {
            lines.push(format!(
                "{name}: characteristic (descending) = {}",
                lin::describe(characteristic)
            ));
        }
        lines.push(format!(
            "Manin to Brandt = {}",
            lin::describe_matrix(&self.forward)
        ));
        lines.push(format!(
            "Brandt to Manin = {}",
            lin::describe_matrix(&self.backward)
        ));
        lines.join("\n")
    }
}

type Operators = [(&'static str, Matrix, Matrix)];

fn find_intertwiner(model: &mut Model, operators: &Operators, d: usize) -> Result<Matrix> {
    if d == 0 {
        return Ok(Vec::new());
    }
    // When a generator is cyclic, the same polynomial in its iterates
    // gives the change of basis. Direct verification still certifies
    // simultaneous intertwining with the other operators.
    let (_, source, target) = &operators[0];
    for coefficient in [0, 1, -1, 2] {
        let mut a = source.clone();
        let mut b = target.clone();
        if coefficient != 0 {
            for ((left, right), (c, e)) in a
                .iter_mut()
                .zip(&mut b)
                .zip(operators[1].1.iter().zip(&operators[1].2))
            {
                lin::add_scaled(model, left, c, Q::integer(coefficient))?;
                lin::add_scaled(model, right, e, Q::integer(coefficient))?;
            }
        }
        if let Some(cyclic_a) = cyclic(model, &a)? {
            if let Some(cyclic_b) = cyclic(model, &b)? {
                let (inverse, _) = lin::inverse(model, &cyclic_a)?
                    .ok_or_else(|| certificate("cyclic module basis"))?;
                let candidate = lin::multiply(model, &cyclic_b, &inverse)?;
                if simultaneous(model, &candidate, operators)? {
                    return Ok(candidate);
                }
            }
        }
    }
    // Solve all homogeneous equations for the entries of the intertwiner
    // when the operators do not admit a cyclic separating generator.
    let mut equations = Vec::new();
    for (_, source, target) in operators {
        for i in 0..d {
            for j in 0..d {
                let mut row = vec![ZERO; d * d];
                for k in 0..d {
                    row[i * d + k] = row[i * d + k].add(source[k][j])?;
                    row[k * d + j] = lin::sub(row[k * d + j], target[i][k])?;
                }
                equations.push(row);
            }
        }
    }
    let frame = Frame::nullspace(model, d * d, equations)?;
    // Search actual invertible combinations. A failed bounded search makes
    // no assertion about abstract equivalence.
    let mut coefficients = vec![0usize; frame.vectors.len()];
    for _ in 0..4096 {
        let mut place = 0;
        while place < coefficients.len() && coefficients[place] == d {
            coefficients[place] = 0;
            place += 1;
        }
        if place == coefficients.len() {
            break;
        }
        coefficients[place] += 1;
        let weights = coefficients
            .iter()
            .map(|&c| Q::integer(c as i128))
            .collect::<Vector>();
        let candidate = frame.expand(model, &weights)?;
        let matrix = candidate.chunks(d).map(|row| row.to_vec()).collect();
        if lin::inverse(model, &matrix)?.is_some() && simultaneous(model, &matrix, operators)? {
            return Ok(matrix);
        }
    }
    Err(Error::plain(
        "bounded search did not find an invertible simultaneous Manin–Brandt intertwiner",
    ))
}

fn cyclic(model: &mut Model, matrix: &Matrix) -> Result<Option<Matrix>> {
    let n = matrix.len();
    let mut seeds = lin::identity(n);
    seeds.push((1..=n).map(|j| Q::integer(j as i128)).collect());
    for mut seed in seeds {
        let mut powers = Vec::with_capacity(n);
        for _ in 0..n {
            let next = lin::apply(model, matrix, &seed)?;
            powers.push(seed);
            seed = next;
        }
        let frame = lin::columns(&powers, n);
        if lin::inverse(model, &frame)?.is_some() {
            return Ok(Some(frame));
        }
    }
    Ok(None)
}

fn simultaneous(model: &mut Model, candidate: &Matrix, operators: &Operators) -> Result<bool> {
    for (_, source, target) in operators {
        if !intertwines(model, candidate, source, target)? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn intertwines(model: &mut Model, x: &Matrix, source: &Matrix, target: &Matrix) -> Result<bool> {
    Ok(lin::multiply(model, x, source)? == lin::multiply(model, target, x)?)
}
