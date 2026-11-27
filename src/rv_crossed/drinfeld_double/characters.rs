use crate::error::{Error, Result};
use crate::model::Model;

use super::super::path_integral::Quotient;
use super::super::postnikov::{Cochain, Postnikov};
use super::Simple;

pub(super) fn exponent(model: &mut Model, group: &Quotient) -> Result<usize> {
    let mut exponent = 1;
    for x in 0..group.order() {
        let mut order = 1;
        let mut power = x;
        while power != 0 {
            model.spend(1)?;
            power = group.mul(power, x);
            order += 1;
        }
        let common = crate::rv_lens::gcd(exponent as u128, order as u128) as usize;
        exponent = exponent / common * order;
    }
    Ok(exponent)
}

pub(super) fn transgress(
    model: &mut Model,
    group: &Quotient,
    post: &Postnikov,
    k: &Cochain,
    g: usize,
    exponent: usize,
) -> Result<Vec<Vec<usize>>> {
    let n = group.order();
    let p = post.prime;
    let mut alpha = vec![vec![0; n]; n];
    for (x, row) in alpha.iter_mut().enumerate().skip(1) {
        for (y, value) in row.iter_mut().enumerate().skip(1) {
            let a = post.coefficient(model, k, &[g, x, y])?;
            let b = post.coefficient(model, k, &[x, y, g])?;
            let c = post.coefficient(model, k, &[x, g, y])?;
            *value = exponent * p.add(p.add(a, b), p.negate(c)) as usize;
        }
    }
    verify_multiplier(model, group, &alpha, exponent * p.get() as usize)?;
    Ok(alpha)
}

pub(super) fn verify_multiplier(
    model: &mut Model,
    group: &Quotient,
    alpha: &[Vec<usize>],
    conductor: usize,
) -> Result<()> {
    let n = group.order();
    model.spend(n * n * n)?;
    if (0..n).any(|i| alpha[i][0] != 0 || alpha[0][i] != 0) {
        return Err(Error::plain("projective multiplier must be normalized"));
    }
    for x in 0..n {
        for y in 0..n {
            let xy = group.mul(x, y);
            for z in 0..n {
                if (alpha[x][y] + alpha[xy][z]) % conductor
                    != (alpha[y][z] + alpha[x][group.mul(y, z)]) % conductor
                {
                    return Err(Error::plain(
                        "projective multiplier fails its 2-cocycle equation",
                    ));
                }
            }
        }
    }
    Ok(())
}

pub(super) fn radical(model: &mut Model, alpha: &[Vec<usize>]) -> Result<Vec<usize>> {
    let n = alpha.len();
    model.spend(n * n)?;
    Ok((0..n)
        .filter(|&x| (0..n).all(|y| alpha[x][y] == alpha[y][x]))
        .collect())
}

pub(super) fn extend(
    model: &mut Model,
    group: &Quotient,
    flux: usize,
    alpha: &[Vec<usize>],
    radical: &[usize],
    conductor: usize,
) -> Result<Vec<Simple>> {
    let n = group.order();
    let square = n / radical.len();
    let dimension = square.isqrt();
    if dimension * dimension * radical.len() != n {
        return Err(Error::plain(
            "double commutator radical has nonsquare index",
        ));
    }

    let mut known = vec![0];
    let mut initial = vec![None; n];
    initial[0] = Some(0);
    let mut characters = vec![initial];
    for &r in radical.iter().skip(1) {
        if !known.contains(&r) {
            let extension = Extension::new(group, alpha, &known, r, conductor);
            model.spend(1 + characters.len() * extension.powers.len() * (1 + n))?;
            characters = extension.characters(alpha, &known, characters, conductor)?;
            known = extension.cosets.into_iter().flatten().collect();
        }
    }
    if characters.len() != radical.len() {
        return Err(Error::plain(
            "double character count disagrees with the radical order",
        ));
    }
    model.spend(characters.len() * radical.len() * radical.len())?;
    for phases in &characters {
        for &x in radical {
            for &y in radical {
                let a = phases[x].expect("radical character support");
                let b = phases[y].expect("radical character support");
                let xy = phases[group.mul(x, y)].expect("radical is a subgroup");
                if (a + b) % conductor != (xy + alpha[x][y]) % conductor {
                    return Err(Error::plain(
                        "double central character fails a projective product",
                    ));
                }
            }
        }
    }
    Ok(characters
        .into_iter()
        .map(|phases| Simple {
            flux,
            dimension,
            phases,
        })
        .collect())
}

// Adjoin a generator r to the subgroup already covered. If r^m is the
// first power back in that subgroup, each character has precisely m lifts.
pub(super) struct Extension {
    powers: Vec<usize>,
    carries: Vec<usize>,
    cosets: Vec<Vec<usize>>,
    return_to: usize,
    total_carry: usize,
}

impl Extension {
    pub(super) fn new(
        group: &Quotient,
        alpha: &[Vec<usize>],
        known: &[usize],
        r: usize,
        conductor: usize,
    ) -> Self {
        let mut powers = vec![0];
        let mut carries = vec![0];
        let mut current = 0;
        let mut carry = 0;
        loop {
            carry = (carry + alpha[current][r]) % conductor;
            current = group.mul(current, r);
            if known.contains(&current) {
                break;
            }
            powers.push(current);
            carries.push(carry);
        }
        let cosets = known
            .iter()
            .map(|&x| powers.iter().map(|&y| group.mul(x, y)).collect())
            .collect();
        Self {
            powers,
            carries,
            cosets,
            return_to: current,
            total_carry: carry,
        }
    }

    pub(super) fn subgroup(&self) -> Vec<usize> {
        self.cosets.iter().flatten().copied().collect()
    }

    fn characters(
        &self,
        alpha: &[Vec<usize>],
        known: &[usize],
        previous: Vec<Vec<Option<usize>>>,
        conductor: usize,
    ) -> Result<Vec<Vec<Option<usize>>>> {
        let m = self.powers.len();
        let mut result = Vec::with_capacity(m * previous.len());
        for phases in previous {
            for choice in 0..m {
                result.push(self.lift(alpha, known, &phases, conductor, choice)?);
            }
        }
        Ok(result)
    }

    pub(super) fn lift(
        &self,
        alpha: &[Vec<usize>],
        known: &[usize],
        phases: &[Option<usize>],
        conductor: usize,
        choice: usize,
    ) -> Result<Vec<Option<usize>>> {
        let m = self.powers.len();
        let rhs = (phases[self.return_to].expect("power lies in known subgroup")
            + self.total_carry)
            % conductor;
        if !conductor.is_multiple_of(m) || !rhs.is_multiple_of(m) {
            return Err(Error::plain(
                "double projective character power has no lift at its conductor",
            ));
        }
        let generator_phase = rhs / m + choice * (conductor / m);
        let mut lifted = phases.to_vec();
        for (&x, coset) in known.iter().zip(&self.cosets) {
            let phase = phases[x].expect("known subgroup");
            for (j, &xy) in coset.iter().enumerate() {
                let correction = self.carries[j] + alpha[x][self.powers[j]];
                lifted[xy] =
                    Some((phase + j * generator_phase + conductor * 2 - correction) % conductor);
            }
        }
        Ok(lifted)
    }
}
