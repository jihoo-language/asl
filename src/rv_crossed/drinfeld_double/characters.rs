use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;

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
) -> Result<Vec<Vec<u32>>> {
    let n = group.order();
    let p = post.prime;
    let mut alpha = vec![vec![0; n]; n];
    for (x, row) in alpha.iter_mut().enumerate().skip(1) {
        for (y, value) in row.iter_mut().enumerate().skip(1) {
            let a = post.coefficient(model, k, &[g, x, y])?;
            let b = post.coefficient(model, k, &[x, y, g])?;
            let c = post.coefficient(model, k, &[x, g, y])?;
            *value = p.add(p.add(a, b), p.negate(c));
        }
    }
    model.spend(n * n * n)?;
    for x in 0..n {
        for y in 0..n {
            let xy = group.mul(x, y);
            for z in 0..n {
                if p.add(alpha[x][y], alpha[xy][z]) != p.add(alpha[y][z], alpha[x][group.mul(y, z)])
                {
                    return Err(Error::plain(
                        "transgressed double multiplier fails its 2-cocycle equation",
                    ));
                }
            }
        }
    }
    Ok(alpha)
}

pub(super) fn radical(model: &mut Model, alpha: &[Vec<u32>]) -> Result<Vec<usize>> {
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
    alpha: &[Vec<u32>],
    radical: &[usize],
    prime: Prime,
    exponent: usize,
) -> Result<Vec<Simple>> {
    let conductor = exponent * prime.get() as usize;
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
            let extension = Extension::new(group, alpha, &known, r, prime);
            model.spend(1 + characters.len() * extension.powers.len() * (1 + n))?;
            characters = extension.characters(alpha, &known, characters, exponent, conductor)?;
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
                if (a + b) % conductor != (xy + exponent * alpha[x][y] as usize) % conductor {
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
struct Extension {
    powers: Vec<usize>,
    carries: Vec<u32>,
    cosets: Vec<Vec<usize>>,
    return_to: usize,
    total_carry: u32,
}

impl Extension {
    fn new(group: &Quotient, alpha: &[Vec<u32>], known: &[usize], r: usize, p: Prime) -> Self {
        let mut powers = vec![0];
        let mut carries = vec![0];
        let mut current = 0;
        let mut carry = 0;
        loop {
            carry = p.add(carry, alpha[current][r]);
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

    fn characters(
        &self,
        alpha: &[Vec<u32>],
        known: &[usize],
        previous: Vec<Vec<Option<usize>>>,
        exponent: usize,
        conductor: usize,
    ) -> Result<Vec<Vec<Option<usize>>>> {
        let m = self.powers.len();
        let mut result = Vec::with_capacity(m * previous.len());
        for phases in previous {
            let rhs = (phases[self.return_to].expect("power lies in known subgroup")
                + exponent * self.total_carry as usize)
                % conductor;
            if !conductor.is_multiple_of(m) || !rhs.is_multiple_of(m) {
                return Err(Error::plain(
                    "double projective character power has no lift at its conductor",
                ));
            }
            for choice in 0..m {
                let generator_phase = rhs / m + choice * (conductor / m);
                let mut lifted = phases.clone();
                for (&x, coset) in known.iter().zip(&self.cosets) {
                    let phase = phases[x].expect("known subgroup");
                    for (j, &xy) in coset.iter().enumerate() {
                        let correction =
                            exponent * (self.carries[j] + alpha[x][self.powers[j]]) as usize;
                        lifted[xy] = Some(
                            (phase + j * generator_phase + conductor * 2 - correction) % conductor,
                        );
                    }
                }
                result.push(lifted);
            }
        }
        Ok(result)
    }
}
