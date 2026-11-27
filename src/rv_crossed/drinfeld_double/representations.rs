use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_cyclotomic;

use super::characters::Extension;
use super::{negate_phase, Quotient, Simple};

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Image {
    pub(crate) target: usize,
    pub(crate) phase: usize,
}

pub(super) struct Irrep {
    basis: Rc<InductionBasis>,
    pub(super) action: Vec<Vec<Image>>,
}

impl Irrep {
    pub(super) fn data(&self, name: &str, i: usize) -> String {
        format!(
            "rvDoubleRepData[{name},{i}] = isotropic={:?}, coset_representatives={:?}, dimension={}",
            self.basis.subgroup, self.basis.representatives, self.basis.representatives.len(),
        )
    }

    pub(super) fn inverse_image(&self, x: usize, column: usize, conductor: usize) -> Image {
        let (target, image) = self.action[x]
            .iter()
            .enumerate()
            .find(|(_, image)| image.target == column)
            .expect("verified monomial permutation");
        Image {
            target,
            phase: negate_phase(conductor, image.phase),
        }
    }
}

struct InductionBasis {
    subgroup: Vec<usize>,
    representatives: Vec<usize>,
}

pub(super) fn induce_all(
    model: &mut Model,
    group: &Quotient,
    alpha: &[Vec<u32>],
    radical: &[usize],
    simples: &[Simple],
    p: Prime,
    exponent: usize,
) -> Result<Vec<Irrep>> {
    let n = group.order();
    let conductor = p.get() as usize * exponent;
    let mut subgroup = radical.to_vec();
    let mut extensions = Vec::new();
    model.spend(n * n)?;
    for x in 1..n {
        if !subgroup.contains(&x) && subgroup.iter().all(|&y| alpha[x][y] == alpha[y][x]) {
            let step = Extension::new(group, alpha, &subgroup, x, p);
            let next = step.subgroup();
            extensions.push((subgroup, step));
            subgroup = next;
        }
    }
    subgroup.sort_unstable();
    if subgroup.len() * subgroup.len() != n * radical.len() {
        return Err(Error::plain(
            "commutator isotropic subgroup has incorrect order",
        ));
    }
    let mut representatives = Vec::new();
    let mut decomposition = vec![(usize::MAX, usize::MAX); n];
    for t in 0..n {
        if decomposition[t].0 != usize::MAX {
            continue;
        }
        let index = representatives.len();
        representatives.push(t);
        for &l in &subgroup {
            decomposition[group.mul(t, l)] = (index, l);
        }
    }
    let basis = Rc::new(InductionBasis {
        subgroup,
        representatives,
    });
    model.allocate(basis.subgroup.len() + basis.representatives.len() + 9)?;
    let mut answer = Vec::with_capacity(simples.len());
    for simple in simples {
        if basis.representatives.len() != simple.dimension {
            return Err(Error::plain(
                "induced projective representation has incorrect dimension",
            ));
        }
        let mut eta = simple.phases.clone();
        for (known, step) in &extensions {
            model.spend(n)?;
            eta = step.lift(alpha, known, &eta, exponent, conductor, 0)?;
        }
        model.spend(basis.subgroup.len() * basis.subgroup.len())?;
        for &x in &basis.subgroup {
            for &y in &basis.subgroup {
                if (eta[x].unwrap() + eta[y].unwrap()) % conductor
                    != (eta[group.mul(x, y)].unwrap() + exponent * alpha[x][y] as usize) % conductor
                {
                    return Err(Error::plain(
                        "maximal isotropic character fails a projective product",
                    ));
                }
            }
        }
        let action = (0..n)
            .map(|x| {
                basis
                    .representatives
                    .iter()
                    .map(|&t| {
                        let (target, l) = decomposition[group.mul(x, t)];
                        let next_rep = basis.representatives[target];
                        let phase = (eta[l].expect("induced character is defined on L")
                            + conductor
                            + exponent * alpha[x][t] as usize
                            - exponent * alpha[next_rep][l] as usize)
                            % conductor;
                        Image { target, phase }
                    })
                    .collect()
            })
            .collect();
        let irrep = Irrep {
            basis: basis.clone(),
            action,
        };
        verify(model, group, alpha, simple, &irrep, exponent, conductor)?;
        model.allocate(5 + 3 * n + 2 * n * simple.dimension)?;
        answer.push(irrep);
    }
    Ok(answer)
}

fn verify(
    model: &mut Model,
    group: &Quotient,
    alpha: &[Vec<u32>],
    simple: &Simple,
    irrep: &Irrep,
    exponent: usize,
    conductor: usize,
) -> Result<()> {
    let n = group.order();
    model.spend(n * n * simple.dimension + n * simple.dimension)?;
    for (x, row) in irrep.action.iter().enumerate() {
        for (y, second) in irrep.action.iter().enumerate() {
            let product = &irrep.action[group.mul(x, y)];
            let multiplier = exponent * alpha[x][y] as usize;
            for (b, c) in second.iter().zip(product) {
                let a = row[b.target];
                if a.target != c.target
                    || (a.phase + b.phase) % conductor != (c.phase + multiplier) % conductor
                {
                    return Err(Error::plain(
                        "induced matrices fail the twisted projective representation law",
                    ));
                }
            }
        }
        if let Some(phase) = simple.phases[x] {
            if row
                .iter()
                .enumerate()
                .any(|(i, image)| image.target != i || image.phase != phase)
            {
                return Err(Error::plain(
                    "induced matrix disagrees with its radical character",
                ));
            }
        } else {
            let mut trace = vec![0; conductor];
            for (i, image) in row.iter().enumerate() {
                if image.target == i {
                    trace[image.phase] += 1;
                }
            }
            if !rv_cyclotomic::reduced(model, conductor, trace)?.is_empty() {
                return Err(Error::plain(
                    "induced matrix has nonzero trace outside the radical",
                ));
            }
        }
    }
    Ok(())
}
