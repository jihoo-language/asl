use std::collections::VecDeque;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
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
    pub(super) fn isotropic(&self) -> &[usize] {
        &self.basis.subgroup
    }

    pub(super) fn representatives(&self) -> &[usize] {
        &self.basis.representatives
    }

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
    alpha: &[Vec<usize>],
    radical: &[usize],
    simples: &[Simple],
    conductor: usize,
) -> Result<Vec<Irrep>> {
    let n = group.order();
    let mut subgroup = radical.to_vec();
    let mut extensions = Vec::new();
    model.spend(n * n)?;
    for x in 1..n {
        if !subgroup.contains(&x) && subgroup.iter().all(|&y| alpha[x][y] == alpha[y][x]) {
            let step = Extension::new(group, alpha, &subgroup, x, conductor);
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
            eta = step.lift(alpha, known, &eta, conductor, 0)?;
        }
        model.spend(basis.subgroup.len() * basis.subgroup.len())?;
        for &x in &basis.subgroup {
            for &y in &basis.subgroup {
                if (eta[x].unwrap() + eta[y].unwrap()) % conductor
                    != (eta[group.mul(x, y)].unwrap() + alpha[x][y]) % conductor
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
                            + alpha[x][t]
                            - alpha[next_rep][l])
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
        verify(model, group, alpha, simple, &irrep, conductor)?;
        model.allocate(5 + 3 * n + 2 * n * simple.dimension)?;
        answer.push(irrep);
    }
    Ok(answer)
}

fn verify(
    model: &mut Model,
    group: &Quotient,
    alpha: &[Vec<usize>],
    simple: &Simple,
    irrep: &Irrep,
    conductor: usize,
) -> Result<()> {
    verify_action(model, group, alpha, &irrep.action, conductor)?;
    model.spend(group.order() * simple.dimension)?;
    for (x, row) in irrep.action.iter().enumerate() {
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

pub(super) fn verify_action(
    model: &mut Model,
    group: &Quotient,
    alpha: &[Vec<usize>],
    action: &[Vec<Image>],
    conductor: usize,
) -> Result<()> {
    let n = group.order();
    model.spend(n * n * action[0].len())?;
    for (x, row) in action.iter().enumerate() {
        for (y, second) in action.iter().enumerate() {
            let product = &action[group.mul(x, y)];
            let multiplier = alpha[x][y];
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
    }
    Ok(())
}

pub(super) fn tensor_action(
    model: &mut Model,
    left: &[Vec<Image>],
    right: &[Vec<Image>],
    phases: &[usize],
    conductor: usize,
) -> Result<Vec<Vec<Image>>> {
    let width = right[0].len();
    model.spend(left.len() * left[0].len() * width)?;
    Ok(left
        .iter()
        .zip(right)
        .zip(phases)
        .map(|((a, b), &phase)| {
            a.iter()
                .flat_map(|a| {
                    b.iter().map(move |b| Image {
                        target: a.target * width + b.target,
                        phase: (phase + a.phase + b.phase) % conductor,
                    })
                })
                .collect()
        })
        .collect())
}

pub(super) type IntertwinerOrbit = Vec<(usize, usize)>;

pub(super) fn intertwiner_orbits(
    model: &mut Model,
    tensor: &[Vec<Image>],
    source: &[Vec<Image>],
    conductor: usize,
) -> Result<Vec<IntertwinerOrbit>> {
    let rows = tensor[0].len();
    let cols = source[0].len();
    let size = rows * cols;
    model.spend(3 * tensor.len() * size)?;
    let mut phases = vec![None; size];
    let mut orbits = Vec::new();
    for root in 0..size {
        if phases[root].is_some() {
            continue;
        }
        let mut orbit = Vec::new();
        let mut pending = VecDeque::from([root]);
        phases[root] = Some(0);
        let mut consistent = true;
        while let Some(i) = pending.pop_front() {
            let (row, column) = (i / cols, i % cols);
            let v = phases[i].expect("known intertwiner phase");
            orbit.push((i, v));
            for (a, b) in tensor.iter().zip(source) {
                let (a, b) = (a[row], b[column]);
                let j = a.target * cols + b.target;
                let w = (v + a.phase + conductor - b.phase) % conductor;
                if let Some(existing) = phases[j] {
                    consistent &= existing == w;
                } else {
                    phases[j] = Some(w);
                    pending.push_back(j);
                }
            }
        }
        if consistent {
            orbits.push(orbit);
        }
    }
    Ok(orbits)
}
