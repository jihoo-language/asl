use std::cmp::Ordering;
use std::mem;
use std::rc::Rc;

use crate::bp;
use crate::bp_classes::{self, Grade, Support};
use crate::bp_ext;
use crate::bp_polynomial::Polynomial;
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};

type Inputs<'a> = [(TypeId, &'a Support); 3];

pub(crate) struct Layout {
    inputs: [Grade; 3],
    ab: Grade,
    bc: Grade,
    u: Grade,
    v: Grade,
    pub(crate) target: Grade,
}

impl Layout {
    pub(crate) fn new(inputs: [Grade; 3]) -> Result<Self> {
        let [a, b, c] = inputs;
        if inputs.iter().any(|grade| grade.shape.arrows == 0) {
            return Err(Error::plain("bpTriple requires positive cobar degrees"));
        }
        let ab = a.cup(b)?;
        let bc = b.cup(c)?;
        let u = Grade {
            shape: bp_ext::previous_shape(ab.shape)?,
            ..ab
        };
        let v = Grade {
            shape: bp_ext::previous_shape(bc.shape)?,
            ..bc
        };
        let target = a.cup(v)?;
        Ok(Self {
            inputs,
            ab,
            bc,
            u,
            v,
            target,
        })
    }

    pub(crate) fn defined(&self, model: &mut Model, inputs: Inputs<'_>) -> Result<bool> {
        for (left, right, grade) in [
            (inputs[0], inputs[1], self.ab),
            (inputs[1], inputs[2], self.bc),
        ] {
            let ty = model.bp_ext_type(grade)?;
            let Value::BPExt(product) = bp_classes::cup(model, left, right, ty)? else {
                unreachable!()
            };
            if !product.is_empty() {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub(crate) fn build(&self, model: &mut Model, inputs: Inputs<'_>) -> Result<Coset> {
        let [a_shape, b_shape, c_shape] = self.inputs.map(|grade| grade.shape);
        let lift = |model: &mut Model, (ty, support): (TypeId, &Support)| {
            if support.is_empty() {
                return Ok(Polynomial::default());
            }
            model.bp_ext(ty)?.representative(model, support)
        };
        let a = lift(model, inputs[0])?;
        let b = lift(model, inputs[1])?;
        let ab = bp::cup(model, a_shape, b_shape, &a, &b)?;
        require_boundary(model, self.ab, &ab, "A cup B")?;
        let c = lift(model, inputs[2])?;
        let bc = bp::cup(model, b_shape, c_shape, &b, &c)?;
        require_boundary(model, self.bc, &bc, "B cup C")?;
        let u = bp_ext::null(model, self.ab.shape, self.ab.t, &ab)?;
        let v = bp_ext::null(model, self.bc.shape, self.bc.t, &bc)?;
        let layout = bp::Massey::new([a_shape, b_shape, c_shape, self.u.shape, self.v.shape])?;
        let representative = layout.apply(model, [&a, &b, &c, &u, &v])?;
        let ty = model.bp_ext_type(self.target)?;
        let target = model.bp_ext(ty)?;
        debug_assert_eq!(target.invariants.free, 0);
        let relations = &target.invariants.torsion;
        let mut lattice = Lattice::new(model, relations)?;
        if !relations.is_empty() {
            for (grade, fixed, on_left) in [(self.v, &a, true), (self.u, &c, false)] {
                if fixed.is_zero() {
                    continue;
                }
                let ty = model.bp_ext_type(grade)?;
                let group = model.bp_ext(ty)?;
                for index in 0..group.rank() {
                    let generator = group.representative(model, &[(index, 1)])?;
                    let image = if on_left {
                        bp::cup(model, a_shape, grade.shape, fixed, &generator)?
                    } else {
                        bp::cup(model, grade.shape, c_shape, &generator, fixed)?
                    };
                    if !image.is_zero() {
                        let coordinates = target.project(model, &image)?;
                        if !coordinates.is_empty() {
                            let vector = dense(model, relations.len(), &coordinates)?;
                            lattice.insert(model, vector)?;
                        }
                    }
                }
            }
        }
        lattice.canonicalize(model)?;
        let base = if representative.is_zero() {
            Support::from([])
        } else {
            target.project(model, &representative)?
        };
        let base = lattice.residue(model, &base)?;
        let prime = self.target.shape.prime;
        let exponent = lattice.subgroup_exponent(model, prime, relations)?;
        model.allocate(lattice.cells() + 3 * base.len())?;
        Ok(Coset {
            lattice,
            base,
            prime,
            exponent,
        })
    }
}

fn require_boundary(
    model: &mut Model,
    grade: Grade,
    cochain: &Polynomial,
    product: &str,
) -> Result<()> {
    if cochain.is_zero() {
        return Ok(());
    }
    let ty = model.bp_ext_type(grade)?;
    if !model.bp_ext(ty)?.project(model, cochain)?.is_empty() {
        return Err(Error::plain(format!(
            "bpTriple is undefined because {product} is nonzero"
        )));
    }
    Ok(())
}

pub(crate) struct Coset {
    lattice: Lattice,
    base: Support,
    prime: u32,
    exponent: usize,
}

impl Coset {
    pub(crate) fn base(&self) -> Value {
        Value::BPExt(self.base.clone())
    }

    pub(crate) fn size(&self) -> String {
        match self.exponent {
            0 => "1".to_owned(),
            exponent => format!("{}^{exponent}", self.prime),
        }
    }

    pub(crate) fn contains(&self, model: &mut Model, candidate: &Support) -> Result<bool> {
        let reduced = self.lattice.residue(model, candidate)?;
        Ok(reduced == self.base)
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        model.spend(
            self.lattice.cells() + other.lattice.cells() + self.base.len() + other.base.len(),
        )?;
        Ok((&self.lattice.columns, &self.base).cmp(&(&other.lattice.columns, &other.base)))
    }
}

// Full-rank integral column Hermite basis. Starting with the cyclic
// relations ensures every pivot stays positive. The leading entry of
// column i is at row i, and later entries are reduced by later columns.
struct Lattice {
    columns: Vec<Vec<i128>>,
}

impl Lattice {
    fn new(model: &mut Model, relations: &[i128]) -> Result<Self> {
        let n = relations.len();
        let cells = n
            .checked_mul(3 + 2 * n)
            .filter(|&size| size <= MAX_PIVOT_CELLS)
            .ok_or_else(|| {
                Error::plain("BP triple indeterminacy lattice exceeds the storage limit")
            })?;
        model.spend(cells)?;
        let mut columns = vec![vec![0; n]; n];
        for (i, &relation) in relations.iter().enumerate() {
            columns[i][i] = relation;
        }
        Ok(Self { columns })
    }

    fn cells(&self) -> usize {
        let n = self.columns.len();
        n * (3 + 2 * n)
    }

    fn insert(&mut self, model: &mut Model, mut vector: Vec<i128>) -> Result<()> {
        model.spend(self.columns.len())?;
        for (i, column) in self.columns.iter_mut().enumerate() {
            while vector[i] != 0 {
                reduce_at(model, &mut vector[i..], &column[i..])?;
                if vector[i] != 0 {
                    mem::swap(&mut vector, column);
                }
            }
        }
        Ok(())
    }

    fn canonicalize(&mut self, model: &mut Model) -> Result<()> {
        for i in (0..self.columns.len()).rev() {
            let (prefix, suffix) = self.columns.split_at_mut(i + 1);
            let current = &mut prefix[i];
            for (j, column) in ((i + 1)..).zip(suffix) {
                reduce_at(model, &mut current[j..], &column[j..])?;
            }
        }
        Ok(())
    }

    fn residue(&self, model: &mut Model, vector: &Support) -> Result<Support> {
        let mut reduced = dense(model, self.columns.len(), vector)?;
        for (i, column) in self.columns.iter().enumerate() {
            reduce_at(model, &mut reduced[i..], &column[i..])?;
        }
        Ok(Rc::from(
            reduced
                .into_iter()
                .enumerate()
                .filter(|&(_, x)| x != 0)
                .collect::<Vec<_>>(),
        ))
    }

    fn subgroup_exponent(
        &self,
        model: &mut Model,
        prime: u32,
        relations: &[i128],
    ) -> Result<usize> {
        let prime = i128::from(prime);
        let mut exponent = 0usize;
        for (i, (column, &relation)) in self.columns.iter().zip(relations).enumerate() {
            let pivot = column[i];
            if relation % pivot != 0 {
                return Err(Error::plain("BP triple lattice lost a cyclic relation"));
            }
            let mut orders = relation / pivot;
            while orders > 1 {
                model.spend(1)?;
                if orders % prime != 0 {
                    return Err(Error::plain("BP triple indeterminacy is not p-primary"));
                }
                orders /= prime;
                exponent += 1;
            }
        }
        Ok(exponent)
    }
}

fn dense(model: &mut Model, length: usize, support: &Support) -> Result<Vec<i128>> {
    model.spend(length + support.len())?;
    let mut vector = vec![0; length];
    for &(index, coefficient) in support.iter() {
        vector[index] = coefficient;
    }
    Ok(vector)
}

fn reduce_at(model: &mut Model, vector: &mut [i128], column: &[i128]) -> Result<()> {
    model.spend(1)?;
    let q = vector[0].div_euclid(column[0]);
    if q == 0 {
        return Ok(());
    }
    vector[0] = vector[0].rem_euclid(column[0]);
    model.spend(vector.len() - 1)?;
    let overflow = || Error::plain("BP triple lattice arithmetic exceeds the signed i128 range");
    for (value, &entry) in vector[1..].iter_mut().zip(&column[1..]) {
        *value = value
            .checked_sub(q.checked_mul(entry).ok_or_else(overflow)?)
            .ok_or_else(overflow)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hermite_indeterminacy_handles_mixed_cyclic_orders() {
        let mut model = Model::new();
        let mut verify = || -> Result<()> {
            let relations = [4, 8];
            let mut one = Lattice::new(&mut model, &relations)?;
            one.insert(&mut model, vec![2, 2])?;
            one.canonicalize(&mut model)?;
            let mut two = Lattice::new(&mut model, &relations)?;
            two.insert(&mut model, vec![0, 4])?;
            two.insert(&mut model, vec![2, 6])?;
            two.canonicalize(&mut model)?;
            assert_eq!(one.columns, two.columns);
            assert_eq!(one.subgroup_exponent(&mut model, 2, &relations)?, 2);
            let origin = one.residue(&mut model, &[(0, 1)].into())?;
            assert_eq!(one.residue(&mut model, &[(0, 3), (1, 2)].into())?, origin);
            assert_ne!(one.residue(&mut model, &[(0, 3)].into())?, origin);
            Ok(())
        };
        verify().unwrap_or_else(|error| panic!("{error}"));
    }
}
