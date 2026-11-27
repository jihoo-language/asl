use std::cmp::Ordering;
use std::rc::Rc;

use crate::bp;
use crate::bp_classes::{self, Grade, Support};
use crate::bp_ext;
use crate::bp_polynomial::Polynomial;
use crate::error::{Error, Result};
use crate::integer_smith::ColumnHermite;
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
        let mut lattice = ColumnHermite::from_diagonal(model, relations)?;
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
        let base = residue(model, &lattice, &base)?;
        let prime = self.target.shape.prime;
        let exponent = subgroup_exponent(model, &lattice, prime, relations)?;
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
    lattice: ColumnHermite,
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
        let reduced = residue(model, &self.lattice, candidate)?;
        Ok(reduced == self.base)
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        model.spend(
            self.lattice.cells() + other.lattice.cells() + self.base.len() + other.base.len(),
        )?;
        Ok((self.lattice.columns(), &self.base).cmp(&(other.lattice.columns(), &other.base)))
    }
}

fn residue(model: &mut Model, lattice: &ColumnHermite, vector: &Support) -> Result<Support> {
    let mut reduced = dense(model, lattice.columns().len(), vector)?;
    lattice.reduce(model, &mut reduced)?;
    Ok(Rc::from(
        reduced
            .into_iter()
            .enumerate()
            .filter(|&(_, x)| x != 0)
            .collect::<Vec<_>>(),
    ))
}

fn subgroup_exponent(
    model: &mut Model,
    lattice: &ColumnHermite,
    prime: u32,
    relations: &[i128],
) -> Result<usize> {
    let prime = i128::from(prime);
    let mut exponent = 0usize;
    for (i, (column, &relation)) in lattice.columns().iter().zip(relations).enumerate() {
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
fn dense(model: &mut Model, length: usize, support: &Support) -> Result<Vec<i128>> {
    model.spend(length + support.len())?;
    let mut vector = vec![0; length];
    for &(index, coefficient) in support.iter() {
        vector[index] = coefficient;
    }
    Ok(vector)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hermite_indeterminacy_handles_mixed_cyclic_orders() {
        let mut model = Model::new();
        let mut verify = || -> Result<()> {
            let relations = [4, 8];
            let mut one = ColumnHermite::from_diagonal(&mut model, &relations)?;
            one.insert(&mut model, vec![2, 2])?;
            one.canonicalize(&mut model)?;
            let mut two = ColumnHermite::from_diagonal(&mut model, &relations)?;
            two.insert(&mut model, vec![0, 4])?;
            two.insert(&mut model, vec![2, 6])?;
            two.canonicalize(&mut model)?;
            assert_eq!(one.columns(), two.columns());
            assert_eq!(subgroup_exponent(&mut model, &one, 2, &relations)?, 2);
            let origin = residue(&mut model, &one, &[(0, 1)].into())?;
            assert_eq!(residue(&mut model, &one, &[(0, 3), (1, 2)].into())?, origin);
            assert_ne!(residue(&mut model, &one, &[(0, 3)].into())?, origin);
            Ok(())
        };
        verify().unwrap_or_else(|error| panic!("{error}"));
    }
}
