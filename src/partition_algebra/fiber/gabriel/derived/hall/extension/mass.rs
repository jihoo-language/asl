use std::collections::BTreeMap;

use crate::elliptic_modular::decimal;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::wide_integer::Natural;

use super::super::{Hall, Record, Shape};
use super::{digits, Groupoid, CLASS_LIMIT};

type Integer = Natural<64>;

fn arithmetic() -> Error {
    Error::plain("partition extension groupoid cardinality exceeds the exact 2048-bit bound")
}
fn from(x: u128) -> Integer {
    Integer::from_small(x)
}
fn product(a: Integer, b: Integer) -> Result<Integer> {
    a.checked_product(b).ok_or_else(arithmetic)
}
fn show(x: Integer) -> String {
    decimal(Wide {
        magnitude: x,
        negative: false,
    })
}

fn power(model: &mut Model, q: u32, dim: usize) -> Result<Integer> {
    let mut out = from(1);
    for _ in 0..dim {
        model.spend(64)?;
        out = out.multiply_u128(u128::from(q)).ok_or_else(arithmetic)?;
    }
    Ok(out)
}

fn fraction(model: &mut Model, a: Integer, b: Integer) -> Result<String> {
    let d = a.gcd(model, b)?;
    let a = a.divide_exact(model, d)?.ok_or_else(arithmetic)?;
    let b = b.divide_exact(model, d)?.ok_or_else(arithmetic)?;
    if b == from(1) {
        Ok(show(a))
    } else {
        Ok(format!("{}/{}", show(a), show(b)))
    }
}

impl Hall {
    pub(super) fn intertwiner_dimension(&self, m: &Shape, n: &Shape) -> usize {
        let mut dimension = 0;
        for &(i, l, a) in &m.0 {
            for &(j, k, b) in &n.0 {
                let paths = &self.path_heads[usize::from(j)][..usize::from(k)];
                let count = paths
                    .iter()
                    .enumerate()
                    .filter(|&(t, &head)| {
                        head == usize::from(i) && t + usize::from(l) >= usize::from(k)
                    })
                    .count();
                dimension += count * usize::from(a) * usize::from(b);
            }
        }
        dimension
    }

    pub(super) fn projective_ext_dimension(&self, m: &Shape, n: &Shape) -> usize {
        let mut dimension = 0;
        for &(head, length, copies) in &m.0 {
            let head = usize::from(head);
            let length = usize::from(length);
            let projective = &self.path_heads[head];
            if length == projective.len() {
                continue;
            }
            let first = projective[length];
            let syzygy_length = projective.len() - length;
            for &(target, target_length, count) in &n.0 {
                let paths = &self.path_heads[usize::from(target)][..usize::from(target_length)];
                for (layer, &vertex) in paths.iter().enumerate() {
                    if vertex != first || layer + syzygy_length < paths.len() {
                        continue;
                    }
                    // Restriction of a projective map to the syzygy is
                    // right multiplication by the path at this cutoff.
                    let restricts = layer.checked_sub(length).is_some_and(|start| {
                        paths[start] == head && start + projective.len() >= paths.len()
                    });
                    if !restricts {
                        dimension += usize::from(copies) * usize::from(count);
                    }
                }
            }
        }
        dimension
    }

    pub(crate) fn hom_report(&self, model: &mut Model, m: &Record, n: &Record) -> Result<String> {
        let dim = self.intertwiner_dimension(m.module(), n.module());
        let q = self.field.order();
        let cardinality = if (q.ilog2() as usize + 1) * dim + 1 >= 2048 {
            format!("{q}^{dim}")
        } else {
            show(power(model, q, dim)?)
        };
        Ok(format!(
            "rvPartitionHallHom[q={q}] = dimension={dim}, cardinality={cardinality}"
        ))
    }

    pub(crate) fn automorphism_report(&self, model: &mut Model, m: &Record) -> Result<String> {
        let shape = m.module();
        let dim = self.intertwiner_dimension(shape, shape);
        let q = self.field.order();
        let semisimple = semisimple_dimension(shape);
        let radical = dim - semisimple;
        let size = if (q.ilog2() as usize + 1) * dim + 1 >= 2048 {
            let factors = shape
                .0
                .iter()
                .map(|&(_, _, m)| format!("GL({m},{q})"))
                .collect::<Vec<_>>()
                .join("*");
            format!("{q}^{radical}*{factors}")
        } else {
            show(self.aut(model, shape)?)
        };
        Ok(format!("rvPartitionHallAut[q={q}] = endDimension={dim}, radicalDimension={radical}, cardinality={size}"))
    }

    fn aut(&self, model: &mut Model, module: &Shape) -> Result<Integer> {
        let dim = self.intertwiner_dimension(module, module);
        let radical = dim
            .checked_sub(semisimple_dimension(module))
            .ok_or_else(|| {
                Error::plain("partition Hall endomorphism radical dimension was negative")
            })?;
        let q = self.field.order();
        let mut value = power(model, q, radical)?;
        for &(_, _, m) in &module.0 {
            let order = power(model, q, usize::from(m))?;
            let mut smaller = from(1);
            for _ in 0..m {
                model.spend(64)?;
                value = product(value, order.difference(smaller))?;
                smaller = smaller
                    .multiply_u128(u128::from(q))
                    .ok_or_else(arithmetic)?;
            }
        }
        Ok(value)
    }
}

fn semisimple_dimension(shape: &Shape) -> usize {
    shape.0.iter().map(|&(_, _, m)| usize::from(m).pow(2)).sum()
}

impl Groupoid {
    pub(crate) fn groupoid_report(&self, model: &mut Model) -> Result<String> {
        let q = self.hall.field.order();
        let cocycles = self.setup.variables - self.equations.rank;
        let boundaries = self.boundaries.rank;
        let ext = self.dimension();
        let gauge = power(model, q, self.setup.gauges)?;
        let stabilizer = power(model, q, self.hom)?;
        let classes = power(model, q, ext)?;
        let cocycle_size = power(model, q, cocycles)?;
        if product(cocycle_size, stabilizer)? != product(classes, gauge)? {
            return Err(Error::plain(
                "partition extension groupoid failed the gauge orbit–stabilizer identity",
            ));
        }
        let arrows: Vec<_> = self
            .setup
            .arrows
            .iter()
            .enumerate()
            .filter_map(|(source, range)| {
                range
                    .as_ref()
                    .map(|r| (source, self.hall.successor[source].unwrap(), r.len()))
            })
            .collect();
        let splittings: Vec<_> = self
            .setup
            .splittings
            .iter()
            .map(std::ops::Range::len)
            .collect();
        let mass = fraction(model, classes, stabilizer)?;
        Ok(format!("rvPartitionExtensionGroupoid[{}] = q={q}, arrows={arrows:?}, vertexSplittings={splittings:?}, Z1={cocycles}, B1={boundaries}, Ext1={ext}, Hom={}, classes={}, stabilizer={}, mass={mass}",self.name,self.hom,show(classes),show(stabilizer)))
    }

    pub(crate) fn distribution_report(&self, model: &mut Model) -> Result<String> {
        if let Some(s) = self.distribution.borrow().as_ref() {
            return Ok(s.clone());
        }
        let q = self.hall.field.order();
        let count = self.cardinality();
        if count > CLASS_LIMIT {
            return Err(Error::plain(
                "partition extension class enumeration exceeds 50,000 representatives",
            ));
        }
        let mut actual = BTreeMap::<Shape, u128>::new();
        for j in 0..count {
            model.spend(1 + self.dimension())?;
            let class = self.class(digits(j, q, self.dimension()));
            let shape = self.middle_shape(model, &class)?;
            *actual.entry(shape).or_default() += 1;
        }
        let total: u128 = actual.values().sum();
        if from(total) != power(model, q, self.dimension())? {
            return Err(Error::plain(
                "partition extension class count did not equal q to the Ext dimension",
            ));
        }
        let stabilizer = power(model, q, self.hom)?;
        let end_aut = product(
            self.hall.aut(model, &self.quotient)?,
            self.hall.aut(model, &self.submodule)?,
        )?;
        let weight = product(stabilizer, end_aut)?;
        let dims = self
            .hall
            .dim(&self.quotient)
            .into_iter()
            .zip(self.hall.dim(&self.submodule))
            .map(|(a, b)| a + b)
            .collect();
        let candidates = self.hall.candidate_shapes(model, dims)?;
        let mut verified = Vec::new();
        for l in &candidates {
            let enumerated = actual.remove(l).unwrap_or(0);
            let flags = self
                .hall
                .coefficient(model, l, &self.quotient, &self.submodule)?;
            let auto = self.hall.aut(model, l)?;
            if product(from(enumerated), auto)? != product(from(flags), weight)? {
                return Err(Error::plain(
                    "partition Hall Riedtmann–Peng identity failed for a middle isotype",
                ));
            }
            if enumerated != 0 {
                verified.push(format!(
                    "({:?}, classes={enumerated}, hall={flags}, aut={}, mass={})",
                    l.0,
                    show(auto),
                    fraction(model, from(enumerated), stabilizer)?
                ));
            }
        }
        if !actual.is_empty() {
            return Err(Error::plain(
                "partition extension middle escaped the candidate isotypes",
            ));
        }
        let text=format!("rvPartitionExtensionDistribution[{}] = q={q}, total={total}, checkedIsotypes={}, middles=[{}]",self.name,candidates.len(),verified.join(", "));
        model.allocate(text.len())?;
        *self.distribution.borrow_mut() = Some(text.clone());
        Ok(text)
    }
}
