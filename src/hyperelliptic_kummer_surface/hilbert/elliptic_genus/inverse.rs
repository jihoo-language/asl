use std::ops::Bound::{Excluded, Included};
use std::rc::Rc;

use crate::error::Result;
use crate::model::Model;

use super::super::Magnitude;
use super::{
    decimal, eta_product, multiply, overflow, series, Cutoff, Expansion, Integer, Product,
};

#[derive(Default)]
pub(crate) struct Cache {
    rectangle: Option<(Cutoff, Rc<Product>)>,
    eta_inverse: Vec<Integer>,
}

#[derive(Clone, Copy)]
pub(crate) enum Chamber {
    Inside,
    Outside,
}

pub(crate) enum Query {
    Coefficient(i128, i128, i128, Chamber),
    WallJump(i128, i128, i128),
    Residue(i128, i128),
}

impl Query {
    pub(super) fn report(self, model: &mut Model) -> Result<String> {
        let (label, value) = match self {
            Self::Coefficient(m, n, r, chamber) => {
                let name = match chamber {
                    Chamber::Inside => "inside",
                    Chamber::Outside => "outside",
                };
                (
                    format!("rvIgusaInverse[{m}, {n}, {r}, {name}]"),
                    coefficient(model, m, n, Integer::from_i128(r), chamber)?,
                )
            }
            Self::WallJump(m, n, r) => (
                format!("rvIgusaWallJump[{m}, {n}, {r}]"),
                wall_jump(model, m, n, Integer::from_i128(r))?,
            ),
            Self::Residue(m, n) => (format!("rvIgusaResidue[{m}, {n}]"), residue(model, m, n)?),
        };
        Ok(format!("{label} = {}", decimal(value)))
    }
}

fn shifted_degree(index: i128) -> Option<u128> {
    match index {
        -1 => Some(0),
        n if n >= 0 => Some(n as u128 + 1),
        _ => None,
    }
}

pub(super) fn coefficient(
    model: &mut Model,
    m: i128,
    n: i128,
    r: Integer,
    chamber: Chamber,
) -> Result<Integer> {
    let mut total = Integer::from_i128(0);
    let Some((p, q)) = shifted_degree(m).zip(shifted_degree(n)) else {
        return Ok(total);
    };
    // Every root with a+b>0 satisfies |s| ≤ a+b, so a full row at
    // (p,q) has no Laurent exponent outside [−p−q,p+q].
    let on_empty_side = match chamber {
        Chamber::Inside => r.negative || r.is_zero(),
        Chamber::Outside => !r.negative || r.is_zero(),
    };
    if on_empty_side
        && !r
            .magnitude
            .cmp(Magnitude::from_small(p.saturating_add(q)))
            .is_lt()
    {
        return Ok(total);
    }
    let product = rectangle(model, p, q)?;
    // The cutoff or a containing cached rectangle checked both indices.
    let (p, q) = (p as usize, q as usize);
    let first = Included((p, q, i128::MIN));
    let last = Included((p, q, i128::MAX));
    let bounds = match (chamber, r.signed()) {
        (Chamber::Inside, Some(r)) => (first, Excluded((p, q, r))),
        (Chamber::Outside, Some(r)) => (Excluded((p, q, r)), last),
        // The degree cutoff bounds every Laurent exponent tightly. A
        // wider r lies beyond the whole row on the contributing side.
        (_, None) => (first, last),
    };
    for (&(_, _, s), &value) in product.range(bounds) {
        model.spend(1)?;
        let mut scale = r
            .add(Integer::from_i128(s).negate())
            .map_err(|_| overflow())?;
        scale.negative = false;
        total = total
            .add(multiply(model, value, scale)?)
            .map_err(|_| overflow())?;
    }
    Ok(total)
}

pub(super) fn wall_jump(model: &mut Model, m: i128, n: i128, r: Integer) -> Result<Integer> {
    if r.is_zero() {
        return Ok(Integer::from_i128(0));
    }
    let residue = residue(model, m, n)?;
    multiply(model, residue, r)
}

fn rectangle(model: &mut Model, p: u128, q: u128) -> Result<Rc<Product>> {
    if let Some((cutoff, product)) = &model.igusa_inverse_cache.rectangle {
        if p <= cutoff.p as u128 && q <= cutoff.q as u128 {
            return Ok(product.clone());
        }
    }
    let cutoff = Cutoff::new(model, p, q)?;
    let product = Rc::new(series(model, cutoff, Expansion::IgusaInverse)?);
    // One bounded rectangle warms every Laurent slice below its cutoff.
    model.igusa_inverse_cache.rectangle = Some((cutoff, product.clone()));
    Ok(product)
}

fn residue(model: &mut Model, m: i128, n: i128) -> Result<Integer> {
    let Some((p, q)) = shifted_degree(m).zip(shifted_degree(n)) else {
        return Ok(Integer::from_i128(0));
    };
    let bound = p.max(q);
    if bound >= model.igusa_inverse_cache.eta_inverse.len() as u128 {
        let cutoff = Cutoff::new(model, 0, bound)?;
        model.igusa_inverse_cache.eta_inverse = eta_product(model, cutoff.q, -24)?;
    }
    let eta = &model.igusa_inverse_cache.eta_inverse;
    let (first, second) = (eta[p as usize], eta[q as usize]);
    multiply(model, first, second)
}
