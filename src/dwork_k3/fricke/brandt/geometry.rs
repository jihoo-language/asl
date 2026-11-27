use std::cell::OnceCell;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use super::{certificate, Basis};
use crate::dwork_k3::fricke::{ENDPOINT_PRODUCT, ENDPOINT_SUM};
use crate::ec_poly::{self as poly, Field, Scalar};
use crate::error::{Error, Result};
use crate::finite_field_factor;
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear;
use crate::rv_unram::{self, Ring};

pub(super) const MAX_PRIME: u64 = 101;

impl Basis {
    pub(super) fn geometry(model: &mut Model, field_name: Rc<str>, ring: Rc<Ring>) -> Result<Self> {
        let p = ring.prime();
        if ring.degree() != 2 || ring.precision() != 1 || !(5..=MAX_PRIME).contains(&p) {
            return Err(Error::plain(format!(
                "supersingular Fricke geometry requires a named quadratic residue field with 5 ≤ p ≤ {MAX_PRIME}"
            )));
        }
        let field = Field { ring: &ring };
        let js = supersingular_js(model, field)?;
        let product = field.polynomial(&ENDPOINT_PRODUCT);
        let sum = field.polynomial(&ENDPOINT_SUM);
        let mut vertices = BTreeMap::new();
        for j in &js {
            let mut norm = product.clone();
            for (c, s) in norm.iter_mut().zip(&sum) {
                let term = field.multiply(model, j, s)?;
                *c = field.subtract(model, c, &term)?;
            }
            let j_squared = field.multiply(model, j, j)?;
            norm[0] = field.add(model, &norm[0], &j_squared)?;
            for (h, _) in split_roots(model, field, &norm)? {
                let endpoint_sum = poly::evaluate(model, field, &sum, &h)?;
                let other = field.subtract(model, &endpoint_sum, j)?;
                let endpoint_product = poly::evaluate(model, field, &product, &h)?;
                if !js.contains(&other) || field.multiply(model, j, &other)? != endpoint_product {
                    return Err(certificate("both supersingular endpoints"));
                }
                let pair = if *j <= other {
                    [j.clone(), other]
                } else {
                    [other, j.clone()]
                };
                if let Some(previous) = vertices.insert(h, pair.clone()) {
                    if previous != pair {
                        return Err(certificate("consistent unordered endpoints"));
                    }
                }
            }
        }
        let (vertices, endpoints): (Vec<_>, Vec<_>) = vertices.into_iter().unzip();
        let mut frobenius = Vec::with_capacity(vertices.len());
        for h in &vertices {
            let image = ring.frobenius(model, h, 1)?;
            frobenius.push(
                vertices
                    .binary_search(&image)
                    .map_err(|_| certificate("Frobenius preserves the Fricke divisor"))?,
            );
        }
        if frobenius
            .iter()
            .enumerate()
            .any(|(i, &j)| frobenius[j] != i)
        {
            return Err(certificate("quadratic Frobenius involution"));
        }
        let j_polynomial = prime_polynomial(model, field, js.iter())?;
        let h_polynomial = prime_polynomial(model, field, vertices.iter())?;
        let js = js.into_iter().collect::<Vec<_>>();
        let n = vertices.len();
        model.allocate(28 + 6 * js.len() + 15 * n)?;
        Ok(Basis {
            field_name,
            ring,
            operators: std::array::from_fn(|_| OnceCell::new()),
            p,
            js,
            j_polynomial,
            h_polynomial,
            vertices,
            endpoints,
            frobenius,
        })
    }
}

fn supersingular_js(model: &mut Model, field: Field<'_>) -> Result<BTreeSet<Scalar>> {
    let p = field.ring.prime() as u32;
    let prime = Prime::validated(p);
    let m = (p - 1) / 2;
    let mut binomial = 1u32;
    let mut deuring = vec![field.one()];
    for k in 1..=m {
        let inverse = prime_linear::invert(model, prime, k)?;
        binomial = binomial * (m + 1 - k) % p * inverse % p;
        deuring.push(field.ring.scalar(u64::from(binomial * binomial % p)));
    }
    let mut js = BTreeSet::new();
    for (lambda, multiplicity) in split_roots(model, field, &deuring)? {
        let one_minus = field.subtract(model, &field.one(), &lambda)?;
        if multiplicity != 1 || field.is_zero(&lambda) || field.is_zero(&one_minus) {
            return Err(certificate(
                "simple Legendre supersingular roots away from 0 and 1",
            ));
        }
        let lambda_squared = field.multiply(model, &lambda, &lambda)?;
        let invariant = field.add(model, &one_minus, &lambda_squared)?;
        let denominator = field.multiply(model, &lambda, &one_minus)?;
        let denominator = field.multiply(model, &denominator, &denominator)?;
        let inverse = field.inverse_unit(model, &denominator)?;
        let numerator = field.ring.power(model, &invariant, 3)?;
        let numerator = field.ring.scale(model, &numerator, u64::from(256 % p))?;
        js.insert(field.multiply(model, &numerator, &inverse)?);
    }
    let correction = match p % 12 {
        1 => 0,
        5 | 7 => 1,
        11 => 2,
        _ => return Err(certificate("prime characteristic ≥ 5")),
    };
    if js.len() != (p / 12 + correction) as usize
        || js.contains(&field.zero()) != (p % 3 == 2)
        || js.contains(&field.signed(1728)) != (p % 4 == 3)
    {
        return Err(certificate(
            "Deuring's supersingular j count and exceptional j values",
        ));
    }
    Ok(js)
}

fn split_roots(
    model: &mut Model,
    field: Field<'_>,
    polynomial: &[Scalar],
) -> Result<Vec<(Scalar, usize)>> {
    let factors = finite_field_factor::factor_with_multiplicity(model, field, polynomial)?;
    let mut roots = Vec::with_capacity(factors.len());
    let mut degree = 0;
    for (d, multiplicity, linear) in factors {
        if d != 1 {
            return Err(certificate(
                "complete splitting over the quadratic residue field",
            ));
        }
        degree += multiplicity;
        roots.push((
            field.subtract(model, &field.zero(), &linear[0])?,
            multiplicity,
        ));
    }
    if degree + 1 != polynomial.len() {
        return Err(certificate("complete factor multiplicities"));
    }
    Ok(roots)
}

fn prime_polynomial<'a>(
    model: &mut Model,
    field: Field<'_>,
    roots: impl IntoIterator<Item = &'a Scalar>,
) -> Result<Vec<u64>> {
    let mut polynomial = vec![field.one()];
    for root in roots {
        let constant = field.subtract(model, &field.zero(), root)?;
        polynomial = poly::multiply(model, field, &polynomial, &[constant, field.one()])?;
    }
    polynomial
        .iter()
        .map(|c| rv_unram::require_scalar(c))
        .collect()
}
