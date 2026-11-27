use crate::dwork_k3::fricke::{ENDPOINT_PRODUCT, ENDPOINT_SUM};
use crate::ec_poly::{self as poly, Field, Polynomial as FieldPoly};
use crate::error::{Error, Result};
use crate::finite_field_factor::factor_with_multiplicity;
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;
use crate::rv_unram::{self, Ring};

use super::super::{residue, Poly};
use super::{certificate, Kind, Polynomial};

pub(super) fn eval(
    model: &mut Model,
    ring: &Ring,
    polynomial: &Poly,
    x: &[u64],
) -> Result<Vec<u64>> {
    let mut value = ring.scalar(0);
    for &coefficient in polynomial.0.iter().rev() {
        value = ring.multiply(model, &value, x)?;
        let c = residue(model, coefficient, ring.modulus())?;
        value[0] = rv_unram::add_mod(value[0], c, ring.modulus());
    }
    Ok(value)
}

fn field_rows(model: &mut Model, ring: &Ring, polynomial: &Poly) -> Result<FieldPoly> {
    polynomial
        .0
        .iter()
        .map(|&c| Ok(ring.scalar(residue(model, c, ring.modulus())?)))
        .collect()
}

fn endpoint_polynomial(model: &mut Model, ring: &Ring, h: &[u64]) -> Result<FieldPoly> {
    let evaluate = |model: &mut Model, coefficients: &[i128]| -> Result<Vec<u64>> {
        let coefficients = coefficients
            .iter()
            .map(|&c| ring.scalar(rv_unram::reduce(c, ring.modulus())))
            .collect::<Vec<_>>();
        ring.evaluate_polynomial(model, &coefficients, h)
    };
    let sum = evaluate(model, &ENDPOINT_SUM)?;
    let product = evaluate(model, &ENDPOINT_PRODUCT)?;
    let middle = ring.negate(model, &sum)?;
    Ok(vec![product, middle, ring.scalar(1)])
}

pub(crate) fn pair_coefficient(
    model: &mut Model,
    ring: &Ring,
    h: &[u64],
    i: usize,
) -> Result<Vec<u64>> {
    Ok(endpoint_polynomial(model, ring, h)?
        .get(i)
        .cloned()
        .unwrap_or_else(|| ring.scalar(0)))
}

type Factors = Vec<(usize, usize, FieldPoly)>;

fn factors(model: &mut Model, f: Field<'_>, p: &FieldPoly) -> Result<Factors> {
    let factors = factor_with_multiplicity(model, f, p)?;
    let mut reconstructed = vec![f.one()];
    for (_, multiplicity, factor) in &factors {
        for _ in 0..*multiplicity {
            reconstructed = poly::multiply(model, f, &reconstructed, factor)?;
        }
    }
    if *p != reconstructed {
        return Err(certificate());
    }
    Ok(factors)
}

fn field(ring: &Ring) -> Result<Field<'_>> {
    if ring.precision() != 1 {
        return Err(Error::plain(
            "Fricke CM factorization requires a residue field at precision one",
        ));
    }
    Ok(Field { ring })
}

fn describe_pair(model: &mut Model, f: Field<'_>, p: &FieldPoly) -> Result<String> {
    let factors = factors(model, f, p)?;
    let status = if factors.len() == 1 && factors[0].0 == 2 {
        let mut conjugate = vec![f.zero(), f.one()];
        for _ in 0..f.ring.degree() {
            conjugate = poly::power_mod(model, f, &conjugate, f.ring.prime(), p)?;
        }
        let sum = f.ring.negate(model, &p[1])?;
        if conjugate != vec![sum, f.signed(-1)] {
            return Err(certificate());
        }
        "nonsplit, Frobenius exchanges the roots"
    } else if factors.len() == 1 && factors[0].1 == 2 {
        "repeated endpoint"
    } else {
        "two distinct split endpoints"
    };
    let characteristic = if f.ring.prime() == 2 {
        ", characteristic-two polynomial reduction"
    } else {
        ""
    };
    Ok(format!("endpoint polynomial={p:?}, factors (degree, multiplicity, coefficients)={factors:?}, {status}{characteristic}"))
}

pub(crate) fn endpoint_report(model: &mut Model, ring: &Ring, h: &[u64]) -> Result<String> {
    let f = field(ring)?;
    let p = endpoint_polynomial(model, ring, h)?;
    Ok(format!(
        "rvFrickeCMEndpoints = H={h:?}, {}",
        describe_pair(model, f, &p)?
    ))
}

pub(super) fn factors_report(
    model: &mut Model,
    value: &Polynomial,
    name: &str,
    ring: &Ring,
) -> Result<String> {
    let f = field(ring)?;
    let main = field_rows(model, f.ring, &value.support)?;
    let main = factors(model, f, &main)?;
    let mut directed = Vec::new();
    match &value.kind {
        Kind::Norm { strata, .. } => {
            for (m, factor) in strata {
                let reduced = field_rows(model, f.ring, factor)?;
                let reduced = factors(model, f, &reduced)?;
                directed.push(format!("({m}, {reduced:?})"));
            }
        }
        Kind::Edge { strata, .. } => {
            for (a, b, factor) in strata {
                let reduced = field_rows(model, f.ring, factor)?;
                let reduced = factors(model, f, &reduced)?;
                directed.push(format!("({a}, {b}, {reduced:?})"));
            }
        }
    }
    let characteristic = if ring.prime() == 2 {
        ", characteristic-two polynomial reduction"
    } else {
        ""
    };
    Ok(format!(
        "rvFrickeCMPolyFactors[{name}] = factors (degree, multiplicity, coefficients)={main:?}, reduced directed strata=[{}]{characteristic}",
        directed.join(", ")
    ))
}

pub(super) fn edge_endpoints(
    model: &mut Model,
    value: &Polynomial,
    name: &str,
    ring: &Ring,
    h: &[u64],
) -> Result<String> {
    let f = field(ring)?;
    let Kind::Edge { sources, .. } = &value.kind else {
        return Err(Error::plain(
            "Fricke CM edge endpoint recovery requires a declared edge polynomial",
        ));
    };
    if !poly::zero(&eval(model, f.ring, &value.support, h)?) {
        return Err(Error::plain(
            "H is outside the reduced characteristic-zero Fricke edge support",
        ));
    }
    let p = endpoint_polynomial(model, ring, h)?;
    let mut compatibilities = Vec::new();
    for source in sources {
        let class = source.unramified_coefficients(model, f.ring)?;
        let common = poly::monic_gcd(model, f, p.clone(), class)?;
        if common.len() < 2 {
            return Err(certificate());
        }
        compatibilities.push(common);
    }
    let overlap = poly::monic_gcd(
        model,
        f,
        compatibilities[0].clone(),
        compatibilities[1].clone(),
    )?;
    let order_overlap = if sources[0].discriminant != sources[1].discriminant && overlap.len() > 1 {
        ", residue order labels overlap"
    } else {
        ""
    };
    Ok(format!(
        "rvFrickeCMEdgeEndpoints[{name}] = H={h:?}, {}, order-compatible factors ({}, {})={compatibilities:?}{order_overlap}",
        describe_pair(model, f, &p)?, sources[0].discriminant, sources[1].discriminant
    ))
}
