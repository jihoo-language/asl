use crate::ec_poly::{self as poly, Field};
use crate::elliptic_modular::{self, Integer, Polynomial};
use crate::error::{Error, Result};
use crate::finite_field_factor::factor_with_multiplicity;
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_unram::Ring;

use super::eta::eisenstein_difference;

const AUDIT: usize = 8;

pub(crate) fn level(ell: usize) -> Result<Prime> {
    if matches!(ell, 3 | 5) {
        Ok(Prime::validated(ell as u32))
    } else {
        Err(Error::plain(
            "Dwork Fricke Hecke equations support good levels 3 and 5",
        ))
    }
}

pub(crate) fn construct(model: &mut Model, prime: Prime) -> Result<Polynomial> {
    let ell = prime.get();
    let polynomial = Polynomial::from_hauptmodul(model, prime, AUDIT, shifted_hauptmodul)?;
    // Check the integral polynomial against both Frobenius components
    // of its reduction at the isogeny characteristic.
    let d = ell as usize;
    for i in 0..=d + 1 {
        for j in 0..=d + 1 {
            let coefficient = polynomial.wide_coefficient::<40>(i, j)?;
            let residue = coefficient.magnitude.residue(ell);
            let residue = if coefficient.negative {
                prime.negate(residue)
            } else {
                residue
            };
            let expected = match (i, j) {
                (0, b) | (b, 0) if b == d + 1 => 1,
                (1, 1) => ell - 1,
                (a, b) if a == d && b == d => ell - 1,
                _ => 0,
            };
            if residue != expected {
                return Err(Error::plain(
                    "Dwork Fricke Hecke equation failed the Kronecker congruence",
                ));
            }
        }
    }
    Ok(polynomial)
}

pub(crate) fn equation_report(model: &mut Model, ell: usize) -> Result<String> {
    let d = elliptic_modular::degree(level(ell)?);
    let polynomial = model.dwork_hecke_equation(ell)?;
    let mut terms = Vec::new();
    for j in 0..=d {
        for i in 0..=d {
            let coefficient = polynomial.coefficient_text(i, j);
            if coefficient != "0" {
                terms.push(format!("({i}, {j}): {coefficient}"));
            }
        }
    }
    Ok(format!(
        "rvDworkHeckeEquation[{ell}] = bidegree=({d}, {d}), terms={} [{}], certified symmetry, q coefficients through {AUDIT}, and Kronecker congruence",
        terms.len(), terms.join(", ")
    ))
}

pub(crate) fn coefficient_report(
    model: &mut Model,
    ell: usize,
    i: usize,
    j: usize,
) -> Result<String> {
    let p = model.dwork_hecke_equation(ell)?;
    Ok(format!(
        "rvDworkHeckeCoeff[{ell}, {i}, {j}] = {}",
        p.coefficient_text(i, j)
    ))
}

fn shifted_hauptmodul(model: &mut Model, cutoff: usize) -> Result<Vec<Integer>> {
    let weights = eisenstein_difference(cutoff);
    model.spend(cutoff * cutoff)?;
    // For T=qt, theta(T)/T = -(2E2(2tau)-E2(tau)-1).
    // Integrate this equation for T and its reciprocal over Z.
    let mut t = vec![Integer::from_i128(1)];
    let mut inv = t.clone();
    for n in 1..=cutoff {
        let mut a = Integer::from_i128(0);
        let mut b = Integer::from_i128(0);
        for (i, weight) in weights.iter().enumerate().take(n + 1).skip(1) {
            a = a.add(t[n - i].multiply(-weight.numerator)?)?;
            b = b.add(inv[n - i].multiply(weight.numerator)?)?;
        }
        t.push(elliptic_modular::divide_small(a, n)?);
        inv.push(elliptic_modular::divide_small(b, n)?);
    }
    // qH=T+24q+4096q^2/T.
    if cutoff > 0 {
        t[1] = t[1].add(Integer::from_i128(24))?;
    }
    for (coefficient, value) in t.iter_mut().skip(2).zip(inv) {
        *coefficient = coefficient.add(value.multiply(4096)?)?;
    }
    if cutoff >= 3
        && !t
            .iter()
            .take(4)
            .zip([1, 0, 4372, 96256])
            .all(|(c, expected)| c.signed() == Some(expected))
    {
        return Err(Error::plain("Dwork Fricke Hauptmodul normalization failed"));
    }
    Ok(t)
}

#[derive(Clone, Copy)]
pub(crate) enum Coordinate {
    Z,
    T,
}

impl Coordinate {
    pub(crate) fn apply(self, model: &mut Model, ring: &Ring, x: &[u64]) -> Result<Vec<u64>> {
        require_field(ring)?;
        if poly::zero(x) {
            return Err(Error::plain(
                "Dwork Hauptmodul coordinate requires a nonzero input",
            ));
        }
        let field = Field { ring };
        let inverse = ring.inverse(model, x)?;
        match self {
            Self::Z => field.subtract(model, &inverse, &field.signed(104)),
            Self::T => {
                let partner = ring.scale(model, &inverse, 4096 % ring.prime())?;
                let direct = ring.add(model, x, &field.signed(24))?;
                ring.add(model, &direct, &partner)
            }
        }
    }
}

pub(crate) fn neighbors(model: &mut Model, ell: usize, ring: &Ring, x: &[u64]) -> Result<String> {
    require_field(ring)?;
    let d = elliptic_modular::degree(level(ell)?);
    let equation = model.dwork_hecke_equation(ell)?;
    let f = Field { ring };
    let specialized = (0..=d)
        .map(|j| equation.specialization(model, ring, x, j))
        .collect::<Result<Vec<_>>>()?;
    let factors = factor_with_multiplicity(model, f, &specialized)?;
    let mut roots = Vec::new();
    let mut product = vec![f.one()];
    for (degree, multiplicity, factor) in &factors {
        if *degree == 1 {
            roots.push((ring.negate(model, &factor[0])?, *multiplicity));
        }
        for _ in 0..*multiplicity {
            product = poly::multiply(model, f, &product, factor)?;
        }
    }
    if product != specialized {
        return Err(Error::plain(
            "Dwork Fricke Hecke factors failed their product check",
        ));
    }
    roots.sort_unstable();
    let characteristic = if ring.prime() == ell as u64 {
        ", characteristic equals level, Frobenius correspondence"
    } else {
        ""
    };
    Ok(format!(
        "rvDworkHeckeNeighbors[{ell}] = H={x:?}, factors (degree, multiplicity, coefficients)={factors:?}, rational roots (value, multiplicity)={roots:?}{characteristic}"
    ))
}

fn require_field(ring: &Ring) -> Result<()> {
    if ring.prime() == 2 || ring.precision() != 1 {
        return Err(Error::plain(
            "Dwork Fricke Hecke neighbors require an odd residue field at precision one",
        ));
    }
    Ok(())
}
