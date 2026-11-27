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

pub(crate) mod composition;

pub(crate) fn degree(level: usize) -> Result<usize> {
    if level == 9 {
        return Ok(12);
    }
    Ok(elliptic_modular::degree(prime_level(level)?))
}

pub(crate) fn prime_level(ell: usize) -> Result<Prime> {
    if matches!(ell, 3 | 5) {
        Ok(Prime::validated(ell as u32))
    } else {
        Err(Error::plain(
            "Dwork Fricke prime Hecke equations support good levels 3 and 5",
        ))
    }
}

pub(crate) fn construct(model: &mut Model, level: usize) -> Result<Polynomial> {
    let d = degree(level)?;
    let prime = prime_level(if level == 9 { 3 } else { level })?;
    let polynomial = if level == 9 {
        Polynomial::from_cyclic_square(model, prime, AUDIT, shifted_hauptmodul)?
    } else {
        Polynomial::from_hauptmodul(model, prime, AUDIT, shifted_hauptmodul)?
    };
    let ell = prime.get();
    let mut expected = vec![vec![0i32; d + 1]; d + 1];
    let k = if level == 9 { 9 } else { ell as usize };
    let kronecker = [(0, k + 1, 1), (k + 1, 0, 1), (1, 1, -1), (k, k, -1)];
    if level == 9 {
        // (Y-X^9)(Y^9-X)(Y-X)^2 in the level characteristic.
        for (i, j, coefficient) in kronecker {
            for (x, y, c) in [(2, 0, 1), (1, 1, -2), (0, 2, 1)] {
                expected[i + x][j + y] += coefficient * c;
            }
        }
    } else {
        for (i, j, c) in kronecker {
            expected[i][j] = c;
        }
    }
    for (i, row) in expected.iter().enumerate() {
        for (j, &c) in row.iter().enumerate() {
            let coefficient = polynomial.wide_coefficient::<40>(i, j)?;
            let residue = coefficient.magnitude.residue(ell);
            let residue = if coefficient.negative {
                prime.negate(residue)
            } else {
                residue
            };
            if residue != c.rem_euclid(ell as i32) as u32 {
                return Err(Error::plain(
                    "Dwork Fricke Hecke equation failed the Kronecker congruence",
                ));
            }
        }
    }
    if level == 9 {
        composition::certify_nine(model, &polynomial)?;
    }
    Ok(polynomial)
}

pub(crate) fn equation_report(model: &mut Model, ell: usize) -> Result<String> {
    let d = degree(ell)?;
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
    let resultant = if ell == 9 {
        ", certified Res(Phi_3,Phi_3)=(X-Y)^4 Phi_9 over Z"
    } else {
        ""
    };
    Ok(format!(
        "rvDworkHeckeEquation[{ell}] = bidegree=({d}, {d}), terms={} [{}], certified symmetry, q coefficients through {AUDIT}, and Kronecker congruence{resultant}",
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

fn eta_pair(model: &mut Model, cutoff: usize) -> Result<(Vec<Integer>, Vec<Integer>)> {
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
    Ok((t, inv))
}

fn shifted_hauptmodul(model: &mut Model, cutoff: usize) -> Result<Vec<Integer>> {
    let (mut t, inv) = eta_pair(model, cutoff)?;
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

pub(crate) fn oriented_equation(model: &mut Model, level: usize) -> Result<Polynomial> {
    let prime = prime_level(level)?;
    let polynomial = Polynomial::from_hauptmodul(model, prime, AUDIT, |model, cutoff| {
        let (t, _) = eta_pair(model, cutoff)?;
        Ok(t)
    })?;
    let ell = level;
    // The Kronecker congruence uses the same determinant-ell
    // correspondences on the ordered eta Hauptmodul.
    for j in 0..=ell + 1 {
        for i in 0..=ell + 1 {
            let coefficient = polynomial.wide_coefficient::<20>(i, j)?;
            let residue = coefficient.magnitude.residue(prime.get());
            let residue = if coefficient.negative {
                prime.negate(residue)
            } else {
                residue
            };
            let expected = match (i, j) {
                (a, b) if a == ell && b == ell || a == 1 && b == 1 => ell - 1,
                (a, 0) | (0, a) if a == ell + 1 => 1,
                _ => 0,
            };
            if residue as usize != expected {
                return Err(Error::plain(
                    "ordered Γ0(2) Hecke equation failed the Kronecker congruence",
                ));
            }
        }
    }
    Ok(polynomial)
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
    let d = degree(ell)?;
    let equation = model.dwork_hecke_equation(ell)?;
    let specialized = composition::specialize(model, &equation, ring, x, d)?;
    let data = factors(model, ring, &specialized)?;
    let characteristic = if ring.prime() == 3 && ell == 9 {
        ", cyclic level in characteristic three, schematic multiplicities"
    } else if ring.prime() == ell as u64 {
        ", characteristic equals level, Frobenius correspondence"
    } else {
        ""
    };
    Ok(format!(
        "rvDworkHeckeNeighbors[{ell}] = H={x:?}, {}{characteristic}",
        data.describe()
    ))
}

type Factor = (usize, usize, poly::Polynomial);

struct Factorization {
    factors: Vec<Factor>,
    roots: Vec<(poly::Scalar, usize)>,
}

impl Factorization {
    fn describe(&self) -> String {
        format!(
            "factors (degree, multiplicity, coefficients)={:?}, rational roots (value, multiplicity)={:?}",
            self.factors, self.roots
        )
    }
}

fn factors(model: &mut Model, ring: &Ring, polynomial: &[poly::Scalar]) -> Result<Factorization> {
    let f = Field { ring };
    let factors = factor_with_multiplicity(model, f, polynomial)?;
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
    if product != polynomial {
        return Err(Error::plain(
            "Dwork Fricke Hecke factors failed their product check",
        ));
    }
    roots.sort_unstable();
    Ok(Factorization { factors, roots })
}

fn require_field(ring: &Ring) -> Result<()> {
    if ring.prime() == 2 || ring.precision() != 1 {
        return Err(Error::plain(
            "Dwork Fricke Hecke neighbors require an odd residue field at precision one",
        ));
    }
    Ok(())
}
