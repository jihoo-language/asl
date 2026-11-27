use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::model::Model;
use crate::rv_qp_quad::Wide;

use super::RANK;

// A Berkowitz term of degree ≤17 contains at most 17 matrix entries. The
// extra limbs allow exact sums of such terms even with full i128 entries.
const LIMBS: usize = 72;

#[derive(Clone, Copy)]
struct Integer(Wide<LIMBS>);

impl Integer {
    fn from_i128(value: i128) -> Self {
        Self(Wide::from_i128(value))
    }
}

impl PartialEq for Integer {
    fn eq(&self, other: &Self) -> bool {
        self.0.magnitude == other.0.magnitude
            && (self.0.negative == other.0.negative || self.0.is_zero())
    }
}

impl Eq for Integer {}

struct Integers;

impl CoefficientRing for Integers {
    type Element = Integer;

    fn cells(&self) -> usize {
        LIMBS + 1
    }

    fn zero(&self) -> Integer {
        Integer::from_i128(0)
    }

    fn one(&self) -> Integer {
        Integer::from_i128(1)
    }

    fn is_zero(&self, value: &Integer) -> bool {
        value.0.is_zero()
    }

    fn is_unit(&self, value: &Integer) -> bool {
        matches!(value.0.signed(), Some(1 | -1))
    }

    fn add(&self, model: &mut Model, a: &Integer, b: &Integer) -> Result<Integer> {
        model.spend(LIMBS)?;
        a.0.add(b.0).map(Integer).map_err(|_| arithmetic_range())
    }

    fn subtract(&self, model: &mut Model, a: &Integer, b: &Integer) -> Result<Integer> {
        self.add(model, a, &Integer(b.0.negate()))
    }

    fn multiply(&self, model: &mut Model, a: &Integer, b: &Integer) -> Result<Integer> {
        model
            .spend(LIMBS + a.0.magnitude.bits().div_ceil(32) * b.0.magnitude.bits().div_ceil(32))?;
        a.0.multiply_wide(b.0)
            .map(Integer)
            .map_err(|_| arithmetic_range())
    }

    fn inverse_unit(&self, _model: &mut Model, value: &Integer) -> Result<Integer> {
        if self.is_unit(value) {
            Ok(*value)
        } else {
            Err(Error::plain(
                "K3 characteristic calculation expected an integral unit",
            ))
        }
    }
}

fn arithmetic_range() -> Error {
    Error::plain("K3 lattice spectrum exceeds the exact arithmetic range")
}

fn characteristic(model: &mut Model, matrix: [[i128; RANK]; RANK]) -> Result<Vec<Integer>> {
    hw_algebra::storage(RANK * RANK, Integers.cells())?;
    model.spend(RANK * RANK * Integers.cells())?;
    let entries = matrix
        .into_iter()
        .flatten()
        .map(Integer::from_i128)
        .collect();
    let matrix = Matrix::from_entries(RANK, entries)?;
    hw_algebra::characteristic(model, &Integers, &matrix)
}

pub(super) fn characteristic_coefficients(
    model: &mut Model,
    matrix: [[i128; RANK]; RANK],
) -> Result<Vec<i128>> {
    coefficients(&characteristic(model, matrix)?)
}

pub(super) fn power_trace(
    model: &mut Model,
    matrix: [[i128; RANK]; RANK],
    exponent: u32,
) -> Result<Wide<LIMBS>> {
    // The coefficients are exact before computing power sums. Newton's
    // identities here use only integer multiplication and addition.
    let trace = if exponent == 1 {
        matrix
            .iter()
            .enumerate()
            .try_fold(Integers.zero(), |trace, (i, row)| {
                Integers.add(model, &trace, &Integer::from_i128(row[i]))
            })?
    } else {
        let polynomial = characteristic(model, matrix)?;
        let exponent = exponent as usize;
        let mut traces = Vec::with_capacity(exponent + 1);
        traces.push(Integers.zero());
        for power in 1..=exponent {
            let mut trace = Integers.zero();
            for (i, coefficient) in polynomial.iter().enumerate().skip(1).take(power) {
                let factor = if i == power {
                    Integer::from_i128(power as i128)
                } else {
                    traces[power - i]
                };
                let term = Integers.multiply(model, coefficient, &factor)?;
                trace = Integers.subtract(model, &trace, &term)?;
            }
            traces.push(trace);
        }
        traces[exponent]
    };
    Ok(trace.0)
}

pub(super) fn entropy(
    model: &mut Model,
    name: &str,
    matrix: [[i128; RANK]; RANK],
) -> Result<String> {
    let polynomial = characteristic(model, matrix)?;
    let characteristic = coefficients(&polynomial)?;
    let (cyclotomic, residual) = strip_cyclotomic(model, polynomial)?;
    let description = if residual == [Integers.one()] {
        "dynamical degree=1".into()
    } else {
        // An integral, positive-cone-preserving isometry of signature (1,n)
        // has at most one reciprocal pair away from the unit circle. Every
        // irreducible factor supported on the unit circle is cyclotomic.
        let at_one = residual.iter().try_fold(Integers.zero(), |sum, entry| {
            Integers.add(model, &sum, entry)
        })?;
        if residual.len() < 3
            || residual.len() % 2 == 0
            || !residual.iter().eq(residual.iter().rev())
            || at_one.0.is_zero()
            || !at_one.0.negative
        {
            return Err(Error::plain(
                "K3 automorphism residual failed the reciprocal and real-root checks",
            ));
        }
        let kind = if residual.len() == 3 {
            "reciprocal quadratic"
        } else {
            "Salem factor"
        };
        let residual = coefficients(&residual)?;
        format!("{kind}={residual:?}, dynamical degree=largest real root greater than 1 of {residual:?}")
    };
    Ok(format!("rvK3Entropy[{name}] = characteristic={characteristic:?}, cyclotomic={cyclotomic:?}, {description}"))
}

fn coefficients(polynomial: &[Integer]) -> Result<Vec<i128>> {
    polynomial
        .iter()
        .map(|entry| {
            entry.0.signed().ok_or_else(|| {
                Error::plain("K3 characteristic coefficient exceeds signed 128-bit range")
            })
        })
        .collect()
}

fn strip_cyclotomic(
    model: &mut Model,
    mut polynomial: Vec<Integer>,
) -> Result<(Vec<(usize, usize)>, Vec<Integer>)> {
    // In φ(m)^2/m = ∏ p^(a−2)(p−1)^2 only the factor for 2 can be
    // smaller than 1, by at most 1/2. Thus φ(m)≤17 implies m≤578.
    const MAX_ORDER: usize = 2 * RANK * RANK;
    model.spend(MAX_ORDER)?;
    let mut totient: Vec<_> = (0..=MAX_ORDER).collect();
    for p in 2..=MAX_ORDER {
        if totient[p] == p {
            model.spend(MAX_ORDER / p)?;
            for multiple in (p..=MAX_ORDER).step_by(p) {
                totient[multiple] -= totient[multiple] / p;
            }
        }
    }
    let mut factors = Vec::new();
    for (order, degree) in totient.into_iter().enumerate().skip(1) {
        if polynomial.len() == 1 {
            break;
        }
        if degree >= polynomial.len() {
            continue;
        }
        let modulus = model.cyclotomic(order)?;
        let divisor: Vec<_> = modulus
            .iter()
            .rev()
            .copied()
            .map(Integer::from_i128)
            .collect();
        let mut multiplicity = 0;
        while divisor.len() <= polynomial.len() {
            let (quotient, remainder) =
                hw_algebra::divide_monic(model, &Integers, &polynomial, &divisor)?;
            if remainder.iter().any(|entry| !Integers.is_zero(entry)) {
                break;
            }
            polynomial = quotient;
            multiplicity += 1;
        }
        if multiplicity != 0 {
            factors.push((order, multiplicity));
        }
    }
    Ok((factors, polynomial))
}
