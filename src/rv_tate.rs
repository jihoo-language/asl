use std::rc::Rc;

use crate::ec_poly::{self as poly, Field};
use crate::error::{Error, Result};
use crate::finite_field_factor::factor_with_multiplicity;
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;
use crate::rv_ore::gcd;
use crate::rv_unram::{add_mod, mul_mod, reduce, sub_mod, Ring};
use crate::syntax::ParsedExpr;

const MAX_DEGREE: usize = 256;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub(crate) struct Precision {
    pub(crate) prime: Prime,
    pub(crate) digits: usize,
    pub(crate) modulus: u64,
}

impl Precision {
    pub(crate) fn new(model: &mut Model, p: usize, digits: usize) -> Result<Self> {
        let prime = match p {
            2 => Prime::validated(2),
            3.. if !p.is_multiple_of(2) => Prime::new(model, p)?,
            _ => {
                return Err(Error::plain(
                    "Tate coefficient characteristic must be prime",
                ))
            }
        };
        if digits == 0 {
            return Err(Error::plain("Tate coefficient precision must be positive"));
        }
        let modulus = u32::try_from(digits)
            .ok()
            .and_then(|digits| u64::from(prime.get()).checked_pow(digits))
            .ok_or_else(|| Error::plain("Tate coefficient modulus p^s exceeds u64"))?;
        Ok(Self {
            prime,
            digits,
            modulus,
        })
    }

    fn p(self) -> u64 {
        u64::from(self.prime.get())
    }
}

pub(crate) enum Expression {
    Polynomial(usize, usize, Vec<i128>),
    Multiply(Box<[ParsedExpr; 2]>),
}

pub(crate) enum Report {
    Coefficient(usize),
    Prepare,
    Factors,
    Newton,
    Unramified(usize),
}

impl Report {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            Self::Coefficient(_) => "rvTateCoeff",
            Self::Prepare => "rvTatePrepare",
            Self::Factors => "rvTateFactors",
            Self::Newton => "rvTateNewton",
            Self::Unramified(_) => "rvTateUnramified",
        }
    }

    pub(crate) fn write(
        self,
        model: &mut Model,
        ring: Precision,
        coefficients: &[u64],
        coordinate: &str,
    ) -> Result<String> {
        let name = self.name();
        if let Self::Coefficient(j) = self {
            let entry = coefficients.get(j).copied().unwrap_or(0);
            return Ok(format!(
                "{name}[{j}] = coordinate={coordinate}, precision={}, modulus={}, residue={entry}",
                ring.digits, ring.modulus
            ));
        }
        let Some(prepared) = Preparation::new(model, ring, coefficients)? else {
            return Ok(format!(
                "{name} = coordinate={coordinate}, precision={}, unresolved",
                ring.digits
            ));
        };
        let header = format!(
            "{name} = coordinate={coordinate}, mu={}, precision={}, modulus={}",
            prepared.valuation, prepared.ring.digits, prepared.ring.modulus
        );
        match self {
            Self::Prepare => Ok(format!(
                "{header}, degree={}, polynomial={:?}, unit={:?}",
                prepared.polynomial.len() - 1,
                prepared.polynomial,
                prepared.unit
            )),
            Self::Factors => {
                let factors = prepared.factors(model)?;
                let triples: Vec<_> = factors
                    .iter()
                    .map(|f| (&f.residual, f.multiplicity, &f.lift))
                    .collect();
                Ok(format!("{header}, factors={triples:?}"))
            }
            Self::Newton => {
                model.spend(prepared.polynomial.len().pow(2))?;
                let (slopes, remaining) = newton_polygon(prepared.ring, &prepared.polynomial);
                Ok(format!(
                    "{header}, valuations={slopes:?}, unresolved={remaining:?}"
                ))
            }
            Self::Unramified(index) => {
                let factors = prepared.factors(model)?;
                let factor = factors.get(index).ok_or_else(|| {
                    Error::plain(format!(
                        "Weierstrass primary factor index {index} is out of range"
                    ))
                })?;
                if factor.multiplicity > 1 {
                    return Ok(format!(
                        "{header}, factor={index}, residual={:?}, multiplicity={}, unresolved",
                        factor.residual, factor.multiplicity
                    ));
                }
                let roots = prepared.unramified_roots(model, factor)?;
                Ok(format!(
                    "{header}, factor={index}, defining_polynomial={:?}, roots={roots:?}",
                    factor.residual
                ))
            }
            Self::Coefficient(_) => unreachable!(),
        }
    }
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Series(Vec<u64>);

impl Series {
    pub(crate) fn literal(ring: Precision, coefficients: Vec<i128>) -> Result<Self> {
        check_degree(coefficients.len())?;
        Ok(Self(trim(
            coefficients
                .into_iter()
                .map(|a| reduce(a, ring.modulus))
                .collect(),
        )))
    }

    pub(crate) fn zero() -> Self {
        Self(Vec::new())
    }

    pub(crate) fn coefficients(&self) -> &[u64] {
        &self.0
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + self.0.len())?;
        Ok(Value::Tate(Rc::new(self)))
    }

    pub(crate) fn add(&self, model: &mut Model, ring: Precision, other: &Self) -> Result<Self> {
        model.spend(self.0.len().max(other.0.len()))?;
        let values = (0..self.0.len().max(other.0.len()))
            .map(|j| {
                add_mod(
                    self.0.get(j).copied().unwrap_or(0),
                    other.0.get(j).copied().unwrap_or(0),
                    ring.modulus,
                )
            })
            .collect();
        Ok(Self(trim(values)))
    }

    pub(crate) fn negate(&self, model: &mut Model, ring: Precision) -> Result<Self> {
        model.spend(self.0.len())?;
        Ok(Self(
            self.0
                .iter()
                .map(|&a| sub_mod(0, a, ring.modulus))
                .collect(),
        ))
    }

    pub(crate) fn multiply(
        &self,
        model: &mut Model,
        ring: Precision,
        other: &Self,
    ) -> Result<Self> {
        if self.0.is_empty() || other.0.is_empty() {
            return Ok(Self::zero());
        }
        check_degree(self.0.len() + other.0.len() - 1)?;
        Ok(Self(product(model, &self.0, &other.0, ring.modulus)?))
    }
}

struct Preparation {
    ring: Precision,
    valuation: usize,
    polynomial: Vec<u64>,
    unit: Vec<u64>,
}

impl Preparation {
    fn new(model: &mut Model, original: Precision, input: &[u64]) -> Result<Option<Self>> {
        let p = original.p();
        let valuation = input
            .iter()
            .map(|&a| order(a, original))
            .min()
            .unwrap_or(original.digits);
        if valuation == original.digits {
            return Ok(None);
        }
        let scale = p.pow(valuation as u32);
        let ring = Precision {
            prime: original.prime,
            digits: original.digits - valuation,
            modulus: original.modulus / scale,
        };
        let h = trim(input.iter().map(|a| a / scale).collect());
        let residue = trim(h.iter().map(|a| a % p).collect());
        let u = *residue
            .last()
            .expect("unit Tate norm has nonzero reduction");
        let inverse_u = u64::from(invert(model, ring.prime, u as u32)?);
        let initial: Vec<_> = residue.iter().map(|&a| mul_mod(a, inverse_u, p)).collect();
        let mut polynomial = initial.clone();
        let mut unit = vec![u];
        let mut step = p;
        while step < ring.modulus {
            let error = defect(model, &h, &unit, &polynomial, step, p)?;
            let (quotient, mut remainder) = divide_monic(model, error, &initial, p)?;
            for coefficient in &mut remainder {
                *coefficient = mul_mod(*coefficient, inverse_u, p);
            }
            correct(&mut unit, &quotient, step);
            correct(&mut polynomial, &remainder, step);
            step *= p;
        }
        if product(model, &unit, &polynomial, ring.modulus)? != h {
            return Err(inconsistent());
        }
        Ok(Some(Self {
            ring,
            valuation,
            polynomial,
            unit,
        }))
    }

    fn factors(&self, model: &mut Model) -> Result<Vec<Primary>> {
        let p = self.ring.p();
        let residue = Ring::new(model, p as usize, 1, vec![0, 1])?;
        let field = Field { ring: &residue };
        let target = self
            .polynomial
            .iter()
            .map(|a| residue.scalar(a % p))
            .collect::<Vec<_>>();
        let factors = factor_with_multiplicity(model, field, &target)?;
        let mut result = Vec::with_capacity(factors.len());
        let mut remaining = self.polynomial.clone();
        let count = factors.len();
        for (index, (_, multiplicity, factor)) in factors.into_iter().enumerate() {
            let irreducible: Vec<_> = factor.iter().map(|a| a[0]).collect();
            let lift = if index + 1 == count {
                remaining.clone()
            } else {
                let mut primary = vec![1];
                for _ in 0..multiplicity {
                    primary = product(model, &primary, &irreducible, p)?;
                }
                let (lift, complement) =
                    split_primary(model, self.ring, field, &remaining, primary)?;
                remaining = complement;
                lift
            };
            result.push(Primary {
                residual: irreducible,
                multiplicity,
                lift,
            });
        }
        Ok(result)
    }

    fn unramified_roots(&self, model: &mut Model, factor: &Primary) -> Result<Vec<Vec<u64>>> {
        let ring = Ring::new(
            model,
            self.ring.p() as usize,
            self.ring.digits,
            factor.residual.iter().map(|&a| i128::from(a)).collect(),
        )?;
        let degree = ring.degree();
        let mut root = if degree == 1 {
            ring.scalar(sub_mod(0, factor.residual[0], ring.modulus()))
        } else {
            ring.literal(model, vec![0, 1])?
        };
        let modulus = ring.modulus();
        let lifted: Vec<_> = factor.lift.iter().map(|&a| ring.scalar(a)).collect();
        let polynomial: Vec<_> = self.polynomial.iter().map(|&a| ring.scalar(a)).collect();
        let derivative: Vec<_> = factor
            .lift
            .iter()
            .enumerate()
            .skip(1)
            .map(|(j, &a)| ring.scalar(mul_mod(j as u64, a, modulus)))
            .collect();
        let mut digits = 1;
        while digits < self.ring.digits {
            let value = ring.evaluate_polynomial(model, &lifted, &root)?;
            let slope = ring.evaluate_polynomial(model, &derivative, &root)?;
            let inverse = ring.inverse(model, &slope)?;
            let correction = ring.multiply(model, &inverse, &value)?;
            let negative = ring.negate(model, &correction)?;
            root = ring.add(model, &root, &negative)?;
            digits *= 2;
        }
        let mut conjugates = Vec::with_capacity(degree);
        for k in 0..degree {
            let conjugate = ring.frobenius(model, &root, k)?;
            let value = ring.evaluate_polynomial(model, &polynomial, &conjugate)?;
            if !poly::zero(&value) {
                return Err(inconsistent());
            }
            conjugates.push(conjugate);
        }
        Ok(conjugates)
    }
}

struct Primary {
    residual: Vec<u64>,
    multiplicity: usize,
    lift: Vec<u64>,
}

fn split_primary(
    model: &mut Model,
    ring: Precision,
    field: Field<'_>,
    target: &[u64],
    residual: Vec<u64>,
) -> Result<(Vec<u64>, Vec<u64>)> {
    let p = ring.p();
    let reduced: Vec<_> = target.iter().map(|&a| a % p).collect();
    let (complement, remainder) = divide_monic(model, reduced, &residual, p)?;
    if !remainder.is_empty() {
        return Err(inconsistent());
    }
    let field_residual = residual
        .iter()
        .map(|&a| field.ring.scalar(a))
        .collect::<Vec<_>>();
    let field_complement = complement
        .iter()
        .map(|&a| field.ring.scalar(a))
        .collect::<Vec<_>>();
    let inverse = poly::inverse_mod(model, field, &field_complement, &field_residual)?;
    let inverse = inverse.into_iter().map(|a| a[0]).collect::<Vec<_>>();
    let (mut left, mut right) = (residual.clone(), complement.clone());
    let mut step = p;
    while step < ring.modulus {
        let error = defect(model, target, &left, &right, step, p)?;
        let adjusted = product(model, &error, &inverse, p)?;
        let (_, d_left) = divide_monic(model, adjusted, &residual, p)?;
        let removed = product(model, &d_left, &complement, p)?;
        let difference = (0..error.len().max(removed.len()))
            .map(|j| {
                sub_mod(
                    error.get(j).copied().unwrap_or(0),
                    removed.get(j).copied().unwrap_or(0),
                    p,
                )
            })
            .collect();
        let (d_right, remainder) = divide_monic(model, difference, &residual, p)?;
        if !remainder.is_empty() || d_left.len() >= left.len() || d_right.len() >= right.len() {
            return Err(inconsistent());
        }
        correct(&mut left, &d_left, step);
        correct(&mut right, &d_right, step);
        step *= p;
    }
    if product(model, &left, &right, ring.modulus)? != target {
        return Err(inconsistent());
    }
    Ok((left, right))
}

type Slope = (usize, usize, usize);

fn newton_polygon(ring: Precision, coefficients: &[u64]) -> (Vec<Slope>, Option<Slope>) {
    let points: Vec<_> = coefficients.iter().map(|&a| order(a, ring)).collect();
    let mut end = points.len() - 1;
    let mut slopes = Vec::new();
    while end != 0 {
        let mut known = None;
        let mut unknown = None;
        for (j, &v) in points[..end].iter().enumerate() {
            let slot = if v == ring.digits {
                &mut unknown
            } else {
                &mut known
            };
            let slope = (v - points[end], end - j);
            if slot.is_none_or(|old| less(slope, old)) {
                *slot = Some(slope);
            }
        }
        // Starting from the monic right vertex, a strict comparison to the
        // valuation floors proves that the next entire edge is fixed.
        match (known, unknown) {
            (Some((rise, span)), bound) if bound.is_none_or(|b| less((rise, span), b)) => {
                slopes.push(reduced_slope(rise, span, span));
                end -= span;
            }
            (_, Some((rise, span))) => return (slopes, Some(reduced_slope(rise, span, end))),
            _ => unreachable!("ordered Newton vertices always have a left candidate"),
        }
    }
    (slopes, None)
}

fn reduced_slope(n: usize, d: usize, multiplicity: usize) -> Slope {
    let divisor = gcd(n as u128, d as u128) as usize;
    (n / divisor, d / divisor, multiplicity)
}

fn less((n, d): (usize, usize), (a, b): (usize, usize)) -> bool {
    n * b < a * d
}

fn order(mut value: u64, ring: Precision) -> usize {
    if value == 0 {
        return ring.digits;
    }
    let mut valuation = 0;
    while value.is_multiple_of(ring.p()) {
        value /= ring.p();
        valuation += 1;
    }
    valuation
}

fn defect(
    model: &mut Model,
    target: &[u64],
    left: &[u64],
    right: &[u64],
    step: u64,
    p: u64,
) -> Result<Vec<u64>> {
    let modulus = step * p;
    let current = product(model, left, right, modulus)?;
    (0..target.len().max(current.len()))
        .map(|j| {
            let a = target.get(j).copied().unwrap_or(0) % modulus;
            let b = current.get(j).copied().unwrap_or(0);
            let delta = sub_mod(a, b, modulus);
            if !delta.is_multiple_of(step) {
                return Err(inconsistent());
            }
            Ok(delta / step)
        })
        .collect()
}

fn correct(values: &mut Vec<u64>, correction: &[u64], step: u64) {
    values.resize(values.len().max(correction.len()), 0);
    for (coefficient, digit) in values.iter_mut().zip(correction) {
        *coefficient += step * digit;
    }
}

fn divide_monic(
    model: &mut Model,
    values: Vec<u64>,
    divisor: &[u64],
    p: u64,
) -> Result<(Vec<u64>, Vec<u64>)> {
    let mut remaining = trim(values);
    let count = remaining.len().saturating_sub(divisor.len() - 1);
    model.spend(count.saturating_mul(divisor.len()))?;
    let mut quotient = vec![0; count];
    for j in (0..count).rev() {
        let leading = remaining.pop().expect("coefficient of current degree");
        quotient[j] = leading;
        for (k, &coefficient) in divisor[..divisor.len() - 1].iter().enumerate() {
            remaining[j + k] = sub_mod(remaining[j + k], mul_mod(leading, coefficient, p), p);
        }
    }
    Ok((trim(quotient), trim(remaining)))
}

fn product(model: &mut Model, a: &[u64], b: &[u64], modulus: u64) -> Result<Vec<u64>> {
    if a.is_empty() || b.is_empty() {
        return Ok(Vec::new());
    }
    model.spend(a.len().saturating_mul(b.len()))?;
    let mut output = vec![0; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            output[i + j] = add_mod(output[i + j], mul_mod(x, y, modulus), modulus);
        }
    }
    Ok(trim(output))
}

fn trim(mut values: Vec<u64>) -> Vec<u64> {
    while values.last() == Some(&0) {
        values.pop();
    }
    values
}

fn check_degree(length: usize) -> Result<()> {
    if length > MAX_DEGREE + 1 {
        return Err(Error::plain(format!(
            "Tate polynomial degree exceeds {MAX_DEGREE}"
        )));
    }
    Ok(())
}

fn inconsistent() -> Error {
    Error::plain("inconsistent Weierstrass lifting identity")
}
