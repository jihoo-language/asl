use std::cmp::Ordering;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear;
use crate::rv_tate;
use crate::rv_unram::{add_mod, mul_mod, reduce};
use crate::syntax::ParsedExpr;

use super::super::{residual_hecke::lift_simple_root, Integer};
use super::amice::{self, Amice};
use super::mazur_tate::prime_power;

#[derive(Clone, Copy)]
pub(crate) struct Chart(pub(crate) i128, pub(crate) usize);

pub(crate) struct Spec(pub(crate) amice::Spec, pub(crate) Chart);

pub(crate) enum Expression {
    Disc(Chart, Box<ParsedExpr>),
    Restrict(Chart, Box<ParsedExpr>),
    Root(usize, Box<ParsedExpr>),
}

pub(crate) enum Report {
    Coefficient(usize),
    Zeros,
}

impl Report {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            Self::Coefficient(_) => "rvMazurDiscCoeff",
            Self::Zeros => "rvMazurDiscZeros",
        }
    }
}

#[derive(Debug, Eq, PartialEq, Hash)]
pub(crate) struct Grade {
    source: amice::Grade,
    center: i128,
    depth: usize,
    scale: u64,
    precision: usize,
    modulus: u64,
}

impl Grade {
    pub(crate) fn new(source: amice::Grade, Chart(center, depth): Chart) -> Result<Self> {
        let prime = source.prime;
        let p = u128::from(prime.get());
        if !center.unsigned_abs().is_multiple_of(p) {
            return Err(Error::plain("Amice weight disc center must belong to pZ_p"));
        }
        if depth == 0 {
            return Err(Error::plain("Amice weight disc depth must be at least one"));
        }
        let scale = prime_power(prime, depth).map_err(|_| radius_limit())?;
        let v = valuation(center.unsigned_abs(), p, depth);
        let precision = source.precision.min(source.principal_level() + v);
        let modulus = prime_power(prime, precision)?;
        Ok(Self {
            source,
            center,
            depth,
            scale,
            precision,
            modulus,
        })
    }

    pub(crate) fn restrict(&self, Chart(offset, refinement): Chart) -> Result<Self> {
        if refinement == 0 {
            return Err(Error::plain(
                "Amice weight disc refinement must be at least one",
            ));
        }
        let depth = self
            .depth
            .checked_add(refinement)
            .ok_or_else(radius_limit)?;
        let center = Integer::from_i128(offset)
            .multiply_wide(Integer::from_i128(i128::from(self.scale)))
            .and_then(|offset| offset.add(Integer::from_i128(self.center)))
            .ok()
            .and_then(Integer::signed)
            .ok_or_else(|| {
                Error::plain("Amice weight disc center exceeds the signed 128-bit range")
            })?;
        Self::new(self.source.clone(), Chart(center, depth))
    }

    pub(crate) fn tate_precision(&self) -> rv_tate::Precision {
        rv_tate::Precision {
            prime: self.source.prime,
            digits: self.precision,
            modulus: self.modulus,
        }
    }

    pub(crate) fn root_possible(&self) -> bool {
        // Every nonconstant coefficient is divisible by the chart scale.
        self.precision > self.depth
    }

    pub(crate) fn checked_residue(&self, residue: usize) -> Result<u32> {
        let p = self.source.prime.get();
        if residue >= p as usize {
            return Err(Error::plain(format!(
                "Amice root residue must be between 0 and {}",
                p - 1
            )));
        }
        Ok(residue as u32)
    }

    pub(crate) fn type_name(&self, root: bool) -> String {
        let name = if root {
            "RVMazurDiscRoot"
        } else {
            "RVMazurDisc"
        };
        format!(
            "{name}({},{},{},{},{},{})",
            self.source.ordinary,
            self.source.level,
            self.source.index,
            self.source.character,
            self.center,
            self.depth,
        )
    }
}

pub(crate) struct Disc {
    source: Rc<Amice>,
    coefficients: Vec<u64>,
}

impl Disc {
    pub(crate) fn new(model: &mut Model, grade: &Grade, source: Rc<Amice>) -> Result<Self> {
        let width = 1 + (grade.precision - 1) / grade.depth;
        let coefficients = pullback(model, source.coefficients(), grade, width)?;
        Ok(Self {
            source,
            coefficients,
        })
    }

    pub(crate) fn restrict(&self, model: &mut Model, target: &Grade) -> Result<Self> {
        // Reuse the full finite-level representative. Its finer-disc precision can
        // exceed the precision of the original chart's truncated pullback.
        Self::new(model, target, self.source.clone())
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(2 + self.coefficients.len())?;
        Ok(Value::MazurDisc(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        self.source.compare(model, &rhs.source)
    }

    pub(crate) fn coefficients(&self) -> &[u64] {
        &self.coefficients
    }

    pub(crate) fn data(&self, grade: &Grade) -> String {
        format!(
            "precision={}, modulus={}, coefficients={:?}",
            grade.precision, grade.modulus, self.coefficients,
        )
    }

    pub(crate) fn report(
        &self,
        model: &mut Model,
        grade: &Grade,
        report: Report,
    ) -> Result<String> {
        let name = report.name();
        match report {
            Report::Coefficient(j) => Ok(format!(
                "{name}[{j}] = precision={}, modulus={}, residue={}",
                grade.precision,
                grade.modulus,
                self.coefficients.get(j).copied().unwrap_or(0),
            )),
            Report::Zeros => {
                let Some(normalized) = self.normalize(grade) else {
                    return Ok(format!(
                        "{name} = precision={}, unresolved",
                        grade.precision
                    ));
                };
                let reduction = &normalized.reduction;
                let degree = reduction.len() - 1;
                let roots = rational_residue_roots(model, grade.source.prime, reduction)?;
                let nonrational = degree
                    - roots
                        .iter()
                        .map(|(_, multiplicity)| multiplicity)
                        .sum::<usize>();
                Ok(format!(
                    "{name} = precision={}, mu={}, lambda={degree}, reduction={reduction:?}, residual_roots={roots:?}, roots_in_nonrational_residues={nonrational}",
                    grade.precision, normalized.valuation,
                ))
            }
        }
    }

    pub(crate) fn root(&self, model: &mut Model, grade: &Grade, residue: u32) -> Result<Root> {
        let normalized = self.normalize(grade).ok_or_else(|| {
            Error::plain("Amice weight disc is unresolved at the available precision")
        })?;
        let prime = grade.source.prime;
        let p = u64::from(prime.get());
        let (value, derivative) =
            evaluate_modular_with_derivative(&normalized.reduction, prime, residue);
        if value != 0 {
            return Err(Error::plain(format!(
                "Amice weight disc residue {residue} has no root in the normalized reduction"
            )));
        }
        if derivative == 0 {
            let multiplicity = residue_multiplicity(&normalized.reduction, prime, residue);
            return Err(Error::plain(format!(
                "Amice weight disc residue {residue} has multiplicity {multiplicity}. A simple Hensel root is required"
            )));
        }

        let divisor = prime_power(prime, normalized.valuation)?;
        let modulus = grade.modulus / divisor;
        let inverse_derivative = prime_linear::invert(model, prime, derivative)?;
        let x = lift_simple_root(
            model,
            prime,
            residue,
            inverse_derivative,
            modulus,
            |model, x, power| {
                model.spend(self.coefficients.len())?;
                Ok(self
                    .coefficients
                    .iter()
                    .rev()
                    .fold(0, |value, &coefficient| {
                        add_mod(
                            mul_mod(value, x, power),
                            (coefficient / divisor) % power,
                            power,
                        )
                    }))
            },
        )?;
        let t_modulus = u128::from(grade.scale) * u128::from(modulus);
        let t = translated_residue(grade.center, grade.scale, x, t_modulus);
        debug_assert_eq!(x % p, u64::from(residue));
        Ok(Root {
            x,
            x_precision: grade.precision - normalized.valuation,
            x_modulus: modulus,
            t,
            t_precision: grade.depth + grade.precision - normalized.valuation,
            t_modulus,
        })
    }

    fn normalize(&self, grade: &Grade) -> Option<Normalized> {
        let p = u128::from(grade.source.prime.get());
        let order = self
            .coefficients
            .iter()
            .map(|&c| valuation(u128::from(c), p, grade.precision))
            .min()?;
        if order == grade.precision {
            return None;
        }
        let divisor = (p as u64).pow(order as u32);
        let mut reduction: Vec<_> = self
            .coefficients
            .iter()
            .map(|c| ((c / divisor) % p as u64) as u32)
            .collect();
        while reduction.last() == Some(&0) {
            reduction.pop();
        }
        Some(Normalized {
            valuation: order,
            reduction,
        })
    }
}

struct Normalized {
    valuation: usize,
    reduction: Vec<u32>,
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Root {
    x: u64,
    x_precision: usize,
    x_modulus: u64,
    t: u128,
    t_precision: usize,
    t_modulus: u128,
}

impl Root {
    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(9)?;
        Ok(Value::MazurDiscRoot(Rc::new(self)))
    }

    pub(crate) fn data(&self) -> String {
        format!(
            "x={}, x_precision={}, x_modulus={}, t={}, t_precision={}, t_modulus={}",
            self.x, self.x_precision, self.x_modulus, self.t, self.t_precision, self.t_modulus,
        )
    }
}

fn pullback(
    model: &mut Model,
    coefficients: &[u64],
    grade: &Grade,
    width: usize,
) -> Result<Vec<u64>> {
    model.spend(coefficients.len().saturating_mul(width))?;
    let modulus = grade.modulus;
    let constant = add_mod(1, reduce(grade.center, modulus), modulus);
    let linear = grade.scale % modulus;
    let mut result = vec![0; width];
    // The existing group-ring variable is 1+T. Substitute 1+a+p^e X.
    for &coefficient in coefficients.iter().rev() {
        for j in (1..width).rev() {
            result[j] = add_mod(
                mul_mod(constant, result[j], modulus),
                mul_mod(linear, result[j - 1], modulus),
                modulus,
            );
        }
        result[0] = add_mod(mul_mod(constant, result[0], modulus), coefficient, modulus);
    }
    Ok(result)
}

fn rational_residue_roots(
    model: &mut Model,
    prime: Prime,
    coefficients: &[u32],
) -> Result<Vec<(u32, usize)>> {
    if coefficients.len() == 1 {
        return Ok(Vec::new());
    }
    model.spend((prime.get() as usize).saturating_mul(coefficients.len()))?;
    let mut roots = Vec::new();
    for residue in 0..prime.get() {
        let (value, derivative) = evaluate_modular_with_derivative(coefficients, prime, residue);
        if value == 0 {
            let multiplicity = if derivative != 0 {
                1
            } else {
                model.spend(coefficients.len().saturating_mul(coefficients.len()))?;
                residue_multiplicity(coefficients, prime, residue)
            };
            roots.push((residue, multiplicity));
        }
    }
    Ok(roots)
}

fn evaluate_modular_with_derivative(coefficients: &[u32], prime: Prime, x: u32) -> (u32, u32) {
    let (mut value, mut derivative) = (0, 0);
    for &coefficient in coefficients.iter().rev() {
        derivative = prime.add(prime.mul(derivative, x), value);
        value = prime.add(prime.mul(value, x), coefficient);
    }
    (value, derivative)
}

fn residue_multiplicity(coefficients: &[u32], prime: Prime, x: u32) -> usize {
    let mut current = coefficients.to_vec();
    let mut multiplicity = 0;
    while current.len() > 1 {
        let mut quotient_term = 0;
        for coefficient in current.iter_mut().rev() {
            let term = prime.add(*coefficient, prime.mul(quotient_term, x));
            *coefficient = quotient_term;
            quotient_term = term;
        }
        if quotient_term != 0 {
            break;
        }
        // Synthetic division left the quotient in the lower coefficients.
        current.pop();
        multiplicity += 1;
    }
    multiplicity
}

fn valuation(mut value: u128, p: u128, cap: usize) -> usize {
    if value == 0 {
        return cap;
    }
    let mut valuation = 0;
    while valuation < cap && value.is_multiple_of(p) {
        value /= p;
        valuation += 1;
    }
    valuation
}

fn translated_residue(center: i128, scale: u64, x: u64, modulus: u128) -> u128 {
    let magnitude = center.unsigned_abs() % modulus;
    let a = if center < 0 && magnitude != 0 {
        modulus - magnitude
    } else {
        magnitude
    };
    let step = u128::from(scale) * u128::from(x);
    // Both terms are already reduced, so avoid overflowing even when the
    // final modulus uses the full unsigned 128-bit range.
    if a >= modulus - step {
        a - (modulus - step)
    } else {
        a + step
    }
}

fn radius_limit() -> Error {
    Error::plain("Amice weight disc radius exceeds the unsigned 64-bit prime-power range")
}
