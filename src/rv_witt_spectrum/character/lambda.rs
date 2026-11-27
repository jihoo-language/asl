use crate::error::{Error, Result};
use crate::model::Model;

use super::super::{add, multiply, one, zero, Integer, SchurShape};
use super::stable::{add_scaled, binomial_product, choose, Expansions, StableCharacter};
use super::table::{accumulate, Coefficients, Partition, MAX_WEIGHT};

pub(crate) enum Functor {
    Adams(u128),
    Exterior(u128),
    Symmetric(u128),
    Schur(SchurShape),
}

impl Functor {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            Self::Adams(_) => "rvStableSchurAdams",
            Self::Exterior(_) => "rvStableSchurExterior",
            Self::Symmetric(_) => "rvStableSchurSymmetric",
            Self::Schur(_) => "rvStableInnerSchur",
        }
    }

    fn weight(&self) -> u128 {
        match self {
            Self::Adams(n) | Self::Exterior(n) | Self::Symmetric(n) => *n,
            Self::Schur(shape) => shape.weight(),
        }
    }

    fn elementary_scalar(&self, value: i128) -> Option<i128> {
        if value == 0 {
            return Some(0);
        }
        let (height, width) = match self {
            Self::Exterior(n) => (*n, 1),
            Self::Symmetric(n) => (1, *n),
            Self::Schur(shape) => (shape.rows().len() as u128, shape.rows()[0]),
            Self::Adams(_) => return Some(value),
        };
        // A negative scalar uses the transposed shape with the sign of
        // its weight. Vanishing beyond the scalar rank comes first.
        let rank = value.unsigned_abs();
        let required = if value > 0 { height } else { width };
        if required > rank {
            return Some(0);
        }
        if rank == 1 {
            return Some(if value < 0 && !self.weight().is_multiple_of(2) {
                -1
            } else {
                1
            });
        }
        None
    }
}

impl StableCharacter {
    pub(crate) fn functor(&self, model: &mut Model, operation: &Functor) -> Result<Self> {
        let weight = operation.weight();
        if let Functor::Adams(d) = operation {
            if *d == 0 {
                return Err(Error::plain("stable Adams index must be positive"));
            }
        }
        let unit = || Ok(Self::seed(SchurShape::new(Vec::new())?));
        if weight == 0 {
            return unit();
        }
        let degree = self.highest_weight();
        if weight == 1 || (degree == 0 && matches!(operation, Functor::Adams(_))) {
            return self.scale(model, 1);
        }
        if degree == 0 {
            let value = self.0.values().next().copied().unwrap_or(0);
            if let Some(result) = operation.elementary_scalar(value) {
                return unit()?.scale(model, result);
            }
        }
        // A cycle of length i in an input monomial can have come from
        // one of length at most d i under powering. The Newton and
        // Jacobi–Trudi degrees respect this weight bound.
        let bound = weight
            .checked_mul(degree.max(1))
            .filter(|&n| n <= u128::from(MAX_WEIGHT))
            .ok_or_else(|| {
                Error::plain("stable Schur operation exceeds the weighted cycle limit of 16")
            })? as u8;
        let actual = if degree == 0 { 0 } else { bound };
        let mut expansions = Expansions::new(model, actual)?;
        let polynomial = expansions.character(model, &self.0)?;
        let result = match operation {
            Functor::Adams(d) => adams(model, &expansions, &polynomial, *d as u8, actual)?,
            _ => {
                let plan = Plan::new(operation, weight as usize);
                let series = newton_series(model, &expansions, &polynomial, degree as u8, &plan)?;
                plan.determinant(model, &series)?
            }
        };
        expansions.reduced_basis(model, actual, result)
    }
}

type CountVector = [u8; MAX_WEIGHT as usize + 1];

fn powered_counts(cycles: &[u8], exponent: u8) -> CountVector {
    let mut counts = [0; MAX_WEIGHT as usize + 1];
    for &cycle in cycles {
        let (mut a, mut b) = (cycle, exponent);
        while b != 0 {
            (a, b) = (b, a % b);
        }
        counts[usize::from(cycle / a)] += a;
    }
    counts
}

fn counts(cycles: &[u8]) -> Vec<(u8, u8)> {
    cycles
        .chunk_by(|a, b| a == b)
        .map(|run| (run[0], run.len() as u8))
        .collect()
}

fn adams(
    model: &mut Model,
    expansions: &Expansions,
    polynomial: &Coefficients,
    exponent: u8,
    bound: u8,
) -> Result<Coefficients> {
    if exponent == 1 || bound == 0 {
        return Ok(polynomial.clone());
    }
    let terms: Vec<_> = polynomial
        .iter()
        .map(|(cycles, &value)| (counts(cycles), value))
        .collect();
    let mut evaluations = Coefficients::new();
    for cycle_type in expansions.cycle_types(bound) {
        let powered = powered_counts(cycle_type, exponent);
        let mut value = zero();
        for (knots, coefficient) in &terms {
            model.spend(knots.len() + 1)?;
            let mut multiplicity = 1i128;
            for &(length, required) in knots {
                let available = powered[usize::from(length)];
                if available < required {
                    multiplicity = 0;
                    break;
                }
                multiplicity *= choose(available, required);
            }
            if multiplicity != 0 {
                let term = multiply(model, *coefficient, Integer::from_i128(multiplicity))?;
                value = add(value, term)?;
            }
        }
        if !value.is_zero() {
            evaluations.insert(cycle_type.clone(), value);
        }
    }
    let mut result = Coefficients::new();
    for cycle_type in expansions.cycle_types(bound) {
        let value =
            finite_difference(model, &counts(cycle_type), &mut Vec::new(), 1, &evaluations)?;
        if !value.is_zero() {
            result.insert(cycle_type.clone(), value);
        }
    }
    Ok(result)
}

fn finite_difference(
    model: &mut Model,
    counts: &[(u8, u8)],
    chosen: &mut Partition,
    factor: i128,
    evaluations: &Coefficients,
) -> Result<Integer> {
    let Some((&(length, count), rest)) = counts.split_first() else {
        return match evaluations.get(chosen) {
            Some(&value) => multiply(model, value, Integer::from_i128(factor)),
            None => Ok(zero()),
        };
    };
    let mut total = zero();
    let start = chosen.len();
    for j in 0..=count {
        model.spend(1 + usize::from(j))?;
        chosen.resize(start + usize::from(j), length);
        let ways = choose(count, j)
            * if (count - j).is_multiple_of(2) {
                factor
            } else {
                -factor
            };
        let term = finite_difference(model, rest, chosen, ways, evaluations)?;
        total = add(total, term)?;
        chosen.truncate(start);
    }
    Ok(total)
}

fn unit_polynomial() -> Coefficients {
    Coefficients::from([(Vec::new(), one())])
}

fn newton_series(
    model: &mut Model,
    expansions: &Expansions,
    polynomial: &Coefficients,
    degree: u8,
    plan: &Plan,
) -> Result<Vec<Coefficients>> {
    let ghosts = (1..=plan.largest)
        .map(|j| adams(model, expansions, polynomial, j as u8, degree * j as u8))
        .collect::<Result<Vec<_>>>()?;
    let mut values = vec![unit_polynomial()];
    for m in 1..=plan.largest {
        let mut numerator = Coefficients::new();
        for (j, ghost) in ghosts.iter().take(m).enumerate() {
            let product = binomial_product(model, &values[m - j - 1], ghost)?;
            let negative = plan.exterior && !j.is_multiple_of(2);
            add_scaled(
                model,
                &mut numerator,
                &product,
                Integer::from_i128(if negative { -1 } else { 1 }),
            )?;
        }
        for coefficient in numerator.values_mut() {
            model.spend(1)?;
            let (magnitude, remainder) = coefficient.magnitude.divide_u32(m as u32);
            if remainder != 0 {
                return Err(Error::plain("stable Schur Newton division is inexact"));
            }
            coefficient.magnitude = magnitude;
        }
        values.push(numerator);
    }
    Ok(values)
}

struct Plan {
    rows: Vec<usize>,
    size: usize,
    largest: usize,
    exterior: bool,
}

impl Plan {
    fn new(operation: &Functor, size: usize) -> Self {
        let (rows, exterior) = match operation {
            Functor::Exterior(_) => (vec![size], true),
            Functor::Symmetric(_) => (vec![size], false),
            Functor::Schur(shape) => {
                let rows = shape.rows();
                if rows[0] < rows.len() as u128 {
                    (
                        (1..=rows[0])
                            .map(|column| rows.partition_point(|&width| width >= column))
                            .collect(),
                        true,
                    )
                } else {
                    (rows.iter().map(|&width| width as usize).collect(), false)
                }
            }
            Functor::Adams(_) => unreachable!(),
        };
        let largest = rows[0] + rows.len() - 1;
        Self {
            rows,
            size,
            largest,
            exterior,
        }
    }

    fn determinant(&self, model: &mut Model, series: &[Coefficients]) -> Result<Coefficients> {
        let order = self.rows.len();
        if order == 1 {
            return Ok(series[self.size].clone());
        }
        // Each minor uses the first k rows and the columns in mask. Its
        // polynomial degree depends only on that mask. A partial term
        // of order above the shape weight cannot complete to the determinant.
        let end = (1usize << order) - 1;
        let mut minors = vec![None; end + 1];
        minors[0] = Some((0usize, unit_polynomial()));
        for mask in 0..end {
            let Some((degree, minor)) = minors[mask].take() else {
                continue;
            };
            let row = mask.count_ones() as usize;
            for column in 0..order {
                let bit = 1usize << column;
                if mask & bit != 0 {
                    continue;
                }
                let Some(index) = (self.rows[row] + column).checked_sub(row) else {
                    continue;
                };
                if degree + index > self.size || series[index].is_empty() {
                    continue;
                }
                let negative = !(mask >> (column + 1)).count_ones().is_multiple_of(2);
                let product = binomial_product(model, &minor, &series[index])?;
                let (_, sum) =
                    minors[mask | bit].get_or_insert_with(|| (degree + index, Coefficients::new()));
                for (cycle_type, value) in product {
                    let value = if negative { value.negate() } else { value };
                    accumulate(sum, cycle_type, value)?;
                }
            }
        }
        Ok(minors[end]
            .take()
            .map(|(_, result)| result)
            .unwrap_or_default())
    }
}
