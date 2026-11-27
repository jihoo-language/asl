use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::elliptic_modular::decimal;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_witt_spectrum::norm_power;

use super::{
    coefficient_overflow, evaluate_polynomial, intermediate_overflow, natural, polynomial_report,
    Algebra, Binary, Diagram, Integer, Polynomial, Query, Unary,
};

const DEGREE_LIMIT: u128 = 128;
const DENOMINATOR_LIMIT: u8 = 64;
const FACTORS: [&str; 4] = ["2", "t", "(t-1)", "(t-2)"];
type Dense = Vec<Integer>;

#[derive(Clone, Copy, Default, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Denominator([u8; 4]);

impl Denominator {
    pub(crate) fn new(powers: [u128; 4]) -> Result<Self> {
        let limit = u128::from(DENOMINATOR_LIMIT);
        if powers.iter().any(|&x| x > limit) || powers.iter().sum::<u128>() > limit {
            return Err(denominator_limit());
        }
        Ok(Self(powers.map(|x| x as u8)))
    }

    fn product(self, rhs: Self) -> Result<Self> {
        Self::new(std::array::from_fn(|i| {
            u128::from(self.0[i]) + u128::from(rhs.0[i])
        }))
    }

    fn lcm(self, rhs: Self) -> Result<Self> {
        Self::new(std::array::from_fn(|i| u128::from(self.0[i].max(rhs.0[i]))))
    }

    fn render(self) -> String {
        let factors = self
            .0
            .iter()
            .zip(FACTORS)
            .filter_map(|(&exponent, factor)| match exponent {
                0 => None,
                1 => Some(factor.to_string()),
                _ => Some(format!("{factor}^{exponent}")),
            })
            .collect::<Vec<_>>();
        if factors.is_empty() {
            "1".into()
        } else {
            factors.join("*")
        }
    }

    fn require_regular(self, n: u128) -> Result<()> {
        if n < 3 && self.0[n as usize + 1] != 0 {
            return Err(Error::plain(format!(
                "partition rational coefficient has a pole at t={n}"
            )));
        }
        Ok(())
    }

    fn at(self, model: &mut Model, n: u128) -> Result<Integer> {
        self.require_regular(n)?;
        let mut result = Integer::from_i128(1);
        for (i, &exponent) in self.0.iter().enumerate().filter(|(_, e)| **e != 0) {
            let value = if i == 0 {
                natural(2)
            } else {
                natural(n)
                    .add(natural((i - 1) as u128).negate())
                    .map_err(|_| intermediate_overflow())?
            };
            let factor = norm_power(model, value, u128::from(exponent))?;
            result = result
                .multiply_wide(factor)
                .map_err(|_| intermediate_overflow())?;
        }
        Ok(result)
    }
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct LocalAlgebra {
    numerator: Algebra,
    denominator: Denominator,
}

impl LocalAlgebra {
    pub(crate) fn zero(rank: usize) -> Self {
        Self {
            numerator: Algebra::zero(rank),
            denominator: Denominator::default(),
        }
    }

    pub(crate) fn identity(rank: usize) -> Self {
        Self {
            numerator: Algebra::identity(rank),
            denominator: Denominator::default(),
        }
    }

    pub(crate) fn embed(numerator: Algebra) -> Result<Self> {
        if numerator.terms.keys().any(|key| key.1 > DEGREE_LIMIT) {
            return Err(degree_limit());
        }
        Ok(Self {
            numerator,
            denominator: Denominator::default(),
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.numerator.rank
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.numerator.terms.is_empty()
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        self.charge(model)?;
        Ok(Value::RVPartitionLocal(Rc::new(self)))
    }

    pub(crate) fn charge(&self, model: &mut Model) -> Result<()> {
        model.allocate(8 + (8 + 2 * self.rank()) * self.numerator.terms.len())
    }

    pub(crate) fn render(&self) -> String {
        format!(
            "RVPartitionLocal numerator={}, denominator={}",
            self.numerator.render(),
            self.denominator.render()
        )
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        let compare = self.numerator.compare(model, &rhs.numerator)?;
        Ok(compare.then_with(|| self.denominator.cmp(&rhs.denominator)))
    }

    pub(crate) fn divide(&self, model: &mut Model, denominator: Denominator) -> Result<Self> {
        if self.is_zero() {
            return Ok(self.clone());
        }
        Self::normalized(
            model,
            self.numerator.clone(),
            self.denominator.product(denominator)?,
        )
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let common = self.denominator.lcm(rhs.denominator)?;
        let a = self.lift(model, common)?;
        let b = rhs.lift(model, common)?;
        let numerator = a.add(model, &b)?;
        Self::normalized(model, numerator, common)
    }

    pub(crate) fn binary(&self, model: &mut Model, operation: Binary, rhs: &Self) -> Result<Self> {
        if self.is_zero() || rhs.is_zero() {
            return Ok(Self::zero(operation.rank(self.rank(), rhs.rank())?));
        }
        let denominator = self.denominator.product(rhs.denominator)?;
        let numerator = self.numerator.binary(model, operation, &rhs.numerator)?;
        Self::normalized(model, numerator, denominator)
    }

    pub(crate) fn unary(&self, model: &mut Model, operation: Unary) -> Result<Self> {
        let numerator = self.numerator.unary(model, operation)?;
        Self::normalized(model, numerator, self.denominator)
    }

    fn lift(&self, model: &mut Model, common: Denominator) -> Result<Algebra> {
        let mut numerator = self.numerator.clone();
        let (scale, degree) = (
            common.0[0] - self.denominator.0[0],
            common.0[1] - self.denominator.0[1],
        );
        if scale != 0 || degree != 0 {
            numerator = numerator.unary(model, Unary::Scale(1i128 << scale, u128::from(degree)))?;
        }
        for factor in 2..=3 {
            for _ in self.denominator.0[factor]..common.0[factor] {
                let leading = numerator.unary(model, Unary::Scale(1, 1))?;
                let constant = numerator.unary(model, Unary::Scale(-((factor - 1) as i128), 0))?;
                numerator = leading.add(model, &constant)?;
            }
        }
        Ok(numerator)
    }

    fn normalized(model: &mut Model, numerator: Algebra, denominator: Denominator) -> Result<Self> {
        if numerator.terms.is_empty() {
            return Ok(Self::zero(numerator.rank));
        }
        if denominator == Denominator::default() {
            return Self::embed(numerator);
        }
        let rank = numerator.rank;
        let mut polys = BTreeMap::<Diagram, Dense>::new();
        for ((diagram, degree), coefficient) in numerator.terms {
            insert_dense(
                model,
                polys.entry(diagram).or_default(),
                degree,
                Integer::from_i128(coefficient),
            )?;
        }
        let (diagrams, mut coefficients): (Vec<_>, Vec<_>) = polys.into_iter().unzip();
        let denominator = reduce(model, &mut coefficients, denominator)?;
        let mut terms = BTreeMap::new();
        for (diagram, coefficients) in diagrams.into_iter().zip(coefficients) {
            if coefficients.len() as u128 > DEGREE_LIMIT + 1 {
                return Err(degree_limit());
            }
            for (degree, coefficient) in coefficients
                .into_iter()
                .enumerate()
                .filter(|(_, c)| !c.is_zero())
            {
                let coefficient = coefficient.signed().ok_or_else(coefficient_overflow)?;
                super::require_terms(terms.len() + 1)?;
                terms.insert((diagram.clone(), degree as u128), coefficient);
            }
        }
        Ok(Self {
            numerator: Algebra { rank, terms },
            denominator,
        })
    }

    fn trace_polynomial(&self, model: &mut Model) -> Result<(Polynomial, Denominator)> {
        let mut poly = Vec::new();
        for (degree, coefficient) in self.numerator.trace_polynomial(model)? {
            insert_dense(model, &mut poly, degree, coefficient)?;
        }
        let denominator = reduce(model, std::slice::from_mut(&mut poly), self.denominator)?;
        let poly = poly
            .into_iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|(d, c)| (d as u128, c))
            .collect();
        Ok((poly, denominator))
    }

    pub(crate) fn trace_report(&self, model: &mut Model, name: &str) -> Result<String> {
        let (polynomial, denominator) = self.trace_polynomial(model)?;
        Ok(format!(
            "{name} = numerator={}, denominator={}",
            polynomial_report(polynomial),
            denominator.render()
        ))
    }

    pub(crate) fn query(&self, model: &mut Model, query: Query) -> Result<String> {
        match query {
            Query::Coefficient(..) => Ok(format!(
                "{} (numerator), denominator={}",
                self.numerator.query(model, query)?,
                self.denominator.render()
            )),
            Query::Trace => self.trace_report(model, "rvPartitionTrace"),
            Query::TraceAt(n) => {
                let (polynomial, denominator) = self.trace_polynomial(model)?;
                let den = denominator.at(model, n)?;
                let value = evaluate_polynomial(model, polynomial, n)?;
                Ok(format!(
                    "rvPartitionTraceAt[{n}] = {}",
                    rational(model, value, den)?
                ))
            }
            Query::ZeroAt(n) => {
                self.denominator.require_regular(n)?;
                self.numerator.query(model, query)
            }
            Query::Entry(n, upper, lower) => {
                let den = self.denominator.at(model, n)?;
                let value = self.numerator.entry(model, n, upper, lower)?;
                Ok(format!(
                    "rvPartitionEntry[{n}] = {}",
                    rational(model, value, den)?
                ))
            }
        }
    }
}

fn insert_dense(
    model: &mut Model,
    polynomial: &mut Dense,
    degree: u128,
    coefficient: Integer,
) -> Result<()> {
    if degree > 2 * DEGREE_LIMIT + super::MAX_RANK as u128 {
        return Err(degree_limit());
    }
    let size = (degree + 1) as usize;
    if polynomial.len() < size {
        model.spend(size - polynomial.len())?;
        polynomial.resize(size, Integer::from_i128(0));
    }
    polynomial[size - 1] = coefficient;
    Ok(())
}

fn reduce(
    model: &mut Model,
    polys: &mut [Dense],
    mut denominator: Denominator,
) -> Result<Denominator> {
    if polys.iter().all(Vec::is_empty) {
        return Ok(Denominator::default());
    }
    for factor in 0..4 {
        while denominator.0[factor] > 0 {
            let mut quotients = Vec::with_capacity(polys.len());
            for poly in polys.iter() {
                model.spend(poly.len())?;
                let Some(quotient) = divide_factor(poly, factor)? else {
                    break;
                };
                quotients.push(quotient);
            }
            if quotients.len() != polys.len() {
                break;
            }
            for (poly, quotient) in polys.iter_mut().zip(quotients) {
                *poly = quotient;
            }
            denominator.0[factor] -= 1;
        }
    }
    Ok(denominator)
}

fn divide_factor(polynomial: &[Integer], factor: usize) -> Result<Option<Dense>> {
    if polynomial.is_empty() {
        return Ok(Some(Vec::new()));
    }
    if factor == 0 {
        let mut quotient = Vec::with_capacity(polynomial.len());
        for value in polynomial {
            let (magnitude, remainder) = value.magnitude.divide_u32(2);
            if remainder != 0 {
                return Ok(None);
            }
            quotient.push(Integer {
                magnitude,
                negative: value.negative,
            });
        }
        return Ok(Some(quotient));
    }
    let root = (factor - 1) as i128;
    let mut quotient = vec![Integer::from_i128(0); polynomial.len() - 1];
    let mut top = Integer::from_i128(0);
    for (index, &value) in polynomial.iter().enumerate().rev() {
        top = top
            .multiply(root)
            .and_then(|carry| value.add(carry))
            .map_err(|_| intermediate_overflow())?;
        if index != 0 {
            quotient[index - 1] = top;
        }
    }
    Ok(top.is_zero().then_some(quotient))
}

fn rational(model: &mut Model, value: Integer, den: Integer) -> Result<String> {
    if value.is_zero() {
        return Ok("0".into());
    }
    let gcd = value.magnitude.gcd(model, den.magnitude)?;
    let magnitude = value
        .magnitude
        .divide_exact(model, gcd)?
        .ok_or_else(intermediate_overflow)?;
    let numerator = Integer {
        magnitude,
        negative: value.negative != den.negative,
    };
    let denominator = den
        .magnitude
        .divide_exact(model, gcd)?
        .ok_or_else(intermediate_overflow)?;
    if denominator.small() == Some(1) {
        return Ok(decimal(numerator));
    }
    Ok(format!(
        "{}/{}",
        decimal(numerator),
        decimal(Integer {
            magnitude: denominator,
            negative: false
        })
    ))
}

fn denominator_limit() -> Error {
    Error::plain("partition localization denominator exceeds 64 factors")
}
fn degree_limit() -> Error {
    Error::plain("partition localization numerator exceeds t degree 128")
}
