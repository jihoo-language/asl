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
const FACTORS: [(&str, Factor); 7] = [
    ("2", Factor::Prime(2)),
    ("t", Factor::Root(0)),
    ("(t-1)", Factor::Root(1)),
    ("(t-2)", Factor::Root(2)),
    ("3", Factor::Prime(3)),
    ("(t-3)", Factor::Root(3)),
    ("(t-4)", Factor::Root(4)),
];

#[derive(Clone, Copy, Eq, PartialEq)]
pub(super) enum Factor {
    Prime(u8),
    Root(u8),
}
type Dense = Vec<Integer>;

#[derive(Clone, Copy, Default, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Denominator([u8; FACTORS.len()]);

impl Denominator {
    pub(crate) fn new(powers: [u128; 4]) -> Result<Self> {
        Self::extended([powers[0], powers[1], powers[2], powers[3], 0, 0, 0])
    }

    pub(crate) fn from_literal(powers: Vec<u128>) -> Result<Self> {
        match powers.as_slice() {
            [a, b, c, d] => Self::new([*a, *b, *c, *d]),
            [a, b, c, d, e, f, g] => Self::extended([*a, *b, *c, *d, *e, *f, *g]),
            _ => Err(Error::plain(
                "partition denominator requires four or seven exponents",
            )),
        }
    }

    pub(super) fn factor(value: Factor) -> Result<Self> {
        let index = FACTORS
            .iter()
            .position(|&(_, factor)| factor == value)
            .ok_or_else(|| {
                Error::plain("path eigenvalue difference lies outside the diagram localization")
            })?;
        let mut powers = [0; FACTORS.len()];
        powers[index] = 1;
        Ok(Self(powers))
    }

    fn extended(powers: [u128; FACTORS.len()]) -> Result<Self> {
        let limit = u128::from(DENOMINATOR_LIMIT);
        if powers.iter().any(|&x| x > limit) || powers.iter().sum::<u128>() > limit {
            return Err(denominator_limit());
        }
        Ok(Self(powers.map(|x| x as u8)))
    }

    fn product(self, rhs: Self) -> Result<Self> {
        Self::extended(std::array::from_fn(|i| {
            u128::from(self.0[i]) + u128::from(rhs.0[i])
        }))
    }

    fn lcm(self, rhs: Self) -> Result<Self> {
        Self::extended(std::array::from_fn(|i| u128::from(self.0[i].max(rhs.0[i]))))
    }

    fn render(self) -> String {
        let factors = self
            .0
            .iter()
            .zip(FACTORS)
            .filter_map(|(&exponent, (factor, _))| match exponent {
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
        if self.0.iter().zip(FACTORS).any(
            |(&e, (_, factor))| matches!(factor, Factor::Root(root) if e!=0 && u128::from(root)==n),
        ) {
            return Err(Error::plain(format!(
                "partition rational coefficient has a pole at t={n}"
            )));
        }
        Ok(())
    }

    fn at(self, model: &mut Model, n: u128) -> Result<Integer> {
        self.require_regular(n)?;
        let mut result = Integer::from_i128(1);
        for (&exponent, (_, factor)) in self.0.iter().zip(FACTORS).filter(|(e, _)| **e != 0) {
            let value = match factor {
                Factor::Prime(p) => natural(u128::from(p)),
                Factor::Root(root) => natural(n)
                    .add(natural(u128::from(root)).negate())
                    .map_err(|_| intermediate_overflow())?,
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
        for (index, (_, factor)) in FACTORS.into_iter().enumerate() {
            let missing = common.0[index] - self.denominator.0[index];
            if missing == 0 {
                continue;
            }
            match factor {
                Factor::Prime(p) => {
                    let scale = i128::from(p)
                        .checked_pow(u32::from(missing))
                        .ok_or_else(coefficient_overflow)?;
                    numerator = numerator.unary(model, Unary::Scale(scale, 0))?;
                }
                Factor::Root(0) => {
                    numerator = numerator.unary(model, Unary::Scale(1, u128::from(missing)))?
                }
                Factor::Root(root) => {
                    for _ in 0..missing {
                        let leading = numerator.unary(model, Unary::Scale(1, 1))?;
                        let constant =
                            numerator.unary(model, Unary::Scale(-i128::from(root), 0))?;
                        numerator = leading.add(model, &constant)?;
                    }
                }
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
        let (diagrams, mut coefficients): (Vec<_>, Vec<_>) =
            split_dense(model, numerator)?.into_iter().unzip();
        let denominator = reduce(model, &mut coefficients, denominator)?;
        let numerator = pack_dense(rank, diagrams.into_iter().zip(coefficients))?;
        Ok(Self {
            numerator,
            denominator,
        })
    }

    pub(super) fn trace_scalar(&self, model: &mut Model) -> Result<Self> {
        let polynomial = self.numerator.trace_polynomial(model)?;
        Self::scalar(model, polynomial, self.denominator)
    }

    pub(super) fn trace_pair(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        debug_assert_eq!(self.rank(), rhs.rank());
        if self.is_zero() || rhs.is_zero() {
            return Ok(Self::zero(0));
        }
        let rank = self.rank();
        let left = self.numerator.diagram_groups(model)?;
        let right: Vec<_> = rhs
            .numerator
            .diagram_groups(model)?
            .into_iter()
            .map(|(diagram, polynomial)| (diagram.reflect(rank), polynomial))
            .collect();
        let mut polynomial = Polynomial::new();
        for (a, p) in &left {
            for (b, q) in &right {
                model.spend(1 + 2 * rank + p.len() * q.len())?;
                let loops = a.joint_components(b, rank) as u128;
                for &(i, x) in p {
                    for &(j, y) in q {
                        super::accumulate(&mut polynomial, i + j + loops, Integer::product(x, y))?;
                    }
                }
            }
        }
        Self::scalar(
            model,
            polynomial,
            self.denominator.product(rhs.denominator)?,
        )
    }

    fn scalar(model: &mut Model, polynomial: Polynomial, denominator: Denominator) -> Result<Self> {
        let mut dense = Vec::new();
        for (degree, value) in polynomial {
            insert_dense(model, &mut dense, degree, value)?;
        }
        let denominator = reduce(model, std::slice::from_mut(&mut dense), denominator)?;
        let numerator = pack_dense(0, [(Diagram::identity(0), dense)])?;
        Ok(Self {
            numerator,
            denominator,
        })
    }

    pub(super) fn divide_by_scalar(&self, model: &mut Model, divisor: &Self) -> Result<Self> {
        debug_assert_eq!(divisor.rank(), 0);
        if divisor.is_zero() {
            return Err(Error::plain("division by a zero partition scalar"));
        }
        if self.is_zero() {
            return Ok(self.clone());
        }
        let mut factor = split_dense(model, divisor.numerator.clone())?
            .into_values()
            .next()
            .unwrap();
        let mut powers = [0; FACTORS.len()];
        for (index, (_, kind)) in FACTORS.into_iter().enumerate() {
            model.spend(factor.len())?;
            while let Some(quotient) = divide_factor(&factor, kind)? {
                powers[index] += 1;
                if powers[index] > u128::from(DENOMINATOR_LIMIT) {
                    return Err(denominator_limit());
                }
                factor = quotient;
                model.spend(factor.len())?;
            }
        }
        // The divisor denominator is an integral polynomial in the
        // localization generators. Multiplying it back uses the same lift
        // as ordinary rational addition. Cancel any remaining factor exactly.
        let numerator = self.lift(model, self.denominator.product(divisor.denominator)?)?;
        let numerator = if factor.len() == 1 && factor[0].signed() == Some(1) {
            numerator
        } else {
            let mut polys = split_dense(model, numerator)?;
            for polynomial in polys.values_mut() {
                *polynomial = exact_quotient(model, polynomial, &factor)?;
            }
            pack_dense(self.rank(), polys)?
        };
        let denominator = self.denominator.product(Denominator::extended(powers)?)?;
        Self::normalized(model, numerator, denominator)
    }

    pub(super) fn rotate(&self, model: &mut Model) -> Result<Self> {
        let rank = self.rank();
        model.spend((1 + 2 * rank) * self.numerator.terms.len())?;
        let terms = self
            .numerator
            .terms
            .iter()
            .map(|((diagram, degree), &value)| ((diagram.rotate(rank), *degree), value))
            .collect();
        Ok(Self {
            numerator: Algebra { rank, terms },
            denominator: self.denominator,
        })
    }

    pub(super) fn require_regular_at(&self, n: u128) -> Result<()> {
        self.denominator.require_regular(n)
    }

    pub(super) fn scalar_display(
        &self,
        model: &mut Model,
        dimension: Option<u128>,
    ) -> Result<String> {
        debug_assert_eq!(self.rank(), 0);
        let polynomial = self.numerator.trace_polynomial(model)?;
        match dimension {
            Some(n) => {
                let denominator = self.denominator.at(model, n)?;
                let value = evaluate_polynomial(model, polynomial, n)?;
                rational(model, value, denominator)
            }
            None if polynomial.keys().all(|&degree| degree == 0)
                && self
                    .denominator
                    .0
                    .iter()
                    .zip(FACTORS)
                    .all(|(&power, (_, factor))| {
                        power == 0 || matches!(factor, Factor::Prime(_))
                    }) =>
            {
                let denominator = self.denominator.at(model, 0)?;
                let value = polynomial
                    .get(&0)
                    .copied()
                    .unwrap_or_else(|| Integer::from_i128(0));
                rational(model, value, denominator)
            }
            None if self.denominator == Denominator::default() => Ok(polynomial_report(polynomial)),
            None => Ok(format!(
                "({})/{}",
                polynomial_report(polynomial),
                self.denominator.render()
            )),
        }
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

fn split_dense(model: &mut Model, numerator: Algebra) -> Result<BTreeMap<Diagram, Dense>> {
    let mut polys = BTreeMap::<Diagram, Dense>::new();
    for ((diagram, degree), coefficient) in numerator.terms {
        insert_dense(
            model,
            polys.entry(diagram).or_default(),
            degree,
            Integer::from_i128(coefficient),
        )?;
    }
    Ok(polys)
}

fn pack_dense(
    rank: usize,
    diagrams: impl IntoIterator<Item = (Diagram, Dense)>,
) -> Result<Algebra> {
    let mut terms = BTreeMap::new();
    for (diagram, coefficients) in diagrams {
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
    Ok(Algebra { rank, terms })
}

fn exact_quotient(model: &mut Model, dividend: &[Integer], divisor: &[Integer]) -> Result<Dense> {
    let nonlocal = || Error::plain("Deligne scalar quotient lies outside the diagram localization");
    if dividend.len() < divisor.len() {
        return Err(nonlocal());
    }
    let mut remainder = dividend.to_vec();
    let mut quotient = vec![Integer::from_i128(0); dividend.len() - divisor.len() + 1];
    let leading = *divisor.last().unwrap();
    for (degree, value) in quotient.iter_mut().enumerate().rev() {
        let top = remainder[degree + divisor.len() - 1];
        if top.is_zero() {
            continue;
        }
        *value = top.divide_exact(model, leading)?.ok_or_else(nonlocal)?;
        model.spend(divisor.len())?;
        for (coefficient, &factor) in remainder[degree..].iter_mut().zip(divisor) {
            *coefficient = coefficient.add(value.multiply_wide(factor)?.negate())?;
        }
    }
    if remainder.iter().any(|coefficient| !coefficient.is_zero()) {
        return Err(nonlocal());
    }
    Ok(quotient)
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
    for (index, (_, factor)) in FACTORS.into_iter().enumerate() {
        while denominator.0[index] > 0 {
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
            denominator.0[index] -= 1;
        }
    }
    Ok(denominator)
}

fn divide_factor(polynomial: &[Integer], factor: Factor) -> Result<Option<Dense>> {
    if polynomial.is_empty() {
        return Ok(Some(Vec::new()));
    }
    let Factor::Root(root) = factor else {
        let Factor::Prime(p) = factor else {
            unreachable!()
        };
        let mut quotient = Vec::with_capacity(polynomial.len());
        for value in polynomial {
            let (magnitude, remainder) = value.magnitude.divide_u32(u32::from(p));
            if remainder != 0 {
                return Ok(None);
            }
            quotient.push(Integer {
                magnitude,
                negative: value.negative,
            });
        }
        return Ok(Some(quotient));
    };
    let root = i128::from(root);
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
