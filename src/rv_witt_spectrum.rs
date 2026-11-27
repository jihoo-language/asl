use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};
use crate::wide_integer::Natural;

const LIMBS: usize = 16;
pub(crate) const MAX_DEGREE: usize = 256;
type Integer = Wide<LIMBS>;
type Magnitude = Natural<LIMBS>;
type Polynomial = Vec<Integer>;

mod cyclic;
mod necklace;
mod norm;
mod schur;
pub(crate) use cyclic::CyclicPower;
pub(crate) use necklace::coordinate as primitive_coordinate;
pub(crate) use norm::{cycles as norm_cycles, power as norm_power};
pub(crate) use schur::{Branch as SchurBranch, Shape as SchurShape};

pub(crate) enum Expression {
    Literal(Vec<i128>, Option<Vec<i128>>),
    Line(i128),
    Tensor(Box<[ParsedExpr; 2]>),
    Adams(u128, Box<ParsedExpr>),
    Verschiebung(u128, Box<ParsedExpr>),
    Norm(u128, Box<ParsedExpr>),
    CyclicPower(CyclicPower, Box<ParsedExpr>),
    Schur(SchurShape, Box<ParsedExpr>),
    Exterior(u128, Box<ParsedExpr>),
    Bridge(Bridge),
    HilbertBridge(HilbertBridge),
}

pub(crate) struct HilbertBridge {
    pub(crate) word: Label,
    pub(crate) base: Label,
    pub(crate) length: u128,
    pub(crate) r: u128,
}

pub(crate) struct Bridge {
    pub(crate) word: Label,
    pub(crate) base: Label,
    pub(crate) stratum: Label,
    pub(crate) r: u128,
    pub(crate) i: u128,
    pub(crate) specht: bool,
}

pub(crate) enum Query {
    Data,
    Ghost(u128),
    Factor(u128),
    Necklace(u128),
    NormGhost(u128, u128),
    CyclicGhost(CyclicPower, u128),
    SchurGhost(SchurShape, u128),
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Spectrum {
    numerator: Box<[i128]>,
    denominator: Box<[i128]>,
}

impl Spectrum {
    pub(crate) fn literal(
        model: &mut Model,
        numerator: Vec<i128>,
        denominator: Option<Vec<i128>>,
    ) -> Result<Self> {
        let denominator = denominator.unwrap_or_else(|| vec![1]);
        model.spend(numerator.len() + denominator.len())?;
        if numerator.first() != Some(&1) || denominator.first() != Some(&1) {
            return Err(Error::plain(
                "Witt determinant factors must have constant coefficient 1",
            ));
        }
        let trim = |mut coefficients: Vec<i128>| {
            while coefficients.last() == Some(&0) {
                coefficients.pop();
            }
            coefficients
                .into_iter()
                .map(Integer::from_i128)
                .collect::<Polynomial>()
        };
        let (numerator, denominator) = (trim(numerator), trim(denominator));
        require_degree(numerator.len() - 1)?;
        require_degree(denominator.len() - 1)?;
        Self::canonical(model, numerator, denominator)
    }

    pub(crate) fn zero() -> Self {
        Self {
            numerator: Box::new([1]),
            denominator: Box::new([1]),
        }
    }

    pub(crate) fn line(eigenvalue: i128) -> Result<Self> {
        if eigenvalue == 0 {
            return Err(Error::plain(
                "Witt spectrum lines require a nonzero eigenvalue",
            ));
        }
        Ok(Self {
            numerator: Box::new([
                1,
                eigenvalue.checked_neg().ok_or_else(coefficient_overflow)?,
            ]),
            denominator: Box::new([1]),
        })
    }

    pub(crate) fn effective(coefficients: Polynomial) -> Result<Self> {
        if coefficients.first().and_then(|value| value.signed()) != Some(1)
            || coefficients.last().is_some_and(|value| value.is_zero())
        {
            return Err(integrity());
        }
        require_degree(coefficients.len() - 1)?;
        Ok(Self {
            numerator: narrow(coefficients)?,
            denominator: Box::new([1]),
        })
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(4 + 2 * (self.numerator.len() + self.denominator.len()))?;
        Ok(Value::RVWittSpectrum(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<std::cmp::Ordering> {
        model.spend(
            self.numerator.len().min(other.numerator.len())
                + self.denominator.len().min(other.denominator.len()),
        )?;
        Ok(self.cmp(other))
    }

    pub(crate) fn render(&self) -> String {
        format!(
            "RVWittSpectrum{{{:?}, {:?}}}",
            self.numerator, self.denominator
        )
    }

    pub(crate) fn negate(&self) -> Self {
        Self {
            numerator: self.denominator.clone(),
            denominator: self.numerator.clone(),
        }
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let numerator = product(model, &lift(&self.numerator), &lift(&rhs.numerator))?;
        let denominator = product(model, &lift(&self.denominator), &lift(&rhs.denominator))?;
        Self::canonical(model, numerator, denominator)
    }

    pub(crate) fn tensor(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let left = tensor_factors(model, &self.numerator, &rhs.numerator)?;
        let right = tensor_factors(model, &self.denominator, &rhs.denominator)?;
        let numerator = product(model, &left, &right)?;
        let left = tensor_factors(model, &self.numerator, &rhs.denominator)?;
        let right = tensor_factors(model, &self.denominator, &rhs.numerator)?;
        let denominator = product(model, &left, &right)?;
        Self::canonical(model, numerator, denominator)
    }

    pub(crate) fn adams(&self, model: &mut Model, r: u128) -> Result<Self> {
        if r == 0 {
            return Err(Error::plain("Witt spectrum Adams index must be positive"));
        }
        if r == 1 {
            return Ok(self.clone());
        }
        let numerator = adams_factor(model, &self.numerator, r)?;
        let denominator = adams_factor(model, &self.denominator, r)?;
        Self::canonical(model, numerator, denominator)
    }

    pub(crate) fn verschiebung(&self, model: &mut Model, n: u128) -> Result<Self> {
        if n == 0 {
            return Err(Error::plain(
                "Witt spectrum Verschiebung index must be positive",
            ));
        }
        let substitute = |model: &mut Model, coefficients: &[i128]| -> Result<Box<[i128]>> {
            let degree = coefficients.len() - 1;
            if degree == 0 || n == 1 {
                return Ok(coefficients.into());
            }
            if n > (MAX_DEGREE / degree) as u128 {
                return Err(degree_limit());
            }
            let step = n as usize;
            let length = degree * step + 1;
            model.spend(length)?;
            let mut output = vec![0; length];
            for (i, &coefficient) in coefficients.iter().enumerate() {
                output[i * step] = coefficient;
            }
            Ok(output.into_boxed_slice())
        };
        // A Bezout identity over Q remains an identity after substitution.
        // The two substituted factors are already coprime.
        Ok(Self {
            numerator: substitute(model, &self.numerator)?,
            denominator: substitute(model, &self.denominator)?,
        })
    }

    pub(crate) fn exterior(&self, model: &mut Model, m: u128) -> Result<Self> {
        if m == 0 {
            return Self::line(1);
        }
        if self.denominator.len() != 1 {
            return Err(Error::plain(
                "Witt spectrum exterior powers require an effective spectrum",
            ));
        }
        let n = self.numerator.len() - 1;
        if m == 1 {
            return Ok(self.clone());
        }
        if m > n as u128 {
            return Ok(Self::zero());
        }
        let m = m as usize;
        if m == n {
            let last = Integer::from_i128(self.numerator[n]);
            let coefficient = if m.is_multiple_of(2) {
                last.negate()
            } else {
                last
            };
            return Self::effective(vec![one(), coefficient]);
        }
        let mut rank = 1;
        for i in 1..=m.min(n - m) {
            rank = rank * (n - i + 1) / i;
            require_degree(rank)?;
        }
        let ghosts = power_sums(model, &self.numerator, m * rank)?;
        let mut exterior_ghosts = Vec::with_capacity(rank);
        for k in 1..=rank {
            let sample: Vec<_> = (1..=m).map(|i| ghosts[i * k]).collect();
            let value = *determinant_from_ghosts(model, &sample)?.last().unwrap();
            exterior_ghosts.push(if m.is_multiple_of(2) {
                value
            } else {
                value.negate()
            });
        }
        Self::effective(determinant_from_ghosts(model, &exterior_ghosts)?)
    }

    pub(crate) fn query(&self, model: &mut Model, query: Query) -> Result<String> {
        match query {
            Query::Data => {
                let p = self.numerator.len() - 1;
                let q = self.denominator.len() - 1;
                let rank = p as isize - q as isize;
                Ok(format!("rvWittSpectrumData = numerator degree={p}, denominator degree={q}, virtual rank={rank}"))
            }
            Query::Ghost(k) => {
                if k == 0 {
                    return Err(Error::plain("Witt spectrum ghost index must be positive"));
                }
                let ghost = self
                    .ghost(model, k)?
                    .signed()
                    .ok_or_else(coefficient_overflow)?;
                Ok(format!("rvWittSpectrumGhost[{k}] = {ghost}"))
            }
            Query::Factor(j) => {
                let get = |coefficients: &[i128]| {
                    usize::try_from(j)
                        .ok()
                        .and_then(|index| coefficients.get(index))
                        .copied()
                        .unwrap_or(0)
                };
                Ok(format!(
                    "rvWittSpectrumCoeff[{j}] = numerator={}, denominator={}",
                    get(&self.numerator),
                    get(&self.denominator),
                ))
            }
            Query::NormGhost(n, k) => {
                let trace = self
                    .norm_ghost(model, n, k)?
                    .signed()
                    .ok_or_else(coefficient_overflow)?;
                Ok(format!("rvWittSpectrumNormGhost[{n}, {k}] = {trace}"))
            }
            Query::CyclicGhost(projection, k) => {
                let ghost = self
                    .cyclic_ghost(model, projection, k)?
                    .signed()
                    .ok_or_else(coefficient_overflow)?;
                Ok(format!("{} = {ghost}", projection.ghost_label(k)))
            }
            Query::SchurGhost(shape, k) => {
                let ghost = self
                    .schur_ghost(model, &shape, k)?
                    .signed()
                    .ok_or_else(coefficient_overflow)?;
                Ok(format!(
                    "rvWittSpectrumSchurGhost[{:?}, {k}] = {ghost}",
                    shape.rows()
                ))
            }
            Query::Necklace(n) => {
                if n == 0 {
                    return Err(Error::plain(
                        "Witt spectrum necklace index must be positive",
                    ));
                }
                let elementary = |coefficients: &[i128]| match coefficients {
                    [1] => Some(0),
                    [1, middle @ .., last] if middle.iter().all(|&c| c == 0) => {
                        let degree = (coefficients.len() - 1) as u128;
                        match last {
                            -1 => Some(i128::from(n == degree)),
                            1 => Some(i128::from(n == 2 * degree) - i128::from(n == degree)),
                            _ => None,
                        }
                    }
                    _ => None,
                };
                let value = match (elementary(&self.numerator), elementary(&self.denominator)) {
                    (Some(a), Some(b)) => a - b,
                    _ => primitive_coordinate(model, n, |model, k| self.ghost(model, k))?
                        .signed()
                        .ok_or_else(coefficient_overflow)?,
                };
                Ok(format!("rvWittSpectrumNecklace[{n}] = {value}"))
            }
        }
    }

    fn ghost(&self, model: &mut Model, k: u128) -> Result<Integer> {
        let pos = power_trace(model, &self.numerator, k)?;
        let neg = power_trace(model, &self.denominator, k)?;
        add(pos, neg.negate())
    }

    fn canonical(model: &mut Model, mut p: Polynomial, mut q: Polynomial) -> Result<Self> {
        // Reciprocal polynomials are monic. Their primitive gcd over Q
        // is monic integral, so it divides both by exact monic division.
        if p.len() > 1 && q.len() > 1 {
            p.reverse();
            q.reverse();
            let factor = polynomial_gcd(model, p.clone(), q.clone())?;
            if factor.len() > 1 {
                let (pn, pr) = divide_monic(model, p, &factor)?;
                let (qn, qr) = divide_monic(model, q, &factor)?;
                if !pr.is_empty() || !qr.is_empty() {
                    return Err(integrity());
                }
                (p, q) = (pn, qn);
            }
            p.reverse();
            q.reverse();
        }
        require_degree(p.len() - 1)?;
        require_degree(q.len() - 1)?;
        Ok(Self {
            numerator: narrow(p)?,
            denominator: narrow(q)?,
        })
    }
}

fn tensor_factors(model: &mut Model, p: &[i128], q: &[i128]) -> Result<Polynomial> {
    let rank = (p.len() - 1) * (q.len() - 1);
    require_degree(rank)?;
    if rank == 0 {
        return Ok(vec![one()]);
    }
    let a = power_sums(model, p, rank)?;
    let b = power_sums(model, q, rank)?;
    let ghosts = a
        .into_iter()
        .zip(b)
        .skip(1)
        .map(|(left, right)| multiply(model, left, right))
        .collect::<Result<Vec<_>>>()?;
    determinant_from_ghosts(model, &ghosts)
}

fn adams_factor(model: &mut Model, p: &[i128], r: u128) -> Result<Polynomial> {
    let rank = p.len() - 1;
    if rank == 0 || r == 1 {
        return Ok(lift(p));
    }
    let characteristic: Polynomial = p.iter().rev().copied().map(Integer::from_i128).collect();
    let root = power_mod(model, &characteristic, r)?;
    let initial = power_sums(model, p, rank - 1)?;
    let mut power = root.clone();
    let mut ghosts = Vec::with_capacity(rank);
    for k in 1..=rank {
        ghosts.push(trace_residue(model, &power, &initial)?);
        if k < rank {
            power = multiply_mod(model, &power, &root, &characteristic)?;
        }
    }
    determinant_from_ghosts(model, &ghosts)
}

fn power_trace(model: &mut Model, p: &[i128], k: u128) -> Result<Integer> {
    let rank = p.len() - 1;
    if rank == 0 {
        return Ok(zero());
    }
    if k <= rank as u128 {
        return Ok(*power_sums(model, p, k as usize)?.last().unwrap());
    }
    let characteristic: Polynomial = p.iter().rev().copied().map(Integer::from_i128).collect();
    let power = power_mod(model, &characteristic, k)?;
    let initial = power_sums(model, p, rank - 1)?;
    trace_residue(model, &power, &initial)
}

fn trace_at_product(model: &mut Model, p: &[i128], a: u128, b: u128) -> Result<Integer> {
    if let Some(exponent) = a.checked_mul(b) {
        return power_trace(model, p, exponent);
    }
    let rank = p.len() - 1;
    if rank == 0 {
        return Ok(zero());
    }
    let characteristic: Polynomial = p.iter().rev().copied().map(Integer::from_i128).collect();
    let root = power_mod(model, &characteristic, a)?;
    let residue = residue_power(model, &characteristic, root, b)?;
    let sums = power_sums(model, p, rank - 1)?;
    trace_residue(model, &residue, &sums)
}

fn trace_residue(model: &mut Model, residue: &[Integer], sums: &[Integer]) -> Result<Integer> {
    residue.iter().zip(sums).try_fold(zero(), |sum, (&a, &b)| {
        let term = multiply(model, a, b)?;
        add(sum, term)
    })
}

fn power_sums(model: &mut Model, p: &[i128], through: usize) -> Result<Polynomial> {
    model.spend(through + 1)?;
    let mut sums = Vec::with_capacity(through + 1);
    sums.push(Integer::from_i128((p.len() - 1) as i128));
    for k in 1..=through {
        let mut result = zero();
        for (j, &coefficient) in p.iter().enumerate().skip(1).take(k) {
            let factor = if j == k {
                Integer::from_i128(k as i128)
            } else {
                sums[k - j]
            };
            let term = multiply(model, factor, Integer::from_i128(coefficient))?;
            result = add(result, term.negate())?;
        }
        sums.push(result);
    }
    Ok(sums)
}

pub(crate) fn determinant_from_ghosts(model: &mut Model, ghosts: &[Integer]) -> Result<Polynomial> {
    model.spend(ghosts.len() + 1)?;
    let mut coefficients = Vec::with_capacity(ghosts.len() + 1);
    coefficients.push(one());
    for k in 1..=ghosts.len() {
        let mut numerator = zero();
        for (&trace, &coefficient) in ghosts.iter().zip(coefficients.iter().rev()) {
            let term = multiply(model, trace, coefficient)?;
            numerator = add(numerator, term)?;
        }
        let (magnitude, remainder) = numerator.magnitude.divide(model, k as u128)?;
        if remainder != 0 {
            return Err(Error::plain("Witt spectrum Newton division is inexact"));
        }
        coefficients.push(Integer {
            magnitude,
            negative: !numerator.negative,
        });
    }
    Ok(coefficients)
}

fn power_mod(model: &mut Model, characteristic: &[Integer], exponent: u128) -> Result<Polynomial> {
    let root = divide_monic(model, vec![zero(), one()], characteristic)?.1;
    residue_power(model, characteristic, root, exponent)
}

fn residue_power(
    model: &mut Model,
    characteristic: &[Integer],
    mut root: Polynomial,
    mut exponent: u128,
) -> Result<Polynomial> {
    let mut result = vec![one()];
    while exponent != 0 {
        if !exponent.is_multiple_of(2) {
            result = multiply_mod(model, &result, &root, characteristic)?;
        }
        exponent >>= 1;
        if exponent != 0 {
            root = multiply_mod(model, &root, &root, characteristic)?;
        }
    }
    Ok(result)
}

fn multiply_mod(
    model: &mut Model,
    left: &[Integer],
    right: &[Integer],
    divisor: &[Integer],
) -> Result<Polynomial> {
    let value = product(model, left, right)?;
    Ok(divide_monic(model, value, divisor)?.1)
}

fn divide_monic(
    model: &mut Model,
    mut remainder: Polynomial,
    divisor: &[Integer],
) -> Result<(Polynomial, Polynomial)> {
    let mut quotient = vec![zero(); (remainder.len() + 1).saturating_sub(divisor.len())];
    for (shift, coefficient) in quotient.iter_mut().enumerate().rev() {
        *coefficient = remainder[shift + divisor.len() - 1];
        subtract_multiple(model, &mut remainder[shift..], divisor, *coefficient)?;
    }
    trim(&mut remainder);
    Ok((quotient, remainder))
}

fn polynomial_gcd(model: &mut Model, mut a: Polynomial, mut b: Polynomial) -> Result<Polynomial> {
    while b.len() > 1 {
        let leading = *b.last().unwrap();
        while a.len() >= b.len() {
            let top = *a.last().unwrap();
            let common = top.magnitude.gcd(model, leading.magnitude)?;
            let scale = exact_divide(model, leading, common)?;
            let multiple = exact_divide(model, top, common)?;
            for value in &mut a {
                *value = multiply(model, *value, scale)?;
            }
            let shift = a.len() - b.len();
            subtract_multiple(model, &mut a[shift..], &b, multiple)?;
            primitive(model, &mut a)?;
        }
        (a, b) = (b, a);
    }
    if !b.is_empty() {
        return Ok(vec![one()]);
    }
    if a.last().and_then(|c| c.signed()) != Some(1) {
        return Err(integrity());
    }
    Ok(a)
}

fn primitive(model: &mut Model, value: &mut Polynomial) -> Result<()> {
    trim(value);
    let Some(leading) = value.last() else {
        return Ok(());
    };
    let negative = leading.negative;
    let mut content = leading.magnitude;
    for coefficient in value.iter().rev().skip(1) {
        if content.small() == Some(1) {
            break;
        }
        content = content.gcd(model, coefficient.magnitude)?;
    }
    if content.small() != Some(1) || negative {
        for coefficient in value {
            let result = exact_divide(model, *coefficient, content)?;
            *coefficient = if negative { result.negate() } else { result };
        }
    }
    Ok(())
}

fn exact_divide(model: &mut Model, value: Integer, divisor: Magnitude) -> Result<Integer> {
    let magnitude = value
        .magnitude
        .divide_exact(model, divisor)?
        .ok_or_else(integrity)?;
    Ok(Integer {
        magnitude,
        negative: value.negative,
    })
}

fn subtract_multiple(
    model: &mut Model,
    into: &mut [Integer],
    source: &[Integer],
    multiple: Integer,
) -> Result<()> {
    for (target, &value) in into.iter_mut().zip(source) {
        let term = multiply(model, value, multiple)?;
        *target = add(*target, term.negate())?;
    }
    Ok(())
}

fn product(model: &mut Model, a: &[Integer], b: &[Integer]) -> Result<Polynomial> {
    if a.is_empty() || b.is_empty() {
        return Ok(Vec::new());
    }
    let width = a.len() + b.len() - 1;
    model.spend(width)?;
    let mut result = vec![zero(); width];
    for (i, &a) in a.iter().enumerate().filter(|(_, a)| !a.is_zero()) {
        for (target, &b) in result[i..].iter_mut().zip(b) {
            let term = multiply(model, a, b)?;
            *target = add(*target, term)?;
        }
    }
    trim(&mut result);
    Ok(result)
}

fn trim(polynomial: &mut Polynomial) {
    while polynomial.last().is_some_and(|c| c.is_zero()) {
        polynomial.pop();
    }
}

fn lift(p: &[i128]) -> Polynomial {
    p.iter().copied().map(Integer::from_i128).collect()
}

fn narrow(p: Polynomial) -> Result<Box<[i128]>> {
    p.into_iter()
        .map(|c| c.signed().ok_or_else(coefficient_overflow))
        .collect()
}

fn multiply(model: &mut Model, a: Integer, b: Integer) -> Result<Integer> {
    model.spend(1 + a.magnitude.bits().div_ceil(32) * b.magnitude.bits().div_ceil(32))?;
    a.multiply_wide(b).map_err(|_| overflow())
}

fn add(a: Integer, b: Integer) -> Result<Integer> {
    a.add(b).map_err(|_| overflow())
}

pub(crate) fn require_degree(degree: usize) -> Result<()> {
    if degree > MAX_DEGREE {
        return Err(degree_limit());
    }
    Ok(())
}

fn require_rank(rank: Integer) -> Result<usize> {
    if rank.negative && !rank.is_zero() {
        return Err(Error::plain("Witt spectrum parity rank is negative"));
    }
    let degree = rank
        .magnitude
        .small()
        .and_then(|value| usize::try_from(value).ok())
        .ok_or_else(degree_limit)?;
    require_degree(degree)?;
    Ok(degree)
}

fn degree_limit() -> Error {
    Error::plain("Witt spectrum exceeds the determinant degree limit")
}

fn zero() -> Integer {
    Integer::from_i128(0)
}
fn one() -> Integer {
    Integer::from_i128(1)
}

fn integrity() -> Error {
    Error::plain("Witt spectrum polynomial normalization failed an exactness check")
}

fn overflow() -> Error {
    Error::plain("Witt spectrum exact intermediate exceeds 512 bits")
}

fn coefficient_overflow() -> Error {
    Error::plain("Witt spectrum coefficient or ghost exceeds signed 128-bit range")
}
