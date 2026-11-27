use std::collections::{btree_map::Entry, BTreeMap};
use std::fmt::{self, Write as _};
use std::mem::size_of;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::siegel::Weight as SiegelWeight;
use crate::syntax::Label;

use super::{Integer, Magnitude};
use crate::hyperelliptic_kummer_surface::{enriques, nikulin::basis_for};

mod cyclic;
pub(crate) mod eigenpacket;
pub(crate) mod farey;
pub(crate) mod hecke;
pub(crate) mod hurwitz;
pub(crate) mod inverse;
pub(crate) mod maass;
mod mock;
pub(crate) mod period_polynomial;

mod generalized;

type Seed = BTreeMap<(usize, i128), Integer>;
type Weight = (usize, usize, i128);
type Product = BTreeMap<Weight, Integer>;

const INTEGER_WORDS: usize = size_of::<Integer>().div_ceil(size_of::<usize>());
const TERM_WORDS: usize = size_of::<(Weight, Integer)>().div_ceil(size_of::<usize>()) + 4;
const TERM_LIMIT: usize = MAX_PIVOT_CELLS / TERM_WORDS;

pub(crate) enum Query {
    EllipticGenus(Index),
    Generalized(Index),
    Igusa(Index),
    IgusaInverse(inverse::Query),
    Delta5(SiegelWeight),
    Symplectic(usize, Index),
    Translate(Label, usize, Index),
    Enriques(Label, Label, usize, Index),
}

#[derive(Clone, Copy)]
pub(crate) struct Index(pub(crate) u128, pub(crate) u128, pub(crate) i128);

impl fmt::Display for Index {
    fn fmt(&self, output: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self(p, q, r) = *self;
        write!(output, "{p}, {q}, {r}")
    }
}

#[derive(Clone, Copy)]
enum Expansion {
    Dmvv,
    Symplectic,
    HalfDmvv,
    Cyclic(usize),
    Igusa,
    IgusaInverse,
    Delta5,
}

impl Query {
    pub(super) fn report(self, model: &mut Model) -> Result<String> {
        let (name, coefficient) = match self {
            Self::EllipticGenus(index) => (
                format!("rvK3EllipticGenus[{index}]"),
                dmvv(model, index, Expansion::Dmvv)?,
            ),
            Self::Generalized(index @ Index(length, degree, r)) => (
                format!("rvGenKummerElliptic[{index}]"),
                generalized::coefficient(model, length, degree, r)?,
            ),
            Self::Igusa(index) => (format!("rvIgusaPhi10[{index}]"), igusa(model, index)?),
            Self::IgusaInverse(query) => return query.report(model),
            Self::Delta5(weight) => return delta5_report(model, weight),
            Self::Symplectic(order, index) => {
                let expansion = if order == 1 {
                    Expansion::Dmvv
                } else {
                    Expansion::Cyclic(order)
                };
                (
                    format!("rvK3SymplecticElliptic[{order}, {index}]"),
                    dmvv(model, index, expansion)?,
                )
            }
            Self::Translate(basis, mask, index) => {
                basis_for(
                    model,
                    &basis,
                    "Kummer translation elliptic genera require a certified genus-two level-2 basis",
                )?;
                let expansion = if mask == 0 {
                    Expansion::Dmvv
                } else {
                    Expansion::Symplectic
                };
                let value = dmvv(model, index, expansion)?;
                (
                    format!("rvKummerTranslateElliptic[{}, {mask}, {index}]", basis.text),
                    value,
                )
            }
            Self::Enriques(basis, base, mask, index) => {
                enriques::Action::resolve(model, &basis, &base, mask)?;
                let value = free_involution(model, index)?;
                (
                    format!(
                        "rvKummerEnriquesElliptic[{}, {}, {mask}, {index}]",
                        basis.text, base.text
                    ),
                    value,
                )
            }
        };
        Ok(format!("{name} = {}", decimal(coefficient)))
    }
}

fn dmvv(
    model: &mut Model,
    Index(length, degree, r): Index,
    expansion: Expansion,
) -> Result<Integer> {
    if length == 0 {
        return Ok(Integer::from_i128(i128::from(degree == 0 && r == 0)));
    }
    let bound = length
        .checked_mul(degree)
        .and_then(|nd| nd.checked_mul(4))
        .and_then(|nd| length.checked_mul(length).and_then(|nn| nd.checked_add(nn)));
    if bound.is_some_and(|square| r.unsigned_abs() > square.isqrt()) {
        return Ok(Integer::from_i128(0));
    }
    let cutoff = Cutoff::new(model, length, degree)?;
    let product = series(model, cutoff, expansion)?;
    Ok(coefficient(&product, (cutoff.p, cutoff.q, r)))
}

fn free_involution(model: &mut Model, Index(length, degree, r): Index) -> Result<Integer> {
    // The primitive trace vanishes. Both eigenspaces have multiplicity c/2,
    // so the two factors combine as (1 - X²)^(-c/2).
    if !length.is_multiple_of(2) || !degree.is_multiple_of(2) || r % 2 != 0 {
        return Ok(Integer::from_i128(0));
    }
    dmvv(
        model,
        Index(length / 2, degree / 2, r / 2),
        Expansion::HalfDmvv,
    )
}

fn igusa(model: &mut Model, Index(m, n, r): Index) -> Result<Integer> {
    if m == 0 || n == 0 {
        return Ok(Integer::from_i128(0));
    }
    let cutoff = Cutoff::new(model, m - 1, n - 1)?;
    let product = series(model, cutoff, Expansion::Igusa)?;
    // Multiply by the Weyl factor P q (y - 2 + y⁻¹).
    let mut answer = Integer::from_i128(0);
    for (shift, scale) in [(-1, 1), (0, -2), (1, 1)] {
        if let Some(index) = r.checked_add(shift) {
            let term = coefficient(&product, (cutoff.p, cutoff.q, index));
            answer = answer
                .add(term.multiply(scale).map_err(|_| overflow())?)
                .map_err(|_| overflow())?;
        }
    }
    Ok(answer)
}

fn delta5(model: &mut Model, SiegelWeight([a, b, r]): SiegelWeight) -> Result<Integer> {
    if a < 1 || b < 1 || [a, b, r].iter().any(|coordinate| coordinate % 2 == 0) {
        return Ok(Integer::from_i128(0));
    }
    // Every product root apart from the isolated axis factor in the
    // prefactor satisfies |r| ≤ a+b in ordinary coordinates.
    if r.unsigned_abs() >= (a as u128) + (b as u128) {
        return Ok(Integer::from_i128(0));
    }
    let cutoff = Cutoff::new(model, (a / 2) as u128, (b / 2) as u128)?;
    let product = series(model, cutoff, Expansion::Delta5)?;
    // Multiply by (Pq)^(1/2) (y^(1/2) - y^(-1/2)).
    let r = r.div_euclid(2);
    let upper = coefficient(&product, (cutoff.p, cutoff.q, r));
    let lower = coefficient(&product, (cutoff.p, cutoff.q, r + 1));
    upper.add(lower.negate()).map_err(|_| overflow())
}

pub(crate) fn delta5_report(model: &mut Model, weight: SiegelWeight) -> Result<String> {
    Ok(format!(
        "rvSiegelDelta5[{weight}] = {}",
        decimal(delta5(model, weight)?)
    ))
}

pub(crate) fn half_k3_coefficient(model: &mut Model, degree: u128, r: i128) -> Result<Integer> {
    let cutoff = Cutoff::new(model, degree, 1)?;
    let seed = jacobi_zero(model, cutoff.p)?;
    Ok(seed
        .get(&(cutoff.p, r))
        .copied()
        .unwrap_or(Integer::from_i128(0)))
}

fn jacobi_zero(model: &mut Model, degree: usize) -> Result<Seed> {
    let mut seed = Seed::new();
    add_theta(model, &mut seed, degree, false)?;
    add_theta(model, &mut seed, degree, true)?;
    halve_seed(model, &mut seed)?;
    Ok(seed)
}

fn jacobi_minus_two(model: &mut Model, cutoff: Cutoff) -> Result<Product> {
    let mut seed = Product::from([
        ((0, 0, -1), Integer::from_i128(1)),
        ((0, 0, 0), Integer::from_i128(-2)),
        ((0, 0, 1), Integer::from_i128(1)),
    ]);
    // (y - 2 + y⁻¹) ∏ (1 - qʲy)^2 (1 - qʲy⁻¹)^2 / (1 - qʲ)^4.
    for j in 1..=cutoff.q {
        for (r, exponent) in [(-1, 2), (0, -4), (1, 2)] {
            multiply_factor(
                model,
                &mut seed,
                cutoff,
                (0, j, r),
                Integer::from_i128(exponent),
                false,
            )?;
        }
    }
    Ok(seed)
}

#[derive(Clone, Copy)]
struct Cutoff {
    p: usize,
    q: usize,
}

impl Cutoff {
    fn new(model: &mut Model, p: u128, q: u128) -> Result<Self> {
        let p = usize::try_from(p).map_err(|_| storage())?;
        let q = usize::try_from(q).map_err(|_| storage())?;
        let slices = p.saturating_add(1).saturating_mul(q.saturating_add(1));
        let theta = p.saturating_mul(q).saturating_mul(2).saturating_add(1);
        if slices
            .saturating_mul(INTEGER_WORDS)
            .saturating_add(theta.saturating_mul(INTEGER_WORDS + 2))
            > MAX_PIVOT_CELLS
        {
            return Err(storage());
        }
        model.spend(slices)?;
        Ok(Self { p, q })
    }

    fn order(self, (p, q, _): Weight) -> usize {
        // Every factor has at least one positive series degree.
        self.p
            .checked_div(p)
            .unwrap_or(usize::MAX)
            .min(self.q.checked_div(q).unwrap_or(usize::MAX))
    }
}

fn series(model: &mut Model, cutoff: Cutoff, expansion: Expansion) -> Result<Product> {
    if let Expansion::Cyclic(order) = expansion {
        return cyclic::series(model, cutoff, order);
    }
    let mut positive = Seed::new();
    let primitive_degree = cutoff.p * cutoff.q;
    add_theta(model, &mut positive, primitive_degree, false)?;
    let negative = if matches!(expansion, Expansion::Symplectic) {
        // The θ₂ component is c_g. Thus c_-=(c-c_g)/2 is half
        // of the unshifted part and c_+=c_-+c_g adds θ₂.
        halve_seed(model, &mut positive)?;
        model.spend(positive.len() * INTEGER_WORDS)?;
        positive.clone()
    } else {
        Seed::new()
    };
    add_theta(model, &mut positive, primitive_degree, true)?;
    if matches!(expansion, Expansion::HalfDmvv | Expansion::Delta5) {
        halve_seed(model, &mut positive)?;
    }
    let reciprocal = !matches!(expansion, Expansion::Igusa | Expansion::Delta5);
    let mut product = Product::from([((0, 0, 0), Integer::from_i128(1))]);
    let first_p = usize::from(reciprocal && !matches!(expansion, Expansion::IgusaInverse));
    for p in first_p..=cutoff.p {
        for q in 0..=cutoff.q {
            if p == 0 && q == 0 {
                continue;
            }
            // The Jacobi bound r² ≤ 4pq+1 follows directly from the two
            // theta numerators. Retain every resulting Laurent exponent.
            let degree = p * q;
            for (negative_eigenvalue, seed) in [(false, &positive), (true, &negative)] {
                for (&(_, r), &c) in seed.range((degree, i128::MIN)..=(degree, i128::MAX)) {
                    let exponent = if reciprocal { c.negate() } else { c };
                    multiply_factor(
                        model,
                        &mut product,
                        cutoff,
                        (p, q, r),
                        exponent,
                        negative_eigenvalue,
                    )?;
                }
            }
        }
    }
    Ok(product)
}

fn eta_product(model: &mut Model, bound: usize, exponent: i128) -> Result<Vec<Integer>> {
    // Logarithmic derivative of ∏ (1-qʲ)^exponent.
    let mut sigma = vec![0i128; bound + 1];
    for d in 1..=bound {
        for n in (d..=bound).step_by(d) {
            model.spend(1)?;
            sigma[n] += d as i128;
        }
    }
    let mut series = Vec::with_capacity(bound + 1);
    series.push(Integer::from_i128(1));
    for n in 1..=bound {
        let mut total = Integer::from_i128(0);
        for (&coefficient, &scale) in series.iter().rev().zip(&sigma[1..=n]) {
            let term = multiply(model, coefficient, Integer::from_i128(-exponent * scale))?;
            total = total.add(term).map_err(|_| overflow())?;
        }
        let (magnitude, remainder) = total.magnitude.divide(model, n as u128)?;
        if remainder != 0 {
            return Err(Error::plain("eta product coefficient division is inexact"));
        }
        series.push(Integer { magnitude, ..total });
    }
    Ok(series)
}

fn halve_seed(model: &mut Model, seed: &mut Seed) -> Result<()> {
    for coefficient in seed.values_mut() {
        let (magnitude, remainder) = coefficient.magnitude.divide(model, 2)?;
        if remainder != 0 {
            return Err(Error::plain(
                "K3 elliptic genus coefficient does not halve exactly",
            ));
        }
        coefficient.magnitude = magnitude;
    }
    Ok(())
}

fn add_theta(model: &mut Model, seed: &mut Seed, degree: usize, shifted: bool) -> Result<()> {
    // θ₃² and θ₄² contribute 16 times the even powers of t=q^(1/2).
    // For θ₂² the common initial monomial cancels, leaving degree
    // k(k+1)/2 in q with a y shift of one after squaring.
    let (stride, y_shift, amplitude) = if shifted { (1, 1, 2) } else { (2, 0, 16) };
    let top = stride * degree;
    let mut numerator = Vec::new();
    for k in 0.. {
        let power = if shifted { k * (k + 1) / 2 } else { k * k };
        if power > top {
            break;
        }
        numerator.push((power, k as i128));
        if shifted || k != 0 {
            numerator.push((power, -(k as i128) - y_shift));
        }
    }
    let denominator = numerator.iter().filter(|&&(_, r)| !shifted || r >= 0);
    let mut squared = vec![0i128; top + 1];
    for &(a, _) in denominator.clone() {
        for &(b, _) in denominator.clone() {
            model.spend(1)?;
            if a + b <= top {
                squared[a + b] += 1;
            }
        }
    }
    // In θ₂ the denominator sums exactly the nonnegative k. Both squared
    // denominators have constant term one, so no rational division is needed.
    let denominator = squared
        .into_iter()
        .enumerate()
        .skip(1)
        .filter(|&(_, coefficient)| coefficient != 0)
        .collect::<Vec<_>>();
    let mut inverse = Vec::with_capacity(top + 1);
    inverse.push(Integer::from_i128(amplitude));
    for d in 1..=top {
        let mut sum = Integer::from_i128(0);
        for &(shift, weight) in denominator.iter().take_while(|&&(shift, _)| shift <= d) {
            let term = multiply(model, inverse[d - shift], Integer::from_i128(weight))?;
            sum = sum.add(term).map_err(|_| overflow())?;
        }
        model.spend(INTEGER_WORDS)?;
        inverse.push(sum.negate());
    }
    for &(a, k) in &numerator {
        for &(b, l) in &numerator {
            model.spend(1)?;
            let initial = a + b;
            for d in (initial.next_multiple_of(stride)..=top).step_by(stride) {
                add(
                    model,
                    seed,
                    (d / stride, k + l + y_shift),
                    inverse[d - initial],
                )?;
            }
        }
    }
    Ok(())
}

fn multiply_factor(
    model: &mut Model,
    product: &mut Product,
    cutoff: Cutoff,
    (p, q, r): Weight,
    exponent: Integer,
    negative_eigenvalue: bool,
) -> Result<()> {
    // Generalized binomial recurrence for (1 - X)^e. For eigenvalue
    // −1 replace X by −X, which negates each odd-degree coefficient.
    let mut binomial = Integer::from_i128(1);
    let mut factors = Vec::new();
    for index in 1..=cutoff.order((p, q, r)) {
        let numerator = Integer::from_i128(index as i128 - 1)
            .add(exponent.negate())
            .map_err(|_| overflow())?;
        if numerator.is_zero() {
            break;
        }
        let term = multiply(model, binomial, numerator)?;
        let (magnitude, remainder) = term.magnitude.divide(model, index as u128)?;
        if remainder != 0 {
            return Err(Error::plain("elliptic genus binomial division is inexact"));
        }
        binomial = Integer { magnitude, ..term };
        factors.push(if negative_eigenvalue && !index.is_multiple_of(2) {
            binomial.negate()
        } else {
            binomial
        });
    }
    // Snapshot only terms that can contribute. A positive P or q degree
    // ensures successive choices from a factor are all finite.
    model.spend(product.len())?;
    let inputs = product
        .iter()
        .filter(|&(&(a, b, _), _)| a + p <= cutoff.p && b + q <= cutoff.q)
        .map(|(&index, &coefficient)| (index, coefficient))
        .collect::<Vec<_>>();
    for ((a, b, s), value) in inputs {
        for (index, &factor) in factors.iter().enumerate() {
            let order = index + 1;
            let target = (a + order * p, b + order * q, s + order as i128 * r);
            if target.0 > cutoff.p || target.1 > cutoff.q {
                break;
            }
            let contribution = multiply(model, value, factor)?;
            add(model, product, target, contribution)?;
        }
    }
    Ok(())
}

fn coefficient(product: &Product, index: Weight) -> Integer {
    product
        .get(&index)
        .copied()
        .unwrap_or_else(|| Integer::from_i128(0))
}

fn add<K: Ord>(
    model: &mut Model,
    series: &mut BTreeMap<K, Integer>,
    index: K,
    value: Integer,
) -> Result<()> {
    model.spend(INTEGER_WORDS)?;
    if value.is_zero() {
        return Ok(());
    }
    let full = series.len() >= TERM_LIMIT;
    match series.entry(index) {
        Entry::Occupied(mut entry) => {
            let sum = entry.get().add(value).map_err(|_| overflow())?;
            if sum.is_zero() {
                entry.remove();
            } else {
                *entry.get_mut() = sum;
            }
        }
        Entry::Vacant(entry) => {
            if full {
                return Err(storage());
            }
            entry.insert(value);
        }
    }
    Ok(())
}

fn multiply(model: &mut Model, a: Integer, b: Integer) -> Result<Integer> {
    model.spend(1 + a.magnitude.bits().div_ceil(32) * b.magnitude.bits().div_ceil(32))?;
    a.multiply_wide(b).map_err(|_| overflow())
}

fn scale_by_power(
    model: &mut Model,
    mut value: Integer,
    base: u128,
    mut exponent: u128,
) -> Result<Integer> {
    if base == 1 {
        return Ok(value);
    }
    let mut power = Integer {
        magnitude: Magnitude::from_small(base),
        negative: false,
    };
    while exponent != 0 {
        if exponent % 2 == 1 {
            value = multiply(model, value, power)?;
        }
        exponent /= 2;
        if exponent != 0 {
            power = multiply(model, power, power)?;
        }
    }
    Ok(value)
}

pub(crate) fn decimal(value: Integer) -> String {
    if let Some(small) = value.signed() {
        return small.to_string();
    }
    let mut magnitude = value.magnitude;
    let mut digits = Vec::new();
    while !magnitude.is_zero() {
        let (quotient, remainder) = magnitude.divide_u32(1_000_000_000);
        digits.push(remainder);
        magnitude = quotient;
    }
    let mut result = if value.negative {
        "-".to_owned()
    } else {
        String::new()
    };
    let first = digits.pop().expect("nonzero wide integer");
    write!(result, "{first}").unwrap();
    for digit in digits.into_iter().rev() {
        write!(result, "{digit:09}").unwrap();
    }
    result
}

fn storage() -> Error {
    Error::plain("elliptic genus calculation exceeds the coefficient storage limit")
}

fn overflow() -> Error {
    Error::plain("elliptic genus exact intermediate exceeds 512 bits")
}
