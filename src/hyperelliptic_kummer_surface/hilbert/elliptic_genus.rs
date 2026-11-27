use std::collections::{btree_map::Entry, BTreeMap};
use std::fmt::Write as _;
use std::mem::size_of;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;

use super::Integer;

type Seed = BTreeMap<(usize, i128), Integer>;
type Weight = (usize, usize, i128);
type Product = BTreeMap<Weight, Integer>;

const INTEGER_WORDS: usize = size_of::<Integer>().div_ceil(size_of::<usize>());
const TERM_WORDS: usize = size_of::<(Weight, Integer)>().div_ceil(size_of::<usize>()) + 4;
const TERM_LIMIT: usize = MAX_PIVOT_CELLS / TERM_WORDS;

pub(crate) enum Query {
    EllipticGenus(u128, u128, i128),
    Igusa(u128, u128, i128),
}

impl Query {
    pub(super) fn report(self, model: &mut Model) -> Result<String> {
        let (name, p, q, r, coefficient) = match self {
            Self::EllipticGenus(length, degree, r) => (
                "rvK3EllipticGenus",
                length,
                degree,
                r,
                dmvv(model, length, degree, r)?,
            ),
            Self::Igusa(m, n, r) => ("rvIgusaPhi10", m, n, r, igusa(model, m, n, r)?),
        };
        Ok(format!("{name}[{p}, {q}, {r}] = {}", decimal(coefficient)))
    }
}

fn dmvv(model: &mut Model, length: u128, degree: u128, r: i128) -> Result<Integer> {
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
    let product = series(model, cutoff, true)?;
    Ok(coefficient(&product, (cutoff.p, cutoff.q, r)))
}

fn igusa(model: &mut Model, m: u128, n: u128, r: i128) -> Result<Integer> {
    if m == 0 || n == 0 {
        return Ok(Integer::from_i128(0));
    }
    let cutoff = Cutoff::new(model, m - 1, n - 1)?;
    let product = series(model, cutoff, false)?;
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

fn series(model: &mut Model, cutoff: Cutoff, reciprocal: bool) -> Result<Product> {
    let mut seed = Seed::new();
    let primitive_degree = cutoff.p * cutoff.q;
    add_theta(model, &mut seed, primitive_degree, false)?;
    add_theta(model, &mut seed, primitive_degree, true)?;
    let mut product = Product::from([((0, 0, 0), Integer::from_i128(1))]);
    for p in usize::from(reciprocal)..=cutoff.p {
        for q in 0..=cutoff.q {
            if p == 0 && q == 0 {
                continue;
            }
            // The Jacobi bound r² ≤ 4pq+1 follows directly from the two
            // theta numerators. Retain every resulting Laurent exponent.
            let degree = p * q;
            for (&(_, r), &c) in seed.range((degree, i128::MIN)..=(degree, i128::MAX)) {
                let exponent = if reciprocal { c.negate() } else { c };
                multiply_factor(model, &mut product, cutoff, (p, q, r), exponent)?;
            }
        }
    }
    Ok(product)
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
) -> Result<()> {
    // Generalized binomial recurrence for (1 - X)^e.
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
        factors.push(binomial);
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

fn decimal(value: Integer) -> String {
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
