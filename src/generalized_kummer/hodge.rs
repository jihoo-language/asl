use std::mem::size_of;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::hilbert::{
    checked_product, intermediate_overflow, unsigned, Integer, SurfaceKind,
};
use crate::model::Model;

use super::cohomology::{exact_order_sums, torsion_characters, torsion_divisors};
use super::require_translation_order;

pub(crate) struct Command(pub Operation, pub u128, pub u128, pub u128);

pub(crate) enum Operation {
    K3Hilbert,
    JacobianHilbert,
    Generalized,
    Translate(u128),
    Invert(u128),
    Spectrum(u128),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<String> {
        let Self(operation, n, p, q) = self;
        match operation {
            Operation::K3Hilbert => {
                let value = hilbert_hodge(model, SurfaceKind::K3, n, p, q)?;
                Ok(format!("rvK3HilbHodge[{n}, {p}, {q}] = {value}"))
            }
            Operation::JacobianHilbert => {
                let value = hilbert_hodge(model, SurfaceKind::Abelian, n, p, q)?;
                Ok(format!("rvJacHilbHodge[{n}, {p}, {q}] = {value}"))
            }
            Operation::Generalized => {
                require_translation_order(n, 1)?;
                let value = Character::new(model, n, p, q, false)?.evaluate(model, 1)?;
                let count = hodge_dimension(value)?;
                Ok(format!("rvGenKummerHodge[{n}, {p}, {q}] = {count}"))
            }
            Operation::Translate(h) => {
                require_translation_order(n, h)?;
                let value = Character::new(model, n, p, q, false)?.evaluate(model, h)?;
                let trace = value
                    .signed()
                    .ok_or_else(|| Error::plain("Hodge character exceeds signed 128-bit range"))?;
                Ok(format!(
                    "rvGenKummerTranslateHodge[{n}, {h}, {p}, {q}] = {trace}"
                ))
            }
            Operation::Invert(h) => {
                require_translation_order(n, h)?;
                let value = Character::new(model, n, p, q, true)?.evaluate(model, h)?;
                let trace = remove_parity_sign(value, p, q)
                    .signed()
                    .ok_or_else(|| Error::plain("Hodge character exceeds signed 128-bit range"))?;
                Ok(format!(
                    "rvGenKummerInversionHodge[{n}, {h}, {p}, {q}] = {trace}"
                ))
            }
            Operation::Spectrum(h) => {
                require_translation_order(n, h)?;
                let spectrum = Character::new(model, n, p, q, false)?.spectrum(model, h)?;
                Ok(format!(
                    "rvGenKummerHodgeSpectrum[{n}, {h}, {p}, {q}] = {spectrum:?}"
                ))
            }
        }
    }
}

const ABELIAN: [[usize; 3]; 3] = [[1, 2, 1], [2, 4, 2], [1, 2, 1]];
const K3: [[usize; 3]; 3] = [[1, 0, 1], [0, 20, 0], [1, 0, 1]];

fn hilbert_hodge(model: &mut Model, kind: SurfaceKind, n: u128, p: u128, q: u128) -> Result<u128> {
    let Some((p, q)) = hodge_indices(n, p, q) else {
        return Ok(0);
    };
    let surface = match kind {
        SurfaceKind::K3 => {
            if odd(p, q) {
                return Ok(0);
            }
            if p == 0 {
                return Ok(1);
            }
            &K3
        }
        SurfaceKind::Abelian => {
            if p == 0 {
                // Use only the length-one holomorphic generators: one even
                // (0,2) class and two independent odd (0,1) classes.
                return Ok(if q.is_multiple_of(2) {
                    1 + u128::from(q != 0 && q / 2 < n)
                } else {
                    2
                });
            }
            &ABELIAN
        }
    };
    let table = Fock::new(model, n, p, q, surface)?;
    hodge_dimension(remove_parity_sign(table.at(n, table.p, table.q), p, q))
}

// Hodge conjugation interchanges p and q. Poincaré duality reflects both
// indices, taking the inverse automorphism, which has the same order.
fn hodge_indices(length: u128, mut p: u128, mut q: u128) -> Option<(u128, u128)> {
    if p.div_ceil(2) > length || q.div_ceil(2) > length {
        return None;
    }
    let half_sum = p / 2 + q / 2;
    if half_sum > length || (half_sum == length && (p | q) & 1 != 0) {
        // Here 2*length is less than p+q. Both complements fit u128,
        // even when the dimension or p+q itself does not.
        let complement = |index| {
            if index > length {
                length - (index - length)
            } else {
                length + (length - index)
            }
        };
        (p, q) = (complement(p), complement(q));
    }
    Some((p.min(q), p.max(q)))
}

struct Fock {
    length: usize,
    p: usize,
    q: usize,
    stride: usize,
    coefficients: Vec<Integer>,
}

impl Fock {
    fn new(
        model: &mut Model,
        n: u128,
        p: u128,
        q: u128,
        surface: &[[usize; 3]; 3],
    ) -> Result<Self> {
        // Apart from the length-one (0,0) class, every generator has
        // total Hodge degree at least its length.
        let length = usize::try_from(n.min(p.saturating_add(q))).map_err(|_| storage_limit())?;
        let p = usize::try_from(p).map_err(|_| storage_limit())?;
        let q = usize::try_from(q).map_err(|_| storage_limit())?;
        let stride = p.saturating_add(1).saturating_mul(q.saturating_add(1));
        let cells = length.saturating_add(1).saturating_mul(stride);
        if cells.saturating_mul(size_of::<Integer>().div_ceil(size_of::<usize>())) > MAX_PIVOT_CELLS
        {
            return Err(storage_limit());
        }
        model.spend(cells)?;
        let mut table = Self {
            length,
            p,
            q,
            stride,
            coefficients: vec![Integer::from_i128(0); cells],
        };
        table.coefficients[0] = Integer::from_i128(1);
        for m in 1..=length.min(p.min(q) + 1) {
            let shift = m - 1;
            for (a, row) in surface.iter().enumerate() {
                for (b, &dimension) in row.iter().enumerate() {
                    if a + shift > p || b + shift > q || (m == 1 && a + b == 0) {
                        continue;
                    }
                    for _ in 0..dimension {
                        table.add_generator(model, m, a + shift, b + shift, (a + b) % 2 != 0)?;
                    }
                }
            }
        }
        // Restore the length-one (0,0) symmetric factor by prefix summation
        // in length, rather than taking a row of arbitrarily large length.
        model.spend(cells - stride)?;
        for i in stride..cells {
            let previous = table.coefficients[i - stride];
            add(&mut table.coefficients[i], previous)?;
        }
        Ok(table)
    }

    fn add_generator(
        &mut self,
        model: &mut Model,
        m: usize,
        a: usize,
        b: usize,
        exterior: bool,
    ) -> Result<()> {
        let last = self.length - m;
        model.spend((last + 1) * (self.p + 1 - a) * (self.q + 1 - b))?;
        let width = self.q + 1;
        let shift = m * self.stride + a * width + b;
        for index in 0..=last {
            let length = if exterior { last - index } else { index };
            for p in 0..=self.p - a {
                let start = length * self.stride + p * width;
                for from in start..=start + self.q - b {
                    let source = self.coefficients[from];
                    let value = if exterior { source.negate() } else { source };
                    add(&mut self.coefficients[from + shift], value)?;
                }
            }
        }
        Ok(())
    }

    fn at(&self, n: u128, p: usize, q: usize) -> Integer {
        let row = n.min(self.length as u128) as usize;
        self.coefficients[row * self.stride + p * (self.q + 1) + q]
    }

    fn divided(&self, model: &mut Model, n: u128, p: usize, q: usize) -> Result<Integer> {
        // The inverse of E(A)=(1-x)^2(1-y)^2 has coefficient (p+1)(q+1).
        model.spend((p + 1) * (q + 1))?;
        let mut result = Integer::from_i128(0);
        for a in 0..=p {
            for b in 0..=q {
                let source = self.at(n, a, b);
                if !source.is_zero() {
                    let factor = ((p + 1 - a) * (q + 1 - b)) as i128;
                    let term = checked_product(model, source, Integer::from_i128(factor))?;
                    add(&mut result, term)?;
                }
            }
        }
        Ok(result)
    }
}

struct Character {
    n: u128,
    divisors: Vec<usize>,
    contributions: Vec<Integer>,
}

impl Character {
    fn new(model: &mut Model, n: u128, p: u128, q: u128, inversion: bool) -> Result<Self> {
        let mut divisors = vec![1];
        let contributions = match hodge_indices(n - 1, p, q) {
            None => vec![Integer::from_i128(0)],
            // Holomorphic forms are generated by the translation-invariant
            // symplectic form, with one power in each even degree.
            Some((0, q)) => vec![Integer::from_i128(i128::from(q.is_multiple_of(2)))],
            Some((p, q)) => {
                let table = Fock::new(model, n, p, q, &ABELIAN)?;
                if inversion {
                    // Inversion interchanges χ and χ^{-1}. Only the
                    // character of order 1 and the order-2 terms contribute.
                    if n.is_multiple_of(2) && p >= n / 2 {
                        divisors.push(2);
                    }
                } else if p >= n.div_ceil(2) {
                    // Each nontrivial divisor shifts both indices by at least
                    // n/2. If any term can contribute, n fits the table bound.
                    divisors = torsion_divisors(model, n as usize, table.p)?;
                }
                divisors
                    .iter()
                    .map(|&d| {
                        let length = n / d as u128;
                        let shift = (n - length) as usize;
                        let value =
                            table.divided(model, length, table.p - shift, table.q - shift)?;
                        Ok(remove_parity_sign(value, p, q))
                    })
                    .collect::<Result<_>>()?
            }
        };
        Ok(Self {
            n,
            divisors,
            contributions,
        })
    }

    fn evaluate(&self, model: &mut Model, h: u128) -> Result<Integer> {
        if self.divisors.len() == 1 {
            return Ok(self.contributions[0]);
        }
        let counts = torsion_characters(model, self.n, h, &self.divisors)?;
        let mut trace = Integer::from_i128(0);
        for (count, &contribution) in counts.into_iter().zip(&self.contributions) {
            let term = checked_product(model, Integer::from_i128(count), contribution)?;
            add(&mut trace, term)?;
        }
        Ok(trace)
    }

    fn spectrum(&self, model: &mut Model, h: u128) -> Result<Vec<(usize, u128)>> {
        if h == 1 || self.divisors.len() == 1 {
            let dimension = hodge_dimension(self.evaluate(model, 1)?)?;
            return Ok((dimension != 0)
                .then_some((1, dimension))
                .into_iter()
                .collect());
        }
        // This branch requires a nontrivial divisor contribution, which
        // bounds h <= n by the coefficient table's length.
        let h = h as usize;
        model.spend(h)?;
        let orders: Vec<_> = (1..=h).filter(|d| h.is_multiple_of(*d)).collect();
        let mut traces = Vec::with_capacity(orders.len());
        for &order in &orders {
            traces.push(self.evaluate(model, order as u128)?);
        }
        let dimension = hodge_dimension(traces[0])?;
        let totients = exact_order_sums(model, &orders, |e| Ok(e as i128))?;
        let mut spectrum = Vec::new();
        let mut rank = 0u128;
        for (&d, phi) in orders.iter().zip(totients) {
            // Sums of powers of primitive u-th roots are obtained by
            // subtracting proper orders from the sum over all e-th roots.
            let power = h / d;
            let ramanujan = exact_order_sums(model, &orders, |e| {
                Ok(if power.is_multiple_of(e) {
                    e as i128
                } else {
                    0
                })
            })?;
            let mut sum = Integer::from_i128(0);
            for (weight, &trace) in ramanujan.into_iter().zip(&traces) {
                let term = checked_product(model, Integer::from_i128(weight), trace)?;
                add(&mut sum, term)?;
            }
            let (quotient, remainder) = sum.magnitude.divide(model, h as u128)?;
            if (sum.negative && !sum.is_zero()) || remainder != 0 {
                return Err(spectrum_integrity());
            }
            let count = quotient.small().ok_or_else(dimension_overflow)?;
            if count != 0 {
                let term = count
                    .checked_mul(phi as u128)
                    .ok_or_else(spectrum_integrity)?;
                rank = rank.checked_add(term).ok_or_else(spectrum_integrity)?;
                spectrum.push((d, count));
            }
        }
        if rank != dimension {
            return Err(spectrum_integrity());
        }
        Ok(spectrum)
    }
}

fn odd(p: u128, q: u128) -> bool {
    (p ^ q) & 1 != 0
}

fn remove_parity_sign(value: Integer, p: u128, q: u128) -> Integer {
    if odd(p, q) {
        value.negate()
    } else {
        value
    }
}

fn add(target: &mut Integer, term: Integer) -> Result<()> {
    if !term.is_zero() {
        *target = target.add(term).map_err(|_| intermediate_overflow())?;
    }
    Ok(())
}

fn hodge_dimension(value: Integer) -> Result<u128> {
    if value.negative && !value.is_zero() {
        return Err(Error::plain("Hodge number must be nonnegative"));
    }
    unsigned(value).ok_or_else(dimension_overflow)
}

fn dimension_overflow() -> Error {
    Error::plain("Hodge number exceeds unsigned 128-bit range")
}

fn spectrum_integrity() -> Error {
    Error::plain("Hodge character failed the cyclotomic multiplicity check")
}

fn storage_limit() -> Error {
    Error::plain("Hodge coefficient table exceeds the storage limit")
}
