use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::hw_closed::Curve;
use crate::hyperelliptic::Residues;
use crate::model::{Model, MAX_SHOWN_BYTES};
use crate::odd_cobar::Prime;
use crate::rv_ore::gcd;
use crate::rv_unram::{add_mod, mul_mod, sub_mod};
use crate::syntax::Label;

use super::{graded_index, node_orbits, storage_limit};

pub(crate) struct Query {
    pub(crate) curve: Label,
    pub(crate) length: u128,
    pub(crate) extension: u128,
    pub(crate) degree: u128,
    pub(crate) prime: usize,
    pub(crate) order: usize,
}

impl Query {
    pub(crate) fn report(&self, model: &mut Model) -> Result<String> {
        let (prime, curve) = self.parameters(model)?;
        let cells = self.order.saturating_add(1);
        if cells.saturating_mul(3).saturating_sub(1) > MAX_SHOWN_BYTES {
            return Err(output_limit());
        }
        let coefficients = self.calculate(model, &curve, prime)?;
        let report = format!(
            "rvKummerHilbChar[{}, {}, {}, {}, {}, {}] = {coefficients:?}",
            self.curve.text, self.length, self.extension, self.degree, self.prime, self.order,
        );
        if report.len() > MAX_SHOWN_BYTES {
            return Err(output_limit());
        }
        Ok(report)
    }

    pub(crate) fn coefficient(&self, model: &mut Model) -> Result<(Prime, u32)> {
        let (prime, curve) = self.parameters(model)?;
        if let Some(value) = self.direct(model, &curve, prime, self.order)? {
            return Ok((prime, value));
        }
        Ok((
            prime,
            self.calculate(model, &curve, prime)?[self.order] as u32,
        ))
    }

    fn parameters(&self, model: &mut Model) -> Result<(Prime, Rc<Curve>)> {
        let curve = model
            .kummer_curve(&self.curve.text)
            .map_err(|error| error.at(self.curve.span))?;
        let prime = if self.prime == 2 {
            Prime::validated(2)
        } else if self.prime < 3 || self.prime.is_multiple_of(2) {
            return Err(Error::plain(
                "Kummer Hilbert-scheme characteristic modulus must be prime",
            ));
        } else {
            Prime::new(model, self.prime)?
        };
        Ok((prime, curve))
    }

    fn direct(
        &self,
        model: &mut Model,
        curve: &Curve,
        prime: Prime,
        order: usize,
    ) -> Result<Option<u32>> {
        if order == 0 {
            return Ok(Some(1));
        }
        match graded_index(self.length, self.degree) {
            None => Ok(Some(0)),
            Some((_, 0, _)) if order > 1 => Ok(Some(0)),
            Some((_, 0, 0)) => Ok(Some(prime.negate(1))),
            Some((_, 0, twist)) => {
                let ell = u64::from(prime.get());
                let q = power(
                    model,
                    curve.ring.prime() % ell,
                    curve.ring.degree() as u128,
                    ell,
                )?;
                let relative_q = power(model, q, self.extension, ell)?;
                let eigenvalue = power(model, relative_q, twist, ell)? as u32;
                Ok(Some(prime.negate(eigenvalue)))
            }
            Some(_) => Ok(None),
        }
    }

    fn calculate(&self, model: &mut Model, curve: &Curve, prime: Prime) -> Result<Vec<u64>> {
        let entries = self.order.saturating_add(1);
        if entries.saturating_mul(2) > MAX_PIVOT_CELLS {
            return Err(storage_limit());
        }
        model.spend(entries)?;
        let (length, degree, twist) = graded_index(self.length, self.degree).unwrap_or((0, 0, 0));
        if self.order == 0 || degree == 0 {
            let mut coefficients = vec![0; entries];
            coefficients[0] = 1;
            if self.order > 0 {
                coefficients[1] = u64::from(
                    self.direct(model, curve, prime, 1)?
                        .expect("constant or Tate factor"),
                );
            }
            return Ok(coefficients);
        }
        let length = usize::try_from(length).map_err(|_| storage_limit())?;
        let degree = usize::try_from(degree).map_err(|_| storage_limit())?;
        let cells = length
            .saturating_add(1)
            .saturating_mul(degree.saturating_add(1));
        // The coefficient table, Witt corrections, and symmetric powers are temporary.
        let total = cells
            .saturating_add(entries * 2)
            .saturating_add(length.saturating_add(1).saturating_mul(4))
            .saturating_add(32);
        if total > MAX_PIVOT_CELLS {
            return Err(storage_limit());
        }
        model.spend(total)?;
        let precision = self.order.ilog(prime.get() as usize) as usize + 1;
        let ring = Residues::new(prime, precision)?;
        let mut surface = SurfacePowers::new(model, curve, self.extension, twist, ring)?;
        recover(model, ring, self.order, |model, j| {
            surface.trace(model, j, length, degree)
        })
    }
}

fn output_limit() -> Error {
    Error::plain("Kummer Hilbert-scheme characteristic report exceeds the output limit")
}

fn power(model: &mut Model, mut base: u64, mut exponent: u128, modulus: u64) -> Result<u64> {
    let mut result = 1;
    while exponent != 0 {
        model.spend(1)?;
        if exponent & 1 != 0 {
            result = mul_mod(result, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exponent >>= 1;
    }
    Ok(result)
}

// Big Witt coordinates of det(1-TF) solve the integral ghost equations.
fn recover(
    model: &mut Model,
    ring: Residues,
    order: usize,
    mut trace: impl FnMut(&mut Model, usize) -> Result<u64>,
) -> Result<Vec<u64>> {
    let modulus = ring.modulus;
    let ell = u64::from(ring.prime.get());
    let mut corrections = vec![0; order + 1];
    let mut polynomial = vec![0; order + 1];
    polynomial[0] = 1;
    for j in 1..=order {
        let residue = sub_mod(trace(model, j)?, corrections[j], modulus);
        let mut unit = j as u64;
        let mut ell_part = 1;
        let mut valuation = 0;
        while unit.is_multiple_of(ell) {
            unit /= ell;
            ell_part *= ell;
            valuation += 1;
        }
        if !residue.is_multiple_of(ell_part) {
            return Err(Error::plain(
                "Kummer Hilbert-scheme Witt ghost division is inexact",
            ));
        }
        let reduced = ring.reduce_to(ring.precision - valuation)?;
        let inverse = reduced.inverse_unit(model, &(unit % reduced.modulus))?;
        let coordinate = mul_mod(residue / ell_part, inverse, reduced.modulus);
        if coordinate == 0 {
            continue;
        }
        // An uncertainty of ell^(precision-valuation) in the coordinate vanishes
        // modulo the original precision after multiplying its contribution by j.
        let mut power = coordinate;
        for target in (2 * j..=order).step_by(j) {
            model.spend(1)?;
            power = mul_mod(power, coordinate, modulus);
            let term = mul_mod(j as u64, power, modulus);
            corrections[target] = add_mod(corrections[target], term, modulus);
        }
        let factor = coordinate % ell;
        if factor != 0 {
            multiply_factor(model, &mut polynomial, j, factor, ell)?;
        }
    }
    Ok(polynomial)
}

// Multiply by 1-aT^shift, truncating at the size of the supplied polynomial.
fn multiply_factor(
    model: &mut Model,
    polynomial: &mut [u64],
    shift: usize,
    a: u64,
    modulus: u64,
) -> Result<()> {
    model.spend(polynomial.len() - shift)?;
    for j in (shift..polynomial.len()).rev() {
        polynomial[j] = sub_mod(
            polynomial[j],
            mul_mod(a, polynomial[j - shift], modulus),
            modulus,
        );
    }
    Ok(())
}

struct SurfacePowers {
    ring: Residues,
    extension: u128,
    orbits: Vec<usize>,
    frobenius: Matrix<u64>,
    matrix: Matrix<u64>,
    field_order: u64,
    q: u64,
    twist: u64,
    twist_power: u64,
}

impl SurfacePowers {
    fn new(
        model: &mut Model,
        curve: &Curve,
        extension: u128,
        twist: u128,
        ring: Residues,
    ) -> Result<Self> {
        let q = curve.ring.size_of_residue_field()?;
        let orbits = node_orbits(model, curve, q)?;
        let weil = curve.weil(model)?;
        let last = |i: usize| sub_mod(0, ring.coefficient(weil[i]), ring.modulus);
        let companion = Matrix::from_entries(
            4,
            vec![
                0,
                0,
                0,
                last(4),
                1,
                0,
                0,
                last(3),
                0,
                1,
                0,
                last(2),
                0,
                0,
                1,
                last(1),
            ],
        )?;
        let mut frobenius = hw_algebra::identity(model, &ring, 4)?;
        let mut factor = companion;
        let mut exponent = extension;
        while exponent != 0 {
            if exponent & 1 != 0 {
                frobenius = hw_algebra::multiply(model, &ring, &frobenius, &factor)?;
            }
            exponent >>= 1;
            if exponent != 0 {
                factor = hw_algebra::multiply(model, &ring, &factor, &factor)?;
            }
        }
        let matrix = hw_algebra::identity(model, &ring, 4)?;
        let field_order = power(model, q % ring.modulus, extension, ring.modulus)?;
        let twist = power(model, field_order, twist, ring.modulus)?;
        Ok(Self {
            ring,
            extension,
            orbits,
            frobenius,
            matrix,
            field_order,
            q: 1,
            twist,
            twist_power: 1,
        })
    }

    fn trace(&mut self, model: &mut Model, j: usize, length: usize, degree: usize) -> Result<u64> {
        let modulus = self.ring.modulus;
        self.matrix = hw_algebra::multiply(model, &self.ring, &self.matrix, &self.frobenius)?;
        self.q = mul_mod(self.q, self.field_order, modulus);
        self.twist_power = mul_mod(self.twist_power, self.twist, modulus);
        let h2 = self.middle_polynomial(model, j)?;
        let trace = fock_trace(model, length, degree, self.q, modulus, &h2)?;
        Ok(mul_mod(trace, self.twist_power, modulus))
    }

    fn middle_polynomial(&self, model: &mut Model, j: usize) -> Result<Vec<u64>> {
        let q = self.q;
        let modulus = self.ring.modulus;
        let a = |a, b| add_mod(a, b, modulus);
        let s = |a, b| sub_mod(a, b, modulus);
        let p = |a, b| mul_mod(a, b, modulus);
        let entry = |r, c| *self.matrix.entry(r, c);
        let (mut trace, mut pairs) = (0, 0);
        model.spend(16)?;
        for r in 0..4 {
            trace = a(trace, entry(r, r));
            for c in r + 1..4 {
                pairs = a(
                    pairs,
                    s(p(entry(r, r), entry(c, c)), p(entry(r, c), entry(c, r))),
                );
            }
        }
        // The exterior square has two Q-eigenvalues and the cross quartic.
        let q2 = p(q, q);
        let linear = s(a(q, q), pairs);
        let quadratic = p(q, s(a(p(trace, trace), linear), pairs));
        let mut polynomial = vec![1, linear, quadratic, p(q2, linear), p(q2, q2)];
        polynomial.resize(7, 0);
        multiply_factor(model, &mut polynomial[..6], 1, q, modulus)?;
        multiply_factor(model, &mut polynomial, 1, q, modulus)?;
        for &orbit in &self.orbits {
            let d = orbit as u128;
            let degree_mod_d = (self.extension % d) * (j as u128 % d) % d;
            let count = gcd(d, degree_mod_d) as usize;
            let length = orbit / count;
            let factor = power(model, q, length as u128, modulus)?;
            for _ in 0..count {
                polynomial.resize(polynomial.len() + length, 0);
                multiply_factor(model, &mut polynomial, length, factor, modulus)?;
            }
        }
        Ok(polynomial)
    }
}

fn reciprocal(
    model: &mut Model,
    polynomial: &[u64],
    limit: usize,
    modulus: u64,
) -> Result<Vec<u64>> {
    let mut values = Vec::with_capacity(limit + 1);
    values.push(1);
    for i in 1..=limit {
        let terms = i.min(polynomial.len() - 1);
        model.spend(terms)?;
        let mut value = 0;
        for j in 1..=terms {
            value = sub_mod(
                value,
                mul_mod(polynomial[j], values[i - j], modulus),
                modulus,
            );
        }
        values.push(value);
    }
    Ok(values)
}

fn fock_trace(
    model: &mut Model,
    length: usize,
    degree: usize,
    q: u64,
    modulus: u64,
    h2: &[u64],
) -> Result<u64> {
    let width = degree + 1;
    let symmetry = [
        vec![1; 1 + length.min(degree)],
        reciprocal(model, h2, length.min(degree), modulus)?,
        reciprocal(
            model,
            &[1, sub_mod(0, mul_mod(q, q, modulus), modulus)],
            length.min(degree / 2),
            modulus,
        )?,
    ];
    model.spend((length + 1) * width)?;
    let mut coefficients = vec![0; (length + 1) * width];
    coefficients[0] = 1;
    let mut tate = 1;
    for m in 1..=length.min(width) {
        for (b, complete) in symmetry.iter().enumerate() {
            let weight = m - 1 + b;
            // The (m,b)=(1,0) factor is 1/(1-x). Summing over lengths at
            // the end accounts for it without recurring on the degree-zero factor.
            if weight == 0 {
                continue;
            }
            let limit = (length / m).min(degree / weight);
            if limit == 0 {
                continue;
            }
            model.spend(limit)?;
            let mut factor = 1;
            let scaled: Vec<_> = complete[1..=limit]
                .iter()
                .map(|&h| {
                    factor = mul_mod(factor, tate, modulus);
                    mul_mod(factor, h, modulus)
                })
                .collect();
            for n in (0..=length - m).rev() {
                let max_degree = (2 * n).min(degree - weight);
                model.spend(max_degree + 1)?;
                for d in 0..=max_degree {
                    let old = coefficients[n * width + d];
                    if old == 0 {
                        continue;
                    }
                    let terms = ((length - n) / m).min((degree - d) / weight);
                    model.spend(terms)?;
                    for (i, &h) in scaled[..terms].iter().enumerate() {
                        let j = i + 1;
                        let index = (n + j * m) * width + d + j * weight;
                        coefficients[index] =
                            add_mod(coefficients[index], mul_mod(old, h, modulus), modulus);
                    }
                }
            }
        }
        tate = mul_mod(tate, q, modulus);
    }
    Ok(coefficients
        .chunks_exact(width)
        .fold(0, |sum, row| add_mod(sum, row[degree], modulus)))
}
