use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear;
use crate::rv_unram::{add_mod, mul_mod, sub_mod};
use crate::syntax::{Label, ParsedExpr};

use super::eigenpacket::{Element, Packet};
use super::hecke::{self, Operator};
use super::Integer;

pub(crate) enum Statement {
    Declare(Label, Label, usize, i128),
    Sturm(Label),
}

pub(crate) enum Expression {
    Reduce(Label, Box<ParsedExpr>),
    Fourier(Label, u128),
    Eisenstein(u128, usize, u128),
}

pub(crate) struct Specialization {
    pub(crate) source: Label,
    packet: Rc<Packet>,
    pub(crate) prime: Prime,
    root: u32,
    inverse_derivative: u32,
}

impl Specialization {
    pub(crate) fn new(
        model: &mut Model,
        source: Label,
        packet: Rc<Packet>,
        p: usize,
        root: i128,
    ) -> Result<Self> {
        let prime = checked_prime(model, p)?;
        let root = root.rem_euclid(i128::from(prime.get())) as u32;
        let (value, derivative) =
            evaluate_with_derivative(model, packet.polynomial(), prime, root)?;
        if value != 0 {
            return Err(Error::plain(
                "residual Hecke value is not a root of the generator polynomial",
            ));
        }
        if derivative == 0 {
            return Err(Error::plain("residual Hecke generator root must be simple"));
        }
        let inverse_derivative = prime_linear::invert(model, prime, derivative)?;
        Ok(Self {
            source,
            packet,
            prime,
            root,
            inverse_derivative,
        })
    }

    pub(crate) fn reduce(&self, model: &mut Model, value: &Element) -> Result<u32> {
        let (numerator, mut unit) = value.fraction_parts();
        let p = self.prime.get();
        let prime = u64::from(p);
        let mut modulus = prime;
        loop {
            model.spend(1)?;
            let (quotient, remainder) = unit.divide_u32(p);
            if remainder != 0 {
                break;
            }
            unit = quotient;
            modulus = modulus.checked_mul(prime).ok_or_else(precision)?;
        }
        let clearing = modulus / prime;
        let root = self.lift_root(model, modulus)?;
        let numerator = evaluate(model, numerator.iter().rev().copied(), root, modulus)?;
        if !numerator.is_multiple_of(clearing) {
            return Err(Error::plain(
                "Hecke eigenpacket value is nonintegral at the specified prime and root",
            ));
        }
        let inverse = prime_linear::invert(model, self.prime, unit.residue(p))?;
        Ok(self.prime.mul((numerator / clearing) as u32, inverse))
    }

    pub(crate) fn fourier(&self, model: &mut Model, n: u128) -> Result<u32> {
        if n == 0 {
            return Ok(0);
        }
        let operator = Operator::basis(n)?;
        let value = self.packet.image_element(model, &operator)?;
        // The normalization a_1=1 makes a_n the T_n eigenvalue. Its
        // characteristic polynomial is integral, so it is integral at the
        // chosen Hensel embedding, even if these coordinates have a pole.
        self.reduce(model, &value)
    }

    pub(crate) fn sturm(&self, model: &mut Model, name: &str) -> Result<String> {
        let weight = self.packet.weight;
        let bound = weight / 12;
        let p = self.prime.get();
        let label = format!("rvHeckeEisensteinSturm[{name}]");
        let Some(constant) = eisenstein_constant(model, weight, self.prime)? else {
            return Ok(format!(
                "{label} = unavailable modulo {p}: Eisenstein coefficient q^0 is nonintegral"
            ));
        };
        for n in 0..=bound {
            let eisenstein = if n == 0 {
                constant
            } else {
                divisor_sum(model, weight, n, self.prime)?
            };
            let eigenform = self.fourier(model, n)?;
            if eisenstein != eigenform {
                return Ok(format!(
                    "{label} = mismatch at q^{n} modulo {p} (eigenform {eigenform}, Eisenstein {eisenstein}; Sturm bound {bound})"
                ));
            }
        }
        Ok(format!(
            "{label} = congruent modulo {p} (Sturm bound {bound})"
        ))
    }

    fn lift_root(&self, model: &mut Model, target: u64) -> Result<u64> {
        let p = u64::from(self.prime.get());
        let mut root = u64::from(self.root);
        let mut step = p;
        while step < target {
            let next = step * p;
            let value = evaluate(model, self.packet.polynomial().iter().copied(), root, next)?;
            if !value.is_multiple_of(step) {
                return Err(Error::plain(
                    "inconsistent Hensel lift of the residual Hecke root",
                ));
            }
            let digit = self.prime.negate(
                self.prime
                    .mul((value / step) as u32, self.inverse_derivative),
            );
            root += u64::from(digit) * step;
            step = next;
        }
        Ok(root)
    }
}

pub(crate) fn eisenstein(
    model: &mut Model,
    weight: u128,
    p: usize,
    n: u128,
) -> Result<(Prime, u32)> {
    hecke::validate_weight(weight)?;
    if weight == 2 {
        return Err(Error::plain(
            "level-one Eisenstein coefficient requires an even weight at least 4",
        ));
    }
    let prime = checked_prime(model, p)?;
    let value = if n == 0 {
        eisenstein_constant(model, weight, prime)?.ok_or_else(|| {
            Error::plain("Eisenstein coefficient q^0 is nonintegral at the specified prime")
        })?
    } else {
        divisor_sum(model, weight, n, prime)?
    };
    Ok((prime, value))
}

fn checked_prime(model: &mut Model, p: usize) -> Result<Prime> {
    if p < 5 || p.is_multiple_of(2) {
        return Err(Error::plain(
            "residual Hecke specialization requires a prime at least 5",
        ));
    }
    Prime::new(model, p)
}

fn evaluate_with_derivative(
    model: &mut Model,
    coefficients: &[Integer],
    prime: Prime,
    root: u32,
) -> Result<(u32, u32)> {
    let (mut value, mut derivative) = (0, 0);
    for &c in coefficients {
        model.spend(1)?;
        derivative = prime.add(prime.mul(derivative, root), value);
        let signed = c.magnitude.residue(prime.get());
        value = prime.add(
            prime.mul(value, root),
            if c.negative {
                prime.negate(signed)
            } else {
                signed
            },
        );
    }
    Ok((value, derivative))
}

fn evaluate(
    model: &mut Model,
    coefficients: impl IntoIterator<Item = Integer>,
    root: u64,
    modulus: u64,
) -> Result<u64> {
    let mut value = 0;
    for c in coefficients {
        let (_, residue) = c.magnitude.divide(model, u128::from(modulus))?;
        let leading = mul_mod(value, root, modulus);
        value = if c.negative {
            sub_mod(leading, residue as u64, modulus)
        } else {
            add_mod(leading, residue as u64, modulus)
        };
    }
    Ok(value)
}

fn divisor_sum(model: &mut Model, weight: u128, n: u128, prime: Prime) -> Result<u32> {
    let p = u64::from(prime.get());
    let mut sum = 0;
    let limit = n.isqrt();
    let work = usize::try_from(limit)
        .map_err(|_| Error::plain("Eisenstein divisor sum exceeds the computation limit"))?;
    model.spend(work)?;
    for d in 1..=limit {
        if n.is_multiple_of(d) {
            let other = n / d;
            let first = power_mod(model, (d % u128::from(p)) as u64, weight - 1, p)?;
            sum = add_mod(sum, first, p);
            if d != other {
                let second = power_mod(model, (other % u128::from(p)) as u64, weight - 1, p)?;
                sum = add_mod(sum, second, p);
            }
        }
    }
    Ok(sum as u32)
}

fn eisenstein_constant(model: &mut Model, weight: u128, prime: Prime) -> Result<Option<u32>> {
    let weight = usize::try_from(weight).map_err(|_| bernoulli_limit())?;
    let rows = weight.checked_add(1).ok_or_else(bernoulli_limit)?;
    let work = weight.checked_mul(rows).ok_or_else(bernoulli_limit)? / 2;
    model.spend(work).map_err(|_| bernoulli_limit())?;

    // B_k = sum_m (-1)^m m! S(k,m)/(m+1). Clear just the p-part of
    // each possible denominator (m+1)*2k, then retain one extra digit.
    // This determines both integrality of -B_k/(2k) and its residue.
    let p = u64::from(prime.get());
    let two_k = 2 * weight as u64;
    let mut maximal = 1;
    while maximal <= rows as u64 / p {
        maximal *= p;
    }
    let clearing = maximal
        .checked_mul(p_part(two_k, p))
        .ok_or_else(precision)?;
    let modulus = clearing.checked_mul(p).ok_or_else(precision)?;

    let mut stirling = vec![0; rows];
    stirling[0] = 1;
    for k in 1..=weight {
        for m in (1..=k).rev() {
            stirling[m] = mul_mod(
                m as u64,
                add_mod(stirling[m], stirling[m - 1], modulus),
                modulus,
            );
        }
        stirling[0] = 0;
    }
    let phi = (modulus / p) * (p - 1);
    let mut sum = 0;
    for (m, value) in stirling
        .into_iter()
        .enumerate()
        .filter(|&(_, value)| value != 0)
    {
        let denominator = two_k * (m as u64 + 1);
        let part = p_part(denominator, p);
        let unit = denominator / part;
        let inverse = power_mod(model, unit % modulus, u128::from(phi - 1), modulus)?;
        let term = mul_mod(mul_mod(value, clearing / part, modulus), inverse, modulus);
        sum = if m.is_multiple_of(2) {
            sub_mod(sum, term, modulus)
        } else {
            add_mod(sum, term, modulus)
        };
    }
    Ok(sum
        .is_multiple_of(clearing)
        .then_some((sum / clearing) as u32))
}

fn p_part(mut n: u64, p: u64) -> u64 {
    let original = n;
    while n.is_multiple_of(p) {
        n /= p;
    }
    original / n
}

pub(super) fn power_mod(
    model: &mut Model,
    mut base: u64,
    mut exponent: u128,
    modulus: u64,
) -> Result<u64> {
    let mut value = 1;
    while exponent != 0 {
        model.spend(1)?;
        if exponent & 1 == 1 {
            value = mul_mod(value, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exponent /= 2;
    }
    Ok(value)
}

fn precision() -> Error {
    Error::plain("residual Hecke calculation needs prime-power precision beyond the unsigned 64-bit modulus limit")
}

fn bernoulli_limit() -> Error {
    Error::plain("Eisenstein Bernoulli coefficient exceeds the exact computation limit")
}
