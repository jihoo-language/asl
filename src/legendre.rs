use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;

pub(crate) enum Query {
    Trace {
        prime: usize,
        lambda: usize,
    },
    Hensel {
        prime: usize,
        lambda: usize,
        precision: usize,
    },
    Dwork {
        prime: usize,
        lambda: usize,
        precision: usize,
    },
    Period {
        prime: usize,
        z: i128,
        order: usize,
        precision: usize,
    },
    Congruence {
        prime: usize,
        precision: usize,
        z: Option<i128>,
    },
}

impl Query {
    pub(crate) fn run(self, model: &mut Model) -> Result<String> {
        let line = match self {
            Self::Trace { prime, lambda } => {
                let curve = Curve::new(model, prime, lambda)?;
                let trace = curve.trace(model)?;
                format!(
                    "legendreTrace[{prime}, {lambda}] = {trace}, points={}",
                    curve.point_count(trace)
                )
            }
            Self::Hensel {
                prime,
                lambda,
                precision,
            } => {
                let (modulus, trace, unit) =
                    Curve::new(model, prime, lambda)?.hensel(model, precision)?;
                let root =
                    unit.map_or_else(|| "supersingular".into(), |root| format!("unit={root}"));
                format!("legendreHensel[{prime}, {lambda}, {precision}] mod {modulus} = trace={trace}, {root}")
            }
            Self::Dwork {
                prime,
                lambda,
                precision,
            } => {
                let result = Curve::new(model, prime, lambda)?.dwork(model, precision)?;
                let Dwork {
                    modulus,
                    teich,
                    previous,
                    current,
                    raw,
                    unit,
                    trace,
                } = result;
                let root = raw.zip(unit).map_or_else(
                    || "supersingular".into(),
                    |(raw, root)| format!("ratio={raw}, unit={root}"),
                );
                format!(
                    "legendreDwork[{prime}, {lambda}, {precision}] mod {modulus} = z={teich}, B_{precision}={current}, B_{}={previous}, {root}, trace={trace}",
                    precision - 1
                )
            }
            Self::Period {
                prime,
                z,
                order,
                precision,
            } => {
                let (modulus, value) = period(model, prime, z, order, precision)?;
                format!("dworkPeriod[{prime}, {z}, {order}, {precision}] mod {modulus} = {value}")
            }
            Self::Congruence {
                prime,
                precision,
                z: Some(z),
            } => {
                let (modulus, left, right) = evaluated_congruence(model, prime, precision, z)?;
                format!("dworkCongruence[{prime}, {precision}, {z}] mod {modulus} = left={left}, right={right}")
            }
            Self::Congruence {
                prime,
                precision,
                z: None,
            } => {
                let (modulus, count, equal) = polynomial_congruence(model, prime, precision)?;
                let result = if equal { "yes" } else { "no" };
                format!("dworkCongruence[{prime}, {precision}] mod {modulus} = {result}, coefficients={count}")
            }
        };
        Ok(line)
    }
}

fn odd_prime(model: &mut Model, prime: usize) -> Result<u64> {
    if prime < 3 || prime.is_multiple_of(2) {
        return Err(Error::plain("Dwork computations require an odd prime"));
    }
    Ok(u64::from(Prime::new(model, prime)?.get()))
}

struct Curve {
    prime: u64,
    lambda: u64,
}

struct Dwork {
    modulus: u64,
    teich: u64,
    previous: u64,
    current: u64,
    raw: Option<u64>,
    unit: Option<u64>,
    trace: i64,
}

impl Curve {
    fn new(model: &mut Model, prime: usize, lambda: usize) -> Result<Self> {
        let prime = odd_prime(model, prime)?;
        let lambda = lambda as u64;
        if lambda >= prime {
            return Err(Error::plain(
                "Legendre parameter must be represented by an integer in 0..p",
            ));
        }
        if lambda < 2 {
            return Err(Error::plain("Legendre parameters 0 and 1 are singular"));
        }
        Ok(Self { prime, lambda })
    }

    fn trace(&self, model: &mut Model) -> Result<i64> {
        let p = self.prime;
        model.spend(p as usize)?;
        let mut sum = 0;
        for x in 0..p {
            let rhs = multiply(
                multiply(x, subtract(x, 1, p), p),
                subtract(x, self.lambda, p),
                p,
            );
            if rhs != 0 {
                let character = modular_power(model, rhs, (p - 1) / 2, p)?;
                sum += if character == 1 { 1 } else { -1 };
            }
        }
        Ok(-sum)
    }

    fn point_count(&self, trace: i64) -> u64 {
        (self.prime as i64 + 1 - trace) as u64
    }

    fn hensel(&self, model: &mut Model, precision: usize) -> Result<(u64, i64, Option<u64>)> {
        let modulus = Modulus::new(self.prime, precision)?;
        let trace = self.trace(model)?;
        let unit = hensel(model, &modulus, trace)?;
        Ok((modulus.value(), trace, unit))
    }

    fn dwork(&self, model: &mut Model, precision: usize) -> Result<Dwork> {
        let modulus = Modulus::new(self.prime, precision)?;
        let p = self.prime;
        let value = modulus.value();
        let terms = value;
        let prefix = modulus.powers[precision - 1];
        let divided = multiply(self.lambda, inverse(model, 16 % p, p)?, p);
        let teich = modular_power(model, divided, prefix, value)?;
        let (current, previous) = periods(model, &modulus, teich, terms, prefix)?;
        let trace = self.trace(model)?;
        let mut raw = None;
        let mut unit = None;
        if !current.is_multiple_of(p) && !previous.is_multiple_of(p) {
            let ratio = multiply(current, inverse(model, previous, value)?, value);
            raw = Some(ratio);
            unit = Some(if p % 4 == 3 {
                subtract(0, ratio, value)
            } else {
                ratio
            });
        }
        if hensel(model, &modulus, trace)? != unit {
            return Err(Error::plain(
                "Dwork period ratio and point-count Hensel root disagree",
            ));
        }
        Ok(Dwork {
            modulus: value,
            teich,
            previous,
            current,
            raw,
            unit,
            trace,
        })
    }
}

struct Modulus {
    prime: u64,
    powers: Vec<u64>,
}

impl Modulus {
    fn new(prime: u64, precision: usize) -> Result<Self> {
        if precision == 0 {
            return Err(Error::plain("p-adic precision must be positive"));
        }
        let mut powers = vec![1u64];
        for _ in 0..precision {
            powers.push(
                powers
                    .last()
                    .unwrap()
                    .checked_mul(prime)
                    .ok_or_else(|| Error::plain("p-adic modulus is out of range for u64"))?,
            );
        }
        Ok(Self { prime, powers })
    }

    fn value(&self) -> u64 {
        *self.powers.last().unwrap()
    }
}

fn power_of_prime(p: u64, exponent: usize) -> Result<u64> {
    let mut power = 1u64;
    for _ in 0..exponent {
        power = power
            .checked_mul(p)
            .ok_or_else(|| Error::plain("Dwork truncation term count is out of range"))?;
    }
    Ok(power)
}

fn admit_terms(model: &mut Model, terms: u64) -> Result<()> {
    let terms = usize::try_from(terms)
        .map_err(|_| Error::plain("Dwork truncation term count is out of range"))?;
    model.spend(terms)
}

fn add(a: u64, b: u64, modulus: u64) -> u64 {
    ((u128::from(a) + u128::from(b)) % u128::from(modulus)) as u64
}

fn subtract(a: u64, b: u64, modulus: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        modulus - (b - a)
    }
}

fn multiply(a: u64, b: u64, modulus: u64) -> u64 {
    ((u128::from(a) * u128::from(b)) % u128::from(modulus)) as u64
}

fn modular_power(model: &mut Model, mut a: u64, mut exponent: u64, modulus: u64) -> Result<u64> {
    model.spend(2 * (u64::BITS - exponent.leading_zeros()) as usize)?;
    let mut power = 1;
    while exponent > 0 {
        if exponent & 1 == 1 {
            power = multiply(power, a, modulus);
        }
        a = multiply(a, a, modulus);
        exponent >>= 1;
    }
    Ok(power)
}

fn inverse(model: &mut Model, value: u64, modulus: u64) -> Result<u64> {
    let (mut old, mut next) = (modulus, value % modulus);
    let (mut a, mut b) = (0, 1);
    while next != 0 {
        model.spend(1)?;
        let quotient = old / next;
        (old, next) = (next, old % next);
        (a, b) = (b, subtract(a, multiply(quotient, b, modulus), modulus));
    }
    if old != 1 {
        return Err(Error::plain("Dwork arithmetic requires a p-adic unit"));
    }
    Ok(a)
}

fn hensel(model: &mut Model, modulus: &Modulus, trace: i64) -> Result<Option<u64>> {
    unit_hensel(model, modulus, trace, modulus.prime)
}

pub(crate) fn closed_point_hensel(
    model: &mut Model,
    prime: u64,
    precision: usize,
    trace: i64,
    q: u64,
) -> Result<Option<u64>> {
    unit_hensel(model, &Modulus::new(prime, precision)?, trace, q)
}

fn unit_hensel(model: &mut Model, modulus: &Modulus, trace: i64, norm: u64) -> Result<Option<u64>> {
    let p = modulus.prime;
    let a = trace.rem_euclid(p as i64) as u64;
    if a == 0 {
        return Ok(None);
    }
    let slope_inverse = inverse(model, a, p)?;
    let mut root = a;
    for powers in modulus.powers[1..].windows(2) {
        model.spend(1)?;
        let (old, next) = (powers[0], powers[1]);
        let trace = i128::from(trace).rem_euclid(i128::from(next)) as u64;
        let square = multiply(root, root, next);
        let linear = multiply(trace, root, next);
        let defect = add(subtract(square, linear, next), norm % next, next);
        let digit = multiply(subtract(0, defect / old, p), slope_inverse, p);
        root += digit * old;
    }
    Ok(Some(root))
}

struct Binomial {
    unit: u64,
    valuation: usize,
}

impl Binomial {
    fn new() -> Self {
        Self {
            unit: 1,
            valuation: 0,
        }
    }

    fn square(&self, modulus: &Modulus) -> u64 {
        let Some(&factor) = modulus.powers.get(2 * self.valuation) else {
            return 0;
        };
        let value = modulus.value();
        multiply(factor, multiply(self.unit, self.unit, value), value)
    }

    fn advance(&mut self, model: &mut Model, modulus: &Modulus, n_next: u64) -> Result<()> {
        // C_n=binom(2n,n). Divide out powers of p before inverting the
        // denominator, including when C_n² currently vanishes.
        let (numerator, numerator_v) = split_prime(4 * n_next - 2, modulus.prime);
        let (denominator, denominator_v) = split_prime(n_next, modulus.prime);
        model.spend(1 + numerator_v + denominator_v)?;
        self.valuation = (self.valuation + numerator_v)
            .checked_sub(denominator_v)
            .ok_or_else(|| Error::plain("central binomial valuation became negative"))?;
        let value = modulus.value();
        let denominator_inv = inverse(model, denominator, value)?;
        self.unit = multiply(
            multiply(self.unit, numerator, value),
            denominator_inv,
            value,
        );
        Ok(())
    }
}

fn split_prime(mut value: u64, prime: u64) -> (u64, usize) {
    let mut count = 0;
    while value.is_multiple_of(prime) {
        value /= prime;
        count += 1;
    }
    (value, count)
}

pub(crate) fn period_coefficients(
    model: &mut Model,
    prime: u64,
    precision: usize,
    apply: impl FnMut(&mut Model, u64, u64) -> Result<()>,
) -> Result<()> {
    let modulus = Modulus::new(prime, precision)?;
    for_each_coefficient(model, &modulus, modulus.value(), apply)
}

fn for_each_coefficient(
    model: &mut Model,
    modulus: &Modulus,
    terms: u64,
    mut apply: impl FnMut(&mut Model, u64, u64) -> Result<()>,
) -> Result<()> {
    admit_terms(model, terms)?;
    let mut binomial = Binomial::new();
    for next in 1..=terms {
        apply(model, next, binomial.square(modulus))?;
        if next != terms {
            binomial.advance(model, modulus, next)?;
        }
    }
    Ok(())
}

fn periods(
    model: &mut Model,
    modulus: &Modulus,
    z: u64,
    terms: u64,
    prefix: u64,
) -> Result<(u64, u64)> {
    let value = modulus.value();
    let mut total = 0;
    let mut partial = 0;
    let mut x = 1;
    for_each_coefficient(model, modulus, terms, |_, next, coefficient| {
        total = add(total, multiply(coefficient, x, value), value);
        if next == prefix {
            partial = total;
        }
        if next != terms {
            x = multiply(x, z, value);
        }
        Ok(())
    })?;
    Ok((total, partial))
}

fn period(
    model: &mut Model,
    prime: usize,
    z: i128,
    order: usize,
    precision: usize,
) -> Result<(u64, u64)> {
    let p = odd_prime(model, prime)?;
    let modulus = Modulus::new(p, precision)?;
    let terms = power_of_prime(p, order)?;
    let (value, _) = periods(
        model,
        &modulus,
        z.rem_euclid(i128::from(modulus.value())) as u64,
        terms,
        terms,
    )?;
    Ok((modulus.value(), value))
}

fn evaluated_congruence(
    model: &mut Model,
    prime: usize,
    precision: usize,
    z: i128,
) -> Result<(u64, u64, u64)> {
    let p = odd_prime(model, prime)?;
    let modulus = Modulus::new(p, precision)?;
    let value = modulus.value();
    let longest = value
        .checked_mul(p)
        .ok_or_else(|| Error::plain("Dwork truncation term count is out of range"))?;
    let z = z.rem_euclid(i128::from(value)) as u64;
    let zp = modular_power(model, z, p, value)?;
    let (f_next, f) = periods(model, &modulus, z, longest, value)?;
    let (g, g_previous) = periods(model, &modulus, zp, value, value / p)?;
    Ok((
        value,
        multiply(f_next, g_previous, value),
        multiply(f, g, value),
    ))
}

fn polynomial_congruence(
    model: &mut Model,
    prime: usize,
    precision: usize,
) -> Result<(u64, usize, bool)> {
    let p = odd_prime(model, prime)?;
    let modulus = Modulus::new(p, precision)?;
    let value = modulus.value();
    let stride = p as usize;
    let terms = usize::try_from(value)
        .ok()
        .and_then(|n| n.checked_mul(stride))
        .ok_or_else(polynomial_limit)?;
    let middle = terms / stride;
    let previous = middle / stride;
    let product = terms
        .checked_add(middle - stride)
        .ok_or_else(polynomial_limit)?;
    if terms
        .checked_add(product)
        .is_none_or(|cells| cells > MAX_PIVOT_CELLS)
    {
        return Err(polynomial_limit());
    }
    model.spend(terms + product)?;
    let mut coefficients = Vec::with_capacity(terms);
    let mut binomial = Binomial::new();
    for next in 1..=terms {
        coefficients.push(binomial.square(&modulus));
        if next != terms {
            binomial.advance(model, &modulus, next as u64)?;
        }
    }
    let mut difference = vec![0; product];
    for (outer, inner, subtract_product) in [
        (&coefficients[..], &coefficients[..previous], false),
        (&coefficients[..middle], &coefficients[..middle], true),
    ] {
        for (i, &a) in outer.iter().enumerate() {
            if a == 0 {
                continue;
            }
            model.spend(inner.len())?;
            for (j, &b) in inner.iter().enumerate() {
                let coefficient = multiply(a, b, value);
                let entry = &mut difference[i + stride * j];
                *entry = if subtract_product {
                    subtract(*entry, coefficient, value)
                } else {
                    add(*entry, coefficient, value)
                };
            }
        }
    }
    Ok((value, product, difference.iter().all(|&term| term == 0)))
}

fn polynomial_limit() -> Error {
    Error::plain("Dwork coefficientwise congruence exceeds the polynomial storage limit")
}
