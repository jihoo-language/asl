use super::{level::Level, Certificate};
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::syntax::Label;

pub(crate) enum Statement {
    Check(Label, Vec<Label>),
    Define(Label, Label, Vec<Label>),
    Scalar(Label, usize),
}

struct Depth {
    prime: u64,
    cap: u32,
    lower: u32,
    upper: u32,
    scalar: u64,
}

pub(crate) struct Bounds {
    q: u64,
    trace: i128,
    delta: i128,
    fundamental: i128,
    conductor: u64,
    depths: Vec<Depth>,
}

struct Observation {
    available: u32,
    depth: u32,
    power: u64,
    scalar: u64,
}

impl Observation {
    fn read(level: &Level, prime: u64) -> Self {
        let available = valuation(level.modulus, prime);
        let [[a, b], [c, d]] = level.frobenius.0;
        let mut depth = 0;
        let mut power = 1;
        for _ in 0..available {
            let next = power * prime;
            if [b, c, a.abs_diff(d)]
                .iter()
                .any(|value| !value.is_multiple_of(next))
            {
                break;
            }
            power = next;
            depth += 1;
        }
        Self {
            available,
            depth,
            power,
            scalar: a % power,
        }
    }
}

impl Bounds {
    pub(crate) fn read(model: &mut Model, base: &Label, levels: &[Label]) -> Result<Self> {
        let curve = model.ec_curve(&base.text)?;
        let (q, trace) = curve.trace(model)?;
        let trace = i128::from(trace);
        if trace.rem_euclid(i128::from(curve.ring.prime())) == 0 {
            return Err(Error::new(
                base.span,
                "Frobenius scalar-depth certification requires an ordinary elliptic curve",
            ));
        }
        let delta = trace * trace - 4 * i128::from(q);
        let (fundamental, conductor, mut depths) = quadratic_data(model, delta)?;
        for label in levels {
            model.spend(1)?;
            let level = model.ec_torsion_level(&label.text)?;
            if !level.is_plain() || level.base_name.as_ref() != base.text {
                return Err(Error::new(label.span, "endomorphism-order evidence requires a plain level on the same named base curve"));
            }
            for bound in &mut depths {
                let observed = Observation::read(&level, bound.prime);
                if observed.depth > bound.cap {
                    return Err(inconsistent(label));
                }
                let common = bound.prime.pow(observed.depth.min(bound.lower));
                if bound.scalar % common != observed.scalar % common {
                    return Err(inconsistent(label));
                }
                if observed.depth > bound.lower {
                    bound.lower = observed.depth;
                    bound.scalar = observed.scalar;
                }
                // A nonscalar action at the next available power pins
                // down the index. Reaching the conductor cap does too.
                if observed.depth < observed.available {
                    bound.upper = bound.upper.min(observed.depth);
                }
                if bound.lower > bound.upper {
                    return Err(inconsistent(label));
                }
            }
        }
        Ok(Self {
            q,
            trace,
            delta,
            fundamental,
            conductor,
            depths,
        })
    }

    fn indices(&self) -> (u64, u64) {
        self.depths.iter().fold((1, 1), |(lower, upper), depth| {
            (
                lower * depth.prime.pow(depth.lower),
                upper * depth.prime.pow(depth.upper),
            )
        })
    }

    fn discriminant(&self) -> Option<i128> {
        let (lower, upper) = self.indices();
        (lower == upper).then(|| self.fundamental * i128::from(self.conductor / lower).pow(2))
    }

    pub(crate) fn certify(&self, model: &mut Model, curve: Label) -> Result<Certificate> {
        let discriminant = self.discriminant().ok_or_else(|| Error::new(
            curve.span, "insufficient Frobenius scalar-depth evidence; use rvECEndTorsionCheck for the missing prime powers",
        ))?;
        let model_curve = model.ec_curve(&curve.text)?;
        let j = model_curve.j_invariant(model)?;
        Certificate::from_discriminant(model, curve.text, self.q, self.trace, discriminant, j)
    }

    pub(crate) fn report(&self, base: &str) -> String {
        let (lower, upper) = self.indices();
        let depths: Vec<_> = self
            .depths
            .iter()
            .map(|depth| (depth.prime, depth.lower, depth.upper))
            .collect();
        let missing: Vec<_> = self
            .depths
            .iter()
            .filter(|depth| depth.lower < depth.upper)
            .map(|depth| depth.prime.pow(depth.lower + 1))
            .collect();
        let discriminant = self
            .discriminant()
            .map_or_else(|| "none".into(), |d| d.to_string());
        format!(
            "rvECEndTorsionCheck[{base}] = q={}, trace={}, frobenius_discriminant={}, fundamental_discriminant={}, frobenius_conductor={}, index_bounds=[{lower}, {upper}], conductor_bounds=[{}, {}], discriminant={discriminant}, depths={depths:?}, missing_levels={missing:?}",
            self.q, self.trace, self.delta, self.fundamental, self.conductor, self.conductor / upper, self.conductor / lower,
        )
    }
}

pub(crate) fn scalar_report(
    model: &mut Model,
    level: &Level,
    label: &Label,
    p: usize,
) -> Result<String> {
    if !level.is_plain() {
        return Err(Error::new(
            label.span,
            "scalar-depth query requires a plain elliptic torsion level",
        ));
    }
    let prime = if p == 2 {
        2
    } else {
        u64::from(Prime::new(model, p)?.get())
    };
    let obs = Observation::read(level, prime);
    if obs.available == 0 {
        return Err(Error::new(
            label.span,
            "scalar-depth prime must divide the named level",
        ));
    }
    let Observation {
        available,
        depth,
        power,
        scalar,
    } = obs;
    let mut report = format!(
        "rvECPlainLevelScalar[{},{}] = available_depth={available}, scalar_depth={depth}, scalar={scalar}, scalar_modulus={power}",
        label.text, prime,
    );
    let remaining = level.modulus / power;
    if depth != 0 && remaining > 1 {
        let [[a, b], [c, d]] = level.frobenius.0;
        let signed = |x| i128::from(x) - i128::from(scalar);
        let diagonal = |x| (signed(x) / i128::from(power)).rem_euclid(i128::from(remaining)) as u64;
        let columns = [[diagonal(a), b / power], [c / power, diagonal(d)]];
        let base = model.ec_curve(&level.base_name)?;
        let (q, trace) = base.trace(model)?;
        let scale = i128::from(power);
        let scalar = i128::from(scalar);
        let trace = i128::from(trace);
        let trace_numerator = trace - 2 * scalar;
        let norm_numerator = i128::from(q) - trace * scalar + scalar * scalar;
        let matrix = super::level::Matrix(columns);
        if trace_numerator % scale != 0
            || norm_numerator % (scale * scale) != 0
            || i128::from(matrix.trace(remaining))
                != (trace_numerator / scale).rem_euclid(i128::from(remaining))
            || i128::from(matrix.det(remaining))
                != (norm_numerator / (scale * scale)).rem_euclid(i128::from(remaining))
        {
            return Err(Error::plain(
                "divided Frobenius failed its integral trace or determinant check",
            ));
        }
        report.push_str(&format!(
            ", divided_level={remaining}, divided_columns={columns:?}, divided_trace={}, divided_norm={}",
            trace_numerator / scale, norm_numerator / (scale * scale),
        ));
    }
    Ok(report)
}

fn valuation(mut n: u64, prime: u64) -> u32 {
    let mut depth = 0;
    while n.is_multiple_of(prime) {
        n /= prime;
        depth += 1;
    }
    depth
}

pub(crate) fn frobenius_conductor(model: &mut Model, delta: i128) -> Result<u64> {
    let (_, conductor, _) = quadratic_data(model, delta)?;
    Ok(conductor)
}

fn quadratic_data(model: &mut Model, delta: i128) -> Result<(i128, u64, Vec<Depth>)> {
    if delta >= 0 || !matches!(delta.rem_euclid(4), 0 | 1) {
        return Err(Error::plain(
            "ordinary Frobenius does not have a negative quadratic discriminant",
        ));
    }
    let factors = factor(model, delta.unsigned_abs())?;
    let squarefree = -factors
        .iter()
        .filter(|(_, e)| e % 2 != 0)
        .map(|&(p, _)| p as i128)
        .product::<i128>();
    let even_discriminant = squarefree.rem_euclid(4) != 1;
    let fundamental = if even_discriminant {
        4 * squarefree
    } else {
        squarefree
    };
    let mut conductor = 1u64;
    let mut depths = Vec::new();
    for (prime, multiplicity) in factors {
        let adjust = u32::from(prime == 2 && even_discriminant);
        let cap = (multiplicity / 2).checked_sub(adjust).ok_or_else(|| {
            Error::plain("Frobenius discriminant has inconsistent two-adic valuation")
        })?;
        if cap > 0 {
            let prime = u64::try_from(prime).map_err(|_| quadratic_range())?;
            conductor = conductor
                .checked_mul(prime.pow(cap))
                .ok_or_else(quadratic_range)?;
            depths.push(Depth {
                prime,
                cap,
                lower: 0,
                upper: cap,
                scalar: 0,
            });
        }
    }
    let reconstructed = fundamental
        .checked_mul(i128::from(conductor).pow(2))
        .ok_or_else(quadratic_range)?;
    if reconstructed != delta {
        return Err(Error::plain(
            "Frobenius conductor failed discriminant reconstruction",
        ));
    }
    Ok((fundamental, conductor, depths))
}

fn factor(model: &mut Model, mut n: u128) -> Result<Vec<(u128, u32)>> {
    let mut factors = Vec::new();
    let mut candidate = 2;
    while candidate <= n / candidate {
        model.spend(1)?;
        let mut exponent = 0;
        while n.is_multiple_of(candidate) {
            n /= candidate;
            exponent += 1;
        }
        if exponent != 0 {
            factors.push((candidate, exponent));
        }
        candidate += if candidate == 2 { 1 } else { 2 };
    }
    if n > 1 {
        factors.push((n, 1));
    }
    Ok(factors)
}

fn inconsistent(label: &Label) -> Error {
    Error::new(label.span, "inconsistent Frobenius scalar-depth evidence")
}

fn quadratic_range() -> Error {
    Error::plain("Frobenius conductor exceeds the supported integer range")
}
