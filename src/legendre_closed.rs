use crate::error::{Error, Result};
use crate::legendre;
use crate::model::Model;
use crate::rv_unram::{require_scalar, Ring};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Parameter,
    Factor(usize),
    Unit,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Parameter => "rvLegendreParameter",
            Self::Factor(_) => "rvLegendreFactor",
            Self::Unit => "rvLegendreUnit",
        }
    }

    pub(crate) fn apply(self, model: &mut Model, ring: &Ring, lambda: &[u64]) -> Result<Vec<u64>> {
        validate(ring, lambda)?;
        let z = parameter(model, ring, lambda)?;
        if matches!(self, Self::Parameter) {
            return Ok(z);
        }
        let periods = Periods::compute(model, ring, &z)?;
        let first = periods.signed_ratio(model, ring)?.ok_or_else(|| {
            Error::plain("Legendre closed point is supersingular and has no unit period factor")
        })?;
        match self {
            Self::Factor(index) => ring.frobenius(model, &first, index),
            Self::Unit => ring.conjugate_sum_or_product(model, &first, true),
            Self::Parameter => unreachable!(),
        }
    }
}

pub(crate) enum Query {
    Data,
    Trace,
    Hensel,
}

impl Query {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            Self::Data => "rvLegendreData",
            Self::Trace => "rvLegendreTrace",
            Self::Hensel => "rvLegendreHensel",
        }
    }

    pub(crate) fn report(self, model: &mut Model, ring: &Ring, lambda: &[u64]) -> Result<String> {
        validate(ring, lambda)?;
        match self {
            Self::Data => {
                let z = parameter(model, ring, lambda)?;
                let periods = Periods::compute(model, ring, &z)?;
                let current = ring.valuation(model, &periods.current)?;
                let previous = ring.valuation(model, &periods.previous)?;
                let status = if current == 0 && previous == 0 {
                    "ordinary"
                } else {
                    "supersingular"
                };
                let describe = |v: usize| {
                    if v == ring.precision() {
                        format!(">={v}")
                    } else {
                        v.to_string()
                    }
                };
                Ok(format!(
                    "rvLegendreData mod {} = z={z:?}, B_{}={:?}, B_{}={:?}, valuations=({}, {}), {status}",
                    ring.modulus(), ring.precision(), periods.current, ring.precision()-1, periods.previous,
                    describe(current), describe(previous),
                ))
            }
            Self::Trace => {
                let q = ring.size_of_residue_field()?;
                let trace = ring.legendre_trace(model, lambda)?;
                let points = i128::from(q) + 1 - i128::from(trace);
                Ok(format!(
                    "rvLegendreTrace = q={q}, trace={trace}, points={points}"
                ))
            }
            Self::Hensel => {
                let q = ring.size_of_residue_field()?;
                let trace = ring.legendre_trace(model, lambda)?;
                let detail =
                    legendre::closed_point_hensel(model, ring.prime(), ring.precision(), trace, q)?
                        .map_or_else(|| "supersingular".to_owned(), |root| format!("unit={root}"));
                Ok(format!(
                    "rvLegendreHensel mod {} = q={q}, trace={trace}, {detail}",
                    ring.modulus()
                ))
            }
        }
    }
}

pub(crate) fn validate(ring: &Ring, lambda: &[u64]) -> Result<()> {
    if ring.prime() == 2 {
        return Err(Error::plain(
            "Legendre closed point requires an odd residue characteristic",
        ));
    }
    let residue = ring.residue(lambda);
    if residue[1..].iter().all(|&c| c == 0) && residue[0] <= 1 {
        return Err(Error::plain(
            "Legendre parameters 0 and 1 are singular in the residue field",
        ));
    }
    Ok(())
}

fn parameter(model: &mut Model, ring: &Ring, lambda: &[u64]) -> Result<Vec<u64>> {
    let inverse = ring.inverse(model, &ring.scalar(16))?;
    let divided = ring.scale(model, lambda, require_scalar(&inverse)?)?;
    ring.teichmuller(model, &divided)
}

struct Periods {
    current: Vec<u64>,
    previous: Vec<u64>,
}

impl Periods {
    fn compute(model: &mut Model, ring: &Ring, z: &[u64]) -> Result<Self> {
        let mut current = ring.scalar(0);
        let mut previous = ring.scalar(0);
        let mut power = ring.scalar(1);
        let terms = ring.modulus();
        let prefix = terms / ring.prime();
        legendre::period_coefficients(
            model,
            ring.prime(),
            ring.precision(),
            |model, next, coefficient| {
                if coefficient != 0 {
                    let term = ring.scale(model, &power, coefficient)?;
                    current = ring.add(model, &current, &term)?;
                }
                if next == prefix {
                    previous.clone_from(&current);
                }
                if next < terms {
                    power = ring.multiply(model, &power, z)?;
                }
                Ok(())
            },
        )?;
        Ok(Self { current, previous })
    }

    fn signed_ratio(&self, model: &mut Model, ring: &Ring) -> Result<Option<Vec<u64>>> {
        if ring.valuation(model, &self.current)? != 0 || ring.valuation(model, &self.previous)? != 0
        {
            return Ok(None);
        }
        // Integer polynomial coefficients are fixed by sigma. Thus this
        // denominator really is B_(s-1)(sigma z), and Frobenius permutes
        // all the signed ratios from this first one.
        let shifted = ring.frobenius(model, &self.previous, 1)?;
        let inverse = ring.inverse(model, &shifted)?;
        let ratio = ring.multiply(model, &self.current, &inverse)?;
        let signed = if ring.prime() % 4 == 3 {
            ring.negate(model, &ratio)?
        } else {
            ratio
        };
        Ok(Some(signed))
    }
}
