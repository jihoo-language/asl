use crate::elliptic_modular::decimal;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::rv_ore::{Powers, Sum};
use crate::syntax::Label;

use super::{linear as lin, Operator, Space};
use reciprocal::Integer;

mod geometry;
mod reciprocal;

const MAX_EXPONENT: usize = 64;

pub(crate) enum Command {
    Zeta(Label, usize, usize),
    Geometry(Label, Label),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<String> {
        match self {
            Self::Zeta(name, ell, exponent) => {
                let space = model.manin_space(&name)?;
                let data = space.frob_data(model, ell, exponent)?;
                Ok(data.report(space.geometry.level, ell, exponent))
            }
            Self::Geometry(name, field) => {
                let space = model.manin_space(&name)?;
                let ring = model.rv_unram(&field)?;
                space.geometric(model, field.text.into(), ring)
            }
        }
    }
}

struct Data {
    q: Integer,
    trace: Integer,
    points: Integer,
    jacobian: Integer,
    numerator: Vec<Integer>,
}

impl Data {
    fn report(&self, level: usize, ell: usize, exponent: usize) -> String {
        let coefficients = self
            .numerator
            .iter()
            .copied()
            .map(decimal)
            .collect::<Vec<_>>()
            .join(", ");
        let q = decimal(self.q);
        format!(
            "X0({level}) at ℓ={ell}, extension degree {exponent}, q={q}, genus {}\n\
             zeta numerator (ascending) = [{coefficients}]\n\
             zeta denominator = (1-T)(1-{q}T)\n\
             Frobenius trace = {}; #X0 = {}; #Jac = {}",
            (self.numerator.len() - 1) / 2,
            decimal(self.trace),
            decimal(self.points),
            decimal(self.jacobian)
        )
    }
}

impl Space {
    fn frob_data(&self, model: &mut Model, ell: usize, exponent: usize) -> Result<Data> {
        require_good(model, self.geometry.prime, ell)?;
        if !(1..=MAX_EXPONENT).contains(&exponent) {
            return Err(Error::plain(format!(
                "Manin Frobenius extension degree must lie in 1..={MAX_EXPONENT}"
            )));
        }
        let hecke = self.operator(model, Operator::Hecke(ell))?;
        let cut = self.plus.restrict(model, &hecke)?;
        let characteristic = lin::characteristic(model, &cut)?;
        let genus = cut.len();
        let mut numerator = vec![0i128; 2 * genus + 1];
        let mut factor = vec![1i128];
        // T^g χ(T^-1 + ℓT) = sum_i c_i T^i (1 + ℓ T²)^(g-i).
        for (i, c) in characteristic.into_iter().enumerate().rev() {
            if c.denominator != 1 {
                return Err(certificate(
                    "the full cuspidal Hecke polynomial is integral",
                ));
            }
            for (entry, &term) in numerator[i..].iter_mut().zip(&factor) {
                *entry = checked_add(*entry, checked_mul(c.numerator, term)?)?;
            }
            if i != 0 {
                let mut next = vec![0i128; factor.len() + 2];
                for (j, &term) in factor.iter().enumerate() {
                    next[j] = checked_add(next[j], term)?;
                    next[j + 2] = checked_add(next[j + 2], checked_mul(ell as i128, term)?)?;
                }
                factor = next;
            }
        }
        let coefficients = numerator
            .into_iter()
            .map(Integer::from_i128)
            .collect::<Vec<_>>();
        reciprocal::reciprocity(model, &coefficients, Integer::from_i128(ell as i128))?;
        let (numerator, trace) = reciprocal::power(model, &coefficients, exponent)?;
        let q = reciprocal::q_power(model, ell, exponent)?;
        if exponent != 1 {
            reciprocal::reciprocity(model, &numerator, q)?;
        }
        let points = reciprocal::add(reciprocal::add(q, Integer::from_i128(1))?, trace.negate())?;
        let jacobian = numerator
            .iter()
            .copied()
            .try_fold(Integer::from_i128(0), reciprocal::add)?;
        if points.negative
            || points.signed().is_some_and(|value| value < 4)
            || jacobian.negative
            || jacobian.is_zero()
        {
            return Err(certificate(
                "the point count includes all four cusps and the Jacobian has positive order",
            ));
        }
        Ok(Data {
            q,
            trace,
            points,
            jacobian,
            numerator,
        })
    }

    pub(crate) fn frob_polynomial(
        &self,
        model: &mut Model,
        ell: usize,
        exponent: usize,
    ) -> Result<Value> {
        let data = self.frob_data(model, ell, exponent)?;
        let mut polynomial = Sum::default();
        for (degree, coefficient) in data.numerator.into_iter().enumerate() {
            let coefficient=coefficient.signed().ok_or_else(|| Error::plain(
                "Manin Frobenius polynomial exceeds RVOre's signed coefficient range; use rvManinZeta for wide exact coefficients"
            ))?;
            polynomial.add(
                model,
                Powers {
                    q: degree,
                    theta: 0,
                },
                coefficient,
            )?;
        }
        polynomial.finish(model)
    }
}

fn require_good(model: &mut Model, p: usize, ell: usize) -> Result<()> {
    if !(3..=31).contains(&ell) || ell == p {
        return Err(Error::plain(
            "Manin Frobenius requires a good prime ℓ ≤ 31, coprime to 2p",
        ));
    }
    Prime::new(model, ell)?;
    Ok(())
}

fn checked_add(a: i128, b: i128) -> Result<i128> {
    a.checked_add(b).ok_or_else(arithmetic)
}

fn checked_mul(a: i128, b: i128) -> Result<i128> {
    a.checked_mul(b).ok_or_else(arithmetic)
}

fn arithmetic() -> Error {
    Error::plain("Manin Frobenius arithmetic exceeds the exact coefficient range")
}

fn certificate(identity: &str) -> Error {
    Error::plain(format!(
        "Manin Eichler–Shimura certificate failed: {identity}"
    ))
}
