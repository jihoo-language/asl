use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_ore::{Powers, Sum, Support};

use super::period;

mod connection;
mod series;

#[derive(Clone, Copy)]
pub(crate) enum Operator {
    PicardFuchs,
    Clausen,
}

impl Operator {
    fn rank(self) -> usize {
        match self {
            Self::PicardFuchs => 3,
            Self::Clausen => 2,
        }
    }

    fn with_rank(rank: usize) -> Result<Self> {
        match rank {
            2 => Ok(Self::Clausen),
            3 => Ok(Self::PicardFuchs),
            _ => Err(Error::plain("Dwork mirror operator rank must be 2 or 3")),
        }
    }

    fn factors(self) -> (i128, &'static [i128]) {
        match self {
            Self::PicardFuchs => (4, &[1, 2, 3]),
            Self::Clausen => (8, &[1, 3]),
        }
    }

    // Coefficients of P(theta) in theta^rank - z P(theta).
    fn coefficients(self) -> [i128; 4] {
        let (scale, factors) = self.factors();
        let mut coefficients = [4, 0, 0, 0];
        for &offset in factors {
            coefficients = std::array::from_fn(|i| {
                let shifted = i.checked_sub(1).map_or(0, |j| coefficients[j]);
                offset * coefficients[i] + scale * shifted
            });
        }
        coefficients
    }

    fn terms(self) -> impl Iterator<Item = (Powers, i128)> {
        std::iter::once((
            Powers {
                q: 0,
                theta: self.rank(),
            },
            1,
        ))
        .chain(
            self.coefficients()
                .into_iter()
                .enumerate()
                .take(self.rank() + 1)
                .map(|(theta, coefficient)| (Powers { q: 1, theta }, -coefficient)),
        )
    }

    pub(crate) fn value(self, model: &mut Model) -> Result<Value> {
        let mut sum = Sum::default();
        for (powers, coefficient) in self.terms() {
            sum.add(model, powers, coefficient)?;
        }
        sum.finish(model)
    }
}

pub(crate) enum Report {
    Periods(usize),
    MirrorMap(usize),
    Monodromy(usize, i128),
    Connection,
    Congruence(usize, i128, usize),
}

impl Report {
    pub(crate) fn run(self, model: &mut Model) -> Result<String> {
        match self {
            Self::Periods(degree) => series::period_report(model, degree),
            Self::MirrorMap(degree) => series::mirror_map(model, degree),
            Self::Monodromy(degree, shift) => series::monodromy(model, degree, shift),
            Self::Connection => connection::report(),
            Self::Congruence(prime, psi, precision) => {
                let p = super::quartic_prime(model, prime, psi)?;
                period::mirror_congruence(model, p, psi, precision)
            }
        }
    }
}

pub(crate) fn apply(
    model: &mut Model,
    rank: usize,
    solution: usize,
    degree: usize,
    operator: &Support,
) -> Result<String> {
    let family = Operator::with_rank(rank)?;
    let jet = series::apply(model, family, solution, degree, operator)?;
    Ok(format!(
        "rvDworkMirrorApply[{rank}, {solution}, {degree}] = {jet}"
    ))
}

fn integrity() -> Error {
    Error::plain("Dwork mirror periods failed an exact operator or series identity")
}
