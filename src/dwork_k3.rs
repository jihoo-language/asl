use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;

pub(crate) mod mirror;
mod period;
mod series;

pub(crate) use period::closed_unit;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Log,
    Group,
    Inverse,
    PSeries,
    Unit,
}

impl Operation {
    pub(crate) fn spelling(self) -> &'static str {
        match self {
            Self::Log => "rvDworkK3Log",
            Self::Group => "rvDworkK3Group",
            Self::Inverse => "rvDworkK3Inverse",
            Self::PSeries => "rvDworkK3PSeries",
            Self::Unit => "rvDworkK3Unit",
        }
    }
}

pub(crate) struct Query {
    pub(crate) operation: Operation,
    pub(crate) prime: usize,
    pub(crate) psi: i128,
    pub(crate) cutoff: usize,
}

impl Query {
    pub(crate) fn run(self, model: &mut Model) -> Result<String> {
        let Self {
            operation,
            prime,
            psi,
            cutoff,
        } = self;
        let p = quartic_prime(model, prime, psi)?;
        let residue = psi.rem_euclid(i128::from(p.get())) as u32;
        let report = match operation {
            Operation::Unit => period::unit_root(model, p, psi, cutoff)?,
            _ => series::report(model, p, psi, residue, cutoff, operation)?,
        };
        Ok(format!(
            "{}[{prime}, {psi}, {cutoff}] = {report}",
            operation.spelling()
        ))
    }
}

fn quartic_prime(model: &mut Model, prime: usize, psi: i128) -> Result<Prime> {
    if prime < 3 || prime.is_multiple_of(2) {
        return Err(odd_characteristic());
    }
    let p = Prime::new(model, prime)?;
    let residue = psi.rem_euclid(i128::from(p.get())) as u32;
    let square = p.mul(residue, residue);
    if p.mul(square, square) == 1 {
        return Err(singular());
    }
    Ok(p)
}

fn numerator_factors(n: u64, j: u64) -> [u64; 4] {
    let k = n - 4 * j;
    [k, k - 1, k - 2, k - 3]
}

fn odd_characteristic() -> Error {
    Error::plain("Dwork quartics require an odd prime characteristic")
}

fn singular() -> Error {
    Error::plain("singular Dwork quartic: psi^4 is 1 in the residue field")
}
