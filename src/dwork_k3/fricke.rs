use crate::error::{Error, Result};
use crate::model::Model;

mod eta;
mod fiber;
mod identity;

pub(crate) use fiber::{base_report, define, fiber_report, map, Invariant};

pub(crate) enum Report {
    Series(usize),
    Identity,
    Compare(usize, i128, usize),
}

impl Report {
    pub(crate) fn run(self, model: &mut Model) -> Result<String> {
        match self {
            Self::Series(degree) => eta::report(model, degree),
            Self::Identity => identity::report(model),
            Self::Compare(prime, psi, precision) => fiber::compare(model, prime, psi, precision),
        }
    }
}

fn certificate() -> Error {
    Error::plain("Dwork Fricke correspondence failed an exact identity")
}

fn arithmetic() -> Error {
    Error::plain("Dwork Fricke exact rational arithmetic exceeds the signed 128-bit range")
}
