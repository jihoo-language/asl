use std::collections::BTreeMap;
use std::fmt::Write;

use crate::error::{Error, Result};
use crate::generalized_kummer::hodge;
use crate::model::{Model, MAX_SHOWN_BYTES};
use crate::rv_ore::gcd;
use crate::syntax::Label;
use crate::wide_integer::Natural;

use super::{hilbert::SurfaceKind, rational_sparse::decimal};

pub(crate) enum Command {
    Declare(Label, Vec<(u128, u128)>),
    Data(Label),
    SurfaceHodge(Label, SurfaceKind, u128, u128),
    GeneralizedHodge(Label, u128, u128, u128),
}

pub(crate) struct Partition {
    parts: Box<[(u128, u128)]>,
    pub(crate) weight: u128,
    pub(crate) length: u128,
    pub(crate) gcd: u128,
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<Option<String>> {
        let report = match self {
            Self::Declare(name, parts) => {
                let shape = Partition::new(model, parts)?;
                model.define_hilbert_stratum(name, shape)?;
                return Ok(None);
            }
            Self::Data(name) => {
                let shape = model.hilbert_stratum(&name)?;
                shape.report(model, &name.text)?
            }
            Self::SurfaceHodge(name, kind, p, q) => {
                let shape = model.hilbert_stratum(&name)?;
                let value = hodge::stratum(model, &shape, kind, None, p, q)?;
                let operation = match kind {
                    SurfaceKind::K3 => "rvK3StratumHodge",
                    SurfaceKind::Abelian => "rvJacStratumHodge",
                };
                format!("{operation}[{}, {p}, {q}] = {value}", name.text)
            }
            Self::GeneralizedHodge(name, h, p, q) => {
                let shape = model.hilbert_stratum(&name)?;
                let value = hodge::stratum(model, &shape, SurfaceKind::Abelian, Some(h), p, q)?;
                format!(
                    "rvGenKummerStratumHodge[{}, {h}, {p}, {q}] = {value}",
                    name.text
                )
            }
        };
        Ok(Some(report))
    }
}

impl Partition {
    fn new(model: &mut Model, parts: Vec<(u128, u128)>) -> Result<Self> {
        let mut combined = BTreeMap::<u128, u128>::new();
        let (mut weight, mut length, mut divisor) = (0u128, 0u128, 0);
        for (size, count) in parts {
            model.spend(1 + combined.len().max(1).ilog2() as usize)?;
            if size == 0 {
                return Err(Error::plain(
                    "Hilbert–Chow partition part size must be positive",
                ));
            }
            if count == 0 {
                continue;
            }
            weight = size
                .checked_mul(count)
                .and_then(|part| weight.checked_add(part))
                .ok_or_else(weight_overflow)?;
            length = length.checked_add(count).ok_or_else(weight_overflow)?;
            divisor = gcd(divisor, size);
            let entry = combined.entry(size).or_default();
            *entry = entry.checked_add(count).ok_or_else(weight_overflow)?;
        }
        Ok(Self {
            parts: combined.into_iter().collect(),
            weight,
            length,
            gcd: divisor,
        })
    }

    pub(crate) fn parts(&self) -> &[(u128, u128)] {
        &self.parts
    }

    pub(crate) fn shift(&self) -> u128 {
        self.weight - self.length
    }

    fn report(&self, model: &mut Model, name: &str) -> Result<String> {
        let shift = self.shift();
        let codimension = Natural::<5>::from_small(shift)
            .multiply_u128(2)
            .expect("twice an unsigned 128-bit integer fits in 129 bits");
        let mut report = format!(
            "rvHilbStratumData[{name}] = n={}, L={}, e={}, fiber dimension={shift}, codimension={}, intersection multiplier={}",
            self.weight,
            self.length,
            self.gcd,
            decimal(codimension),
            if shift.is_multiple_of(2) { "" } else { "-" },
        );
        model.spend(self.parts.len())?;
        let exact = self.parts.iter().filter(|&&(size, _)| size > 1).try_fold(
            1u128,
            |product, &(size, count)| {
                let power = u32::try_from(count).ok()?;
                product.checked_mul(size.checked_pow(power)?)
            },
        );
        if let Some(value) = exact {
            write!(report, "{value}").unwrap();
        } else {
            let mut factors = self.parts.iter().filter(|&&(size, _)| size > 1).peekable();
            while let Some(&(size, count)) = factors.next() {
                if count == 1 {
                    write!(report, "{size}").unwrap();
                } else {
                    write!(report, "{size}^{count}").unwrap();
                }
                if factors.peek().is_some() {
                    report.push('*');
                }
                if report.len() > MAX_SHOWN_BYTES {
                    return Err(Error::plain(
                        "Hilbert–Chow stratum report exceeds the rendering limit",
                    ));
                }
            }
        }
        Ok(report)
    }
}

fn weight_overflow() -> Error {
    Error::plain("Hilbert–Chow partition weight exceeds unsigned 128-bit range")
}
