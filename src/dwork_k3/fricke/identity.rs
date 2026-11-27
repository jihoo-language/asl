use std::collections::BTreeMap;

use crate::elliptic_modular::Polynomial;
use crate::error::Result;
use crate::model::Model;
use crate::odd_cobar::Prime;

use super::{arithmetic, certificate};

type Coefficients = [(i32, i128)];

#[derive(Clone, Eq, PartialEq)]
struct Laurent(BTreeMap<i32, i128>);

impl Laurent {
    fn new(terms: &Coefficients) -> Self {
        Self(terms.iter().copied().collect())
    }

    fn add_into(&mut self, other: &Self, scale: i128) -> Result<()> {
        for (&exponent, &coefficient) in &other.0 {
            let entry = self.0.entry(exponent).or_default();
            *entry = coefficient
                .checked_mul(scale)
                .and_then(|term| entry.checked_add(term))
                .ok_or_else(arithmetic)?;
        }
        self.0.retain(|_, c| *c != 0);
        Ok(())
    }

    fn product(&self, other: &Self) -> Result<Self> {
        let mut product = Self::new(&[]);
        for (&i, &u) in &self.0 {
            for (&j, &v) in &other.0 {
                let slot = product.0.entry(i + j).or_default();
                *slot = u
                    .checked_mul(v)
                    .and_then(|c| slot.checked_add(c))
                    .ok_or_else(arithmetic)?;
            }
        }
        product.0.retain(|_, c| *c != 0);
        Ok(product)
    }

    fn powers(&self) -> Result<[Self; 4]> {
        let square = self.product(self)?;
        let cube = square.product(self)?;
        Ok([Self::new(&[(0, 1)]), self.clone(), square, cube])
    }

    fn involution(&self) -> Result<Self> {
        let mut transformed = BTreeMap::new();
        for (&exponent, &coefficient) in &self.0 {
            let power = 4096i128
                .checked_pow(exponent.unsigned_abs())
                .ok_or_else(arithmetic)?;
            let value = if exponent >= 0 {
                coefficient.checked_mul(power).ok_or_else(arithmetic)?
            } else {
                if coefficient % power != 0 {
                    return Err(certificate());
                }
                coefficient / power
            };
            transformed.insert(-exponent, value);
        }
        Ok(Self(transformed))
    }
}

pub(super) fn report(model: &mut Model) -> Result<String> {
    let j1 = Laurent::new(&[
        (-2, 256i128.pow(3)),
        (-1, 3 * 256i128.pow(2)),
        (0, 768),
        (1, 1),
    ]);
    let j2 = Laurent::new(&[(-1, 16i128.pow(3)), (0, 3 * 16i128.pow(2)), (1, 48), (2, 1)]);
    let v = Laurent::new(&[(-1, 4096), (0, 128), (1, 1)]);
    let constant = Laurent::new(&[(0, 1)]);
    let mut sum = j1.clone();
    sum.add_into(&j2, 1)?;
    let mut expected_sum = v.product(&v)?;
    expected_sum.add_into(&v, -207)?;
    expected_sum.add_into(&constant, 3456)?;
    let mut expected_product = v.clone();
    expected_product.add_into(&constant, 144)?;
    if sum != expected_sum
        || j1.product(&j2)? != expected_product.powers()?[3]
        || j1.involution()? != j2
        || j2.involution()? != j1
        || v.involution()? != v
    {
        return Err(certificate());
    }

    let phi = Polynomial::new(model, Prime::validated(2))?;
    let (jp1, jp2) = (j1.powers()?, j2.powers()?);
    let mut value = Laurent::new(&[]);
    for (i, p1) in jp1.iter().enumerate() {
        for (j, p2) in jp2.iter().enumerate() {
            let term = phi
                .wide_coefficient::<5>(i, j)?
                .signed()
                .ok_or_else(arithmetic)?;
            if term != 0 {
                value.add_into(&p1.product(p2)?, term)?;
            }
        }
    }
    if !value.0.is_empty() {
        return Err(certificate());
    }
    Ok("rvDworkFrickeIdentity = j1=(t+256)^3/t^2, j2=(t+16)^3/t, v=t+128+4096/t, J^2-(v^2-207v+3456)J+(v+144)^3, certified W2 exchange and Phi_2(j1,j2)=0 over Q(t)".into())
}
