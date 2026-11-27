use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::wide_integer::Natural;

use super::distribution::{bound, Distribution};
use super::Field;

pub(super) fn absolute_trace(model: &mut Model, field: &Field, value: u32) -> Result<u32> {
    let p = field.characteristic();
    let mut q = field.order();
    let mut current = value;
    let mut trace = 0;
    while q > 1 {
        trace = field.add(trace, current);
        q /= p;
        if q > 1 {
            current = frobenius(model, field, current)?;
        }
    }
    if trace >= p {
        return Err(Error::plain(
            "partition quiver absolute trace did not descend to the prime field",
        ));
    }
    Ok(trace)
}

fn frobenius(model: &mut Model, field: &Field, mut value: u32) -> Result<u32> {
    let mut power = field.characteristic();
    let mut out = 1;
    while power != 0 {
        model.spend(2)?;
        if power & 1 != 0 {
            out = field.mul(out, value);
        }
        power >>= 1;
        value = field.mul(value, value);
    }
    Ok(out)
}

impl Distribution {
    fn phases(model: &mut Model, field: &Field) -> Result<Vec<usize>> {
        (0..field.order())
            .map(|x| Ok(absolute_trace(model, field, x)? as usize))
            .collect()
    }

    fn frequency(&self, field: &Field, phases: &[usize], lambda: u32) -> Result<Vec<u128>> {
        let mut frequency = vec![0u128; field.characteristic() as usize];
        for (value, &count) in self.framed.iter().enumerate() {
            let phase = phases[field.mul(lambda, value as u32) as usize];
            frequency[phase] = frequency[phase].checked_add(count).ok_or_else(bound)?;
        }
        Ok(frequency)
    }

    pub(super) fn fourier(&self, model: &mut Model, field: &Field, lambda: u32) -> Result<Value> {
        let q = field.order() as usize;
        model.spend(q)?;
        let phases = Self::phases(model, field)?;
        let mut frequency = self.frequency(field, &phases, lambda)?;
        let last = frequency.pop().unwrap();
        // Phi_p = 1+z+...+z^(p−1). Subtract first in u128 so this tests
        // the exact signed output range after cyclotomic cancellation.
        let coefficients = frequency
            .into_iter()
            .map(|n| {
                Natural::<8>::from_small(n.abs_diff(last))
                    .signed(n < last)
                    .ok_or_else(|| {
                        Error::plain("partition quiver Fourier coefficient exceeds signed 128 bits")
                    })
            })
            .collect::<Result<Vec<_>>>()?;
        crate::rv_cyclotomic::persist(model, coefficients)
    }

    pub(super) fn inversion(&self, model: &mut Model, field: &Field, name: &str) -> Result<String> {
        let q = field.order();
        let p = field.characteristic() as usize;
        model.spend((q as usize).pow(2) * (p + 1))?;
        let phases = Self::phases(model, field)?;
        let frequencies = (0..q)
            .map(|lambda| self.frequency(field, &phases, lambda))
            .collect::<Result<Vec<_>>>()?;
        for (value, &count) in self.framed.iter().enumerate() {
            let mut output = vec![Natural::<8>::ZERO; p];
            for (lambda, frequency) in frequencies.iter().enumerate() {
                let phase = phases[field.mul(lambda as u32, value as u32) as usize];
                for (i, &n) in frequency.iter().enumerate() {
                    output[(i + p - phase) % p]
                        .add_multiple(Natural::from_small(n), 1)
                        .ok_or_else(bound)?;
                }
            }
            let mut expected = output[1];
            expected
                .add_multiple(Natural::from_small(count), q)
                .ok_or_else(bound)?;
            if output[0] != expected || output[1..].iter().any(|&n| n != output[1]) {
                return Err(Error::plain(
                    "partition quiver Fourier inversion failed in the prime cyclotomic ring",
                ));
            }
        }
        Ok(format!(
            "rvPartitionQuiverFourierCheck[{name}] = q={q}, recovered={:?}",
            self.framed
        ))
    }
}
