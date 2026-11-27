use crate::error::{Error, Result};
use crate::hw_algebra::{self as algebra, CoefficientRing, Matrix};
use crate::model::Model;
use crate::rv_newton::valuation_slopes;
use crate::rv_ore::gcd;
use crate::rv_qp_quad::mul_mod;

pub(super) struct Newton {
    pub(super) determinant: i64,
    pub(super) positive: Vec<Slope>,
}

#[derive(Clone, Copy)]
pub(super) struct Slope {
    pub(super) numerator: u128,
    denominator: u128,
    pub(super) count: usize,
}

impl Slope {
    pub(super) fn past_threshold(self, prime: u128) -> bool {
        self.numerator * (prime - 1) > self.denominator
    }

    pub(super) fn report(self, offset: usize) -> String {
        let numerator = self.numerator + (offset as u128) * self.denominator;
        if self.denominator == 1 {
            format!("{numerator}: {}", self.count)
        } else {
            format!("{numerator}/{}: {}", self.denominator, self.count)
        }
    }
}

pub(super) fn certify(
    model: &mut Model,
    q: u64,
    polynomial: &[i128],
    degree: u128,
    prime: u128,
    level: usize,
) -> Result<Newton> {
    if polynomial == [1] {
        return Ok(Newton {
            determinant: 0,
            positive: Vec::new(),
        });
    }
    let max = max_precision(prime);
    let mut precision = 1;
    loop {
        let modulus = prime.pow(precision);
        let ring = Residues { prime, modulus };
        let characteristic = shifted_characteristic(model, &ring, q, polynomial, degree, level)?;
        if let Some(determinant) = valuation(prime, *characteristic.last().unwrap()) {
            // The companion eigenvalues are q^r times the normalized
            // eigenvalues. Since q is an ell-unit, their shifted valuations
            // are unchanged. Unresolved coefficients are >= precision,
            // strictly above the chord between the now-known endpoints.
            let valuations = characteristic
                .into_iter()
                .rev()
                .map(|c| valuation(prime, c));
            let mut positive = Vec::new();
            let mut total = 0;
            for (numerator, count) in valuation_slopes(model, valuations)? {
                let v = u128::try_from(numerator).map_err(|_| inconsistent())?;
                total += v;
                if v != 0 {
                    let divisor = gcd(v, count as u128);
                    positive.push(Slope {
                        numerator: v / divisor,
                        denominator: count as u128 / divisor,
                        count,
                    });
                }
            }
            if total != determinant as u128 {
                return Err(inconsistent());
            }
            return Ok(Newton {
                determinant,
                positive,
            });
        }
        if precision == max {
            return Err(Error::plain(
                "Kummer K3 Iwasawa Newton polygon exceeds the checked ell-adic precision",
            ));
        }
        precision = (precision * 2).min(max);
    }
}

pub(super) fn valuation(prime: u128, mut value: u128) -> Option<i64> {
    if value == 0 {
        return None;
    }
    let mut order = 0;
    while value.is_multiple_of(prime) {
        value /= prime;
        order += 1;
    }
    Some(order)
}

fn max_precision(prime: u128) -> u32 {
    let maximum = i128::MAX as u128 / prime;
    let mut value = 1;
    let mut precision = 0;
    while value <= maximum {
        value *= prime;
        precision += 1;
    }
    precision
}

fn shifted_characteristic(
    model: &mut Model,
    ring: &Residues,
    q: u64,
    polynomial: &[i128],
    degree: u128,
    level: usize,
) -> Result<Vec<u128>> {
    let n = polynomial.len() - 1;
    let mut entries = vec![0; n * n];
    for i in 0..n {
        let coefficient = polynomial[n - i].rem_euclid(ring.modulus as i128) as u128;
        entries[i * n + n - 1] = ring.subtract(model, &0, &coefficient)?;
        if i != 0 {
            entries[i * n + i - 1] = 1;
        }
    }
    let companion = Matrix::from_entries(n, entries)?;
    let mut powered = matrix_power(model, ring, companion, degree)?;
    let mut q_power = ring.power(model, u128::from(q) % ring.modulus, degree)?;
    for _ in 0..level {
        powered = matrix_power(model, ring, powered, ring.prime)?;
        q_power = ring.power(model, q_power, ring.prime)?;
    }
    let mut shifted = powered.entries().to_vec();
    for i in 0..n {
        let entry = &mut shifted[i * (n + 1)];
        *entry = ring.subtract(model, entry, &q_power)?;
    }
    algebra::characteristic(model, ring, &Matrix::from_entries(n, shifted)?)
}

fn matrix_power(
    model: &mut Model,
    ring: &Residues,
    mut matrix: Matrix<u128>,
    mut exponent: u128,
) -> Result<Matrix<u128>> {
    let mut output = algebra::identity(model, ring, matrix.size())?;
    while exponent != 0 {
        if exponent & 1 != 0 {
            output = algebra::multiply(model, ring, &output, &matrix)?;
        }
        exponent >>= 1;
        if exponent != 0 {
            matrix = algebra::multiply(model, ring, &matrix, &matrix)?;
        }
    }
    Ok(output)
}

struct Residues {
    prime: u128,
    modulus: u128,
}

impl Residues {
    fn power(&self, model: &mut Model, mut term: u128, mut exponent: u128) -> Result<u128> {
        let mut output = 1;
        while exponent != 0 {
            if exponent & 1 != 0 {
                output = self.multiply(model, &output, &term)?;
            }
            exponent >>= 1;
            if exponent != 0 {
                term = self.multiply(model, &term, &term)?;
            }
        }
        Ok(output)
    }
}

impl CoefficientRing for Residues {
    type Element = u128;

    fn cells(&self) -> usize {
        4
    }
    fn zero(&self) -> u128 {
        0
    }
    fn one(&self) -> u128 {
        1
    }
    fn is_zero(&self, value: &u128) -> bool {
        *value == 0
    }
    fn is_unit(&self, value: &u128) -> bool {
        !value.is_multiple_of(self.prime)
    }

    fn add(&self, model: &mut Model, a: &u128, b: &u128) -> Result<u128> {
        model.spend(1)?;
        Ok((a + b) % self.modulus)
    }

    fn subtract(&self, model: &mut Model, a: &u128, b: &u128) -> Result<u128> {
        model.spend(1)?;
        Ok((self.modulus + a - b) % self.modulus)
    }

    fn multiply(&self, model: &mut Model, a: &u128, b: &u128) -> Result<u128> {
        mul_mod(model, *a, *b, self.modulus)
    }

    fn inverse_unit(&self, model: &mut Model, value: &u128) -> Result<u128> {
        if !self.is_unit(value) {
            return Err(Error::plain(
                "Kummer K3 Iwasawa coefficient is not an ell-adic unit",
            ));
        }
        self.power(model, *value, self.modulus - self.modulus / self.prime - 1)
    }
}

fn inconsistent() -> Error {
    Error::plain("Kummer K3 Iwasawa Newton slopes failed their determinant check")
}
