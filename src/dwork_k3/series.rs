use crate::elliptic_modular::decimal;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_ore::gcd;
use crate::rv_qp_quad::Wide;

use super::{numerator_factors, Operation};

const MAX_DEGREE: usize = 32;
type Integer = Wide<40>;

pub(super) fn report(
    model: &mut Model,
    prime: Prime,
    psi: i128,
    residue: u32,
    degree: usize,
    operation: Operation,
) -> Result<String> {
    if !(1..=MAX_DEGREE).contains(&degree) {
        return Err(Error::plain(format!(
            "Dwork formal series degree must be in 1..={MAX_DEGREE}"
        )));
    }
    let periods = periods(model, psi, degree)?;
    if matches!(operation, Operation::Log) {
        return Ok(logarithm(&periods));
    }

    let target = match operation {
        Operation::Group => Target::Addition,
        Operation::Inverse => Target::Endomorphism(-1),
        Operation::PSeries => Target::Endomorphism(i128::from(prime.get())),
        _ => unreachable!(),
    };
    let law = solve(model, &periods, target)?;
    let mut result = law.report();
    if matches!(operation, Operation::PSeries) {
        let supersingular_fermat = residue == 0 && prime.get() % 4 == 3;
        result.push_str(&format!(
            ", mod {}={}, {}",
            prime.get(),
            law.reduce(prime).report(),
            law.height(prime, supersingular_fermat)?
        ));
    }
    Ok(result)
}

fn periods(model: &mut Model, psi: i128, degree: usize) -> Result<Vec<Integer>> {
    let z = Integer::from_i128(psi).multiply(-4).map_err(wide_error)?;
    let square = z.multiply_wide(z).map_err(wide_error)?;
    let fourth = square.multiply_wide(square).map_err(wide_error)?;
    let tails = [
        Integer::from_i128(1),
        z,
        square,
        square.multiply_wide(z).map_err(wide_error)?,
    ];
    let mut periods = vec![Integer::from_i128(0)];
    for n in 0..degree as u64 {
        let mut sum = Integer::from_i128(0);
        let mut coefficient = Integer::from_i128(1);
        for j in 0..=n / 4 {
            model.spend(1)?;
            sum = sum
                .multiply_wide(fourth)
                .and_then(|value| value.add(coefficient))
                .map_err(wide_error)?;
            if j != n / 4 {
                for factor in numerator_factors(n, j) {
                    coefficient = coefficient
                        .multiply(i128::from(factor))
                        .map_err(wide_error)?;
                }
                for _ in 0..4 {
                    let (value, remainder) = coefficient.magnitude.divide_u32(j as u32 + 1);
                    if remainder != 0 {
                        return Err(integrity());
                    }
                    coefficient.magnitude = value;
                }
            }
        }
        periods.push(
            sum.multiply_wide(tails[(n % 4) as usize])
                .map_err(wide_error)?,
        );
    }
    Ok(periods)
}

fn logarithm(periods: &[Integer]) -> String {
    let terms = periods
        .iter()
        .enumerate()
        .skip(1)
        .filter(|(_, value)| !value.is_zero())
        .map(|(degree, value)| {
            let divisor = gcd(
                u128::from(value.magnitude.residue(degree as u32)),
                degree as u128,
            ) as u32;
            let (magnitude, _) = value.magnitude.divide_u32(divisor);
            let numerator = decimal(Integer {
                magnitude,
                ..*value
            });
            let denominator = degree / divisor as usize;
            if denominator == 1 {
                format!("{degree}: {numerator}")
            } else {
                format!("{degree}: {numerator}/{denominator}")
            }
        });
    format!("[{}]", terms.collect::<Vec<_>>().join(", "))
}

#[derive(Clone, Copy)]
enum Target {
    Addition,
    Endomorphism(i128),
}

impl Target {
    fn coefficient(self, period: Integer, degree: usize, y: usize) -> Result<Integer> {
        match self {
            Self::Addition if y == 0 || y == degree => Ok(period),
            Self::Addition => Ok(Integer::from_i128(0)),
            Self::Endomorphism(scale) => period.multiply(scale).map_err(wide_error),
        }
    }
}

struct Series {
    degree: usize,
    bivariate: bool,
    coefficients: Vec<Integer>,
}

impl Series {
    fn new(degree: usize, bivariate: bool) -> Self {
        let size = if bivariate {
            (degree + 1) * (degree + 2) / 2
        } else {
            degree + 1
        };
        Self {
            degree,
            bivariate,
            coefficients: vec![Integer::from_i128(0); size],
        }
    }

    fn index(&self, degree: usize, y: usize) -> usize {
        if self.bivariate {
            degree * (degree + 1) / 2 + y
        } else {
            degree
        }
    }

    fn y_max(&self, degree: usize) -> usize {
        if self.bivariate {
            degree
        } else {
            0
        }
    }

    fn get(&self, degree: usize, y: usize) -> Integer {
        self.coefficients[self.index(degree, y)]
    }

    fn set(&mut self, degree: usize, y: usize, value: Integer) {
        let index = self.index(degree, y);
        self.coefficients[index] = value;
    }

    fn report(&self) -> String {
        let mut terms = Vec::new();
        for degree in 1..=self.degree {
            for y in 0..=self.y_max(degree) {
                let value = self.get(degree, y);
                if value.is_zero() {
                    continue;
                }
                let value = decimal(value);
                terms.push(if self.bivariate {
                    format!("({}, {}): {value}", degree - y, y)
                } else {
                    format!("{degree}: {value}")
                });
            }
        }
        format!("[{}]", terms.join(", "))
    }

    fn reduce(&self, p: Prime) -> Self {
        let mut reduced = Self::new(self.degree, self.bivariate);
        for (&coefficient, residue) in self.coefficients.iter().zip(&mut reduced.coefficients) {
            let value = coefficient.magnitude.residue(p.get());
            let value = if coefficient.negative {
                p.negate(value)
            } else {
                value
            };
            *residue = Integer::from_i128(i128::from(value));
        }
        reduced
    }

    fn height(&self, p: Prime, supersingular_fermat: bool) -> Result<String> {
        let prime = p.get() as usize;
        let mut degree = prime;
        let mut height = 1;
        for (n, coefficient) in self.coefficients.iter().enumerate().skip(1) {
            if n == degree && coefficient.magnitude.residue(p.get()) != 0 {
                if supersingular_fermat {
                    return Err(integrity());
                }
                return Ok(format!("height={height}"));
            }
            if coefficient.magnitude.residue(p.get()) != 0 {
                return Err(integrity());
            }
            if n == degree {
                degree = degree.saturating_mul(prime);
                height += 1;
            }
        }
        Ok(if supersingular_fermat {
            "height=infinite (supersingular Fermat quartic)".into()
        } else {
            format!("height lower bound={height}")
        })
    }
}

// Clearing the denominators in the logarithm gives D*l(X). Its linear
// term D*X means each homogeneous degree is solved by one exact division.
// Stienstra's formal group is over Z for this flat integral hypersurface.
fn solve(model: &mut Model, periods: &[Integer], target: Target) -> Result<Series> {
    let degree = periods.len() - 1;
    let denominator = (2..=degree as u128).fold(1, |d, n| d / gcd(d, n) * n);
    let log = periods
        .iter()
        .enumerate()
        .skip(1)
        .map(|(n, c)| {
            c.multiply((denominator / n as u128) as i128)
                .map_err(wide_error)
        })
        .collect::<Result<Vec<_>>>()?;
    let bivariate = matches!(target, Target::Addition);
    let mut powers = (0..=degree)
        .map(|_| Series::new(degree, bivariate))
        .collect::<Vec<_>>();
    for d in 1..=degree {
        for k in 2..=d {
            let (previous, rest) = powers.split_at_mut(k);
            for y in 0..=previous[1].y_max(d) {
                let coefficient =
                    product_coefficient(model, &previous[k - 1], &previous[1], k - 1, d, y)?;
                rest[0].set(d, y, coefficient);
            }
        }
        for y in 0..=powers[1].y_max(d) {
            let mut coefficient = target.coefficient(log[d - 1], d, y)?;
            for k in 2..=d {
                let product = powers[k]
                    .get(d, y)
                    .multiply_wide(log[k - 1])
                    .map_err(wide_error)?;
                coefficient = coefficient.add(product.negate()).map_err(wide_error)?;
            }
            let (magnitude, remainder) = coefficient.magnitude.divide(model, denominator)?;
            if remainder != 0 {
                return Err(integrity());
            }
            powers[1].set(
                d,
                y,
                Integer {
                    magnitude,
                    ..coefficient
                },
            );
        }
    }
    Ok(powers.swap_remove(1))
}

fn product_coefficient(
    model: &mut Model,
    power: &Series,
    law: &Series,
    order: usize,
    degree: usize,
    y: usize,
) -> Result<Integer> {
    let mut sum = Integer::from_i128(0);
    for d in order..degree {
        let other_degree = degree - d;
        let low = y.saturating_sub(other_degree);
        let high = y.min(d);
        model.spend(high + 1 - low)?;
        for j in low..=high {
            let (a, b) = (power.get(d, j), law.get(other_degree, y - j));
            if !a.is_zero() && !b.is_zero() {
                sum = a
                    .multiply_wide(b)
                    .and_then(|term| sum.add(term))
                    .map_err(wide_error)?;
            }
        }
    }
    Ok(sum)
}

fn wide_error(_: Error) -> Error {
    Error::plain("Dwork formal series exceeds the checked exact coefficient range")
}

fn integrity() -> Error {
    Error::plain("Dwork formal series failed its exact integrality or height check")
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;

    const CUTOFF: usize = 7;
    type Polynomial = BTreeMap<[usize; 3], i128>;

    fn symbol(axis: usize) -> Polynomial {
        let mut exponent = [0; 3];
        exponent[axis] = 1;
        BTreeMap::from([(exponent, 1)])
    }

    fn times(a: &Polynomial, b: &Polynomial) -> Polynomial {
        let mut product = Polynomial::new();
        for (i, a) in a {
            for (j, b) in b {
                let exponent = std::array::from_fn(|k| i[k] + j[k]);
                if exponent.iter().sum::<usize>() <= CUTOFF {
                    *product.entry(exponent).or_default() += a * b;
                }
            }
        }
        product.retain(|_, c| *c != 0);
        product
    }

    fn powers(a: &Polynomial) -> Vec<Polynomial> {
        let mut powers = vec![BTreeMap::from([([0; 3], 1)])];
        for _ in 0..CUTOFF {
            powers.push(times(powers.last().unwrap(), a));
        }
        powers
    }

    fn substitute(law: &Series, a: &Polynomial, b: &Polynomial) -> Polynomial {
        let ap = powers(a);
        let bp = powers(b);
        let mut result = Polynomial::new();
        for n in 1..=law.degree {
            for j in 0..=law.y_max(n) {
                let coefficient = law.get(n, j).signed().unwrap();
                if coefficient != 0 {
                    for (exponents, product) in times(&ap[n - j], &bp[j]) {
                        *result.entry(exponents).or_default() += product * coefficient;
                    }
                }
            }
        }
        result.retain(|_, c| *c != 0);
        result
    }

    #[test]
    fn exact_three_variable_associativity_inverse_and_repeated_addition() {
        let mut model = Model::new();
        for (p, psi) in [(5, 0), (7, 3), (13, 2)] {
            let periods =
                periods(&mut model, psi, CUTOFF).unwrap_or_else(|error| panic!("{error}"));
            let law = solve(&mut model, &periods, Target::Addition)
                .unwrap_or_else(|error| panic!("{error}"));
            let inverse = solve(&mut model, &periods, Target::Endomorphism(-1))
                .unwrap_or_else(|error| panic!("{error}"));
            let mul_p = solve(&mut model, &periods, Target::Endomorphism(p))
                .unwrap_or_else(|error| panic!("{error}"));
            let (a, b, c) = (symbol(0), symbol(1), symbol(2));
            let zero = Polynomial::new();
            assert_eq!(substitute(&law, &a, &zero), a);
            assert_eq!(substitute(&law, &a, &b), substitute(&law, &b, &a));
            let left = substitute(&law, &substitute(&law, &a, &b), &c);
            let right = substitute(&law, &a, &substitute(&law, &b, &c));
            assert_eq!(left, right);
            assert_eq!(substitute(&law, &a, &substitute(&inverse, &a, &zero)), zero);
            let mut repeated = zero.clone();
            for _ in 0..p {
                repeated = substitute(&law, &repeated, &a);
            }
            assert_eq!(repeated, substitute(&mul_p, &a, &zero));
        }
    }
}
