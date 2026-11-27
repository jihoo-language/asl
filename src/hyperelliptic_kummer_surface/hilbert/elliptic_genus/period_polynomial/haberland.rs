use crate::error::{Error, Result};
use crate::integer_smith::ColumnHermite;
use crate::lee_linear::Rational;
use crate::model::Model;
use crate::rv_ore::gcd;

use super::super::eigenpacket::{determinant, fraction};
use super::super::{decimal, Integer, Magnitude};
use super::{arithmetic, binomial_table, multiply, plus, ratio, Grade, Polynomial, Sign};

pub(crate) struct IntegralBasis(Vec<Polynomial>);

impl IntegralBasis {
    pub(crate) fn new(model: &mut Model, grade: Grade) -> Result<Self> {
        let rational = model.period_basis(grade)?;
        let dimension = rational.dimension();
        let w = grade.degree();
        let mut denominator = 1;
        for p in &rational.vectors {
            for coefficient in &p.0 {
                denominator = lcm(denominator, coefficient.denominator)?;
            }
        }
        // In free coordinates x, integrality means that each rational
        // coefficient row R_i pairs integrally with x. Compute the column
        // Hermite basis H of the lattice spanned by the integral rows D R_i.
        // The free rows of R are identity, so D*Z^d is already in this lattice.
        let mut lattice = ColumnHermite::from_diagonal(model, &vec![denominator; dimension])?;
        for i in 0..=w {
            let row = rational
                .vectors
                .iter()
                .map(|p| {
                    let coefficient = p.0[i];
                    let scale = denominator / coefficient.denominator;
                    let entry = coefficient
                        .numerator
                        .checked_mul(scale)
                        .ok_or_else(arithmetic)?;
                    Ok(entry.rem_euclid(denominator))
                })
                .collect::<Result<Vec<_>>>()?;
            if row.iter().any(|&entry| entry != 0) {
                lattice.insert(model, row)?;
                // Keep the intermediate columns small before the next row.
                lattice.canonicalize(model)?;
            }
        }

        // The dual condition H^t x in D*Z^d gives all integral Manin
        // polynomials, including those missed by clearing vectors separately.
        let mut basis = Vec::with_capacity(dimension);
        for position in 0..dimension {
            let coordinates = dual_column(model, lattice.columns(), denominator, position)?;
            let mut polynomial = Polynomial::zero(grade);
            for (p, coordinate) in rational.vectors.iter().zip(coordinates) {
                if coordinate != 0 {
                    let term = p.scale(model, Rational::integer(coordinate))?;
                    polynomial = polynomial.add(model, &term)?;
                }
            }
            if polynomial
                .0
                .iter()
                .any(|coefficient| coefficient.denominator != 1)
            {
                return Err(not_integral());
            }
            polynomial.check(model, grade)?;
            basis.push(polynomial);
        }
        model.allocate(1 + basis.len() * (1 + 4 * (w + 1)))?;
        Ok(Self(basis))
    }

    pub(crate) fn vector(&self, index: usize) -> Result<Polynomial> {
        self.0.get(index).cloned().ok_or_else(|| {
            Error::plain("integral period-polynomial basis index is outside this cusp space")
        })
    }
}

fn dual_column(
    model: &mut Model,
    columns: &[Vec<i128>],
    denominator: i128,
    index: usize,
) -> Result<Vec<i128>> {
    let mut solution = vec![0_i128; columns.len()];
    for (i, column) in columns.iter().enumerate().rev() {
        let mut rhs = if i == index { denominator } else { 0 };
        model.spend(1 + column.len() - i)?;
        for (&entry, &known) in column.iter().zip(&solution).skip(i + 1) {
            let term = entry.checked_mul(known).ok_or_else(arithmetic)?;
            rhs = rhs.checked_sub(term).ok_or_else(arithmetic)?;
        }
        if rhs % column[i] != 0 {
            return Err(not_integral());
        }
        solution[i] = rhs / column[i];
    }
    Ok(solution)
}

fn not_integral() -> Error {
    Error::plain("the parabolic Manin lattice dual has a nonintegral coordinate")
}

fn lcm(a: i128, b: i128) -> Result<i128> {
    let common = gcd(a as u128, b as u128) as i128;
    (a / common).checked_mul(b).ok_or_else(arithmetic)
}

pub(crate) fn pair(model: &mut Model, p: &Polynomial, q: &Polynomial) -> Result<Rational> {
    let w = p.0.len() - 1;
    let binomial = binomial_table(model, w)?;
    let translated = translation_difference(model, p, &binomial)?;
    contract(model, &translated, q, &binomial[w])
}

pub(crate) fn report(model: &mut Model, grade: Grade) -> Result<String> {
    let plus = model.integral_period_basis(grade)?;
    let minus = model.integral_period_basis(Grade {
        sign: Sign::Minus,
        ..grade
    })?;
    let dimension = plus.0.len();
    let binomial = binomial_table(model, grade.degree())?;
    let mut entries = Vec::with_capacity(dimension);
    let mut denominator = 1_i128;
    for p in &plus.0 {
        let translated = translation_difference(model, p, &binomial)?;
        let mut row = Vec::with_capacity(dimension);
        for q in &minus.0 {
            let entry = contract(model, &translated, q, &binomial[grade.degree()])?;
            denominator = lcm(denominator, entry.denominator)?;
            row.push(entry);
        }
        entries.push(row);
    }
    let cleared = entries
        .iter()
        .map(|row| {
            row.iter()
                .map(|entry| {
                    let value = entry
                        .numerator
                        .checked_mul(denominator / entry.denominator)
                        .ok_or_else(arithmetic)?;
                    Ok(Integer::from_i128(value))
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let factors = smith_factors(model, &cleared)?;
    let matrix: Vec<_> = entries
        .iter()
        .map(|row| format!("[{}]", display(row)))
        .collect();
    let integer_factors = factors.iter().copied().map(decimal).collect::<Vec<_>>();
    let elementary = factors
        .into_iter()
        .map(|value| fraction(model, value, Magnitude::from_small(denominator as u128)))
        .collect::<Result<Vec<_>>>()?;
    Ok(format!(
        "rvPeriodPairing[{}] = denominator={denominator}, matrix=[{}], smith=[{}], elementary=[{}]",
        grade.weight,
        matrix.join(", "),
        display(&integer_factors),
        display(&elementary)
    ))
}

fn smith_factors(model: &mut Model, matrix: &[Vec<Integer>]) -> Result<Vec<Integer>> {
    let n = matrix.len();
    let mut factors = Vec::with_capacity(n);
    let mut previous = Magnitude::from_small(1);
    for size in 1..=n {
        let mut divisor = Magnitude::ZERO;
        // The gcd of the size-i minors is the product of the first i
        // Smith factors. Bareiss avoids growing unimodular transforms.
        'minors: for rows in subsets(n, size) {
            for columns in subsets(n, size) {
                model.spend(1 + size * size)?;
                let minor = rows
                    .iter()
                    .map(|&i| columns.iter().map(|&j| matrix[i][j]).collect())
                    .collect();
                let determinant = determinant(model, minor)?;
                divisor = divisor.gcd(model, determinant.magnitude)?;
                // Each current minor is divisible by the previous
                // determinantal divisor, so we have reached the gcd.
                if divisor == previous {
                    break 'minors;
                }
            }
        }
        if divisor.is_zero() {
            return Err(Error::plain(
                "the parabolic Haberland cross-pairing is degenerate",
            ));
        }
        let magnitude = divisor.divide_exact(model, previous)?.ok_or_else(|| {
            Error::plain("Haberland Smith factors have inconsistent determinantal divisors")
        })?;
        factors.push(Integer {
            magnitude,
            negative: false,
        });
        previous = divisor;
    }
    Ok(factors)
}

fn subsets(n: usize, size: usize) -> impl Iterator<Item = Vec<usize>> {
    std::iter::successors(Some((0..size).collect::<Vec<_>>()), move |items| {
        let i = (0..size).rfind(|&i| items[i] < n - size + i)?;
        let mut next = items.clone();
        let start = next[i] + 1;
        for (j, item) in next[i..].iter_mut().enumerate() {
            *item = start + j;
        }
        Some(next)
    })
}

fn translation_difference(
    model: &mut Model,
    p: &Polynomial,
    binomial: &[Vec<i128>],
) -> Result<Vec<Rational>> {
    // For T and T^-1, X goes to X+Y and X-Y respectively.
    // The translation difference retains the odd binomial terms only.
    let mut result = vec![Rational::ZERO; p.0.len()];
    for (j, &coefficient) in p.0.iter().enumerate() {
        if coefficient.numerator == 0 || j == 0 {
            continue;
        }
        model.spend(1 + j / 2)?;
        let doubled = multiply(coefficient, Rational::integer(2))?;
        for (i, term) in result[..j]
            .iter_mut()
            .enumerate()
            .skip(1 - j % 2)
            .step_by(2)
        {
            let contribution = multiply(doubled, Rational::integer(binomial[j][i]))?;
            *term = plus(*term, contribution)?;
        }
    }
    Ok(result)
}

fn contract(
    model: &mut Model,
    translated: &[Rational],
    q: &Polynomial,
    binomial: &[i128],
) -> Result<Rational> {
    let mut result = Rational::ZERO;
    model.spend(translated.len())?;
    for (i, (&a, &b)) in translated.iter().zip(q.0.iter().rev()).enumerate() {
        if a.numerator == 0 || b.numerator == 0 {
            continue;
        }
        let inverse = ratio(if i.is_multiple_of(2) { 1 } else { -1 }, binomial[i])?;
        result = plus(result, multiply(multiply(a, inverse)?, b)?)?;
    }
    Ok(result)
}

fn display<T: std::fmt::Display>(entries: &[T]) -> String {
    entries
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boundary_pairs_to_zero_before_taking_the_quotient() {
        fn verify() -> Result<()> {
            let model = &mut Model::new();
            let grade = Grade::new(12, Sign::Plus)?;
            let mut boundary = Polynomial::zero(grade);
            boundary.0[0] = Rational::integer(-1);
            boundary.0[10] = Rational::integer(1);
            let odd = Polynomial::new(
                model,
                Grade {
                    sign: Sign::Minus,
                    ..grade
                },
                [0, 4, 0, -25, 0, 42, 0, -25, 0, 4, 0]
                    .into_iter()
                    .map(Rational::integer)
                    .collect(),
            )?;
            assert_eq!(pair(model, &boundary, &odd)?, Rational::ZERO);
            assert_eq!(pair(model, &odd, &boundary)?, Rational::ZERO);
            Ok(())
        }
        verify().unwrap_or_else(|error| panic!("{error}"));
    }
}
