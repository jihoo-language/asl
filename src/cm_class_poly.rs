use std::fmt::Write;

use crate::elliptic_modular;
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::{cm_neighbors, cm_ring::Order};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_qp_quad::Wide;
use crate::rv_unram::{self, Ring};
use crate::syntax::{Label, ParsedExpr};

const LIMBS: usize = 128;
type Integer = Wide<LIMBS>;

pub(crate) enum Statement {
    Seed(Label, i128),
    Lift(Label, Label, usize),
    Data(Label),
    Coefficient(Label, usize),
}

pub(crate) enum Expression {
    Evaluate(Label, Box<ParsedExpr>),
    ModularCoefficient(Label, Label, usize),
}

pub(crate) struct Polynomial {
    pub(crate) discriminant: i128,
    coefficients: Poly,
}

impl Polynomial {
    pub(crate) fn seed(model: &mut Model, discriminant: i128) -> Result<Self> {
        let singular_modulus =
            match discriminant {
                -3 => 0,
                -4 => 1728,
                -7 => -3375,
                -8 => 8000,
                -11 => -32768,
                -19 => -884736,
                -43 => -884736000,
                -67 => -147197952000,
                -163 => -262537412640768000,
                _ => return Err(Error::plain(
                    "CM class polynomial seed requires a class-number-one fundamental discriminant",
                )),
            };
        Self::persist(
            model,
            discriminant,
            Poly(vec![
                Integer::from_i128(-singular_modulus),
                Integer::from_i128(1),
            ]),
        )
    }

    pub(crate) fn lift(
        model: &mut Model,
        source: &Self,
        ascending: Option<&Self>,
        level: Prime,
    ) -> Result<Self> {
        let ell = i128::from(level.get());
        let discriminant = source.discriminant.checked_mul(ell * ell).ok_or_else(|| {
            Error::plain("descending CM class polynomial discriminant exceeds signed i128")
        })?;
        // Load the exact modular polynomial before iterating cyclic kernels.
        let modular = model.elliptic_modular_polynomial(level)?;
        let source_order = Order::new(source.discriminant)?;
        let classes = source_order.elements(model, source.degree())?;
        if classes.len() != source.degree() {
            return Err(Error::plain(
                "source CM polynomial degree differs from its class number",
            ));
        }
        let counts = cm_neighbors::counts(model, source_order.principal(), level)?;
        for class in classes {
            if cm_neighbors::counts(model, class, level)? != counts {
                return Err(Error::plain(
                    "CM cyclic neighbor counts disagree within one order",
                ));
            }
        }
        let target_degree = Order::new(discriminant)?.class_number(model)?;
        let down_power = multiplicity(source.degree(), counts.descending, target_degree)?;
        let horizontal = &source.coefficients;
        let upper = if counts.ascending == 0 {
            None
        } else {
            let d = source.discriminant / (ell * ell);
            let upper = ascending
                .filter(|value| value.discriminant == d)
                .ok_or_else(|| {
                    Error::plain(format!(
                        "ascending CM class polynomial for discriminant {d} must be registered"
                    ))
                })?;
            let power = multiplicity(source.degree(), counts.ascending, upper.degree())?;
            Some((&upper.coefficients, power))
        };

        let norm = resultant(
            model,
            &source.coefficients,
            &modular,
            elliptic_modular::degree(level),
        )?;
        let horizontal_factor = horizontal.power(model, counts.horizontal)?;
        let upper_factor = if let Some((polynomial, power)) = upper {
            polynomial.power(model, power)?
        } else {
            Poly::one()
        };
        let remaining = norm
            .clone()
            .exact_divide(model, &horizontal_factor)?
            .exact_divide(model, &upper_factor)?;
        let recovered = remaining.monic_root(model, target_degree, down_power)?;
        let reconstructed = recovered
            .power(model, down_power)?
            .multiply(model, &horizontal_factor)?
            .multiply(model, &upper_factor)?;
        if !norm.same(&reconstructed) {
            return Err(Error::plain(
                "CM class polynomial resultant factorization failed verification",
            ));
        }
        Self::persist(model, discriminant, recovered)
    }

    fn persist(model: &mut Model, discriminant: i128, coefficients: Poly) -> Result<Self> {
        model.allocate(2 + coefficients.0.len() * (LIMBS + 1))?;
        Ok(Self {
            discriminant,
            coefficients,
        })
    }

    pub(crate) fn same(&self, other: &Self) -> bool {
        self.discriminant == other.discriminant && self.coefficients.same(&other.coefficients)
    }

    pub(crate) fn degree(&self) -> usize {
        self.coefficients.0.len() - 1
    }

    pub(crate) fn coefficient_text(&self, index: usize) -> String {
        decimal(self.coefficients.coefficient(index))
    }

    pub(crate) fn modular_coefficient(
        &self,
        model: &mut Model,
        ring: &Ring,
        index: usize,
    ) -> Result<Vec<u64>> {
        Ok(ring.scalar(residue(
            model,
            self.coefficients.coefficient(index),
            ring.prime(),
        )?))
    }

    pub(crate) fn evaluate(
        &self,
        model: &mut Model,
        ring: &Ring,
        input: &[u64],
    ) -> Result<Vec<u64>> {
        let mut value = ring.scalar(0);
        let p = ring.prime();
        for &coefficient in self.coefficients.0.iter().rev() {
            value = ring.multiply(model, &value, input)?;
            value[0] = rv_unram::add_mod(value[0], residue(model, coefficient, p)?, p);
        }
        Ok(value)
    }

    pub(crate) fn unramified_coefficients(
        &self,
        model: &mut Model,
        ring: &Ring,
    ) -> Result<Vec<Vec<u64>>> {
        self.coefficients
            .0
            .iter()
            .map(|&coefficient| Ok(ring.scalar(residue(model, coefficient, ring.modulus())?)))
            .collect()
    }
}

fn multiplicity(source: usize, edges: usize, target: usize) -> Result<usize> {
    let total = source.checked_mul(edges).ok_or_else(limit)?;
    if target == 0 || !total.is_multiple_of(target) || total == 0 {
        return Err(Error::plain(
            "CM cyclic kernel multiplicity is inconsistent with class numbers",
        ));
    }
    Ok(total / target)
}

fn resultant(
    model: &mut Model,
    source: &Poly,
    modular: &elliptic_modular::Polynomial,
    degree: usize,
) -> Result<Poly> {
    let h = source.0.len() - 1;
    model
        .spend(
            h.saturating_mul(h)
                .saturating_mul(h.saturating_mul(degree).saturating_add(1))
                .saturating_mul(LIMBS),
        )
        .map_err(|_| limit())?;
    let zero = Integer::from_i128(0);
    let mut matrix = vec![vec![Poly(vec![zero; degree + 1]); h]; h];
    // The j-th coefficient in Y acts by multiplication in Z[X]/(H).
    for j in 0..=degree {
        let mut coefficients = (0..=degree)
            .map(|i| modular.wide_coefficient(i, j))
            .collect::<Result<Vec<Integer>>>()?;
        for i in (h..coefficients.len()).rev() {
            let value = coefficients[i];
            if !value.is_zero() {
                for k in 0..h {
                    let index = i - h + k;
                    coefficients[index] =
                        coefficients[index].add(value.multiply_wide(source.0[k])?.negate())?;
                }
            }
        }
        coefficients.resize(h, zero);
        // Store multiplication columns as rows. The transposed matrix has
        // the same determinant and monic leading principal minors.
        for (index, column) in matrix.iter_mut().enumerate() {
            for (entry, &value) in column.iter_mut().zip(&coefficients) {
                entry.0[j] = value;
            }
            if index + 1 < h {
                let last = coefficients[h - 1];
                for row in (0..h).rev() {
                    let shifted = if row == 0 {
                        zero
                    } else {
                        coefficients[row - 1]
                    };
                    coefficients[row] = shifted.add(last.multiply_wide(source.0[row])?.negate())?;
                }
            }
        }
    }
    for row in &mut matrix {
        for entry in row {
            entry.trim();
        }
    }

    // M(Y) starts with Y^degree I. Every leading principal minor is
    // monic, so Bareiss elimination needs no pivot swaps and each
    // division is exact monic division in Z[Y].
    let mut predecessor = Poly::one();
    for k in 0..h - 1 {
        let pivot = matrix[k][k].clone();
        for i in k + 1..h {
            for j in k + 1..h {
                let diagonal = pivot.multiply(model, &matrix[i][j])?;
                let cross = matrix[i][k].multiply(model, &matrix[k][j])?;
                matrix[i][j] = diagonal
                    .subtract(&cross)?
                    .exact_divide(model, &predecessor)?;
            }
        }
        predecessor = pivot;
    }
    let determinant = matrix
        .pop()
        .and_then(|mut row| row.pop())
        .ok_or_else(arithmetic)?;
    if !determinant.is_monic() || determinant.0.len() != h * degree + 1 {
        return Err(Error::plain(
            "CM modular norm has an inconsistent degree or leading coefficient",
        ));
    }
    Ok(determinant)
}

#[derive(Clone)]
struct Poly(Vec<Integer>);

impl Poly {
    fn one() -> Self {
        Self(vec![Integer::from_i128(1)])
    }

    fn coefficient(&self, i: usize) -> Integer {
        self.0
            .get(i)
            .copied()
            .unwrap_or_else(|| Integer::from_i128(0))
    }

    fn is_monic(&self) -> bool {
        self.0
            .last()
            .is_some_and(|c| !c.negative && c.magnitude.small() == Some(1))
    }

    fn trim(&mut self) {
        while self.0.last().is_some_and(|c| c.is_zero()) {
            self.0.pop();
        }
    }

    fn same(&self, other: &Self) -> bool {
        self.0.len() == other.0.len()
            && self.0.iter().zip(&other.0).all(|(a, b)| {
                a.magnitude == b.magnitude && (a.is_zero() || a.negative == b.negative)
            })
    }

    fn multiply(&self, model: &mut Model, other: &Self) -> Result<Self> {
        if self.0.is_empty() || other.0.is_empty() {
            return Ok(Self(Vec::new()));
        }
        model
            .spend(
                LIMBS
                    .saturating_mul(self.0.len())
                    .saturating_mul(other.0.len()),
            )
            .map_err(|_| limit())?;
        let mut coefficients = vec![Integer::from_i128(0); self.0.len() + other.0.len() - 1];
        for (i, &a) in self.0.iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for (j, &b) in other.0.iter().enumerate() {
                if !b.is_zero() {
                    coefficients[i + j] = coefficients[i + j].add(a.multiply_wide(b)?)?;
                }
            }
        }
        Ok(Self(coefficients))
    }

    fn subtract(mut self, other: &Self) -> Result<Self> {
        self.0
            .resize(self.0.len().max(other.0.len()), Integer::from_i128(0));
        for (target, &value) in self.0.iter_mut().zip(&other.0) {
            *target = target.add(value.negate())?;
        }
        self.trim();
        Ok(self)
    }

    fn exact_divide(mut self, model: &mut Model, divisor: &Self) -> Result<Self> {
        if !divisor.is_monic() {
            return Err(Error::plain("CM modular norm requires a monic divisor"));
        }
        let degree = divisor.0.len() - 1;
        let terms = self.0.len().saturating_sub(degree);
        model
            .spend(LIMBS.saturating_mul(terms).saturating_mul(divisor.0.len()))
            .map_err(|_| limit())?;
        let mut quotient = vec![Integer::from_i128(0); terms];
        for (offset, target) in quotient.iter_mut().enumerate().rev() {
            let coefficient = self.0[offset + degree];
            *target = coefficient;
            if !coefficient.is_zero() {
                for (i, &value) in divisor.0.iter().take(degree).enumerate() {
                    self.0[offset + i] =
                        self.0[offset + i].add(coefficient.multiply_wide(value)?.negate())?;
                }
            }
        }
        if self.0.iter().take(degree).any(|c| !c.is_zero()) {
            return Err(Error::plain(
                "CM modular norm has a nonzero remainder on monic division",
            ));
        }
        Ok(Self(quotient))
    }

    fn power(&self, model: &mut Model, exponent: usize) -> Result<Self> {
        let mut value = Self::one();
        for _ in 0..exponent {
            value = value.multiply(model, self)?;
        }
        Ok(value)
    }

    fn monic_root(&self, model: &mut Model, degree: usize, power: usize) -> Result<Self> {
        if !self.is_monic() || self.0.len() != degree * power + 1 {
            return Err(Error::plain(
                "CM resultant factor degree does not match the descending class number",
            ));
        }
        if power == 1 {
            return Ok(self.clone());
        }
        let mut root = Self(vec![Integer::from_i128(0); degree + 1]);
        root.0[degree] = Integer::from_i128(1);
        for i in 1..=degree {
            let partial = root.power(model, power)?;
            let k = degree * power - i;
            let difference = self.0[k].add(partial.0[k].negate())?;
            let (magnitude, remainder) = difference
                .magnitude
                .divide_u32(u32::try_from(power).map_err(|_| limit())?);
            if remainder != 0 {
                return Err(Error::plain(
                    "CM resultant factor has no integral monic root",
                ));
            }
            root.0[degree - i] = Integer {
                magnitude,
                negative: difference.negative,
            };
        }
        Ok(root)
    }
}

fn decimal(value: Integer) -> String {
    if value.is_zero() {
        return "0".into();
    }
    let mut magnitude = value.magnitude;
    let mut digits = Vec::new();
    while !magnitude.is_zero() {
        let (rest, digit) = magnitude.divide_u32(1_000_000_000);
        digits.push(digit);
        magnitude = rest;
    }
    let sign = if value.negative { "-" } else { "" };
    let mut output = format!("{sign}{}", digits.pop().expect("nonzero CM coefficient"));
    for digit in digits.into_iter().rev() {
        write!(output, "{digit:09}").expect("writing to a String");
    }
    output
}

fn residue(model: &mut Model, value: Integer, p: u64) -> Result<u64> {
    let (_, remainder) = value.magnitude.divide(model, u128::from(p))?;
    let remainder = remainder as u64;
    Ok(if value.negative && remainder != 0 {
        p - remainder
    } else {
        remainder
    })
}

fn limit() -> Error {
    Error::plain("CM class polynomial norm exceeds the computation limit")
}

fn arithmetic() -> Error {
    Error::plain("CM class polynomial norm exceeds the supported exact arithmetic range")
}
