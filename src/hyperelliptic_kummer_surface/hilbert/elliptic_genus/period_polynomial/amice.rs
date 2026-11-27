use std::cmp::Ordering;
use std::collections::HashMap;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result, Span};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::rv_unram::{add_mod, mul_mod, sub_mod};
use crate::syntax::{Label, ParsedExpr};

use super::super::residual_hecke::power_mod;
use super::mazur_tate::{prime_power, Ordinary};

pub(crate) struct Spec(
    pub(crate) Label,
    pub(crate) usize,
    pub(crate) usize,
    pub(crate) usize,
);

pub(crate) enum Expression {
    Transform(Spec),
    Project(Box<ParsedExpr>),
    GroupCoefficient(usize, Box<ParsedExpr>),
}

pub(crate) enum Report {
    Coefficient(usize),
    Iwasawa(usize),
    Character(usize),
}

impl Report {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            Self::Coefficient(_) => "rvMazurAmiceCoeff",
            Self::Iwasawa(_) => "rvMazurAmiceIwasawa",
            Self::Character(_) => "rvMazurAmiceCharacter",
        }
    }
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) ordinary: Box<str>,
    pub(crate) level: usize,
    pub(crate) index: usize,
    pub(crate) character: usize,
    pub(super) prime: Prime,
    pub(super) precision: usize,
    modulus: u64,
    length: usize,
}

impl Grade {
    pub(crate) fn resolve(model: &mut Model, spec: Spec) -> Result<(Self, Rc<Ordinary>)> {
        let Spec(name, level, index, character) = spec;
        let measure = model.mazur_ordinary(&name)?;
        measure.check_index(index)?;
        let prime = measure.prime;
        let p = prime.get();
        if level < base_level(p) {
            return Err(Error::plain(
                "Amice level must be at least one for odd primes and two for p=2",
            ));
        }
        let tame_order = if p == 2 { 2 } else { p - 1 } as usize;
        if character >= tame_order {
            return Err(Error::plain(format!(
                "Amice tame character index must be between 0 and {}",
                tame_order - 1
            )));
        }
        let cyclotomic_modulus = prime_power(prime, level)?;
        let tame_modulus = if p == 2 { 4 } else { u64::from(p) };
        let length = bounded_level_size(cyclotomic_modulus / tame_modulus)?;
        Ok((
            Self {
                ordinary: name.text.into(),
                level,
                index,
                character,
                prime,
                precision: measure.precision,
                modulus: measure.modulus,
                length,
            },
            measure,
        ))
    }

    pub(super) fn principal_level(&self) -> usize {
        self.level - base_level(self.prime.get())
    }

    pub(crate) fn project(&self) -> Result<Self> {
        if self.principal_level() == 0 {
            return Err(Error::plain("cannot project below the first Amice level"));
        }
        Ok(Self {
            level: self.level - 1,
            length: self.length / self.prime.get() as usize,
            ..self.clone()
        })
    }

    pub(crate) fn label(&self, span: Span) -> Label {
        Label {
            text: self.ordinary.to_string(),
            span,
        }
    }

    pub(crate) fn cardinality(&self) -> usize {
        let modulus = usize::try_from(self.modulus).unwrap_or(usize::MAX);
        modulus.saturating_pow(self.length as u32)
    }
}

pub(crate) struct Amice(Vec<u64>);

impl Amice {
    pub(crate) fn zero(grade: &Grade) -> Self {
        Self(vec![0; grade.length])
    }

    pub(crate) fn transform(model: &mut Model, grade: &Grade, measure: &Ordinary) -> Result<Self> {
        let level = measure.prepare(model, grade.level, grade.index)?;
        let size = bounded_level_size(level.modulus)?;
        model.spend(size + grade.length)?;
        let principal_logs = principal_logs(grade.prime.get(), level.modulus, grade.length)?;
        let tame = Tame::new(model, grade, level.modulus)?;
        let mut result = Self::zero(grade);
        let p = u64::from(grade.prime.get());
        for a in 1..level.modulus {
            if !a.is_multiple_of(p) {
                let (principal, weight) = tame.decompose(a);
                let j = *principal_logs.get(&principal).ok_or_else(|| {
                    Error::plain("Amice principal unit is outside its cyclic group")
                })?;
                let value = measure.measure_prepared(model, &level, a)?;
                let term = mul_mod(weight, value, grade.modulus);
                result.0[j] = add_mod(result.0[j], term, grade.modulus);
            }
        }
        Ok(result)
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + self.0.len())?;
        Ok(Value::MazurAmice(Rc::new(self)))
    }

    pub(crate) fn coefficients(&self) -> &[u64] {
        &self.0
    }

    pub(crate) fn group_coefficient(&self, j: usize) -> u64 {
        self.0[j % self.0.len()]
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(self.0.len())?;
        Ok(self.0.cmp(&rhs.0))
    }

    pub(crate) fn add(&self, model: &mut Model, grade: &Grade, rhs: &Self) -> Result<Self> {
        model.spend(grade.length)?;
        Ok(Self(
            self.0
                .iter()
                .zip(&rhs.0)
                .map(|(&left, &right)| add_mod(left, right, grade.modulus))
                .collect(),
        ))
    }

    pub(crate) fn negate(&self, model: &mut Model, grade: &Grade) -> Result<Self> {
        model.spend(grade.length)?;
        Ok(Self(
            self.0
                .iter()
                .map(|&value| sub_mod(0, value, grade.modulus))
                .collect(),
        ))
    }

    pub(crate) fn project(&self, model: &mut Model, target: &Grade) -> Result<Self> {
        model.spend(self.0.len())?;
        let mut result = Self::zero(target);
        for (j, &value) in self.0.iter().enumerate() {
            let slot = &mut result.0[j % target.length];
            *slot = add_mod(*slot, value, target.modulus);
        }
        Ok(result)
    }

    pub(crate) fn report(
        &self,
        model: &mut Model,
        grade: &Grade,
        report: Report,
    ) -> Result<String> {
        let name = report.name();
        match report {
            Report::Coefficient(k) => {
                let digits = coefficient_precision(grade, k);
                if digits == 0 {
                    return Ok(format!("{name}[{k}] = precision=0, unknown"));
                }
                let modulus = prime_power(grade.prime, digits)?;
                let value = self.coefficient(model, k, modulus)?;
                Ok(format!(
                    "{name}[{k}] = precision={digits}, modulus={modulus}, residue={value}"
                ))
            }
            Report::Iwasawa(limit) => {
                let modulus = u64::from(grade.prime.get());
                // The first L coefficients are all known modulo p, and higher
                // coefficients can vary modulo p between compatible lifts.
                let checked = limit.min(grade.length - 1);
                for k in 0..=checked {
                    if self.coefficient(model, k, modulus)? != 0 {
                        return Ok(format!("{name}[{limit}] = mu=0, lambda={k}"));
                    }
                }
                Ok(format!(
                    "{name}[{limit}] = mu=unknown, lambda=unknown, first_unit_degree_at_least={}",
                    checked + 1
                ))
            }
            Report::Character(depth) => {
                let (order, polynomial) = self.character(model, grade, depth)?;
                Ok(format!(
                    "{name}[{depth}] = order={order}, modulus={}, polynomial={polynomial:?}",
                    grade.modulus
                ))
            }
        }
    }

    fn coefficient(&self, model: &mut Model, k: usize, modulus: u64) -> Result<u64> {
        if k >= self.0.len() {
            return Ok(0);
        }
        model.spend(self.0.len().saturating_mul(k + 1))?;
        // Begin Pascal's recurrence at row k, where the requested binomial is one.
        // The preceding entries are row k and can be constructed at the same cost.
        let mut row = vec![0; k + 1];
        row[0] = 1;
        for end in 1..=k {
            advance_binomial(&mut row[..=end], modulus);
        }
        let mut value = 0;
        for &group_coefficient in &self.0[k..] {
            value = add_mod(value, mul_mod(group_coefficient, row[k], modulus), modulus);
            advance_binomial(&mut row, modulus);
        }
        Ok(value)
    }

    fn character(
        &self,
        model: &mut Model,
        grade: &Grade,
        depth: usize,
    ) -> Result<(usize, Vec<u64>)> {
        if depth > grade.principal_level() {
            return Err(Error::plain(
                "Amice character conductor exceeds the available principal-unit level",
            ));
        }
        let order = prime_power(grade.prime, depth)? as usize;
        model.spend(self.0.len() + order)?;
        let mut values = vec![0; order];
        for (j, &value) in self.0.iter().enumerate() {
            let slot = &mut values[j % order];
            *slot = add_mod(*slot, value, grade.modulus);
        }
        if depth > 0 {
            let block = order / grade.prime.get() as usize;
            let degree = order - block;
            for j in 0..degree {
                values[j] = sub_mod(values[j], values[degree + j % block], grade.modulus);
            }
            values.truncate(degree);
        }
        Ok((order, values))
    }
}

fn advance_binomial(row: &mut [u64], modulus: u64) {
    for j in (1..row.len()).rev() {
        row[j] = add_mod(row[j], row[j - 1], modulus);
    }
}

fn coefficient_precision(grade: &Grade, mut k: usize) -> usize {
    if k == 0 {
        return grade.precision;
    }
    let p = grade.prime.get() as usize;
    let mut log = 0;
    while k >= p {
        k /= p;
        log += 1;
    }
    grade
        .precision
        .min(grade.principal_level().saturating_sub(log))
}

fn base_level(p: u32) -> usize {
    if p == 2 {
        2
    } else {
        1
    }
}

fn bounded_level_size(size: u64) -> Result<usize> {
    usize::try_from(size)
        .ok()
        .filter(|&size| size <= MAX_PIVOT_CELLS)
        .ok_or_else(|| Error::plain("Amice cyclotomic level exceeds the finite table limit"))
}

fn principal_logs(p: u32, modulus: u64, length: usize) -> Result<HashMap<u64, usize>> {
    let generator = if p == 2 { 5 } else { 1 + u64::from(p) };
    let mut table = HashMap::with_capacity(length);
    let mut power = 1;
    for j in 0..length {
        if table.insert(power, j).is_some() {
            return Err(Error::plain(
                "Amice principal generator has incorrect order",
            ));
        }
        power = mul_mod(power, generator, modulus);
    }
    if power != 1 {
        return Err(Error::plain(
            "Amice principal generator has incorrect order",
        ));
    }
    Ok(table)
}

struct Tame {
    prime: u64,
    level_modulus: u64,
    weights: Vec<(u64, u64)>,
}

impl Tame {
    fn new(model: &mut Model, grade: &Grade, level_modulus: u64) -> Result<Self> {
        let p = u64::from(grade.prime.get());
        let weights = if p == 2 {
            let sign = if grade.character == 0 {
                1
            } else {
                grade.modulus - 1
            };
            vec![(1, 1), (level_modulus - 1, sign)]
        } else {
            let mut weights = vec![(0, 0); p as usize];
            for a in 1..p {
                let omega = power_mod(model, a, u128::from(level_modulus / p), level_modulus)?;
                let inverse = power_mod(model, omega, u128::from(p - 2), level_modulus)?;
                let lift = power_mod(model, a, u128::from(grade.modulus / p), grade.modulus)?;
                let weight = power_mod(model, lift, grade.character as u128, grade.modulus)?;
                weights[a as usize] = (inverse, weight);
            }
            weights
        };
        Ok(Self {
            prime: p,
            level_modulus,
            weights,
        })
    }

    fn decompose(&self, a: u64) -> (u64, u64) {
        let index = if self.prime == 2 {
            (a % 4) / 2
        } else {
            a % self.prime
        };
        let (inverse, weight) = self.weights[index as usize];
        (mul_mod(a, inverse, self.level_modulus), weight)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn checked<T>(result: Result<T>) -> T {
        result.unwrap_or_else(|error| panic!("{error}"))
    }

    #[test]
    fn two_adic_sign_and_principal_components() {
        let mut model = Model::new();
        for level in 2..=6 {
            let length = 1 << (level - 2);
            let level_modulus = 1 << level;
            let mut grade = Grade {
                ordinary: "O".into(),
                level,
                index: 0,
                character: 0,
                prime: Prime::validated(2),
                precision: 4,
                modulus: 16,
                length,
            };
            let trivial = checked(Tame::new(&mut model, &grade, level_modulus));
            grade.character = 1;
            let sign = checked(Tame::new(&mut model, &grade, level_modulus));
            let logs = checked(principal_logs(2, level_modulus, length));
            for a in (1..level_modulus).step_by(2) {
                let (principal, weight) = sign.decompose(a);
                assert!(logs.contains_key(&principal));
                assert_eq!(principal % 4, 1);
                assert_eq!(trivial.decompose(a), (principal, 1));
                assert_eq!(sign.decompose(level_modulus - a), (principal, 16 - weight));
                assert_eq!(weight, if a % 4 == 1 { 1 } else { 15 });
            }
        }
    }

    #[test]
    fn two_adic_projection_and_primitive_evaluation() {
        let mut model = Model::new();
        let grade = Grade {
            ordinary: "O".into(),
            level: 4,
            index: 0,
            character: 1,
            prime: Prime::validated(2),
            precision: 4,
            modulus: 16,
            length: 4,
        };
        let a = Amice(vec![3, 7, 4, 11]);
        assert_eq!(
            checked(a.character(&mut model, &grade, 2)),
            (4, vec![15, 12])
        );
        let low = checked(grade.project());
        let projected = checked(a.project(&mut model, &low));
        assert_eq!(projected.0, [7, 2]);
        assert_eq!(
            checked(projected.character(&mut model, &low, 1)),
            (2, vec![5])
        );
        assert_eq!(checked(a.character(&mut model, &grade, 1)), (2, vec![5]));
        assert_eq!(coefficient_precision(&grade, 0), 4);
        assert_eq!(coefficient_precision(&grade, 1), 2);
        assert_eq!(coefficient_precision(&grade, 2), 1);
        assert_eq!(coefficient_precision(&grade, 3), 1);
        assert_eq!(coefficient_precision(&grade, 4), 0);
    }

    fn disc_fixture(prime: u32, level: usize, precision: usize) -> Grade {
        let prime = Prime::validated(prime);
        Grade {
            ordinary: "O".into(),
            level,
            index: 0,
            character: 0,
            prime,
            precision,
            modulus: checked(prime_power(prime, precision)),
            length: checked(prime_power(prime, level - base_level(prime.get()))) as usize,
        }
    }

    #[test]
    fn weight_disc_counts_multiple_and_nonrational_reductions() {
        use super::super::disc;

        let mut model = Model::new();
        let source = disc_fixture(5, 3, 3);
        let grade = checked(disc::Grade::new(source.clone(), disc::Chart(0, 1)));
        for (row, expected, expected_error) in [
            (
                [36, 113, 1], // (T-5)^2 in the group variable 1+T
                "reduction=[1, 3, 1], residual_roots=[(1, 2)], roots_in_nonrational_residues=0",
                "has multiplicity 2",
            ),
            (
                [51, 123, 1], // T^2 + 2*5^2
                "reduction=[2, 0, 1], residual_roots=[], roots_in_nonrational_residues=2",
                "has no root",
            ),
        ] {
            let mut amplitude = Amice::zero(&source);
            amplitude.0[..row.len()].copy_from_slice(&row);
            let disc = checked(disc::Disc::new(&mut model, &grade, Rc::new(amplitude)));
            let report = checked(disc.report(&mut model, &grade, disc::Report::Zeros));
            assert!(report.starts_with("rvMazurDiscZeros = precision=3, mu=2, lambda=2"));
            assert!(report.ends_with(expected), "{report}");
            let Err(error) = disc.root(&mut model, &grade, 1) else {
                panic!("residue 1 is not simple");
            };
            assert!(error.message.contains(expected_error));
        }
    }

    #[test]
    fn finite_representatives_agree_with_the_dyadic_and_odd_chart_precision() {
        use super::super::disc;

        let mut model = Model::new();
        for (p, level, center, depth, precision) in [
            (2, 5, 2, 2, 4),
            (2, 5, 0, 1, 4),
            (3, 3, -9, 2, 4),
            (5, 2, 0, 3, 4),
        ] {
            let source = disc_fixture(p, level, 6);
            let grade = checked(disc::Grade::new(source.clone(), disc::Chart(center, depth)));
            let mut relation = vec![0; source.length + 1];
            relation[0] = source.modulus - 1;
            relation[source.length] = 1;
            // Deliberately pass an unreduced representative Z^L - 1.
            let disc = checked(disc::Disc::new(
                &mut model,
                &grade,
                Rc::new(Amice(relation)),
            ));
            assert_eq!(
                checked(disc.report(&mut model, &grade, disc::Report::Zeros)),
                format!("rvMazurDiscZeros = precision={precision}, unresolved")
            );
        }
    }
}
