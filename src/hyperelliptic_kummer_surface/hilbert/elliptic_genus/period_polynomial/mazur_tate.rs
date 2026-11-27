use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result, Span};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear;
use crate::rv_unram::{add_mod, mul_mod, sub_mod};
use crate::syntax::{Label, ParsedExpr};

use super::super::eigenpacket::{Element, Packet};
use super::super::hecke::Operator;
use super::super::residual_hecke::{invert_unit, lift_simple_root, power_mod, Specialization};
use super::super::{scale_by_power, Integer, Magnitude, INTEGER_WORDS};
use super::{binomial_table, eigenperiod::Eigenperiod, modular_symbol::scaled_moment};

pub(crate) struct Spec(
    pub(crate) Label,
    pub(crate) usize,
    pub(crate) usize,
    pub(crate) usize,
);

pub(crate) enum Expression {
    Moment(Spec, i128),
    Theta(Spec),
    Coefficient(i128, Box<ParsedExpr>),
    Change(Direction, Box<ParsedExpr>),
    Scale(Box<[ParsedExpr; 2]>),
    Measure(Label, usize, usize, i128),
}

pub(crate) enum Statement {
    Ordinary(Label, Label, Embedding, usize),
    Data(Label),
}

pub(crate) enum Embedding {
    Rational(usize),
    Residual(Label),
}

#[derive(Clone, Copy)]
pub(crate) enum Direction {
    Trace,
    Transfer,
}

impl Direction {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Trace => "rvMazurTateTrace",
            Self::Transfer => "rvMazurTateTransfer",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub(crate) struct Grade {
    pub(crate) period: Box<str>,
    pub(crate) prime: Prime,
    pub(crate) level: usize,
    pub(crate) index: usize,
    pub(crate) modulus: u64,
}

impl Grade {
    pub(crate) fn resolve(model: &mut Model, spec: Spec) -> Result<(Self, Rc<Eigenperiod>)> {
        let Spec(name, input_prime, level, index) = spec;
        let period = model.hecke_eigenperiod(&name)?;
        check_index(index, &period)?;
        let prime = checked_prime(model, input_prime)?;
        let modulus = prime_power(prime, level)?;
        Ok((
            Self {
                period: name.text.into(),
                prime,
                level,
                index,
                modulus,
            },
            period,
        ))
    }

    pub(crate) fn label(&self, span: Span) -> Label {
        Label {
            text: self.period.to_string(),
            span,
        }
    }

    pub(crate) fn change(&self, direction: Direction) -> Result<Self> {
        let level = match direction {
            Direction::Trace => self
                .level
                .checked_sub(1)
                .ok_or_else(|| Error::plain("Mazur–Tate trace requires level at least one"))?,
            Direction::Transfer => self.level.checked_add(1).ok_or_else(level_limit)?,
        };
        let modulus = prime_power(self.prime, level)?;
        Ok(Self {
            level,
            modulus,
            ..self.clone()
        })
    }

    fn unit(&self, residue: u64) -> bool {
        self.modulus == 1 || !residue.is_multiple_of(u64::from(self.prime.get()))
    }

    pub(crate) fn moment(
        &self,
        model: &mut Model,
        period: &Eigenperiod,
        representative: i128,
    ) -> Result<Value> {
        let binomial = binomial_table(model, period.grade.degree())?;
        let a = representative.rem_euclid(i128::from(self.modulus)) as u64;
        scaled_moment(model, period, self.modulus, self.index, a, &binomial)?.persist(model)
    }
}

type Support = BTreeMap<u64, Element>;

pub(crate) struct Distribution {
    degree: usize,
    coefficients: Support,
}

impl Distribution {
    pub(crate) fn zero_for(period: &Eigenperiod) -> Self {
        Self::zero(period.packet.degree())
    }

    fn zero(degree: usize) -> Self {
        Self {
            degree,
            coefficients: Support::new(),
        }
    }

    pub(crate) fn theta(model: &mut Model, grade: &Grade, period: &Eigenperiod) -> Result<Self> {
        let mut result = Self::zero(period.packet.degree());
        charge_residues(model, grade.modulus)?;
        let binomial = binomial_table(model, period.grade.degree())?;
        for a in 0..grade.modulus {
            if grade.unit(a) {
                let value = scaled_moment(model, period, grade.modulus, grade.index, a, &binomial)?;
                result.insert(model, a, value)?;
            }
        }
        Ok(result)
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + self.coefficients.len() * self.term_cells())?;
        Ok(Value::MazurTate(Rc::new(self)))
    }

    pub(crate) fn terms(&self) -> impl Iterator<Item = (u64, &Element)> + '_ {
        self.coefficients.iter().map(|(&a, value)| (a, value))
    }

    pub(crate) fn coefficient(&self, model: &mut Model, grade: &Grade, a: i128) -> Result<Value> {
        let a = a.rem_euclid(i128::from(grade.modulus)) as u64;
        self.coefficients
            .get(&a)
            .cloned()
            .unwrap_or_else(|| Element::zero(self.degree))
            .persist(model)
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        model.spend(1 + self.coefficients.len() * self.term_cells())?;
        let mut result = Self {
            degree: self.degree,
            coefficients: self.coefficients.clone(),
        };
        for (a, value) in rhs.terms() {
            result.insert(model, a, value.clone())?;
        }
        Ok(result)
    }

    pub(crate) fn negate(&self, model: &mut Model) -> Result<Self> {
        let coefficients = self
            .terms()
            .map(|(a, value)| Ok((a, value.negate(model)?)))
            .collect::<Result<_>>()?;
        Ok(Self {
            degree: self.degree,
            coefficients,
        })
    }

    pub(crate) fn scale(
        &self,
        model: &mut Model,
        packet: &Packet,
        scalar: &Element,
    ) -> Result<Self> {
        let mut result = Self::zero(self.degree);
        if scalar.is_zero() {
            return Ok(result);
        }
        for (a, value) in self.terms() {
            let product = packet.multiply_elements(model, value, scalar)?;
            result.insert(model, a, product)?;
        }
        Ok(result)
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        for ((a, x), (b, y)) in self.terms().zip(rhs.terms()) {
            model.spend(1)?;
            let order = a.cmp(&b);
            if !order.is_eq() {
                return Ok(order);
            }
            let order = x.compare(model, y)?;
            if !order.is_eq() {
                return Ok(order);
            }
        }
        Ok(self.coefficients.len().cmp(&rhs.coefficients.len()))
    }

    pub(crate) fn change(
        &self,
        model: &mut Model,
        source: &Grade,
        target: &Grade,
        direction: Direction,
    ) -> Result<Self> {
        let mut result = Self::zero(self.degree);
        match direction {
            Direction::Trace => {
                for (a, value) in self.terms() {
                    result.insert(model, a % target.modulus, value.clone())?;
                }
            }
            Direction::Transfer => {
                let p = u64::from(source.prime.get());
                for (a, value) in self.terms() {
                    charge_residues(model, p)?;
                    for b in 0..p {
                        let lift = a + b * source.modulus;
                        if target.unit(lift) {
                            result.insert(model, lift, value.clone())?;
                        }
                    }
                }
            }
        }
        Ok(result)
    }

    fn term_cells(&self) -> usize {
        (self.degree + 1) * INTEGER_WORDS + 6
    }

    fn insert(&mut self, model: &mut Model, a: u64, value: Element) -> Result<()> {
        model.spend(1 + self.coefficients.len().max(1).ilog2() as usize)?;
        if value.is_zero() {
            return Ok(());
        }
        if let Some(old) = self.coefficients.get_mut(&a) {
            let sum = old.add(model, &value)?;
            if sum.is_zero() {
                self.coefficients.remove(&a);
            } else {
                *old = sum;
            }
        } else {
            if self.coefficients.len() >= MAX_PIVOT_CELLS / self.term_cells() {
                return Err(Error::plain(
                    "Mazur–Tate group-ring support exceeds the storage limit",
                ));
            }
            self.coefficients.insert(a, value);
        }
        Ok(())
    }
}

pub(crate) struct Ordinary {
    period: Rc<Eigenperiod>,
    embedding: Rc<Specialization>,
    precision: usize,
    pub(crate) modulus: u64,
    a_p: u64,
    alpha: u64,
    inverse_alpha: u64,
    lower_factor: u64,
    denominator_exponent: u32,
    denominator_scale: Integer,
}

impl Ordinary {
    pub(crate) fn new(
        model: &mut Model,
        name: &Label,
        embedding: Embedding,
        precision: usize,
    ) -> Result<Self> {
        let period = model.hecke_eigenperiod(name)?;
        if precision == 0 {
            return Err(Error::plain(
                "ordinary Mazur–Tate precision must be positive",
            ));
        }
        let embedding = match embedding {
            Embedding::Rational(p) => {
                let prime = checked_prime(model, p)?;
                Rc::new(Specialization::rational(
                    model,
                    period.source.clone(),
                    period.packet.clone(),
                    prime,
                )?)
            }
            Embedding::Residual(name) => {
                let embedding = model.hecke_residual(&name)?;
                if embedding.source.text != period.source.text {
                    return Err(Error::new(
                        name.span,
                        "ordinary Mazur–Tate embedding must use the eigenperiod's named packet",
                    ));
                }
                embedding
            }
        };
        let prime = embedding.prime;
        let p = u64::from(prime.get());
        let modulus = prime_power(prime, precision)?;
        let t_p = Operator::basis(u128::from(p))?;
        let eigenvalue = period.packet.image_element(model, &t_p)?;
        let a_p = embedding.reduce_to_modulus(model, &eigenvalue, modulus)?;
        let residue = (a_p % p) as u32;
        if residue == 0 {
            return Err(Error::plain(
                "ordinary Mazur–Tate measure requires a p-adic unit a_p",
            ));
        }
        let constant = power_mod(model, p % modulus, period.packet.weight - 1, modulus)?;
        let inverse_derivative = prime_linear::invert(model, prime, residue)?;
        let alpha = lift_simple_root(
            model,
            prime,
            residue,
            inverse_derivative,
            modulus,
            |_, alpha, modulus| {
                let square = mul_mod(alpha, alpha, modulus);
                let linear = mul_mod(a_p % modulus, alpha, modulus);
                Ok(add_mod(
                    sub_mod(square, linear, modulus),
                    constant % modulus,
                    modulus,
                ))
            },
        )?;
        let inverse_alpha = invert_unit(model, prime.get(), alpha, modulus)?;
        let lower_factor = power_mod(model, p % modulus, period.grade.degree() as u128, modulus)?;
        let lower_factor = mul_mod(lower_factor, inverse_alpha, modulus);

        // After scaling the full period by p^e, integral matrices and
        // the integral Horner formula keep all S_n(a) integral.
        let mut denominator_exponent = 0;
        for coefficient in &period.coefficients {
            denominator_exponent =
                denominator_exponent.max(embedding.integrality_deficit(model, coefficient)?);
        }
        let denominator_scale = scale_by_power(
            model,
            Integer::from_i128(1),
            u128::from(p),
            u128::from(denominator_exponent),
        )?;
        Ok(Self {
            period,
            embedding,
            precision,
            modulus,
            a_p,
            alpha,
            inverse_alpha,
            lower_factor,
            denominator_exponent,
            denominator_scale,
        })
    }

    pub(crate) fn data(&self, name: &str) -> String {
        let p = self.embedding.prime.get();
        format!(
            "rvMazurTateOrdinaryData[{name}] = p={p}, precision={}, modulus={}, a_p={}, alpha={}, period_scale={p}^{}",
            self.precision, self.modulus, self.a_p, self.alpha, self.denominator_exponent
        )
    }

    pub(crate) fn measure(
        &self,
        model: &mut Model,
        level: usize,
        index: usize,
        representative: i128,
    ) -> Result<u64> {
        if level == 0 {
            return Err(Error::plain(
                "ordinary Mazur–Tate measure requires level at least one",
            ));
        }
        check_index(index, &self.period)?;
        let p = self.embedding.prime;
        let level_modulus = prime_power(p, level)?;
        let lower_modulus = level_modulus / u64::from(p.get());
        let binomial = binomial_table(model, self.period.grade.degree())?;
        let a = representative.rem_euclid(i128::from(level_modulus)) as u64;
        let a_lower = a % lower_modulus;
        let moment = self.reduced_moment(model, level_modulus, index, a, &binomial)?;
        let lower = self.reduced_moment(model, lower_modulus, index, a_lower, &binomial)?;
        let correction = mul_mod(self.lower_factor, lower, self.modulus);
        let scale = power_mod(model, self.inverse_alpha, level as u128, self.modulus)?;
        Ok(mul_mod(
            scale,
            sub_mod(moment, correction, self.modulus),
            self.modulus,
        ))
    }

    fn reduced_moment(
        &self,
        model: &mut Model,
        modulus: u64,
        index: usize,
        a: u64,
        binomial: &[Vec<i128>],
    ) -> Result<u64> {
        let value = scaled_moment(model, &self.period, modulus, index, a, binomial)?;
        let scaled = value.scale(model, self.denominator_scale, Magnitude::from_small(1))?;
        self.embedding
            .reduce_to_modulus(model, &scaled, self.modulus)
    }
}

fn check_index(index: usize, period: &Eigenperiod) -> Result<()> {
    if index > period.grade.degree() {
        return Err(Error::plain(
            "Mazur–Tate moment index must be between 0 and k−2",
        ));
    }
    Ok(())
}

fn checked_prime(model: &mut Model, p: usize) -> Result<Prime> {
    if p == 2 {
        Ok(Prime::validated(2))
    } else if p < 2 || p.is_multiple_of(2) {
        Err(Error::plain(format!("Mazur–Tate prime {p} is not prime")))
    } else {
        Prime::new(model, p)
    }
}

fn prime_power(p: Prime, level: usize) -> Result<u64> {
    let exponent = u32::try_from(level).map_err(|_| level_limit())?;
    u64::from(p.get())
        .checked_pow(exponent)
        .ok_or_else(level_limit)
}

fn charge_residues(model: &mut Model, count: u64) -> Result<()> {
    let work = usize::try_from(count).map_err(|_| level_limit())?;
    model.spend(work)
}

fn level_limit() -> Error {
    Error::plain("Mazur–Tate prime-power level or precision exceeds the unsigned 64-bit range")
}
