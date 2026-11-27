use std::collections::BTreeMap;
use std::rc::Rc;

use super::{even_index, even_mask, Action, Algebra, Element as Spin, EVEN_RANK};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::mukai::abelian::og6::cusp::nilpotent::Limit;
use crate::hyperelliptic_kummer_surface::rational_sparse::{
    Coefficients, Domain, Integer, Magnitude, WORK_LIMBS,
};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

mod smith;

pub(crate) enum Expression {
    Basis(Label, u64),
    Embed(Box<ParsedExpr>),
    ToSpin(Box<ParsedExpr>),
    Cup(Box<[ParsedExpr; 2]>),
    Scale(i128, Box<ParsedExpr>),
    Apply(Label, Operation, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Power(usize),
    Divided(usize),
    Exp(i128),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Power(_) => "rvKSExteriorPower",
            Self::Divided(_) => "rvKSExteriorDivided",
            Self::Exp(_) => "rvKSExteriorExp",
        }
    }
}

pub(crate) enum Command {
    Data(Label, usize),
    Smith(Label, usize, usize),
    Coefficient(u64, Box<ParsedExpr>),
}

pub(crate) enum Exterior {}

impl Domain for Exterior {
    const NAME: &'static str = "Kuga–Satake exterior";
}

pub(crate) type Form = Coefficients<Exterior, u64>;
type Terms = BTreeMap<u64, Integer>;

impl Form {
    pub(crate) fn blade(mask: u64) -> Self {
        Self::singleton(mask, 1)
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        self.charge_storage(model)?;
        Ok(Value::KSExterior(Rc::new(self)))
    }

    fn integral(model: &mut Model, terms: &Terms) -> Result<Self> {
        Self::normalized(
            model,
            terms.iter().map(|(&mask, value)| (mask, value)),
            Magnitude::from_small(1),
        )
    }

    pub(crate) fn embed(model: &mut Model, spin: &Spin) -> Result<Self> {
        let terms = spin
            .terms()
            .iter()
            .map(|&(mask, coefficient)| Ok((1 << even_index(mask), coefficient.resize()?)))
            .collect::<Result<Terms>>()?;
        Self::integral(model, &terms)
    }

    pub(crate) fn to_spin(&self, model: &mut Model) -> Result<Spin> {
        let terms = self
            .terms()
            .iter()
            .map(|&(bit, coefficient)| {
                debug_assert_eq!(bit.count_ones(), 1);
                Ok((
                    even_mask(bit.trailing_zeros() as usize),
                    coefficient.resize()?,
                ))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;
        Spin::integral(model, &terms)
    }

    pub(crate) fn cup(&self, model: &mut Model, right: &Self) -> Result<Self> {
        let (Some((a, _)), Some((b, _))) = (self.terms().first(), right.terms().first()) else {
            return Ok(Self::zero());
        };
        if (a.count_ones() + b.count_ones()) as usize > EVEN_RANK {
            return Ok(Self::zero());
        }
        let mut terms = Terms::new();
        for &(a, left_coefficient) in self.terms() {
            model.spend(right.terms().len() * (WORK_LIMBS + EVEN_RANK))?;
            let left_coefficient: Integer = left_coefficient.resize()?;
            for &(b, right_coefficient) in right.terms() {
                if a & b == 0 {
                    let mut coefficient =
                        left_coefficient.multiply_wide(right_coefficient.resize()?)?;
                    if wedge_negative(a, b) {
                        coefficient = coefficient.negate();
                    }
                    Self::accumulate(&mut terms, a | b, coefficient)?;
                }
            }
        }
        Self::integral(model, &terms)
    }

    fn widen(&self) -> Result<Terms> {
        self.terms()
            .iter()
            .map(|&(mask, coefficient)| Ok((mask, coefficient.resize()?)))
            .collect()
    }
}

impl Algebra {
    pub(crate) fn exterior_apply(
        &self,
        model: &mut Model,
        limit: &Limit,
        operation: Operation,
        input: &Form,
    ) -> Result<Form> {
        let Some((mask, _)) = input.terms().first() else {
            return Ok(Form::zero());
        };
        let degree = mask.count_ones() as usize;
        let bound = degree.min(EVEN_RANK - degree);
        match operation {
            Operation::Power(j) | Operation::Divided(j) if j > bound => Ok(Form::zero()),
            Operation::Power(0) | Operation::Divided(0) | Operation::Exp(0) => {
                input.signed_map(model, |_| false)
            }
            Operation::Exp(_) if bound == 0 => input.signed_map(model, |_| false),
            Operation::Power(j) => {
                let action = self.action(model, limit)?;
                let mut terms = input.widen()?;
                for _ in 0..j {
                    terms = action.derivative(model, &terms)?;
                    if terms.is_empty() {
                        return Ok(Form::zero());
                    }
                }
                Form::integral(model, &terms)
            }
            Operation::Divided(j) => {
                self.action(model, limit)?
                    .exterior_divided(model, j, input, degree > EVEN_RANK / 2)
            }
            Operation::Exp(m) => {
                self.action(model, limit)?
                    .exterior_exp(model, m, input, degree > EVEN_RANK / 2)
            }
        }
    }

    pub(crate) fn exterior_report(
        &self,
        model: &mut Model,
        label: &str,
        limit: &Limit,
        degree: usize,
    ) -> Result<String> {
        if degree > EVEN_RANK {
            return Err(Error::plain(
                "Kuga–Satake cohomological degree must be between 0 and 64",
            ));
        }
        let rank = self.action(model, limit)?.smith(model)?.len();
        let character = character(model, rank, degree)?;
        let dimension: u128 = character.iter().sum();
        let graded = character
            .iter()
            .enumerate()
            .filter(|(_, count)| **count != 0)
            .map(|(i, &count)| (degree as isize + i as isize - rank as isize, count))
            .collect::<Vec<_>>();
        // A highest weight j contributes once in weights j,j-2,...,-j.
        let jordan = (0..=rank)
            .filter_map(|j| {
                let count = character[rank + j] - character.get(rank + j + 2).unwrap_or(&0);
                (count != 0).then_some((j + 1, count))
            })
            .collect::<Vec<_>>();
        let index = jordan
            .last()
            .expect("an exterior degree from 0 through 64 is nonzero")
            .0;
        Ok(format!(
            "rvKSExteriorData[{label},{degree}] = dimension={dimension}, nilpotency_index={index}, graded={graded:?}, jordan={jordan:?}"
        ))
    }
}

impl Action {
    fn derivative(&self, model: &mut Model, input: &Terms) -> Result<Terms> {
        let mut terms = Terms::new();
        for (&mask, &coefficient) in input {
            for source in bits(mask) {
                let column = &self.columns[source];
                model.spend((1 + column.len()) * WORK_LIMBS)?;
                let bit = 1 << source;
                let reduced = mask ^ bit;
                let position = (mask & (bit - 1)).count_ones();
                for &(target, entry) in column {
                    let bit = 1 << target;
                    if reduced & bit == 0 {
                        let swap_count = position + (reduced & (bit - 1)).count_ones();
                        let mut value = coefficient.multiply_wide(entry)?;
                        if !swap_count.is_multiple_of(2) {
                            value = value.negate();
                        }
                        Form::accumulate(&mut terms, reduced | bit, value)?;
                    }
                }
            }
        }
        Ok(terms)
    }

    fn exterior_columns(
        &self,
        model: &mut Model,
        m: i128,
        dual: bool,
    ) -> Result<super::RegularColumns> {
        // The complementary wedge pairing identifies Λ^k(A) with
        // Λ^(64−k)(A^−transpose) for det(A)=1. Here A=1+mS and S²=0.
        let mut linear = (0..EVEN_RANK).map(|_| Vec::new()).collect::<Vec<_>>();
        for (source, column) in self.columns.iter().enumerate() {
            model.spend((1 + column.len()) * WORK_LIMBS)?;
            for &(target, entry) in column {
                let scaled = entry.multiply(m)?;
                let (i, j, value) = if dual {
                    (target, source, scaled.negate())
                } else {
                    (source, target, scaled)
                };
                linear[i].push((j, value));
            }
        }
        Ok(linear)
    }

    fn exterior_exp(&self, model: &mut Model, m: i128, input: &Form, dual: bool) -> Result<Form> {
        let mut linear = self.exterior_columns(model, m, dual)?;
        let one = Integer::from_i128(1);
        for (i, column) in linear.iter_mut().enumerate() {
            if let Some((_, entry)) = column.iter_mut().find(|(j, _)| *j == i) {
                *entry = entry.add(one)?;
            } else {
                column.push((i, one));
            }
            column.retain(|(_, coefficient)| !coefficient.is_zero());
        }
        exterior_transform(model, input, dual, |model, mask, coefficient| {
            let mut partial = Terms::from([(0, coefficient)]);
            for source in bits(mask) {
                let mut next = Terms::new();
                wedge_column(model, &mut next, &partial, &linear[source])?;
                partial = next;
            }
            Ok(partial)
        })
    }

    fn exterior_divided(
        &self,
        model: &mut Model,
        order: usize,
        input: &Form,
        dual: bool,
    ) -> Result<Form> {
        let linear = self.exterior_columns(model, 1, dual)?;
        exterior_transform(model, input, dual, |model, mask, coefficient| {
            let mut possible = bits(mask).filter(|&i| !linear[i].is_empty()).count();
            if possible < order {
                return Ok(Terms::new());
            }
            let mut powers = vec![Terms::new(); order + 1];
            powers[0].insert(0, coefficient);
            for (step, source) in bits(mask).enumerate() {
                let mut next = vec![Terms::new(); order + 1];
                possible -= usize::from(!linear[source].is_empty());
                let minimum = order.saturating_sub(possible);
                let identity = [(source, Integer::from_i128(1))];
                for (j, partial) in powers.iter().enumerate().take(step.min(order) + 1) {
                    if j >= minimum {
                        wedge_column(model, &mut next[j], partial, &identity)?;
                    }
                    if j < order && j + 1 >= minimum {
                        wedge_column(model, &mut next[j + 1], partial, &linear[source])?;
                    }
                }
                let entries = next.iter().map(Terms::len).sum::<usize>();
                if entries * (WORK_LIMBS + 3) > crate::complex::MAX_PIVOT_CELLS {
                    return Err(Form::work_limit());
                }
                powers = next;
            }
            Ok(powers
                .pop()
                .expect("the requested divided power is present"))
        })
    }
}

fn exterior_transform(
    model: &mut Model,
    input: &Form,
    dual: bool,
    blade: impl Fn(&mut Model, u64, Integer) -> Result<Terms>,
) -> Result<Form> {
    let mut output = Terms::new();
    for &(mask, coefficient) in input.terms() {
        let coefficient: Integer = coefficient.resize()?;
        let (mask, coefficient) = if dual {
            (!mask, complement_coefficient(mask, coefficient))
        } else {
            (mask, coefficient)
        };
        for (key, value) in blade(model, mask, coefficient)? {
            let (key, value) = if dual {
                (!key, complement_coefficient(!key, value))
            } else {
                (key, value)
            };
            Form::accumulate(&mut output, key, value)?;
        }
    }
    Form::integral(model, &output)
}

fn wedge_column(
    model: &mut Model,
    result: &mut Terms,
    partial: &Terms,
    column: &[(usize, Integer)],
) -> Result<()> {
    if partial.is_empty() || column.is_empty() {
        return Ok(());
    }
    model.spend((1 + column.len()) * (1 + partial.len()) * WORK_LIMBS)?;
    for (&subset, &coefficient) in partial {
        for &(target, entry) in column {
            let bit = 1 << target;
            if subset & bit == 0 {
                let mut value = coefficient.multiply_wide(entry)?;
                if wedge_negative(subset, bit) {
                    value = value.negate();
                }
                Form::accumulate(result, subset | bit, value)?;
            }
        }
    }
    Ok(())
}

fn complement_coefficient(mask: u64, coefficient: Integer) -> Integer {
    if wedge_negative(mask, !mask) {
        coefficient.negate()
    } else {
        coefficient
    }
}

fn bits(mut mask: u64) -> impl Iterator<Item = usize> {
    std::iter::from_fn(move || {
        if mask == 0 {
            return None;
        }
        let index = mask.trailing_zeros() as usize;
        mask &= mask - 1;
        Some(index)
    })
}

fn wedge_negative(left: u64, right: u64) -> bool {
    let inversions: u32 = bits(right).map(|i| ((left >> i) & !1).count_ones()).sum();
    !inversions.is_multiple_of(2)
}

fn character(model: &mut Model, rank: usize, degree: usize) -> Result<Vec<u128>> {
    let width = 2 * rank + 1;
    model.spend(EVEN_RANK * (degree + 1) * width)?;
    let mut table = vec![vec![0; width]; degree + 1];
    table[0][rank] = 1;
    for i in 0..EVEN_RANK {
        let shift = if i < rank {
            1
        } else if i < 2 * rank {
            -1
        } else {
            0
        };
        for k in (1..=degree.min(i + 1)).rev() {
            for j in 0..width {
                if let Some(target) = j.checked_add_signed(shift).filter(|&j| j < width) {
                    // Every entry counts subsets of at most 64 basis elements.
                    table[k][target] += table[k - 1][j];
                }
            }
        }
    }
    Ok(table
        .pop()
        .expect("the requested exterior degree is present"))
}
