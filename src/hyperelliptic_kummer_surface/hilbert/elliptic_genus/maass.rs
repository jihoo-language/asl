use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_ore::gcd;
use crate::siegel::{discriminant, reduced_jacobi_degree};
use crate::syntax::Label;

use super::super::Magnitude;
use super::{
    add, coefficient, decimal, eta_product, jacobi_minus_two, jacobi_zero, multiply, overflow,
    storage, Cutoff, Index, Integer, INTEGER_WORDS, MAX_PIVOT_CELLS,
};

const ZERO: Integer = Integer {
    magnitude: Magnitude::ZERO,
    negative: false,
};
type Series = Vec<Integer>;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Family {
    MinusTwo,
    Zero,
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Monomial(pub(crate) Family, pub(crate) [u128; 3]);

pub(crate) enum Command {
    Declare(Label, u128, Vec<(Monomial, i128)>),
    Jacobi(Label, u128, i128),
    Lift(Label, Index),
    Chi12(Index),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<Option<String>> {
        let (label, value) = match self {
            Self::Declare(name, weight, terms) => {
                model.register_siegel_maass(name, weight, terms)?;
                return Ok(None);
            }
            Self::Jacobi(name, n, r) => (
                format!("rvSiegelMaassJacobi[{}, {n}, {r}]", name.text),
                model.siegel_maass(&name)?.jacobi(model, n, r)?,
            ),
            Self::Lift(name, index) => (
                format!("rvSiegelMaassCoeff[{}, {index}]", name.text),
                model.siegel_maass(&name)?.lift(model, index)?,
            ),
            Self::Chi12(index) => {
                let term = (Monomial(Family::Zero, [0, 0, 1]), Integer::from_i128(1));
                let form = Form {
                    weight: 12,
                    terms: BTreeMap::from([term]),
                };
                (format!("rvIgusaChi12[{index}]"), form.lift(model, index)?)
            }
        };
        Ok(Some(format!("{label} = {}", decimal(value))))
    }
}

pub(crate) struct Form {
    weight: u128,
    terms: BTreeMap<Monomial, Integer>,
}

impl Form {
    pub(crate) fn new(
        model: &mut Model,
        weight: u128,
        terms: Vec<(Monomial, i128)>,
    ) -> Result<Self> {
        if weight < 10 || !weight.is_multiple_of(2) {
            return Err(Error::plain("Maaß lift requires even weight at least 10"));
        }
        let mut canonical = BTreeMap::new();
        for (term @ Monomial(family, [a, b, d]), scale) in terms {
            if d == 0 {
                return Err(Error::plain(
                    "Jacobi cusp monomials require a positive Delta exponent",
                ));
            }
            let first_weight = match family {
                Family::MinusTwo => 10,
                Family::Zero => 12,
            };
            let total = (d - 1)
                .checked_mul(12)
                .and_then(|total| total.checked_add(first_weight))
                .and_then(|total| total.checked_add(a.checked_mul(4)?))
                .and_then(|total| total.checked_add(b.checked_mul(6)?));
            if total != Some(weight) {
                return Err(Error::plain(
                    "Jacobi cusp monomial does not have the declared weight",
                ));
            }
            add(model, &mut canonical, term, Integer::from_i128(scale))?;
        }
        model.allocate(canonical.len() * (INTEGER_WORDS + 12))?;
        Ok(Self {
            weight,
            terms: canonical,
        })
    }

    fn jacobi(&self, model: &mut Model, n: u128, r: i128) -> Result<Integer> {
        let Some((degree, parity)) = self.reduced_bound(discriminant(n, 1, r)?, r)? else {
            return Ok(ZERO);
        };
        Ok(self.reduced_series(model, degree)?[degree][parity])
    }

    fn lift(&self, model: &mut Model, Index(m, n, r): Index) -> Result<Integer> {
        let discriminant = discriminant(m, n, r)?;
        let Some((degree, _)) = self.reduced_bound(discriminant, r)? else {
            return Ok(ZERO);
        };
        let reduced = self.reduced_series(model, degree)?;
        let content = gcd(gcd(m, n), r.unsigned_abs());
        let mut value = ZERO;
        for a in 1..=content.isqrt() {
            model.spend(1)?;
            if !content.is_multiple_of(a) {
                continue;
            }
            let b = content / a;
            for d in [Some(a), (a != b).then_some(b)].into_iter().flatten() {
                let (quotient, first_remainder) = discriminant.magnitude.divide(model, d)?;
                let (magnitude, second_remainder) = quotient.divide(model, d)?;
                debug_assert_eq!((first_remainder, second_remainder), (0, 0));
                let parity = (r.unsigned_abs() / d % 2) as i128;
                let degree = reduced_jacobi_degree(
                    Integer {
                        magnitude,
                        negative: false,
                    },
                    parity,
                )?
                .ok_or_else(storage)? as usize;
                let coefficient = reduced[degree][parity as usize];
                if !coefficient.is_zero() {
                    let term = scale_by_power(model, coefficient, d, self.weight - 1)?;
                    value = value.add(term).map_err(|_| overflow())?;
                }
            }
        }
        Ok(value)
    }

    fn reduced_bound(&self, discriminant: Integer, r: i128) -> Result<Option<(usize, usize)>> {
        if self.terms.is_empty() || discriminant.is_zero() || discriminant.negative {
            return Ok(None);
        }
        let parity = r.rem_euclid(2);
        let degree = reduced_jacobi_degree(discriminant, parity)?.ok_or_else(storage)?;
        if self.terms.keys().all(|Monomial(_, [_, _, d])| *d > degree) {
            return Ok(None);
        }
        let degree = usize::try_from(degree).map_err(|_| storage())?;
        // Bound the dense polynomial buffers alongside the sparse Jacobi expansion.
        if degree.saturating_add(1).saturating_mul(16 * INTEGER_WORDS) > MAX_PIVOT_CELLS {
            return Err(storage());
        }
        Ok(Some((degree, parity as usize)))
    }

    fn reduced_series(&self, model: &mut Model, degree: usize) -> Result<Vec<[Integer; 2]>> {
        let terms = self
            .terms
            .iter()
            .filter(|(Monomial(_, [_, _, d]), _)| *d <= degree as u128);
        let minimum = terms
            .clone()
            .map(|(Monomial(_, [_, _, d]), _)| *d as usize)
            .min()
            .unwrap();
        let cutoff = Cutoff::new(model, 0, (degree - minimum) as u128)?;
        let max = cutoff.q;
        let e4 = if terms.clone().any(|(Monomial(_, [a, _, _]), _)| *a != 0) {
            eisenstein(model, max, 3, 240)?
        } else {
            Vec::new()
        };
        let e6 = if terms.clone().any(|(Monomial(_, [_, b, _]), _)| *b != 0) {
            eisenstein(model, max, 5, -504)?
        } else {
            Vec::new()
        };
        let delta = eta_product(model, max, 24)?;
        let mut result = vec![[ZERO; 2]; degree + 1];
        for family in [Family::MinusTwo, Family::Zero] {
            let mut selected = terms
                .clone()
                .filter(|(Monomial(kind, _), _)| *kind == family)
                .peekable();
            if selected.peek().is_none() {
                continue;
            }
            let weak = weak_jacobi_rows(model, cutoff, family)?;
            for (&Monomial(_, [a, b, d]), &scale) in selected {
                let d = d as usize;
                let bound = degree - d;
                let mut modular = vec![scale];
                for (base, exponent) in [(&e4, a), (&e6, b), (&delta, d as u128)] {
                    modular = multiply_power(model, modular, base, exponent, bound)?;
                }
                for (lag, &scalar) in modular.iter().enumerate() {
                    if scalar.is_zero() {
                        continue;
                    }
                    for (j, row) in weak.iter().take(bound - lag + 1).enumerate() {
                        for (output, &entry) in result[d + lag + j].iter_mut().zip(row) {
                            let term = multiply(model, scalar, entry)?;
                            *output = output.add(term).map_err(|_| overflow())?;
                        }
                    }
                }
            }
        }
        Ok(result)
    }
}

fn weak_jacobi_rows(
    model: &mut Model,
    cutoff: Cutoff,
    family: Family,
) -> Result<Vec<[Integer; 2]>> {
    match family {
        Family::MinusTwo => {
            let product = jacobi_minus_two(model, cutoff)?;
            Ok((0..=cutoff.q)
                .map(|q| {
                    [
                        coefficient(&product, (0, q, 0)),
                        coefficient(&product, (0, q, 1)),
                    ]
                })
                .collect())
        }
        Family::Zero => {
            let seed = jacobi_zero(model, cutoff.q)?;
            Ok((0..=cutoff.q)
                .map(|q| {
                    [
                        *seed.get(&(q, 0)).unwrap_or(&ZERO),
                        *seed.get(&(q, 1)).unwrap_or(&ZERO),
                    ]
                })
                .collect())
        }
    }
}

fn eisenstein(model: &mut Model, bound: usize, power: u32, scale: i128) -> Result<Series> {
    let mut series = vec![ZERO; bound + 1];
    series[0] = Integer::from_i128(1);
    for d in 1..=bound {
        let term = Integer::from_i128((d as i128).pow(power) * scale);
        for n in (d..=bound).step_by(d) {
            model.spend(1)?;
            series[n] = series[n].add(term).map_err(|_| overflow())?;
        }
    }
    Ok(series)
}

fn multiply_power(
    model: &mut Model,
    mut out: Series,
    input: &[Integer],
    mut exponent: u128,
    bound: usize,
) -> Result<Series> {
    if exponent == 0 {
        return Ok(out);
    }
    model.spend(bound + 1)?;
    let mut power = input[..=bound].to_vec();
    while exponent != 0 {
        if exponent % 2 == 1 {
            out = convolve(model, &out, &power, bound)?;
        }
        exponent /= 2;
        if exponent != 0 {
            power = convolve(model, &power, &power, bound)?;
        }
    }
    Ok(out)
}

fn convolve(
    model: &mut Model,
    left: &[Integer],
    right: &[Integer],
    bound: usize,
) -> Result<Series> {
    let mut out = vec![ZERO; (left.len() + right.len() - 1).min(bound + 1)];
    for (i, &a) in left.iter().enumerate() {
        if a.is_zero() {
            continue;
        }
        for (j, &b) in right.iter().take(out.len() - i).enumerate() {
            let term = multiply(model, a, b)?;
            out[i + j] = out[i + j].add(term).map_err(|_| overflow())?;
        }
    }
    Ok(out)
}

fn scale_by_power(
    model: &mut Model,
    mut value: Integer,
    base: u128,
    mut exponent: u128,
) -> Result<Integer> {
    if base == 1 {
        return Ok(value);
    }
    let mut power = Integer {
        magnitude: Magnitude::from_small(base),
        negative: false,
    };
    while exponent != 0 {
        if exponent % 2 == 1 {
            value = multiply(model, value, power)?;
        }
        exponent /= 2;
        if exponent != 0 {
            power = multiply(model, power, power)?;
        }
    }
    Ok(value)
}
