use std::cmp::Ordering;
use std::collections::HashMap;
use std::rc::Rc;

use crate::cobar::{self, Monomial, Support, Word};
use crate::complex::xor_sorted;
use crate::error::{Error, Result, Span};
use crate::model::{Model, Value};
use crate::motivic::Grade;
use crate::profile::Profile;
use crate::real_motivic;

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct Factor {
    u: usize,
    xi: Monomial,
}

impl Factor {
    fn unit() -> Self {
        Self {
            u: 0,
            xi: [].into(),
        }
    }

    fn is_unit(&self) -> bool {
        self.u == 0 && self.xi.is_empty()
    }
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Term {
    vertex: usize,
    factors: Rc<[Factor]>,
}

pub(crate) type Cochain = Rc<[Term]>;
pub(crate) type FactorLiteral = (usize, Vec<usize>, Span);

impl Term {
    fn scalar(s: usize, vertex: usize) -> Self {
        Self {
            vertex,
            factors: vec![Factor::unit(); s].into(),
        }
    }

    pub(crate) fn vertex(&self) -> usize {
        self.vertex
    }

    pub(crate) fn factors(&self) -> impl Iterator<Item = (usize, &[usize])> {
        self.factors
            .iter()
            .map(|factor| (factor.u, factor.xi.as_ref()))
    }

    pub(crate) fn cells(&self) -> usize {
        2 + self
            .factors
            .iter()
            .map(|factor| 2 + factor.xi.len())
            .sum::<usize>()
    }
}

pub(crate) fn classical_degree(grade: Grade) -> Result<usize> {
    let k = grade
        .t
        .checked_sub(i128::from(grade.w))
        .ok_or_else(classical_degree_overflow)?;
    usize::try_from(k).map_err(|_| classical_degree_overflow())
}

fn classical_degree_overflow() -> Error {
    Error::plain("rho-localized degree t-w is outside the supported classical range")
}

fn checked_add(left: usize, right: usize) -> Result<usize> {
    left.checked_add(right)
        .ok_or_else(|| Error::plain("rho-localized exponent is out of range"))
}

pub(crate) fn section_grade(s: usize, t: i128, k: usize) -> Result<Grade> {
    let w = t
        .checked_sub(k as i128)
        .and_then(|w| i64::try_from(w).ok())
        .ok_or_else(|| Error::plain("rho-localized weight is out of range"))?;
    Ok(Grade { s, t, w })
}

pub(crate) fn shift(grade: Grade, power: i128) -> Result<Grade> {
    let t = grade
        .t
        .checked_sub(power)
        .ok_or_else(|| Error::plain("rhoShift topological degree is out of range"))?;
    let w = i128::from(grade.w)
        .checked_sub(power)
        .and_then(|w| i64::try_from(w).ok())
        .ok_or_else(|| Error::plain("rhoShift weight is out of range"))?;
    Ok(Grade { t, w, ..grade })
}

pub(crate) fn literal(
    model: &mut Model,
    t: i128,
    vertex: usize,
    input: Vec<FactorLiteral>,
) -> Result<(Grade, Value)> {
    cobar::check_storage(input.len().saturating_mul(2).saturating_add(2))?;
    let mut k = vertex;
    let mut factors = Vec::with_capacity(input.len());
    for (u, mut xi, span) in input {
        model.spend(1 + xi.len())?;
        xi.truncate(
            xi.iter()
                .rposition(|&power| power != 0)
                .map_or(0, |i| i + 1),
        );
        let factor = Factor { u, xi: xi.into() };
        if factor.is_unit() {
            return Err(Error::new(
                span,
                "rho-localized word factors must have positive degree",
            ));
        }
        k = checked_add(k, u)?;
        factors.push(factor);
    }
    let xi_word: Word = factors.iter().map(|factor| factor.xi.clone()).collect();
    let (s, classical) = cobar::word_metadata(&xi_word)?;
    let grade = section_grade(s, t, checked_add(k, classical)?)?;
    let cochain = persist(
        model,
        vec![Term {
            vertex,
            factors: factors.into(),
        }],
    )?;
    Ok((grade, cochain))
}

pub(crate) fn zero() -> Value {
    Value::RhoLoc([].into())
}

pub(crate) fn persist(model: &mut Model, terms: Vec<Term>) -> Result<Value> {
    let cells = terms.iter().map(Term::cells).sum();
    cobar::check_storage(cells)?;
    model.allocate(cells)?;
    Ok(Value::RhoLoc(terms.into()))
}

pub(crate) fn add(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Value> {
    model.spend(left.iter().chain(right).map(Term::cells).sum())?;
    persist(model, xor_sorted(left, right))
}

pub(crate) fn compare(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Ordering> {
    for (a, b) in left.iter().rev().zip(right.iter().rev()) {
        model.spend(a.cells() + b.cells())?;
        let order = a.cmp(b);
        if order != Ordering::Equal {
            return Ok(order);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

fn toggle(model: &mut Model, support: &mut Support<Term>, term: Term) -> Result<()> {
    let cells = term.cells();
    support.toggle(model, term, cells)
}

pub(crate) fn differential(model: &mut Model, terms: &[Term]) -> Result<Value> {
    let mut image = Support::new();
    for term in terms {
        model.spend(term.cells())?;
        // The first coface sends T_0 to T_0 + U_1. Its term with U_1^0
        // and the last coface both insert a unit arrow.
        for power in submasks(term.vertex).take_while(|&power| power > 0) {
            let first = Factor {
                u: power,
                xi: [].into(),
            };
            let factors = std::iter::once(first)
                .chain(term.factors.iter().cloned())
                .collect();
            toggle(
                model,
                &mut image,
                Term {
                    vertex: term.vertex - power,
                    factors,
                },
            )?;
        }
        for (slot, factor) in term.factors.iter().enumerate() {
            for (left, right) in cobar::coproduct(model, Profile::Full, &factor.xi)? {
                for power in submasks(factor.u) {
                    model.spend(1)?;
                    if power == 0 && left.is_empty() || power == factor.u && right.is_empty() {
                        continue;
                    }
                    let left = Factor {
                        u: power,
                        xi: left.clone(),
                    };
                    let right = Factor {
                        u: factor.u - power,
                        xi: right.clone(),
                    };
                    let factors = term.factors[..slot]
                        .iter()
                        .cloned()
                        .chain([left, right])
                        .chain(term.factors[slot + 1..].iter().cloned())
                        .collect();
                    toggle(
                        model,
                        &mut image,
                        Term {
                            vertex: term.vertex,
                            factors,
                        },
                    )?;
                }
            }
        }
    }
    persist(model, image.into_vec())
}

fn submasks(exponent: usize) -> impl Iterator<Item = usize> {
    std::iter::successors(Some(exponent), move |&power| {
        power.checked_sub(1).map(|previous| previous & exponent)
    })
}

pub(crate) fn homotopy(model: &mut Model, terms: &[Term]) -> Result<Value> {
    let mut image = Support::new();
    for term in terms {
        model.spend(term.cells())?;
        if term.vertex != 0 {
            continue;
        }
        // An eligible unit group arrow has positive U exponent since
        // the word is normalized. Only the first positive U can be used.
        let Some(slot) = term.factors.iter().position(|factor| factor.u != 0) else {
            continue;
        };
        let factor = &term.factors[slot];
        if !factor.xi.is_empty() {
            continue;
        }
        let rest = term.factors[..slot]
            .iter()
            .chain(&term.factors[slot + 1..])
            .cloned()
            .collect();
        let base = Term {
            vertex: 0,
            factors: rest,
        };
        for lift in transport(model, base, slot, factor.u)? {
            toggle(model, &mut image, lift)?;
        }
    }
    persist(model, image.into_vec())
}

pub(crate) fn concat(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Value> {
    let mut product = Support::new();
    for a in left {
        for b in right {
            let cells = a.cells() + b.cells();
            cobar::check_storage(cells)?;
            model.spend(cells)?;
            let factors = a.factors.iter().chain(&*b.factors).cloned().collect();
            let base = Term {
                vertex: a.vertex,
                factors,
            };
            for term in transport(model, base, a.factors.len(), b.vertex)? {
                toggle(model, &mut product, term)?;
            }
        }
    }
    persist(model, product.into_vec())
}

fn transport(model: &mut Model, mut base: Term, prefix: usize, power: usize) -> Result<Vec<Term>> {
    if prefix == 0 {
        base.vertex = checked_add(base.vertex, power)?;
        return Ok(vec![base]);
    }
    let mut terms = vec![base];
    let mut bits = power;
    // Frobenius sends each bit to one vertex or one preceding increment.
    while bits != 0 {
        let bit = 1usize << bits.trailing_zeros();
        bits ^= bit;
        let mut expanded = Support::new();
        for term in terms {
            let vertex = checked_add(term.vertex, bit)?;
            toggle(
                model,
                &mut expanded,
                Term {
                    vertex,
                    factors: term.factors.clone(),
                },
            )?;
            for slot in 0..prefix {
                model.spend(term.cells())?;
                let mut factors = term.factors.to_vec();
                factors[slot].u = checked_add(factors[slot].u, bit)?;
                toggle(
                    model,
                    &mut expanded,
                    Term {
                        vertex: term.vertex,
                        factors: factors.into(),
                    },
                )?;
            }
        }
        terms = expanded.into_vec();
    }
    Ok(terms)
}

pub(crate) fn project(model: &mut Model, terms: &[Term]) -> Result<Value> {
    let mut words = Vec::new();
    let mut cells = 0;
    for term in terms {
        model.spend(term.cells())?;
        if term.vertex == 0 && term.factors.iter().all(|factor| factor.u == 0) {
            let word: Word = term
                .factors
                .iter()
                .map(|factor| factor.xi.clone())
                .collect();
            cells += 1 + word.len();
            cobar::check_storage(cells)?;
            words.push(word);
        }
    }
    model.allocate(cells)?;
    Ok(Value::Cobar(words.into()))
}

pub(crate) fn section(model: &mut Model, words: &[Word]) -> Result<Value> {
    let mut terms = Vec::with_capacity(words.len());
    for word in words {
        model.spend(cobar::word_cells(word))?;
        let factors = word
            .iter()
            .map(|xi| Factor {
                u: 0,
                xi: xi.clone(),
            })
            .collect();
        terms.push(Term { vertex: 0, factors });
    }
    persist(model, terms)
}

pub(crate) fn localize(model: &mut Model, input: &[real_motivic::Term]) -> Result<Value> {
    let mut localized = Support::new();
    for term in input {
        model.spend(term.cells())?;
        let s = term.word.len();
        let factors = term
            .word
            .iter()
            .map(|factor| {
                let last = factor.iter().rposition(|&exponent| exponent >= 2);
                let xi = last.map_or_else(
                    || [].into(),
                    |i| factor[..=i].iter().map(|power| power / 2).collect(),
                );
                Factor { u: 0, xi }
            })
            .collect();
        let mut product = vec![Term {
            vertex: term.tau,
            factors,
        }];
        for (slot, factor) in term.word.iter().enumerate() {
            let Some(last_odd) = factor.iter().rposition(|&power| power & 1 != 0) else {
                continue;
            };
            let mut tau = vec![at_slot(
                s,
                slot,
                Factor {
                    u: 1,
                    xi: [].into(),
                },
            )];
            for (index, &power) in factor[..=last_odd].iter().enumerate() {
                if power & 1 != 0 {
                    product = multiply(model, &product, &tau)?;
                }
                if index < last_odd {
                    tau = next_tau(model, &tau, s, slot, index)?;
                }
            }
        }
        for term in product {
            toggle(model, &mut localized, term)?;
        }
    }
    persist(model, localized.into_vec())
}

fn at_slot(s: usize, slot: usize, factor: Factor) -> Term {
    let mut term = Term::scalar(s, 0);
    Rc::make_mut(&mut term.factors)[slot] = factor;
    term
}

fn next_tau(
    model: &mut Model,
    previous: &[Term],
    s: usize,
    slot: usize,
    index: usize,
) -> Result<Vec<Term>> {
    let mut terms = Support::new();
    // The rho exponent is forced by the tridegree. In these coordinates
    // tau_(i+1) is the square of tau_i plus T_j times xi_(i+1).
    for term in previous {
        let square = multiply_term(model, term, term)?;
        toggle(model, &mut terms, square)?;
    }
    let mut generator = vec![0; index + 1];
    generator[index] = 1;
    let xi = at_slot(
        s,
        slot,
        Factor {
            u: 0,
            xi: generator.into(),
        },
    );
    for term in transport(model, xi, slot + 1, 1)? {
        toggle(model, &mut terms, term)?;
    }
    Ok(terms.into_vec())
}

fn multiply(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Vec<Term>> {
    let mut result = Support::new();
    for a in left {
        for b in right {
            let term = multiply_term(model, a, b)?;
            toggle(model, &mut result, term)?;
        }
    }
    Ok(result.into_vec())
}

fn multiply_term(model: &mut Model, a: &Term, b: &Term) -> Result<Term> {
    debug_assert_eq!(a.factors.len(), b.factors.len());
    let cells = a.cells() + b.cells();
    cobar::check_storage(cells)?;
    model.spend(cells)?;
    let factors = a
        .factors
        .iter()
        .zip(&*b.factors)
        .map(|(a, b)| {
            let u = checked_add(a.u, b.u)?;
            let xi = if a.xi.is_empty() {
                b.xi.clone()
            } else if b.xi.is_empty() {
                a.xi.clone()
            } else {
                (0..a.xi.len().max(b.xi.len()))
                    .map(|index| {
                        checked_add(
                            a.xi.get(index).copied().unwrap_or(0),
                            b.xi.get(index).copied().unwrap_or(0),
                        )
                    })
                    .collect::<Result<Monomial>>()?
            };
            Ok(Factor { u, xi })
        })
        .collect::<Result<Rc<[Factor]>>>()?;
    Ok(Term {
        vertex: checked_add(a.vertex, b.vertex)?,
        factors,
    })
}

pub(crate) fn basis(model: &mut Model, grade: Grade, limit: usize) -> Result<Vec<Term>> {
    if grade.t < i128::from(grade.w) {
        return Ok(Vec::new());
    }
    let k = classical_degree(grade)?;
    let s = grade.s;
    let Some(excess) = k.checked_sub(s) else {
        return Ok(Vec::new());
    };
    if limit == 0 {
        return Err(cobar::basis_limit(limit));
    }
    if s == 0 {
        let term = Term::scalar(0, k);
        model.allocate(term.cells())?;
        return Ok(vec![term]);
    }
    // Even at the minimum degree, each slot can hold U or xi_1.
    // Varying the vertex exponent also gives excess + 1 choices.
    if s > limit.ilog2() as usize || excess >= limit {
        return Err(cobar::basis_limit(limit));
    }
    let mut degrees = vec![1; s + 1];
    degrees[s] += excess;
    let mut words = Vec::new();
    let mut classical = HashMap::new();
    let mut catalog = HashMap::new();
    let mut cells = 0;
    let mut catalog_cells = 0;
    loop {
        model.spend(s + 1)?;
        for &degree in &degrees[1..] {
            if let std::collections::hash_map::Entry::Vacant(entry) = catalog.entry(degree) {
                let mut factors = Vec::new();
                for q in 0..=degree {
                    if let std::collections::hash_map::Entry::Vacant(entry) = classical.entry(q) {
                        let powers = if q == 0 {
                            vec![Monomial::from([])]
                        } else {
                            cobar::basis(model, 1, q, limit)?
                                .into_iter()
                                .map(|word| word[0].clone())
                                .collect()
                        };
                        entry.insert(powers);
                    }
                    for xi in &classical[&q] {
                        model.spend(2 + xi.len())?;
                        catalog_cells += 2 + xi.len();
                        cobar::check_storage(catalog_cells + cells)?;
                        factors.push(Factor {
                            u: degree - q,
                            xi: xi.clone(),
                        });
                    }
                }
                entry.insert(factors);
            }
        }
        let options: Vec<_> = degrees[1..].iter().map(|d| catalog[d].as_slice()).collect();
        let mut digits = vec![0; s];
        loop {
            if words.len() == limit {
                return Err(cobar::basis_limit(limit));
            }
            let factors = options
                .iter()
                .zip(&digits)
                .map(|(choices, &i)| choices[i].clone())
                .collect();
            let word = Term {
                vertex: degrees[0] - 1,
                factors,
            };
            model.spend(word.cells())?;
            cells += word.cells();
            cobar::check_storage(catalog_cells + cells)?;
            words.push(word);
            if !cobar::next_choices(&options, &mut digits) {
                break;
            }
        }
        if !cobar::next_degrees(&mut degrees) {
            break;
        }
    }
    model.spend(cells * (1 + words.len().ilog2() as usize))?;
    model.allocate(cells)?;
    words.sort_unstable();
    Ok(words)
}
