use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use crate::complex::{xor_sorted, MAX_PIVOT_CELLS};
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::rv_hochster::{self, Summand};
use crate::rv_ideal::{self, Ideal};
use crate::rv_koszul::{self, Grade as TorGrade};
use crate::rv_poly::{self, Monomial};
use crate::rv_t_transfer::{self, Operation};

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) ideal: Rc<str>,
    pub(crate) b: i64,
    pub(crate) weight: Monomial,
}

impl Grade {
    pub(crate) fn from_tor(grade: TorGrade, b: i64) -> Self {
        Self {
            ideal: grade.ideal,
            b,
            weight: grade.multidegree.unwrap(),
        }
    }

    pub(crate) fn differential(&self) -> Result<Self> {
        let b = self.b.checked_sub(1).ok_or_else(degree_limit)?;
        Ok(Self { b, ..self.clone() })
    }

    pub(crate) fn connes(&self) -> Result<Self> {
        let b = self.b.checked_add(1).ok_or_else(degree_limit)?;
        Ok(Self { b, ..self.clone() })
    }
}

#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Factor {
    pub(crate) weight: Monomial,
    pub(crate) s: usize,
    pub(crate) generator: usize,
}

impl Factor {
    fn unit(&self) -> bool {
        self.s == 0 && self.weight.iter().all(|&exponent| exponent == 0)
    }

    fn tor_grade(&self, ideal: Rc<str>) -> TorGrade {
        TorGrade {
            ideal,
            s: self.s,
            t: self.weight.iter().sum(),
            multidegree: Some(self.weight.clone()),
        }
    }
}

pub(crate) type Word = Rc<[Factor]>;
pub(crate) type Support = Rc<[Word]>;

fn degree_limit() -> Error {
    Error::plain("rank Hochschild degree is out of range")
}

fn basis_limit() -> Error {
    Error::plain("rank Hochschild basis exceeds the enumeration limit")
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "rank Hochschild calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

pub(crate) fn word_cells(word: &[Factor]) -> usize {
    1 + word
        .iter()
        .map(|factor| 3 + factor.weight.len())
        .sum::<usize>()
}

fn support_cells(words: &[Word]) -> usize {
    1 + words.iter().map(|word| word_cells(word)).sum::<usize>()
}

pub(crate) fn persist(model: &mut Model, words: Vec<Word>) -> Result<Support> {
    model.allocate(support_cells(&words))?;
    Ok(words.into())
}

pub(crate) fn add(model: &mut Model, a: &[Word], b: &[Word]) -> Result<Support> {
    model.spend(support_cells(a) + support_cells(b))?;
    persist(model, xor_sorted(a, b))
}

pub(crate) fn compare(model: &mut Model, a: &[Word], b: &[Word]) -> Result<Ordering> {
    model.spend(support_cells(a) + support_cells(b))?;
    Ok(a.iter().rev().cmp(b.iter().rev()))
}

pub(crate) fn length_part(model: &mut Model, len: usize, words: &[Word]) -> Result<Support> {
    model.spend(words.len())?;
    persist(
        model,
        words
            .iter()
            .filter(|word| word.len() - 1 == len)
            .cloned()
            .collect(),
    )
}

pub(crate) fn word_type(
    model: &Model,
    inputs: impl ExactSizeIterator<Item = TypeId>,
) -> Result<Grade> {
    let length = inputs
        .len()
        .checked_sub(1)
        .ok_or_else(|| Error::plain("rvHHWord requires a coefficient class"))?;
    let mut combined: Option<TorGrade> = None;
    for (i, ty) in inputs.enumerate() {
        let Kind::RVTor(grade) = model.kind(ty) else {
            return Err(Error::plain("rvHHWord requires RVMTor classes"));
        };
        if grade.multidegree.is_none() {
            return Err(Error::plain("rvHHWord requires RVMTor classes"));
        }
        if i > 0 && grade.t == 0 {
            return Err(Error::plain(
                "rvHHWord requires positive multidegree in every reduced input",
            ));
        }
        combined = Some(if let Some(sum) = combined {
            sum.product(grade)?
        } else {
            grade.clone()
        });
    }
    let combined = combined.unwrap();
    let b = combined
        .s
        .checked_add(length)
        .and_then(|b| b.try_into().ok())
        .ok_or_else(degree_limit)?;
    Ok(Grade::from_tor(combined, b))
}

pub(crate) fn word(model: &mut Model, inputs: &[(TypeId, Rc<[usize]>)]) -> Result<Support> {
    if inputs.iter().any(|(_, class)| class.is_empty()) {
        return Ok([].into());
    }
    let factors = inputs
        .iter()
        .map(|(ty, _)| {
            let Kind::RVTor(grade) = model.kind(*ty) else {
                unreachable!()
            };
            Factor {
                weight: grade.multidegree.clone().unwrap(),
                s: grade.s,
                generator: 0,
            }
        })
        .collect::<Vec<_>>();
    let count = inputs
        .iter()
        .try_fold(1usize, |count, (_, class)| count.checked_mul(class.len()))
        .ok_or_else(basis_limit)?;
    let cells = count.saturating_mul(word_cells(&factors));
    check_storage(cells)?;
    model.spend(cells)?;
    let mut terms = Vec::with_capacity(count);
    for mut code in 0..count {
        let mut word = factors.clone();
        for (factor, (_, class)) in word.iter_mut().zip(inputs).rev() {
            factor.generator = class[code % class.len()];
            code /= class.len();
        }
        terms.push(word.into());
    }
    persist(model, terms)
}

#[derive(Default)]
struct Sum {
    terms: BTreeSet<Word>,
    cells: usize,
}

impl Sum {
    fn toggle(&mut self, model: &mut Model, term: Vec<Factor>) -> Result<()> {
        let cost = word_cells(&term);
        model.spend(cost * (1 + self.terms.len().max(1).ilog2() as usize))?;
        let word: Word = term.into();
        if self.terms.remove(&word) {
            self.cells -= cost;
        } else {
            self.cells += cost;
            check_storage(self.cells)?;
            self.terms.insert(word);
        }
        Ok(())
    }

    fn finish(self) -> Vec<Word> {
        self.terms.into_iter().collect()
    }
}

pub(crate) struct Differential<'a> {
    ideal: Rc<str>,
    squarefree: Option<&'a [usize]>,
    width: usize,
    ranks: BTreeMap<Monomial, Option<Vec<Summand>>>,
    cache: BTreeMap<Word, Word>,
    cells: usize,
}

impl<'a> Differential<'a> {
    pub(crate) fn new(model: &mut Model, grade: &Grade, ideal: &'a Ideal) -> Result<Self> {
        Ok(Self {
            ideal: grade.ideal.clone(),
            squarefree: ideal.squarefree_masks(model)?,
            width: grade.weight.len(),
            ranks: BTreeMap::new(),
            cache: BTreeMap::new(),
            cells: 0,
        })
    }

    fn allows(&mut self, model: &mut Model, grade: &TorGrade) -> Result<bool> {
        if grade.s > self.width.min(grade.t) {
            return Ok(false);
        }
        let Some(masks) = self.squarefree else {
            return Ok(true);
        };
        let weight = grade.multidegree.as_ref().unwrap();
        model.spend(1 + weight.len())?;
        if weight.iter().any(|&power| power > 1) {
            return Ok(false);
        }
        if let Some(cached) = self.ranks.get(weight) {
            return Ok(cached
                .as_ref()
                .is_none_or(|ranks| ranks.iter().any(|summand| summand.s == grade.s)));
        }
        // This is a vanishing shortcut for small collapsed complexes. When
        // it is inconclusive, compute the product directly on its support.
        let ranks =
            rv_hochster::squarefree_ranks(model, masks, rv_ideal::support_mask(weight, 1), 12)?;
        let allowed = ranks
            .as_ref()
            .is_none_or(|ranks| ranks.iter().any(|summand| summand.s == grade.s));
        self.cells += weight.len() + 2 * ranks.as_ref().map_or(0, Vec::len);
        check_storage(self.cells)?;
        self.ranks.insert(weight.clone(), ranks);
        Ok(allowed)
    }

    fn evaluate(&mut self, model: &mut Model, factors: &[Factor]) -> Result<Word> {
        let cost = word_cells(factors);
        model.spend(cost * (1 + self.cache.len().max(1).ilog2() as usize))?;
        if let Some(output) = self.cache.get(factors) {
            return Ok(output.clone());
        }
        let mut grade = factors[0].tor_grade(self.ideal.clone());
        for factor in &factors[1..] {
            grade = grade.product(&factor.tor_grade(self.ideal.clone()))?;
        }
        grade.s = grade
            .s
            .checked_add(factors.len() - 2)
            .ok_or_else(degree_limit)?;
        let result: Word = if !self.allows(model, &grade)? {
            [].into()
        } else {
            let target = model.rv_complex_type(grade.clone(), true)?;
            let inputs = factors
                .iter()
                .map(|factor| {
                    let ty = model.rv_complex_type(factor.tor_grade(self.ideal.clone()), true)?;
                    Ok((ty, [factor.generator].into()))
                })
                .collect::<Result<Vec<_>>>()?;
            let Value::Class(class) =
                rv_t_transfer::evaluate(model, Operation::Higher, target, &inputs)?
            else {
                unreachable!()
            };
            class
                .iter()
                .map(|&generator| Factor {
                    weight: grade.multidegree.clone().unwrap(),
                    s: grade.s,
                    generator,
                })
                .collect()
        };
        self.cells += cost + word_cells(&result);
        check_storage(self.cells)?;
        self.cache.insert(factors.into(), result.clone());
        Ok(result)
    }

    pub(crate) fn apply(&mut self, model: &mut Model, words: &[Word]) -> Result<Vec<Word>> {
        let mut result = Sum::default();
        for word in words {
            let len = word.len();
            model.spend(1 + self.width)?;
            let coefficient_is_unit = word[0].unit();
            for start in 0..len {
                let mut block = Vec::new();
                block.push(word[start].clone());
                for size in 2..=len {
                    model.spend(1)?;
                    let end = start + size;
                    block.push(word[(end - 1) % len].clone());
                    let wraps = end > len;
                    let contains_unit = coefficient_is_unit && (start == 0 || wraps);
                    let factors: Word = if contains_unit {
                        if size > 2 {
                            continue;
                        }
                        let other = if start == 0 { &word[1] } else { &word[start] };
                        [other.clone()].into()
                    } else {
                        self.evaluate(model, &block)?
                    };
                    for factor in factors.iter() {
                        let mut term = Vec::with_capacity(len + 1 - size);
                        if start == 0 {
                            term.push(factor.clone());
                            term.extend_from_slice(&word[end..]);
                        } else if wraps {
                            term.push(factor.clone());
                            term.extend_from_slice(&word[end - len..start]);
                        } else {
                            term.extend_from_slice(&word[..start]);
                            term.push(factor.clone());
                            term.extend_from_slice(&word[end..]);
                        }
                        result.toggle(model, term)?;
                    }
                }
            }
        }
        Ok(result.finish())
    }
}

pub(crate) fn differential(model: &mut Model, grade: &Grade, words: &[Word]) -> Result<Vec<Word>> {
    if words.is_empty() {
        return Ok(Vec::new());
    }
    let ideal = model.rv_ideal_named(&grade.ideal).unwrap();
    Differential::new(model, grade, &ideal)?.apply(model, words)
}

pub(crate) fn connes(model: &mut Model, grade: &Grade, words: &[Word]) -> Result<Vec<Word>> {
    let mut result = Sum::default();
    let unit = Factor {
        weight: vec![0; grade.weight.len()].into(),
        s: 0,
        generator: 0,
    };
    for word in words {
        model.spend(word_cells(word))?;
        if !word[0].unit() {
            for start in 0..word.len() {
                let mut term = Vec::with_capacity(word.len() + 1);
                term.push(unit.clone());
                term.extend_from_slice(&word[start..]);
                term.extend_from_slice(&word[..start]);
                result.toggle(model, term)?;
            }
        }
    }
    Ok(result.finish())
}

fn eligible_factors(model: &mut Model, grade: &Grade, b: usize) -> Result<Vec<Factor>> {
    let ideal = model.rv_ideal_named(&grade.ideal).unwrap();
    let generators = rv_koszul::monomial_initial(model, &ideal)?
        .ok_or_else(|| Error::plain("rank Hochschild chains require a monomial rank ideal"))?;
    let zero: Monomial = vec![0; grade.weight.len()].into();
    if generators.iter().any(|&powers| powers == zero.as_ref()) {
        return Ok(Vec::new());
    }
    let mut multiples = BTreeSet::from([zero.clone()]);
    let mut cells = 1 + zero.len();
    for powers in generators {
        model.spend(1 + powers.len())?;
        if powers.iter().zip(grade.weight.iter()).any(|(a, b)| a > b) {
            continue;
        }
        let seen = multiples.iter().cloned().collect::<Vec<_>>();
        for other in seen {
            model.spend((1 + powers.len()) * (1 + multiples.len().max(1).ilog2() as usize))?;
            let lcm: Monomial = powers
                .iter()
                .zip(other.iter())
                .map(|(&a, &b)| a.max(b))
                .collect();
            if multiples.insert(lcm) {
                cells += 1 + powers.len();
                check_storage(cells)?;
            }
        }
    }

    // The Taylor resolution has shifts at lcms of subsets of monomial
    // generators. Only these weights can support Tor.
    let mut factors = vec![Factor {
        weight: zero,
        s: 0,
        generator: 0,
    }];
    for weight in multiples.into_iter().skip(1) {
        let t = rv_poly::degree(&weight)?;
        for s in 1..=b.min(t).min(weight.len()) {
            model.spend(1)?;
            let ty = model.rv_complex_type(
                TorGrade {
                    ideal: grade.ideal.clone(),
                    s,
                    t,
                    multidegree: Some(weight.clone()),
                },
                true,
            )?;
            let rank = model.rv_tor(ty)?.rank();
            cells += (3 + weight.len()) * rank;
            check_storage(cells)?;
            factors.extend((0..rank).map(|generator| Factor {
                weight: weight.clone(),
                s,
                generator,
            }));
        }
    }
    Ok(factors)
}

pub(crate) fn basis(model: &mut Model, grade: &Grade, limit: usize) -> Result<Vec<Word>> {
    basis_at(model, grade, grade.b.into(), limit)
}

pub(crate) fn basis_at(
    model: &mut Model,
    grade: &Grade,
    degree: i128,
    limit: usize,
) -> Result<Vec<Word>> {
    let bound = 2 * rv_poly::degree(&grade.weight)? as i128;
    if degree < 0 || degree > bound {
        return Ok(Vec::new());
    }
    let b = usize::try_from(degree).map_err(|_| degree_limit())?;
    let factors = eligible_factors(model, grade, b)?;
    let mut remaining = grade.weight.to_vec();
    let mut degree = b;
    let mut active: Vec<usize> = Vec::new();
    let mut next = 0;
    let mut basis = Vec::new();
    let mut cells = 0;
    loop {
        if next == factors.len() {
            let Some(previous) = active.pop() else {
                break;
            };
            let factor = &factors[previous];
            degree += factor.s + usize::from(!active.is_empty());
            for (remaining, value) in remaining.iter_mut().zip(factor.weight.iter()) {
                *remaining += value;
            }
            next = previous + 1;
            continue;
        }
        let factor = &factors[next];
        let reduced = !active.is_empty();
        let cost = factor.s + usize::from(reduced);
        model.spend(1 + remaining.len())?;
        if cost > degree || factor.weight.iter().zip(&remaining).any(|(a, b)| a > b) {
            next += 1;
            continue;
        }
        for (remaining, value) in remaining.iter_mut().zip(factor.weight.iter()) {
            *remaining -= value;
        }
        degree -= cost;
        active.push(next);
        check_storage(active.len())?;
        if remaining.iter().all(|&v| v == 0) {
            if degree == 0 {
                if basis.len() == limit {
                    return Err(basis_limit());
                }
                let word: Word = active.iter().map(|&i| factors[i].clone()).collect();
                cells += word_cells(&word);
                check_storage(cells)?;
                basis.push(word);
            }
            next = factors.len();
        } else if degree < 2 {
            next = factors.len();
        } else {
            // The initial (unit) factor is unavailable in a reduced slot.
            next = 1;
        }
    }
    Ok(basis)
}
