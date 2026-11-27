use std::cmp::Ordering;
use std::collections::{hash_map::Entry, HashMap};
use std::ops::RangeInclusive;
use std::rc::Rc;

use crate::cobar::{self, Monomial, Support, Word};
use crate::complex::xor_sorted;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::motivic::{self, Grade};
use crate::odd_cobar::FactorLiteral;
use crate::profile::Profile;

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Term {
    pub(crate) rho: usize,
    pub(crate) tau: usize,
    pub(crate) word: Word,
}

pub(crate) type Cochain = Rc<[Term]>;

// At a fixed tridegree the word determines both scalar exponents. Encode
// its factors once, and keep only their indices in each sparse row key.
#[derive(Default)]
pub(crate) struct Rows {
    factors: HashMap<Monomial, usize>,
    words: HashMap<Box<[usize]>, usize>,
    cells: usize,
}

impl Rows {
    pub(crate) fn column(&mut self, model: &mut Model, terms: &[Term]) -> Result<Vec<usize>> {
        let mut support = Vec::with_capacity(terms.len());
        for term in terms {
            support.push(self.coordinate(model, &term.word)?);
        }
        crate::may::sort_support(model, &mut support)?;
        Ok(support)
    }

    fn coordinate(&mut self, model: &mut Model, word: &Word) -> Result<usize> {
        let size = 1 + word.len();
        model.spend(size)?;
        let mut key = Vec::with_capacity(word.len());
        for factor in word.iter() {
            model.spend(1 + factor.len())?;
            let next = self.factors.len();
            let index = match self.factors.entry(factor.clone()) {
                Entry::Occupied(entry) => *entry.get(),
                Entry::Vacant(entry) => {
                    let cells = self.cells + 2 + factor.len();
                    cobar::check_storage(cells)?;
                    self.cells = cells;
                    *entry.insert(next)
                }
            };
            key.push(index);
        }
        let next = self.words.len();
        match self.words.entry(key.into_boxed_slice()) {
            Entry::Occupied(entry) => Ok(*entry.get()),
            Entry::Vacant(entry) => {
                let cells = self.cells + size;
                cobar::check_storage(cells)?;
                self.cells = cells;
                Ok(*entry.insert(next))
            }
        }
    }
}

impl Term {
    fn scalar(s: usize, rho: usize, tau: usize) -> Self {
        Self {
            rho,
            tau,
            word: vec![[].into(); s].into(),
        }
    }

    fn arrow(s: usize, slot: usize, index: usize) -> Self {
        let mut factors = vec![[].into(); s];
        factors[slot] = generator(index);
        Self {
            rho: 0,
            tau: 0,
            word: factors.into(),
        }
    }

    pub(crate) fn cells(&self) -> usize {
        2 + cobar::word_cells(&self.word)
    }

    fn normalized(&self) -> bool {
        self.word.iter().all(|factor| !factor.is_empty())
    }
}

fn generator(index: usize) -> Monomial {
    let mut factor = vec![0; index + 1];
    factor[index] = 1;
    factor.into()
}

fn checked_sum(left: usize, right: usize) -> Result<usize> {
    left.checked_add(right)
        .ok_or_else(|| Error::plain("real-motivic exponent is out of range"))
}

pub(crate) fn shifted(grade: Grade, power: usize) -> Result<Grade> {
    let t = grade
        .t
        .checked_sub(power as i128)
        .ok_or_else(|| Error::plain("real-motivic topological degree is out of range"))?;
    Ok(Grade {
        t,
        ..grade.shift(power)?
    })
}

pub(crate) fn literal(
    model: &mut Model,
    rho: usize,
    tau: usize,
    factors: Vec<FactorLiteral>,
) -> Result<(Grade, Value)> {
    cobar::check_storage(factors.len().saturating_mul(2).saturating_add(3))?;
    let mut intrinsic = Vec::with_capacity(factors.len());
    let mut squarefree = Vec::with_capacity(factors.len());
    let mut powers = Vec::new();
    let mut odd_count = 0;
    for (slot, (xi, tau_powers, span)) in factors.into_iter().enumerate() {
        let width = xi.len().max(tau_powers.len());
        model.spend(1 + width)?;
        let mut raw = Vec::with_capacity(width);
        let mut canonical = Vec::with_capacity(width);
        for index in 0..width {
            let a = xi.get(index).copied().unwrap_or(0);
            let b = tau_powers.get(index).copied().unwrap_or(0);
            let two_a = checked_sum(a, a)?;
            raw.push(checked_sum(two_a, b)?);
            canonical.push(checked_sum(two_a, b & 1)?);
            odd_count = checked_sum(odd_count, b)?;
            if b > 1 {
                powers.push((slot, index, b & !1));
            }
        }
        trim(&mut raw);
        trim(&mut canonical);
        if raw.is_empty() {
            return Err(Error::new(
                span,
                "real-motivic word factors must have positive degree",
            ));
        }
        intrinsic.push(Monomial::from(raw));
        squarefree.push(Monomial::from(canonical));
    }
    let intrinsic: Word = intrinsic.into();
    cobar::check_storage(cobar::word_cells(&intrinsic))?;
    let (s, top) = cobar::word_metadata(&intrinsic)?;
    let w = i64::try_from(((top - odd_count) / 2) as i128 - rho as i128 - tau as i128)
        .map_err(|_| Error::plain("real-motivic weight is out of range"))?;
    let grade = Grade {
        s,
        t: top as i128 - rho as i128,
        w,
    };
    let mut product = vec![Term {
        rho,
        tau,
        word: squarefree.into(),
    }];
    for (slot, index, power) in powers {
        product = multiply_power(
            model,
            product,
            vec![Term::arrow(s, slot, index)],
            power,
            u128::MAX,
        )?;
    }
    // Each term introduced by an arrow generator has positive arrow degree.
    let value = persist(model, product)?;
    Ok((grade, value))
}

pub(crate) fn zero() -> Value {
    Value::RMCobar([].into())
}

pub(crate) fn persist(model: &mut Model, terms: Vec<Term>) -> Result<Value> {
    let cells = terms.iter().map(Term::cells).sum::<usize>();
    cobar::check_storage(cells)?;
    model.allocate(cells)?;
    Ok(Value::RMCobar(terms.into()))
}

pub(crate) fn add(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Vec<Term>> {
    model.spend(left.iter().chain(right).map(Term::cells).sum())?;
    let result = xor_sorted(left, right);
    cobar::check_storage(result.iter().map(Term::cells).sum())?;
    Ok(result)
}

pub(crate) fn compare(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Ordering> {
    for (left, right) in left.iter().rev().zip(right.iter().rev()) {
        model.spend(left.cells() + right.cells())?;
        let order = left.cmp(right);
        if order != Ordering::Equal {
            return Ok(order);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

pub(crate) fn rho_shift(model: &mut Model, terms: &Cochain, power: usize) -> Result<Value> {
    if power == 0 {
        return Ok(Value::RMCobar(terms.clone()));
    }
    let shifted = shift_support(model, terms, power)?;
    persist(model, shifted)
}

pub(crate) fn shift_support(model: &mut Model, terms: &[Term], power: usize) -> Result<Vec<Term>> {
    model.spend(terms.iter().map(Term::cells).sum())?;
    terms
        .iter()
        .map(|term| {
            Ok(Term {
                rho: checked_sum(term.rho, power)?,
                ..term.clone()
            })
        })
        .collect()
}

pub(crate) fn specialize(model: &mut Model, terms: &[Term]) -> Result<Value> {
    let words = complex_support(model, terms)?;
    model.allocate(words.len())?;
    Ok(Value::Cobar(words.into()))
}

pub(crate) fn complex_support(model: &mut Model, terms: &[Term]) -> Result<Vec<Word>> {
    let mut words: Vec<_> = terms
        .iter()
        .take_while(|term| term.rho == 0)
        .map(|term| term.word.clone())
        .collect();
    // At fixed tridegree the scalar tau power is forced by the word.
    let cells = words.iter().map(cobar::word_cells).sum::<usize>();
    model.spend(cells * (1 + words.len().max(1).ilog2() as usize))?;
    words.sort_unstable();
    Ok(words)
}

pub(crate) fn basis(model: &mut Model, grade: Grade, limit: usize) -> Result<Vec<Term>> {
    basis_window(model, grade, 0..=usize::MAX, limit)
}

pub(crate) fn basis_window(
    model: &mut Model,
    grade: Grade,
    filtration: RangeInclusive<usize>,
    limit: usize,
) -> Result<Vec<Term>> {
    let Grade { s, t, w } = grade;
    let w = i128::from(w);
    if t < w || (s == 0 && t > 0) {
        return Ok(Vec::new());
    }
    if s == 0 {
        let rho = usize::try_from(-t).map_err(|_| basis_overflow())?;
        if !filtration.contains(&rho) {
            return Ok(Vec::new());
        }
        if limit == 0 {
            return Err(cobar::basis_limit(limit));
        }
        let tau = usize::try_from(t - w).map_err(|_| basis_overflow())?;
        let unit = Term::scalar(0, rho, tau);
        model.allocate(unit.cells())?;
        return Ok(vec![unit]);
    }
    let start = t.saturating_add(*filtration.start() as i128).max(s as i128);
    let stop = t
        .checked_sub(w)
        .and_then(|bound| bound.checked_mul(2))
        .ok_or_else(basis_overflow)?
        .min(t.saturating_add(*filtration.end() as i128));
    let mut terms = Vec::new();
    let mut cells = 0;
    for top in start..=stop {
        model.spend(1)?;
        let actual_top = usize::try_from(top).map_err(|_| basis_overflow())?;
        let rho = usize::try_from(top - t).map_err(|_| basis_overflow())?;
        let required = w + rho as i128;
        let weight = i64::try_from(required).map_err(|_| basis_overflow())?;
        let slice = Grade {
            s,
            t: top,
            w: weight,
        };
        for word in motivic::basis(model, slice, limit - terms.len())? {
            let intrinsic = motivic::intrinsic_weight(&word, actual_top) as i128;
            let tau = usize::try_from(intrinsic - required).map_err(|_| basis_overflow())?;
            let term = Term { rho, tau, word };
            model.spend(term.cells())?;
            cells += term.cells();
            cobar::check_storage(cells)?;
            terms.push(term);
        }
    }
    model.spend(cells * (1 + terms.len().max(1).ilog2() as usize))?;
    model.allocate(2 * terms.len())?;
    terms.sort_unstable();
    Ok(terms)
}

fn basis_overflow() -> Error {
    Error::plain("real-motivic basis degree exceeds the supported range")
}

pub(crate) fn differential_support(model: &mut Model, terms: &[Term]) -> Result<Vec<Term>> {
    differential_below(model, terms, u128::MAX)
}

pub(crate) fn differential_below(
    model: &mut Model,
    terms: &[Term],
    cutoff: u128,
) -> Result<Vec<Term>> {
    let mut image = Support::new();
    for term in terms {
        model.spend(term.cells())?;
        if term.rho as u128 >= cutoff {
            continue;
        }
        if term.rho as u128 + 1 == cutoff {
            // At a single filtration level, specialization to rho=0
            // identifies the column with a complex-motivic differential.
            let odd = motivic::odd_exponents(&term.word);
            for word in cobar::differential_word(model, &term.word)? {
                let next_odd = motivic::odd_exponents(&word);
                let tau = checked_sum(term.tau, (odd - next_odd) / 2)?;
                toggle(
                    model,
                    &mut image,
                    Term {
                        tau,
                        word,
                        ..term.clone()
                    },
                )?;
            }
            continue;
        }
        if term.tau != 0 {
            // Transport affects only the new first slot. Reduce it
            // before adjoining the old word.
            let first = Term::scalar(1, term.rho, 0);
            for mut coface in transport(model, first, 1, term.tau, cutoff)? {
                if coface.normalized() {
                    coface.word = cobar::join(&coface.word, &term.word);
                    toggle(model, &mut image, coface)?;
                }
            }
        }
        for (slot, factor) in term.word.iter().enumerate() {
            for coface in split(model, term, slot, factor, cutoff)? {
                if coface.normalized() {
                    toggle(model, &mut image, coface)?;
                }
            }
        }
        // The last coface inserts a unit arrow and has zero normalized image.
    }
    Ok(image.into_vec())
}

fn split(
    model: &mut Model,
    term: &Term,
    slot: usize,
    factor: &Monomial,
    cutoff: u128,
) -> Result<Vec<Term>> {
    let split_word = |left, right| {
        term.word[..slot]
            .iter()
            .cloned()
            .chain([left, right])
            .chain(term.word[slot + 1..].iter().cloned())
            .collect()
    };
    // With at most one odd generator no coproduct factor can
    // contain a tau square. The ordinary Milnor expansion is exact.
    if cobar::odd_exponents(factor) <= 1 {
        let mut coface = Vec::new();
        let mut cells = 0;
        for (left, right) in cobar::coproduct(model, Profile::Full, factor)? {
            if left.is_empty() || right.is_empty() {
                continue;
            }
            let split = Term {
                word: split_word(left, right),
                ..term.clone()
            };
            model.spend(split.cells())?;
            cells += split.cells();
            cobar::check_storage(cells)?;
            coface.push(split);
        }
        return Ok(coface);
    }

    let width = term.word.len() + 1;
    let even_factor = factor.iter().map(|power| power & !1).collect();
    let mut coface = Support::new();
    for (left, right) in cobar::coproduct(model, Profile::Full, &even_factor)? {
        let word = split_word(left, right);
        toggle(
            model,
            &mut coface,
            Term {
                word,
                ..term.clone()
            },
        )?;
    }
    let mut coface = coface.into_vec();
    for (index, _) in factor
        .iter()
        .enumerate()
        .filter(|(_, power)| *power & 1 != 0)
    {
        let mut factors = Support::new();
        for (left, right) in cobar::coproduct(model, Profile::Full, &generator(index))? {
            let mut word: Vec<Monomial> = vec![[].into(); width];
            word[slot] = left;
            word[slot + 1] = right;
            toggle(
                model,
                &mut factors,
                Term {
                    rho: 0,
                    tau: 0,
                    word: word.into(),
                },
            )?;
        }
        coface = multiply(model, &coface, &factors.into_vec(), cutoff)?;
    }
    Ok(coface)
}

pub(crate) fn concat_support(
    model: &mut Model,
    left: &[Term],
    right: &[Term],
) -> Result<Vec<Term>> {
    concat_below(model, left, right, u128::MAX)
}

pub(crate) fn concat_below(
    model: &mut Model,
    left: &[Term],
    right: &[Term],
    cutoff: u128,
) -> Result<Vec<Term>> {
    let mut result = Support::new();
    for a in left {
        for b in right {
            if a.rho as u128 + b.rho as u128 >= cutoff {
                model.spend(1)?;
                break;
            }
            let cells = a.cells() + b.cells();
            cobar::check_storage(cells)?;
            model.spend(cells)?;
            let word = cobar::join(&a.word, &b.word);
            let term = Term {
                rho: checked_sum(a.rho, b.rho)?,
                tau: a.tau,
                word,
            };
            for term in transport(model, term, a.word.len(), b.tau, cutoff)? {
                toggle(model, &mut result, term)?;
            }
        }
    }
    Ok(result.into_vec())
}

// Transport a scalar through the preceding arrows via the right unit.
fn transport(
    model: &mut Model,
    mut term: Term,
    slots: usize,
    power: usize,
    cutoff: u128,
) -> Result<Vec<Term>> {
    if slots == 0 || power == 0 {
        term.tau = checked_sum(term.tau, power)?;
        return Ok(vec![term]);
    }
    let s = term.word.len();
    let mut scalar = Support::new();
    toggle(model, &mut scalar, Term::scalar(s, 0, 1))?;
    for slot in 0..slots {
        let term = Term {
            rho: 1,
            ..Term::arrow(s, slot, 0)
        };
        toggle(model, &mut scalar, term)?;
    }
    multiply_power(model, vec![term], scalar.into_vec(), power, cutoff)
}

fn multiply_power(
    model: &mut Model,
    mut result: Vec<Term>,
    mut factor: Vec<Term>,
    mut exponent: usize,
    cutoff: u128,
) -> Result<Vec<Term>> {
    let Some(first) = result.first() else {
        return Ok(result);
    };
    let factor_cutoff = cutoff.saturating_sub(first.rho as u128);
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = multiply(model, &result, &factor, cutoff)?;
            if result.is_empty() {
                break;
            }
        }
        exponent >>= 1;
        if exponent != 0 {
            let mut square = Support::new();
            // Frobenius in the commutative tensor algebra over F2.
            for term in &factor {
                multiply_terms(model, term, term, &mut square, factor_cutoff)?;
            }
            factor = square.into_vec();
        }
    }
    Ok(result)
}

fn multiply(model: &mut Model, left: &[Term], right: &[Term], cutoff: u128) -> Result<Vec<Term>> {
    let mut result = Support::new();
    for a in left {
        for b in right {
            if a.rho as u128 + b.rho as u128 >= cutoff {
                model.spend(1)?;
                break;
            }
            multiply_terms(model, a, b, &mut result, cutoff)?;
        }
    }
    Ok(result.into_vec())
}

fn toggle(model: &mut Model, support: &mut Support<Term>, term: Term) -> Result<()> {
    let cells = term.cells();
    support.toggle(model, term, cells)
}

fn trim(factor: &mut Vec<usize>) {
    if let Some(last) = factor.iter().rposition(|&power| power != 0) {
        factor.truncate(last + 1);
    } else {
        factor.clear();
    }
}

#[derive(Clone)]
struct Reduction {
    rho: usize,
    tau: usize,
    factors: Vec<Vec<usize>>,
    squares: Vec<(usize, usize)>,
}

impl Reduction {
    fn product(left: &Term, right: &Term) -> Result<Self> {
        let mut squares = Vec::new();
        let factors = left
            .word
            .iter()
            .zip(&*right.word)
            .enumerate()
            .map(|(slot, (a, b))| {
                (0..a.len().max(b.len()))
                    .map(|index| {
                        let x = a.get(index).copied().unwrap_or(0);
                        let y = b.get(index).copied().unwrap_or(0);
                        if x & y & 1 != 0 {
                            squares.push((slot, index));
                            checked_sum(x - 1, y - 1)
                        } else {
                            checked_sum(x, y)
                        }
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            rho: checked_sum(left.rho, right.rho)?,
            tau: checked_sum(left.tau, right.tau)?,
            factors,
            squares,
        })
    }

    fn exponent(&mut self, slot: usize, index: usize) -> Result<&mut usize> {
        if index >= usize::BITS as usize {
            return Err(basis_overflow());
        }
        let factor = &mut self.factors[slot];
        factor.resize(factor.len().max(index + 1), 0);
        Ok(&mut factor[index])
    }

    fn add_tau(&mut self, slot: usize, index: usize) -> Result<()> {
        let exponent = self.exponent(slot, index)?;
        if *exponent & 1 != 0 {
            *exponent -= 1;
            self.squares.push((slot, index));
        } else {
            *exponent += 1;
        }
        Ok(())
    }

    fn cells(&self) -> usize {
        3 + self
            .factors
            .iter()
            .map(|factor| 1 + factor.len())
            .sum::<usize>()
            + 2 * self.squares.len()
    }

    fn finish(self) -> Term {
        let word = self
            .factors
            .into_iter()
            .map(|mut factor| {
                trim(&mut factor);
                factor.into()
            })
            .collect();
        Term {
            rho: self.rho,
            tau: self.tau,
            word,
        }
    }
}

fn multiply_terms(
    model: &mut Model,
    left: &Term,
    right: &Term,
    output: &mut Support<Term>,
    cutoff: u128,
) -> Result<()> {
    debug_assert_eq!(left.word.len(), right.word.len());
    model.spend(left.cells() + right.cells())?;
    if left.rho as u128 + right.rho as u128 >= cutoff {
        return Ok(());
    }
    let first = Reduction::product(left, right)?;
    let mut queued = first.cells();
    cobar::check_storage(queued)?;
    let mut stack = vec![first];
    while let Some(mut term) = stack.pop() {
        let cells = term.cells();
        queued -= cells;
        model.spend(cells)?;
        let Some((slot, index)) = term.squares.pop() else {
            toggle(model, output, term.finish())?;
            continue;
        };

        // At slot j, tau_i^2 equals
        // (tau + rho * sum_{k<=j} tau_{0,k}) * xi_{i+1,j} + rho * tau_{i+1,j}.
        let mut xi = term.clone();
        let exponent = xi.exponent(slot, index)?;
        *exponent = checked_sum(*exponent, 2)?;
        let mut scalar = xi.clone();
        scalar.tau = checked_sum(scalar.tau, 1)?;
        queue(model, &mut stack, &mut queued, scalar)?;
        if term.rho as u128 + 1 < cutoff {
            xi.rho = checked_sum(xi.rho, 1)?;
            for preceding in 0..=slot {
                let mut carry = xi.clone();
                carry.add_tau(preceding, 0)?;
                queue(model, &mut stack, &mut queued, carry)?;
            }
            term.rho = checked_sum(term.rho, 1)?;
            term.add_tau(slot, index + 1)?;
            queue(model, &mut stack, &mut queued, term)?;
        }
    }
    Ok(())
}

fn queue(
    model: &mut Model,
    stack: &mut Vec<Reduction>,
    queued: &mut usize,
    term: Reduction,
) -> Result<()> {
    let cells = term.cells();
    model.spend(cells)?;
    cobar::check_storage(*queued + cells)?;
    *queued += cells;
    stack.push(term);
    Ok(())
}
