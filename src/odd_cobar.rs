use std::cmp::Ordering;
use std::collections::{btree_map::Entry, BTreeMap, HashMap};
use std::rc::Rc;

use crate::error::{Error, Result, Span};
use crate::model::{Model, Value};
use crate::odd_profile;
use crate::profile::Profile;

const MAX_WORD_CELLS: usize = 300_000;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Prime(u32);

impl Prime {
    pub(crate) fn new(model: &mut Model, p: usize) -> Result<Self> {
        let p = u32::try_from(p).map_err(|_| Error::plain("prime must fit in u32"))?;
        if p < 3 || p % 2 == 0 {
            return Err(Error::plain("PCobar requires an odd prime"));
        }
        let mut divisor = 3;
        while divisor <= p / divisor {
            model.spend(1)?;
            if p % divisor == 0 {
                return Err(Error::plain(format!("{p} is not prime")));
            }
            divisor += 2;
        }
        Ok(Self(p))
    }

    pub(crate) fn get(self) -> u32 {
        self.0
    }

    pub(crate) fn reduce(self, value: usize) -> u32 {
        (value % self.0 as usize) as u32
    }

    pub(crate) fn add(self, a: u32, b: u32) -> u32 {
        ((u64::from(a) + u64::from(b)) % u64::from(self.0)) as u32
    }

    pub(crate) fn mul(self, a: u32, b: u32) -> u32 {
        ((u64::from(a) * u64::from(b)) % u64::from(self.0)) as u32
    }

    pub(crate) fn negate(self, value: u32) -> u32 {
        if value == 0 {
            0
        } else {
            self.0 - value
        }
    }

    fn signed(self, value: u32, negative: bool) -> u32 {
        if negative {
            self.negate(value)
        } else {
            value
        }
    }
}

// Polynomial exponents are trimmed, and exterior indices are increasing.
#[derive(Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Monomial {
    pub(crate) exponents: Box<[usize]>,
    pub(crate) taus: Box<[usize]>,
    pub(crate) degree: usize,
}

pub(crate) type Factor = Rc<Monomial>;
pub(crate) type Word = Rc<[Factor]>;
pub(crate) type Term = (Word, u32);
pub(crate) type Support = Rc<[Term]>;
pub(crate) type FactorLiteral = (Vec<usize>, Vec<usize>, Span);
type Tensor = (Factor, Factor);

impl Monomial {
    pub(crate) fn cells(&self) -> usize {
        1 + self.exponents.len() + self.taus.len()
    }

    fn odd(&self) -> bool {
        self.degree % 2 == 1
    }
}

pub(crate) fn factor(mut exponents: Vec<usize>, taus: Vec<usize>, degree: usize) -> Factor {
    let width = exponents
        .iter()
        .rposition(|&exponent| exponent != 0)
        .map_or(0, |index| index + 1);
    exponents.truncate(width);
    Rc::new(Monomial {
        exponents: exponents.into_boxed_slice(),
        taus: taus.into_boxed_slice(),
        degree,
    })
}

fn unit() -> Factor {
    factor(Vec::new(), Vec::new(), 0)
}

fn degree_overflow() -> Error {
    Error::plain("odd-prime cobar internal degree is out of range")
}

pub(crate) fn power(prime: Prime, index: usize) -> Result<usize> {
    if index >= usize::BITS as usize {
        return Err(degree_overflow());
    }
    (prime.get() as usize)
        .checked_pow(index as u32)
        .ok_or_else(degree_overflow)
}

pub(crate) fn xi_degree(prime: Prime, index: usize) -> Result<usize> {
    (power(prime, index)? - 1)
        .checked_mul(2)
        .ok_or_else(degree_overflow)
}

fn tau_degree(prime: Prime, index: usize) -> Result<usize> {
    xi_degree(prime, index)?
        .checked_add(1)
        .ok_or_else(degree_overflow)
}

fn monomial_degree(prime: Prime, exponents: &[usize], taus: &[usize]) -> Result<usize> {
    let mut degree = 0usize;
    for (index, &exponent) in exponents.iter().enumerate() {
        if exponent != 0 {
            let contribution = xi_degree(prime, index + 1)?
                .checked_mul(exponent)
                .ok_or_else(degree_overflow)?;
            degree = degree
                .checked_add(contribution)
                .ok_or_else(degree_overflow)?;
        }
    }
    for &index in taus {
        degree = degree
            .checked_add(tau_degree(prime, index)?)
            .ok_or_else(degree_overflow)?;
    }
    Ok(degree)
}

fn xi(prime: Prime, index: usize, exponent: usize) -> Result<Factor> {
    if index == 0 {
        return Ok(unit());
    }
    let degree = xi_degree(prime, index)?
        .checked_mul(exponent)
        .ok_or_else(degree_overflow)?;
    let mut exponents = vec![0; index];
    exponents[index - 1] = exponent;
    Ok(factor(exponents, Vec::new(), degree))
}

fn tau(prime: Prime, index: usize) -> Result<Factor> {
    Ok(factor(Vec::new(), vec![index], tau_degree(prime, index)?))
}

pub(crate) fn word_cells(word: &Word) -> usize {
    1 + word.iter().map(|factor| factor.cells()).sum::<usize>()
}

fn word_degree(word: &Word) -> Result<usize> {
    word.iter().try_fold(0usize, |degree, factor| {
        degree
            .checked_add(factor.degree)
            .ok_or_else(degree_overflow)
    })
}

pub(crate) fn word_metadata(word: &Word) -> Result<(usize, usize)> {
    Ok((word.len(), word_degree(word)?))
}

pub(crate) fn term_cells((word, _): &Term) -> usize {
    1 + word_cells(word)
}

pub(crate) fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_WORD_CELLS {
        return Err(Error::plain(
            "odd-prime cobar calculation exceeds the word storage limit",
        ));
    }
    Ok(())
}

pub(crate) fn zero() -> Value {
    Value::PCobar([].into())
}

pub(crate) fn literal(
    model: &mut Model,
    prime: Prime,
    factors: Vec<FactorLiteral>,
) -> Result<((usize, usize), Value)> {
    check_storage(factors.len().saturating_mul(2).saturating_add(2))?;
    let mut word = Vec::with_capacity(factors.len());
    let mut cells = 2;
    for (exponents, taus, span) in factors {
        model.spend(1 + exponents.len() + taus.len())?;
        if taus.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(Error::new(span, "tau indices must be strictly increasing"));
        }
        let degree = monomial_degree(prime, &exponents, &taus).map_err(|error| error.at(span))?;
        if degree == 0 {
            return Err(Error::new(
                span,
                "cobar word factors must have positive degree",
            ));
        }
        let monomial = factor(exponents, taus, degree);
        cells += monomial.cells();
        check_storage(cells)?;
        word.push(monomial);
    }
    let word: Word = word.into();
    let degree = word_metadata(&word)?;
    model.allocate(cells)?;
    Ok((degree, Value::PCobar([(word, 1)].into())))
}

pub(crate) fn persist(model: &mut Model, terms: Vec<Term>) -> Result<Value> {
    let cells = terms.iter().map(term_cells).sum();
    check_storage(cells)?;
    model.spend(cells)?;
    model.allocate(cells)?;
    Ok(Value::PCobar(terms.into()))
}

pub(crate) fn add(model: &mut Model, prime: Prime, left: &[Term], right: &[Term]) -> Result<Value> {
    let terms = add_support(model, prime, left, right)?;
    persist(model, terms)
}

pub(crate) fn add_support(
    model: &mut Model,
    prime: Prime,
    mut left: &[Term],
    mut right: &[Term],
) -> Result<Vec<Term>> {
    model.spend(left.iter().chain(right).map(term_cells).sum())?;
    let mut terms = Vec::new();
    while let (Some((a, rest_left)), Some((b, rest_right))) =
        (left.split_first(), right.split_first())
    {
        match a.0.cmp(&b.0) {
            Ordering::Less => {
                terms.push(a.clone());
                left = rest_left;
            }
            Ordering::Greater => {
                terms.push(b.clone());
                right = rest_right;
            }
            Ordering::Equal => {
                let coefficient = prime.add(a.1, b.1);
                if coefficient != 0 {
                    terms.push((a.0.clone(), coefficient));
                }
                left = rest_left;
                right = rest_right;
            }
        }
    }
    terms.extend_from_slice(left);
    terms.extend_from_slice(right);
    check_storage(terms.iter().map(term_cells).sum())?;
    Ok(terms)
}

pub(crate) fn scale(
    model: &mut Model,
    prime: Prime,
    scalar: usize,
    input: &[Term],
) -> Result<Value> {
    let scalar = prime.reduce(scalar);
    if scalar == 0 || input.is_empty() {
        return Ok(zero());
    }
    model.spend(input.len())?;
    let terms = input
        .iter()
        .map(|(word, coefficient)| (word.clone(), prime.mul(*coefficient, scalar)))
        .collect();
    persist(model, terms)
}

pub(crate) fn compare(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Ordering> {
    // Compare the highest word first, then its field digit, as in enumeration.
    for (a, b) in left.iter().rev().zip(right.iter().rev()) {
        model.spend(term_cells(a) + term_cells(b))?;
        let order = a.cmp(b);
        if order != Ordering::Equal {
            return Ok(order);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

pub(crate) fn concat(
    model: &mut Model,
    prime: Prime,
    left: &[Term],
    right: &[Term],
) -> Result<Value> {
    let terms = concat_support(model, prime, left, right)?;
    persist(model, terms)
}

pub(crate) fn concat_support(
    model: &mut Model,
    prime: Prime,
    left: &[Term],
    right: &[Term],
) -> Result<Vec<Term>> {
    let mut terms = Vec::new();
    let mut cells = 0;
    // The common left word length makes the split unique and keeps order.
    for (a, ac) in left {
        for (b, bc) in right {
            let size = word_cells(a) + word_cells(b);
            cells += size;
            check_storage(cells)?;
            model.spend(size)?;
            word_degree(a)?
                .checked_add(word_degree(b)?)
                .ok_or_else(degree_overflow)?;
            let word = if b.is_empty() {
                a.clone()
            } else if a.is_empty() {
                b.clone()
            } else {
                a.iter().chain(b.iter()).cloned().collect()
            };
            terms.push((word, prime.mul(*ac, *bc)));
        }
    }
    Ok(terms)
}

pub(crate) fn differential(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    input: &[Term],
) -> Result<Value> {
    let terms = differential_profile_support(model, prime, profile, input)?;
    persist(model, terms)
}

pub(crate) fn differential_word(model: &mut Model, prime: Prime, word: &Word) -> Result<Vec<Term>> {
    differential_support(model, prime, &[(word.clone(), 1)])
}

pub(crate) fn differential_profile_word(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    word: &Word,
) -> Result<Vec<Term>> {
    differential_profile_support(model, prime, profile, &[(word.clone(), 1)])
}

pub(crate) fn differential_support(
    model: &mut Model,
    prime: Prime,
    input: &[Term],
) -> Result<Vec<Term>> {
    differential_profile_support(model, prime, Profile::Full, input)
}

pub(crate) fn differential_profile_support(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    input: &[Term],
) -> Result<Vec<Term>> {
    let mut result = Accumulator::new(prime);
    for (word, coefficient) in input {
        let cells = word_cells(word);
        model.spend(cells)?;
        let mut prefix_odd = false;
        for (position, current) in word.iter().enumerate() {
            for ((left, right), coproduct_coefficient) in
                reduced_coproduct(model, prime, profile, current)?
            {
                let size = 1 + cells + left.cells() + right.cells() - current.cells();
                check_storage(size)?;
                let coefficient = prime.signed(
                    prime.mul(*coefficient, coproduct_coefficient),
                    prefix_odd ^ left.odd(),
                );
                let mut output = Vec::with_capacity(word.len() + 1);
                output.extend_from_slice(&word[..position]);
                output.extend([left, right]);
                output.extend_from_slice(&word[position + 1..]);
                result.insert(model, output.into(), coefficient, size)?;
            }
            // Each desuspended even-degree factor has odd stem.
            prefix_odd ^= !current.odd();
        }
    }
    Ok(result.into_terms())
}

fn reduced_coproduct(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    monomial: &Monomial,
) -> Result<Vec<(Tensor, u32)>> {
    let mut terms = vec![((unit(), unit()), 1)];
    for (index, mut exponent) in monomial.exponents.iter().copied().enumerate() {
        let generator = index + 1;
        let p = prime.get() as usize;
        let mut frobenius = 1usize;
        while exponent != 0 {
            let digit = exponent % p;
            if digit != 0 {
                let mut column = Vec::with_capacity(generator + 1);
                for i in 0..=generator {
                    // Frobenius acts on an entirely even tensor algebra here.
                    let left_exponent = if i == generator {
                        0
                    } else {
                        frobenius
                            .checked_mul(power(prime, i)?)
                            .ok_or_else(degree_overflow)?
                    };
                    column.push((
                        xi(prime, generator - i, left_exponent)?,
                        xi(prime, i, frobenius)?,
                    ));
                }
                for _ in 0..digit {
                    terms = multiply_tensor(model, prime, profile, &terms, &column)?;
                }
            }
            exponent /= p;
            if exponent != 0 {
                frobenius = frobenius.checked_mul(p).ok_or_else(degree_overflow)?;
            }
        }
    }
    for &generator in &monomial.taus {
        let mut column = vec![(tau(prime, generator)?, unit())];
        for i in 0..=generator {
            let left_exponent = if i == generator { 0 } else { power(prime, i)? };
            column.push((xi(prime, generator - i, left_exponent)?, tau(prime, i)?));
        }
        terms = multiply_tensor(model, prime, profile, &terms, &column)?;
    }
    terms.retain(|((left, right), _)| left.degree != 0 && right.degree != 0);
    Ok(terms)
}

fn multiply_tensor(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    terms: &[(Tensor, u32)],
    column: &[Tensor],
) -> Result<Vec<(Tensor, u32)>> {
    let mut result = Accumulator::new(prime);
    for ((a, b), coefficient) in terms {
        for (c, d) in column {
            model.spend(a.cells() + b.cells() + c.cells() + d.cells())?;
            let Some((left, left_sign)) = multiply_monomials(prime, profile, a, c)? else {
                continue;
            };
            let Some((right, right_sign)) = multiply_monomials(prime, profile, b, d)? else {
                continue;
            };
            let size = 1 + left.cells() + right.cells();
            let negative = (b.odd() && c.odd()) ^ left_sign ^ right_sign;
            result.insert(
                model,
                (left, right),
                prime.signed(*coefficient, negative),
                size,
            )?;
        }
    }
    Ok(result.into_terms())
}

fn multiply_monomials(
    prime: Prime,
    profile: Profile,
    a: &Factor,
    b: &Factor,
) -> Result<Option<(Factor, bool)>> {
    if a.degree == 0 {
        return Ok(Some((b.clone(), false)));
    }
    if b.degree == 0 {
        return Ok(Some((a.clone(), false)));
    }
    let mut left = a.taus.as_ref();
    let mut right = b.taus.as_ref();
    let mut taus = Vec::with_capacity(left.len() + right.len());
    let mut odd = false;
    while let (Some((&i, rest_left)), Some((&j, rest_right))) =
        (left.split_first(), right.split_first())
    {
        match i.cmp(&j) {
            Ordering::Less => {
                taus.push(i);
                left = rest_left;
            }
            Ordering::Equal => return Ok(None),
            Ordering::Greater => {
                odd ^= left.len() % 2 != 0;
                taus.push(j);
                right = rest_right;
            }
        }
    }
    taus.extend_from_slice(left);
    taus.extend_from_slice(right);
    let mut exponents = a.exponents.to_vec();
    exponents.resize(exponents.len().max(b.exponents.len()), 0);
    for (index, (exponent, &other)) in exponents.iter_mut().zip(&b.exponents).enumerate() {
        *exponent = exponent.checked_add(other).ok_or_else(degree_overflow)?;
        // Further products cannot reduce an exponent that leaves the profile.
        if !odd_profile::allows_power(prime, profile, index, *exponent) {
            return Ok(None);
        }
    }
    let degree = a.degree.checked_add(b.degree).ok_or_else(degree_overflow)?;
    Ok(Some((factor(exponents, taus, degree), odd)))
}

struct Accumulator<K> {
    prime: Prime,
    terms: BTreeMap<K, u32>,
    cells: usize,
}

impl<K: Ord> Accumulator<K> {
    fn new(prime: Prime) -> Self {
        Self {
            prime,
            terms: BTreeMap::new(),
            cells: 0,
        }
    }

    fn insert(&mut self, model: &mut Model, key: K, coefficient: u32, cells: usize) -> Result<()> {
        model.spend(cells * (1 + self.terms.len().max(1).ilog2() as usize))?;
        if coefficient == 0 {
            return Ok(());
        }
        match self.terms.entry(key) {
            Entry::Occupied(mut entry) => {
                let coefficient = self.prime.add(*entry.get(), coefficient);
                if coefficient == 0 {
                    entry.remove();
                    self.cells -= cells;
                } else {
                    *entry.get_mut() = coefficient;
                }
            }
            Entry::Vacant(entry) => {
                check_storage(self.cells + cells)?;
                self.cells += cells;
                entry.insert(coefficient);
            }
        }
        Ok(())
    }

    fn into_terms(self) -> Vec<(K, u32)> {
        self.terms.into_iter().collect()
    }
}

pub(crate) fn basis_profile(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    s: usize,
    t: usize,
    limit: usize,
) -> Result<Vec<Word>> {
    if let Profile::An(n) = profile {
        if s != 0 && t >= s && odd_profile::changes_basis(prime, n, t - s + 1) {
            return odd_profile::basis(model, prime, n, s, t, limit);
        }
    }
    basis(model, prime, s, t, limit)
}

pub(crate) fn basis(
    model: &mut Model,
    prime: Prime,
    s: usize,
    t: usize,
    limit: usize,
) -> Result<Vec<Word>> {
    model.spend(1)?;
    if s == 0 {
        return if t == 0 {
            if limit == 0 {
                return Err(basis_limit(limit));
            }
            model.allocate(1)?;
            Ok(vec![[].into()])
        } else {
            Ok(Vec::new())
        };
    }
    let Some(excess) = t.checked_sub(s) else {
        return Ok(Vec::new());
    };
    // Tau_0 has degree one. The next permitted degree is |xi_1|.
    if excess != 0 && (excess as u128) < 2 * u128::from(prime.get() - 1) - 1 {
        return Ok(Vec::new());
    }
    check_storage(s.saturating_mul(2).saturating_add(1))?;

    let mut degrees = vec![1; s];
    degrees[s - 1] += excess;
    let mut catalogue = HashMap::new();
    let mut catalogue_cells = 0;
    let mut words = Vec::new();
    let mut cells = 0;
    loop {
        model.spend(s)?;
        for &degree in &degrees {
            if let std::collections::hash_map::Entry::Vacant(entry) = catalogue.entry(degree) {
                catalogue_cells += 1;
                check_storage(cells + catalogue_cells)?;
                entry.insert(monomials(model, prime, degree, &mut catalogue_cells)?);
            }
            if catalogue[&degree].is_empty() {
                break;
            }
        }
        let factors = degrees
            .iter()
            .map(|degree| {
                catalogue
                    .get(degree)
                    .filter(|monomials| !monomials.is_empty())
            })
            .collect::<Option<Vec<_>>>();
        if let Some(factors) = factors {
            let mut digits = vec![0; s];
            loop {
                if words.len() == limit {
                    return Err(basis_limit(limit));
                }
                let word: Word = factors
                    .iter()
                    .zip(&digits)
                    .map(|(choices, &digit)| choices[digit].clone())
                    .collect();
                let size = word_cells(&word);
                cells += size;
                check_storage(cells + catalogue_cells)?;
                model.spend(size)?;
                words.push(word);
                if !next_choices(&factors, &mut digits) {
                    break;
                }
            }
        }
        if !next_degrees(&mut degrees) {
            break;
        }
    }
    model.spend(cells * (1 + words.len().max(1).ilog2() as usize))?;
    model.allocate(cells)?;
    words.sort_unstable();
    Ok(words)
}

pub(crate) fn basis_limit(limit: usize) -> Error {
    Error::plain(format!(
        "odd-prime cobar basis exceeds the limit of {limit} words"
    ))
}

fn next_degrees(degrees: &mut [usize]) -> bool {
    let Some(index) = (1..degrees.len()).rev().find(|&index| degrees[index] > 1) else {
        return false;
    };
    let excess = degrees[index] - 1;
    degrees[index - 1] += 1;
    degrees[index] = 1;
    *degrees.last_mut().unwrap() = excess;
    true
}

fn next_choices(factors: &[&Vec<Factor>], digits: &mut [usize]) -> bool {
    for (choices, digit) in factors.iter().zip(digits).rev() {
        *digit += 1;
        if *digit < choices.len() {
            return true;
        }
        *digit = 0;
    }
    false
}

fn monomials(
    model: &mut Model,
    prime: Prime,
    degree: usize,
    cells: &mut usize,
) -> Result<Vec<Factor>> {
    let xis: Vec<_> = (1..)
        .map_while(|index| {
            xi_degree(prime, index)
                .ok()
                .filter(|&weight| weight <= degree)
        })
        .collect();
    let taus: Vec<_> = (0..)
        .map_while(|index| {
            tau_degree(prime, index)
                .ok()
                .filter(|&weight| weight <= degree)
        })
        .collect();
    model.spend(xis.len() + taus.len())?;
    let mut catalogue = Monomials {
        degree,
        xis,
        taus,
        terms: Vec::new(),
        cells,
    };
    catalogue.exterior(model, degree, 0, &mut Vec::new())?;
    Ok(catalogue.terms)
}

struct Monomials<'a> {
    degree: usize,
    xis: Vec<usize>,
    taus: Vec<usize>,
    terms: Vec<Factor>,
    cells: &'a mut usize,
}

impl Monomials<'_> {
    fn exterior(
        &mut self,
        model: &mut Model,
        remaining: usize,
        start: usize,
        taus: &mut Vec<usize>,
    ) -> Result<()> {
        model.spend(1)?;
        self.polynomial(model, remaining, taus)?;
        // The recursion depth is the number of generators whose degrees fit
        // a machine integer, even when the degree has many exterior subsets.
        for index in start..self.taus.len() {
            let Some(rest) = remaining.checked_sub(self.taus[index]) else {
                break;
            };
            taus.push(index);
            self.exterior(model, rest, index + 1, taus)?;
            taus.pop();
        }
        Ok(())
    }

    fn polynomial(&mut self, model: &mut Model, remaining: usize, taus: &[usize]) -> Result<()> {
        if remaining == 0 {
            return self.push(model, Vec::new(), taus);
        }
        let Some(&base) = self.xis.first() else {
            return Ok(());
        };
        if !remaining.is_multiple_of(base) {
            return Ok(());
        }
        let mut digits = vec![0; self.xis.len()];
        digits[0] = remaining / base;
        loop {
            self.push(model, digits.clone(), taus)?;
            let mut index = 1;
            while index < digits.len() {
                model.spend(1)?;
                let price = self.xis[index] / base;
                if digits[0] >= price {
                    digits[0] -= price;
                    digits[index] += 1;
                    break;
                }
                // Every polynomial degree is a multiple of |xi_1|.
                // Its exponent pays for each increment and absorbs carries.
                digits[0] += digits[index] * price;
                digits[index] = 0;
                index += 1;
            }
            if index == digits.len() {
                return Ok(());
            }
        }
    }

    fn push(&mut self, model: &mut Model, exponents: Vec<usize>, taus: &[usize]) -> Result<()> {
        model.spend(1 + exponents.len() + taus.len())?;
        let term = factor(exponents, taus.to_vec(), self.degree);
        *self.cells += term.cells();
        check_storage(*self.cells)?;
        self.terms.push(term);
        Ok(())
    }
}

pub(crate) fn elements(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    s: usize,
    t: usize,
    limit: usize,
) -> Result<Vec<Value>> {
    if limit == 0 {
        return Err(Error::plain(
            "odd-prime cochain space exceeds the enumeration limit",
        ));
    }
    let p = prime.get() as usize;
    let mut size = 1;
    let mut rank_limit = 0;
    while size <= limit / p {
        size *= p;
        rank_limit += 1;
    }
    let words = basis_profile(model, prime, profile, s, t, rank_limit)?;
    let size = (0..words.len()).fold(1usize, |size, _| size * p);
    model.allocate(size)?;
    model.spend(size.saturating_mul(1 + words.len()))?;
    let mut values = Vec::with_capacity(size);
    for mut index in 0..size {
        let mut terms = Vec::new();
        for word in &words {
            let coefficient = (index % p) as u32;
            if coefficient != 0 {
                terms.push((word.clone(), coefficient));
            }
            index /= p;
        }
        model.allocate(2 * terms.len())?;
        values.push(Value::PCobar(terms.into()));
    }
    Ok(values)
}
