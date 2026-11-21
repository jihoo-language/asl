use std::cell::OnceCell;
use std::collections::BTreeSet;
use std::rc::Rc;

use crate::complex::xor_sorted;
use crate::error::{Error, Result};
use crate::linear::{Basis, QuotientBasis};
use crate::model::Model;
use crate::rv_ideal::Ideal;
use crate::rv_poly::{self, Monomial, Polynomial};

const MAX_BASIS: usize = 32_768;
const MAX_CELLS: usize = 300_000;

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) ideal: Rc<str>,
    pub(crate) s: usize,
    pub(crate) t: usize,
}

impl Grade {
    pub(crate) fn differential(&self) -> Result<Self> {
        let s = self
            .s
            .checked_sub(1)
            .ok_or_else(|| Error::plain("rvkD requires positive Koszul homological degree"))?;
        Ok(Self { s, ..self.clone() })
    }

    pub(crate) fn primitive(&self) -> Result<Self> {
        let s = self.s.checked_add(1).ok_or_else(grade_limit)?;
        Ok(Self { s, ..self.clone() })
    }

    pub(crate) fn product(&self, other: &Self) -> Result<Self> {
        if self.ideal != other.ideal {
            return Err(Error::plain(
                "Koszul multiplication requires the same rank ideal",
            ));
        }
        let s = self.s.checked_add(other.s).ok_or_else(grade_limit)?;
        let t = self.t.checked_add(other.t).ok_or_else(grade_limit)?;
        Ok(Self {
            s,
            t,
            ..self.clone()
        })
    }
}

fn grade_limit() -> Error {
    Error::plain("rank Koszul degree is out of range")
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Term {
    pub(crate) wedge: usize,
    pub(crate) powers: Monomial,
}

pub(crate) type Support = Rc<[Term]>;

pub(crate) fn mask(ideal: &Ideal, primitives: &[usize]) -> Result<usize> {
    if primitives.windows(2).any(|pair| pair[0] >= pair[1]) {
        return Err(Error::plain(
            "Koszul wedge indices must be strictly increasing",
        ));
    }
    let mut mask = 0;
    for primitive in primitives {
        let index = ideal.selected.binary_search(primitive).map_err(|_| {
            Error::plain(format!(
                "primitive {primitive} is outside the rank ideal's selection"
            ))
        })?;
        mask |= 1usize << index;
    }
    Ok(mask)
}

pub(crate) fn term(
    model: &mut Model,
    ideal: &Ideal,
    wedge: usize,
    coefficient: Polynomial,
) -> Result<Support> {
    let normal = ideal.normal(model, coefficient)?;
    model.spend(rv_poly::cells(&normal))?;
    persist(
        model,
        normal
            .iter()
            .cloned()
            .map(|powers| Term { wedge, powers })
            .collect(),
    )
}

pub(crate) fn cells(support: &[Term]) -> usize {
    1 + support
        .first()
        .map_or(0, |term| (2 + term.powers.len()) * support.len())
}

pub(crate) fn persist(model: &mut Model, terms: Vec<Term>) -> Result<Support> {
    model.allocate(cells(&terms))?;
    Ok(terms.into())
}

pub(crate) fn compare(
    model: &mut Model,
    left: &[Term],
    right: &[Term],
) -> Result<std::cmp::Ordering> {
    model.spend(cells(left) + cells(right))?;
    Ok(left.iter().rev().cmp(right.iter().rev()))
}

pub(crate) fn add(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Support> {
    let result = sum_terms(model, left, right)?;
    persist(model, result)
}

pub(crate) fn sum_terms(model: &mut Model, left: &[Term], right: &[Term]) -> Result<Vec<Term>> {
    model.spend(cells(left) + cells(right))?;
    Ok(xor_sorted(left, right))
}

pub(crate) fn differential(model: &mut Model, ideal: &Ideal, support: &[Term]) -> Result<Support> {
    let terms = differential_terms(model, ideal, support)?;
    persist(model, terms)
}

pub(crate) fn differential_terms(
    model: &mut Model,
    ideal: &Ideal,
    support: &[Term],
) -> Result<Vec<Term>> {
    model.spend(support.len())?;
    let mut terms = BTreeSet::new();
    for term in support {
        let mut remaining = term.wedge;
        while remaining != 0 {
            let index = remaining.trailing_zeros() as usize;
            let bit = 1usize << index;
            remaining ^= bit;
            let mut powers = term.powers.to_vec();
            powers[index] = powers[index].checked_add(1).ok_or_else(grade_limit)?;
            toggle(
                model,
                &mut terms,
                Term {
                    wedge: term.wedge ^ bit,
                    powers: powers.into(),
                },
            )?;
        }
    }
    normal_terms(model, ideal, terms)
}

pub(crate) fn product(
    model: &mut Model,
    ideal: &Ideal,
    left: &[Term],
    right: &[Term],
) -> Result<Support> {
    let terms = product_terms(model, ideal, left, right)?;
    persist(model, terms)
}

pub(crate) fn product_terms(
    model: &mut Model,
    ideal: &Ideal,
    left: &[Term],
    right: &[Term],
) -> Result<Vec<Term>> {
    let mut terms = BTreeSet::new();
    for a in left {
        for b in right {
            model.spend(1)?;
            if a.wedge & b.wedge == 0 {
                toggle(
                    model,
                    &mut terms,
                    Term {
                        wedge: a.wedge | b.wedge,
                        powers: rv_poly::product_term(&a.powers, &b.powers)?,
                    },
                )?;
            }
        }
    }
    normal_terms(model, ideal, terms)
}

fn toggle(model: &mut Model, terms: &mut BTreeSet<Term>, term: Term) -> Result<()> {
    model.spend((2 + term.powers.len()) * (1 + terms.len().max(1).ilog2() as usize))?;
    if !terms.remove(&term) {
        check_storage((2 + term.powers.len()) * (terms.len() + 1))?;
        terms.insert(term);
    }
    Ok(())
}

fn normal_terms(model: &mut Model, ideal: &Ideal, terms: BTreeSet<Term>) -> Result<Vec<Term>> {
    let mut result = Vec::new();
    let mut terms = terms.into_iter().peekable();
    while let Some(Term { wedge, powers }) = terms.next() {
        let mut group = vec![powers];
        while terms.peek().is_some_and(|term| term.wedge == wedge) {
            group.push(terms.next().unwrap().powers);
        }
        model.spend(group.len())?;
        let normal = ideal.normal(model, group.into())?;
        check_storage((2 + ideal.selected.len()) * (result.len() + normal.len()))?;
        result.extend(normal.iter().cloned().map(|powers| Term { wedge, powers }));
    }
    Ok(result)
}

pub(crate) fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain(
            "rank Koszul calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

fn basis_limit() -> Error {
    Error::plain("rank Koszul basis exceeds the enumeration limit")
}

#[derive(Default)]
pub(crate) struct Slice {
    wedges: Vec<usize>,
    monomials: Vec<Monomial>,
}

impl Slice {
    pub(crate) fn build(
        model: &mut Model,
        ideal: &Ideal,
        s: usize,
        t: usize,
        limit: usize,
    ) -> Result<Self> {
        let width = ideal.selected.len();
        if s > width || s > t {
            return Ok(Self::default());
        }
        let wedge_count = choose_bounded(width, s, limit);
        let monomials = standard_monomials(model, ideal, t - s, limit / wedge_count)?;
        if monomials.is_empty() {
            return Ok(Self::default());
        }
        check_storage(wedge_count + rv_poly::cells(&monomials))?;
        let mut wedges = Vec::with_capacity(wedge_count);
        wedge_masks(model, width, s, 0, &mut wedges)?;
        Ok(Self { wedges, monomials })
    }

    pub(crate) fn len(&self) -> usize {
        self.wedges.len() * self.monomials.len()
    }

    pub(crate) fn cells(&self) -> usize {
        self.wedges.len() + rv_poly::cells(&self.monomials)
    }

    pub(crate) fn term(&self, index: usize) -> Term {
        let width = self.monomials.len();
        Term {
            wedge: self.wedges[index / width],
            powers: self.monomials[index % width].clone(),
        }
    }

    fn encode(&self, model: &mut Model, support: &[Term]) -> Result<Vec<usize>> {
        let width = self.monomials.len();
        model.spend(
            support.len()
                * (1 + self.wedges.len().max(1).ilog2() as usize
                    + self
                        .monomials
                        .first()
                        .map_or(0, |p| (1 + p.len()) * (1 + width.ilog2() as usize))),
        )?;
        support
            .iter()
            .map(|term| {
                let wedge = self.wedges.binary_search(&term.wedge).ok();
                let coefficient = self.monomials.binary_search(&term.powers).ok();
                match (wedge, coefficient) {
                    (Some(wedge), Some(coefficient)) => Ok(wedge * width + coefficient),
                    _ => Err(Error::plain(
                        "rank Koszul term is outside its standard basis",
                    )),
                }
            })
            .collect()
    }
}

fn choose_bounded(n: usize, k: usize, limit: usize) -> usize {
    let k = k.min(n - k);
    let mut value = 1u128;
    for j in 1..=k {
        value = value * (n - k + j) as u128 / j as u128;
        if value > limit as u128 {
            return limit + 1;
        }
    }
    value as usize
}

fn wedge_masks(
    model: &mut Model,
    upper: usize,
    remaining: usize,
    mask: usize,
    out: &mut Vec<usize>,
) -> Result<()> {
    model.spend(1)?;
    if remaining == 0 {
        out.push(mask);
    } else {
        for bit in remaining - 1..upper {
            wedge_masks(model, bit, remaining - 1, mask | (1usize << bit), out)?;
        }
    }
    Ok(())
}

fn standard_monomials(
    model: &mut Model,
    ideal: &Ideal,
    degree: usize,
    limit: usize,
) -> Result<Vec<Monomial>> {
    let leading: Vec<_> = ideal
        .polynomials()
        .map(|p| p.last().unwrap().as_ref())
        .collect();
    model.spend(1 + leading.len())?;
    let mut powers = vec![0; ideal.selected.len()];
    let mut output = Vec::new();
    standard_search(model, &mut powers, 0, degree, &leading, &mut output, limit)?;
    Ok(output)
}

fn standard_search(
    model: &mut Model,
    powers: &mut [usize],
    position: usize,
    degree: usize,
    leading: &[&[usize]],
    output: &mut Vec<Monomial>,
    limit: usize,
) -> Result<()> {
    model.spend(1)?;
    let free = powers.len() - position;
    if free <= 1 || degree == 0 {
        model.spend(1 + powers.len() + leading.len() * free)?;
        let tail = &mut powers[position..];
        tail.fill(0);
        if let Some(first) = tail.first_mut() {
            *first = degree;
        } else if degree != 0 {
            return Ok(());
        }
        if leading
            .iter()
            .any(|row| row[position..].iter().zip(tail.iter()).all(|(a, b)| a <= b))
        {
            return Ok(());
        }
        if output.len() == limit {
            return Err(basis_limit());
        }
        check_storage((2 + powers.len()) * (output.len() + 1))?;
        output.push(Rc::from(&*powers));
        return Ok(());
    }

    // Over an interval between consecutive leading exponents, exactly the
    // same suffix monomials are forbidden. Pure suffix powers bound the total
    // remaining degree, so a large interval can be skipped at once.
    model.spend((1 + leading.len()) * (1 + leading.len().max(1).ilog2() as usize))?;
    let mut ordered = leading.to_vec();
    ordered.sort_unstable_by_key(|row| row[position]);
    let mut pure_bounds = vec![usize::MAX; free - 1];
    let mut end = 0;
    let mut lower = 0;
    loop {
        while let Some(row) = ordered.get(end).filter(|row| row[position] <= lower) {
            let tail = &row[position + 1..];
            model.spend(1 + tail.len())?;
            let mut support = tail
                .iter()
                .enumerate()
                .filter(|(_, exponent)| **exponent != 0);
            let Some((index, &exponent)) = support.next() else {
                return Ok(());
            };
            if support.next().is_none() {
                pure_bounds[index] = pure_bounds[index].min(exponent - 1);
            }
            end += 1;
        }
        let upper = ordered
            .get(end)
            .map_or(degree, |row| degree.min(row[position] - 1));
        model.spend(free)?;
        let capacity = pure_bounds.iter().copied().fold(0, usize::saturating_add);
        for exponent in lower.max(degree.saturating_sub(capacity))..=upper {
            powers[position] = exponent;
            standard_search(
                model,
                powers,
                position + 1,
                degree - exponent,
                &ordered[..end],
                output,
                limit,
            )?;
        }
        if upper == degree {
            return Ok(());
        }
        lower = upper + 1;
    }
}

#[derive(Default)]
pub(crate) struct Tor {
    slice: Slice,
    quotient: QuotientBasis,
    preimages: OnceCell<BoundaryPreimages>,
}

struct BoundaryPreimages {
    source: Slice,
    images: Basis,
}

fn boundary_basis(
    model: &mut Model,
    ideal: &Ideal,
    source: &Slice,
    target: &Slice,
    track: bool,
) -> Result<Basis> {
    let mut basis = Basis::default();
    for index in 0..source.len() {
        let support = differential_terms(model, ideal, &[source.term(index)])?;
        let column = target.encode(model, &support)?;
        let witness = if track { vec![index] } else { Vec::new() };
        basis.insert_with_witness(model, column, witness)?;
    }
    Ok(basis)
}

impl Tor {
    pub(crate) fn build(model: &mut Model, ideal: &Ideal, grade: &Grade) -> Result<Self> {
        let Grade { s, t, .. } = *grade;
        if (s != 0 || t != 0) && ideal.polynomials().next().is_none() {
            // The Koszul resolution of the ground field is exact here.
            return Ok(Self::default());
        }
        let slice = Slice::build(model, ideal, s, t, MAX_BASIS)?;
        if slice.len() == 0 {
            return Ok(Self::default());
        }
        let incoming = Slice::build(model, ideal, s + 1, t, MAX_BASIS)?;
        let outgoing = if s == 0 {
            Slice::default()
        } else {
            Slice::build(model, ideal, s - 1, t, MAX_BASIS)?
        };
        check_storage(slice.cells() + incoming.cells() + outgoing.cells())?;

        let boundaries = boundary_basis(model, ideal, &incoming, &slice, false)?;
        let quotient = QuotientBasis::build(model, boundaries, slice.len(), |model, index| {
            let support = differential_terms(model, ideal, &[slice.term(index)])?;
            outgoing.encode(model, &support)
        })?;
        model.allocate(slice.cells())?;
        Ok(Self {
            slice,
            quotient,
            preimages: OnceCell::new(),
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.quotient.rank()
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        ideal: &Ideal,
        support: &[Term],
    ) -> Result<Vec<usize>> {
        if self.slice.len() == 0 {
            if differential_terms(model, ideal, support)?.is_empty() {
                return Ok(Vec::new());
            }
        } else {
            let coordinates = self.slice.encode(model, support)?;
            if let Some(class) = self.quotient.project(model, coordinates)? {
                return Ok(class);
            }
        }
        Err(Error::plain("rvTorClass requires a Koszul cycle"))
    }

    pub(crate) fn is_boundary(
        &self,
        model: &mut Model,
        ideal: &Ideal,
        support: &[Term],
    ) -> Result<bool> {
        if support.is_empty() {
            return Ok(true);
        }
        if self.slice.len() == 0 {
            // The only nonempty input in this case is in the acyclic part of
            // the Koszul complex over the polynomial ring itself.
            return Ok(differential_terms(model, ideal, support)?.is_empty());
        }
        let coordinates = self.slice.encode(model, support)?;
        self.quotient.is_boundary(model, coordinates)
    }

    pub(crate) fn primitive(
        &self,
        model: &mut Model,
        ideal: &Ideal,
        grade: &Grade,
        support: &[Term],
    ) -> Result<Option<Vec<Term>>> {
        if support.is_empty() {
            return Ok(Some(Vec::new()));
        }
        if self.slice.len() == 0 {
            if !self.is_boundary(model, ideal, support)? {
                return Ok(None);
            }
            return Ok(Some(polynomial_primitive(model, support)?));
        }
        let coordinates = self.slice.encode(model, support)?;
        let preimages = if let Some(preimages) = self.preimages.get() {
            preimages
        } else {
            if !self.quotient.is_boundary(model, coordinates.clone())? {
                return Ok(None);
            }
            let source_grade = grade.primitive()?;
            let source = Slice::build(model, ideal, source_grade.s, source_grade.t, MAX_BASIS)?;
            let images = boundary_basis(model, ideal, &source, &self.slice, true)?;
            model.allocate(source.cells() + images.cells())?;
            self.preimages
                .get_or_init(|| BoundaryPreimages { source, images })
        };
        let (remainder, witness) =
            preimages
                .images
                .reduce_with_witness(model, coordinates, Vec::new())?;
        Ok(remainder.is_empty().then(|| {
            witness
                .into_iter()
                .map(|index| preimages.source.term(index))
                .collect()
        }))
    }

    pub(crate) fn representative_terms(
        &self,
        model: &mut Model,
        class: &[usize],
    ) -> Result<Vec<Term>> {
        let coordinates = self.quotient.representative(model, class)?;
        Ok(coordinates
            .into_iter()
            .map(|index| self.slice.term(index))
            .collect())
    }

    pub(crate) fn representative(&self, model: &mut Model, class: &[usize]) -> Result<Support> {
        let support = self.representative_terms(model, class)?;
        persist(model, support)
    }
}

fn polynomial_primitive(model: &mut Model, support: &[Term]) -> Result<Vec<Term>> {
    // Each multigraded summand of positive weight is contracted using its
    // first active variable: divide by x_i and insert e_i if it is absent.
    let mut preimage = BTreeSet::new();
    for term in support {
        model.spend(1 + term.powers.len())?;
        for (index, &power) in term.powers.iter().enumerate() {
            let bit = 1usize << index;
            if term.wedge & bit != 0 {
                break;
            }
            if power != 0 {
                let mut powers = term.powers.to_vec();
                powers[index] -= 1;
                toggle(
                    model,
                    &mut preimage,
                    Term {
                        wedge: term.wedge | bit,
                        powers: powers.into(),
                    },
                )?;
                break;
            }
        }
    }
    Ok(preimage.into_iter().collect())
}
