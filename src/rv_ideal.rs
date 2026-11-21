use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::q_koszul::Indices;
use crate::rv_poly::{self, Monomial, Polynomial};

const MAX_CELLS: usize = 300_000;

pub(crate) struct Ideal {
    pub(crate) selected: Indices,
    degrees: Vec<usize>,
    basis: Vec<Row>,
}

struct Row {
    // Elimination rows use pure lex without an ordinary homogeneous degree.
    degree: Option<usize>,
    polynomial: Polynomial,
    certificate: Certificate,
}

#[derive(Clone, Default)]
struct Certificate {
    coefficients: BTreeMap<usize, Polynomial>,
    cells: usize,
}

impl Certificate {
    fn generator(index: usize, unit: &Polynomial) -> Self {
        Self {
            coefficients: BTreeMap::from([(index, unit.clone())]),
            cells: 1 + rv_poly::cells(unit),
        }
    }

    fn add_shifted(&mut self, model: &mut Model, other: &Self, powers: &[usize]) -> Result<()> {
        for (&index, coefficient) in &other.coefficients {
            model.spend(1 + self.coefficients.len().max(1).ilog2() as usize)?;
            let shifted = shift(model, coefficient, powers)?;
            let value = if let Some(old) = self.coefficients.remove(&index) {
                self.cells -= 1 + rv_poly::cells(&old);
                rv_poly::add(model, &old, &shifted)?
            } else {
                shifted
            };
            if !value.is_empty() {
                self.cells += 1 + rv_poly::cells(&value);
                check_storage(self.cells)?;
                self.coefficients.insert(index, value);
            }
        }
        Ok(())
    }
}

impl Row {
    fn leading(&self) -> &Monomial {
        self.polynomial.last().expect("nonzero basis row")
    }

    fn cells(&self) -> usize {
        2 + rv_poly::cells(&self.polynomial) + self.certificate.cells
    }
}

impl Ideal {
    pub(crate) fn build(
        model: &mut Model,
        selected: Indices,
        generators: Vec<(usize, Polynomial)>,
    ) -> Result<Self> {
        model.spend(1 + generators.len())?;
        check_storage(1 + generators.len())?;
        let degrees: Vec<_> = generators.iter().map(|&(degree, _)| degree).collect();
        let unit = rv_poly::monomial(vec![0; selected.len()]);
        let basis = if let Some(index) = generators
            .iter()
            .position(|(degree, p)| *degree == 0 && !p.is_empty())
        {
            vec![Row {
                degree: Some(0),
                polynomial: unit.clone(),
                certificate: Certificate::generator(index, &unit),
            }]
        } else {
            let rows = generators
                .into_iter()
                .enumerate()
                .map(|(index, (degree, polynomial))| Row {
                    degree: Some(degree),
                    polynomial,
                    certificate: Certificate::generator(index, &unit),
                });
            groebner(model, rows)?
        };
        Self::from_basis(model, selected, degrees, basis)
    }

    pub(crate) fn derived(
        model: &mut Model,
        selected: Indices,
        generators: Vec<Polynomial>,
    ) -> Result<Self> {
        let mut rows = Vec::with_capacity(generators.len());
        for polynomial in generators {
            if let Some(leading) = polynomial.last() {
                let degree = rv_poly::degree(leading)?;
                model.spend(rv_poly::cells(&polynomial))?;
                for term in polynomial.iter() {
                    if rv_poly::degree(term)? != degree {
                        return Err(Error::plain("derived ideal generator is not homogeneous"));
                    }
                }
                rows.push(Row {
                    degree: Some(degree),
                    polynomial,
                    certificate: Certificate::default(),
                });
            }
        }
        let mut basis = groebner(model, rows)?;
        let unit = rv_poly::monomial(vec![0; selected.len()]);
        let cells = basis.iter().map(Row::cells).sum::<usize>();
        check_storage(cells + basis.len() * (2 + rv_poly::cells(&unit)))?;
        model.spend(1 + basis.len())?;
        let mut degrees = Vec::with_capacity(basis.len());
        // The generating list of a derived ideal is its reduced basis.
        for (index, row) in basis.iter_mut().enumerate() {
            degrees.push(row.degree.expect("homogeneous basis row"));
            row.certificate = Certificate::generator(index, &unit);
        }
        Self::from_basis(model, selected, degrees, basis)
    }

    fn from_basis(
        model: &mut Model,
        selected: Indices,
        degrees: Vec<usize>,
        basis: Vec<Row>,
    ) -> Result<Self> {
        model.allocate(1 + degrees.len() + basis.iter().map(Row::cells).sum::<usize>())?;
        Ok(Self {
            selected,
            degrees,
            basis,
        })
    }

    pub(crate) fn polynomials(&self) -> impl Iterator<Item = &Polynomial> {
        self.basis.iter().map(|row| &row.polynomial)
    }

    pub(crate) fn generator(&self, index: usize) -> Option<(usize, Polynomial)> {
        self.basis.get(index).map(|row| {
            (
                row.degree.expect("homogeneous basis row"),
                row.polynomial.clone(),
            )
        })
    }

    pub(crate) fn normal(&self, model: &mut Model, polynomial: Polynomial) -> Result<Polynomial> {
        reduce(model, polynomial, &self.basis, None, None)
    }

    pub(crate) fn witness_degree(&self, index: usize, degree: usize) -> Result<usize> {
        let &generator_degree = self.degrees.get(index).ok_or_else(|| {
            Error::plain(format!("rvWitness generator index {index} is out of range"))
        })?;
        degree.checked_sub(generator_degree).ok_or_else(|| {
            Error::plain(
                "rvWitness cannot use a generator of degree greater than the polynomial degree",
            )
        })
    }

    pub(crate) fn witness(
        &self,
        model: &mut Model,
        index: usize,
        polynomial: Polynomial,
    ) -> Result<Polynomial> {
        let mut certificate = Certificate::default();
        let remainder = reduce(model, polynomial, &self.basis, None, Some(&mut certificate))?;
        if !remainder.is_empty() {
            return Err(Error::plain("rvWitness requires a polynomial in the ideal"));
        }
        Ok(certificate
            .coefficients
            .remove(&index)
            .unwrap_or_else(|| [].into()))
    }

    pub(crate) fn hilbert(&self, model: &mut Model, degree: usize) -> Result<String> {
        hilbert(model, self.selected.len(), degree, &self.basis)
    }
}

pub(crate) fn eliminate(model: &mut Model, generators: Vec<Polynomial>) -> Result<Vec<Polynomial>> {
    eliminate_prefix(model, generators, 1)
}

// Pure lex eliminates any initial block of variables. A row with zero
// leading exponents on the block has no term supported there.
pub(crate) fn eliminate_prefix(
    model: &mut Model,
    generators: Vec<Polynomial>,
    prefix: usize,
) -> Result<Vec<Polynomial>> {
    let rows = generators.into_iter().map(|polynomial| Row {
        degree: None,
        polynomial,
        certificate: Certificate::default(),
    });
    let mut result = Vec::new();
    let mut cells = 0;
    for row in groebner(model, rows)? {
        if row.leading()[..prefix].iter().all(|&power| power == 0) {
            let cost = rv_poly::cells(&row.polynomial);
            cells += cost;
            check_storage(cells)?;
            model.spend(cost)?;
            result.push(
                row.polynomial
                    .iter()
                    .map(|term| Rc::from(&term[prefix..]))
                    .collect(),
            );
        }
    }
    Ok(result)
}

pub(crate) fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain(
            "symbolic ideal calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

fn divides(divisor: &[usize], target: &[usize]) -> bool {
    divisor.iter().zip(target).all(|(&a, &b)| a <= b)
}

fn quotient(target: &[usize], divisor: &[usize]) -> Monomial {
    target.iter().zip(divisor).map(|(&a, &b)| a - b).collect()
}

fn join(left: &[usize], right: &[usize]) -> Monomial {
    left.iter().zip(right).map(|(&a, &b)| a.max(b)).collect()
}

fn shift(model: &mut Model, polynomial: &[Monomial], powers: &[usize]) -> Result<Polynomial> {
    model.spend(rv_poly::cells(polynomial))?;
    polynomial
        .iter()
        .map(|term| rv_poly::product_term(term, powers))
        .collect()
}

fn reduce(
    model: &mut Model,
    polynomial: Polynomial,
    basis: &[Row],
    omit: Option<usize>,
    mut certificate: Option<&mut Certificate>,
) -> Result<Polynomial> {
    if polynomial.is_empty() || basis.is_empty() {
        return Ok(polynomial);
    }
    let width = polynomial[0].len();
    model.spend(rv_poly::cells(&polynomial) * (1 + polynomial.len().ilog2() as usize))?;
    let mut pending: BTreeSet<_> = polynomial.iter().cloned().collect();
    let mut remainder = Vec::new();
    while let Some(leading) = pending.pop_last() {
        model.spend(1 + width)?;
        let mut divisor = None;
        for (index, row) in basis.iter().enumerate() {
            model.spend(1 + width)?;
            if Some(index) != omit && divides(row.leading(), &leading) {
                divisor = Some(row);
                break;
            }
        }
        let Some(row) = divisor else {
            remainder.push(leading);
            continue;
        };
        let powers = quotient(&leading, row.leading());
        // The leading term cancels. A fixed monomial shift preserves lex order,
        // so every replacement term is smaller than the cancelled term.
        for term in &row.polynomial[..row.polynomial.len() - 1] {
            model.spend((1 + width) * (1 + pending.len().max(1).ilog2() as usize))?;
            let term = rv_poly::product_term(term, &powers)?;
            if !pending.remove(&term) {
                pending.insert(term);
                check_storage((1 + width) * (remainder.len() + pending.len()))?;
            }
        }
        if let Some(certificate) = certificate.as_deref_mut() {
            certificate.add_shifted(model, &row.certificate, &powers)?;
        }
    }
    remainder.reverse();
    Ok(remainder.into())
}

fn groebner(model: &mut Model, rows: impl IntoIterator<Item = Row>) -> Result<Vec<Row>> {
    let mut basis = Vec::new();
    let mut storage = 0;
    for row in rows {
        if !row.polynomial.is_empty() {
            append_reduced(model, &mut basis, &mut storage, row)?;
        }
    }

    // Enumerate all pairs in insertion order, including pairs with every newly
    // appended remainder, without retaining a quadratic queue of pairs.
    let mut j = 1;
    while j < basis.len() {
        for i in 0..j {
            if let Some(row) = s_polynomial(model, &basis[i], &basis[j])? {
                append_reduced(model, &mut basis, &mut storage, row)?;
            }
        }
        j += 1;
    }
    minimal_reduced(model, basis)
}

fn append_reduced(
    model: &mut Model,
    basis: &mut Vec<Row>,
    storage: &mut usize,
    mut row: Row,
) -> Result<()> {
    row.polynomial = reduce(
        model,
        row.polynomial,
        basis,
        None,
        Some(&mut row.certificate),
    )?;
    if !row.polynomial.is_empty() {
        *storage += row.cells();
        check_storage(*storage)?;
        basis.push(row);
    }
    Ok(())
}

fn s_polynomial(model: &mut Model, left: &Row, right: &Row) -> Result<Option<Row>> {
    let a = left.leading();
    let b = right.leading();
    model.spend(1 + a.len())?;
    // Relatively prime initial monomials pass Buchberger's product criterion.
    if (left.polynomial.len() == 1 && right.polynomial.len() == 1)
        || !a.iter().zip(b.iter()).any(|(&x, &y)| x != 0 && y != 0)
    {
        return Ok(None);
    }
    let lcm = join(a, b);
    let degree = left.degree.map(|_| rv_poly::degree(&lcm)).transpose()?;
    let left_powers = quotient(&lcm, a);
    let right_powers = quotient(&lcm, b);
    let left_tail = shift(
        model,
        &left.polynomial[..left.polynomial.len() - 1],
        &left_powers,
    )?;
    let right_tail = shift(
        model,
        &right.polynomial[..right.polynomial.len() - 1],
        &right_powers,
    )?;
    let polynomial = rv_poly::add(model, &left_tail, &right_tail)?;
    if polynomial.is_empty() {
        return Ok(None);
    }
    let mut certificate = Certificate::default();
    certificate.add_shifted(model, &left.certificate, &left_powers)?;
    certificate.add_shifted(model, &right.certificate, &right_powers)?;
    Ok(Some(Row {
        degree,
        polynomial,
        certificate,
    }))
}

fn minimal_reduced(model: &mut Model, basis: Vec<Row>) -> Result<Vec<Row>> {
    let mut keep = vec![true; basis.len()];
    for (index, row) in basis.iter().enumerate() {
        let leading = row.leading();
        for (other, divisor) in basis.iter().enumerate() {
            model.spend(1 + leading.len())?;
            if index != other
                && (divisor.leading() != leading || other < index)
                && divides(divisor.leading(), leading)
            {
                keep[index] = false;
                break;
            }
        }
    }
    let mut basis: Vec<_> = basis
        .into_iter()
        .zip(keep)
        .filter_map(|(row, keep)| keep.then_some(row))
        .collect();
    let mut storage = basis.iter().map(Row::cells).sum::<usize>();
    for index in 0..basis.len() {
        let row = &basis[index];
        model.spend(1 + row.certificate.coefficients.len())?;
        let mut certificate = row.certificate.clone();
        let polynomial = reduce(
            model,
            row.polynomial.clone(),
            &basis,
            Some(index),
            Some(&mut certificate),
        )?;
        let reduced = Row {
            degree: row.degree,
            polynomial,
            certificate,
        };
        storage = storage - row.cells() + reduced.cells();
        check_storage(storage)?;
        basis[index] = reduced;
    }
    let width = basis.first().map_or(0, |row| row.leading().len());
    model.spend(basis.len() * (1 + width) * (1 + basis.len().max(1).ilog2() as usize))?;
    basis.sort_unstable_by(|a, b| {
        a.degree
            .cmp(&b.degree)
            .then_with(|| a.leading().cmp(b.leading()))
    });
    Ok(basis)
}

fn hilbert(model: &mut Model, width: usize, degree: usize, basis: &[Row]) -> Result<String> {
    if basis.first().is_some_and(|row| row.degree == Some(0)) {
        return Ok("0".into());
    }
    if width == 0 {
        return Ok(if degree == 0 { "1" } else { "0" }.into());
    }
    // Inclusion–exclusion over multiples of the minimal initial monomials.
    // Merge equal LCMs after each generator, discarding degrees above the query.
    let empty: Monomial = vec![0; width].into();
    let mut lcms = BTreeMap::from([(empty, 1i128)]);
    for row in basis
        .iter()
        .take_while(|row| row.degree.expect("homogeneous basis row") <= degree)
    {
        let mut changes = BTreeMap::new();
        for (term, &coefficient) in &lcms {
            model.spend((1 + width) * (1 + lcms.len().ilog2() as usize))?;
            let lcm = join(term, row.leading());
            if rv_poly::degree(&lcm).is_ok_and(|d| d <= degree) {
                change_coefficient(
                    &mut changes,
                    lcm,
                    coefficient.checked_neg().ok_or_else(count_limit)?,
                )?;
                check_storage((width + 3) * (lcms.len() + changes.len()))?;
            }
        }
        for (term, difference) in changes {
            change_coefficient(&mut lcms, term, difference)?;
        }
    }
    // The number of multiples in degree d depends only on the LCM's degree.
    let mut by_degree = BTreeMap::new();
    for (term, coefficient) in lcms {
        model.spend((1 + width) * (1 + by_degree.len().max(1).ilog2() as usize))?;
        change_coefficient(&mut by_degree, rv_poly::degree(&term)?, coefficient)?;
    }
    let mut positive = Natural::default();
    let mut negative = Natural::default();
    for (initial_degree, coefficient) in by_degree {
        let free_degree = degree - initial_degree;
        let count = Natural::compositions(model, free_degree, width)?;
        let total = if coefficient > 0 {
            &mut positive
        } else {
            &mut negative
        };
        total.add_multiple(model, &count, coefficient.unsigned_abs())?;
    }
    positive.subtract(&negative)?;
    Ok(positive.decimal())
}

fn count_limit() -> Error {
    Error::plain("Hilbert inclusion–exclusion coefficient is out of range")
}

fn change_coefficient<Key: Ord>(
    terms: &mut BTreeMap<Key, i128>,
    term: Key,
    difference: i128,
) -> Result<()> {
    use std::collections::btree_map::Entry;
    match terms.entry(term) {
        Entry::Occupied(mut entry) => {
            let coefficient = entry
                .get()
                .checked_add(difference)
                .ok_or_else(count_limit)?;
            if coefficient == 0 {
                entry.remove();
            } else {
                *entry.get_mut() = coefficient;
            }
        }
        Entry::Vacant(entry) => {
            entry.insert(difference);
        }
    }
    Ok(())
}

// Binomial counts may exceed machine integers even when their difference is
// small. A count uses at most O(width * word_bits) digits in this calculation.
const RADIX: u128 = 1_000_000_000;

#[derive(Default)]
struct Natural(Vec<u32>);

impl Natural {
    fn compositions(model: &mut Model, degree: usize, width: usize) -> Result<Self> {
        let mut value = Self(vec![1]);
        // C(degree + width - 1, width - 1), with exact division at each step.
        for i in 1..width as u128 {
            model.spend(1 + 2 * value.0.len())?;
            let mut carry = 0;
            for digit in &mut value.0 {
                carry += u128::from(*digit) * (degree as u128 + i);
                *digit = (carry % RADIX) as u32;
                carry /= RADIX;
            }
            while carry != 0 {
                value.0.push((carry % RADIX) as u32);
                carry /= RADIX;
            }
            let mut remainder = 0;
            for digit in value.0.iter_mut().rev() {
                let numerator = remainder * RADIX + u128::from(*digit);
                *digit = (numerator / i) as u32;
                remainder = numerator % i;
            }
            debug_assert_eq!(remainder, 0);
            value.normalize();
        }
        Ok(value)
    }

    fn add_multiple(&mut self, model: &mut Model, value: &Self, mut factor: u128) -> Result<()> {
        let mut offset = 0;
        while factor != 0 {
            model.spend(1 + value.0.len())?;
            let multiplier = (factor % RADIX) as u64;
            factor /= RADIX;
            let mut carry = 0;
            let size = offset + value.0.len();
            self.0.resize(size.max(self.0.len()), 0);
            for (index, &digit) in value.0.iter().enumerate() {
                let target = &mut self.0[offset + index];
                carry += u64::from(*target) + u64::from(digit) * multiplier;
                *target = (carry % RADIX as u64) as u32;
                carry /= RADIX as u64;
            }
            let mut index = size;
            while carry != 0 {
                model.spend(1)?;
                if index == self.0.len() {
                    self.0.push(0);
                }
                carry += u64::from(self.0[index]);
                self.0[index] = (carry % RADIX as u64) as u32;
                carry /= RADIX as u64;
                index += 1;
            }
            offset += 1;
        }
        Ok(())
    }

    fn subtract(&mut self, other: &Self) -> Result<()> {
        self.normalize();
        let magnitude = self
            .0
            .len()
            .cmp(&other.0.len())
            .then_with(|| self.0.iter().rev().cmp(other.0.iter().rev()));
        if magnitude == Ordering::Less {
            return Err(Error::plain("Hilbert count is unexpectedly negative"));
        }
        let mut borrow = 0;
        for (index, digit) in self.0.iter_mut().enumerate() {
            let difference =
                i64::from(*digit) - i64::from(other.0.get(index).copied().unwrap_or(0)) - borrow;
            borrow = i64::from(difference < 0);
            *digit = (difference + borrow * RADIX as i64) as u32;
        }
        debug_assert_eq!(borrow, 0);
        self.normalize();
        Ok(())
    }

    fn normalize(&mut self) {
        while self.0.last() == Some(&0) {
            self.0.pop();
        }
    }

    fn decimal(self) -> String {
        use std::fmt::Write;
        let Some((first, digits)) = self.0.split_last() else {
            return "0".into();
        };
        let mut output = first.to_string();
        for digit in digits.iter().rev() {
            let _ = write!(output, "{digit:09}");
        }
        output
    }
}
