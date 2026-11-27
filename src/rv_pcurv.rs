use std::cmp::Ordering;
use std::collections::{btree_map::Entry, BTreeMap, HashMap};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::rv_qconnection;
use crate::rv_qgr::{self, Grade, Schubert, Shape};

const MAX_CELLS: usize = 300_000;
const MAX_PRIME: usize = 4096;
const MAX_WEDGE_AMBIENT: usize = 256;
const MAX_CHAR_RANK: usize = 64;

pub(crate) type Support = Rc<[(Schubert, u32)]>;
type Polynomial = Vec<(usize, u32)>;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct ModGrade {
    pub(crate) prime: Prime,
    pub(crate) grassmannian: Grade,
}

impl ModGrade {
    pub(crate) fn checked(model: &mut Model, p: usize, k: usize, n: usize) -> Result<Self> {
        let grassmannian = Grade::checked(k, n)?;
        let prime = checked_prime(model, p)?;
        Ok(Self {
            prime,
            grassmannian,
        })
    }
}

pub(crate) fn checked_prime(model: &mut Model, number: usize) -> Result<Prime> {
    if number == 2 {
        return Ok(Prime::validated(2));
    }
    if number < 2 || number.is_multiple_of(2) {
        return Err(Error::plain(format!(
            "quantum connection characteristic {number} is not prime"
        )));
    }
    Prime::new(model, number)
}

fn steps(prime: Prime) -> Result<usize> {
    let p = prime.get() as usize;
    if p > MAX_PRIME {
        return Err(Error::plain(
            "quantum p-curvature characteristic exceeds the limit of 4096",
        ));
    }
    Ok(p)
}

fn storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain(
            "quantum p-curvature calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

fn overflow() -> Error {
    Error::plain("quantum p-curvature polynomial exponent is out of range")
}

fn cells(term: &Schubert) -> usize {
    4 + term.shape.len()
}

fn finite(prime: Prime, integer: i128) -> u32 {
    integer.rem_euclid(i128::from(prime.get())) as u32
}

struct Sum<K> {
    terms: BTreeMap<K, u32>,
    cells: usize,
}

impl<K: Ord> Default for Sum<K> {
    fn default() -> Self {
        Self {
            terms: BTreeMap::new(),
            cells: 0,
        }
    }
}

impl<K: Ord> Sum<K> {
    fn add(&mut self, model: &mut Model, prime: Prime, key: K, count: usize, c: u32) -> Result<()> {
        model.spend(count * (1 + self.terms.len().max(1).ilog2() as usize))?;
        if c == 0 {
            return Ok(());
        }
        match self.terms.entry(key) {
            Entry::Vacant(slot) => {
                self.cells += count;
                storage(self.cells)?;
                slot.insert(c);
            }
            Entry::Occupied(mut slot) => {
                let value = prime.add(*slot.get(), c);
                if value == 0 {
                    slot.remove();
                    self.cells -= count;
                } else {
                    *slot.get_mut() = value;
                }
            }
        }
        Ok(())
    }

    fn finish(self) -> Vec<(K, u32)> {
        self.terms.into_iter().collect()
    }
}

impl Sum<Schubert> {
    fn term(&mut self, model: &mut Model, prime: Prime, term: Schubert, c: u32) -> Result<()> {
        let count = cells(&term);
        self.add(model, prime, term, count, c)
    }
}

fn persist(model: &mut Model, terms: Vec<(Schubert, u32)>) -> Result<Value> {
    let mut count = 0;
    for (term, _) in &terms {
        count += cells(term);
        storage(count)?;
    }
    model.allocate(count)?;
    Ok(Value::RVQGrP(terms.into()))
}

pub(crate) fn zero() -> Value {
    Value::RVQGrP([].into())
}

pub(crate) fn reduction(
    model: &mut Model,
    grade: ModGrade,
    input: &rv_qgr::Support,
) -> Result<Value> {
    model.spend(input.len())?;
    let terms = input
        .iter()
        .filter_map(|(term, c)| {
            let coefficient = finite(grade.prime, *c);
            (coefficient != 0).then(|| (term.clone(), coefficient))
        })
        .collect();
    persist(model, terms)
}

pub(crate) fn basis(model: &mut Model, grade: ModGrade, parts: Vec<usize>) -> Result<Value> {
    let shape = grade.grassmannian.shape(model, parts)?;
    persist(model, vec![(Schubert { quantum: 0, shape }, 1)])
}

pub(crate) fn compare(model: &mut Model, left: &Support, right: &Support) -> Result<Ordering> {
    if Rc::ptr_eq(left, right) {
        return Ok(Ordering::Equal);
    }
    for (a, b) in left.iter().zip(right.iter()) {
        model.spend(cells(&a.0))?;
        let comparison = a.cmp(b);
        if !comparison.is_eq() {
            return Ok(comparison);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

pub(crate) fn add(model: &mut Model, p: Prime, a: &Support, b: &Support) -> Result<Value> {
    let mut sum = Sum::default();
    for (term, c) in a.iter().chain(b.iter()) {
        sum.term(model, p, term.clone(), *c)?;
    }
    persist(model, sum.finish())
}

pub(crate) fn scale(model: &mut Model, p: Prime, input: &Support, scalar: i128) -> Result<Value> {
    let scalar = finite(p, scalar);
    if scalar == 0 {
        return Ok(zero());
    }
    model.spend(input.len())?;
    persist(
        model,
        input
            .iter()
            .map(|(term, c)| (term.clone(), p.mul(*c, scalar)))
            .collect(),
    )
}

pub(crate) fn shift(model: &mut Model, input: &Support, power: usize) -> Result<Value> {
    model.spend(input.len())?;
    let terms = input
        .iter()
        .map(|(term, c)| {
            let quantum = term.quantum.checked_add(power).ok_or_else(overflow)?;
            Ok((
                Schubert {
                    quantum,
                    shape: term.shape.clone(),
                },
                *c,
            ))
        })
        .collect::<Result<_>>()?;
    persist(model, terms)
}

pub(crate) fn transpose(model: &mut Model, p: Prime, input: &Support) -> Result<Value> {
    let mut sum = Sum::default();
    for (term, c) in input.iter() {
        let shape = rv_qgr::transpose_shape(model, &term.shape, sum.cells)?;
        sum.term(
            model,
            p,
            Schubert {
                quantum: term.quantum,
                shape,
            },
            *c,
        )?;
    }
    persist(model, sum.finish())
}

pub(crate) fn classical(model: &mut Model, input: &Support) -> Result<Value> {
    let length = input.partition_point(|(term, _)| term.quantum == 0);
    model.spend(length)?;
    persist(model, input[..length].to_vec())
}

pub(crate) fn coefficient(
    model: &mut Model,
    grade: ModGrade,
    terms: &Support,
    degree: usize,
    parts: Vec<usize>,
) -> Result<u32> {
    let shape = grade.grassmannian.shape(model, parts)?;
    let term = Schubert {
        quantum: degree,
        shape,
    };
    model.spend(cells(&term) * (1 + terms.len().max(1).ilog2() as usize))?;
    Ok(terms
        .binary_search_by(|(entry, _)| entry.cmp(&term))
        .map(|i| terms[i].1)
        .unwrap_or(0))
}

fn nabla_support(
    model: &mut Model,
    grade: ModGrade,
    terms: &[(Schubert, u32)],
) -> Result<Vec<(Schubert, u32)>> {
    let mut sum = Sum::default();
    let p = grade.prime;
    for (term, c) in terms {
        sum.term(model, p, term.clone(), p.mul(*c, p.reduce(term.quantum)))?;
        for mut image in rv_qgr::hyperplane(model, grade.grassmannian, &term.shape, sum.cells)? {
            image.quantum = image
                .quantum
                .checked_add(term.quantum)
                .ok_or_else(overflow)?;
            sum.term(model, p, image, *c)?;
        }
    }
    Ok(sum.finish())
}

pub(crate) fn nabla(model: &mut Model, grade: ModGrade, input: &Support) -> Result<Value> {
    let output = nabla_support(model, grade, input)?;
    persist(model, output)
}

fn direct_column(
    model: &mut Model,
    grade: ModGrade,
    shape: &Shape,
) -> Result<Vec<(Schubert, u32)>> {
    let steps = steps(grade.prime)?;
    model.spend(steps)?;
    let basis = [(
        Schubert {
            quantum: 0,
            shape: shape.clone(),
        },
        1,
    )];
    let first = nabla_support(model, grade, &basis)?;
    let mut power = first.clone();
    for _ in 1..steps {
        power = nabla_support(model, grade, &power)?;
        if power.is_empty() {
            break;
        }
    }
    let mut sum = Sum::default();
    let p = grade.prime;
    for (term, c) in power {
        sum.term(model, p, term, c)?;
    }
    for (term, c) in first {
        sum.term(model, p, term, p.negate(c))?;
    }
    Ok(sum.finish())
}

pub(crate) fn column(model: &mut Model, grade: ModGrade, parts: Vec<usize>) -> Result<Value> {
    let shape = grade.grassmannian.shape(model, parts)?;
    let output = direct_column(model, grade, &shape)?;
    persist(model, output)
}

fn projective_step(
    model: &mut Model,
    p: Prime,
    n: usize,
    wrap: u32,
    terms: &[((usize, usize), u32)],
) -> Result<Vec<((usize, usize), u32)>> {
    let mut sum = Sum::default();
    for &((q, beta), c) in terms {
        sum.add(model, p, (q, beta), 3, p.mul(c, p.reduce(q)))?;
        if beta + 1 == n {
            sum.add(model, p, (q + 1, 0), 3, p.mul(c, wrap))?;
        } else {
            sum.add(model, p, (q, beta + 1), 3, c)?;
        }
    }
    Ok(sum.finish())
}

fn projective_column(
    model: &mut Model,
    p: Prime,
    n: usize,
    wrap: u32,
    beta: usize,
) -> Result<Vec<((usize, usize), u32)>> {
    let steps = steps(p)?;
    model.spend(steps)?;
    let first = projective_step(model, p, n, wrap, &[((0, beta), 1)])?;
    let mut power = first.clone();
    for _ in 1..steps {
        power = projective_step(model, p, n, wrap, &power)?;
        if power.is_empty() {
            break;
        }
    }
    let mut sum = Sum::default();
    for (term, c) in power {
        sum.add(model, p, term, 3, c)?;
    }
    for (term, c) in first {
        sum.add(model, p, term, 3, p.negate(c))?;
    }
    Ok(sum.finish())
}

fn wedge_column(model: &mut Model, grade: ModGrade, shape: &Shape) -> Result<Vec<(Schubert, u32)>> {
    let Grade { rank, ambient } = grade.grassmannian;
    if ambient > MAX_WEDGE_AMBIENT {
        return Err(Error::plain(
            "quantum Satake ambient dimension exceeds the limit of 256",
        ));
    }
    let p = grade.prime;
    let wrap = if rank.is_multiple_of(2) {
        p.negate(1)
    } else {
        1
    };
    let beta: Vec<_> = (0..rank)
        .rev()
        .map(|i| shape.get(i).copied().unwrap_or(0) + rank - 1 - i)
        .collect();
    let mut sum = Sum::default();
    // The p-curvature of an exterior-power connection is the sum of the
    // actions on each factor. The tensor-factor operators commute in the
    // relative p-power formula, so there are no mixed terms in char p.
    for (i, &factor) in beta.iter().enumerate() {
        for ((quantum, replacement), mut c) in projective_column(model, p, ambient, wrap, factor)? {
            model.spend(rank)?;
            let mut next = beta.clone();
            next.remove(i);
            let position = match next.binary_search(&replacement) {
                Ok(_) => continue,
                Err(position) => position,
            };
            if !i.abs_diff(position).is_multiple_of(2) {
                c = p.negate(c);
            }
            next.insert(position, replacement);
            next.reverse();
            let shape = rv_qgr::beta_shape(&next);
            sum.term(model, p, Schubert { quantum, shape }, c)?;
        }
    }
    Ok(sum.finish())
}

pub(crate) fn apply(
    model: &mut Model,
    grade: ModGrade,
    input: &Support,
    satake: bool,
) -> Result<Value> {
    if input.is_empty() {
        return Ok(zero());
    }
    steps(grade.prime)?;
    let mut sum = Sum::default();
    let mut columns = BTreeMap::new();
    let mut kept = 0;
    for (term, c) in input.iter() {
        model.spend(cells(term) * (1 + columns.len().max(1).ilog2() as usize))?;
        let column = match columns.entry(term.shape.clone()) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => {
                let result = if satake {
                    wedge_column(model, grade, &term.shape)?
                } else {
                    direct_column(model, grade, &term.shape)?
                };
                kept += 2 + term.shape.len();
                for (entry, _) in &result {
                    kept += cells(entry);
                    storage(kept)?;
                }
                storage(kept)?;
                entry.insert(result)
            }
        };
        for (image, scalar) in column {
            let quantum = image
                .quantum
                .checked_add(term.quantum)
                .ok_or_else(overflow)?;
            let image = Schubert {
                quantum,
                shape: image.shape.clone(),
            };
            sum.term(model, grade.prime, image, grade.prime.mul(*c, *scalar))?;
        }
    }
    persist(model, sum.finish())
}

fn matrix(model: &mut Model, mut grade: ModGrade) -> Result<Vec<Vec<Polynomial>>> {
    let step_count = steps(grade.prime)?;
    // Transposition preserves the connection and its characteristic polynomial.
    if grade.grassmannian.rank > grade.grassmannian.width() {
        grade.grassmannian = grade.grassmannian.complement();
    }
    let shapes = rv_qconnection::schubert_basis(model, grade.grassmannian, MAX_CHAR_RANK)?;
    let size = shapes.len();
    model.spend(size * step_count)?;
    let positions: HashMap<_, _> = shapes
        .iter()
        .enumerate()
        .map(|(i, shape)| (shape, i))
        .collect();
    let mut matrix = vec![vec![Vec::new(); size]; size];
    let mut count = size * size;
    for (i, shape) in shapes.iter().enumerate() {
        for (term, coefficient) in direct_column(model, grade, shape)? {
            count += 2;
            storage(count)?;
            matrix[positions[&term.shape]][i].push((term.quantum, coefficient));
        }
    }
    Ok(matrix)
}

fn accumulate(
    model: &mut Model,
    prime: Prime,
    polynomial: &mut Sum<usize>,
    a: &Polynomial,
    b: &Polynomial,
    sign: u32,
) -> Result<()> {
    for &(i, c) in a {
        for &(j, d) in b {
            let q = i.checked_add(j).ok_or_else(overflow)?;
            polynomial.add(model, prime, q, 2, prime.mul(sign, prime.mul(c, d)))?;
        }
    }
    Ok(())
}

fn push(poly: Polynomial, vector: &mut Vec<Polynomial>, kept: &mut usize) -> Result<()> {
    *kept += 1 + 2 * poly.len();
    storage(*kept)?;
    vector.push(poly);
    Ok(())
}

fn characteristic(model: &mut Model, grade: ModGrade) -> Result<Vec<Polynomial>> {
    let matrix = matrix(model, grade)?;
    let p = grade.prime;
    let negative = p.negate(1);
    let one = vec![(0, 1)];
    let mut coefficients = vec![one.clone()];
    // Division-free Berkowitz recurrence. For [a R; S M], convolve the
    // descending coefficients for M with [1, -a, -RS, -RMS, ...].
    for start in (0..matrix.len()).rev() {
        let degree = matrix.len() - start;
        let tail = start + 1;
        let mut series = vec![one.clone()];
        let mut kept = 3;
        let scalar = matrix[start][start]
            .iter()
            .map(|&(q, c)| (q, p.negate(c)))
            .collect();
        push(scalar, &mut series, &mut kept)?;
        let mut image: Vec<_> = matrix[tail..]
            .iter()
            .map(|row| row[start].clone())
            .collect();
        for power in 2..=degree {
            let mut dot = Sum::default();
            model.spend(1 + image.len())?;
            for (a, b) in matrix[start][tail..].iter().zip(&image) {
                accumulate(model, p, &mut dot, a, b, negative)?;
            }
            push(dot.finish(), &mut series, &mut kept)?;
            if power < degree {
                let mut next = Vec::with_capacity(image.len());
                let mut image_cells = 0;
                for row in &matrix[tail..] {
                    let mut dot = Sum::default();
                    model.spend(1 + image.len())?;
                    for (a, b) in row[tail..].iter().zip(&image) {
                        accumulate(model, p, &mut dot, a, b, 1)?;
                    }
                    push(dot.finish(), &mut next, &mut image_cells)?;
                }
                image = next;
            }
        }
        let mut next = Vec::with_capacity(degree + 1);
        let mut coefficients_cells = 0;
        for i in 0..=degree {
            let mut coefficient = Sum::default();
            model.spend(1 + (i + 1).min(coefficients.len()))?;
            for (j, value) in coefficients.iter().enumerate().take(i + 1) {
                accumulate(model, p, &mut coefficient, &series[i - j], value, 1)?;
            }
            push(coefficient.finish(), &mut next, &mut coefficients_cells)?;
        }
        coefficients = next;
    }
    coefficients.reverse();
    Ok(coefficients)
}

#[derive(Clone, Copy)]
pub(crate) enum Report {
    Characteristic,
    Coefficient(usize),
    Horizontal,
}

pub(crate) fn report(model: &mut Model, grade: ModGrade, report: Report) -> Result<String> {
    let mut coefficients = characteristic(model, grade)?;
    if matches!(report, Report::Horizontal) {
        let prime = grade.prime.get() as usize;
        for (q, _) in coefficients.iter_mut().flatten() {
            if !q.is_multiple_of(prime) {
                return Err(Error::plain(
                    "p-curvature characteristic polynomial violated horizontality",
                ));
            }
            *q /= prime;
        }
    }
    let result = match report {
        Report::Characteristic => format!("rvQGrCurvatureChar = {coefficients:?}"),
        Report::Coefficient(j) => {
            let value = coefficients.get(j).map_or(&[][..], Vec::as_slice);
            format!("rvQGrCurvatureCoeff = {value:?}")
        }
        Report::Horizontal => format!("rvQGrCurvatureHorizontal = {coefficients:?}"),
    };
    if result.len() > 65_536 {
        return Err(Error::plain(
            "shown p-curvature polynomial exceeds the output limit",
        ));
    }
    Ok(result)
}
