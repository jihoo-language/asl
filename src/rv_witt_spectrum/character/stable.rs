use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::wide_integer::Natural;

use super::super::{add, multiply, zero, Integer, SchurShape};
use super::table::{
    accumulate, beta_numbers, factorial, from_power_basis, pack, partitions, power_basis,
    CharacterOracle, Coefficients, Partition, MAX_WEIGHT,
};
use super::{
    add_terms, charge_compare, coefficient_report, render_terms, scale_terms, store_terms,
    Character, Terms,
};

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct StableCharacter(pub(super) Terms);

impl StableCharacter {
    pub(crate) fn seed(tail: SchurShape) -> Self {
        Self(Terms::from([(tail, 1)]))
    }

    pub(crate) fn zero() -> Self {
        Self(Terms::new())
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        store_terms(model, &self.0)?;
        Ok(Value::RVStableSchurCharacter(Rc::new(self)))
    }

    pub(crate) fn render(&self) -> String {
        format!("RVStableSchurCharacter{{{}}}", render_terms(&self.0))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        charge_compare(model, &self.0, &other.0)?;
        Ok(self.cmp(other))
    }

    pub(crate) fn coefficient(&self, shape: &SchurShape) -> String {
        coefficient_report("rvStableSchurCharacterCoeff", &self.0, shape)
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        Ok(Self(add_terms(model, &self.0, &rhs.0)?))
    }

    pub(crate) fn scale(&self, model: &mut Model, factor: i128) -> Result<Self> {
        Ok(Self(scale_terms(model, &self.0, factor)?))
    }

    pub(super) fn highest_weight(&self) -> u128 {
        self.0.keys().map(SchurShape::weight).max().unwrap_or(0)
    }

    pub(crate) fn product(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        if self.0.is_empty() || rhs.0.is_empty() {
            return Ok(Self::zero());
        }
        let a = self.highest_weight();
        let b = rhs.highest_weight();
        if a == 0 {
            return rhs.scale(model, *self.0.values().next().unwrap());
        }
        if b == 0 {
            return self.scale(model, *rhs.0.values().next().unwrap());
        }
        let degree = a
            .checked_add(b)
            .filter(|&n| n <= u128::from(MAX_WEIGHT))
            .ok_or_else(|| {
                Error::plain(
                    "reduced Kronecker product exceeds the combined tail limit of 16 boxes",
                )
            })? as u8;
        let mut expansions = Expansions::new(model, degree)?;
        let left = expansions.character(model, &self.0)?;
        let right = expansions.character(model, &rhs.0)?;
        let polynomial = binomial_product(model, &left, &right)?;
        expansions.reduced_basis(model, degree, polynomial)
    }

    pub(crate) fn specialize(&self, model: &mut Model, rank: u128) -> Result<Character> {
        let mut wide: BTreeMap<SchurShape, Integer> = BTreeMap::new();
        for (tail, &coefficient) in &self.0 {
            if let Some((shape, negative)) = straighten(model, tail, rank)? {
                let value = Integer::from_i128(coefficient);
                let value = if negative { value.negate() } else { value };
                let entry = wide.entry(shape).or_insert_with(zero);
                *entry = add(*entry, value)?;
            }
        }
        let terms = wide
            .into_iter()
            .filter(|(_, value)| !value.is_zero())
            .map(|(shape, value)| Ok((shape, value.signed().ok_or_else(super::overflow)?)))
            .collect::<Result<Terms>>()?;
        Ok(Character {
            weight: rank,
            terms,
        })
    }
}

pub(super) struct Expansions {
    finite: CharacterOracle,
    classes: Vec<Vec<Partition>>,
    stable: BTreeMap<Partition, Coefficients>,
}

impl Expansions {
    pub(super) fn new(model: &mut Model, degree: u8) -> Result<Self> {
        let classes = (0..=degree)
            .map(|n| partitions(model, n))
            .collect::<Result<_>>()?;
        Ok(Self {
            finite: CharacterOracle::default(),
            classes,
            stable: BTreeMap::new(),
        })
    }

    pub(super) fn character(&mut self, model: &mut Model, terms: &Terms) -> Result<Coefficients> {
        let mut polynomial = BTreeMap::new();
        for (shape, &coefficient) in terms {
            let tail = small_tail(shape);
            let expansion = self.expansion(model, &tail)?;
            add_scaled(
                model,
                &mut polynomial,
                expansion,
                Integer::from_i128(coefficient),
            )?;
        }
        Ok(polynomial)
    }

    pub(super) fn cycle_types(&self, degree: u8) -> impl Iterator<Item = &Partition> {
        self.classes[..=usize::from(degree)].iter().flatten()
    }

    fn expansion(&mut self, model: &mut Model, tail: &Partition) -> Result<&Coefficients> {
        if !self.stable.contains_key(tail) {
            let runs: Vec<_> = tail
                .chunk_by(|a, b| a == b)
                .map(|run| (run[0], run.len()))
                .collect();
            let mut result = BTreeMap::new();
            // First-row Jacobi–Trudi cofactors are s_{tail/(1^j)}.
            // Vertical Pieri removes j boxes, at most one from each row.
            vertical_strips(
                model,
                &runs,
                &mut Vec::new(),
                false,
                &mut |model, shape, negative| {
                    let degree: u8 = shape.iter().sum();
                    let beta = beta_numbers(shape.iter().copied());
                    for cycle_type in &self.classes[degree as usize] {
                        let value = self.finite.value(model, beta, pack(cycle_type))?;
                        let value = Integer::from_i128(if negative { -value } else { value });
                        accumulate(&mut result, cycle_type.clone(), value)?;
                    }
                    Ok(())
                },
            )?;
            self.stable.insert(tail.clone(), result);
        }
        Ok(&self.stable[tail])
    }

    pub(super) fn reduced_basis(
        &mut self,
        model: &mut Model,
        degree: u8,
        mut polynomial: Coefficients,
    ) -> Result<StableCharacter> {
        let mut terms = Terms::new();
        for n in (0..=degree).rev() {
            let top = polynomial
                .iter()
                .filter(|(class, _)| class.iter().sum::<u8>() == n)
                .map(|(class, &value)| (class.clone(), value))
                .collect::<Coefficients>();
            if top.is_empty() {
                continue;
            }
            // The top binomial layer is the ordinary irreducible character
            // table. Subtract its stable lifts before descending in weight.
            let basis = power_basis(model, n, top)?;
            let top = from_power_basis(model, n, factorial(n), &basis, &mut self.finite)?;
            for (shape, coefficient) in top.terms {
                let tail = small_tail(&shape);
                let expansion = self.expansion(model, &tail)?;
                add_scaled(
                    model,
                    &mut polynomial,
                    expansion,
                    Integer::from_i128(coefficient).negate(),
                )?;
                terms.insert(shape, coefficient);
            }
        }
        if !polynomial.is_empty() {
            return Err(Error::plain(
                "stable Kronecker basis conversion did not exhaust the character polynomial",
            ));
        }
        Ok(StableCharacter(terms))
    }
}

fn small_tail(shape: &SchurShape) -> Partition {
    shape.rows().iter().map(|&row| row as u8).collect()
}

pub(super) fn add_scaled(
    model: &mut Model,
    result: &mut Coefficients,
    source: &Coefficients,
    scale: Integer,
) -> Result<()> {
    for (class, &value) in source {
        let value = multiply(model, value, scale)?;
        accumulate(result, class.clone(), value)?;
    }
    Ok(())
}

fn vertical_strips(
    model: &mut Model,
    runs: &[(u8, usize)],
    rows: &mut Partition,
    negative: bool,
    visit: &mut impl FnMut(&mut Model, &[u8], bool) -> Result<()>,
) -> Result<()> {
    let Some((&(width, count), rest)) = runs.split_first() else {
        return visit(model, rows, negative);
    };
    let start = rows.len();
    for removed in 0..=count {
        model.spend(count + 1)?;
        rows.resize(start + count - removed, width);
        if width > 1 {
            rows.resize(start + count, width - 1);
        }
        vertical_strips(
            model,
            rest,
            rows,
            negative ^ !removed.is_multiple_of(2),
            visit,
        )?;
        rows.truncate(start);
    }
    Ok(())
}

pub(super) fn binomial_product(
    model: &mut Model,
    left: &Coefficients,
    right: &Coefficients,
) -> Result<Coefficients> {
    let mut product = Coefficients::new();
    for (left, &a) in left {
        for (right, &b) in right {
            let value = multiply(model, a, b)?;
            let groups = counts(left, right);
            overlaps(model, &groups, &mut Vec::new(), value, &mut product)?;
        }
    }
    Ok(product)
}

fn counts(left: &[u8], right: &[u8]) -> Vec<(u8, u8, u8)> {
    let mut counts = [[0u8; MAX_WEIGHT as usize + 1]; 2];
    for (j, partition) in [left, right].into_iter().enumerate() {
        for &part in partition {
            counts[j][part as usize] += 1;
        }
    }
    (1..=MAX_WEIGHT)
        .rev()
        .filter_map(|i| {
            let (a, b) = (counts[0][i as usize], counts[1][i as usize]);
            (a != 0 || b != 0).then_some((i, a, b))
        })
        .collect()
}

fn overlaps(
    model: &mut Model,
    groups: &[(u8, u8, u8)],
    cycle_type: &mut Partition,
    coefficient: Integer,
    result: &mut Coefficients,
) -> Result<()> {
    let Some((&(length, a, b), rest)) = groups.split_first() else {
        return accumulate(result, cycle_type.clone(), coefficient);
    };
    let start = cycle_type.len();
    // Count pairs of subsets of sizes a and b with a fixed union and
    // intersection of size j, for cycles of one length at a time.
    for shared in 0..=a.min(b) {
        let count = a + b - shared;
        model.spend(count as usize + 1)?;
        cycle_type.resize(start + count as usize, length);
        let ways = choose(count, a) * choose(a, shared);
        let scale = multiply(model, coefficient, Integer::from_i128(ways))?;
        overlaps(model, rest, cycle_type, scale, result)?;
        cycle_type.truncate(start);
    }
    Ok(())
}

pub(super) fn choose(n: u8, k: u8) -> i128 {
    (1..=k).fold(1, |value, i| value * i128::from(n - i + 1) / i128::from(i))
}

fn straighten(
    model: &mut Model,
    tail: &SchurShape,
    rank: u128,
) -> Result<Option<(SchurShape, bool)>> {
    let base = natural(rank);
    let weight = natural(tail.weight());
    // The shifted tail entries are strictly decreasing. Only the first
    // entry of the completed sequence can move. Compare the two sides
    // of N+j versus |tail|+tail_j in wide arithmetic.
    let mut position = 0;
    for &row in tail.rows() {
        model.spend(1)?;
        let left = add(base, natural((position + 1) as u128))?;
        let right = add(weight, natural(row))?;
        match left.magnitude.cmp(right.magnitude) {
            Ordering::Equal => return Ok(None),
            Ordering::Greater => break,
            Ordering::Less => position += 1,
        }
    }
    let inserted = add(add(base, natural(position as u128))?, weight.negate())?;
    if inserted.negative && !inserted.is_zero() {
        return Ok(None);
    }
    let inserted = inserted.magnitude.small().ok_or_else(super::overflow)?;
    let (before, after) = tail.rows().split_at(position);
    let mut rows = before
        .iter()
        .map(|&x| x - 1)
        .chain([inserted])
        .chain(after.iter().copied())
        .collect::<Vec<_>>();
    while rows.last() == Some(&0) {
        rows.pop();
    }
    model.spend(rows.len() + 1)?;
    let shape = SchurShape::new(rows)?;
    if shape.weight() != rank {
        return Err(Error::plain(
            "stable Schur specialization has the wrong weight",
        ));
    }
    Ok(Some((shape, !position.is_multiple_of(2))))
}

fn natural(value: u128) -> Integer {
    Integer {
        magnitude: Natural::from_small(value),
        negative: false,
    }
}
