use std::rc::Rc;

use crate::adams;
use crate::cobar::{self, Word};
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::profile::Profile;

pub(crate) type Bidegree = (usize, usize);
pub(crate) type Input = (Bidegree, Rc<[usize]>);
type Cochain = Vec<Word>;

pub(crate) fn odd_stem((s, t): Bidegree) -> bool {
    (s ^ t) & 1 != 0
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Product,
    Lift,
}

impl Operation {
    pub(crate) fn from_name(name: &str) -> Option<Self> {
        match name {
            "ainfty" => Some(Self::Product),
            "alift" => Some(Self::Lift),
            _ => None,
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Product => "ainfty",
            Self::Lift => "alift",
        }
    }

    pub(crate) fn require_arity(self, arity: usize) -> Result<()> {
        if arity < 2 {
            let name = self.name();
            return Err(Error::plain(format!(
                "{name} requires at least two Ext classes"
            )));
        }
        Ok(())
    }

    pub(crate) fn bidegree(
        self,
        inputs: impl IntoIterator<Item = Bidegree>,
        arity: usize,
    ) -> Result<Bidegree> {
        self.require_arity(arity)?;
        let shift = match self {
            Self::Product => arity - 2,
            Self::Lift => arity - 1,
        };
        sum_degrees(inputs, shift, self.name())
    }

    fn zero(self) -> Value {
        match self {
            Self::Product => Value::Class([].into()),
            Self::Lift => cobar::zero(),
        }
    }
}

pub(crate) fn sum_degrees(
    inputs: impl IntoIterator<Item = Bidegree>,
    shift: usize,
    name: &str,
) -> Result<Bidegree> {
    let mut s = 0u128;
    let mut t = 0usize;
    for (length, degree) in inputs {
        s += length as u128;
        t = t
            .checked_add(degree)
            .ok_or_else(|| Error::plain(format!("{name} internal degree is out of range")))?;
    }
    let s = s
        .checked_sub(shift as u128)
        .ok_or_else(|| Error::plain(format!("{name} produces a negative cohomological degree")))?;
    let s = usize::try_from(s)
        .map_err(|_| Error::plain(format!("{name} cohomological degree is out of range")))?;
    Ok((s, t))
}

pub(crate) fn evaluate(
    model: &mut Model,
    profile: Profile,
    inputs: &[Input],
    operation: Operation,
) -> Result<Value> {
    let arity = inputs.len();
    let (result_s, result_t) =
        operation.bidegree(inputs.iter().map(|&(degree, _)| degree), arity)?;
    model.spend(arity)?;
    if inputs.iter().any(|(_, coordinates)| coordinates.is_empty()) {
        return Ok(operation.zero());
    }
    // Ext(0,0) is spanned by the strict unit, and other Ext(0,t) vanish.
    if let Some(unit) = inputs.iter().position(|((s, _), _)| *s == 0) {
        return Ok(match (operation, arity) {
            (Operation::Product, 2) => Value::Class(inputs[1 - unit].1.clone()),
            _ => operation.zero(),
        });
    }
    if let (Operation::Product, [((s, t), left), ((u, v), right)]) = (operation, inputs) {
        let left_ty = model.profile_ext_type(profile, *s, *t)?;
        let right_ty = model.profile_ext_type(profile, *u, *v)?;
        let ty = model.profile_ext_type(profile, result_s, result_t)?;
        return adams::yoneda(model, left_ty, left, right_ty, right, ty);
    }

    let intervals = Intervals::build(model, profile, inputs, arity - 1)?;
    let grade = intervals.obstruction_grade(0, arity)?;
    let source = intervals.obstruction(model, 0, arity)?;
    let (closed, coordinates) = project(model, profile, grade, &source)?;
    match operation {
        Operation::Product => {
            model.allocate(coordinates.len())?;
            Ok(Value::Class(coordinates.into()))
        }
        Operation::Lift => {
            let (_, primitive) = split_closed(model, profile, grade, &closed, &coordinates)?;
            model.allocate(primitive.len())?;
            Ok(Value::Cobar(primitive.into()))
        }
    }
}

pub(crate) struct Intervals {
    prefixes: Vec<Bidegree>,
    // Row k stores the lifts of the consecutive intervals of length k + 1.
    lifts: Vec<Vec<Cochain>>,
}

impl Intervals {
    pub(crate) fn build(
        model: &mut Model,
        profile: Profile,
        inputs: &[Input],
        through: usize,
    ) -> Result<Self> {
        let arity = inputs.len();
        let prefix_count = arity.saturating_add(1);
        let interval_count = arity.saturating_mul(prefix_count) / 2;
        model.spend(interval_count)?;
        let mut cells = interval_count.saturating_add(prefix_count);
        check_lifts(cells)?;
        let mut prefixes: Vec<Bidegree> = Vec::with_capacity(prefix_count);
        prefixes.push((0, 0));
        let mut singletons = Vec::with_capacity(arity);
        for &((s, t), ref coordinates) in inputs {
            let &(prev_s, prev_t) = prefixes.last().unwrap();
            let sum_s = prev_s
                .checked_add(s - 1)
                .ok_or_else(intermediate_overflow)?;
            let sum_t = prev_t.checked_add(t).ok_or_else(intermediate_overflow)?;
            prefixes.push((sum_s, sum_t));
            let ty = model.profile_ext_type(profile, s, t)?;
            let words = model.adams(ty)?.representative(model, coordinates)?;
            cells += words.len();
            check_lifts(cells)?;
            singletons.push(words);
        }
        let mut intervals = Self {
            prefixes,
            lifts: vec![singletons],
        };
        for length in 2..=through {
            let mut row = Vec::with_capacity(arity - length + 1);
            for start in 0..=arity - length {
                let end = start + length;
                let grade = intervals.obstruction_grade(start, end)?;
                let source = intervals.obstruction(model, start, end)?;
                let (closed, coordinates) = project(model, profile, grade, &source)?;
                let (_, primitive) = split_closed(model, profile, grade, &closed, &coordinates)?;
                cells += primitive.len();
                check_lifts(cells)?;
                row.push(primitive);
            }
            intervals.lifts.push(row);
        }
        Ok(intervals)
    }

    pub(crate) fn grade(&self, start: usize, end: usize) -> Result<Bidegree> {
        self.shifted_grade(start, end, 1)
    }

    pub(crate) fn lift(&self, start: usize, end: usize) -> &[Word] {
        &self.lifts[end - start - 1][start]
    }

    fn obstruction_grade(&self, start: usize, end: usize) -> Result<Bidegree> {
        self.shifted_grade(start, end, 2)
    }

    fn shifted_grade(&self, start: usize, end: usize, shift: usize) -> Result<Bidegree> {
        let (before_s, before_t) = self.prefixes[start];
        let (after_s, after_t) = self.prefixes[end];
        let s = (after_s - before_s)
            .checked_add(shift)
            .ok_or_else(intermediate_overflow)?;
        Ok((s, after_t - before_t))
    }

    fn obstruction(&self, model: &mut Model, start: usize, end: usize) -> Result<Cochain> {
        sum_splits(
            model,
            &self.lifts,
            start,
            end,
            cobar::concat_temporary,
            cobar::add_support,
        )
    }
}

pub(crate) fn sum_splits<T>(
    model: &mut Model,
    lifts: &[Vec<Vec<T>>],
    start: usize,
    end: usize,
    product: impl Fn(&mut Model, &[T], &[T]) -> Result<Vec<T>>,
    add: impl Fn(&mut Model, &[T], &[T]) -> Result<Vec<T>>,
) -> Result<Vec<T>> {
    let mut sum = Vec::new();
    for middle in start + 1..end {
        model.spend(1)?;
        let left = &lifts[middle - start - 1][start];
        let right = &lifts[end - middle - 1][middle];
        if left.is_empty() || right.is_empty() {
            continue;
        }
        let term = product(model, left, right)?;
        sum = if sum.is_empty() {
            term
        } else {
            add(model, &sum, &term)?
        };
    }
    Ok(sum)
}

fn check_lifts(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "transfer lift table exceeds the storage limit",
        ));
    }
    Ok(())
}

fn intermediate_overflow() -> Error {
    Error::plain("transfer intermediate cobar bidegree is out of range")
}

// The contraction projects an arbitrary cochain by first closing it with
// the canonical primitive of its differential.
pub(crate) fn project(
    model: &mut Model,
    profile: Profile,
    (s, t): Bidegree,
    source: &[Word],
) -> Result<(Cochain, Vec<usize>)> {
    if source.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    let differential = profile.differential_support(model, source)?;
    let closed = if differential.is_empty() {
        source.to_vec()
    } else {
        let next_s = s
            .checked_add(1)
            .ok_or_else(|| Error::plain("cobar differential degree is out of range"))?;
        let ty = model.profile_ext_type(profile, next_s, t)?;
        let next = model.adams(ty)?;
        let primitive = next
            .primitive(model, &differential)?
            .ok_or_else(|| Error::plain("cobar differential has no boundary preimage"))?;
        cobar::add_support(model, source, &primitive)?
    };
    if closed.is_empty() {
        return Ok((closed, Vec::new()));
    }
    let ty = model.profile_ext_type(profile, s, t)?;
    let coordinates = model.adams(ty)?.project(model, &closed)?;
    Ok((closed, coordinates))
}

pub(crate) fn split_closed(
    model: &mut Model,
    profile: Profile,
    (s, t): Bidegree,
    closed: &[Word],
    coordinates: &[usize],
) -> Result<(Cochain, Cochain)> {
    if closed.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    let ty = model.profile_ext_type(profile, s, t)?;
    let adams = model.adams(ty)?;
    let representative = adams.representative(model, coordinates)?;
    let boundary = cobar::add_support(model, closed, &representative)?;
    let primitive = adams.primitive(model, &boundary)?.ok_or_else(|| {
        Error::plain("cobar cycle differs from its Ext representative by a nonboundary")
    })?;
    Ok((representative, primitive))
}
