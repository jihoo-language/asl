use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::{self, Prime, Term};
use crate::odd_ext;
use crate::prime_linear::Coordinate;
use crate::profile::Profile;
use crate::transfer::{odd_stem, Bidegree, Operation};

pub(crate) type Input = (Bidegree, Rc<[Coordinate]>);
type Cochain = Vec<Term>;

pub(crate) fn evaluate(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    inputs: &[Input],
    operation: Operation,
) -> Result<Value> {
    let arity = inputs.len();
    let (result_s, result_t) = operation.bidegree(inputs.iter().map(|&(grade, _)| grade), arity)?;
    model.spend(arity)?;
    if inputs.iter().any(|(_, class)| class.is_empty()) {
        return Ok(zero(operation));
    }
    if let (Operation::Product, [((s, t), left), ((u, v), right)]) = (operation, inputs) {
        let left_ty = model.profile_odd_ext_type(prime, profile, *s, *t)?;
        let right_ty = model.profile_odd_ext_type(prime, profile, *u, *v)?;
        let ty = model.profile_odd_ext_type(prime, profile, result_s, result_t)?;
        return odd_ext::yoneda(model, left_ty, left, right_ty, right, ty);
    }
    // The only nonzero classes of filtration zero are scalar multiples of
    // the strict unit. Binary multiplication above keeps their coefficient.
    if inputs.iter().any(|((s, _), _)| *s == 0) {
        return Ok(zero(operation));
    }

    let intervals = Intervals::build(model, prime, profile, inputs, arity - 1)?;
    let grade = intervals.grade(0, arity, 2)?;
    let source = intervals.obstruction(model, 0, arity)?;
    let negative = unsuspension_sign(inputs);
    match operation {
        Operation::Product => {
            let (_, mut class) = project(model, prime, profile, grade, &source)?;
            if negative {
                negate(model, prime, &mut class)?;
            }
            odd_ext::persist(model, class)
        }
        Operation::Lift => {
            let mut terms = finish_lift(model, prime, profile, grade, &source)?;
            if negative {
                negate(model, prime, &mut terms)?;
            }
            odd_cobar::persist(model, terms)
        }
    }
}

pub(crate) struct Intervals {
    prime: Prime,
    prefixes: Vec<Bidegree>,
    // Row k stores the suspended lift of each interval of length k + 1.
    lifts: Vec<Vec<Cochain>>,
}

impl Intervals {
    pub(crate) fn build(
        model: &mut Model,
        prime: Prime,
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
        let mut prefixes = Vec::<Bidegree>::with_capacity(prefix_count);
        prefixes.push((0, 0));
        let mut singletons = Vec::with_capacity(arity);
        for &((s, t), ref class) in inputs {
            let &(previous_s, previous_t) = prefixes.last().unwrap();
            let sum_s = previous_s
                .checked_add(s - 1)
                .ok_or_else(intermediate_overflow)?;
            let sum_t = previous_t
                .checked_add(t)
                .ok_or_else(intermediate_overflow)?;
            prefixes.push((sum_s, sum_t));
            let ty = model.profile_odd_ext_type(prime, profile, s, t)?;
            let terms = model.odd_ext(ty)?.representative(model, class)?;
            cells += terms.iter().map(odd_cobar::term_cells).sum::<usize>();
            check_lifts(cells)?;
            singletons.push(terms);
        }
        let mut intervals = Self {
            prime,
            prefixes,
            lifts: vec![singletons],
        };
        for length in 2..=through {
            let mut row = Vec::with_capacity(arity - length + 1);
            for start in 0..=arity - length {
                let end = start + length;
                let grade = intervals.grade(start, end, 2)?;
                let source = intervals.obstruction(model, start, end)?;
                let lift = finish_lift(model, prime, profile, grade, &source)?;
                cells += lift.iter().map(odd_cobar::term_cells).sum::<usize>();
                check_lifts(cells)?;
                row.push(lift);
            }
            intervals.lifts.push(row);
        }
        Ok(intervals)
    }

    pub(crate) fn grade(&self, start: usize, end: usize, shift: usize) -> Result<Bidegree> {
        let (before_s, before_t) = self.prefixes[start];
        let (after_s, after_t) = self.prefixes[end];
        let s = (after_s - before_s)
            .checked_add(shift)
            .ok_or_else(intermediate_overflow)?;
        Ok((s, after_t - before_t))
    }

    pub(crate) fn lift(&self, start: usize, end: usize) -> &[Term] {
        &self.lifts[end - start - 1][start]
    }

    fn obstruction(&self, model: &mut Model, start: usize, end: usize) -> Result<Cochain> {
        let mut sum = Vec::new();
        for middle in start + 1..end {
            model.spend(1)?;
            let left = self.lift(start, middle);
            let right = self.lift(middle, end);
            if left.is_empty() || right.is_empty() {
                continue;
            }
            let mut term = odd_cobar::concat_support(model, self.prime, left, right)?;
            // b_2(sigma a, sigma b) = (-1)^|a| sigma(ab).
            if odd_stem(self.grade(start, middle, 1)?) {
                negate(model, self.prime, &mut term)?;
            }
            sum = if sum.is_empty() {
                term
            } else {
                odd_cobar::add_support(model, self.prime, &sum, &term)?
            };
        }
        Ok(sum)
    }
}

pub(crate) fn project(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    grade: Bidegree,
    source: &[Term],
) -> Result<(Cochain, Vec<Coordinate>)> {
    let image = odd_cobar::differential_profile_support(model, prime, profile, source)?;
    project_with_differential(model, prime, profile, grade, source, &image)
}

fn project_with_differential(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    (s, t): Bidegree,
    source: &[Term],
    image: &[Term],
) -> Result<(Cochain, Vec<Coordinate>)> {
    let closed = if image.is_empty() {
        source.to_vec()
    } else {
        // c - P(dc) lies in the cycles.
        let later = s.checked_add(1).ok_or_else(intermediate_overflow)?;
        let mut primitive = odd_ext::canonical_primitive(model, prime, profile, (later, t), image)?
            .ok_or_else(|| Error::plain("odd-prime cobar differential has no boundary preimage"))?;
        negate(model, prime, &mut primitive)?;
        odd_cobar::add_support(model, prime, source, &primitive)?
    };
    if closed.is_empty() {
        return Ok((closed, Vec::new()));
    }
    let ty = model.profile_odd_ext_type(prime, profile, s, t)?;
    let class = model.odd_ext(ty)?.project(model, &closed)?;
    Ok((closed, class))
}

fn finish_lift(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    grade: Bidegree,
    source: &[Term],
) -> Result<Cochain> {
    let image = odd_cobar::differential_profile_support(model, prime, profile, source)?;
    let exact = if image.is_empty() {
        odd_ext::canonical_primitive(model, prime, profile, grade, source)?
    } else {
        None
    };
    let mut primitive = if let Some(primitive) = exact {
        primitive
    } else {
        let (closed, class) =
            project_with_differential(model, prime, profile, grade, source, &image)?;
        split_closed(model, prime, profile, grade, &closed, &class)?.1
    };
    // The suspended recursion uses -h.
    negate(model, prime, &mut primitive)?;
    Ok(primitive)
}

pub(crate) fn split_closed(
    model: &mut Model,
    prime: Prime,
    profile: Profile,
    (s, t): Bidegree,
    closed: &[Term],
    class: &[Coordinate],
) -> Result<(Cochain, Cochain)> {
    if closed.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    let ty = model.profile_odd_ext_type(prime, profile, s, t)?;
    let group = model.odd_ext(ty)?;
    let representative = group.representative(model, class)?;
    let mut correction = representative.clone();
    negate(model, prime, &mut correction)?;
    let boundary = odd_cobar::add_support(model, prime, closed, &correction)?;
    let primitive = group.primitive(model, &boundary)?.ok_or_else(|| {
        Error::plain("odd-prime cobar cycle differs from its representative by a nonboundary")
    })?;
    Ok((representative, primitive))
}

pub(crate) fn unsuspension_sign(inputs: &[Input]) -> bool {
    // Undo suspension with exponent sum_j (n-j)|x_j|, with j starting at one.
    inputs
        .iter()
        .rev()
        .skip(1)
        .step_by(2)
        .fold(false, |parity, &(grade, _)| parity ^ odd_stem(grade))
}

pub(crate) fn negate<T>(model: &mut Model, prime: Prime, support: &mut [(T, u32)]) -> Result<()> {
    model.spend(support.len())?;
    for (_, coefficient) in support {
        *coefficient = prime.negate(*coefficient);
    }
    Ok(())
}

fn zero(operation: Operation) -> Value {
    match operation {
        Operation::Product => odd_ext::zero(),
        Operation::Lift => odd_cobar::zero(),
    }
}

fn check_lifts(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "odd-prime transfer lift table exceeds the storage limit",
        ));
    }
    Ok(())
}

fn intermediate_overflow() -> Error {
    Error::plain("odd-prime transfer intermediate cobar bidegree is out of range")
}
