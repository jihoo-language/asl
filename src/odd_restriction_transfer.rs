use crate::error::Result;
use crate::model::{Model, Value};
use crate::odd_cobar::{self, Prime, Term};
use crate::odd_ext;
use crate::odd_profile;
use crate::odd_transfer::{self, Input, Intervals};
use crate::prime_linear::{self, Coordinate};
use crate::profile::Profile;
use crate::restriction_transfer::{self, TensorContraction, TensorProjection};
use crate::transfer::Bidegree;

pub(crate) fn evaluate(
    model: &mut Model,
    prime: Prime,
    source: Profile,
    target: Profile,
    inputs: &[Input],
) -> Result<Value> {
    let arity = inputs.len();
    let (s, t) = restriction_transfer::bidegree(inputs.iter().map(|&(degree, _)| degree), arity)?;
    model.spend(arity)?;
    if inputs.iter().any(|(_, coordinates)| coordinates.is_empty()) {
        return Ok(odd_ext::zero());
    }
    if let [(_, class)] = inputs {
        if source == target {
            return Ok(Value::OddClass(class.clone()));
        }
        let from = model.profile_odd_ext_type(prime, source, s, t)?;
        let to = model.profile_odd_ext_type(prime, target, s, t)?;
        return odd_ext::restrict(model, from, class, to);
    }
    if source == target || inputs.iter().any(|((s, _), _)| *s == 0) {
        return Ok(odd_ext::zero());
    }

    restriction_transfer::check_partitions(model, arity)?;
    let Profile::An(n) = target else {
        unreachable!()
    };
    let lifts = Intervals::build(model, prime, source, inputs, arity)?;
    let mut projection = TensorProjection::new(Target {
        prime,
        profile: target,
    });
    let mut coordinates = projection.project_blocks(model, arity, |model, start, end| {
        let terms = odd_profile::restrict(model, prime, n, lifts.lift(start, end))?;
        Ok((lifts.grade(start, end, 1)?, terms))
    })?;
    if odd_transfer::unsuspension_sign(inputs) {
        odd_transfer::negate(model, prime, &mut coordinates)?;
    }
    odd_ext::persist(model, coordinates)
}

struct Target {
    prime: Prime,
    profile: Profile,
}

impl TensorContraction for Target {
    type Term = Term;
    type Coordinate = Coordinate;
    const COORDINATE_CELLS: usize = 2;

    fn cells(term: &Term) -> usize {
        odd_cobar::term_cells(term)
    }

    fn project(
        &self,
        model: &mut Model,
        grade: Bidegree,
        cochain: &[Term],
    ) -> Result<(Vec<Term>, Vec<Coordinate>)> {
        odd_transfer::project(model, self.prime, self.profile, grade, cochain)
    }

    fn split_closed(
        &self,
        model: &mut Model,
        grade: Bidegree,
        closed: &[Term],
        coordinates: &[Coordinate],
    ) -> Result<(Vec<Term>, Vec<Term>)> {
        odd_transfer::split_closed(model, self.prime, self.profile, grade, closed, coordinates)
    }

    fn multiply(&self, model: &mut Model, left: &[Term], right: &[Term]) -> Result<Vec<Term>> {
        odd_cobar::concat_support(model, self.prime, left, right)
    }

    fn add(
        &self,
        model: &mut Model,
        sum: &[Coordinate],
        term: &[Coordinate],
        negative: bool,
    ) -> Result<Vec<Coordinate>> {
        let coefficient = if negative { self.prime.negate(1) } else { 1 };
        prime_linear::add_scaled(model, self.prime, sum, term, coefficient)
    }
}
