use std::collections::BTreeMap;
use std::rc::Rc;

use crate::cohomology;
use crate::error::{Error, Result};
use crate::linear::Basis;
use crate::model::{Kind, Model, TypeId, Value};
use crate::rv_ideal::Ideal;
use crate::rv_koszul::{self, Grade, Term};

type Input = (TypeId, Rc<[usize]>);
type Chain = Vec<Term>;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Higher,
    Lift,
}

impl Operation {
    pub(crate) fn from_name(name: &str) -> Option<Self> {
        match name {
            "rvTHigher" => Some(Self::Higher),
            "rvTLift" => Some(Self::Lift),
            _ => None,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Higher => "rvTHigher",
            Self::Lift => "rvTLift",
        }
    }

    pub(crate) fn require_arity(self, arity: usize) -> Result<()> {
        match self {
            Self::Higher if arity < 2 => Err(Error::plain(
                "rvTHigher requires at least two RVTor classes",
            )),
            Self::Lift if arity == 0 => {
                Err(Error::plain("rvTLift requires at least one RVTor class"))
            }
            _ => Ok(()),
        }
    }

    fn shift(self, arity: usize) -> usize {
        match self {
            Self::Higher => arity - 2,
            Self::Lift => arity - 1,
        }
    }

    fn zero(self) -> Value {
        match self {
            Self::Higher => cohomology::zero(),
            Self::Lift => Value::RVKoszul([].into()),
        }
    }
}

pub(crate) fn result_type(
    model: &mut Model,
    operation: Operation,
    inputs: impl ExactSizeIterator<Item = TypeId>,
) -> Result<TypeId> {
    let arity = inputs.len();
    operation.require_arity(arity)?;
    let name = operation.name();
    let mut ideal = None;
    let mut s = operation.shift(arity) as u128;
    let mut t = 0usize;
    for input in inputs {
        let Kind::RVTor(grade) = model.kind(input) else {
            return Err(Error::plain(format!(
                "{name} requires RVTor classes, found {}",
                model.type_name(input)
            )));
        };
        if grade.multidegree.is_some() {
            return Err(Error::plain(format!("{name} requires total RVTor classes")));
        }
        if ideal.get_or_insert_with(|| grade.ideal.clone()) != &grade.ideal {
            return Err(Error::plain(format!("{name} requires the same rank ideal")));
        }
        s += grade.s as u128;
        t = t.checked_add(grade.t).ok_or_else(degree_limit)?;
    }
    let s = usize::try_from(s).map_err(|_| degree_limit())?;
    let grade = Grade::total(ideal.unwrap(), s, t);
    model.rv_complex_type(grade, matches!(operation, Operation::Higher))
}

fn degree_limit() -> Error {
    Error::plain("Koszul transfer bidegree is out of range")
}

pub(crate) fn evaluate(
    model: &mut Model,
    operation: Operation,
    target: TypeId,
    inputs: &[Input],
) -> Result<Value> {
    model.spend(inputs.len())?;
    if inputs.iter().any(|(_, class)| class.is_empty()) {
        return Ok(operation.zero());
    }
    if let [(ty, class)] = inputs {
        return Ok(Value::RVKoszul(
            model.rv_tor(*ty)?.representative(model, class)?,
        ));
    }

    let (Kind::RVTor(grade) | Kind::RVKoszul(grade)) = model.kind(target) else {
        unreachable!();
    };
    let name = grade.ideal.clone();
    let ideal = model.rv_ideal_named(&name).unwrap();
    let contraction = Contraction::new(model, &ideal)?;
    let intervals = Intervals::build(model, &name, &contraction, inputs)?;
    let source = intervals.obstruction(model, &ideal, 0, inputs.len())?;
    let grade = intervals.obstruction_grade(&name, 0, inputs.len())?;
    match operation {
        Operation::Higher => {
            let class = contraction.project(model, &grade, source)?;
            model.allocate(class.len())?;
            Ok(Value::Class(class.into()))
        }
        Operation::Lift => {
            let lift = contraction.homotopy(model, &grade, source)?;
            Ok(Value::RVKoszul(rv_koszul::persist(model, lift)?))
        }
    }
}

struct Intervals {
    prefixes: Vec<(u128, usize)>,
    lifts: Vec<Vec<Chain>>,
}

impl Intervals {
    fn build(
        model: &mut Model,
        name: &Rc<str>,
        contraction: &Contraction<'_>,
        inputs: &[Input],
    ) -> Result<Self> {
        let arity = inputs.len();
        let interval_count = arity.saturating_mul(arity.saturating_add(1)) / 2;
        model.spend(interval_count)?;
        let mut cells = interval_count.saturating_add(3 * (arity + 1));
        rv_koszul::check_storage(cells)?;

        let mut prefixes = vec![(0, 0usize)];
        let mut singletons = Vec::with_capacity(arity);
        for (ty, class) in inputs {
            let Kind::RVTor(grade) = model.kind(*ty) else {
                unreachable!();
            };
            let &(s, t) = prefixes.last().unwrap();
            prefixes.push((
                s + grade.s as u128,
                t.checked_add(grade.t).ok_or_else(degree_limit)?,
            ));
            let chain = model.rv_tor(*ty)?.representative_terms(model, class)?;
            cells += rv_koszul::cells(&chain);
            rv_koszul::check_storage(cells)?;
            singletons.push(chain);
        }

        let mut intervals = Self {
            prefixes,
            lifts: vec![singletons],
        };
        for length in 2..arity {
            let mut row = Vec::with_capacity(arity - length + 1);
            for start in 0..=arity - length {
                let end = start + length;
                let source = intervals.obstruction(model, contraction.ideal, start, end)?;
                let grade = intervals.obstruction_grade(name, start, end)?;
                let lift = contraction.homotopy(model, &grade, source)?;
                cells += rv_koszul::cells(&lift);
                rv_koszul::check_storage(cells)?;
                row.push(lift);
            }
            intervals.lifts.push(row);
        }
        Ok(intervals)
    }

    fn obstruction_grade(&self, name: &Rc<str>, start: usize, end: usize) -> Result<Grade> {
        let (before_s, before_t) = self.prefixes[start];
        let (after_s, after_t) = self.prefixes[end];
        let s = usize::try_from(after_s - before_s + (end - start - 2) as u128)
            .map_err(|_| degree_limit())?;
        Ok(Grade::total(name.clone(), s, after_t - before_t))
    }

    fn obstruction(
        &self,
        model: &mut Model,
        ideal: &Ideal,
        start: usize,
        end: usize,
    ) -> Result<Chain> {
        let mut source = Vec::new();
        for split in start + 1..end {
            model.spend(1)?;
            let left = &self.lifts[split - start - 1][start];
            let right = &self.lifts[end - split - 1][split];
            if left.is_empty() || right.is_empty() {
                continue;
            }
            let product = rv_koszul::product_terms(model, ideal, left, right)?;
            source = if source.is_empty() {
                product
            } else {
                rv_koszul::sum_terms(model, &source, &product)?
            };
            rv_koszul::check_storage(rv_koszul::cells(&source))?;
        }
        Ok(source)
    }
}

struct Contraction<'a> {
    ideal: &'a Ideal,
    initial: Option<Vec<&'a [usize]>>,
}

impl<'a> Contraction<'a> {
    fn new(model: &mut Model, ideal: &'a Ideal) -> Result<Self> {
        let mut initial = Vec::new();
        for polynomial in ideal.polynomials() {
            model.spend(1)?;
            if polynomial.len() != 1 {
                return Ok(Self {
                    ideal,
                    initial: None,
                });
            }
            initial.push(polynomial[0].as_ref());
        }
        Ok(Self {
            ideal,
            initial: (!initial.is_empty()).then_some(initial),
        })
    }

    fn close(&self, model: &mut Model, grade: &Grade, source: Chain) -> Result<Chain> {
        if source.is_empty() {
            return Ok(source);
        }
        // The binary sum need not be a cycle. Remove its component
        // complementary to ker d using the same section as the homotopy.
        let differential = rv_koszul::differential_terms(model, self.ideal, &source)?;
        if differential.is_empty() {
            return Ok(source);
        }
        let previous = grade.differential()?;
        let correction = self
            .primitive(model, &previous, &differential)?
            .ok_or_else(|| Error::plain("Koszul differential has no boundary preimage"))?;
        let closed = rv_koszul::sum_terms(model, &source, &correction)?;
        rv_koszul::check_storage(rv_koszul::cells(&closed))?;
        Ok(closed)
    }

    fn project(&self, model: &mut Model, grade: &Grade, source: Chain) -> Result<Vec<usize>> {
        let closed = self.close(model, grade, source)?;
        if closed.is_empty() {
            return Ok(Vec::new());
        }
        if let Some(initial) = &self.initial {
            if monomial_primitive(model, initial, grade, &closed)?.is_some() {
                return Ok(Vec::new());
            }
        }
        let ty = model.rv_complex_type(grade.clone(), true)?;
        model.rv_tor(ty)?.project(model, self.ideal, &closed)
    }

    fn homotopy(&self, model: &mut Model, grade: &Grade, source: Chain) -> Result<Chain> {
        let closed = self.close(model, grade, source)?;
        if closed.is_empty() {
            return Ok(closed);
        }
        if let Some(initial) = &self.initial {
            if let Some(primitive) = monomial_primitive(model, initial, grade, &closed)? {
                return Ok(primitive);
            }
        }
        let ty = model.rv_complex_type(grade.clone(), true)?;
        let group = model.rv_tor(ty)?;
        let class = group.project(model, self.ideal, &closed)?;
        let representative = group.representative_terms(model, &class)?;
        let boundary = rv_koszul::sum_terms(model, &closed, &representative)?;
        self.primitive(model, grade, &boundary)?.ok_or_else(|| {
            Error::plain("Koszul cycle differs from its Tor representative by a nonboundary")
        })
    }

    fn primitive(
        &self,
        model: &mut Model,
        grade: &Grade,
        boundary: &[Term],
    ) -> Result<Option<Chain>> {
        if boundary.is_empty() {
            return Ok(Some(Vec::new()));
        }
        if let Some(initial) = &self.initial {
            return monomial_primitive(model, initial, grade, boundary);
        }
        let ty = model.rv_complex_type(grade.clone(), true)?;
        model
            .rv_tor(ty)?
            .primitive(model, self.ideal, grade, boundary)
    }
}

// Reduction at a monomial ideal preserves coefficient exponents plus wedge
// indicators. Restricting ordered incoming columns to an individual
// multidegree gives exactly the same linear section as the full reduction.
fn monomial_primitive(
    model: &mut Model,
    initial: &[&[usize]],
    grade: &Grade,
    support: &[Term],
) -> Result<Option<Chain>> {
    let width = initial[0].len();
    let source_degree = grade.primitive()?.s;
    let mut blocks = BTreeMap::<Vec<usize>, Vec<usize>>::new();
    for term in support {
        model.spend((1 + width) * (1 + blocks.len().max(1).ilog2() as usize))?;
        let mut multidegree = term.powers.to_vec();
        for (index, exponent) in multidegree.iter_mut().enumerate() {
            if term.wedge & (1usize << index) != 0 {
                *exponent = exponent.checked_add(1).ok_or_else(degree_limit)?;
            }
        }
        blocks.entry(multidegree).or_default().push(term.wedge);
    }
    let mut primitive = Vec::new();
    for (multidegree, column) in blocks {
        let active: Vec<_> = multidegree
            .iter()
            .enumerate()
            .filter_map(|(index, &e)| (e != 0).then_some(1usize << index))
            .collect();
        let mut images = Basis::default();
        visit_wedges(model, &active, source_degree, 0, &mut |model, wedge| {
            let mut powers = remove_wedge(&multidegree, wedge);
            if !standard(model, initial, &powers)? {
                return Ok(());
            }
            let mut image = Vec::new();
            let mut remaining = wedge;
            // Removing the highest index gives the smallest resulting mask.
            while remaining != 0 {
                let index = remaining.ilog2() as usize;
                let bit = 1usize << index;
                remaining ^= bit;
                powers[index] += 1;
                if standard(model, initial, &powers)? {
                    image.push(wedge ^ bit);
                }
                powers[index] -= 1;
            }
            images.insert_with_witness(model, image, vec![wedge])?;
            Ok(())
        })?;
        let (remainder, witness) = images.reduce_with_witness(model, column, Vec::new())?;
        if !remainder.is_empty() {
            return Ok(None);
        }
        let total = primitive.len() + witness.len();
        rv_koszul::check_storage(1 + (2 + width) * total)?;
        model.spend((1 + width) * witness.len())?;
        primitive.extend(witness.into_iter().map(|wedge| Term {
            wedge,
            powers: remove_wedge(&multidegree, wedge).into(),
        }));
    }
    model.spend((1 + width) * primitive.len() * (1 + primitive.len().max(1).ilog2() as usize))?;
    primitive.sort_unstable();
    Ok(Some(primitive))
}

fn remove_wedge(multidegree: &[usize], wedge: usize) -> Vec<usize> {
    multidegree
        .iter()
        .enumerate()
        .map(|(index, &e)| e - ((wedge >> index) & 1))
        .collect()
}

fn standard(model: &mut Model, initial: &[&[usize]], powers: &[usize]) -> Result<bool> {
    for &leading in initial {
        model.spend(1 + powers.len())?;
        if leading.iter().zip(powers).all(|(&a, &b)| a <= b) {
            return Ok(false);
        }
    }
    Ok(true)
}

fn visit_wedges(
    model: &mut Model,
    active: &[usize],
    remaining: usize,
    wedge: usize,
    visit: &mut impl FnMut(&mut Model, usize) -> Result<()>,
) -> Result<()> {
    model.spend(1)?;
    if remaining == 0 {
        visit(model, wedge)
    } else {
        for index in remaining - 1..active.len() {
            visit_wedges(
                model,
                &active[..index],
                remaining - 1,
                wedge | active[index],
                visit,
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::q_koszul;
    use crate::rv_koszul::Slice;
    use crate::rv_poly;

    #[test]
    fn contraction_identity_and_side_conditions_hold_on_arbitrary_chains() {
        exercise_contraction().unwrap_or_else(|error| panic!("{error}"));
    }

    fn exercise_contraction() -> Result<()> {
        let mut model = Model::new();
        let indices = q_koszul::indices(&mut model, vec![0, 1, 2])?;
        let xy = rv_poly::monomial(vec![1, 1, 0]);
        let z2 = rv_poly::monomial(vec![0, 0, 2]);
        let quadric = rv_poly::add(&mut model, &xy, &z2)?;
        for (name, generators) in [
            ("Polynomial", Vec::new()),
            ("Monomial", vec![(2, xy), (2, z2)]),
            ("Quadric", vec![(2, quadric)]),
        ] {
            let ideal = Ideal::build(&mut model, indices.clone(), generators)?;
            model.define_rv_ideal(name.to_owned(), ideal);
            let ideal = model.rv_ideal_named(name).unwrap();
            let contraction = Contraction::new(&mut model, &ideal)?;
            for s in 0..=indices.len() {
                for t in s..=s + 3 {
                    let grade = Grade::total(name.into(), s, t);
                    let slice = Slice::build(&mut model, &ideal, s, t, 128)?;
                    let ty = model.rv_complex_type(grade.clone(), true)?;
                    let tor = model.rv_tor(ty)?;
                    let next = grade.primitive()?;
                    for index in 0..slice.len() {
                        model.reset_work();
                        let w = vec![slice.term(index)];
                        let p = contraction.project(&mut model, &grade, w.clone())?;
                        let ip = tor.representative_terms(&mut model, &p)?;
                        let h = contraction.homotopy(&mut model, &grade, w.clone())?;
                        let dw = rv_koszul::differential_terms(&mut model, &ideal, &w)?;
                        let dh = rv_koszul::differential_terms(&mut model, &ideal, &h)?;
                        let hd = if dw.is_empty() {
                            Vec::new()
                        } else {
                            let previous = grade.differential()?;
                            if let Some(initial) = &contraction.initial {
                                let ty = model.rv_complex_type(previous.clone(), true)?;
                                let group = model.rv_tor(ty)?;
                                let canonical =
                                    group.primitive(&mut model, &ideal, &previous, &dw)?;
                                let multigraded =
                                    monomial_primitive(&mut model, initial, &previous, &dw)?;
                                assert!(
                                    canonical == multigraded,
                                    "{name} boundary section at ({s},{t})"
                                );
                            }
                            contraction.homotopy(&mut model, &previous, dw)?
                        };
                        assert!(
                            rv_koszul::sum_terms(&mut model, &dh, &hd)?
                                == rv_koszul::sum_terms(&mut model, &w, &ip)?,
                            "{name} contraction identity at ({s},{t}), term {index}",
                        );
                        assert!(contraction
                            .homotopy(&mut model, &next, h.clone())?
                            .is_empty());
                        assert!(contraction.project(&mut model, &next, h)?.is_empty());
                        assert!(contraction
                            .homotopy(&mut model, &grade, ip.clone())?
                            .is_empty());
                        assert_eq!(contraction.project(&mut model, &grade, ip)?, p);
                    }
                }
            }
        }
        Ok(())
    }
}
