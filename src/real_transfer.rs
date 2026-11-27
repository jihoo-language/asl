use std::rc::Rc;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::motivic::Grade;
use crate::real_ext;
use crate::real_motivic::{self, Term};
use crate::transfer::{self, Operation};

pub(crate) type Input = (TypeId, Rc<[usize]>);
type Cochain = Vec<Term>;

#[derive(Clone, Copy, Default)]
struct SuspendedSum {
    s: i128,
    t: i128,
    w: i128,
}

impl SuspendedSum {
    fn push(self, grade: Grade) -> Result<Self> {
        Ok(Self {
            s: self
                .s
                .checked_add(grade.s as i128 - 1)
                .ok_or_else(out_of_range)?,
            t: self.t.checked_add(grade.t).ok_or_else(out_of_range)?,
            w: self
                .w
                .checked_add(i128::from(grade.w))
                .ok_or_else(out_of_range)?,
        })
    }

    fn grade(self, before: Self, shift: i128) -> Result<Option<Grade>> {
        let s = self
            .s
            .checked_sub(before.s)
            .and_then(|s| s.checked_add(shift))
            .ok_or_else(out_of_range)?;
        if s < 0 {
            return Ok(None);
        }
        let s = usize::try_from(s).map_err(|_| out_of_range())?;
        let t = self.t.checked_sub(before.t).ok_or_else(out_of_range)?;
        let w = self.w.checked_sub(before.w).ok_or_else(out_of_range)?;
        let w = i64::try_from(w).map_err(|_| out_of_range())?;
        Ok(Some(Grade { s, t, w }))
    }
}

pub(crate) fn result_type(
    model: &mut Model,
    operation: Operation,
    inputs: impl IntoIterator<Item = TypeId>,
) -> Result<TypeId> {
    let mut sum = SuspendedSum::default();
    for ty in inputs {
        sum = sum.push(read_grade(model, ty, operation)?)?;
    }
    let shift = match operation {
        Operation::Product => 2,
        Operation::Lift => 1,
    };
    let grade = sum.grade(SuspendedSum::default(), shift)?.ok_or_else(|| {
        Error::plain(format!(
            "{} produces a negative cohomological degree",
            operation.name()
        ))
    })?;
    match operation {
        Operation::Product => model.rmext_type(grade),
        Operation::Lift => model.real_cobar_type(grade),
    }
}

fn read_grade(model: &Model, ty: TypeId, operation: Operation) -> Result<Grade> {
    match *model.kind(ty) {
        Kind::RMExt(grade) => Ok(grade),
        _ => Err(Error::plain(format!(
            "{} requires RMExt classes, found {}",
            operation.name(),
            model.type_name(ty)
        ))),
    }
}

fn zero(operation: Operation) -> Value {
    match operation {
        Operation::Product => Value::Class([].into()),
        Operation::Lift => real_motivic::zero(),
    }
}

pub(crate) fn evaluate(
    model: &mut Model,
    operation: Operation,
    output: TypeId,
    inputs: &[Input],
) -> Result<Value> {
    let arity = inputs.len();
    model.spend(arity)?;
    if inputs.iter().any(|(_, coords)| coords.is_empty()) {
        return Ok(zero(operation));
    }
    if let [(ty, coords)] = inputs {
        return match operation {
            Operation::Product => Ok(zero(operation)),
            Operation::Lift => real_ext::representative(model, *ty, coords),
        };
    }
    if let (Operation::Product, [(a, x), (b, y)]) = (operation, inputs) {
        return real_ext::yoneda(model, (*a, x), (*b, y), output);
    }
    let unit = Grade { s: 0, t: 0, w: 0 };
    if inputs
        .iter()
        .any(|(ty, _)| *model.kind(*ty) == Kind::RMExt(unit))
    {
        return Ok(zero(operation));
    }

    let count = arity.saturating_add(1);
    let mut cells = arity.saturating_mul(count) / 2;
    cells = cells.saturating_add(4usize.saturating_mul(count));
    check_storage(cells)?;
    model.spend(cells)?;
    let mut sums = Vec::<SuspendedSum>::with_capacity(count);
    sums.push(SuspendedSum::default());
    let mut singletons = Vec::with_capacity(arity);
    for (ty, coords) in inputs {
        sums.push(
            sums.last()
                .unwrap()
                .push(read_grade(model, *ty, operation)?)?,
        );
        let source = model.real_ext(*ty)?.representative(model, coords)?;
        cells += cochain_cells(&source);
        check_storage(cells)?;
        singletons.push(source);
    }

    let mut lifts = vec![singletons];
    for length in 2..arity {
        let mut row = Vec::with_capacity(arity - length + 1);
        for start in 0..=arity - length {
            let end = start + length;
            let lift = match sums[end].grade(sums[start], 1)? {
                Some(lift_grade) => {
                    let source = obstruction(model, &lifts, start, end)?;
                    homotopy(model, lift_grade.differential()?, source)?
                }
                None => Vec::new(),
            };
            cells += cochain_cells(&lift);
            check_storage(cells)?;
            row.push(lift);
        }
        lifts.push(row);
    }

    let grade = sums[arity].grade(sums[0], 2)?.ok_or_else(out_of_range)?;
    let source = obstruction(model, &lifts, 0, arity)?;
    match operation {
        Operation::Product => {
            let (_, coordinates) = project(model, grade, source)?;
            model.allocate(coordinates.len())?;
            Ok(Value::Class(coordinates.into()))
        }
        Operation::Lift => {
            let lift = homotopy(model, grade, source)?;
            real_motivic::persist(model, lift)
        }
    }
}

fn obstruction(
    model: &mut Model,
    lifts: &[Vec<Cochain>],
    start: usize,
    end: usize,
) -> Result<Cochain> {
    transfer::sum_splits(
        model,
        lifts,
        start,
        end,
        real_motivic::concat_support,
        real_motivic::add,
    )
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "real-motivic transfer lift table exceeds the storage limit",
        ));
    }
    Ok(())
}

fn cochain_cells(terms: &[Term]) -> usize {
    terms.iter().map(Term::cells).sum()
}

fn out_of_range() -> Error {
    Error::plain("real-motivic transfer tridegree exceeds the supported range")
}

// For the linear boundary section q, close x by adding q(dx). The
// remaining non-boundary part projects to the usual RMExt coordinates.
fn project(model: &mut Model, grade: Grade, source: Cochain) -> Result<(Cochain, Vec<usize>)> {
    if source.is_empty() {
        return Ok((source, Vec::new()));
    }
    let differential = real_motivic::differential_support(model, &source)?;
    let closed = if differential.is_empty() {
        source
    } else {
        let primitive = real_ext::linear_primitive(model, grade.differential()?, &differential)?
            .ok_or_else(missing_primitive)?;
        real_motivic::add(model, &source, &primitive)?
    };
    if closed.is_empty() {
        return Ok((closed, Vec::new()));
    }
    let ty = model.rmext_type(grade)?;
    let coordinates = model.real_ext(ty)?.project(model, &closed)?;
    Ok((closed, coordinates))
}

fn homotopy(model: &mut Model, grade: Grade, source: Cochain) -> Result<Cochain> {
    if source.is_empty() || grade.s == 0 {
        return Ok(Vec::new());
    }
    let (closed, coordinates) = project(model, grade, source)?;
    let boundary = if coordinates.is_empty() {
        closed
    } else {
        let ty = model.rmext_type(grade)?;
        let representative = model.real_ext(ty)?.representative(model, &coordinates)?;
        real_motivic::add(model, &closed, &representative)?
    };
    real_ext::linear_primitive(model, grade, &boundary)?.ok_or_else(missing_primitive)
}

fn missing_primitive() -> Error {
    Error::plain("real-motivic contraction boundary has no primitive")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn real_contraction_satisfies_the_sdr_and_side_conditions_on_cobar_columns() {
        exercise_contraction().unwrap_or_else(|error| panic!("{error}"));
    }

    fn exercise_contraction() -> Result<()> {
        let mut model = Model::new();
        for k in 0..=3 {
            for t in -1..=2 * k {
                for s in 0..=3 {
                    let grade = Grade {
                        s,
                        t: i128::from(t),
                        w: t - k,
                    };
                    let ty = model.rmext_type(grade)?;
                    for term in real_motivic::basis(&mut model, grade, 128)? {
                        model.reset_work();
                        let x = vec![term];
                        let (_, p) = project(&mut model, grade, x.clone())?;
                        let ip = model.real_ext(ty)?.representative(&mut model, &p)?;
                        let (_, pi) = project(&mut model, grade, ip.clone())?;
                        assert_eq!(pi, p);
                        assert!(homotopy(&mut model, grade, ip.clone())?.is_empty());

                        let dx = real_motivic::differential_support(&mut model, &x)?;
                        let next = grade.differential()?;
                        let (_, pd) = project(&mut model, next, dx.clone())?;
                        assert!(pd.is_empty());
                        let hd = homotopy(&mut model, next, dx)?;
                        let h = homotopy(&mut model, grade, x.clone())?;
                        let dh = real_motivic::differential_support(&mut model, &h)?;
                        assert_eq!(
                            real_motivic::add(&mut model, &dh, &hd)?,
                            real_motivic::add(&mut model, &x, &ip)?,
                        );
                        if let Some(previous) = grade.previous() {
                            let (_, ph) = project(&mut model, previous, h.clone())?;
                            assert!(ph.is_empty());
                            assert!(homotopy(&mut model, previous, h)?.is_empty());
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
