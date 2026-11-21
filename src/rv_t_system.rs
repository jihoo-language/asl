use std::collections::BTreeMap;
use std::rc::Rc;

use crate::cohomology;
use crate::error::{Error, Result, Span};
use crate::model::{Kind, Model, TypeId, Value};
use crate::rv_ideal::Ideal;
use crate::rv_koszul::{self, Grade, Support, Term};
use crate::syntax::validate_witness_interval;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Class,
    Curvature,
}

impl Operation {
    pub(crate) fn from_name(name: &str) -> Option<Self> {
        match name {
            "rvTMassey" => Some(Self::Class),
            "rvTCurvature" => Some(Self::Curvature),
            _ => None,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Class => "rvTMassey",
            Self::Curvature => "rvTCurvature",
        }
    }
}

pub(crate) struct Signature {
    operation: Operation,
    ideal: Rc<str>,
    prefixes: Vec<(u128, u128)>,
}

impl Signature {
    pub(crate) fn new(
        model: &mut Model,
        operation: Operation,
        inputs: impl ExactSizeIterator<Item = (TypeId, Span)>,
        span: Span,
    ) -> Result<Self> {
        let name = operation.name();
        let arity = inputs.len();
        if arity < 2 {
            return Err(Error::new(
                span,
                format!("{name} requires at least two adjacent inputs"),
            ));
        }
        model.spend(arity)?;
        model.allocate(4 * (arity + 1))?;
        let mut prefixes = Vec::with_capacity(arity + 1);
        prefixes.push((0, 0));
        let mut ideal = None;
        for (ty, location) in inputs {
            let grade = match model.kind(ty) {
                Kind::RVKoszul(grade) | Kind::RVTor(grade) => grade,
                _ => {
                    return Err(Error::new(
                        location,
                        format!(
                            "{name} requires RVTor or RVKoszul inputs, found {}",
                            model.type_name(ty)
                        ),
                    ))
                }
            };
            if ideal.get_or_insert_with(|| grade.ideal.clone()) != &grade.ideal {
                return Err(Error::new(
                    location,
                    format!("{name} requires the same rank ideal"),
                ));
            }
            let &(s, t) = prefixes.last().unwrap();
            prefixes.push((s + grade.s as u128, t + grade.t as u128));
        }
        Ok(Self {
            operation,
            ideal: ideal.unwrap(),
            prefixes,
        })
    }

    pub(crate) fn result_type(&self, model: &mut Model, span: Span) -> Result<TypeId> {
        let arity = self.prefixes.len() - 1;
        let grade = self.grade(0, arity, 2, span)?;
        model.rv_complex_type(grade, matches!(self.operation, Operation::Class))
    }

    pub(crate) fn witness_type(
        &self,
        model: &mut Model,
        from: usize,
        to: usize,
        span: Span,
    ) -> Result<TypeId> {
        validate_witness_interval(self.prefixes.len() - 1, from, to, span)?;
        model.rv_complex_type(self.grade(from, to, 1, span)?, false)
    }

    fn grade(&self, from: usize, to: usize, shift: usize, span: Span) -> Result<Grade> {
        let (s0, t0) = self.prefixes[from];
        let (s1, t1) = self.prefixes[to];
        let overflow = || {
            Error::new(
                span,
                format!(
                    "{} bidegree for interval [{from},{to}] is out of range",
                    self.operation.name()
                ),
            )
        };
        let s = usize::try_from(s1 - s0 + (to - from - shift) as u128).map_err(|_| overflow())?;
        let t = usize::try_from(t1 - t0).map_err(|_| overflow())?;
        Ok(Grade {
            ideal: self.ideal.clone(),
            s,
            t,
        })
    }
}

pub(crate) struct Entry {
    pub(crate) from: usize,
    pub(crate) to: usize,
    pub(crate) ty: TypeId,
    pub(crate) span: Span,
    pub(crate) value: Value,
}

struct Active {
    from: usize,
    to: usize,
    span: Span,
    terms: Support,
}

pub(crate) fn evaluate(
    model: &mut Model,
    arity: usize,
    operation: Operation,
    ty: TypeId,
    entries: &[Entry],
) -> Result<Value> {
    let (Kind::RVKoszul(grade) | Kind::RVTor(grade)) = model.kind(ty) else {
        unreachable!();
    };
    let ideal = model.rv_ideal_named(&grade.ideal).unwrap();
    let active = active_entries(model, entries)?;
    let mut top = Vec::new();
    for row in active.chunk_by(|a, b| a.from == b.from) {
        let from = row[0].from;
        let defects = curvature_row(model, &ideal, row, &active)?;
        for (to, terms) in defects.values {
            if from == 0 && to == arity {
                top = terms;
            } else {
                let error = Error::plain(format!(
                    "defining system equation fails on interval [{from},{to}]"
                ));
                model.spend(entries.len())?;
                return Err(
                    match entries
                        .iter()
                        .find(|entry| (entry.from, entry.to) == (from, to))
                    {
                        Some(entry) => error.at(entry.span),
                        None => error,
                    },
                );
            }
        }
    }
    // Once every other entry of dA + A^2 vanishes, its top entry is a cycle.
    match operation {
        Operation::Curvature => Ok(Value::RVKoszul(rv_koszul::persist(model, top)?)),
        Operation::Class => {
            if top.is_empty() {
                return Ok(cohomology::zero());
            }
            let class = model.rv_tor(ty)?.project(model, &ideal, &top)?;
            model.allocate(class.len())?;
            Ok(Value::Class(class.into()))
        }
    }
}

fn active_entries(model: &mut Model, entries: &[Entry]) -> Result<Vec<Active>> {
    let mut result = Vec::new();
    let mut cells = 0;
    for entry in entries {
        model.spend(1)?;
        let terms = match &entry.value {
            Value::Class(class) if class.is_empty() => continue,
            Value::Class(class) => model
                .rv_tor(entry.ty)?
                .representative_terms(model, class)?
                .into(),
            Value::RVKoszul(terms) => terms.clone(),
            _ => unreachable!(),
        };
        if terms.is_empty() {
            continue;
        }
        cells += 4 + rv_koszul::cells(&terms);
        rv_koszul::check_storage(cells)?;
        result.push(Active {
            from: entry.from,
            to: entry.to,
            span: entry.span,
            terms,
        });
    }
    model.spend(result.len() * (1 + result.len().max(1).ilog2() as usize))?;
    result.sort_unstable_by_key(|entry| (entry.from, entry.to));
    Ok(result)
}

fn curvature_row(
    model: &mut Model,
    ideal: &Ideal,
    row: &[Active],
    entries: &[Active],
) -> Result<Row> {
    let mut result = Row::default();
    for left in row {
        let differential = rv_koszul::differential_terms(model, ideal, &left.terms)?;
        if left.to - left.from == 1 && !differential.is_empty() {
            return Err(Error::new(
                left.span,
                format!(
                    "adjacent cochain has nonzero differential on interval [{},{}]",
                    left.from, left.to
                ),
            ));
        }
        result.add(model, left.to, differential)?;

        model.spend(1 + entries.len().ilog2() as usize)?;
        let start = entries.partition_point(|next| next.from < left.to);
        for right in &entries[start..] {
            model.spend(1)?;
            if right.from != left.to {
                break;
            }
            let product = rv_koszul::product_terms(model, ideal, &left.terms, &right.terms)?;
            result.add(model, right.to, product)?;
        }
    }
    Ok(result)
}

#[derive(Default)]
struct Row {
    values: BTreeMap<usize, Vec<Term>>,
    cells: usize,
}

impl Row {
    fn add(&mut self, model: &mut Model, to: usize, mut terms: Vec<Term>) -> Result<()> {
        if terms.is_empty() {
            return Ok(());
        }
        model.spend(1 + self.values.len().max(1).ilog2() as usize)?;
        if let Some(previous) = self.values.remove(&to) {
            self.cells -= rv_koszul::cells(&previous);
            terms = rv_koszul::sum_terms(model, &previous, &terms)?;
        }
        if !terms.is_empty() {
            self.cells += rv_koszul::cells(&terms);
            rv_koszul::check_storage(self.cells)?;
            self.values.insert(to, terms);
        }
        Ok(())
    }
}
