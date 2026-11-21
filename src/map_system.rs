use std::collections::btree_map::Entry as BTreeEntry;
use std::collections::BTreeMap;

use crate::error::{Error, Result, Span};
use crate::map_cochain::{self, Support};
use crate::map_ext::{self, Term};
use crate::model::{Kind, Model, TypeId, Value};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Class,
    Curvature,
}

impl Operation {
    pub(crate) fn from_name(name: &str) -> Option<Self> {
        match name {
            "eMassey" => Some(Self::Class),
            "eCurvature" => Some(Self::Curvature),
            _ => None,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Class => "eMassey",
            Self::Curvature => "eCurvature",
        }
    }
}

pub(crate) struct Signature {
    operation: Operation,
    vertices: Vec<Vertex>,
}

struct Vertex {
    object: TypeId,
    filtration_sum: i128,
    internal_sum: i128,
}

impl Signature {
    pub(crate) fn new(
        model: &mut Model,
        operation: Operation,
        inputs: impl ExactSizeIterator<Item = (TypeId, Span)> + DoubleEndedIterator,
        span: Span,
    ) -> Result<Self> {
        let arity = inputs.len();
        let name = operation.name();
        if arity < 2 {
            return Err(Error::new(
                span,
                format!("{name} requires at least two adjacent maps"),
            ));
        }
        model.spend(arity)?;
        model.allocate(3 * (arity + 1))?;
        let mut vertices = Vec::<Vertex>::with_capacity(arity + 1);
        for (ty, at) in inputs.rev() {
            let (source, target, s, t) = match *model.kind(ty) {
                Kind::MapExt(source, target, s, t) | Kind::MapCobar(source, target, s, t) => {
                    (source, target, s, t)
                }
                _ => {
                    return Err(Error::new(
                        at,
                        format!(
                            "{name} requires MapExt or MapCobar inputs, found {}",
                            model.type_name(ty)
                        ),
                    ))
                }
            };
            if vertices.is_empty() {
                vertices.push(Vertex {
                    object: source,
                    filtration_sum: 0,
                    internal_sum: 0,
                });
            }
            let previous = vertices.last().unwrap();
            if previous.object != source {
                return Err(Error::new(
                    at,
                    format!("{name} requires composable map endpoints"),
                ));
            }
            vertices.push(Vertex {
                object: target,
                filtration_sum: previous.filtration_sum + s as i128,
                internal_sum: previous.internal_sum + i128::from(t),
            });
        }
        Ok(Self {
            operation,
            vertices,
        })
    }

    pub(crate) fn result_type(&self, model: &mut Model, span: Span) -> Result<TypeId> {
        let arity = self.vertices.len() - 1;
        let (s, t) = self.bidegree(0, arity, 2);
        let name = self.operation.name();
        let s = usize::try_from(s).map_err(|_| {
            let message = if s < 0 {
                format!("{name} produces a negative cohomological degree")
            } else {
                format!("{name} cohomological degree is out of range")
            };
            Error::new(span, message)
        })?;
        let t = i64::try_from(t)
            .map_err(|_| Error::new(span, format!("{name} internal degree is out of range")))?;
        let source = self.vertices[0].object;
        let target = self.vertices[arity].object;
        match self.operation {
            Operation::Class => model.map_ext_type(source, target, s, t),
            Operation::Curvature => model.map_cobar_type(source, target, s, t),
        }
    }

    pub(crate) fn witness_type(
        &self,
        model: &mut Model,
        from: usize,
        to: usize,
        span: Span,
    ) -> Result<TypeId> {
        let arity = self.vertices.len() - 1;
        let interval = || format!("interval [{from},{to}]");
        let error = |message| Error::new(span, message);
        if from >= to || to > arity {
            return Err(error(format!("invalid {}", interval())));
        }
        if to - from == 1 {
            return Err(error(format!(
                "adjacent {} already has an input",
                interval()
            )));
        }
        if to - from == arity {
            return Err(error(format!("top {} cannot have a witness", interval())));
        }
        let (s, t) = self.bidegree(from, to, 1);
        let s = usize::try_from(s).map_err(|_| {
            if s < 0 {
                error(format!(
                    "witness for {} has negative cobar length",
                    interval()
                ))
            } else {
                error(format!(
                    "witness cobar length for {} is out of range",
                    interval()
                ))
            }
        })?;
        let t = i64::try_from(t).map_err(|_| {
            error(format!(
                "witness internal degree for {} is out of range",
                interval()
            ))
        })?;
        model.map_cobar_type(self.vertices[from].object, self.vertices[to].object, s, t)
    }

    fn bidegree(&self, from: usize, to: usize, shift: i128) -> (i128, i128) {
        let start = &self.vertices[from];
        let end = &self.vertices[to];
        (
            end.filtration_sum - start.filtration_sum + shift - (to - from) as i128,
            end.internal_sum - start.internal_sum,
        )
    }
}

pub(crate) struct Entry {
    pub(crate) from: usize,
    pub(crate) to: usize,
    pub(crate) ty: TypeId,
    pub(crate) span: Span,
    pub(crate) value: Value,
}

struct ActiveEntry {
    from: usize,
    to: usize,
    source: TypeId,
    target: TypeId,
    span: Span,
    terms: Support,
}

pub(crate) fn evaluate(
    model: &mut Model,
    arity: usize,
    entries: &[Entry],
    ty: TypeId,
    operation: Operation,
) -> Result<Value> {
    let active = nonzero_entries(model, entries)?;
    let mut top = Vec::new();
    for row in active.chunk_by(|a, b| a.from == b.from) {
        let from = row[0].from;
        let defects = curvature_row(model, row, &active)?;
        for (to, terms) in defects.values {
            if from == 0 && to == arity {
                top = terms;
            } else {
                model.spend(entries.len())?;
                let mut error = Error::plain(format!(
                    "defining system equation fails on interval [{from},{to}]"
                ));
                if let Some(entry) = entries
                    .iter()
                    .find(|entry| (entry.from, entry.to) == (from, to))
                {
                    error = error.at(entry.span);
                }
                return Err(error);
            }
        }
    }
    // The Bianchi identity implies that the top entry is closed once all
    // proper curvature entries vanish.
    match operation {
        Operation::Class => map_cochain::class(model, ty, &top),
        Operation::Curvature => Ok(Value::MapCobar(map_cochain::persist(model, top)?)),
    }
}

fn nonzero_entries(model: &mut Model, entries: &[Entry]) -> Result<Vec<ActiveEntry>> {
    let mut active = Vec::new();
    let mut cells = 0;
    for entry in entries {
        model.spend(1)?;
        let (source, target, terms) = match (model.kind(entry.ty), &entry.value) {
            (&Kind::MapExt(source, target, ..), Value::Class(class)) => {
                let terms = map_cochain::representative(model, entry.ty, class)?;
                (source, target, terms)
            }
            (&Kind::MapCobar(source, target, ..), Value::MapCobar(terms)) => {
                (source, target, terms.clone())
            }
            _ => unreachable!(),
        };
        if terms.is_empty() {
            continue;
        }
        cells += 6 + terms.iter().map(Term::cells).sum::<usize>();
        map_ext::check_storage(cells)?;
        active.push(ActiveEntry {
            from: entry.from,
            to: entry.to,
            source,
            target,
            span: entry.span,
            terms,
        });
    }
    model.spend(cells + active.len() * (1 + active.len().max(1).ilog2() as usize))?;
    active.sort_unstable_by_key(|entry| (entry.from, entry.to));
    Ok(active)
}

fn curvature_row(model: &mut Model, row: &[ActiveEntry], entries: &[ActiveEntry]) -> Result<Row> {
    let mut result = Row::default();
    let search = 1 + entries.len().ilog2() as usize;
    for inner in row {
        let differential =
            map_ext::differential_support(model, inner.source, inner.target, &inner.terms)?;
        if inner.to - inner.from == 1 && !differential.is_empty() {
            return Err(Error::new(
                inner.span,
                format!(
                    "adjacent cochain has nonzero differential on interval [{},{}]",
                    inner.from, inner.to
                ),
            ));
        }
        result.add(model, inner.to, differential)?;

        model.spend(search)?;
        let start = entries.partition_point(|entry| entry.from < inner.to);
        for outer in &entries[start..] {
            model.spend(1)?;
            if outer.from != inner.to {
                break;
            }
            let product = map_ext::compose_terms(model, &outer.terms, &inner.terms)?;
            result.add(model, outer.to, product)?;
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
    fn add(&mut self, model: &mut Model, target: usize, terms: Vec<Term>) -> Result<()> {
        if terms.is_empty() {
            return Ok(());
        }
        model.spend(1 + self.values.len().max(1).ilog2() as usize)?;
        match self.values.entry(target) {
            BTreeEntry::Vacant(entry) => {
                self.cells += 1 + terms.iter().map(Term::cells).sum::<usize>();
                map_ext::check_storage(self.cells)?;
                entry.insert(terms);
            }
            BTreeEntry::Occupied(mut entry) => {
                let previous = entry.get();
                let sum = map_ext::add_terms(model, previous, &terms)?;
                self.cells -= previous.iter().map(Term::cells).sum::<usize>();
                self.cells += sum.iter().map(Term::cells).sum::<usize>();
                map_ext::check_storage(self.cells)?;
                if sum.is_empty() {
                    self.cells -= 1;
                    entry.remove();
                } else {
                    entry.insert(sum);
                }
            }
        }
        Ok(())
    }
}
