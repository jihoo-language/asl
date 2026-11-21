use std::collections::{btree_map::Entry as BTreeEntry, BTreeMap};

use crate::error::{Error, Result, Span};
use crate::model::{Kind, Model, TypeId, Value};
use crate::odd_cobar::{self, Prime, Support, Term};
use crate::odd_ext;
use crate::syntax::validate_witness_interval;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Class,
    Curvature,
    Power,
}

impl Operation {
    pub(crate) fn from_name(name: &str) -> Option<Self> {
        match name {
            "pMassey" => Some(Self::Class),
            "pCurvature" => Some(Self::Curvature),
            "pPowerDef" => Some(Self::Power),
            _ => None,
        }
    }

    pub(crate) fn explicit(self) -> bool {
        !matches!(self, Self::Power)
    }

    fn name(self) -> &'static str {
        match self {
            Self::Class => "pMassey",
            Self::Curvature => "pCurvature",
            Self::Power => "pPowerDef",
        }
    }
}

pub(crate) struct Signature {
    operation: Operation,
    arity: usize,
    prime: Prime,
    prefix: Vec<(i128, i128)>,
}

impl Signature {
    pub(crate) fn new(
        model: &mut Model,
        operation: Operation,
        inputs: impl ExactSizeIterator<Item = (TypeId, Span)>,
        span: Span,
    ) -> Result<Self> {
        let count = inputs.len();
        let name = operation.name();
        let arity = if operation.explicit() {
            if count < 2 {
                return Err(Error::new(
                    span,
                    format!("{name} requires at least two inputs"),
                ));
            }
            count
        } else {
            if count == 0 {
                return Err(Error::new(span, "pPowerDef requires at least one input"));
            }
            count + 1
        };
        let adjacent = if operation.explicit() { count } else { 1 };
        model.spend(count)?;
        let mut prime = None;
        let mut prefix = Vec::<(i128, i128)>::with_capacity(adjacent + 1);
        prefix.push((0, 0));
        for (ty, at) in inputs.take(adjacent) {
            let (p, s, t) = match *model.kind(ty) {
                Kind::PCobar(p, s, t) | Kind::PExt(p, s, t) => (p, s, t),
                _ => {
                    return Err(Error::new(
                        at,
                        format!(
                            "{name} requires PCobar or PExt inputs, found {}",
                            model.type_name(ty)
                        ),
                    ));
                }
            };
            if *prime.get_or_insert(p) != p {
                return Err(Error::new(
                    at,
                    format!("{name} requires a common odd prime"),
                ));
            }
            let &(length, degree) = prefix.last().unwrap();
            prefix.push((length + s as i128, degree + t as i128));
        }
        Ok(Self {
            operation,
            arity,
            prime: prime.unwrap(),
            prefix,
        })
    }

    pub(crate) fn arity(&self) -> usize {
        self.arity
    }

    pub(crate) fn result_type(&self, model: &mut Model, span: Span) -> Result<TypeId> {
        let (s, t) = self.bidegree(0, self.arity, 2, self.operation.name(), span)?;
        match self.operation {
            Operation::Curvature => model.pcobar_type(self.prime, s, t),
            _ => model.odd_ext_type(self.prime, s, t),
        }
    }

    pub(crate) fn witness_type(
        &self,
        model: &mut Model,
        from: usize,
        to: usize,
        span: Span,
    ) -> Result<TypeId> {
        validate_witness_interval(self.arity, from, to, span)?;
        let context = format!("witness for interval [{from},{to}]");
        let (s, t) = self.bidegree(from, to, 1, &context, span)?;
        model.pcobar_type(self.prime, s, t)
    }

    fn bidegree(
        &self,
        from: usize,
        to: usize,
        shift: i128,
        context: &str,
        span: Span,
    ) -> Result<(usize, usize)> {
        let gap = (to - from) as i128;
        let (length, degree) = if self.operation.explicit() {
            let (start_s, start_t) = self.prefix[from];
            let (end_s, end_t) = self.prefix[to];
            (end_s - start_s, end_t - start_t)
        } else {
            let (s, t) = self.prefix[1];
            (s * gap, t * gap)
        };
        let length = length + shift - gap;
        let s = usize::try_from(length).map_err(|_| {
            if length < 0 {
                Error::new(span, format!("{context} has negative cobar length"))
            } else {
                Error::new(span, format!("{context} cobar length is out of range"))
            }
        })?;
        let t = usize::try_from(degree)
            .map_err(|_| Error::new(span, format!("{context} internal degree is out of range")))?;
        Ok((s, t))
    }
}

pub(crate) struct Entry {
    pub(crate) from: usize,
    pub(crate) to: usize,
    pub(crate) ty: TypeId,
    pub(crate) span: Span,
    pub(crate) value: Value,
}

struct Lift {
    from: usize,
    to: usize,
    span: Span,
    negative_bar: bool,
    terms: Support,
}

pub(crate) fn evaluate(
    model: &mut Model,
    arity: usize,
    entries: &[Entry],
    ty: TypeId,
    operation: Operation,
) -> Result<Value> {
    let (Kind::PCobar(prime, ..) | Kind::PExt(prime, ..)) = *model.kind(ty) else {
        unreachable!()
    };
    let mut lifts = Vec::with_capacity(entries.len());
    let mut cells = 0;
    for entry in entries {
        model.spend(1)?;
        let (s, t, terms) = match (model.kind(entry.ty), &entry.value) {
            (&Kind::PCobar(_, s, t), Value::PCobar(terms)) => (s, t, terms.clone()),
            (&Kind::PExt(_, s, t), Value::OddClass(class)) => {
                let terms = if class.is_empty() {
                    [].into()
                } else {
                    model
                        .odd_ext(entry.ty)?
                        .representative(model, class)?
                        .into()
                };
                (s, t, terms)
            }
            _ => unreachable!(),
        };
        if terms.is_empty() && operation.explicit() {
            continue;
        }
        let term_cells = support_cells(&terms);
        cells += 5 + term_cells;
        odd_cobar::check_storage(cells)?;
        model.spend(term_cells)?;
        lifts.push(Lift {
            from: entry.from,
            to: entry.to,
            span: entry.span,
            negative_bar: t % 2 == s % 2,
            terms,
        });
    }
    let curvature = if operation.explicit() {
        explicit_curvature(model, prime, arity, entries, &mut lifts)?
    } else {
        power_curvature(model, prime, arity, &lifts)?
    };
    // The proper equations imply that the top curvature is closed by the
    // graded Leibniz rule, with no input needed at (0, arity).
    match operation {
        Operation::Curvature => odd_cobar::persist(model, curvature),
        _ => odd_ext::class(model, ty, &curvature),
    }
}

fn explicit_curvature(
    model: &mut Model,
    prime: Prime,
    arity: usize,
    entries: &[Entry],
    lifts: &mut [Lift],
) -> Result<Vec<Term>> {
    let search = 1 + lifts.len().max(1).ilog2() as usize;
    model.spend(lifts.len() * search)?;
    lifts.sort_unstable_by_key(|lift| (lift.from, lift.to));
    let mut top = Vec::new();
    for row in lifts.chunk_by(|a, b| a.from == b.from) {
        let from = row[0].from;
        let mut defects = Row::default();
        for a in row {
            let differential = odd_cobar::differential_support(model, prime, &a.terms)?;
            if a.to - a.from == 1 && !differential.is_empty() {
                return Err(adjacent_error(a));
            }
            // Store the right-hand side minus the supplied differential.
            let negative_differential = signed(model, prime, differential, true)?;
            defects.add(model, prime, a.to, negative_differential)?;

            model.spend(search)?;
            let start = lifts.partition_point(|b| b.from < a.to);
            for b in &lifts[start..] {
                model.spend(1)?;
                if b.from != a.to {
                    break;
                }
                let product = bar_product(model, prime, a, b)?;
                defects.add(model, prime, b.to, product)?;
            }
        }
        for (to, terms) in defects.values {
            if from == 0 && to == arity {
                top = terms;
            } else {
                model.spend(entries.len())?;
                let mut error = equation_error(from, to);
                if let Some(entry) = entries.iter().find(|e| (e.from, e.to) == (from, to)) {
                    error = error.at(entry.span);
                }
                return Err(error);
            }
        }
    }
    Ok(top)
}

fn power_curvature(
    model: &mut Model,
    prime: Prime,
    arity: usize,
    lifts: &[Lift],
) -> Result<Vec<Term>> {
    model.spend(lifts.len())?;
    let active: Vec<_> = lifts.iter().filter(|lift| !lift.terms.is_empty()).collect();
    let mut curvature = Vec::new();
    for gap in 1..=arity {
        model.spend(1)?;
        let mut product = Vec::new();
        for a in &active {
            model.spend(1)?;
            if a.to >= gap {
                break;
            }
            let b = &lifts[gap - a.to - 1];
            if !b.terms.is_empty() {
                let term = bar_product(model, prime, a, b)?;
                product = odd_cobar::add_support(model, prime, &product, &term)?;
            }
        }
        if let Some(lift) = lifts.get(gap - 1) {
            let differential = odd_cobar::differential_support(model, prime, &lift.terms)?;
            if odd_cobar::compare(model, &differential, &product)?.is_ne() {
                return Err(if gap == 1 {
                    adjacent_error(lift)
                } else {
                    equation_error(0, gap).at(lift.span)
                });
            }
        } else {
            curvature = product;
        }
    }
    Ok(curvature)
}

fn adjacent_error(lift: &Lift) -> Error {
    Error::new(
        lift.span,
        format!(
            "adjacent input on interval [{},{}] is not closed",
            lift.from, lift.to
        ),
    )
}

fn equation_error(from: usize, to: usize) -> Error {
    Error::plain(format!(
        "defining system equation fails on interval [{from},{to}]"
    ))
}

fn bar_product(model: &mut Model, prime: Prime, left: &Lift, right: &Lift) -> Result<Vec<Term>> {
    let terms = odd_cobar::concat_support(model, prime, &left.terms, &right.terms)?;
    signed(model, prime, terms, left.negative_bar)
}

fn signed(
    model: &mut Model,
    prime: Prime,
    mut terms: Vec<Term>,
    negative: bool,
) -> Result<Vec<Term>> {
    if negative {
        model.spend(terms.len())?;
        for (_, coefficient) in &mut terms {
            *coefficient = prime.negate(*coefficient);
        }
    }
    Ok(terms)
}

fn support_cells(terms: &[Term]) -> usize {
    terms.iter().map(odd_cobar::term_cells).sum()
}

#[derive(Default)]
struct Row {
    values: BTreeMap<usize, Vec<Term>>,
    cells: usize,
}

impl Row {
    fn add(&mut self, model: &mut Model, prime: Prime, to: usize, terms: Vec<Term>) -> Result<()> {
        if terms.is_empty() {
            return Ok(());
        }
        model.spend(1 + self.values.len().max(1).ilog2() as usize)?;
        match self.values.entry(to) {
            BTreeEntry::Vacant(entry) => {
                self.cells += 1 + support_cells(&terms);
                odd_cobar::check_storage(self.cells)?;
                entry.insert(terms);
            }
            BTreeEntry::Occupied(mut entry) => {
                let previous = entry.get();
                let sum = odd_cobar::add_support(model, prime, previous, &terms)?;
                self.cells -= support_cells(previous);
                self.cells += support_cells(&sum);
                odd_cobar::check_storage(self.cells)?;
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
