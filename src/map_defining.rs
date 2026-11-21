use std::rc::Rc;

use crate::error::{Error, Result};
use crate::map_ext::{self, MapAdams, Term};
use crate::model::{Kind, Model, TypeId, Value};

type Input = (TypeId, Rc<[usize]>);

pub(crate) fn predicate(model: &mut Model, target: TypeId, inputs: &[Input]) -> Result<Value> {
    let size = model.elements(target)?.len();
    model.allocate(1 + size)?;
    let mut rows = vec![Value::Point(0); size];
    let mut remaining = size;
    search(model, target, inputs, |model, coordinates| {
        let index = model.index(target, &Value::Class(coordinates.into()))?;
        if matches!(rows[index], Value::Point(0)) {
            rows[index] = Value::Point(1);
            remaining -= 1;
        }
        Ok(remaining == 0)
    })?;
    Ok(Value::Map(rows.into()))
}

pub(crate) fn member(
    model: &mut Model,
    target: TypeId,
    inputs: &[Input],
    query: &[usize],
) -> Result<Value> {
    let found = search(model, target, inputs, |model, coordinates| {
        model.spend(query.len() + coordinates.len())?;
        Ok(query == coordinates)
    })?;
    Ok(Value::Point(usize::from(found)))
}

fn search(
    model: &mut Model,
    target: TypeId,
    inputs: &[Input],
    mut visit: impl FnMut(&mut Model, Vec<usize>) -> Result<bool>,
) -> Result<bool> {
    let first = inputs.len();
    let mut search = Search::new(model, inputs)?;
    let top = search.intervals.len() - 1;
    let mut next = first;
    loop {
        model.spend(1)?;
        if next == top {
            let curvature = search.obstruction(model, top)?;
            let coordinates = if curvature.is_empty() {
                Vec::new()
            } else {
                model.map_ext(target)?.project(model, &curvature)?
            };
            if visit(model, coordinates)? {
                return Ok(true);
            }
        } else if search.choose_primitive(model, next)? {
            next += 1;
            continue;
        } else if next < search.offsets[2] {
            // Pair obstructions depend only on the fixed adjacent entries.
            return Ok(false);
        }

        loop {
            if next == first {
                return Ok(false);
            }
            next -= 1;
            if search.advance_choice(model, next)? {
                next += 1;
                break;
            }
        }
    }
}

#[derive(Clone, Copy)]
struct Interval {
    from: usize,
    to: usize,
    s: i128,
    t: i128,
}

#[derive(Default)]
struct Witness {
    terms: Vec<Term>,
    cells: usize,
    choices: Option<Choices>,
}

impl Witness {
    fn replace(&mut self, model: &mut Model, live: &mut usize, terms: Vec<Term>) -> Result<()> {
        let cells = terms.iter().map(Term::cells).sum();
        let next = *live - self.cells + cells;
        map_ext::check_storage(next)?;
        model.spend(cells)?;
        self.cells = cells;
        self.terms = terms;
        *live = next;
        Ok(())
    }
}

struct Choices {
    group: Rc<MapAdams>,
    bits: Vec<bool>,
}

struct Search {
    objects: Vec<TypeId>,
    offsets: Vec<usize>,
    intervals: Vec<Interval>,
    witnesses: Vec<Witness>,
    cells: usize,
}

impl Search {
    fn new(model: &mut Model, inputs: &[Input]) -> Result<Self> {
        let arity = inputs.len();
        let size = arity
            .checked_add(1)
            .and_then(|endpoints| arity.checked_mul(endpoints))
            .ok_or_else(|| Error::plain("Massey interval count exceeds the storage limit"))?
            / 2;
        let interval_cells = (std::mem::size_of::<Interval>() + std::mem::size_of::<Witness>())
            .div_ceil(std::mem::size_of::<usize>());
        let index_cells = arity.saturating_add(1).saturating_mul(5);
        model.allocate(
            size.saturating_mul(interval_cells)
                .saturating_add(index_cells),
        )?;
        model.spend(size)?;

        let mut objects = Vec::with_capacity(arity + 1);
        let mut prefixes = Vec::with_capacity(arity + 1);
        let mut witnesses = Vec::with_capacity(size - 1);
        let mut cells = 0;
        let (mut total_s, mut total_t) = (0i128, 0i128);
        prefixes.push((0, 0));
        // Calls list the outermost map first. Intervals follow source order.
        for (ty, class) in inputs.iter().rev() {
            let Kind::MapExt(source, target, s, t) = *model.kind(*ty) else {
                unreachable!();
            };
            if objects.is_empty() {
                objects.push(source);
            }
            objects.push(target);
            total_s += s as i128;
            total_t += i128::from(t);
            prefixes.push((total_s, total_t));

            let mut witness = Witness::default();
            if !class.is_empty() {
                let terms = model.map_ext(*ty)?.representative(model, class)?;
                witness.replace(model, &mut cells, terms)?;
            }
            witnesses.push(witness);
        }
        witnesses.resize_with(size - 1, Witness::default);

        let mut offsets = Vec::with_capacity(arity);
        let mut intervals = Vec::with_capacity(size);
        for span in 1..=arity {
            offsets.push(intervals.len());
            for from in 0..=arity - span {
                let to = from + span;
                intervals.push(Interval {
                    from,
                    to,
                    s: prefixes[to].0 - prefixes[from].0 + 1 - span as i128,
                    t: prefixes[to].1 - prefixes[from].1,
                });
            }
        }
        Ok(Self {
            objects,
            offsets,
            intervals,
            witnesses,
            cells,
        })
    }

    fn ty(&self, model: &mut Model, interval: Interval, s: i128) -> Result<TypeId> {
        let s = usize::try_from(s).map_err(|_| intermediate_overflow())?;
        let t = i64::try_from(interval.t).map_err(|_| intermediate_overflow())?;
        model.map_ext_type(self.objects[interval.from], self.objects[interval.to], s, t)
    }

    fn choose_primitive(&mut self, model: &mut Model, index: usize) -> Result<bool> {
        let interval = self.intervals[index];
        let obstruction = self.obstruction(model, index)?;
        let primitive = if obstruction.is_empty() {
            Vec::new()
        } else {
            if interval.s < 0 {
                return Ok(false);
            }
            let ty = self.ty(model, interval, interval.s + 1)?;
            let Some(primitive) = model.map_ext(ty)?.primitive(model, &obstruction)? else {
                return Ok(false);
            };
            primitive
        };
        self.witnesses[index].replace(model, &mut self.cells, primitive)?;
        Ok(true)
    }

    fn advance_choice(&mut self, model: &mut Model, index: usize) -> Result<bool> {
        let interval = self.intervals[index];
        if interval.s < 0 {
            return Ok(false);
        }
        if self.witnesses[index].choices.is_none() {
            let ty = self.ty(model, interval, interval.s)?;
            let group = model.map_ext(ty)?;
            model.allocate(group.rank())?;
            let bits = vec![false; group.rank()];
            self.witnesses[index].choices = Some(Choices { group, bits });
        }

        let witness = &mut self.witnesses[index];
        let Choices { group, bits } = witness.choices.as_mut().unwrap();
        let carry = bits.iter().position(|bit| !bit).unwrap_or(bits.len());
        model.spend(1 + carry)?;
        bits[..carry].fill(false);
        if carry == bits.len() {
            witness.replace(model, &mut self.cells, Vec::new())?;
            return Ok(false);
        }
        bits[carry] = true;

        // A carry toggles exactly these basis classes. Exact changes are
        // absorbed by a unitriangular gauge on longer intervals.
        let mut terms = std::mem::take(&mut witness.terms);
        for coordinate in 0..=carry {
            let generator = group.representative(model, &[coordinate])?;
            terms = map_ext::add_terms(model, &terms, &generator)?;
        }
        witness.replace(model, &mut self.cells, terms)?;
        Ok(true)
    }

    fn obstruction(&self, model: &mut Model, index: usize) -> Result<Vec<Term>> {
        let Interval { from, to, .. } = self.intervals[index];
        let mut sum = Vec::new();
        for split in from + 1..to {
            model.spend(1)?;
            let inner = &self.witnesses[self.offsets[split - from - 1] + from].terms;
            let outer = &self.witnesses[self.offsets[to - split - 1] + split].terms;
            if inner.is_empty() || outer.is_empty() {
                continue;
            }
            let product = map_ext::compose_terms(model, outer, inner)?;
            sum = if sum.is_empty() {
                product
            } else {
                map_ext::add_terms(model, &sum, &product)?
            };
        }
        Ok(sum)
    }
}

fn intermediate_overflow() -> Error {
    Error::plain("Massey defining system intermediate bidegree is out of range")
}
