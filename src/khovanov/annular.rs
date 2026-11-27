use super::{
    choose, persist, Braid, Cell, Chain, ChainSum, Context, Grade, Resolution, Slice, State,
};
use crate::error::{Error, Result};
use crate::integer_smith as integral;
use crate::linear::Basis;
use crate::model::{Model, Value};
use crate::rv_lens;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

impl Resolution {
    pub(super) fn visit_labels(
        &self,
        model: &mut Model,
        xs: usize,
        slice: Option<Slice>,
        mut visit: impl FnMut(&mut Model, &[usize]) -> Result<()>,
    ) -> Result<()> {
        let circles = self.anchors.len();
        let essential = self.essential_count;
        let trivial = circles - essential;
        let Some(slice) = slice else {
            return choose(model, circles, xs, visit);
        };
        if matches!(slice, Slice::AtMost(k) if k >= essential as i64) {
            return choose(model, circles, xs, visit);
        }
        let mut first = xs.saturating_sub(trivial);
        let mut last = xs.min(essential);
        match slice {
            Slice::Exact(k) => {
                let difference = essential as i128 - i128::from(k);
                if difference < 0 || difference % 2 != 0 || difference / 2 > last as i128 {
                    return Ok(());
                }
                let count = (difference / 2) as usize;
                first = first.max(count);
                last = count;
            }
            Slice::AtMost(k) => {
                let needed = (essential as i128 - i128::from(k) + 1) / 2;
                if needed > last as i128 {
                    return Ok(());
                }
                first = first.max(needed as usize);
            }
        }
        if first > last {
            return Ok(());
        }
        model.spend(circles)?;
        let (a, b): (Vec<_>, Vec<_>) = (0..circles).partition(|&i| self.essential[i]);
        for count in first..=last {
            choose(model, essential, count, |model, left| {
                choose(model, trivial, xs - count, |model, right| {
                    model.spend(xs * (1 + xs.max(1).ilog2() as usize))?;
                    let mut labels: Vec<_> = left.iter().map(|&i| a[i]).collect();
                    labels.extend(right.iter().map(|&j| b[j]));
                    labels.sort_unstable();
                    visit(model, &labels)
                })
            })?;
        }
        Ok(())
    }
}

#[derive(Default)]
struct BinarySum {
    cells: usize,
    terms: BTreeSet<Cell>,
}

impl BinarySum {
    fn toggle(&mut self, model: &mut Model, cell: Cell) -> Result<()> {
        let width = cell.cells();
        model.spend(width * (1 + self.terms.len().max(1).ilog2() as usize))?;
        if self.terms.remove(&cell) {
            self.cells -= width;
        } else {
            self.cells += width;
            integral::storage(self.cells)?;
            self.terms.insert(cell);
        }
        Ok(())
    }

    fn finish(self, model: &mut Model) -> Result<Chain> {
        persist(
            model,
            self.terms.into_iter().map(|cell| (cell, 1)).collect(),
        )
    }
}

pub(crate) fn add_binary(model: &mut Model, a: &Chain, b: &Chain) -> Result<Chain> {
    let mut sum = BinarySum::default();
    for (cell, _) in a.iter().chain(b.iter()) {
        sum.toggle(model, cell.clone())?;
    }
    sum.finish(model)
}

pub(crate) fn pair(model: &mut Model, a: &Chain, b: &Chain) -> Result<bool> {
    let mut a = a.iter().peekable();
    let mut b = b.iter().peekable();
    let mut parity = false;
    while let (Some((x, _)), Some((y, _))) = (a.peek(), b.peek()) {
        model.spend(x.cells().min(y.cells()))?;
        let order = x.cmp(y);
        if !order.is_gt() {
            a.next();
        }
        if !order.is_lt() {
            b.next();
        }
        parity ^= order.is_eq();
    }
    Ok(parity)
}

pub(crate) struct Transverse {
    pub(crate) grade: Grade,
    pub(crate) strands: i64,
    pub(crate) boundary: Option<(i64, Chain)>,
    pub(crate) lower: i64,
    pub(crate) dual: Chain,
}

fn dual_from(
    model: &mut Model,
    basis: &Basis,
    target: &[Cell],
    chosen: usize,
    omitted: &[usize],
) -> Result<Chain> {
    let mut coefficients = vec![false; target.len()];
    let mut support = Vec::new();
    let mut width = 0;
    for (i, cell) in target.iter().enumerate() {
        model.spend(1 + omitted.len().max(1).ilog2() as usize)?;
        let mut value = i == chosen;
        if omitted.binary_search(&i).is_err() {
            if let Some(row) = basis.vector(i) {
                // Every row has largest index i. Annihilating that
                // row determines the value at i from smaller indices.
                model.spend(row.len())?;
                value = row[..row.len() - 1]
                    .iter()
                    .fold(false, |s, &j| s ^ coefficients[j]);
            }
        }
        coefficients[i] = value;
        if value {
            width += cell.cells();
            integral::storage(width)?;
            support.push((cell.clone(), 1));
        }
    }
    persist(model, support)
}

fn state_count(crossings: usize) -> Result<usize> {
    1usize
        .checked_shl(crossings as u32)
        .ok_or_else(|| Error::plain("annular Khovanov cube exceeds the work limit"))
}

fn binomials(model: &mut Model, n: usize) -> Result<Vec<i128>> {
    model.spend(n + 1)?;
    integral::storage(2 * (n + 1))?;
    let mut factors = Vec::with_capacity(n + 1);
    let mut value = 1i128;
    for i in 0..=n {
        factors.push(value);
        if i < n {
            let a = (n - i) as u128;
            let b = (i + 1) as u128;
            let common = rv_lens::gcd(a, b);
            value = (value / (b / common) as i128)
                .checked_mul((a / common) as i128)
                .ok_or_else(integral::arithmetic)?;
        }
    }
    Ok(factors)
}

impl Braid {
    pub(crate) fn part(&self, model: &mut Model, input: &Chain, k: i64) -> Result<Chain> {
        let mut context = Context::with_annular(self, true);
        let mut terms = Vec::new();
        for (cell, coefficient) in input.iter() {
            if context.resolution(model, cell.state)?.k(model, &cell.xs)? == k {
                terms.push((cell.clone(), *coefficient));
            }
        }
        persist(model, terms)
    }

    pub(crate) fn to_f2(&self, model: &mut Model, input: &Chain, cutoff: i64) -> Result<Chain> {
        let mut context = Context::with_annular(self, true);
        let mut terms = Vec::new();
        model.spend(input.len())?;
        for (cell, coefficient) in input.iter() {
            if coefficient.rem_euclid(2) == 1
                && context.resolution(model, cell.state)?.k(model, &cell.xs)? <= cutoff
            {
                terms.push((cell.clone(), 1));
            }
        }
        persist(model, terms)
    }

    pub(crate) fn binary_differential(
        &self,
        model: &mut Model,
        input: &Chain,
        cutoff: i64,
        reverse: bool,
    ) -> Result<Chain> {
        let mut context = Context::with_annular(self, true);
        let mut sum = BinarySum::default();
        for (cell, _) in input.iter() {
            context.edges(
                model,
                cell,
                reverse,
                Some(Slice::AtMost(cutoff)),
                |model, target, _| sum.toggle(model, target),
            )?;
        }
        sum.finish(model)
    }

    pub(crate) fn binary_elements(
        &self,
        model: &mut Model,
        grade: Grade,
        cutoff: i64,
        limit: usize,
    ) -> Result<Vec<Value>> {
        let basis =
            Context::with_annular(self, true).basis(model, grade, Some(Slice::AtMost(cutoff)))?;
        let size = 1usize
            .checked_shl(basis.len() as u32)
            .filter(|&size| size <= limit)
            .ok_or_else(|| {
                Error::plain("filtered Khovanov space exceeds the finite element limit")
            })?;
        let width: usize = basis.iter().map(Cell::cells).sum();
        model.allocate(size + (size / 2) * width)?;
        model.spend(size * (1 + basis.len()) * (1 + size.ilog2() as usize))?;
        let mut elements: Vec<Chain> = Vec::with_capacity(size);
        for state in 0..size {
            elements.push(
                basis
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| state & (1 << i) != 0)
                    .map(|(_, cell)| (cell.clone(), 1))
                    .collect(),
            );
        }
        elements.sort_unstable();
        Ok(elements.into_iter().map(Value::KhChain).collect())
    }

    fn psi_cell(&self, model: &mut Model) -> Result<(Grade, Cell)> {
        model.spend(self.crossings.len() + self.ends.len())?;
        let state = self
            .crossings
            .iter()
            .enumerate()
            .fold(0u128, |state, (i, edge)| {
                state | (u128::from(edge.negative) << i)
            });
        let strands = self.ends.len();
        let cell = Cell {
            state,
            xs: (0..strands).collect(),
        };
        Ok((self.grade(&cell, strands), cell))
    }

    pub(crate) fn psi(&self, model: &mut Model) -> Result<(Grade, i64, Chain)> {
        let (grade, cell) = self.psi_cell(model)?;
        Ok((
            grade,
            -(self.ends.len() as i64),
            persist(model, vec![(cell, 1)])?,
        ))
    }

    pub(crate) fn transverse(&self, model: &mut Model) -> Result<Rc<Transverse>> {
        if let Some(result) = self.transverse.borrow().as_ref() {
            return Ok(result.clone());
        }
        let result = Rc::new(self.transverse_complex(model)?);
        *self.transverse.borrow_mut() = Some(result.clone());
        Ok(result)
    }

    fn transverse_complex(&self, model: &mut Model) -> Result<Transverse> {
        let (grade, psi) = self.psi_cell(model)?;
        let strands = self.ends.len() as i64;
        let mut context = Context::with_annular(self, true);
        let target = context.basis(model, grade, None)?;
        let psi_index = target
            .binary_search(&psi)
            .map_err(|_| Error::plain("transverse cell is absent from the Khovanov basis"))?;
        context.edges(
            model,
            &psi,
            false,
            Some(Slice::AtMost(strands)),
            |_, _, _| Err(Error::plain("transverse element is not a Khovanov cycle")),
        )?;
        let source = context.basis(model, grade.shift(-1)?, None)?;
        let width = source
            .iter()
            .chain(target.iter())
            .map(Cell::cells)
            .sum::<usize>()
            + 3 * source.len()
            + target.len();
        integral::storage(width)?;
        let mut ordered = Vec::with_capacity(source.len());
        for (i, cell) in source.iter().enumerate() {
            let k = context.resolution(model, cell.state)?.k(model, &cell.xs)?;
            ordered.push((k, i));
        }
        model.spend(ordered.len() * (1 + ordered.len().max(1).ilog2() as usize))?;
        ordered.sort_unstable();

        let mut span = Basis::default();
        let mut last_remainder = vec![psi_index];
        for band in ordered.chunk_by(|a, b| a.0 == b.0) {
            let k = band[0].0;
            let mut new_pivots = Vec::new();
            for &(_, index) in band {
                let mut image = Vec::new();
                context.edges(
                    model,
                    &source[index],
                    false,
                    Some(Slice::AtMost(k)),
                    |model, cell, _| {
                        model.spend(cell.cells() * (1 + target.len().max(1).ilog2() as usize))?;
                        let coordinate = target.binary_search(&cell).map_err(|_| {
                            Error::plain("Khovanov differential has the wrong bidegree")
                        })?;
                        image.push(coordinate);
                        Ok(())
                    },
                )?;
                model.spend(image.len() * (1 + image.len().max(1).ilog2() as usize))?;
                image.sort_unstable();
                debug_assert!(image.windows(2).all(|a| a[0] != a[1]));
                span.insert_traced(model, image, vec![index], |_, reduced, _| {
                    if let Some(&pivot) = reduced.last() {
                        new_pivots.push(pivot);
                    }
                    Ok(())
                })?;
                integral::storage(width + span.cells())?;
            }
            let (remainder, witness) =
                span.reduce_with_witness(model, vec![psi_index], Vec::new())?;
            if remainder.is_empty() {
                model.spend(new_pivots.len() * (1 + new_pivots.len().max(1).ilog2() as usize))?;
                new_pivots.sort_unstable();
                // Existing pivots do not change on insertion. Omitting
                // this band's pivots recovers the preceding filtration.
                debug_assert!(k - 2 >= -strands);
                let dual = dual_from(model, &span, &target, last_remainder[0], &new_pivots)?;
                let primitive = persist(
                    model,
                    witness
                        .into_iter()
                        .map(|i| (source[i].clone(), 1))
                        .collect(),
                )?;
                return Ok(Transverse {
                    grade,
                    strands,
                    boundary: Some((k, primitive)),
                    lower: k - 2,
                    dual,
                });
            }
            last_remainder = remainder;
        }
        let dual = dual_from(model, &span, &target, last_remainder[0], &[])?;
        Ok(Transverse {
            grade,
            strands,
            boundary: None,
            lower: strands,
            dual,
        })
    }

    pub(crate) fn annular_euler(&self, model: &mut Model) -> Result<Vec<((i64, i64), i128)>> {
        let count = state_count(self.crossings.len())?;
        model.spend(count)?;
        let mut euler = BTreeMap::<(i64, i64), i128>::new();
        for state in 0..count {
            let resolution = self.resolve(model, state as State, true)?;
            let circles = resolution.anchors.len();
            let essential = resolution.essential_count;
            let a = binomials(model, essential)?;
            let b = binomials(model, circles - essential)?;
            let ones = state.count_ones() as usize;
            let j = self.shift(ones) + circles as i64;
            let sign = if ones.abs_diff(self.negative).is_multiple_of(2) {
                1
            } else {
                -1
            };
            for (i, &first) in a.iter().enumerate() {
                let k = essential as i64 - 2 * i as i64;
                for (l, &second) in b.iter().enumerate() {
                    model.spend(1 + euler.len().max(1).ilog2() as usize)?;
                    let degree = (j - 2 * (i + l) as i64, k);
                    let coefficient = (sign * first)
                        .checked_mul(second)
                        .ok_or_else(integral::arithmetic)?;
                    let entry = euler.entry(degree).or_default();
                    *entry = entry
                        .checked_add(coefficient)
                        .ok_or_else(integral::arithmetic)?;
                    if *entry == 0 {
                        euler.remove(&degree);
                    }
                    integral::storage(7 * euler.len())?;
                }
            }
        }
        Ok(euler.into_iter().collect())
    }

    pub(crate) fn annular_check(&self, model: &mut Model) -> Result<usize> {
        let states = state_count(self.crossings.len())?;
        model.spend(states)?;
        let mut context = Context::with_annular(self, true);
        let slice = Some(Slice::AtMost(self.ends.len() as i64));
        let mut tested = 0;
        for state in 0..states {
            let resolution = context.resolution(model, state as State)?;
            let circles = resolution.anchors.len();
            let labels = 1usize.checked_shl(circles as u32).ok_or_else(|| {
                Error::plain("annular Khovanov enhanced cube exceeds the work limit")
            })?;
            model.spend(labels)?;
            tested += labels;
            for bits in 0..labels {
                model.spend(circles)?;
                let xs = (0..circles).filter(|&i| (bits & (1 << i)) != 0).collect();
                let cell = Cell {
                    state: state as State,
                    xs,
                };
                let mut first = Vec::new();
                context.edges(model, &cell, false, slice, |_, term, sign| {
                    first.push((term, sign));
                    Ok(())
                })?;
                let mut output = ChainSum::default();
                for (term, sign) in first {
                    context.edges(model, &term, false, slice, |model, target, second| {
                        output.add(model, target, sign * second)
                    })?;
                }
                // The two verified degrees of d make the components
                // of d² precisely d0², d0*d-2+d-2*d0, and d-2².
                if !output.terms.is_empty() {
                    return Err(Error::plain(
                        "annular Khovanov differential does not square to zero",
                    ));
                }
            }
        }
        Ok(tested)
    }
}
