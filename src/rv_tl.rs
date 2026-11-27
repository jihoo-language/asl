use std::cmp::Ordering;
use std::collections::{btree_map::Entry, BTreeMap};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_cyclotomic;
use crate::rv_lens::{self, Rational};

const MAX_CELLS: usize = 300_000;
const MAX_WORD: usize = 4096;
type Path = Rc<[usize]>;
pub(crate) type Support = Rc<[(Path, Rational)]>;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) level: usize,
    pub(crate) strands: usize,
    pub(crate) height: usize,
    pub(crate) order: usize,
}

impl Grade {
    pub(crate) fn checked(level: usize, strands: usize, height: usize) -> Result<Self> {
        if level == 0 {
            return Err(Error::plain("Temperley–Lieb level must be at least 1"));
        }
        if height > level {
            return Err(Error::plain(
                "Temperley–Lieb endpoint height exceeds the level",
            ));
        }
        let order = level
            .checked_add(2)
            .and_then(|n| n.checked_mul(4))
            .ok_or_else(|| Error::plain("Temperley–Lieb root order is out of range"))?;
        rv_cyclotomic::check_order(order)?;
        Ok(Self {
            level,
            strands,
            height,
            order,
        })
    }

    pub(crate) fn possible(self) -> bool {
        self.height <= self.strands && (self.strands - self.height).is_multiple_of(2)
    }

    pub(crate) fn validate_path(self, model: &mut Model, path: &[usize]) -> Result<()> {
        storage(4 + path.len())?;
        model.spend(path.len())?;
        if path.len().checked_sub(1) != Some(self.strands)
            || path.first() != Some(&0)
            || path.last() != Some(&self.height)
        {
            return Err(Error::plain("Temperley–Lieb path must have the indicated length and endpoints, starting at zero"));
        }
        if path.iter().any(|&height| height > self.level)
            || path.windows(2).any(|a| a[0].abs_diff(a[1]) != 1)
        {
            return Err(Error::plain(
                "Temperley–Lieb path requires unit steps within the level bounds",
            ));
        }
        Ok(())
    }
}

fn arithmetic() -> Error {
    Error::plain("Temperley–Lieb integer arithmetic is out of range")
}

fn storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain(
            "Temperley–Lieb calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

fn cells(path: &[usize], coefficient: &Rational) -> usize {
    5 + path.len() + 2 * coefficient.numerator.len()
}

fn support_cells(terms: &Support) -> usize {
    terms
        .iter()
        .map(|(path, coefficient)| cells(path, coefficient))
        .sum()
}

fn persist(model: &mut Model, terms: Support) -> Result<Value> {
    let total = support_cells(&terms);
    storage(total)?;
    model.allocate(total)?;
    Ok(Value::RvTl(terms))
}

pub(crate) fn zero() -> Value {
    Value::RvTl([].into())
}

pub(crate) fn basis(model: &mut Model, grade: Grade, path: Vec<usize>) -> Result<Value> {
    grade.validate_path(model, &path)?;
    persist(model, [(path.into(), Rational::integer(1))].into())
}

pub(crate) fn coordinate(
    model: &mut Model,
    grade: Grade,
    input: &Support,
    path: &[usize],
) -> Result<Value> {
    grade.validate_path(model, path)?;
    model.spend(path.len() * (1 + input.len().max(1).ilog2() as usize))?;
    Ok(Value::RVCycloQ(
        input
            .binary_search_by(|(p, _)| p.as_ref().cmp(path))
            .map(|i| input[i].1.clone())
            .unwrap_or_else(|_| Rational::zero()),
    ))
}

pub(crate) fn compare(model: &mut Model, a: &Support, b: &Support) -> Result<Ordering> {
    if Rc::ptr_eq(a, b) {
        return Ok(Ordering::Equal);
    }
    for (x, y) in a.iter().zip(b.iter()) {
        model.spend(1 + x.0.len() + x.1.numerator.len().min(y.1.numerator.len()))?;
        let order = x.cmp(y);
        if !order.is_eq() {
            return Ok(order);
        }
    }
    Ok(a.len().cmp(&b.len()))
}

#[derive(Default)]
struct Sum {
    terms: BTreeMap<Path, Rational>,
    cells: usize,
}

impl Sum {
    fn insert(&mut self, model: &mut Model, path: Path, value: Rational) -> Result<()> {
        model.spend((path.len() + 1) * (1 + self.terms.len().max(1).ilog2() as usize))?;
        if value.is_zero() {
            return Ok(());
        }
        let width = cells(&path, &value);
        match self.terms.entry(path) {
            Entry::Vacant(slot) => {
                self.cells += width;
                storage(self.cells)?;
                slot.insert(value);
            }
            Entry::Occupied(mut slot) => {
                let final_value = rv_lens::sum(model, slot.get(), &value)?;
                self.cells -= cells(slot.key(), slot.get());
                if final_value.is_zero() {
                    slot.remove();
                } else {
                    self.cells += cells(slot.key(), &final_value);
                    storage(self.cells)?;
                    *slot.get_mut() = final_value;
                }
            }
        }
        Ok(())
    }
}

pub(crate) fn add(model: &mut Model, a: &Support, b: &Support) -> Result<Value> {
    let mut sum = Sum::default();
    for (path, coefficient) in a.iter().chain(b.iter()) {
        sum.insert(model, path.clone(), coefficient.clone())?;
    }
    persist(model, sum.terms.into_iter().collect())
}

pub(crate) fn negate(model: &mut Model, input: &Support) -> Result<Value> {
    let mut terms = Vec::new();
    for (path, c) in input.iter() {
        terms.push((path.clone(), rv_lens::negative(model, c)?));
    }
    persist(model, terms.into())
}

pub(crate) fn scale(
    model: &mut Model,
    grade: Grade,
    input: &Support,
    c: &Rational,
) -> Result<Value> {
    if c.is_zero() {
        return Ok(zero());
    }
    let mut output = Vec::new();
    let mut kept = 0;
    for (path, coefficient) in input.iter() {
        let scaled = rv_lens::multiply(model, grade.order, coefficient, c)?;
        kept += cells(path, &scaled);
        storage(kept)?;
        output.push((path.clone(), scaled));
    }
    persist(model, output.into())
}

#[derive(Clone, Copy)]
enum Generator {
    Tl(usize),
    Braid(usize, bool),
}

fn check_index(strands: usize, index: usize) -> Result<()> {
    if index == 0 || index >= strands {
        return Err(Error::plain(
            "Temperley–Lieb/Artin index must satisfy 1 <= i < strands",
        ));
    }
    Ok(())
}

fn signed_word(model: &mut Model, strands: usize, word: &[i128]) -> Result<Vec<Generator>> {
    if word.len() > MAX_WORD {
        return Err(Error::plain(
            "Temperley–Lieb braid word exceeds the limit of 4096 letters",
        ));
    }
    model.spend(word.len())?;
    word.iter()
        .map(|&letter| {
            let index = usize::try_from(letter.unsigned_abs())
                .map_err(|_| Error::plain("Artin generator index is out of range"))?;
            check_index(strands, index)?;
            Ok(Generator::Braid(index, letter > 0))
        })
        .collect()
}

#[derive(Clone, Copy)]
struct Block {
    start: usize,
    length: usize,
}

impl Block {
    fn checked(grade: Grade, first: usize, length: usize) -> Result<Self> {
        if length > grade.level + 1 {
            return Err(Error::plain(
                "Jones–Wenzl projector length exceeds the level plus one",
            ));
        }
        let Some(start) = first.checked_sub(1) else {
            return Err(Error::plain(
                "Jones–Wenzl block starts at a one-based strand index",
            ));
        };
        if length > grade.strands || start > grade.strands - length {
            return Err(Error::plain(
                "Jones–Wenzl block lies outside the path strands",
            ));
        }
        Ok(Self { start, length })
    }
}

#[derive(Default)]
struct Context {
    quantum: BTreeMap<(usize, bool), Rational>,
    ratios: BTreeMap<(usize, usize), Rational>,
    projectors: BTreeMap<Path, Support>,
    cells: usize,
}

impl Context {
    fn dimension(
        &mut self,
        model: &mut Model,
        grade: Grade,
        height: usize,
        inverse: bool,
    ) -> Result<Rational> {
        if let Some(existing) = self.quantum.get(&(height, inverse)) {
            return Ok(existing.clone());
        }
        let order = grade.order;
        let n = grade.level + 2;
        let m = height + 1;
        let sign = if height.is_multiple_of(2) { 1 } else { -1 };
        let mut polynomial = vec![0; order];
        let denominator = if inverse {
            // With u=q^(2m) of order l, (1-u) sum_j(l-1-j)u^j=l.
            let l = n / rv_lens::gcd(n as u128, m as u128) as usize;
            model.spend(l)?;
            for j in 0..l - 1 {
                let power = (2 * (m - 1) + 4 * m * j) % order;
                let coefficient = sign * (l - 1 - j) as i128;
                polynomial[power] += coefficient;
                polynomial[(power + 4) % order] -= coefficient;
            }
            l as u128
        } else {
            model.spend(m)?;
            for j in 0..m {
                let power = (2 * (m as i128 - 1 - 2 * j as i128)).rem_euclid(order as i128);
                polynomial[power as usize] += sign;
            }
            1
        };
        let numerator = rv_cyclotomic::reduced(model, order, polynomial)?;
        let value = rv_lens::fraction(model, &numerator, denominator)?;
        self.cells += 5 + 2 * value.numerator.len();
        storage(self.cells)?;
        self.quantum.insert((height, inverse), value.clone());
        Ok(value)
    }

    fn ratio(
        &mut self,
        model: &mut Model,
        grade: Grade,
        outer: usize,
        middle: usize,
    ) -> Result<Rational> {
        if let Some(existing) = self.ratios.get(&(outer, middle)) {
            return Ok(existing.clone());
        }
        let a = self.dimension(model, grade, middle, false)?;
        let b = self.dimension(model, grade, outer, true)?;
        let value = rv_lens::multiply(model, grade.order, &a, &b)?;
        self.cells += 5 + 2 * value.numerator.len();
        storage(self.cells)?;
        self.ratios.insert((outer, middle), value.clone());
        Ok(value)
    }

    fn step(
        &mut self,
        model: &mut Model,
        grade: Grade,
        input: &Support,
        generator: Generator,
    ) -> Result<Support> {
        let (index, twist) = match generator {
            Generator::Tl(index) => (index, None),
            Generator::Braid(index, true) => (index, Some(1)),
            Generator::Braid(index, false) => (index, Some(-1)),
        };
        let mut sum = Sum::default();
        for (path, coefficient) in input.iter() {
            if let Some(twist) = twist {
                let diagonal = rv_lens::shifted(model, grade.order, coefficient, twist)?;
                sum.insert(model, path.clone(), diagonal)?;
            }
            let outer = path[index - 1];
            if outer == path[index + 1] {
                let lower = outer.checked_sub(1);
                let upper = (outer < grade.level).then_some(outer + 1);
                for middle in [lower, upper].into_iter().flatten() {
                    let ratio = self.ratio(model, grade, outer, middle)?;
                    let mut value = rv_lens::multiply(model, grade.order, &ratio, coefficient)?;
                    if let Some(twist) = twist {
                        value = rv_lens::shifted(model, grade.order, &value, -twist)?;
                    }
                    let mut updated = path.to_vec();
                    updated[index] = middle;
                    sum.insert(model, updated.into(), value)?;
                }
            }
            storage(sum.cells + self.cells)?;
        }
        Ok(sum.terms.into_iter().collect())
    }

    fn project_column(
        &mut self,
        model: &mut Model,
        grade: Grade,
        segment: Path,
        reserved: usize,
    ) -> Result<Support> {
        model.spend(1 + segment.len() * (1 + self.projectors.len().max(1).ilog2() as usize))?;
        if let Some(known) = self.projectors.get(&segment) {
            return Ok(known.clone());
        }
        let length = segment.len() - 1;
        let unit: Support = [(segment.clone(), Rational::integer(1))].into();
        if segment[0].abs_diff(segment[length]) == length {
            // Every elementary e_i vanishes on a monotone segment.
            return Ok(unit);
        }

        // f_(r+1) v = f_r v - (d_(r-1)/d_r) f_r e_r f_r v.
        let r = length - 1;
        let prefix = Block {
            start: 0,
            length: r,
        };
        let left = self.project(model, grade, &unit, prefix, reserved)?;
        let result = if left.is_empty() {
            left
        } else {
            let kept = reserved + support_cells(&left);
            storage(self.cells + kept)?;
            let bent = self.step(model, grade, &left, Generator::Tl(r))?;
            let right = self.project(model, grade, &bent, prefix, kept)?;
            let ratio = self.ratio(model, grade, r, r - 1)?;
            let factor = rv_lens::negative(model, &ratio)?;
            let kept = kept + support_cells(&right);
            storage(self.cells + kept)?;
            let mut sum = Sum::default();
            for (path, coefficient) in left.iter() {
                sum.insert(model, path.clone(), coefficient.clone())?;
                storage(self.cells + sum.cells + kept)?;
            }
            for (path, coefficient) in right.iter() {
                let value = rv_lens::multiply(model, grade.order, coefficient, &factor)?;
                sum.insert(model, path.clone(), value)?;
                storage(self.cells + sum.cells + kept)?;
            }
            sum.terms.into_iter().collect()
        };
        self.cells += 1 + segment.len() + support_cells(&result);
        storage(self.cells + reserved)?;
        self.projectors.insert(segment, result.clone());
        Ok(result)
    }

    fn project(
        &mut self,
        model: &mut Model,
        grade: Grade,
        input: &Support,
        block: Block,
        reserved: usize,
    ) -> Result<Support> {
        model.spend(1)?;
        if block.length == grade.level + 1 {
            // At level k the simple object labeled k+1 is negligible.
            // The bounded path quotient kills its Jones–Wenzl idempotent.
            return Ok([].into());
        }
        if block.length <= 1 || input.is_empty() {
            return Ok(input.clone());
        }
        let reserved = reserved + support_cells(input);
        storage(self.cells + reserved)?;
        let interval = block.start..=block.start + block.length;
        let mut sum = Sum::default();
        for (path, coefficient) in input.iter() {
            let local = self.project_column(
                model,
                grade,
                path[interval.clone()].into(),
                reserved + sum.cells,
            )?;
            let kept = reserved + support_cells(&local);
            for (segment, factor) in local.iter() {
                let value = rv_lens::multiply(model, grade.order, coefficient, factor)?;
                let mut updated = path.to_vec();
                updated[interval.clone()].copy_from_slice(segment);
                sum.insert(model, updated.into(), value)?;
                storage(self.cells + sum.cells + kept)?;
            }
        }
        Ok(sum.terms.into_iter().collect())
    }

    fn act(
        &mut self,
        model: &mut Model,
        grade: Grade,
        mut input: Support,
        word: &[Generator],
    ) -> Result<Support> {
        for &step in word.iter().rev() {
            if input.is_empty() {
                break;
            }
            input = self.step(model, grade, &input, step)?;
        }
        Ok(input)
    }
}

pub(crate) fn tl(model: &mut Model, grade: Grade, input: &Support, index: usize) -> Result<Value> {
    check_index(grade.strands, index)?;
    let result = Context::default().step(model, grade, input, Generator::Tl(index))?;
    persist(model, result)
}

pub(crate) fn project(
    model: &mut Model,
    grade: Grade,
    input: &Support,
    first: usize,
    length: usize,
) -> Result<Value> {
    let block = Block::checked(grade, first, length)?;
    let result = Context::default().project(model, grade, input, block, 0)?;
    persist(model, result)
}

pub(crate) fn braid(
    model: &mut Model,
    grade: Grade,
    input: &Support,
    word: &[i128],
) -> Result<Value> {
    let word = signed_word(model, grade.strands, word)?;
    let result = Context::default().act(model, grade, input.clone(), &word)?;
    persist(model, result)
}

fn visit_paths(
    model: &mut Model,
    grade: Grade,
    mut visit: impl FnMut(&mut Model, &[usize]) -> Result<()>,
) -> Result<()> {
    if !grade.possible() {
        return Ok(());
    }
    let length = grade.strands.checked_add(1).ok_or_else(arithmetic)?;
    storage(length.saturating_mul(2))?;
    let mut path = vec![0usize];
    let mut choices = vec![0u8];
    while let Some(choice) = choices.last_mut() {
        model.spend(1)?;
        if path.len() == length {
            visit(model, &path)?;
            choices.pop();
            path.pop();
        } else if *choice < 2 {
            let height = *path.last().unwrap();
            let next = if *choice == 0 {
                height.checked_sub(1)
            } else {
                (height < grade.level).then_some(height + 1)
            };
            *choice += 1;
            if let Some(next) = next {
                let remaining = length - path.len() - 1;
                // Parity is forced by unit steps and the chosen endpoint.
                if next.abs_diff(grade.height) <= remaining {
                    path.push(next);
                    choices.push(0);
                }
            }
        } else {
            choices.pop();
            path.pop();
        }
    }
    Ok(())
}

fn path_count(model: &mut Model, grade: Grade) -> Result<i128> {
    if !grade.possible() {
        return Ok(0);
    }
    let width = grade.level.min(grade.strands) + 1;
    model.spend(grade.strands.checked_mul(width).ok_or_else(arithmetic)?)?;
    let mut counts = vec![0; width];
    counts[0] = 1i128;
    for _ in 0..grade.strands {
        let mut next = vec![0; width];
        for (height, value) in next.iter_mut().enumerate() {
            let lower = height.checked_sub(1).map_or(0, |h| counts[h]);
            let upper = counts.get(height + 1).copied().unwrap_or(0);
            *value = lower.checked_add(upper).ok_or_else(arithmetic)?;
        }
        counts = next;
    }
    Ok(counts[grade.height])
}

fn trace_raw(
    model: &mut Model,
    grade: Grade,
    word: &[Generator],
    context: &mut Context,
) -> Result<Rational> {
    if word.is_empty() {
        return Ok(Rational::integer(path_count(model, grade)?));
    }
    let mut trace = Rational::zero();
    visit_paths(model, grade, |model, path| {
        model.spend(path.len())?;
        let initial: Support = [(path.into(), Rational::integer(1))].into();
        let output = context.act(model, grade, initial, word)?;
        if let Ok(index) = output.binary_search_by(|(p, _)| p.as_ref().cmp(path)) {
            trace = rv_lens::sum(model, &trace, &output[index].1)?;
        }
        Ok(())
    })?;
    Ok(trace)
}

pub(crate) fn trace(model: &mut Model, grade: Grade, word: &[i128]) -> Result<Value> {
    let word = signed_word(model, grade.strands, word)?;
    let value = trace_raw(model, grade, &word, &mut Context::default())?;
    rv_lens::rational_value(model, value)
}

fn power(
    model: &mut Model,
    order: usize,
    mut value: Rational,
    mut exponent: usize,
) -> Result<Rational> {
    let mut result = Rational::integer(1);
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = rv_lens::multiply(model, order, &result, &value)?;
        }
        exponent >>= 1;
        if exponent > 0 {
            value = rv_lens::multiply(model, order, &value, &value)?;
        }
    }
    Ok(result)
}

pub(crate) fn jones(model: &mut Model, grade: Grade, word: &[i128]) -> Result<Value> {
    let mut context = Context::default();
    let generators = signed_word(model, grade.strands, word)?;
    if word.is_empty() {
        let delta = context.dimension(model, grade, 1, grade.strands == 0)?;
        let value = if grade.strands == 0 {
            delta
        } else {
            power(model, grade.order, delta, grade.strands - 1)?
        };
        return rv_lens::rational_value(model, value);
    }
    let mut sum = Rational::zero();
    for height in (grade.strands % 2..=grade.strands.min(grade.level)).step_by(2) {
        let sector = Grade { height, ..grade };
        let trace = trace_raw(model, sector, &generators, &mut context)?;
        if !trace.is_zero() {
            let dimension = context.dimension(model, grade, height, false)?;
            let contribution = rv_lens::multiply(model, grade.order, &trace, &dimension)?;
            sum = rv_lens::sum(model, &sum, &contribution)?;
        }
    }
    let inverse = context.dimension(model, grade, 1, true)?;
    sum = rv_lens::multiply(model, grade.order, &sum, &inverse)?;
    let writhe: i128 = word.iter().map(|&letter| letter.signum()).sum();
    sum = rv_lens::shifted(model, grade.order, &sum, -3 * writhe)?;
    if writhe.rem_euclid(2) != 0 {
        sum = rv_lens::negative(model, &sum)?;
    }
    rv_lens::rational_value(model, sum)
}
fn cable(model: &mut Model, color: usize, braid: &[Generator]) -> Result<Vec<Generator>> {
    let crossings = color
        .checked_mul(color)
        .and_then(|block| block.checked_mul(braid.len()))
        .ok_or_else(arithmetic)?;
    if crossings > MAX_WORD {
        return Err(Error::plain(
            "cabled braid exceeds the limit of 4096 elementary crossings",
        ));
    }
    model.spend(crossings)?;
    let mut result = Vec::with_capacity(crossings);
    for &letter in braid {
        let Generator::Braid(index, positive) = letter else {
            unreachable!()
        };
        let start = result.len();
        let base = (index - 1) * color;
        for row in 0..color {
            for shift in (row + 1..=row + color).rev() {
                result.push(Generator::Braid(base + shift, positive));
            }
        }
        if !positive {
            result[start..].reverse();
        }
    }
    Ok(result)
}

pub(crate) fn colored_jones(
    model: &mut Model,
    grade: Grade,
    color: usize,
    word: &[i128],
) -> Result<Value> {
    if color > grade.level {
        return Err(Error::plain("Jones color must be at most the level"));
    }
    if color == 1 {
        return jones(model, grade, word);
    }
    let braid = signed_word(model, grade.strands, word)?;
    if color == 0 {
        return rv_lens::rational_value(model, Rational::integer(1));
    }
    let mut context = Context::default();
    let inverse = context.dimension(model, grade, color, true)?;
    if braid.is_empty() {
        // The identity on n copies of the color has quantum trace d_color^n.
        let value = if grade.strands == 0 {
            inverse
        } else {
            let dimension = context.dimension(model, grade, color, false)?;
            power(model, grade.order, dimension, grade.strands - 1)?
        };
        return rv_lens::rational_value(model, value);
    }
    let strands = grade.strands.checked_mul(color).ok_or_else(arithmetic)?;
    let cabled = cable(model, color, &braid)?;
    let mut sum = Rational::zero();
    for height in (strands % 2..=strands.min(grade.level)).step_by(2) {
        let sector = Grade {
            strands,
            height,
            ..grade
        };
        let mut trace = Rational::zero();
        visit_paths(model, sector, |model, path| {
            model.spend(path.len())?;
            let mut state: Support = [(path.into(), Rational::integer(1))].into();
            for strand in 0..grade.strands {
                let block = Block {
                    start: strand * color,
                    length: color,
                };
                state = context.project(model, sector, &state, block, 0)?;
                if state.is_empty() {
                    return Ok(());
                }
            }
            let image = context.act(model, sector, state, &cabled)?;
            if let Ok(index) = image.binary_search_by(|(p, _)| p.as_ref().cmp(path)) {
                trace = rv_lens::sum(model, &trace, &image[index].1)?;
            }
            Ok(())
        })?;
        if !trace.is_zero() {
            let dimension = context.dimension(model, grade, height, false)?;
            let weighted = rv_lens::multiply(model, grade.order, &trace, &dimension)?;
            sum = rv_lens::sum(model, &sum, &weighted)?;
        }
    }
    sum = rv_lens::multiply(model, grade.order, &sum, &inverse)?;
    let writhe: i128 = word.iter().map(|&letter| letter.signum()).sum();
    let twist = color as i128 * writhe;
    sum = rv_lens::shifted(model, grade.order, &sum, -twist * (color as i128 + 2))?;
    if twist.rem_euclid(2) != 0 {
        sum = rv_lens::negative(model, &sum)?;
    }
    rv_lens::rational_value(model, sum)
}
