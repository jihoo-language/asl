use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use crate::complex::{xor_sorted, MAX_PIVOT_CELLS};
use crate::error::{Error, Result};
use crate::linear::Basis;
use crate::model::{Kind, Model, TypeId, Value};
use crate::rv_ideal::{self, Ideal};
use crate::rv_koszul;
use crate::rv_t_transfer::{self, Operation};
use crate::syntax::Label;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct PathId(pub(crate) usize);

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) path: PathId,
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) b: usize,
}

impl Grade {
    pub(crate) fn differential(self) -> Result<Self> {
        let b = self
            .b
            .checked_sub(1)
            .ok_or_else(|| Error::plain("rvBarD requires positive suspended homological degree"))?;
        Ok(Self { b, ..self })
    }

    pub(crate) fn max_length(self) -> usize {
        (self.end - self.start).min(self.b)
    }
}

#[derive(Clone, Copy, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Factor {
    pub(crate) end: usize,
    pub(crate) s: usize,
    pub(crate) generator: usize,
}

pub(crate) type Word = Rc<[Factor]>;
pub(crate) type Support = Rc<[Word]>;

type Interval = (usize, usize);
type Multiplication = (usize, Word);

#[derive(Clone, Copy)]
struct Summand {
    s: usize,
    rank: usize,
}

pub(crate) struct Path {
    name: Rc<str>,
    ideal: Rc<Ideal>,
    prefixes: Vec<(Vec<u128>, u128)>,
    homology: RefCell<BTreeMap<Interval, Rc<[Summand]>>>,
    products: RefCell<BTreeMap<Multiplication, Rc<[usize]>>>,
}

impl Path {
    pub(crate) fn build(model: &mut Model, label: Label, steps: Vec<Vec<usize>>) -> Result<Self> {
        let ideal = model.rv_ideal(&label)?;
        if rv_koszul::monomial_initial(model, &ideal)?.is_none() {
            return Err(Error::plain("rvbarpath requires a monomial rank ideal"));
        }
        let width = ideal.selected.len();
        let cells = (2 * width + 2) * (steps.len() + 1);
        check_storage(cells)?;
        model.spend(steps.len() * (1 + width))?;
        let mut prefixes = vec![(vec![0u128; width], 0u128)];
        for step in steps {
            if step.len() != width {
                return Err(Error::plain(
                    "rvbarpath requires one exponent per selected primitive",
                ));
            }
            if step.iter().all(|&exponent| exponent == 0) {
                return Err(Error::plain(
                    "rvbarpath steps must have nonzero multidegree",
                ));
            }
            let (prior, t) = prefixes.last().unwrap();
            let total = t + step.iter().map(|&x| x as u128).sum::<u128>();
            let next = prior
                .iter()
                .zip(step)
                .map(|(&x, y)| x + y as u128)
                .collect();
            prefixes.push((next, total));
        }
        model.allocate(cells)?;
        Ok(Self {
            name: label.text.into(),
            ideal,
            prefixes,
            homology: RefCell::new(BTreeMap::new()),
            products: RefCell::new(BTreeMap::new()),
        })
    }

    pub(crate) fn check_interval(&self, start: usize, end: usize) -> Result<()> {
        if start >= end || end >= self.prefixes.len() {
            return Err(Error::plain(format!(
                "invalid ordered bar interval [{start},{end}]"
            )));
        }
        Ok(())
    }

    pub(crate) fn tor_grade(&self, start: usize, end: usize, s: usize) -> Result<rv_koszul::Grade> {
        self.check_interval(start, end)?;
        let (before, t0) = &self.prefixes[start];
        let (after, t1) = &self.prefixes[end];
        let degree_error = || Error::plain("bar interval multidegree is out of range");
        let t = usize::try_from(t1 - t0).map_err(|_| degree_error())?;
        let multidegree = after
            .iter()
            .zip(before)
            .map(|(a, b)| usize::try_from(a - b).map_err(|_| degree_error()))
            .collect::<Result<Rc<[usize]>>>()?;
        Ok(rv_koszul::Grade {
            ideal: self.name.clone(),
            s,
            t,
            multidegree: Some(multidegree),
        })
    }

    fn ranks(&self, model: &mut Model, start: usize, end: usize) -> Result<Rc<[Summand]>> {
        if let Some(result) = self.homology.borrow().get(&(start, end)).cloned() {
            model.spend(1)?;
            return Ok(result);
        }
        let grade = self.tor_grade(start, end, 0)?;
        let multidegree = grade.multidegree.as_ref().unwrap();
        model.spend(1 + multidegree.len())?;
        let result = if let Some(masks) = self.ideal.squarefree_masks(model)? {
            // In a repeated coordinate support is unchanged when its exterior
            // factor is toggled, giving an elementary contraction.
            if multidegree.iter().any(|&power| power > 1) || masks.is_empty() {
                Vec::new()
            } else {
                let vertices = rv_ideal::support_mask(multidegree, 1);
                squarefree_ranks(model, masks, vertices)?
            }
        } else {
            let mut result = Vec::new();
            for s in 0..=multidegree.iter().filter(|&&power| power != 0).count() {
                let mut grade = grade.clone();
                grade.s = s;
                let ty = model.rv_complex_type(grade, true)?;
                let rank = model.rv_tor(ty)?.rank();
                if rank != 0 {
                    result.push(Summand { s, rank });
                }
            }
            result
        };
        model.allocate(3 + 2 * result.len())?;
        let result: Rc<[Summand]> = result.into();
        self.homology
            .borrow_mut()
            .insert((start, end), result.clone());
        Ok(result)
    }

    pub(crate) fn multiply(
        &self,
        model: &mut Model,
        start: usize,
        word: &[Factor],
        s: usize,
    ) -> Result<Rc<[usize]>> {
        model.spend(1 + word.len())?;
        let terms: Word = word.into();
        let key = (start, terms);
        if let Some(result) = self.products.borrow().get(&key).cloned() {
            return Ok(result);
        }
        let end = word.last().unwrap().end;
        let ranks = self.ranks(model, start, end)?;
        if !ranks.iter().any(|rank| rank.s == s) {
            return Ok([].into());
        }
        let target = model.rv_complex_type(self.tor_grade(start, end, s)?, true)?;
        let mut inputs = Vec::with_capacity(word.len());
        let mut left = start;
        for factor in word {
            let ty = model.rv_complex_type(self.tor_grade(left, factor.end, factor.s)?, true)?;
            inputs.push((ty, Rc::from([factor.generator])));
            left = factor.end;
        }
        let Value::Class(result) =
            rv_t_transfer::evaluate(model, Operation::Higher, target, &inputs)?
        else {
            unreachable!();
        };
        model.allocate(1 + word_cells(word))?;
        self.products.borrow_mut().insert(key, result.clone());
        Ok(result)
    }
}

fn squarefree_ranks(model: &mut Model, masks: &[usize], vertices: usize) -> Result<Vec<Summand>> {
    if !rv_ideal::standard_support(model, masks, 0)? {
        return Ok(Vec::new());
    }
    let size = vertices.count_ones() as usize;
    let mut reduced = vertices;
    for &mask in masks {
        model.spend(1)?;
        if mask.count_ones() == 1 {
            reduced &= !mask;
        }
    }
    if masks.iter().all(|mask| mask.count_ones() <= 2) {
        reduced = collapse_flag(model, masks, reduced)?;
    }
    let mut faces = vec![Vec::new(); reduced.count_ones() as usize + 1];
    let mut count = 0;
    collect_faces(model, masks, reduced, 0, 0, &mut faces, &mut count)?;
    let mut incoming = 0;
    let mut result = Vec::new();
    for (q, layer) in faces.into_iter().enumerate() {
        let mut differential = Basis::default();
        let mut outgoing = 0;
        for face in &layer {
            let mut image = Vec::new();
            let mut possible = reduced & !face;
            while possible != 0 {
                let bit = 1usize << possible.trailing_zeros();
                possible ^= bit;
                let extension = face | bit;
                if rv_ideal::standard_support(model, masks, extension)? {
                    image.push(extension);
                }
            }
            outgoing += usize::from(differential.insert(model, image)?);
        }
        let rank = layer
            .len()
            .checked_sub(incoming + outgoing)
            .ok_or_else(|| Error::plain("invalid augmented simplicial ranks"))?;
        if rank != 0 {
            result.push(Summand { s: size - q, rank });
        }
        incoming = outgoing;
    }
    result.reverse();
    Ok(result)
}

fn collapse_flag(model: &mut Model, masks: &[usize], mut vertices: usize) -> Result<usize> {
    let mut neighborhoods = vec![vertices; usize::BITS as usize];
    for &mask in masks {
        model.spend(1)?;
        if mask.count_ones() == 2 && (mask & vertices) == mask {
            let lower = mask.trailing_zeros() as usize;
            let upper = mask.ilog2() as usize;
            neighborhoods[lower] &= !(1usize << upper);
            neighborhoods[upper] &= !(1usize << lower);
        }
    }
    loop {
        let mut candidates = vertices;
        let mut removed = 0;
        'search: while candidates != 0 {
            model.spend(1)?;
            let index = candidates.trailing_zeros() as usize;
            let vertex = 1usize << index;
            candidates ^= vertex;
            let closed = neighborhoods[index] & vertices;
            let mut neighbors = closed & !vertex;
            while neighbors != 0 {
                model.spend(1)?;
                let neighbor = neighbors.trailing_zeros() as usize;
                neighbors ^= 1usize << neighbor;
                if closed & !neighborhoods[neighbor] == 0 {
                    // In the clique complex the vertex is dominated by this
                    // neighbor. Deleting it preserves reduced cohomology.
                    removed = vertex;
                    break 'search;
                }
            }
        }
        if removed == 0 {
            return Ok(vertices);
        }
        vertices ^= removed;
    }
}

fn collect_faces(
    model: &mut Model,
    masks: &[usize],
    mut options: usize,
    face: usize,
    q: usize,
    layers: &mut [Vec<usize>],
    count: &mut usize,
) -> Result<()> {
    model.spend(1)?;
    *count += 1;
    check_storage(*count)?;
    layers[q].push(face);
    while options != 0 {
        let bit = 1usize << options.trailing_zeros();
        options ^= bit;
        let extended = face | bit;
        if rv_ideal::standard_support(model, masks, extended)? {
            collect_faces(model, masks, options, extended, q + 1, layers, count)?;
        }
    }
    Ok(())
}

fn check_storage(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "rank bar calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

pub(crate) fn word_cells(word: &[Factor]) -> usize {
    1 + 3 * word.len()
}

fn support_cells(words: &[Word]) -> usize {
    1 + words.iter().map(|word| word_cells(word)).sum::<usize>()
}

pub(crate) fn persist(model: &mut Model, words: Vec<Word>) -> Result<Support> {
    model.allocate(support_cells(&words))?;
    Ok(words.into())
}

pub(crate) fn compare(model: &mut Model, left: &[Word], right: &[Word]) -> Result<Ordering> {
    model.spend(support_cells(left) + support_cells(right))?;
    Ok(left.iter().rev().cmp(right.iter().rev()))
}

pub(crate) fn add(model: &mut Model, left: &[Word], right: &[Word]) -> Result<Support> {
    model.spend(support_cells(left) + support_cells(right))?;
    persist(model, xor_sorted(left, right))
}

pub(crate) fn length_part(model: &mut Model, length: usize, support: &[Word]) -> Result<Support> {
    model.spend(1 + support.len())?;
    persist(
        model,
        support
            .iter()
            .filter(|word| word.len() == length)
            .cloned()
            .collect(),
    )
}

pub(crate) fn word_type(
    model: &mut Model,
    path_id: PathId,
    cuts: &[usize],
    inputs: impl ExactSizeIterator<Item = TypeId>,
) -> Result<(Grade, Vec<usize>)> {
    if inputs.len() == 0 || cuts.len() != inputs.len() + 1 {
        return Err(Error::plain(
            "rvBarWord requires one class per cut interval",
        ));
    }
    let path = model.rv_bar_path(path_id);
    let mut suspended = 0usize;
    let mut degrees = Vec::with_capacity(inputs.len());
    for (ends, ty) in cuts.windows(2).zip(inputs) {
        let Kind::RVTor(input) = model.kind(ty) else {
            return Err(Error::plain("rvBarWord requires RVMTor classes"));
        };
        let expected = path.tor_grade(ends[0], ends[1], input.s)?;
        if input != &expected {
            return Err(Error::plain(format!(
                "rvBarWord class has the wrong ideal or multidegree on interval [{},{}]",
                ends[0], ends[1]
            )));
        }
        degrees.push(input.s);
        suspended = suspended
            .checked_add(input.s)
            .and_then(|b| b.checked_add(1))
            .ok_or_else(|| Error::plain("bar suspended homological degree is out of range"))?;
    }
    let grade = Grade {
        path: path_id,
        start: cuts[0],
        end: *cuts.last().unwrap(),
        b: suspended,
    };
    Ok((grade, degrees))
}

pub(crate) fn word(
    model: &mut Model,
    ends: &[usize],
    degrees: &[usize],
    classes: &[Rc<[usize]>],
) -> Result<Support> {
    if classes.iter().any(|x| x.is_empty()) {
        return Ok([].into());
    }
    let cost = 1 + 3 * classes.len();
    let count = classes
        .iter()
        .try_fold(1usize, |n, class| n.checked_mul(class.len()))
        .ok_or_else(|| Error::plain("bar word expansion exceeds the storage limit"))?;
    check_storage(count.saturating_mul(cost))?;
    model.spend(count * cost)?;
    let mut support = Vec::with_capacity(count);
    for mut code in 0..count {
        let mut word = Vec::with_capacity(classes.len());
        for ((&end, &s), class) in ends.iter().zip(degrees).zip(classes).rev() {
            word.push(Factor {
                end,
                s,
                generator: class[code % class.len()],
            });
            code /= class.len();
        }
        word.reverse();
        support.push(word.into());
    }
    persist(model, support)
}

pub(crate) fn differential(model: &mut Model, grade: Grade, terms: &[Word]) -> Result<Vec<Word>> {
    let path = model.rv_bar_path(grade.path);
    let mut sum = BTreeSet::new();
    let mut cells = 0;
    for word in terms {
        let mut left = grade.start;
        for i in 0..word.len() {
            let mut s = word[i].s;
            for j in i + 1..word.len() {
                model.spend(1)?;
                s += word[j].s;
                let combined_s = s + j - i - 1;
                let product = path.multiply(model, left, &word[i..=j], combined_s)?;
                for &generator in product.iter() {
                    model.spend(word_cells(word) * (1 + sum.len().max(1).ilog2() as usize))?;
                    let mut combined = Vec::with_capacity(word.len() - (j - i));
                    combined.extend_from_slice(&word[..i]);
                    combined.push(Factor {
                        end: word[j].end,
                        s: combined_s,
                        generator,
                    });
                    combined.extend_from_slice(&word[j + 1..]);
                    let term: Word = combined.into();
                    if sum.remove(&term) {
                        cells -= word_cells(&term);
                    } else {
                        cells += word_cells(&term);
                        check_storage(cells)?;
                        sum.insert(term);
                    }
                }
            }
            left = word[i].end;
        }
    }
    Ok(sum.into_iter().collect())
}

type State = (usize, usize); // suspended degree and length
type Counts = BTreeMap<State, Option<u128>>;

struct BasisPlan {
    grade: Grade,
    suffix: Vec<Counts>,
    path: Rc<Path>,
}

impl BasisPlan {
    fn build(model: &mut Model, grade: Grade, max_length: usize) -> Result<Self> {
        let path = model.rv_bar_path(grade.path);
        let span = grade.end - grade.start;
        let max_length = max_length.min(grade.max_length());
        let internal = path.prefixes[grade.end].1 - path.prefixes[grade.start].1;
        if max_length == 0
            || grade.b > max_length.saturating_mul(path.ideal.selected.len() + 1)
            || grade.b as u128 > internal + max_length as u128
        {
            return Ok(Self {
                grade,
                suffix: Vec::new(),
                path,
            });
        }
        let mut cells = 3 * (span + 1) + 4;
        check_storage(cells)?;
        let mut suffix = vec![Counts::new(); span + 1];
        suffix[span].insert((0, 0), Some(1));
        for position in (0..span).rev() {
            let left = grade.start + position;
            let (head, tail) = suffix.split_at_mut(position + 1);
            let states = &mut head[position];
            for (offset, endings) in tail.iter().enumerate() {
                model.spend(1)?;
                if endings.is_empty() {
                    continue;
                }
                let right = left + offset + 1;
                let ranks = path.ranks(model, left, right)?;
                for &Summand { s, rank } in ranks.iter() {
                    for (&(b, length), &count) in endings {
                        model.spend(1)?;
                        if length >= max_length || s + 1 > grade.b - b {
                            continue;
                        }
                        let key = (b + s + 1, length + 1);
                        // Mark overflow in an intermediate degree; it can
                        // be irrelevant to the final requested slice.
                        let contribution = count.and_then(|count| count.checked_mul(rank as u128));
                        let log = 1 + states.len().max(1).ilog2() as usize;
                        model.spend(log)?;
                        match states.entry(key) {
                            std::collections::btree_map::Entry::Occupied(mut value) => {
                                *value.get_mut() = value
                                    .get()
                                    .zip(contribution)
                                    .and_then(|(a, b)| a.checked_add(b));
                            }
                            std::collections::btree_map::Entry::Vacant(slot) => {
                                cells += 4;
                                check_storage(cells)?;
                                slot.insert(contribution);
                            }
                        }
                    }
                }
            }
        }
        Ok(Self {
            grade,
            suffix,
            path,
        })
    }

    fn starts(&self) -> impl Iterator<Item = (usize, Option<u128>)> + '_ {
        self.suffix.first().into_iter().flat_map(|states| {
            states
                .range((self.grade.b, 0)..=(self.grade.b, usize::MAX))
                .map(|(&(_, length), &count)| (length, count))
        })
    }

    fn count(&self, keep: impl Fn(usize) -> bool) -> Result<u128> {
        self.starts()
            .filter(|&(k, _)| keep(k))
            .try_fold(0u128, |sum, (_, count)| {
                sum.checked_add(count.ok_or_else(count_limit)?)
                    .ok_or_else(count_limit)
            })
    }

    fn words(&self, model: &mut Model, minimum_length: usize, limit: usize) -> Result<Vec<Word>> {
        let count = self.count(|length| length >= minimum_length)?;
        if count > limit as u128 {
            return Err(Error::plain("rank bar basis exceeds the enumeration limit"));
        }
        let grade = self.grade;
        let mut stack = Vec::new();
        for (length, _) in self
            .starts()
            .filter(|&(length, _)| length >= minimum_length)
        {
            stack.push(PartialWord {
                position: grade.start,
                b: grade.b,
                length,
                factors: Vec::new(),
            });
        }
        let mut retained = 4 * stack.len();
        let mut output = Vec::with_capacity(count as usize);
        while let Some(state) = stack.pop() {
            retained -= 4 + 3 * state.factors.len();
            if state.position == grade.end {
                retained += word_cells(&state.factors);
                check_storage(retained)?;
                output.push(state.factors.into());
                continue;
            }
            for next in state.position + 1..=grade.end {
                model.spend(1)?;
                let endings = &self.suffix[next - grade.start];
                if endings.is_empty() {
                    continue;
                }
                let ranks = self.path.ranks(model, state.position, next)?;
                for &Summand { s, rank } in ranks.iter() {
                    let Some(b) = state.b.checked_sub(s + 1) else {
                        continue;
                    };
                    let length = state.length - 1;
                    model.spend(1 + endings.len().ilog2() as usize)?;
                    if !endings.contains_key(&(b, length)) {
                        continue;
                    }
                    for generator in 0..rank {
                        let cost = 4 + 3 * (state.factors.len() + 1);
                        model.spend(cost)?;
                        retained += cost;
                        check_storage(retained)?;
                        let mut factors = state.factors.clone();
                        factors.push(Factor {
                            end: next,
                            s,
                            generator,
                        });
                        stack.push(PartialWord {
                            position: next,
                            b,
                            length,
                            factors,
                        });
                    }
                }
            }
        }
        model.spend(retained * (1 + output.len().max(1).ilog2() as usize))?;
        output.sort_unstable();
        Ok(output)
    }
}

struct PartialWord {
    position: usize,
    b: usize,
    length: usize,
    factors: Vec<Factor>,
}

fn count_limit() -> Error {
    Error::plain("rank bar dimension is out of range")
}

pub(crate) fn dimension(model: &mut Model, grade: Grade, length: Option<usize>) -> Result<u128> {
    let max = length.unwrap_or_else(|| grade.max_length());
    BasisPlan::build(model, grade, max)?.count(|k| length.is_none_or(|length| length == k))
}

pub(crate) fn basis(
    model: &mut Model,
    grade: Grade,
    max_length: usize,
    limit: usize,
) -> Result<Vec<Word>> {
    BasisPlan::build(model, grade, max_length)?.words(model, 1, limit)
}

pub(crate) fn replacements(model: &mut Model, grade: Grade, limit: usize) -> Result<Vec<Word>> {
    BasisPlan::build(model, grade, grade.max_length())?.words(model, 2, limit)
}
