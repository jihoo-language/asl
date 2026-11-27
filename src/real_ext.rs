use std::collections::BTreeSet;
use std::hash::Hash;

use crate::adams;
use crate::cobar::{self, Monomial, Word};
use crate::cohomology;
use crate::error::{Error, Result};
use crate::linear::{Basis, QuotientBasis};
use crate::may::Rows;
use crate::model::{Kind, Model, TypeId, Value};
use crate::motivic::{self, Grade};
use crate::real_bars;
use crate::real_motivic::{self, Term};
use crate::rho_localized;

const MAX_BASIS: usize = 32_768;

#[derive(Default)]
pub(crate) struct Group {
    terms: Vec<Term>,
    quotient: QuotientBasis,
}

impl Group {
    pub(crate) fn build(
        model: &mut Model,
        grade: Grade,
        cutoff: Option<usize>,
        bars: Option<&real_bars::Group>,
    ) -> Result<Self> {
        if let Some(cycles) = degree_one_cycles(model, grade, cutoff)? {
            return Self::from_cycle_supports(model, cycles);
        }
        if bars.is_none() {
            if let Some(cycles) = filtered_degree_one_cycles(model, grade, cutoff)? {
                return Self::from_cycle_supports(model, cycles);
            }
        }
        Self::from_matrix(model, grade, cutoff, bars)
    }

    fn from_matrix(
        model: &mut Model,
        grade: Grade,
        cutoff: Option<usize>,
        bars: Option<&real_bars::Group>,
    ) -> Result<Self> {
        let max = cutoff.map_or(usize::MAX, |exponent| exponent - 1);
        let bound = cutoff.map_or(u128::MAX, |exponent| exponent as u128);
        let terms = real_motivic::basis_window(model, grade, 0..=max, MAX_BASIS)?;
        if terms.is_empty() {
            return Ok(Self::default());
        }
        let mut boundaries = Basis::default();
        if let Some(previous) = grade.previous() {
            for term in real_motivic::basis_window(model, previous, 0..=max, MAX_BASIS)? {
                let image = real_motivic::differential_below(model, &[term], bound)?;
                let column = indices(model, &terms, &image, wrong_degree)?;
                boundaries.insert(model, column)?;
            }
        }
        if let Some(bars) = bars {
            // The cached barcode supplies a basis of ordinary cohomology.
            // Reduce its representatives into the same boundary complement.
            let mut span = Basis::default();
            for index in bars.living(model, grade.t)? {
                let cycle = bars.get(index)?.lift_at(model, grade.t)?;
                let support = indices(model, &terms, &cycle, wrong_degree)?;
                let reduced = boundaries.reduce(model, support)?;
                if !span.insert(model, reduced)? {
                    return Err(Error::plain("living rho bars are linearly dependent"));
                }
            }
            let quotient = QuotientBasis::from_cycle_span(model, boundaries, span)?;
            return Ok(Self { terms, quotient });
        }
        let mut rows = real_motivic::Rows::default();
        let quotient = QuotientBasis::build(model, boundaries, terms.len(), |model, index| {
            let image = real_motivic::differential_below(
                model,
                std::slice::from_ref(&terms[index]),
                bound,
            )?;
            rows.column(model, &image)
        })?;
        Ok(Self { terms, quotient })
    }

    fn from_cycle_supports(model: &mut Model, cycles: Vec<Vec<Term>>) -> Result<Self> {
        let mut terms = BTreeSet::new();
        let mut cells = 0;
        for t in cycles.iter().flatten() {
            model.spend(t.cells() * (1 + terms.len().max(1).ilog2() as usize))?;
            if terms.insert(t.clone()) {
                cells += t.cells();
                cobar::check_storage(cells)?;
            }
        }
        model.allocate(cells)?;
        let terms: Vec<_> = terms.into_iter().collect();
        let mut span = Basis::default();
        for cycle in cycles {
            let column = indices(model, &terms, &cycle, wrong_degree)?;
            span.insert(model, column)?;
        }
        let quotient = QuotientBasis::from_cycle_span(model, Basis::default(), span)?;
        Ok(Self { terms, quotient })
    }

    pub(crate) fn rank(&self) -> usize {
        self.quotient.rank()
    }

    pub(crate) fn representative(&self, model: &mut Model, coords: &[usize]) -> Result<Vec<Term>> {
        Ok(self
            .quotient
            .representative(model, coords)?
            .into_iter()
            .map(|coordinate| self.terms[coordinate].clone())
            .collect())
    }

    pub(crate) fn project(&self, model: &mut Model, terms: &[Term]) -> Result<Vec<usize>> {
        let support = indices(model, &self.terms, terms, not_closed)?;
        self.quotient
            .project(model, support)?
            .ok_or_else(not_closed)
    }
}

fn degree_one_window(grade: Grade, cutoff: Option<usize>) -> Option<(i128, usize, usize)> {
    if grade.s != 1 || grade.t <= 0 {
        // At positive t there are no boundaries from scalar degree.
        return None;
    }
    let k = grade.t.checked_sub(i128::from(grade.w))?;
    let mut upper = k.checked_mul(2)?;
    if let Some(e) = cutoff {
        upper = upper.min(grade.t.saturating_add(e as i128 - 1));
    }
    let lower = usize::try_from(grade.t).ok()?;
    let upper = usize::try_from(upper.max(0)).ok()?;
    Some((k, lower, upper))
}

fn scalar_degree_one_cycle(
    model: &mut Model,
    top: usize,
    lower: usize,
    bound: u128,
) -> Result<Vec<Term>> {
    let scalar = Term {
        rho: 0,
        tau: top,
        word: [].into(),
    };
    let mut cycle = real_motivic::differential_below(model, &[scalar], lower as u128 + bound)?;
    for term in &mut cycle {
        term.rho = term
            .rho
            .checked_sub(lower)
            .ok_or_else(|| Error::plain("rho primitive has insufficient filtration"))?;
    }
    if cycle.first().map(|term| term.rho) != Some(top - lower) {
        return Err(Error::plain("rho primitive has an unexpected leading term"));
    }
    Ok(cycle)
}

fn degree_one_cycles(
    model: &mut Model,
    grade: Grade,
    cutoff: Option<usize>,
) -> Result<Option<Vec<Vec<Term>>>> {
    let Some((k, lower, upper)) = degree_one_window(grade, cutoff) else {
        return Ok(None);
    };
    if upper < lower {
        return Ok(Some(Vec::new()));
    }
    let mut powers = Vec::new();
    let mut top = 1usize;
    while top <= upper {
        model.spend(1)?;
        if top >= lower {
            if top as i128 != k && top as i128 != 2 * k {
                return Ok(None);
            }
            powers.push(top);
        }
        let Some(next) = top.checked_mul(2) else {
            break;
        };
        top = next;
    }

    // In each rho level the leading cycle is a complex-motivic primitive;
    // degree-one primitives in the classical encoding are xi_1^(2^j).
    // When every possible leading term has a closed lift, they form the
    // full cycle space. At top=k use the actual d(tau^k)/rho^t, and at
    // top=2k use a power of the primitive real generator xi_1.
    let mut cycles = Vec::new();
    for top in powers {
        if top as i128 == k {
            let bound = (upper - lower) as u128 + 1;
            cycles.push(scalar_degree_one_cycle(model, top, lower, bound)?);
        } else {
            let factor: Monomial = [top].into();
            cycles.push(vec![Term {
                rho: top - lower,
                tau: 0,
                word: [factor].into(),
            }]);
        }
    }
    Ok(Some(cycles))
}

type FilteredDirection = (Vec<Term>, Vec<Term>);

#[derive(Default)]
struct FilteredKernel {
    cycles: Vec<Vec<Term>>,
    pending: Vec<FilteredDirection>,
    cells: usize,
}

impl FilteredKernel {
    fn retain(&mut self, model: &mut Model, cochain: &[Term], image: &[Term]) -> Result<()> {
        model.spend(cochain.len() + image.len())?;
        self.cells += degree_one_cells(cochain, image);
        cobar::check_storage(self.cells)
    }

    fn insert(&mut self, model: &mut Model, cochain: Vec<Term>, image: Vec<Term>) -> Result<()> {
        self.retain(model, &cochain, &image)?;
        if image.is_empty() {
            self.cycles.push(cochain);
        } else {
            self.pending.push((cochain, image));
        }
        Ok(())
    }
}

fn degree_one_cells(cochain: &[Term], image: &[Term]) -> usize {
    cochain.iter().chain(image).map(Term::cells).sum()
}

fn filtered_degree_one_cycles(
    model: &mut Model,
    grade: Grade,
    cutoff: Option<usize>,
) -> Result<Option<Vec<Vec<Term>>>> {
    let Some((k, lower, upper)) = degree_one_window(grade, cutoff) else {
        return Ok(None);
    };
    if upper < lower {
        return Ok(Some(Vec::new()));
    }
    let bound = (upper - lower) as u128 + 1;
    let mut kernel = FilteredKernel::default();
    let mut top = 1usize;
    while top <= upper {
        model.spend(1)?;
        if top >= lower {
            if top as i128 == k {
                let cycle = scalar_degree_one_cycle(model, top, lower, bound)?;
                kernel.insert(model, cycle, Vec::new())?;
            } else {
                let tau = usize::try_from(k - top.div_ceil(2) as i128)
                    .map_err(|_| Error::plain("degree-one primitive weight exceeds its range"))?;
                let factor: Monomial = [top].into();
                let cochain = vec![Term {
                    rho: top - lower,
                    tau,
                    word: [factor].into(),
                }];
                let image = real_motivic::differential_below(model, &cochain, bound)?;
                kernel.insert(model, cochain, image)?;
            }
        }
        let Some(next) = top.checked_mul(2) else {
            break;
        };
        top = next;
    }

    // Every cycle's first term in the rho filtration is a d0-primitive.
    // Correct the first obstruction by a d0-preimage, then retain precisely
    // the combinations whose obstruction class cancels at that level.
    while let Some(level) = kernel.pending.iter().map(|(_, image)| image[0].rho).min() {
        model.spend(kernel.pending.len())?;
        let (front, later): (Vec<_>, Vec<_>) = kernel
            .pending
            .into_iter()
            .partition(|(_, image)| image[0].rho == level);
        kernel.pending = later;
        resolve_degree_one_layer(model, grade, level, bound, front, &mut kernel)?;
    }
    Ok(Some(kernel.cycles))
}

fn resolve_degree_one_layer(
    model: &mut Model,
    grade: Grade,
    level: usize,
    bound: u128,
    front: Vec<FilteredDirection>,
    kernel: &mut FilteredKernel,
) -> Result<()> {
    let source = real_motivic::basis_window(model, grade, level..=level, MAX_BASIS)?;
    let mut rows = real_motivic::Rows::default();
    let mut boundaries = Basis::default();
    for (index, term) in source.iter().enumerate() {
        let image =
            real_motivic::differential_below(model, std::slice::from_ref(term), level as u128 + 1)?;
        let column = rows.column(model, &image)?;
        if !column.is_empty() {
            boundaries.insert_with_witness(model, column, vec![index])?;
        }
    }

    let mut directions = Vec::new();
    let mut obstructions = Basis::default();
    for (mut cochain, mut image) in front {
        kernel.cells -= degree_one_cells(&cochain, &image);
        let current = &image[..image.partition_point(|term| term.rho == level)];
        let encoded = rows.column(model, current)?;
        // Remove the part in the d0 image using an actual cochain at this level.
        let (remainder, correction) = boundaries.reduce_with_witness(model, encoded, Vec::new())?;
        if !correction.is_empty() {
            let correction: Vec<_> = correction.into_iter().map(|i| source[i].clone()).collect();
            let delta = real_motivic::differential_below(model, &correction, bound)?;
            cochain = real_motivic::add(model, &cochain, &correction)?;
            image = real_motivic::add(model, &image, &delta)?;
        }
        let index = directions.len();
        kernel.retain(model, &cochain, &image)?;
        directions.push((cochain, image));
        // The kernel of the remaining obstruction matrix is exactly the
        // space of combinations that can persist to a higher rho level.
        if let Some(relation) = obstructions.insert_with_witness(model, remainder, vec![index])? {
            let (mut combined, mut image) = (Vec::new(), Vec::new());
            for index in relation {
                let (c, d) = &directions[index];
                combined = real_motivic::add(model, &combined, c)?;
                image = real_motivic::add(model, &image, d)?;
            }
            if image.first().is_some_and(|first| first.rho <= level) {
                return Err(Error::plain(
                    "degree-one obstruction did not cancel at its rho level",
                ));
            }
            kernel.insert(model, combined, image)?;
        }
    }
    for (cochain, image) in directions {
        kernel.cells -= degree_one_cells(&cochain, &image);
    }
    Ok(())
}

fn not_closed() -> Error {
    Error::plain("operation requires a closed real-motivic cobar cochain")
}

fn wrong_degree() -> Error {
    Error::plain("word is outside its real-motivic cobar tridegree")
}

fn indices(
    model: &mut Model,
    basis: &[Term],
    terms: &[Term],
    missing: fn() -> Error,
) -> Result<Vec<usize>> {
    let steps = 1 + basis.len().max(1).ilog2() as usize;
    terms
        .iter()
        .map(|term| {
            model.spend(steps * term.cells())?;
            basis.binary_search(term).map_err(|_| missing())
        })
        .collect()
}

fn column<T: Hash + Eq>(
    model: &mut Model,
    rows: &mut Rows<T>,
    support: impl IntoIterator<Item = T>,
    cells: fn(&T) -> usize,
) -> Result<Vec<usize>> {
    let mut column = support
        .into_iter()
        .map(|term| rows.coordinate(model, cells(&term), term))
        .collect::<Result<Vec<_>>>()?;
    model.spend(column.len() * (1 + column.len().max(1).ilog2() as usize))?;
    column.sort_unstable();
    Ok(column)
}

pub(crate) fn class(model: &mut Model, ty: TypeId, terms: &[Term]) -> Result<Value> {
    if terms.is_empty() {
        return Ok(cohomology::zero());
    }
    let coords = model.real_ext(ty)?.project(model, terms)?;
    model.allocate(coords.len())?;
    Ok(Value::Class(coords.into()))
}

pub(crate) fn representative(model: &mut Model, input: TypeId, coords: &[usize]) -> Result<Value> {
    if coords.is_empty() {
        return Ok(real_motivic::zero());
    }
    let terms = model.real_ext(input)?.representative(model, coords)?;
    real_motivic::persist(model, terms)
}

pub(crate) fn yoneda(
    model: &mut Model,
    left: (TypeId, &[usize]),
    right: (TypeId, &[usize]),
    output: TypeId,
) -> Result<Value> {
    if left.1.is_empty() || right.1.is_empty() {
        return Ok(cohomology::zero());
    }
    let a = model.real_ext(left.0)?.representative(model, left.1)?;
    let b = model.real_ext(right.0)?.representative(model, right.1)?;
    let terms = real_motivic::concat_support(model, &a, &b)?;
    class(model, output, &terms)
}

pub(crate) fn shift(
    model: &mut Model,
    input: TypeId,
    coords: &[usize],
    output: TypeId,
    power: usize,
) -> Result<Value> {
    if coords.is_empty() {
        return Ok(cohomology::zero());
    }
    let lift = model.real_ext(input)?.representative(model, coords)?;
    let terms = real_motivic::shift_support(model, &lift, power)?;
    let Kind::RMExt(grade) = *model.kind(output) else {
        unreachable!();
    };
    if grade.s == 1
        && grade.t <= 0
        && !model.has_real_ext(output)
        && real_preimage(model, grade, &terms, false)?.is_some()
    {
        // A scalar primitive already certifies a zero in the target group.
        return Ok(cohomology::zero());
    }
    class(model, output, &terms)
}

pub(crate) fn specialize(
    model: &mut Model,
    input: TypeId,
    coords: &[usize],
    output: TypeId,
) -> Result<Value> {
    if coords.is_empty() {
        return Ok(cohomology::zero());
    }
    let lift = model.real_ext(input)?.representative(model, coords)?;
    let Value::Cobar(words) = real_motivic::specialize(model, &lift)? else {
        unreachable!();
    };
    adams::class(model, output, &words)
}

pub(crate) fn project_classical(
    model: &mut Model,
    input: TypeId,
    coords: &[usize],
    output: TypeId,
) -> Result<Value> {
    if coords.is_empty() {
        return Ok(cohomology::zero());
    }
    let lift = model.real_ext(input)?.representative(model, coords)?;
    let words = classical_words(model, &lift)?;
    adams::class(model, output, &words)
}

fn classical_words(model: &mut Model, terms: &[Term]) -> Result<Vec<Word>> {
    let mut projected = cobar::Support::new();
    for term in terms {
        model.spend(term.cells())?;
        // Each tau_i has zero image after killing the vertex and the U
        // variables. The remaining factors are the pure classical xi words.
        if term.tau != 0
            || term
                .word
                .iter()
                .flat_map(|factor| factor.iter())
                .any(|&power| power & 1 != 0)
        {
            continue;
        }
        let word: Word = term
            .word
            .iter()
            .map(|factor| factor.iter().map(|power| power / 2).collect::<Monomial>())
            .collect();
        let cells = cobar::word_cells(&word);
        projected.toggle(model, word, cells)?;
    }
    Ok(projected.into_vec())
}

pub(crate) fn primitive(model: &mut Model, grade: Grade, terms: &[Term]) -> Result<Value> {
    let witness = real_preimage(model, grade, terms, true)?
        .ok_or_else(|| Error::plain("extPrimitive requires a real-motivic cobar boundary"))?;
    real_motivic::persist(model, witness)
}

// Restrict d to the span of the first linearly independent basis columns.
// Its inverse on the image is a linear section, independent of the order
// in which the sparse target rows are discovered.
pub(crate) fn linear_primitive(
    model: &mut Model,
    grade: Grade,
    terms: &[Term],
) -> Result<Option<Vec<Term>>> {
    if terms.is_empty() {
        return Ok(Some(Vec::new()));
    }
    let Some(previous) = grade.previous() else {
        return Ok(None);
    };
    ordered_real_preimage(model, previous, terms, true)
}

// Check the actual incoming differential. Stopping once the target
// reduces to zero avoids reducing unrelated later boundary columns.
fn preimage<T: Clone + Eq + Hash>(
    model: &mut Model,
    source: &[T],
    terms: &[T],
    cells: fn(&T) -> usize,
    differential: impl FnMut(&mut Model, &T) -> Result<Vec<T>>,
    witness: bool,
) -> Result<Option<Vec<T>>> {
    let mut rows = Rows::new("cobar boundary");
    sparse_preimage(
        model,
        source,
        terms,
        |model, terms| column(model, &mut rows, terms.iter().cloned(), cells),
        differential,
        witness,
    )
}

fn sparse_preimage<T: Clone>(
    model: &mut Model,
    source: &[T],
    terms: &[T],
    mut encode: impl FnMut(&mut Model, &[T]) -> Result<Vec<usize>>,
    mut differential: impl FnMut(&mut Model, &T) -> Result<Vec<T>>,
    witness: bool,
) -> Result<Option<Vec<T>>> {
    let mut remainder = encode(model, terms)?;
    let mut coefficients = Vec::new();
    let mut span = Basis::default();
    for (index, term) in source.iter().enumerate() {
        let image = differential(model, term)?;
        let image = encode(model, &image)?;
        let incoming = if witness { vec![index] } else { Vec::new() };
        if span.insert_with_witness(model, image, incoming)?.is_none()
            && remainder.last().is_some_and(|&j| span.contains_pivot(j))
        {
            (remainder, coefficients) = span.reduce_with_witness(model, remainder, coefficients)?;
            if remainder.is_empty() {
                return Ok(Some(
                    coefficients
                        .into_iter()
                        .map(|i| source[i].clone())
                        .collect(),
                ));
            }
        }
    }
    Ok(None)
}

pub(crate) fn degree_one_beta_witness(
    model: &mut Model,
    grade: Grade,
    exponent: usize,
    visible: &[Term],
) -> Result<Option<Vec<Term>>> {
    // A closed lift modulo rho^(2e) suffices. These degree-one cycles
    // have distinct leading rho levels, and there are no incoming
    // scalar boundaries at positive topological degree.
    let Some(cycles) = degree_one_cycles(model, grade, exponent.checked_mul(2))? else {
        return Ok(None);
    };
    let mut remainder = visible.to_vec();
    let mut tail = Vec::new();
    for cycle in cycles {
        let split = cycle.partition_point(|term| term.rho < exponent);
        let (visible_part, higher_part) = cycle.split_at(split);
        if !visible_part.is_empty() && visible_part.first() == remainder.first() {
            remainder = real_motivic::add(model, &remainder, visible_part)?;
            tail = real_motivic::add(model, &tail, higher_part)?;
        }
    }
    if !remainder.is_empty() {
        return Ok(None);
    }
    // The visible part was equal term by term, so the difference
    // from the canonical lift is exactly the tail.
    model.spend(tail.len())?;
    for term in &mut tail {
        term.rho -= exponent;
    }
    tail.truncate(tail.partition_point(|term| term.rho < exponent));
    Ok(Some(tail))
}

pub(crate) fn truncated_preimage(
    model: &mut Model,
    grade: Grade,
    exponent: usize,
    terms: &[Term],
    alternate: impl FnOnce(&mut Model) -> Result<Option<Vec<Term>>>,
) -> Result<Option<Vec<Term>>> {
    if terms.is_empty() {
        return Ok(Some(Vec::new()));
    }
    let Some(previous) = grade.previous() else {
        return Ok(None);
    };
    // A normalized real term has rho exponent at most t - 2w.
    let width = grade
        .t
        .saturating_sub(2 * i128::from(grade.w))
        .saturating_add(1);
    let cutoff = exponent.min(usize::try_from(width.max(0)).unwrap_or(usize::MAX));
    if cutoff == 0 {
        return Ok(None);
    }
    let first = terms[0].rho;
    if first > 0 {
        // Solve the first nonzero associated-graded layer cheaply. Its
        // primitive is accepted only if it lifts against the full target.
        let front = &terms[..terms.partition_point(|term| term.rho == first)];
        let source = real_motivic::basis_window(model, previous, first..=first, MAX_BASIS)?;
        if let Some(witness) = truncated_columns(model, source, front, first + 1)? {
            let image = real_motivic::differential_below(model, &witness, cutoff as u128)?;
            if real_motivic::compare(model, &image, terms)?.is_eq() {
                return Ok(Some(witness));
            }
        }
    }
    if let Some(witness) = alternate(model)? {
        return Ok(Some(witness));
    }
    let source = real_motivic::basis_window(model, previous, 0..=cutoff - 1, MAX_BASIS)?;
    truncated_columns(model, source, terms, cutoff)
}

fn truncated_columns(
    model: &mut Model,
    mut source: Vec<Term>,
    terms: &[Term],
    cutoff: usize,
) -> Result<Option<Vec<Term>>> {
    if source.is_empty() {
        return Ok(None);
    }
    // Higher-filtration columns have fewer potentially nonzero rho layers.
    model.spend(source.len())?;
    source.reverse();
    let mut rows = real_motivic::Rows::default();
    let mut witness = sparse_preimage(
        model,
        &source,
        terms,
        |model, terms| rows.column(model, terms),
        |model, term| {
            real_motivic::differential_below(model, std::slice::from_ref(term), cutoff as u128)
        },
        true,
    )?;
    if let Some(witness) = &mut witness {
        model.spend(witness.len())?;
        witness.reverse();
    }
    Ok(witness)
}

pub(crate) fn real_preimage(
    model: &mut Model,
    grade: Grade,
    terms: &[Term],
    witness: bool,
) -> Result<Option<Vec<Term>>> {
    if terms.is_empty() {
        return Ok(Some(Vec::new()));
    }
    let Some(previous) = grade.previous() else {
        return Ok(None);
    };
    if previous.s == 0 {
        let source = real_motivic::basis(model, previous, 1)?;
        if let [scalar] = source.as_slice() {
            let image = real_motivic::differential_support(model, std::slice::from_ref(scalar))?;
            if real_motivic::compare(model, &image, terms)? == std::cmp::Ordering::Equal {
                return Ok(Some(if witness { source } else { Vec::new() }));
            }
        }
        return Ok(None);
    }

    if let Some(primitive) = left_tau_preimage(model, terms)? {
        return Ok(Some(if witness { primitive } else { Vec::new() }));
    }

    // A non-boundary after specialization is already a certificate
    // against a real primitive, without using the other rho levels.
    let complex = real_motivic::complex_support(model, terms)?;
    if !complex.is_empty() {
        let source = motivic::basis(model, previous, MAX_BASIS)?;
        if preimage(
            model,
            &source,
            &complex,
            cobar::word_cells,
            cobar::differential_word,
            false,
        )?
        .is_none()
        {
            return Ok(None);
        }
    }

    ordered_real_preimage(model, previous, terms, witness)
}

fn ordered_real_preimage(
    model: &mut Model,
    previous: Grade,
    terms: &[Term],
    witness: bool,
) -> Result<Option<Vec<Term>>> {
    let source = real_motivic::basis(model, previous, MAX_BASIS)?;
    preimage(
        model,
        &source,
        terms,
        Term::cells,
        |model, term| real_motivic::differential_support(model, std::slice::from_ref(term)),
        witness,
    )
}

fn left_tau_preimage(model: &mut Model, terms: &[Term]) -> Result<Option<Vec<Term>>> {
    let mut candidate = Vec::with_capacity(terms.len());
    for term in terms {
        model.spend(term.cells())?;
        let (Some(rho), Some(tau), Some((factor, tail))) = (
            term.rho.checked_sub(1),
            term.tau.checked_add(1),
            term.word.split_first(),
        ) else {
            return Ok(None);
        };
        if factor.as_ref() != [1] {
            return Ok(None);
        }
        candidate.push(Term {
            rho,
            tau,
            word: tail.into(),
        });
    }
    // d(tau) = rho[tau_0]. Verify the differential of the proposed
    // left product before returning it as a primitive.
    let image = real_motivic::differential_support(model, &candidate)?;
    let exact = real_motivic::compare(model, &image, terms)? == std::cmp::Ordering::Equal;
    Ok(exact.then_some(candidate))
}

fn classical_boundary(model: &mut Model, grade: Grade, terms: &[Term]) -> Result<bool> {
    let words = classical_words(model, terms)?;
    if words.is_empty() {
        return Ok(true);
    }
    let Some(previous) = grade.s.checked_sub(1) else {
        return Ok(false);
    };
    let k = rho_localized::classical_degree(grade)?;
    let source = cobar::basis(model, previous, k, MAX_BASIS)?;
    Ok(preimage(
        model,
        &source,
        &words,
        cobar::word_cells,
        cobar::differential_word,
        false,
    )?
    .is_some())
}

fn shifted_preimage(
    model: &mut Model,
    grade: Grade,
    terms: &[Term],
    power: usize,
    witness: bool,
) -> Result<Option<Vec<Term>>> {
    if power == 0 {
        return real_preimage(model, grade, terms, witness);
    }
    let shifted = real_motivic::shifted(grade, power)?;
    let shifted_terms = real_motivic::shift_support(model, terms, power)?;
    real_preimage(model, shifted, &shifted_terms, witness)
}

fn verify_cycle(model: &mut Model, terms: &[Term], known_cycle: bool) -> Result<()> {
    if !known_cycle && !real_motivic::differential_support(model, terms)?.is_empty() {
        return Err(not_closed());
    }
    Ok(())
}

pub(crate) fn order(
    model: &mut Model,
    grade: Grade,
    terms: &[Term],
    known_cycle: bool,
) -> Result<Option<usize>> {
    if terms.is_empty() {
        return Ok(Some(0));
    }
    if grade.s == 0 {
        verify_cycle(model, terms, known_cycle)?;
        return Ok(None);
    }

    // For t<=0 the real cobar already equals the localized slice.
    if grade.t <= 0 && (known_cycle || grade.s > 1) {
        if known_cycle {
            return Ok(None);
        }
        verify_cycle(model, terms, false)?;
        return Ok(classical_boundary(model, grade, terms)?.then_some(0));
    }

    // In length one there is no possible scalar primitive before t<=0.
    if grade.s == 1 {
        let bound = usize::try_from(grade.t.max(0)).map_err(|_| exponent_overflow())?;
        if shifted_preimage(model, grade, terms, bound, false)?.is_some() {
            return Ok(Some(bound));
        }
        verify_cycle(model, terms, known_cycle)?;
        return Ok(None);
    }

    if !classical_boundary(model, grade, terms)? {
        verify_cycle(model, terms, known_cycle)?;
        return Ok(None);
    }
    let mut hi = usize::try_from(grade.t).map_err(|_| exponent_overflow())?;
    let mut lo = usize::from(known_cycle);
    let mut known_boundary = false;
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if shifted_preimage(model, grade, terms, mid, false)?.is_some() {
            hi = mid;
            known_boundary = true;
        } else {
            lo = mid + 1;
        }
    }
    if !known_boundary && shifted_preimage(model, grade, terms, hi, false)?.is_none() {
        verify_cycle(model, terms, known_cycle)?;
        return Err(Error::plain(
            "real-motivic cobar boundary disagrees with its stable classical image",
        ));
    }
    // An actual boundary rho^hi z certifies that z was closed, since
    // multiplication by rho is injective on the normalized cochains.
    Ok(Some(hi))
}

fn exponent_overflow() -> Error {
    Error::plain("rho exponent is outside the supported range")
}

pub(crate) fn null(
    model: &mut Model,
    grade: Grade,
    terms: &[Term],
    power: usize,
    known_cycle: bool,
) -> Result<Value> {
    let Some(witness) = shifted_preimage(model, grade, terms, power, true)? else {
        verify_cycle(model, terms, known_cycle)?;
        return Err(Error::plain("rhoNull power does not kill the class"));
    };
    if power != 0 && shifted_preimage(model, grade, terms, power - 1, false)?.is_some() {
        return Err(Error::plain("rhoNull power is not the minimal exponent"));
    }
    real_motivic::persist(model, witness)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn degree_one_shortcut_preserves_ordinary_and_truncated_canonical_representatives() {
        check_degree_one(degree_one_cycles).unwrap_or_else(|error| panic!("{error}"));
    }

    #[test]
    fn filtered_degree_one_preserves_ordinary_and_truncated_canonical_representatives() {
        check_degree_one(filtered_degree_one_cycles).unwrap_or_else(|error| panic!("{error}"));
    }

    fn check_degree_one(
        build: impl Fn(&mut Model, Grade, Option<usize>) -> Result<Option<Vec<Vec<Term>>>>,
    ) -> Result<()> {
        let mut compared = 0;
        for cutoff in [None, Some(1), Some(2), Some(3), Some(5)] {
            for t in 1..=8 {
                for w in -4..=4 {
                    let grade = Grade { s: 1, t, w };
                    let mut model = Model::new();
                    if let Some(cycles) = build(&mut model, grade, cutoff)? {
                        let compact = Group::from_cycle_supports(&mut model, cycles)?;
                        let matrix = Group::from_matrix(&mut model, grade, cutoff, None)?;
                        assert_eq!(compact.rank(), matrix.rank(), "{grade:?}, {cutoff:?}");
                        for i in 0..matrix.rank() {
                            let ordinary = matrix.representative(&mut model, &[i])?;
                            let fast = compact.representative(&mut model, &[i])?;
                            assert_eq!(fast, ordinary, "{grade:?}, {cutoff:?}");
                            assert_eq!(compact.project(&mut model, &fast)?, [i]);
                            compared += 1;
                        }
                        if matrix.rank() > 1 {
                            let coords: Vec<_> = (0..matrix.rank()).collect();
                            let sum = matrix.representative(&mut model, &coords)?;
                            assert_eq!(compact.representative(&mut model, &coords)?, sum);
                            assert_eq!(compact.project(&mut model, &sum)?, coords);
                        }
                    }
                }
            }
        }
        assert!(compared >= 20);
        Ok(())
    }
}
