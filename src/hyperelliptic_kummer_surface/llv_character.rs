use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::{Model, MAX_SHOWN_BYTES};
use crate::rv_qp_quad::Wide;
use crate::wide_integer::Natural;

mod hodge;
mod reduced;
mod weyl;

const RANK: usize = 12;
const LARGEST_RHO: usize = RANK - 1;
const LIMBS: usize = 20;
const WEYL_LIMBS: usize = 32;
const SUBSETS: usize = 1 << RANK;
type Integer = Wide<LIMBS>;
type Magnitude = Natural<WEYL_LIMBS>;
type Signed = Wide<WEYL_LIMBS>;

pub(crate) enum Command {
    Multiplicity(u128, Vec<u128>),
    Dimension(Vec<u128>),
    Decompose(u128),
    Hodge(u128, i128, i128, Vec<u128>),
    HodgeDecompose(u128, i128, i128),
    HodgeTable(u128, Option<Vec<u128>>),
    ReducedMultiplicity(u128, i128, Vec<u128>),
    ReducedDecompose(u128, i128),
    PrimitiveMultiplicity(u128, i128, Vec<i128>),
    PrimitiveDecompose(u128, i128),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<String> {
        match self {
            Self::Multiplicity(n, shape) => {
                let shape = validate(shape)?;
                validate_length(n)?;
                let multiplicity = if weight_exceeds(&shape, n) {
                    0
                } else {
                    Oscillators::new(model, n)?.multiplicity(model, &shape)?
                };
                Ok(format!(
                    "rvHilbLLVMultiplicity[{n}]{shape:?} = {multiplicity}"
                ))
            }
            Self::Dimension(shape) => {
                let shape = validate(shape)?;
                let dimension = weyl::dimension(model, &shape)?;
                let dimension = decimal(model, dimension)?;
                Ok(format!("rvHilbLLVDimension{shape:?} = {dimension}"))
            }
            Self::Decompose(n) => {
                validate_length(n)?;
                decomposition(model, n)
            }
            Self::Hodge(n, p, q, shape) => {
                let shape = validate(shape)?;
                validate_length(n)?;
                hodge::coefficient(model, n, p, q, &shape)
            }
            Self::HodgeDecompose(n, p, q) => {
                validate_length(n)?;
                hodge::decompose(model, n, p, q)
            }
            Self::HodgeTable(n, shape) => {
                let shape = shape.map(validate).transpose()?;
                validate_length(n)?;
                hodge::table(model, n, shape.as_deref())
            }
            Self::ReducedMultiplicity(n, d, shape) => reduced::multiplicity(model, n, d, shape),
            Self::ReducedDecompose(n, d) => reduced::decomposition(model, n, d, false),
            Self::PrimitiveMultiplicity(n, d, shape) => {
                reduced::primitive_multiplicity(model, n, d, shape)
            }
            Self::PrimitiveDecompose(n, d) => reduced::decomposition(model, n, d, true),
        }
    }
}

fn validate(mut shape: Vec<u128>) -> Result<Vec<u128>> {
    if !shape.windows(2).all(|pair| pair[0] >= pair[1]) {
        return Err(Error::plain(
            "type B12 highest weight must be weakly decreasing",
        ));
    }
    while shape.last() == Some(&0) {
        shape.pop();
    }
    if shape.len() > RANK {
        return Err(Error::plain(
            "type B12 highest weight has at most twelve nonzero coordinates",
        ));
    }
    Ok(shape)
}

fn validate_length(n: u128) -> Result<()> {
    if n == 1 {
        return Err(Error::plain(
            "type B12 LLV reports exclude length one, whose LLV type is D12",
        ));
    }
    Ok(())
}

fn weight_exceeds(shape: &[u128], n: u128) -> bool {
    shape
        .iter()
        .try_fold(n, |budget, &part| budget.checked_sub(part))
        .is_none()
}

struct Oscillators {
    n: usize,
    coefficients: Vec<i128>,
}

impl Oscillators {
    fn new(model: &mut Model, n: u128) -> Result<Self> {
        let n = usize::try_from(n).map_err(|_| work_limit())?;
        let width = n.saturating_mul(2).saturating_add(1);
        let cells = width.saturating_mul(n.saturating_add(1));
        guard(cells.saturating_add(4 * SUBSETS))?;
        model
            .spend(cells.saturating_mul(n))
            .map_err(|_| work_limit())?;
        let mut coefficients = vec![0i128; cells];
        coefficients[n] = 1;
        for part in 1..=n {
            for positive in [true, false] {
                for degree in part..=n {
                    let source = (degree - part) * width;
                    let target = degree * width;
                    let charges = if positive { 1..width } else { 0..width - 1 };
                    for charge in charges {
                        let predecessor = if positive { charge - 1 } else { charge + 1 };
                        coefficients[target + charge] = coefficients[target + charge]
                            .checked_add(coefficients[source + predecessor])
                            .ok_or_else(arithmetic)?;
                    }
                }
            }
        }
        Ok(Self { n, coefficients })
    }

    fn get(&self, charge: usize, degree: usize) -> i128 {
        if charge > degree || degree > self.n {
            return 0;
        }
        self.coefficients[degree * (2 * self.n + 1) + self.n + charge]
    }

    fn multiplicity(&self, model: &mut Model, shape: &[u128]) -> Result<i128> {
        self.extracted(model, shape, RANK, None)
    }

    fn extracted(
        &self,
        model: &mut Model,
        shape: &[u128],
        rank: usize,
        charge: Option<usize>,
    ) -> Result<i128> {
        let offset = charge.unwrap_or(0);
        if offset > self.n || weight_exceeds(shape, (self.n - offset) as u128) {
            return Ok(0);
        }
        let mut padded = [0i128; RANK];
        for (out, &part) in padded.iter_mut().zip(shape) {
            *out = part as i128;
        }
        let parts = &padded[..rank];
        let excess = self.n - offset - parts.iter().sum::<i128>() as usize;
        let entries = self.matrix(model, parts, excess)?;
        let tail = min_tail(model, &entries, rank, excess)?;
        if tail[0] > excess {
            return Ok(0);
        }
        // The grading variable was not part of the reduced Weyl denominator.
        // Its coefficient contributes a single scalar oscillator series.
        let initial = if let Some(charge) = charge {
            (0..=excess - tail[0])
                .map(|i| Integer::from_i128(self.get(charge, charge + i)))
                .collect()
        } else {
            vec![Integer::from_i128(1)]
        };
        let answer = determinant(model, &entries, &tail, rank, excess, initial)?;
        if answer < 0 {
            return Err(Error::plain(
                "LLV Weyl extraction produced a negative multiplicity",
            ));
        }
        Ok(answer)
    }

    fn visit_summands(
        &self,
        model: &mut Model,
        mut visit: impl FnMut(&mut Model, &[u128], i128) -> Result<()>,
    ) -> Result<()> {
        let mut shape = Vec::new();
        let mut emit = |model: &mut Model, shape: &[u128]| {
            let multiplicity = self.multiplicity(model, shape)?;
            if multiplicity == 0 {
                Ok(())
            } else {
                visit(model, shape, multiplicity)
            }
        };
        for weight in (0..=self.n).rev() {
            partitions(model, weight, weight, &mut shape, &mut emit)?;
        }
        Ok(())
    }

    fn matrix(&self, model: &mut Model, parts: &[i128], excess: usize) -> Result<Vec<Entry>> {
        let rank = parts.len();
        model.spend(rank * rank * (excess + 1))?;
        let mut entries = Vec::with_capacity(rank * rank);
        for (row, &lambda) in parts.iter().enumerate() {
            for column in 0..rank {
                let shift = lambda + column as i128 - row as i128;
                let charge = shift.unsigned_abs() as usize;
                let other = (lambda + (2 * rank - 1 - row - column) as i128) as usize;
                // Dividing each entry by Q^shift divides the determinant by
                // Q^|lambda|, since the sum of the column shifts is fixed.
                let first = (charge as i128 - shift) as usize;
                let mut coefficients = Vec::new();
                if shift <= self.n as i128 {
                    let last = ((self.n as i128 - shift) as usize).min(excess);
                    for degree in first..=last {
                        let original = (degree as i128 + shift) as usize;
                        coefficients.push(self.get(charge, original) - self.get(other, original));
                    }
                }
                entries.push(Entry {
                    first,
                    coefficients,
                });
            }
        }
        Ok(entries)
    }
}

struct Entry {
    first: usize,
    coefficients: Vec<i128>,
}

fn min_tail(
    model: &mut Model,
    entries: &[Entry],
    rank: usize,
    excess: usize,
) -> Result<Vec<usize>> {
    let subsets = 1 << rank;
    model.spend(rank * subsets / 2)?;
    let mut tail = vec![excess + 1; subsets];
    let end = subsets - 1;
    tail[end] = 0;
    for mask in (0..end).rev() {
        let row = mask.count_ones() as usize;
        for column in 0..rank {
            let bit = 1 << column;
            let entry = &entries[row * rank + column];
            if mask & bit == 0 && !entry.coefficients.is_empty() {
                tail[mask] = tail[mask].min(entry.first + tail[mask | bit]);
            }
        }
    }
    Ok(tail)
}

fn determinant(
    model: &mut Model,
    entries: &[Entry],
    tail: &[usize],
    rank: usize,
    excess: usize,
    initial: Vec<Integer>,
) -> Result<i128> {
    let subsets = tail.len();
    let zero = Integer::from_i128(0);
    let mut partials: Vec<Option<Vec<Integer>>> = (0..subsets).map(|_| None).collect();
    let mut cells = 4 * subsets + LIMBS * initial.len();
    guard(cells)?;
    partials[0] = Some(initial);
    let end = subsets - 1;
    for mask in 0..end {
        let Some(prefix) = partials[mask].take() else {
            continue;
        };
        let row = mask.count_ones() as usize;
        model.spend(rank + prefix.len())?;
        if prefix.iter().all(|value| value.is_zero()) {
            cells -= LIMBS * prefix.len();
            continue;
        }
        for column in 0..rank {
            let bit = 1 << column;
            let entry = &entries[row * rank + column];
            let next = mask | bit;
            if mask & bit != 0 || entry.coefficients.is_empty() || entry.first + tail[next] > excess
            {
                continue;
            }
            let width = excess - tail[next] + 1;
            if partials[next].is_none() {
                cells += LIMBS * width;
                guard(cells)?;
                partials[next] = Some(vec![zero; width]);
            }
            let output = partials[next].as_mut().unwrap();
            let negative = !(mask >> (column + 1)).count_ones().is_multiple_of(2);
            add_product(model, output, &prefix, entry, negative)?;
        }
        cells -= LIMBS * prefix.len();
    }
    let answer = partials[end]
        .as_ref()
        .map_or(zero, |coefficients| coefficients[excess]);
    answer.signed().ok_or_else(arithmetic)
}

fn add_product(
    model: &mut Model,
    output: &mut [Integer],
    prefix: &[Integer],
    entry: &Entry,
    negative: bool,
) -> Result<()> {
    for (degree, left) in prefix
        .iter()
        .copied()
        .take(output.len() - entry.first)
        .enumerate()
    {
        if left.is_zero() {
            continue;
        }
        let start = degree + entry.first;
        let width = (output.len() - start).min(entry.coefficients.len());
        model.spend(width)?;
        for (out, &right) in output[start..start + width]
            .iter_mut()
            .zip(&entry.coefficients)
        {
            if right == 0 {
                continue;
            }
            model.spend(LIMBS)?;
            let product = left.multiply(right).map_err(|_| arithmetic())?;
            *out = out
                .add(if negative { product.negate() } else { product })
                .map_err(|_| arithmetic())?;
        }
    }
    Ok(())
}

fn decomposition(model: &mut Model, length: u128) -> Result<String> {
    use std::fmt::Write;
    let oscillators = Oscillators::new(model, length)?;
    let n = oscillators.n;
    let denominator = weyl::specialized_numerator(model, &[])?;
    let mut total = Magnitude::default();
    let mut grades = vec![Magnitude::default(); 2 * n + 1];
    let mut report = format!("rvHilbLLVDecomposition[{length}] = components=[");
    let mut first = true;
    oscillators.visit_summands(model, |model, shape, multiplicity| {
        let dimension = weyl::dimension(model, shape)?;
        let weights = weyl::specialize(model, shape, &denominator)?;
        let mut check = weights[0];
        for (offset, weight) in weights.into_iter().enumerate() {
            if offset > 0 {
                check.add_multiple(weight, 2).ok_or_else(arithmetic)?;
            }
            let contribution = weight
                .multiply_u128(multiplicity as u128)
                .ok_or_else(arithmetic)?;
            grades[n - offset]
                .add_multiple(contribution, 1)
                .ok_or_else(arithmetic)?;
            if offset > 0 {
                grades[n + offset]
                    .add_multiple(contribution, 1)
                    .ok_or_else(arithmetic)?;
            }
        }
        if check != dimension {
            return Err(Error::plain(
                "LLV specialized weights do not sum to the Weyl dimension",
            ));
        }
        total
            .add_multiple(
                dimension
                    .multiply_u128(multiplicity as u128)
                    .ok_or_else(arithmetic)?,
                1,
            )
            .ok_or_else(arithmetic)?;
        if !first {
            report.push_str(", ");
        }
        first = false;
        let dimension = decimal(model, dimension)?;
        write!(report, "({shape:?}, {multiplicity}, {dimension})").unwrap();
        if report.len() > MAX_SHOWN_BYTES {
            return Err(report_limit());
        }
        Ok(())
    })?;
    let total = decimal(model, total)?;
    write!(report, "] total={total} betti=[").unwrap();
    for (index, betti) in grades.into_iter().enumerate() {
        if index > 0 {
            report.push_str(", ");
        }
        report.push_str(&decimal(model, betti)?);
        if report.len() > MAX_SHOWN_BYTES {
            return Err(report_limit());
        }
    }
    report.push(']');
    Ok(report)
}

fn partitions(
    model: &mut Model,
    weight: usize,
    max_part: usize,
    shape: &mut Vec<u128>,
    emit: &mut impl FnMut(&mut Model, &[u128]) -> Result<()>,
) -> Result<()> {
    if weight == 0 {
        return emit(model, shape);
    }
    let available = RANK - shape.len();
    if available == 0 {
        return Ok(());
    }
    for part in (weight.div_ceil(available)..=max_part.min(weight)).rev() {
        model.spend(1)?;
        shape.push(part as u128);
        partitions(model, weight - part, part, shape, emit)?;
        shape.pop();
    }
    Ok(())
}

fn decimal(model: &mut Model, mut value: Magnitude) -> Result<String> {
    use std::fmt::Write;
    if let Some(small) = value.small() {
        return Ok(small.to_string());
    }
    let mut digits = Vec::new();
    while !value.is_zero() {
        model.spend(WEYL_LIMBS)?;
        let (quotient, remainder) = value.divide_u32(1_000_000_000);
        digits.push(remainder);
        value = quotient;
    }
    let mut text = digits.pop().unwrap().to_string();
    for digit in digits.into_iter().rev() {
        write!(text, "{digit:09}").unwrap();
    }
    Ok(text)
}

fn guard(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        Err(work_limit())
    } else {
        Ok(())
    }
}

fn work_limit() -> Error {
    Error::plain("Hilbert LLV character exceeds the polynomial computation limit")
}

fn report_limit() -> Error {
    Error::plain("Hilbert LLV decomposition exceeds the report limit")
}

fn arithmetic() -> Error {
    Error::plain("Hilbert LLV character exceeds the exact integer range")
}
