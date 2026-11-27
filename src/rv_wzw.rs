use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_cyclotomic::{self, Polynomial};
use crate::rv_fusion::{self, Grade as Fusion};

const MAX_CELLS: usize = 300_000;
const MAX_LABELS: usize = 4096;

type Weight = Rc<[i128]>;
pub(crate) type Support = Rc<[(Weight, Polynomial)]>;
type Cyclic = Vec<(usize, i128)>;

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) fusion: Fusion,
    pub(crate) order: usize,
}

impl Grade {
    pub(crate) fn checked(rank: usize, level: usize) -> Result<Self> {
        let fusion = Fusion::checked(rank, level)?;
        // lcm(2n,4) * L equals lcm(2nL,4L).
        let multiple = if rank.is_multiple_of(2) {
            2 * rank
        } else {
            4 * rank
        };
        let order = rank
            .checked_add(level)
            .and_then(|l| l.checked_mul(multiple))
            .ok_or_else(|| Error::plain("WZW cyclotomic order is out of range"))?;
        rv_cyclotomic::check_order(order)?;
        Ok(Self { fusion, order })
    }

    pub(crate) fn period(self) -> usize {
        self.fusion.rank + self.fusion.level
    }
}

fn arithmetic() -> Error {
    Error::plain("WZW cyclotomic arithmetic is out of range")
}

fn storage(cells: usize) -> Result<()> {
    if cells > MAX_CELLS {
        return Err(Error::plain(
            "WZW modular calculation exceeds the storage limit",
        ));
    }
    Ok(())
}

fn cells(weight: &[i128], coefficient: &[i128]) -> usize {
    4 + 2 * weight.len() + 2 * coefficient.len()
}

fn persist(model: &mut Model, terms: Vec<(Weight, Polynomial)>) -> Result<Value> {
    let mut total = 0;
    for (weight, coefficient) in &terms {
        total += cells(weight, coefficient);
        storage(total)?;
    }
    model.allocate(total)?;
    Ok(Value::RVWzw(terms.into()))
}

pub(crate) fn zero() -> Value {
    Value::RVWzw([].into())
}

pub(crate) fn basis(model: &mut Model, grade: Grade, weight: Vec<i128>) -> Result<Value> {
    let weight = grade.fusion.normalize(model, weight)?;
    persist(model, vec![(weight, [1].into())])
}

pub(crate) fn compare(model: &mut Model, a: &Support, b: &Support) -> Result<Ordering> {
    if Rc::ptr_eq(a, b) {
        return Ok(Ordering::Equal);
    }
    for (x, y) in a.iter().zip(b.iter()) {
        model.spend(1 + x.0.len() + x.1.len().min(y.1.len()))?;
        let order = x.cmp(y);
        if !order.is_eq() {
            return Ok(order);
        }
    }
    Ok(a.len().cmp(&b.len()))
}

pub(crate) fn add(model: &mut Model, a: &Support, b: &Support) -> Result<Value> {
    let mut merged = BTreeMap::new();
    let mut kept = 0;
    for (weight, coefficient) in a.iter().chain(b.iter()) {
        model.spend(1 + weight.len())?;
        let accumulated = merged.entry(weight.clone()).or_insert_with(Vec::new);
        let old = accumulated.len();
        rv_cyclotomic::accumulate(model, accumulated, coefficient, false)?;
        kept += 2 * accumulated.len();
        kept -= 2 * old;
        storage(kept + merged.len() * (4 + 2 * weight.len()))?;
    }
    persist(
        model,
        merged
            .into_iter()
            .filter(|(_, c)| !c.is_empty())
            .map(|(w, c)| (w, c.into()))
            .collect(),
    )
}

pub(crate) fn negate(model: &mut Model, input: &Support) -> Result<Value> {
    let terms = input
        .iter()
        .map(|(weight, coefficient)| {
            model.spend(coefficient.len())?;
            let c = coefficient
                .iter()
                .map(|x| x.checked_neg().ok_or_else(arithmetic))
                .collect::<Result<Polynomial>>()?;
            Ok((weight.clone(), c))
        })
        .collect::<Result<Vec<_>>>()?;
    persist(model, terms)
}

pub(crate) fn scale(
    model: &mut Model,
    grade: Grade,
    input: &Support,
    scalar: &Polynomial,
) -> Result<Value> {
    if scalar.is_empty() {
        return Ok(zero());
    }
    let mut result = Vec::new();
    let mut kept = 0;
    for (weight, coefficient) in input.iter() {
        let value = rv_cyclotomic::product_coefficients(model, grade.order, scalar, coefficient)?;
        if !value.is_empty() {
            kept += cells(weight, &value);
            storage(kept)?;
            result.push((weight.clone(), value.into()));
        }
    }
    persist(model, result)
}

pub(crate) fn charge(model: &mut Model, input: &Support) -> Result<Value> {
    let mut output = Vec::with_capacity(input.len());
    for (weight, coefficient) in input.iter() {
        model.spend(weight.len() * (1 + input.len().max(1).ilog2() as usize))?;
        output.push((rv_fusion::dual_weight(weight), coefficient.clone()));
    }
    output.sort_unstable();
    persist(model, output)
}

pub(crate) fn coordinate(
    model: &mut Model,
    grade: Grade,
    input: &Support,
    parts: Vec<i128>,
) -> Result<Value> {
    let weight = grade.fusion.normalize(model, parts)?;
    model.spend(weight.len() * (1 + input.len().max(1).ilog2() as usize))?;
    Ok(Value::RVCyclo(
        input
            .binary_search_by(|(w, _)| w.cmp(&weight))
            .map(|i| input[i].1.clone())
            .unwrap_or_else(|_| [].into()),
    ))
}

fn shifted_weight(weight: &[i128]) -> Vec<i128> {
    weight
        .iter()
        .rev()
        .enumerate()
        .map(|(i, &entry)| entry + (weight.len() - 1 - i) as i128)
        .collect()
}

fn shifted_norm(weight: &[i128]) -> i128 {
    let sum: i128 = weight.iter().sum();
    weight.len() as i128 * weight.iter().map(|&x| x * x).sum::<i128>() - sum * sum
}

fn twist(grade: Grade, weight: &[i128], power: i128) -> i128 {
    let x = shifted_weight(weight);
    let rho: Vec<_> = (0..weight.len()).rev().map(|i| i as i128).collect();
    let scale = (grade.order / (2 * grade.fusion.rank * grade.period())) as i128;
    let order = grade.order as i128;
    // Reduce the arbitrary signed power first, including i128::MIN.
    power.rem_euclid(order) * (scale * (shifted_norm(&x) - shifted_norm(&rho))).rem_euclid(order)
}

pub(crate) fn t_entry(
    model: &mut Model,
    grade: Grade,
    weight: Vec<i128>,
    power: i128,
) -> Result<Value> {
    let weight = grade.fusion.normalize(model, weight)?;
    rv_cyclotomic::root(model, grade.order, twist(grade, &weight, power))
}

pub(crate) fn t_support(
    model: &mut Model,
    grade: Grade,
    input: &Support,
    power: i128,
) -> Result<Vec<(Weight, Polynomial)>> {
    let mut state = Vec::new();
    let mut kept = 0;
    for (weight, coefficient) in input.iter() {
        model.spend(weight.len())?;
        let output =
            rv_cyclotomic::shifted(model, grade.order, coefficient, twist(grade, weight, power))?;
        kept += cells(weight, &output);
        storage(kept)?;
        state.push((weight.clone(), output.into()));
    }
    Ok(state)
}

pub(crate) fn t_action(
    model: &mut Model,
    grade: Grade,
    input: &Support,
    power: i128,
) -> Result<Value> {
    let output = t_support(model, grade, input, power)?;
    persist(model, output)
}

fn cyclic(raw: Vec<i128>) -> Cyclic {
    raw.into_iter()
        .enumerate()
        .filter(|(_, c)| *c != 0)
        .collect()
}

fn shift_add(
    model: &mut Model,
    output: &mut [i128],
    input: &Cyclic,
    shift: usize,
    factor: i128,
) -> Result<()> {
    model.spend(input.len())?;
    let period = output.len();
    for &(power, scalar) in input {
        let entry = &mut output[(power + shift) % period];
        *entry = scalar
            .checked_mul(factor)
            .and_then(|term| entry.checked_add(term))
            .ok_or_else(arithmetic)?;
    }
    Ok(())
}

fn vacuum_alternant(model: &mut Model, period: usize, x: &[i128], sign: i128) -> Result<Cyclic> {
    let mut result = vec![(0, 1)];
    for (i, &a) in x.iter().enumerate() {
        for &b in &x[i + 1..] {
            model.spend(period)?;
            let mut product = vec![0; period];
            shift_add(
                model,
                &mut product,
                &result,
                (sign * a).rem_euclid(period as i128) as usize,
                1,
            )?;
            shift_add(
                model,
                &mut product,
                &result,
                (sign * b).rem_euclid(period as i128) as usize,
                -1,
            )?;
            result = cyclic(product);
        }
    }
    Ok(result)
}

fn alternant(
    model: &mut Model,
    period: usize,
    x: &[i128],
    y: &[i128],
    sign: i128,
) -> Result<Cyclic> {
    let rank = x.len();
    storage(2 * period * rank * rank)?;
    let shifts: Vec<usize> = x
        .iter()
        .flat_map(|&a| {
            y.iter()
                .map(move |&b| (sign * a * b).rem_euclid(period as i128) as usize)
        })
        .collect();
    let mut matrix: Vec<Cyclic> = shifts.iter().map(|&power| vec![(power, 1)]).collect();
    let mut traces = Vec::with_capacity(rank);
    // Compute power traces in the torsion-free integral cyclic group ring.
    // Multiplication by the original matrix only shifts powers of its root.
    for degree in 1..=rank {
        model.spend(period)?;
        let mut trace = vec![0; period];
        for i in 0..rank {
            shift_add(model, &mut trace, &matrix[i * rank + i], 0, 1)?;
        }
        traces.push(cyclic(trace));
        if degree < rank {
            let mut next = Vec::with_capacity(rank * rank);
            for row in matrix.chunks_exact(rank) {
                for column in 0..rank {
                    model.spend(period)?;
                    let mut product = vec![0; period];
                    for (inner, entry) in row.iter().enumerate() {
                        shift_add(model, &mut product, entry, shifts[inner * rank + column], 1)?;
                    }
                    next.push(cyclic(product));
                }
            }
            matrix = next;
        }
    }
    let mut elementary = vec![vec![(0, 1)]];
    for degree in 1..=rank {
        model.spend(period)?;
        let mut result = vec![0; period];
        for (i, trace) in traces.iter().enumerate().take(degree) {
            let factor = if i % 2 == 0 { 1 } else { -1 };
            for &(power, scalar) in trace {
                shift_add(
                    model,
                    &mut result,
                    &elementary[degree - i - 1],
                    power,
                    scalar.checked_mul(factor).ok_or_else(arithmetic)?,
                )?;
            }
        }
        let divisor = degree as i128;
        for scalar in &mut result {
            if *scalar % divisor != 0 {
                return Err(Error::plain(
                    "WZW Newton determinant division was not exact",
                ));
            }
            *scalar /= divisor;
        }
        elementary.push(cyclic(result));
    }
    Ok(elementary.pop().unwrap())
}

fn s_coefficient(
    model: &mut Model,
    grade: Grade,
    left: &[i128],
    right: &[i128],
    adjoint: bool,
) -> Result<Vec<i128>> {
    let period = grade.period();
    let rank = grade.fusion.rank;
    let x = shifted_weight(left);
    let y = shifted_weight(right);
    let sign = if adjoint { 1 } else { -1 };
    let determinant = if left[rank - 1] == 0 {
        vacuum_alternant(model, period, &y, sign)?
    } else if right[rank - 1] == 0 {
        vacuum_alternant(model, period, &x, sign)?
    } else {
        alternant(model, period, &x, &y, sign)?
    };
    let factor = (grade.order / (rank * period)) as i128;
    let phase = -sign * factor * x.iter().sum::<i128>() * y.iter().sum::<i128>();
    let step = (grade.order / period) as i128;
    model.spend(grade.order)?;
    let mut output = vec![0; grade.order];
    for (i, coefficient) in determinant {
        let exponent = (phase + step * i as i128).rem_euclid(grade.order as i128) as usize;
        output[exponent] = coefficient;
    }
    rv_cyclotomic::reduced(model, grade.order, output)
}

pub(crate) fn s_entry(
    model: &mut Model,
    grade: Grade,
    left: Vec<i128>,
    right: Vec<i128>,
) -> Result<Value> {
    let left = grade.fusion.normalize(model, left)?;
    let right = grade.fusion.normalize(model, right)?;
    let coefficient = s_coefficient(model, grade, &left, &right, false)?;
    rv_cyclotomic::persist(model, coefficient)
}

pub(crate) fn s_support(
    model: &mut Model,
    grade: Grade,
    input: &Support,
    adjoint: bool,
) -> Result<Vec<(Weight, Polynomial)>> {
    if input.is_empty() {
        return Ok(Vec::new());
    }
    let mut output = Vec::new();
    let mut kept = 0;
    // The raw S alternant is symmetric, so conjugating its entries gives
    // the conjugate transpose. This does not divide by the normalization.
    rv_fusion::visit_weights(
        model,
        grade.fusion,
        MAX_LABELS,
        "WZW S action exceeds the limit of 4096 integrable weights",
        |model, weight| {
            let mut sum = Vec::new();
            for (source, coefficient) in input.iter() {
                let entry = s_coefficient(model, grade, weight, source, adjoint)?;
                let term =
                    rv_cyclotomic::product_coefficients(model, grade.order, &entry, coefficient)?;
                rv_cyclotomic::accumulate(model, &mut sum, &term, false)?;
            }
            if !sum.is_empty() {
                kept += cells(weight, &sum);
                storage(kept)?;
                output.push((weight.into(), sum.into()));
            }
            Ok(())
        },
    )?;
    Ok(output)
}
pub(crate) fn s_action(
    model: &mut Model,
    grade: Grade,
    input: &Support,
    adjoint: bool,
) -> Result<Value> {
    let output = s_support(model, grade, input, adjoint)?;
    persist(model, output)
}
