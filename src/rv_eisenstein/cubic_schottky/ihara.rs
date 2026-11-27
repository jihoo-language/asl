use crate::error::{Error, Result};
use crate::model::Model;

use super::{least_rotation, period, Letter, ReportKind};

const MAX_EDGES: usize = 256;
const MAX_DEGREE: usize = 256;

#[derive(Clone, Copy)]
struct Arrow {
    tail: usize,
    head: usize,
}

pub(super) struct Core {
    lengths: Vec<usize>,
    chains: Vec<Vec<usize>>,
    arrows: Vec<Arrow>,
    outgoing: Vec<Vec<usize>>,
    next: Vec<Vec<usize>>,
}

impl Core {
    pub(super) fn new(model: &mut Model, lengths: Vec<usize>) -> Result<Self> {
        let size = lengths
            .iter()
            .try_fold(0usize, |a, &b| a.checked_add(b))
            .ok_or_else(core_limit)?;
        if size > MAX_EDGES {
            return Err(core_limit());
        }
        model.spend(2 * size + 1)?;
        let mut chains = Vec::with_capacity(lengths.len());
        let mut arrows = Vec::with_capacity(2 * size);
        let mut vertices = 1;
        for &length in &lengths {
            let chain: Vec<_> = (vertices..vertices + length - 1).collect();
            vertices += length - 1;
            let trail = std::iter::once(0)
                .chain(chain.iter().copied())
                .chain(std::iter::once(0));
            for (tail, head) in trail.clone().zip(trail.skip(1)) {
                arrows.push(Arrow { tail, head });
                arrows.push(Arrow {
                    tail: head,
                    head: tail,
                });
            }
            chains.push(chain);
        }
        let mut outgoing = vec![Vec::new(); vertices];
        for (i, arrow) in arrows.iter().enumerate() {
            outgoing[arrow.tail].push(i);
        }
        let mut next = Vec::with_capacity(arrows.len());
        for (i, arrow) in arrows.iter().enumerate() {
            let options = &outgoing[arrow.head];
            model.spend(options.len())?;
            next.push(options.iter().copied().filter(|&j| j != i ^ 1).collect());
        }
        Ok(Self {
            lengths,
            chains,
            arrows,
            outgoing,
            next,
        })
    }

    pub(super) fn report(&self, model: &mut Model, name: &str, kind: ReportKind) -> Result<String> {
        match kind {
            ReportKind::Marking | ReportKind::Weighted(..) => unreachable!(),
            ReportKind::Core => {
                let edges: Vec<_> = self
                    .arrows
                    .iter()
                    .step_by(2)
                    .map(|a| (a.tail, a.head))
                    .collect();
                let directed = self.arrows.len();
                Ok(format!("rvCubicSchottkyCore[{name}] = vertices={}, edges={edges:?}, directed={directed}, rank={}",
                    self.outgoing.len(),self.lengths.len()))
            }
            ReportKind::Nonbacktracking(i, j) => {
                let row = self
                    .next
                    .get(i)
                    .filter(|_| j < self.arrows.len())
                    .ok_or_else(|| Error::plain("Schottky directed edge index is out of range"))?;
                Ok(format!(
                    "rvCubicSchottkyNB[{name},{i},{j}] = {}",
                    usize::from(row.contains(&j))
                ))
            }
            ReportKind::Ihara => {
                let polynomial = self.determinant(model)?;
                Ok(format!(
                    "rvCubicIhara[{name}] = reciprocal={polynomial:?}, bass_verified=true"
                ))
            }
            ReportKind::Coefficient(n) => {
                check_degree(n)?;
                let p = self.determinant(model)?;
                let c = zeta_coefficient(model, &p, n)?;
                Ok(format!("rvCubicIharaCoeff[{name},{n}] = {c}"))
            }
            ReportKind::Trace(n) | ReportKind::Primitive(n) | ReportKind::Enumerate(n) => {
                check_degree(n)?;
                let traces = self.traces(model, n)?;
                let p = self.determinant(model)?;
                if logarithmic_traces(model, &p, n)? != traces[1..] {
                    return Err(inconsistent());
                }
                if matches!(kind, ReportKind::Trace(_)) {
                    return Ok(format!("rvCubicSchottkyTrace[{name},{n}] = {}", traces[n]));
                }
                let primitive = primitive_count(&traces, n)?;
                if matches!(kind, ReportKind::Enumerate(_)) {
                    let direct = self.enumerate(model, n, &mut Vec::new())?;
                    if direct != primitive {
                        return Err(inconsistent());
                    }
                    return Ok(format!(
                        "rvCubicSchottkyEnumerate[{name},{n}] = {direct}, verified=true"
                    ));
                }
                Ok(format!(
                    "rvCubicSchottkyPrimitive[{name},{n}] = {primitive}"
                ))
            }
        }
    }

    fn determinant(&self, model: &mut Model) -> Result<Polynomial> {
        let p = self.weighted(model)?;
        if p != self.bass(model)? {
            return Err(inconsistent());
        }
        Ok(p)
    }

    fn weighted(&self, model: &mut Model) -> Result<Polynomial> {
        // I-B(u) is a rank-one update of paired blocks with diagonal 1 and
        // off-diagonal u^m. Taking the two-by-two adjugates clears denominators.
        let factors: Vec<_> = self.lengths.iter().map(|&m| one_minus(2 * m)).collect();
        let (prefix, suffix) = products(model, &factors)?;
        let mut result = prefix.last().expect("constant prefix").clone();
        for (i, &m) in self.lengths.iter().enumerate() {
            let cofactor = multiply(model, &prefix[i], &suffix[i + 1])?;
            let correction = multiply(model, &cofactor, &one_minus(m))?;
            add_shift(model, &mut result, &correction, m, -2)?;
        }
        Ok(result)
    }

    fn bass(&self, model: &mut Model) -> Result<Polynomial> {
        if self.lengths.is_empty() {
            // A single vertex gives (1-u^2)/(1-u^2).
            return Ok(vec![1]);
        }
        let loops = self.outgoing[0]
            .iter()
            .filter(|&&i| self.arrows[i].head == 0)
            .count();
        let root = vec![1, -(loops as i128), self.outgoing[0].len() as i128 - 1];

        let mut factors = Vec::new();
        let mut root_corrections = Vec::new();
        for chain in self.chains.iter().filter(|chain| !chain.is_empty()) {
            // Each interior vertex has degree two, so its diagonal entry is
            // 1+u^2. The reduced chain is tridiagonal with -u on the links.
            for (i, &v) in chain.iter().enumerate() {
                let mut actual: Vec<_> = self.outgoing[v]
                    .iter()
                    .map(|&e| self.arrows[e].head)
                    .collect();
                let mut expected = [
                    i.checked_sub(1).map_or(0, |j| chain[j]),
                    chain.get(i + 1).copied().unwrap_or(0),
                ];
                actual.sort_unstable();
                expected.sort_unstable();
                if actual != expected {
                    return Err(inconsistent());
                }
            }
            let length = chain.len();
            let mut prior = vec![1];
            let mut current = vec![1, 0, 1];
            for _ in 1..length {
                let mut next = current.clone();
                add_shift(model, &mut next, &current, 2, 1)?;
                add_shift(model, &mut next, &prior, 2, -1)?;
                prior = current;
                current = next;
            }
            factors.push(current);
            // The two diagonal endpoint cofactors are det(chain minus an
            // endpoint); the opposite endpoint cofactor is u^(length-1).
            add_shift(model, &mut prior, &[1], length - 1, 1)?;
            root_corrections.push(prior);
        }
        let (prefix, suffix) = products(model, &factors)?;
        let mut numerator = multiply(model, &root, prefix.last().expect("constant prefix"))?;
        for (i, correction) in root_corrections.iter().enumerate() {
            let others = multiply(model, &prefix[i], &suffix[i + 1])?;
            let term = multiply(model, correction, &others)?;
            add_shift(model, &mut numerator, &term, 2, -2)?;
        }
        for _ in 1..self.lengths.len() {
            let copy = numerator.clone();
            add_shift(model, &mut numerator, &copy, 2, -1)?;
        }
        Ok(numerator)
    }

    fn traces(&self, model: &mut Model, n: usize) -> Result<Vec<i128>> {
        let size = self.arrows.len();
        let mut traces = vec![0i128; n + 1];
        traces[0] = size as i128;
        for start in 0..size {
            let mut current = vec![0i128; size];
            current[start] = 1;
            for trace in traces.iter_mut().skip(1) {
                model.spend(size)?;
                let mut next = vec![0i128; size];
                for (i, &coefficient) in current.iter().enumerate().filter(|(_, c)| **c != 0) {
                    model.spend(self.next[i].len())?;
                    for &j in &self.next[i] {
                        next[j] = add(next[j], coefficient)?;
                    }
                }
                *trace = add(*trace, next[start])?;
                current = next;
            }
        }
        Ok(traces)
    }

    fn enumerate(&self, model: &mut Model, n: usize, word: &mut Vec<Letter>) -> Result<i128> {
        if n == 0 {
            return Ok(i128::from(
                !word.is_empty()
                    && word.last() != word.first().map(|&a| -a).as_ref()
                    && period(word) == word.len()
                    && least_rotation(word) == 0,
            ));
        }
        let mut total = 0;
        for (i, &length) in self.lengths.iter().enumerate() {
            if length > n {
                continue;
            }
            let a = (i as Letter) + 1;
            for next in [a, -a] {
                model.spend(1)?;
                if word.last() == Some(&-next) {
                    continue;
                }
                word.push(next);
                total = add(total, self.enumerate(model, n - length, word)?)?;
                word.pop();
            }
        }
        Ok(total)
    }
}

type Polynomial = Vec<i128>;

fn check_degree(n: usize) -> Result<()> {
    if n > MAX_DEGREE {
        Err(Error::plain(
            "Ihara coefficient and trace queries support degrees up to 256",
        ))
    } else {
        Ok(())
    }
}

fn one_minus(n: usize) -> Polynomial {
    let mut p = vec![0; n + 1];
    p[0] = 1;
    p[n] = -1;
    p
}

fn products(
    model: &mut Model,
    factors: &[Polynomial],
) -> Result<(Vec<Polynomial>, Vec<Polynomial>)> {
    let mut pre = vec![vec![1]];
    for factor in factors {
        pre.push(multiply(
            model,
            pre.last().expect("constant prefix"),
            factor,
        )?);
    }
    let mut post = vec![vec![1]];
    for factor in factors.iter().rev() {
        post.push(multiply(
            model,
            post.last().expect("constant suffix"),
            factor,
        )?);
    }
    post.reverse();
    Ok((pre, post))
}

fn normalize(polynomial: &mut Polynomial) {
    while polynomial.last() == Some(&0) {
        polynomial.pop();
    }
}

fn multiply(model: &mut Model, a: &[i128], b: &[i128]) -> Result<Polynomial> {
    let mut result = Vec::new();
    for (degree, &factor) in a.iter().enumerate().filter(|(_, x)| **x != 0) {
        add_shift(model, &mut result, b, degree, factor)?;
    }
    Ok(result)
}

fn add_shift(
    model: &mut Model,
    out: &mut Polynomial,
    source: &[i128],
    shift: usize,
    factor: i128,
) -> Result<()> {
    model.spend(source.len())?;
    if !source.is_empty() {
        out.resize(out.len().max(shift + source.len()), 0);
        for (a, &b) in out[shift..].iter_mut().zip(source) {
            *a = add(*a, mul(b, factor)?)?;
        }
        normalize(out);
    }
    Ok(())
}

fn logarithmic_traces(model: &mut Model, p: &[i128], n: usize) -> Result<Vec<i128>> {
    let mut s = vec![0i128; n + 1];
    if p.first() != Some(&1) {
        return Err(inconsistent());
    }
    for i in 1..=n {
        let mut count = mul(-(i as i128), p.get(i).copied().unwrap_or(0))?;
        model.spend(i.min(p.len()))?;
        for (j, &a) in p.iter().enumerate().take(i).skip(1) {
            count = sub(count, mul(a, s[i - j])?)?;
        }
        s[i] = count;
    }
    Ok(s.split_off(1))
}

fn zeta_coefficient(model: &mut Model, p: &[i128], n: usize) -> Result<i128> {
    let mut s = vec![0i128; n + 1];
    s[0] = 1;
    for i in 1..=n {
        model.spend(i.min(p.len()))?;
        let mut c = 0;
        for (j, &a) in p.iter().enumerate().take(i + 1).skip(1) {
            c = sub(c, mul(a, s[i - j])?)?;
        }
        s[i] = c;
    }
    Ok(s[n])
}

fn primitive_count(traces: &[i128], n: usize) -> Result<i128> {
    if n == 0 {
        return Ok(0);
    }
    let mut count = 0;
    for (d, &trace) in traces.iter().enumerate().take(n + 1).skip(1) {
        if n.is_multiple_of(d) {
            count = add(count, mul(mobius(n / d), trace)?)?;
        }
    }
    let n = n as i128;
    if count < 0 || count % n != 0 {
        return Err(inconsistent());
    }
    Ok(count / n)
}

fn mobius(mut n: usize) -> i128 {
    let mut sign = 1;
    let mut p = 2;
    while p * p <= n {
        if n.is_multiple_of(p) {
            n /= p;
            if n.is_multiple_of(p) {
                return 0;
            }
            sign = -sign;
        }
        p += 1;
    }
    if n > 1 {
        -sign
    } else {
        sign
    }
}

fn add(a: i128, b: i128) -> Result<i128> {
    a.checked_add(b).ok_or_else(arithmetic)
}
fn sub(a: i128, b: i128) -> Result<i128> {
    a.checked_sub(b).ok_or_else(arithmetic)
}
fn mul(a: i128, b: i128) -> Result<i128> {
    a.checked_mul(b).ok_or_else(arithmetic)
}

fn arithmetic() -> Error {
    Error::plain("Ihara polynomial or trace arithmetic exceeds the signed 128-bit exact range")
}
fn core_limit() -> Error {
    Error::plain("explicit Schottky quotient graph supports at most 256 edges")
}
fn inconsistent() -> Error {
    Error::plain(
        "weighted nonbacktracking determinant, Bass determinant or primitive count disagrees",
    )
}
