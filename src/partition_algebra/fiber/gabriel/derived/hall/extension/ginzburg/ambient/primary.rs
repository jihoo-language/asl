use crate::error::{Error, Result};
use crate::model::Model;

use super::distribution::{accumulate_convolution, bound, exact_divide, gl, power, unit};
use super::Field;

type Series = Vec<Vec<u128>>;

pub(super) fn counts(model: &mut Model, field: &Field, n: usize, third: u32) -> Result<Vec<u128>> {
    if n == 0 {
        return Ok(unit(field));
    }
    let q = u128::from(field.order());
    power(q, n * n)?;
    let census = irreducible_census(model, field, n, third)?;
    let groups = (0..=n).map(|d| gl(q, d)).collect::<Result<Vec<_>>>()?;
    let decompositions = (0..=n)
        .map(|total| {
            (0..=total)
                .map(|left| {
                    let denominator = groups[left]
                        .checked_mul(groups[total - left])
                        .ok_or_else(bound)?;
                    exact_divide(groups[total], denominator)
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let arithmetic = Primary {
        field,
        groups,
        decompositions,
    };
    let mut product = arithmetic.one();
    for (h, row) in census.iter().enumerate().skip(1) {
        for (trace, &number) in row.iter().enumerate().filter(|(_, n)| **n != 0) {
            let local = arithmetic.local(model, h, trace as u32)?;
            let power = arithmetic.combine_copies(model, local, number)?;
            product = arithmetic.combine(model, &product, &power)?;
        }
    }
    Ok(product.pop().unwrap())
}

struct Primary<'a> {
    field: &'a Field,
    groups: Vec<u128>,
    decompositions: Vec<Vec<u128>>,
}

impl Primary<'_> {
    fn one(&self) -> Series {
        let mut result = vec![vec![0; self.field.order() as usize]; self.groups.len()];
        result[0][0] = 1;
        result
    }

    fn local(&self, model: &mut Model, h: usize, trace: u32) -> Result<Series> {
        let q = u128::from(self.field.order());
        let p = self.field.characteristic() as usize;
        let mut result = self.one();
        for (size, row) in result.iter_mut().enumerate().skip(h).step_by(h) {
            let n = size / h;
            let mut contribution = 0u128;
            partitions(n, n, &mut Vec::new(), &mut |partition| {
                model.spend(1 + n * n)?;
                let centralizer = centralizer(q, h, partition)?;
                let orbit = exact_divide(self.groups[size], centralizer)?;
                contribution = contribution.checked_add(orbit).ok_or_else(bound)?;
                Ok(())
            })?;
            let value = self.field.mul((n % p) as u32, trace);
            row[value as usize] = contribution;
        }
        Ok(result)
    }

    fn combine(&self, model: &mut Model, a: &Series, b: &Series) -> Result<Series> {
        let mut result = vec![vec![0; self.field.order() as usize]; self.groups.len()];
        let active_a: Vec<_> = a
            .iter()
            .enumerate()
            .filter(|(_, row)| row.iter().any(|&c| c != 0))
            .collect();
        let active_b: Vec<_> = b
            .iter()
            .enumerate()
            .filter(|(_, row)| row.iter().any(|&c| c != 0))
            .collect();
        model.spend(self.groups.len() * self.field.order() as usize)?;
        for &(i, x) in &active_a {
            for &(j, y) in &active_b {
                if let Some(out) = result.get_mut(i + j) {
                    // The two primary spectra are disjoint. Each matrix on the
                    // direct sum determines its two primary subspaces uniquely.
                    accumulate_convolution(
                        model,
                        self.field,
                        out,
                        x,
                        y,
                        self.decompositions[i + j][i],
                    )?;
                }
            }
        }
        Ok(result)
    }

    fn combine_copies(
        &self,
        model: &mut Model,
        mut local: Series,
        mut copies: usize,
    ) -> Result<Series> {
        let mut result = self.one();
        while copies != 0 {
            if copies & 1 != 0 {
                result = self.combine(model, &result, &local)?;
            }
            copies >>= 1;
            if copies != 0 {
                local = self.combine(model, &local, &local)?;
            }
        }
        Ok(result)
    }
}

fn partitions(
    n: usize,
    largest: usize,
    parts: &mut Vec<usize>,
    visit: &mut impl FnMut(&[usize]) -> Result<()>,
) -> Result<()> {
    if n == 0 {
        return visit(parts);
    }
    for part in 1..=n.min(largest) {
        parts.push(part);
        partitions(n - part, part, parts, visit)?;
        parts.pop();
    }
    Ok(())
}

fn centralizer(q: u128, h: usize, partition: &[usize]) -> Result<u128> {
    // |End| has exponent h sum (lambda'_i)^2. For each run of equal
    // blocks multiply by prod_{j=1..m}(1−q^(−hj)).
    let dimension = (1..=partition[0])
        .map(|i| partition.iter().filter(|&&p| p >= i).count().pow(2))
        .sum::<usize>();
    let mut denominator_exponent = 0;
    let mut factors = 1u128;
    for run in partition.chunk_by(|a, b| a == b) {
        for j in 1..=run.len() {
            denominator_exponent += j;
            factors = factors
                .checked_mul(power(q, h * j)? - 1)
                .ok_or_else(bound)?;
        }
    }
    factors
        .checked_mul(power(q, h * (dimension - denominator_exponent))?)
        .ok_or_else(bound)
}

fn irreducible_census(
    model: &mut Model,
    field: &Field,
    n: usize,
    third: u32,
) -> Result<Vec<Vec<usize>>> {
    const POLYNOMIAL_LIMIT: u128 = 100_000;
    let q = field.order() as usize;
    if power(q as u128, n)? > POLYNOMIAL_LIMIT {
        return Err(Error::plain("partition quiver primary census exceeds 100,000 monic polynomial candidates at one degree"));
    }
    let mut counts = vec![vec![0; q]; n + 1];
    let mut totals = vec![0; n + 1];
    for degree in 1..=n {
        let mut f = vec![0; degree + 1];
        f[degree] = 1;
        let candidates = q.pow(degree as u32);
        for code in 0..candidates {
            model.spend(degree + 1)?;
            let mut code = code;
            for a in &mut f[..degree] {
                *a = (code % q) as u32;
                code /= q;
            }
            if irreducible(model, field, &f)? {
                let trace = cube_trace(field, &f, third);
                counts[degree][trace as usize] += 1;
                totals[degree] += 1;
            }
        }
        let degree_sum = (1..=degree)
            .filter(|d| degree.is_multiple_of(*d))
            .map(|d| d * totals[d])
            .sum::<usize>();
        if degree_sum != candidates {
            return Err(Error::plain(
                "partition quiver irreducible census failed the finite-field degree identity",
            ));
        }
    }
    Ok(counts)
}

fn cube_trace(field: &Field, f: &[u32], third: u32) -> u32 {
    let degree = f.len() - 1;
    let coefficient = |i| degree.checked_sub(i).map_or(0, |j| f[j]);
    let (a, b, c) = (coefficient(1), coefficient(2), coefficient(3));
    // Newton gives Tr(alpha³/3)=−a³/3+ab−c, including degrees 1 and 2.
    let cube = field.mul(a, field.mul(a, a));
    field.sub(field.sub(field.mul(a, b), c), field.mul(third, cube))
}

fn irreducible(model: &mut Model, field: &Field, f: &[u32]) -> Result<bool> {
    let degree = f.len() - 1;
    if degree == 1 {
        return Ok(true);
    }
    if f[0] == 0 {
        return Ok(false);
    }
    let mut power = vec![0, 1];
    // Any proper factorization contains an irreducible factor of degree
    // at most floor(degree/2); the degree-j factors divide X^(q^j)−X.
    for _ in 1..=degree / 2 {
        power = pow_mod(model, field, &power, field.order(), f)?;
        let mut difference = power.clone();
        difference.resize(difference.len().max(2), 0);
        difference[1] = field.sub(difference[1], 1);
        trim(&mut difference);
        if !coprime(model, field, f, &difference)? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn trim(f: &mut Vec<u32>) {
    while f.last() == Some(&0) {
        f.pop();
    }
}

fn reduce(model: &mut Model, field: &Field, mut a: Vec<u32>, b: &[u32]) -> Result<Vec<u32>> {
    let leading = field.inverse(model, *b.last().unwrap())?;
    while a.len() >= b.len() {
        model.spend(b.len())?;
        let offset = a.len() - b.len();
        let k = field.mul(*a.last().unwrap(), leading);
        for (x, &y) in a[offset..].iter_mut().zip(b) {
            *x = field.sub(*x, field.mul(k, y));
        }
        trim(&mut a);
    }
    Ok(a)
}

fn coprime(model: &mut Model, field: &Field, a: &[u32], b: &[u32]) -> Result<bool> {
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    while b.len() > 1 {
        let r = reduce(model, field, a, &b)?;
        a = b;
        b = r;
    }
    Ok(b.len() == 1)
}

fn multiply_mod(
    model: &mut Model,
    field: &Field,
    a: &[u32],
    b: &[u32],
    modulus: &[u32],
) -> Result<Vec<u32>> {
    if a.is_empty() || b.is_empty() {
        return Ok(Vec::new());
    }
    let mut product = vec![0; a.len() + b.len() - 1];
    model.spend(a.len() * b.len())?;
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            product[i + j] = field.add(product[i + j], field.mul(x, y));
        }
    }
    trim(&mut product);
    reduce(model, field, product, modulus)
}

fn pow_mod(
    model: &mut Model,
    field: &Field,
    a: &[u32],
    mut k: u32,
    modulus: &[u32],
) -> Result<Vec<u32>> {
    let mut current = a.to_vec();
    let mut value = vec![1];
    while k != 0 {
        if k & 1 != 0 {
            value = multiply_mod(model, field, &value, &current, modulus)?;
        }
        k >>= 1;
        if k != 0 {
            current = multiply_mod(model, field, &current, &current, modulus)?;
        }
    }
    Ok(value)
}
