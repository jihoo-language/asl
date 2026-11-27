use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;

pub(crate) trait CoefficientRing {
    type Element: Clone + Eq;

    fn cells(&self) -> usize;
    fn zero(&self) -> Self::Element;
    fn one(&self) -> Self::Element;
    fn is_zero(&self, value: &Self::Element) -> bool;
    fn is_unit(&self, value: &Self::Element) -> bool;
    fn add(&self, model: &mut Model, a: &Self::Element, b: &Self::Element)
        -> Result<Self::Element>;
    fn subtract(
        &self,
        model: &mut Model,
        a: &Self::Element,
        b: &Self::Element,
    ) -> Result<Self::Element>;
    fn multiply(
        &self,
        model: &mut Model,
        a: &Self::Element,
        b: &Self::Element,
    ) -> Result<Self::Element>;
    fn inverse_unit(&self, model: &mut Model, value: &Self::Element) -> Result<Self::Element>;
}

pub(crate) fn storage(entries: usize, width: usize) -> Result<()> {
    if entries.saturating_mul(width) > MAX_PIVOT_CELLS {
        return Err(Error::plain(
            "Hasse–Witt calculation exceeds the coefficient storage limit",
        ));
    }
    Ok(())
}

#[derive(Clone, Eq, PartialEq)]
pub(crate) struct Matrix<T> {
    size: usize,
    entries: Vec<T>,
}

impl<T> Matrix<T> {
    pub(crate) fn from_entries(size: usize, entries: Vec<T>) -> Result<Self> {
        if size.checked_mul(size) != Some(entries.len()) {
            return Err(Error::plain(
                "Hasse–Witt matrix has the wrong number of entries",
            ));
        }
        Ok(Self { size, entries })
    }

    pub(crate) fn size(&self) -> usize {
        self.size
    }

    pub(crate) fn entries(&self) -> &[T] {
        &self.entries
    }

    fn row(&self, i: usize) -> &[T] {
        &self.entries[i * self.size..(i + 1) * self.size]
    }

    pub(crate) fn entry(&self, i: usize, j: usize) -> &T {
        &self.entries[i * self.size + j]
    }
}

pub(crate) fn identity<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    size: usize,
) -> Result<Matrix<R::Element>> {
    let count = size.saturating_mul(size);
    storage(count, ring.cells())?;
    model.spend(count)?;
    let mut entries = vec![ring.zero(); count];
    for i in 0..size {
        entries[i * (size + 1)] = ring.one();
    }
    Matrix::from_entries(size, entries)
}

pub(crate) fn multiply<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    a: &Matrix<R::Element>,
    b: &Matrix<R::Element>,
) -> Result<Matrix<R::Element>> {
    let n = a.size;
    if n != b.size {
        return Err(Error::plain("Hasse–Witt matrices must have the same size"));
    }
    storage(a.entries.len(), ring.cells())?;
    let mut entries = Vec::with_capacity(a.entries.len());
    for i in 0..n {
        for j in 0..n {
            entries.push(dot(model, ring, a.row(i), (0..n).map(|k| b.entry(k, j)))?);
        }
    }
    Matrix::from_entries(n, entries)
}

fn dot<'a, R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    left: &[R::Element],
    right: impl Iterator<Item = &'a R::Element>,
) -> Result<R::Element>
where
    R::Element: 'a,
{
    model.spend(left.len())?;
    let mut sum = ring.zero();
    for (a, b) in left.iter().zip(right) {
        if !ring.is_zero(a) && !ring.is_zero(b) {
            let product = ring.multiply(model, a, b)?;
            sum = ring.add(model, &sum, &product)?;
        }
    }
    Ok(sum)
}

pub(crate) fn inverse<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    input: &Matrix<R::Element>,
) -> Result<Option<Matrix<R::Element>>> {
    let n = input.size;
    storage(2 * input.entries.len(), ring.cells())?;
    model.spend(input.entries.len())?;
    let mut left = input.entries.clone();
    let mut right = identity(model, ring, n)?.entries;
    for i in 0..n {
        model.spend(n - i)?;
        let Some(pivot) = (i..n).find(|&k| ring.is_unit(&left[k * n + i])) else {
            return Ok(None);
        };
        if pivot != i {
            model.spend(n)?;
            for j in 0..n {
                left.swap(i * n + j, pivot * n + j);
                right.swap(i * n + j, pivot * n + j);
            }
        }
        let unit = ring.inverse_unit(model, &left[i * n + i])?;
        model.spend(n)?;
        for j in 0..n {
            left[i * n + j] = ring.multiply(model, &left[i * n + j], &unit)?;
            right[i * n + j] = ring.multiply(model, &right[i * n + j], &unit)?;
        }
        model.spend(n)?;
        for k in 0..n {
            if k == i || ring.is_zero(&left[k * n + i]) {
                continue;
            }
            model.spend(n)?;
            let scalar = left[k * n + i].clone();
            for j in 0..n {
                let l = ring.multiply(model, &left[i * n + j], &scalar)?;
                let r = ring.multiply(model, &right[i * n + j], &scalar)?;
                left[k * n + j] = ring.subtract(model, &left[k * n + j], &l)?;
                right[k * n + j] = ring.subtract(model, &right[k * n + j], &r)?;
            }
        }
    }
    Ok(Some(Matrix::from_entries(n, right)?))
}

pub(crate) fn characteristic<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    matrix: &Matrix<R::Element>,
) -> Result<Vec<R::Element>> {
    let n = matrix.size;
    storage(5 * (n + 1), ring.cells())?;
    let mut characteristic = vec![ring.one()];
    for start in (0..n).rev() {
        let tail = start + 1;
        let degree = n - start;
        // For [a R; S M], convolve char(M) with [1,-a,-RS,-RMS,...].
        // This is the division-free Berkowitz recurrence over the base ring.
        let mut series = vec![
            ring.one(),
            ring.subtract(model, &ring.zero(), matrix.entry(start, start))?,
        ];
        if degree > 1 {
            model.spend(degree - 1)?;
            let mut image: Vec<_> = (tail..n).map(|i| matrix.entry(i, start).clone()).collect();
            let leading = &matrix.row(start)[tail..];
            for power in 2..=degree {
                let term = dot(model, ring, leading, image.iter())?;
                series.push(ring.subtract(model, &ring.zero(), &term)?);
                if power < degree {
                    image = (tail..n)
                        .map(|i| dot(model, ring, &matrix.row(i)[tail..], image.iter()))
                        .collect::<Result<_>>()?;
                }
            }
        }
        characteristic = (0..=degree)
            .map(|i| {
                dot(
                    model,
                    ring,
                    &characteristic[..i.min(degree - 1) + 1],
                    series[..=i].iter().rev(),
                )
            })
            .collect::<Result<_>>()?;
    }
    Ok(characteristic)
}

pub(crate) type PolynomialDivision<T> = (Vec<T>, Vec<T>);

pub(crate) fn divide_monic<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    dividend: &[R::Element],
    divisor: &[R::Element],
) -> Result<PolynomialDivision<R::Element>> {
    if divisor.first() != Some(&ring.one()) {
        return Err(Error::plain(
            "Hasse–Witt division requires a monic descending coefficient list",
        ));
    }
    if dividend.len() < divisor.len() {
        model.spend(dividend.len())?;
        return Ok((Vec::new(), dividend.to_vec()));
    }
    storage(2 * dividend.len(), ring.cells())?;
    model.spend(dividend.len())?;
    let mut remainder = dividend.to_vec();
    let count = dividend.len() + 1 - divisor.len();
    let mut quotient = Vec::with_capacity(count);
    for i in 0..count {
        let leading = std::mem::replace(&mut remainder[i], ring.zero());
        if !ring.is_zero(&leading) {
            model.spend(divisor.len())?;
            for (j, c) in divisor.iter().enumerate().skip(1) {
                let product = ring.multiply(model, &leading, c)?;
                remainder[i + j] = ring.subtract(model, &remainder[i + j], &product)?;
            }
        }
        quotient.push(leading);
    }
    Ok((quotient, remainder.split_off(count)))
}

pub(crate) fn pow_polynomial<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    base: &[R::Element],
    mut exponent: u64,
    bound: usize,
) -> Result<Vec<R::Element>> {
    if bound == 0 {
        return Ok(Vec::new());
    }
    // The recursive product retains several coefficient vectors of this size.
    storage(bound.saturating_mul(10), ring.cells())?;
    let mut answer = vec![ring.one()];
    let mut power = base[..base.len().min(bound)].to_vec();
    model.spend(power.len())?;
    while exponent != 0 {
        model.spend(1)?;
        if !exponent.is_multiple_of(2) {
            answer = product(model, ring, &answer, &power, bound)?;
        }
        exponent >>= 1;
        if exponent != 0 {
            power = product(model, ring, &power, &power, bound)?;
        }
    }
    Ok(answer)
}

pub(crate) fn multiply_polynomials<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    a: &[R::Element],
    b: &[R::Element],
) -> Result<Vec<R::Element>> {
    if a.is_empty() || b.is_empty() {
        return Ok(Vec::new());
    }
    let bound = a
        .len()
        .checked_add(b.len())
        .ok_or_else(|| Error::plain("polynomial product degree is out of range"))?
        - 1;
    storage(bound.saturating_mul(10), ring.cells())?;
    product(model, ring, a, b, bound)
}

fn product<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    a: &[R::Element],
    b: &[R::Element],
    bound: usize,
) -> Result<Vec<R::Element>> {
    let a = &a[..a.len().min(bound)];
    let b = &b[..b.len().min(bound)];
    let short = a.len().min(b.len());
    if short == 0 {
        return Ok(Vec::new());
    }
    let count = (a.len() + b.len() - 1).min(bound);
    model.spend(count)?;
    let long = a.len().max(b.len());
    if short <= 24 || long > 2 * short {
        let mut result = vec![ring.zero(); count];
        for (i, x) in a.iter().enumerate() {
            let needed = b.len().min(count - i);
            model.spend(1 + needed)?;
            if ring.is_zero(x) {
                continue;
            }
            for (j, y) in b[..needed].iter().enumerate() {
                if !ring.is_zero(y) {
                    let term = ring.multiply(model, x, y)?;
                    result[i + j] = ring.add(model, &result[i + j], &term)?;
                }
            }
        }
        return Ok(result);
    }

    let half = long.div_ceil(2);
    let (a0, a1) = a.split_at(half.min(a.len()));
    let (b0, b1) = b.split_at(half.min(b.len()));
    let low = product(model, ring, a0, b0, count)?;
    let mid_count = count - half;
    // The high product is also subtracted in the middle term, so retain all
    // its coefficients below count-half even if its final placement is later.
    let high = product(model, ring, a1, b1, mid_count)?;
    let aa = polynomial_sum(model, ring, a0, a1)?;
    let bb = polynomial_sum(model, ring, b0, b1)?;
    let mut middle = product(model, ring, &aa, &bb, mid_count)?;
    middle.resize(mid_count, ring.zero());
    model.spend(low.len().min(mid_count) + high.len().min(mid_count))?;
    for (t, x) in middle.iter_mut().zip(&low) {
        *t = ring.subtract(model, t, x)?;
    }
    for (t, x) in middle.iter_mut().zip(&high) {
        *t = ring.subtract(model, t, x)?;
    }
    let mut output = low;
    output.resize(count, ring.zero());
    model.spend(mid_count + count.saturating_sub(2 * half))?;
    for (out, t) in output[half..].iter_mut().zip(middle) {
        *out = ring.add(model, out, &t)?;
    }
    if 2 * half < count {
        for (out, t) in output[2 * half..].iter_mut().zip(high) {
            *out = ring.add(model, out, &t)?;
        }
    }
    Ok(output)
}

fn polynomial_sum<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    a: &[R::Element],
    b: &[R::Element],
) -> Result<Vec<R::Element>> {
    let (long, short) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    model.spend(long.len())?;
    let mut result = long.to_vec();
    for (out, t) in result.iter_mut().zip(short) {
        *out = ring.add(model, out, t)?;
    }
    Ok(result)
}

pub(crate) fn hasse_witt<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    f: &[R::Element],
    q: u64,
) -> Result<Matrix<R::Element>> {
    if f.len() < 2 || !ring.is_unit(f.last().unwrap()) {
        return Err(Error::plain(
            "Hasse–Witt equation requires positive degree and a unit leading coefficient",
        ));
    }
    if q.is_multiple_of(2) {
        return Err(Error::plain(
            "Hasse–Witt coefficient power requires an odd residue size",
        ));
    }
    let genus = (f.len() - 2) / 2;
    if genus == 0 || q == 1 {
        return identity(model, ring, genus);
    }
    let entries = genus.saturating_mul(genus);
    storage(entries, ring.cells())?;
    model.spend(entries + f.len())?;
    let exponent = (q - 1) / 2;
    let lowest = f
        .iter()
        .position(|a| !ring.is_zero(a))
        .expect("polynomial with a unit leading coefficient");
    let f = &f[lowest..];
    let degree = (f.len() - 1) as u128 * u128::from(exponent);
    let shift = lowest as u128 * u128::from(exponent);
    let index = |i: usize, j: usize| {
        (u128::from(q) * j as u128)
            .checked_sub(i as u128)
            .and_then(|d| d.checked_sub(shift))
            .filter(|&d| d <= degree)
    };
    let mut min = None;
    let mut max = 0;
    for i in 1..=genus {
        for j in 1..=genus {
            if let Some(k) = index(i, j) {
                min = Some(min.map_or(k, |m: u128| m.min(k)));
                max = max.max(k);
            }
        }
    }
    let Some(min) = min else {
        return Matrix::from_entries(genus, vec![ring.zero(); entries]);
    };
    let reverse = degree - min < max;
    let needed = if reverse { degree - min } else { max };
    let bound = usize::try_from(needed + 1)
        .map_err(|_| Error::plain("Hasse–Witt coefficient index is out of range"))?;
    let reversed: Vec<_>;
    let base = if reverse {
        model.spend(f.len())?;
        reversed = f.iter().rev().take(bound).cloned().collect();
        reversed.as_slice()
    } else {
        f
    };
    let values = pow_polynomial(model, ring, base, exponent, bound)?;
    let mut coefficients = Vec::with_capacity(entries);
    model.spend(entries)?;
    for i in 1..=genus {
        for j in 1..=genus {
            let value = index(i, j)
                .map(|k| if reverse { degree - k } else { k })
                .and_then(|k| values.get(k as usize))
                .cloned()
                .unwrap_or_else(|| ring.zero());
            coefficients.push(value);
        }
    }
    Matrix::from_entries(genus, coefficients)
}
