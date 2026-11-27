use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::rv_lens::{self, Rational as Scalar};

use super::super::factor;
use super::super::linear::{extend, Rect, Vector};
use super::{
    check_report, check_weights, fusion, Convolution, Error, Field, Model, Result, Square,
};

pub(super) struct Part {
    pub(super) shift: usize,
    pub(super) support: Vec<usize>,
}

pub(super) struct Split {
    pub(super) parts: Vec<Part>,
    pub(super) weights: Vec<usize>,
    pub(super) forward: Square,
    pub(super) inverse: Square,
    pub(super) action: Square,
}

impl Split {
    pub(super) fn compute(model: &mut Model, convolution: &Convolution) -> Result<Self> {
        let Some([a, b]) = convolution.intervals else {
            return Err(Error::plain(
                "LG convolution split requires consecutive circular intervals",
            ));
        };
        let field = Field::new(model, convolution.d)?;
        let d = field.d;
        let n = convolution.weights.len();
        let cells = n.saturating_mul(n).saturating_mul(field.cells());
        if cells.saturating_mul(12) > crate::complex::MAX_PIVOT_CELLS {
            return Err(Error::plain(
                "LG convolution split exceeds the coefficient storage limit",
            ));
        }
        let mut supports: Vec<Vec<_>> = fusion::channels(d, a, b)
            .map(|(start, label)| (0..=label).map(|i| (start + i) % d).collect())
            .collect();
        let stable = supports.iter().map(Vec::len).sum::<usize>();
        let free = n
            .checked_sub(stable)
            .filter(|r| r.is_multiple_of(d))
            .ok_or_else(|| {
                Error::plain("LG convolution interval candidates have incorrect total rank")
            })?
            / d;
        supports.extend((0..free).map(|_| (0..d).collect()));
        supports.sort_by_key(|support| std::cmp::Reverse(support.len()));

        let x = &convolution.powers[1];
        let mut x1 = Rect::zero(n, n);
        let mut x0 = Rect::zero(n, n);
        for (i, &w) in convolution.weights.iter().enumerate() {
            for (j, &v) in convolution.weights.iter().enumerate() {
                x1.rows[i][j] = x.entry(i, j).clone();
                if w == v + 1 {
                    x0.rows[i][j] = x.entry(i, j).clone();
                }
            }
        }
        let mut parts = Vec::with_capacity(supports.len());
        let mut weights = Vec::with_capacity(n);
        let mut columns = Vec::with_capacity(n);
        let mut basis = vec![None; n];
        let mut action = vec![Scalar::zero(); n * n];
        for support in supports {
            let length = support.len();
            let (relation, kernel) = annihilator(model, &field, convolution, &support)?;
            let (shift, mut head) = choose_head(
                model,
                &field,
                &convolution.weights,
                &x0,
                kernel.as_ref(),
                &mut basis,
                length,
            )?;
            let start = columns.len();
            for i in 0..length {
                weights.push(shift + i);
                columns.push(head);
                if i + 1 < length {
                    head = x1.apply(model, &field, columns.last().unwrap())?;
                    action[(start + i + 1) * n + start + i] = Scalar::integer(1);
                } else {
                    for (k, value) in relation.iter().enumerate() {
                        action[(start + k) * n + start + i] = rv_lens::negative(model, value)?;
                    }
                    break;
                }
            }
            parts.push(Part { shift, support });
        }

        if columns.len() != n
            || basis.iter().any(Option::is_none)
            || weights.iter().sum::<usize>() != convolution.weights.iter().sum::<usize>()
        {
            return Err(Error::plain(
                "LG convolution split does not have a constant unit determinant",
            ));
        }
        let forward = Matrix::from_entries(
            n,
            (0..n)
                .flat_map(|row| columns.iter().map(move |c| c[row].clone()))
                .collect(),
        )?;
        let action = Matrix::from_entries(n, action)?;
        check_weights(&forward, &convolution.weights, &weights, 0)?;
        check_weights(&action, &weights, &weights, 1)?;
        // Each equality at z=1 is between matrices with the same homogeneous
        // exponents, so it certifies the full equality over K[z].
        if hw_algebra::multiply(model, &field, x, &forward)?
            != hw_algebra::multiply(model, &field, &forward, &action)?
        {
            return Err(Error::plain(
                "LG convolution cyclic generators failed the polynomial intertwining equation",
            ));
        }
        let inverse = hw_algebra::inverse(model, &field, &forward)?
            .ok_or_else(|| Error::plain("LG convolution split basis is singular"))?;
        check_weights(&inverse, &weights, &convolution.weights, 0)?;
        let identity = &convolution.powers[0];
        if hw_algebra::multiply(model, &field, &inverse, &forward)? != *identity
            || hw_algebra::multiply(model, &field, &forward, &inverse)? != *identity
        {
            return Err(Error::plain(
                "LG convolution split inverse failed verification",
            ));
        }
        model.allocate(cells * 3 + n + parts.iter().map(|p| p.support.len() + 4).sum::<usize>())?;
        Ok(Self {
            parts,
            weights,
            forward,
            inverse,
            action,
        })
    }

    pub(super) fn report(&self, d: usize) -> Result<String> {
        let parts: Vec<_> = self
            .parts
            .iter()
            .map(|p| (p.shift, &p.support, p.support.len() == d))
            .collect();
        let report=format!("(shift,support,free)={parts:?} companion_weights={:?} intertwining=verified determinant=constant-unit inverse=verified",
            self.weights);
        check_report(&report)?;
        Ok(report)
    }
}

fn annihilator(
    model: &mut Model,
    field: &Field,
    convolution: &Convolution,
    support: &[usize],
) -> Result<(Vector, Option<Rect>)> {
    if support.len() == field.d {
        let mut coefficients = vec![Scalar::zero(); field.d];
        coefficients[0] = Scalar::integer(-1);
        return Ok((coefficients, None));
    }
    let n = convolution.weights.len();
    let polynomial = factor(model, field.d, support.iter().copied())?;
    let mut coefficients: Vec<_> = polynomial
        .iter()
        .map(|c| rv_lens::fraction(model, c, 1))
        .collect::<Result<_>>()?;
    let mut image = Rect::zero(n, n);
    for (scalar, power) in coefficients.iter().zip(&convolution.powers) {
        if scalar.is_zero() {
            continue;
        }
        model.spend(n * n)?;
        for (value, entry) in image.rows.iter_mut().flatten().zip(power.entries()) {
            let term = field.multiply(model, scalar, entry)?;
            *value = field.add(model, value, &term)?;
        }
    }
    coefficients.pop();
    Ok((coefficients, Some(image)))
}

fn choose_head(
    model: &mut Model,
    field: &Field,
    weights: &[usize],
    x0: &Rect,
    annihilator: Option<&Rect>,
    basis: &mut Vec<Option<Vector>>,
    length: usize,
) -> Result<(usize, Vector)> {
    let n = weights.len();
    for degree in 0..=weights.iter().copied().max().unwrap_or(0) {
        let allowed: Vec<_> = weights
            .iter()
            .enumerate()
            .filter(|&(_, w)| *w <= degree)
            .map(|(i, _)| i)
            .collect();
        if !allowed.iter().any(|&i| weights[i] == degree) {
            continue;
        }
        let rows = annihilator.map_or(0, |m| m.rows.len());
        let mut restricted = Rect::zero(rows, allowed.len());
        if let Some(annihilator) = annihilator {
            model.spend(rows * allowed.len())?;
            for (out, row) in restricted.rows.iter_mut().zip(&annihilator.rows) {
                for (value, &i) in out.iter_mut().zip(&allowed) {
                    *value = row[i].clone();
                }
            }
        }
        let pivots = restricted.echelon(model, field)?;
        for v in restricted.kernel(model, &pivots)? {
            let mut top = vec![Scalar::zero(); n];
            let mut full = top.clone();
            for (&i, c) in allowed.iter().zip(v) {
                if weights[i] == degree {
                    top[i] = c.clone();
                }
                full[i] = c;
            }
            if top.iter().all(Scalar::is_zero) {
                continue;
            }
            model.spend(n * n)?;
            let mut trial = basis.clone();
            let mut size = 0;
            while size < length && extend(model, field, &mut trial, top.clone())? {
                top = x0.apply(model, field, &top)?;
                size += 1;
            }
            if size == length {
                if top.iter().any(|c| !c.is_zero()) {
                    return Err(Error::plain(
                        "LG convolution head has an incorrect chain length at z=0",
                    ));
                }
                *basis = trial;
                return Ok((degree, full));
            }
        }
    }
    Err(Error::plain(
        "LG convolution could not certify an independent homogeneous cyclic head",
    ))
}
