use crate::binary_field::BinaryField;
use crate::error::{Error, Result};
use crate::model::{Model, TypeId};
use crate::rv_hodge::{self, PivotOrder};
use crate::rv_laurent::{self, Term};

#[derive(Clone, Copy)]
pub(crate) enum Boundary {
    Exact,
    Closure,
}

impl Boundary {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Exact => "rvInIwADLV",
            Self::Closure => "rvInIwADLVClosure",
        }
    }
}

pub(crate) fn check_window(model: &mut Model, ty: TypeId, target: &[i128]) -> Result<()> {
    let size = rv_laurent::dimensions(model, ty)?.0;
    if target.len() != size {
        return Err(Error::plain(
            "the affine window must have one entry per Margolis dimension",
        ));
    }
    validate_window(model, target)
}

pub(crate) fn validate_window(model: &mut Model, target: &[i128]) -> Result<()> {
    let size = target.len();
    model.spend(size)?;
    let mut seen = vec![false; size];
    for &value in target {
        let residue = value.rem_euclid(size as i128) as usize;
        if std::mem::replace(&mut seen[residue], true) {
            return Err(Error::plain(
                "the affine window residues must form a permutation",
            ));
        }
    }
    Ok(())
}

fn window(
    model: &mut Model,
    field: BinaryField,
    size: usize,
    terms: &[Term],
    shift: i64,
    name: &str,
) -> Result<Vec<i128>> {
    model.spend(terms.len())?;
    let common = common_power(terms);
    let shifted;
    let terms = if common == 0 {
        terms
    } else {
        shifted = terms
            .iter()
            .map(|&term| Term {
                power: term.power - common,
                ..term
            })
            .collect::<Vec<_>>();
        &shifted
    };
    let offset = i128::from(common) - i128::from(shift);
    let mut window = vec![0; size];
    for pivot in rv_hodge::pivots(model, field, size, terms, name, PivotOrder::Iwahori)? {
        // The pivot minimizes n*valuation+i-j. The permitted row and column
        // operations therefore lie in the lower-triangular Iwahori mod u.
        window[pivot.column] = pivot.row as i128 + size as i128 * (pivot.exponent + offset);
    }
    Ok(window)
}

// Removing a common power of u shrinks the exponents of a central
// translate before Bareiss forms products of entries.
fn common_power(terms: &[Term]) -> i64 {
    let mut minimum = i64::MAX;
    let mut maximum = i64::MIN;
    for term in terms {
        minimum = minimum.min(term.power);
        maximum = maximum.max(term.power);
    }
    if terms.is_empty() {
        0
    } else if minimum > 0 {
        minimum
    } else {
        maximum.min(0)
    }
}

pub(crate) fn standard(
    model: &mut Model,
    ty: TypeId,
    terms: &[Term],
    name: &str,
) -> Result<Vec<i128>> {
    let grade = rv_hodge::square(model, ty, name)?;
    let size = rv_laurent::dimensions(model, ty)?.0;
    window(model, grade.field, size, terms, 0, name)
}

pub(crate) fn position(model: &mut Model, ty: TypeId, x: &[Term], y: &[Term]) -> Result<Vec<i128>> {
    let name = "rvIwPosition";
    let grade = rv_hodge::square(model, ty, name)?;
    let size = rv_laurent::dimensions(model, ty)?.0;
    let (transformed, shift) = rv_hodge::adjugate_action(model, grade.field, size, x, y, name)?;
    window(model, grade.field, size, &transformed, shift, name)
}

pub(crate) fn at(
    model: &mut Model,
    ty: TypeId,
    a: &[Term],
    x: &[Term],
    step: usize,
    name: &str,
) -> Result<Vec<i128>> {
    let grade = rv_hodge::square(model, ty, name)?;
    let size = rv_laurent::dimensions(model, ty)?.0;
    let conjugate = rv_laurent::frobenius(model, grade.field, x, step)?;
    let image = rv_laurent::product(model, grade.field, size, size, a, &conjugate)?;
    let (transformed, shift) =
        rv_hodge::adjugate_action(model, grade.field, size, x, &image, name)?;
    window(model, grade.field, size, &transformed, shift, name)
}

fn overflow() -> Error {
    Error::plain("affine rank calculation exceeds the signed 128-bit range")
}

pub(crate) fn length(model: &mut Model, window: &[i128]) -> Result<u128> {
    let modulus = window.len() as i128;
    let mut length = 0u128;
    for (index, &first) in window.iter().enumerate() {
        model.spend(window.len() - index - 1)?;
        for &second in &window[index + 1..] {
            let difference = second.checked_sub(first).ok_or_else(overflow)?;
            let term = difference.div_euclid(modulus).unsigned_abs();
            length = length.checked_add(term).ok_or_else(overflow)?;
        }
    }
    Ok(length)
}

fn checked_sum(terms: impl IntoIterator<Item = i128>) -> Result<i128> {
    terms.into_iter().try_fold(0i128, |sum, entry| {
        sum.checked_add(entry).ok_or_else(overflow)
    })
}

pub(crate) fn component(window: &[i128]) -> Result<i128> {
    checked_sum(
        window
            .iter()
            .map(|&entry| entry.div_euclid(window.len() as i128)),
    )
}

// For fixed j and i=residue+n*k, the affine rank is sum max(0,b-k).
fn breakpoints(window: &[i128], j: usize, residue: i128) -> Result<Vec<i128>> {
    let modulus = window.len() as i128;
    let mut result = Vec::with_capacity(window.len());
    for (index, &entry) in window.iter().enumerate() {
        let offset = i128::from(index < j) - i128::from(entry.rem_euclid(modulus) < residue);
        let point = entry
            .div_euclid(modulus)
            .checked_add(offset)
            .ok_or_else(overflow)?;
        result.push(point);
    }
    result.sort_unstable();
    Ok(result)
}

struct RankTail<'a> {
    remaining: &'a [i128],
    sum: i128,
}

impl<'a> RankTail<'a> {
    fn new(remaining: &'a [i128]) -> Result<Self> {
        Ok(Self {
            remaining,
            sum: checked_sum(remaining.iter().copied())?,
        })
    }

    fn evaluate(&mut self, point: i128) -> Result<i128> {
        while let Some((&first, tail)) = self.remaining.split_first() {
            if first > point {
                break;
            }
            self.sum = self.sum.checked_sub(first).ok_or_else(overflow)?;
            self.remaining = tail;
        }
        let subtraction = (self.remaining.len() as i128)
            .checked_mul(point)
            .ok_or_else(overflow)?;
        self.sum.checked_sub(subtraction).ok_or_else(overflow)
    }
}

fn below(model: &mut Model, actual: &[i128], target: &[i128]) -> Result<bool> {
    let size = actual.len();
    if size <= 1 {
        return Ok(actual == target);
    }
    if component(actual)? != component(target)? {
        return Ok(false);
    }

    for j in 0..size {
        for residue in 0..size as i128 {
            model.spend(4 * size * (1 + size.ilog2() as usize))?;
            let actual = breakpoints(actual, j, residue)?;
            let bound = breakpoints(target, j, residue)?;
            let mut left = RankTail::new(&actual)?;
            let mut right = RankTail::new(&bound)?;
            // The difference is affine between breakpoints, is constant
            // below the smallest one, and vanishes above the largest.
            while let Some(&point) = left
                .remaining
                .first()
                .into_iter()
                .chain(right.remaining.first())
                .min()
            {
                if left.evaluate(point)? > right.evaluate(point)? {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}

impl Boundary {
    pub(crate) fn contains(
        self,
        model: &mut Model,
        actual: &[i128],
        target: &[i128],
    ) -> Result<bool> {
        model.spend(actual.len())?;
        match self {
            Self::Exact => Ok(actual == target),
            Self::Closure => below(model, actual, target),
        }
    }
}
