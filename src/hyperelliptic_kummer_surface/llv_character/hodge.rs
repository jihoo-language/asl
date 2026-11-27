use std::fmt::Write;

use crate::error::{Error, Result};
use crate::model::{Model, MAX_SHOWN_BYTES};

use super::{
    arithmetic, decimal, guard, report_limit, weight_exceeds, weyl, Magnitude, Oscillators, Signed,
    LARGEST_RHO, RANK, WEYL_LIMBS,
};

const SECOND_RHO: usize = RANK - 2;

type Charge = (u128, u128);

pub(super) fn coefficient(
    model: &mut Model,
    n: u128,
    p: i128,
    q: i128,
    shape: &[u128],
) -> Result<String> {
    let charge = coordinates(n, p, q);
    let irrep = if let Some(charge) = charge.filter(|&weight| supports(shape, weight)) {
        Specialization::new(model)?
            .character(model, shape)?
            .weight(charge)
    } else {
        Magnitude::default()
    };
    let contribution = if !irrep.is_zero()
        && !weight_exceeds(shape, n)
        && charge.is_some_and(|charge| physical(charge, n))
    {
        let multiplicity = Oscillators::new(model, n)?.multiplicity(model, shape)?;
        weighted(irrep, multiplicity)?
    } else {
        Magnitude::default()
    };
    let irrep = decimal(model, irrep)?;
    let contribution = decimal(model, contribution)?;
    Ok(format!(
        "rvHilbLLVHodge[{n}, {p}, {q}]{shape:?} = irrep={irrep} contribution={contribution}"
    ))
}

pub(super) fn decompose(model: &mut Model, n: u128, p: i128, q: i128) -> Result<String> {
    let mut report = format!("rvHilbLLVHodgeDecomposition[{n}, {p}, {q}] = components=[");
    let mut total = Magnitude::default();
    if let Some(charge) = coordinates(n, p, q).filter(|&charge| physical(charge, n)) {
        let oscillators = Oscillators::new(model, n)?;
        let specialization = Specialization::new(model)?;
        let mut first = true;
        oscillators.visit_summands(model, |model, shape, multiplicity| {
            if !supports(shape, charge) {
                return Ok(());
            }
            let irrep = specialization.character(model, shape)?.weight(charge);
            if irrep.is_zero() {
                return Ok(());
            }
            let contribution = weighted(irrep, multiplicity)?;
            total.add_multiple(contribution, 1).ok_or_else(arithmetic)?;
            if !first {
                report.push_str(", ");
            }
            first = false;
            let irrep = decimal(model, irrep)?;
            let contribution = decimal(model, contribution)?;
            write!(
                report,
                "({shape:?}, {multiplicity}, {irrep}, {contribution})"
            )
            .unwrap();
            check_report(&report)
        })?;
    }
    let total = decimal(model, total)?;
    write!(report, "] total={total}").unwrap();
    check_report(&report)?;
    Ok(report)
}

pub(super) fn table(model: &mut Model, length: u128, shape: Option<&[u128]>) -> Result<String> {
    let n = usize::try_from(length).map_err(|_| super::work_limit())?;
    let width = n.saturating_mul(2).saturating_add(1);
    let mut diamond = Square::new(width, Magnitude::default())?;
    let specialization = Specialization::new(model)?;
    let mut report = format!("rvHilbLLVHodgeTable[{length}]");
    if let Some(shape) = shape {
        let irrep = specialization.character(model, shape)?;
        add_diamond(model, n, &mut diamond, &irrep, 1)?;
        let multiplicity = if weight_exceeds(shape, length) {
            0
        } else {
            Oscillators::new(model, length)?.multiplicity(model, shape)?
        };
        write!(report, "{shape:?} = irrep=").unwrap();
        append_table(model, &diamond, 1, &mut report)?;
        report.push_str(" contribution=");
        append_table(model, &diamond, multiplicity, &mut report)?;
    } else {
        Oscillators::new(model, length)?.visit_summands(model, |model, shape, multiplicity| {
            let irrep = specialization.character(model, shape)?;
            add_diamond(model, n, &mut diamond, &irrep, multiplicity)
        })?;
        report.push_str(" = hodge=");
        append_table(model, &diamond, 1, &mut report)?;
    }
    Ok(report)
}

fn add_diamond(
    model: &mut Model,
    n: usize,
    diamond: &mut Square<Magnitude>,
    irrep: &Square<Magnitude>,
    multiplicity: i128,
) -> Result<()> {
    model.spend(diamond.cells.len())?;
    for p in 0..diamond.width {
        for q in (p % 2..diamond.width).step_by(2) {
            let a = ((p + q) / 2).abs_diff(n);
            let b = p.abs_diff(q) / 2;
            let weight = irrep.weight((a as u128, b as u128));
            if weight.is_zero() {
                continue;
            }
            let value = weighted(weight, multiplicity)?;
            diamond
                .get_mut(p, q)
                .add_multiple(value, 1)
                .ok_or_else(arithmetic)?;
        }
    }
    Ok(())
}

fn append_table(
    model: &mut Model,
    table: &Square<Magnitude>,
    multiplicity: i128,
    report: &mut String,
) -> Result<()> {
    model.spend(table.cells.len())?;
    report.push('[');
    for (index, row) in table.cells.chunks(table.width).enumerate() {
        if index != 0 {
            report.push_str(", ");
        }
        report.push('[');
        for (column, &value) in row.iter().enumerate() {
            if column != 0 {
                report.push_str(", ");
            }
            report.push_str(&decimal(model, weighted(value, multiplicity)?)?);
            check_report(report)?;
        }
        report.push(']');
    }
    report.push(']');
    check_report(report)
}

fn coordinates(n: u128, p: i128, q: i128) -> Option<Charge> {
    if (p & 1) != (q & 1) {
        return None;
    }
    let mid = p / 2 + q / 2 + (p % 2 + q % 2) / 2;
    let a = if mid < 0 {
        n.checked_add(mid.unsigned_abs())?
    } else {
        n.abs_diff(mid as u128)
    };
    Some((a, p.abs_diff(q) / 2))
}

fn physical((a, b): Charge, n: u128) -> bool {
    a <= n && b <= n - a
}

fn supports(shape: &[u128], (a, b): Charge) -> bool {
    let first = shape.first().copied().unwrap_or(0);
    let second = shape.get(1).copied().unwrap_or(0);
    a.max(b) <= first && a.saturating_add(b) <= first.saturating_add(second)
}

fn weighted(value: Magnitude, multiplicity: i128) -> Result<Magnitude> {
    if multiplicity == 0 {
        return Ok(Magnitude::default());
    }
    if multiplicity == 1 {
        return Ok(value);
    }
    value
        .multiply_u128(multiplicity as u128)
        .ok_or_else(arithmetic)
}

fn check_report(report: &str) -> Result<()> {
    if report.len() > MAX_SHOWN_BYTES {
        Err(report_limit())
    } else {
        Ok(())
    }
}

struct Square<T> {
    width: usize,
    cells: Vec<T>,
}

impl<T: Copy> Square<T> {
    fn new(width: usize, zero: T) -> Result<Self> {
        guard(width.saturating_pow(2).saturating_mul(WEYL_LIMBS))?;
        Ok(Self {
            width,
            cells: vec![zero; width * width],
        })
    }

    fn get(&self, a: usize, b: usize) -> T {
        self.cells[a * self.width + b]
    }

    fn get_mut(&mut self, a: usize, b: usize) -> &mut T {
        &mut self.cells[a * self.width + b]
    }
}

impl Square<Magnitude> {
    fn weight(&self, (a, b): Charge) -> Magnitude {
        if a.max(b) >= self.width as u128 {
            return Magnitude::default();
        }
        self.get(a as usize, b as usize)
    }
}

struct Specialization {
    leading: Magnitude,
    terms: Vec<(i128, i128, Signed)>,
}

impl Specialization {
    fn new(model: &mut Model) -> Result<Self> {
        let denominator = alternant(model, &[], 0)?;
        let lead = denominator.get(LARGEST_RHO, SECOND_RHO);
        if lead.negative || lead.is_zero() {
            return Err(Error::plain(
                "rank-two LLV Weyl denominator has the wrong leading term",
            ));
        }
        let mut terms = Vec::new();
        let radius = LARGEST_RHO as i128;
        for a in -radius..=radius {
            for b in -radius..=radius {
                let coefficient =
                    denominator.get(a.unsigned_abs() as usize, b.unsigned_abs() as usize);
                if coefficient.is_zero() || (a, b) == (radius, SECOND_RHO as i128) {
                    continue;
                }
                terms.push((a, b, coefficient));
            }
        }
        Ok(Self {
            leading: lead.magnitude,
            terms,
        })
    }

    fn character(&self, model: &mut Model, shape: &[u128]) -> Result<Square<Magnitude>> {
        let Some(&first) = shape.first() else {
            return Square::new(1, Magnitude::from_small(1));
        };
        let last = usize::try_from(first).map_err(|_| super::work_limit())?;
        let numerator_cells = last.saturating_add(RANK).saturating_pow(2);
        let weight_cells = last.saturating_add(1).saturating_pow(2);
        guard(
            numerator_cells
                .saturating_add(weight_cells)
                .saturating_mul(WEYL_LIMBS),
        )?;
        model.spend(weight_cells.saturating_mul(self.terms.len()))?;
        let numerator = alternant(model, shape, last)?;
        let mut weights = Square::new(last + 1, Magnitude::default())?;
        for a in (0..=last).rev() {
            for b in (0..=last).rev() {
                let x = a + LARGEST_RHO;
                let y = b + SECOND_RHO;
                let mut value = numerator.get(x, y);
                // Divide by the leading term x^11 y^10 in descending lex order.
                // The signed Weyl symmetries give weights at negative coordinates
                // from the same quadrant. Every other term is already solved.
                for &(u, v, denominator) in &self.terms {
                    let i = (x as i128 - u) as usize;
                    let j = (y as i128 - v).unsigned_abs() as usize;
                    if i > last || j > last {
                        continue;
                    }
                    debug_assert!(i > a || (i == a && j > b));
                    let weight = weights.get(i, j);
                    if weight.is_zero() {
                        continue;
                    }
                    model.spend(WEYL_LIMBS)?;
                    let term = denominator
                        .multiply_wide(Signed {
                            magnitude: weight,
                            negative: true,
                        })
                        .map_err(|_| arithmetic())?;
                    value = value.add(term).map_err(|_| arithmetic())?;
                }
                if value.negative && !value.is_zero() {
                    return Err(Error::plain(
                        "rank-two LLV specialization produced a negative weight multiplicity",
                    ));
                }
                *weights.get_mut(a, b) = value
                    .magnitude
                    .divide_exact(model, self.leading)?
                    .ok_or_else(|| {
                        Error::plain("rank-two LLV Weyl specialization did not divide exactly")
                    })?;
            }
        }
        verify_character(model, shape, &weights)?;
        Ok(weights)
    }
}

fn alternant(model: &mut Model, shape: &[u128], last: usize) -> Result<Square<Signed>> {
    let highest: [usize; RANK] =
        std::array::from_fn(|i| shape.get(i).copied().unwrap_or(0) as usize + LARGEST_RHO - i);
    let doubled: [u128; RANK] = std::array::from_fn(|i| (2 * highest[i] + 1) as u128);
    let mut coefficients = Square::new(last + RANK, Signed::from_i128(0))?;
    // The final ten confluent rows are odd powers of the shifted weights.
    // Expand along the first two rows; their entries are interval characters.
    for (i, &a) in highest.iter().enumerate() {
        for (j, &b) in highest.iter().enumerate().skip(i + 1) {
            let magnitude = weyl::odd_vandermonde(model, &doubled, &[i, j])?;
            let term = Signed {
                magnitude,
                negative: (i + j).is_multiple_of(2),
            };
            add(coefficients.get_mut(a, b), term)?;
            add(coefficients.get_mut(b, a), term.negate())?;
        }
    }
    // Each cofactor at a corner fills the rectangle below its two bounds.
    model.spend(2 * WEYL_LIMBS * coefficients.cells.len())?;
    for a in (0..coefficients.width).rev() {
        for b in (0..coefficients.width).rev() {
            if a + 1 < coefficients.width {
                let up = coefficients.get(a + 1, b);
                add(coefficients.get_mut(a, b), up)?;
            }
        }
    }
    for a in 0..coefficients.width {
        for b in (0..coefficients.width - 1).rev() {
            let right = coefficients.get(a, b + 1);
            add(coefficients.get_mut(a, b), right)?;
        }
    }
    Ok(coefficients)
}

fn add(target: &mut Signed, value: Signed) -> Result<()> {
    *target = target.add(value).map_err(|_| arithmetic())?;
    Ok(())
}

fn verify_character(model: &mut Model, shape: &[u128], weights: &Square<Magnitude>) -> Result<()> {
    let dimension = weyl::dimension(model, shape)?;
    let mut total = Magnitude::default();
    model.spend(weights.cells.len())?;
    for a in 0..weights.width {
        for b in 0..weights.width {
            let degree = weights.get(a, b);
            if degree != weights.get(b, a) {
                return Err(Error::plain(
                    "rank-two LLV character is not symmetric in the grading torus",
                ));
            }
            let orbit = (1 + u32::from(a != 0)) * (1 + u32::from(b != 0));
            total.add_multiple(degree, orbit).ok_or_else(arithmetic)?;
        }
    }
    if total != dimension {
        return Err(Error::plain(
            "rank-two LLV character does not sum to the Weyl dimension",
        ));
    }
    Ok(())
}
