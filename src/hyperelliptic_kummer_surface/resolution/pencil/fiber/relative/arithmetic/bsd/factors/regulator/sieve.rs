use std::collections::BTreeSet;
use std::rc::Rc;

use crate::ec_morphism::Rational;
use crate::ec_poly::{self as poly, Field, Scalar};
use crate::elliptic::{Curve, Point};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;

use super::super::super::super::{Kodaira, RelativeJacobian};
use super::{field, Pencil, Section};

const MAX_RESIDUE_ORDER: u64 = 256;
const MAX_LOCAL_GROUP: usize = 2 * (MAX_RESIDUE_ORDER as usize + 1);
const MAX_COMBINATIONS: usize = 4096;
const MAX_FIBERS_WITH_RELEVANT_TORSION: usize = 16;

pub(super) struct Evidence {
    pub(super) excluded: bool,
    pub(super) report: String,
}

pub(super) fn test(
    model: &mut Model,
    pencil: &Pencil,
    relative: &RelativeJacobian,
    generators: &[Rc<Section>],
    prime: u128,
) -> Result<Evidence> {
    let rank = generators.len();
    let field = field(&pencil.curve);
    let q = field.ring.size_of_residue_field()?;
    let ell = usize::try_from(prime).ok();
    let count = ell.and_then(|ell| ell.checked_pow(rank as u32));
    let (Some(ell), Some(size)) = (ell, count.filter(|&s| s <= MAX_COMBINATIONS)) else {
        return Ok(Evidence {
            excluded: false,
            report: format!("ell={prime} inconclusive (combination limit)"),
        });
    };
    if q > MAX_RESIDUE_ORDER {
        return Ok(Evidence {
            excluded: false,
            report: format!("ell={prime} inconclusive (residue field limit)"),
        });
    }
    let mut remaining: Vec<_> = (1..size).collect();
    let mut useful = Vec::new();
    let mut tested = 0;
    for i in 0..q as usize {
        if remaining.is_empty() || tested == MAX_FIBERS_WITH_RELEVANT_TORSION {
            break;
        }
        let t = field.ring.coordinate_row(i);
        let delta = poly::evaluate(model, field, &relative.0.discriminant, &t)?;
        if field.is_zero(&delta) {
            continue;
        }
        let uniformizer = [field.ring.negate(model, &t)?, field.one()];
        let local = relative.minimal_at(model, field, Some(&uniformizer))?;
        if !matches!(local.kodaira, Kodaira::Good) || local.scales != 0 {
            return Err(Error::plain(
                "K3 pencil saturation sample failed the good reduction check",
            ));
        }
        let (a, b) =
            relative.short_coefficients(model, field, Some(&uniformizer), &local, false)?;
        let curve = Curve::new(
            model,
            pencil.curve.grade.ring.clone(),
            pencil.curve.ring.clone(),
            a,
            b,
        )?;
        let points = curve.raw_points(model, MAX_LOCAL_GROUP)?;
        if !points.len().is_multiple_of(ell) {
            continue;
        }
        tested += 1;
        let multiples = points
            .iter()
            .map(|point| curve.multiple(model, point, ell as i128))
            .collect::<Result<BTreeSet<_>>>()?;
        let a2 = poly::evaluate(model, field, &relative.0.a2, &t)?;
        let third = field.inverse_unit(model, &field.signed(3))?;
        let shift = field.multiply(model, &a2, &third)?;
        let shifts = generators
            .iter()
            .map(|section| specialize(model, &curve, &t, &shift, section))
            .collect::<Result<Vec<_>>>()?;
        let mut tables = Vec::with_capacity(rank);
        for point in shifts {
            let mut row = Vec::with_capacity(ell);
            row.push(Point::Infinity);
            for j in 1..ell {
                row.push(curve.add(model, &row[j - 1], &point)?);
            }
            tables.push(row);
        }
        let mut kernel = Vec::new();
        for &candidate in &remaining {
            model.spend(rank + 1)?;
            let mut digits = candidate;
            let mut sum = Point::Infinity;
            for table in &tables {
                let digit = digits % ell;
                if digit != 0 {
                    sum = curve.add(model, &sum, &table[digit])?;
                }
                digits /= ell;
            }
            if multiples.contains(&sum) {
                kernel.push(candidate);
            }
        }
        if kernel.len() < remaining.len() {
            let dimension = kernel_dimension(kernel.len() + 1, ell)?;
            useful.push((t, points.len(), dimension));
            remaining = kernel;
        }
    }
    let excluded = remaining.is_empty();
    let conclusion = if excluded {
        "excluded".into()
    } else {
        let dim = kernel_dimension(remaining.len() + 1, ell)?;
        let limit = if tested == MAX_FIBERS_WITH_RELEVANT_TORSION {
            ", fiber limit reached"
        } else {
            ""
        };
        format!("inconclusive (kernel dimension {dim}{limit})")
    };
    Ok(Evidence {
        excluded,
        report: format!("ell={prime}, reductions (parameter, group order, common kernel dimension)={useful:?}, {conclusion}"),
    })
}

fn specialize(
    model: &mut Model,
    curve: &Curve,
    parameter: &Scalar,
    shift: &Scalar,
    point: &Section,
) -> Result<Point> {
    let Section::Affine(x, y) = point else {
        return Ok(Point::Infinity);
    };
    let field = curve.field();
    let Some(x) = evaluate(model, field, x, parameter)? else {
        // A pole of x at a good minimal Weierstrass place meets O.
        return Ok(Point::Infinity);
    };
    let y = evaluate(model, field, y, parameter)?
        .ok_or_else(|| Error::plain("K3 pencil good reduction has finite x with a pole of y"))?;
    let short_x = field.add(model, &x, shift)?;
    curve.point(model, short_x, y)
}

fn evaluate(
    model: &mut Model,
    field: Field<'_>,
    value: &Rational,
    t: &Scalar,
) -> Result<Option<Scalar>> {
    let denominator = poly::evaluate(model, field, &value.denominator, t)?;
    if field.is_zero(&denominator) {
        return Ok(None);
    }
    let numerator = poly::evaluate(model, field, &value.numerator, t)?;
    let inverse = field.inverse_unit(model, &denominator)?;
    field.multiply(model, &numerator, &inverse).map(Some)
}

fn kernel_dimension(mut order: usize, prime: usize) -> Result<usize> {
    let mut dimension = 0;
    while order > 1 && order.is_multiple_of(prime) {
        dimension += 1;
        order /= prime;
    }
    if order != 1 {
        return Err(Error::plain(
            "K3 pencil specialization kernel has invalid cardinality",
        ));
    }
    Ok(dimension)
}
