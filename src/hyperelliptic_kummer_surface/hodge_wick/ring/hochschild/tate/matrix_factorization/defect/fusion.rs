use crate::rv_cyclotomic as cyclo;

use super::field::Field;
use super::jordan::Summand;
use super::{check_report, Defect, Error, Model, Result};

pub(super) fn interval(defect: &Defect) -> Option<(usize, usize)> {
    let d = defect.d;
    let start = defect
        .indices
        .iter()
        .copied()
        .find(|&r| defect.indices.binary_search(&((r + d - 1) % d)).is_err())?;
    (0..defect.indices.len())
        .all(|i| defect.indices.binary_search(&((start + i) % d)).is_ok())
        .then_some((start, defect.indices.len() - 1))
}

pub(super) fn channels(
    d: usize,
    (a, l): (usize, usize),
    (b, m): (usize, usize),
) -> impl Iterator<Item = (usize, usize)> {
    let low = l.abs_diff(m);
    let high = (l + m).min(2 * (d - 2) - l - m);
    (low..=high)
        .step_by(2)
        .map(move |n| ((a + b + (l + m - n) / 2) % d, n))
}

fn action(model: &mut Model, field: &Field, defect: &Defect) -> Result<Vec<Vec<usize>>> {
    let mut rows = vec![Vec::new(); defect.d];
    for (j, row) in rows.iter_mut().enumerate().skip(1) {
        let summand = Summand::new(model, field, defect, j)?;
        row.extend(summand.lengths().filter(|&i| i < defect.d));
    }
    Ok(rows)
}

pub(crate) fn verify_fusion(model: &mut Model, left: &Defect, right: &Defect) -> Result<String> {
    if left.mark != right.mark || left.d != right.d {
        return Err(Error::plain(
            "LG defect fusion requires the same Tate marking",
        ));
    }
    let (Some(a), Some(b)) = (interval(left), interval(right)) else {
        return Err(Error::plain(
            "LG defect fusion rule requires consecutive circular intervals",
        ));
    };
    let d = left.d;
    let field = Field::new(model, d)?;
    let left_action = action(model, &field, left)?;
    let right_action = action(model, &field, right)?;
    let mut result = Vec::new();
    let mut combined = vec![Vec::new(); d];
    let mut bulk = vec![Vec::new(); d - 1];
    let mut reverse = vec![Vec::new(); d - 1];
    for (start, n) in channels(d, a, b) {
        let support = (0..=n).map(|i| (start + i) % d).collect();
        let channel = Defect::new(model, left.mark.clone(), d, support)?;
        for (all, parts) in combined.iter_mut().zip(action(model, &field, &channel)?) {
            all.extend(parts);
        }
        for (k, (forward, backward)) in bulk.iter_mut().zip(&mut reverse).enumerate() {
            let first = channel.transmission(model, k, false)?;
            let second = channel.transmission(model, k, true)?;
            cyclo::accumulate(model, forward, &first, false)?;
            cyclo::accumulate(model, backward, &second, false)?;
        }
        result.push((start, n));
    }
    for (right, predicted) in right_action.iter().zip(&mut combined).skip(1) {
        let mut composed: Vec<_> = right
            .iter()
            .flat_map(|&b| left_action[b].iter().copied())
            .collect();
        composed.sort_unstable();
        predicted.sort_unstable();
        model.spend(composed.len() + predicted.len())?;
        if &composed != predicted {
            return Err(Error::plain(
                "LG defect fusion failed the computed stable boundary action",
            ));
        }
    }
    for (orientation, rows) in [(false, bulk), (true, reverse)] {
        for (k, expected) in rows.iter().enumerate() {
            let a = left.transmission(model, k, orientation)?;
            let b = right.transmission(model, k, orientation)?;
            if &cyclo::product_coefficients(model, d, &a, &b)? != expected {
                return Err(Error::plain(
                    "LG defect fusion failed the computed bulk residue action",
                ));
            }
        }
    }
    let report = format!(
        "channels={result:?} verified_brane_actions={} verified_bulk_degrees={} orientations=2",
        d - 1,
        d - 1
    );
    check_report(&report)?;
    Ok(report)
}
