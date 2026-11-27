use crate::rv_cyclotomic as cyclo;

use super::field::Field;
use super::jordan::Summand;
use super::{check_report, Defect, Error, Model, Result};

fn interval(defect: &Defect) -> Result<(usize, usize)> {
    let d = defect.d;
    let start = defect
        .indices
        .iter()
        .copied()
        .find(|&r| defect.indices.binary_search(&((r + d - 1) % d)).is_err())
        .unwrap();
    for i in 0..defect.indices.len() {
        if defect.indices.binary_search(&((start + i) % d)).is_err() {
            return Err(Error::plain(
                "LG defect fusion rule requires consecutive circular intervals",
            ));
        }
    }
    Ok((start, defect.indices.len() - 1))
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
    let (a, l) = interval(left)?;
    let (b, m) = interval(right)?;
    let d = left.d;
    let low = l.abs_diff(m);
    let high = (l + m).min(2 * (d - 2) - l - m);
    let field = Field::new(model, d)?;
    let left_action = action(model, &field, left)?;
    let right_action = action(model, &field, right)?;
    let mut channels = Vec::new();
    let mut combined = vec![Vec::new(); d];
    let mut bulk = vec![Vec::new(); d - 1];
    let mut reverse = vec![Vec::new(); d - 1];
    for n in (low..=high).step_by(2) {
        let start = (a + b + (l + m - n) / 2) % d;
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
        channels.push((start, n));
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
        "channels={channels:?} verified_brane_actions={} verified_bulk_degrees={} orientations=2",
        d - 1,
        d - 1
    );
    check_report(&report)?;
    Ok(report)
}
