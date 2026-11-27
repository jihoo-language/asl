use std::collections::BTreeSet;

use crate::complex::{self, OrderComplex, Simplex};
use crate::error::{Error, Result};
use crate::model::{Model, Poset};
use crate::odd_cobar::Prime;
use crate::prime_linear::{Basis, Coordinate, QuotientBasis};
use crate::sheaf;
use crate::sheaf_base_change::Square;

type Chain = Vec<(u32, Vec<usize>)>;

struct Group {
    degree: i64,
    rank: usize,
    cycle: Chain,
}

#[derive(Clone, Copy)]
enum Request {
    First,
    Ranks,
    Cycle(i64, usize),
}

pub(crate) fn universal(model: &mut Model, square: &Square, prime: Prime) -> Result<String> {
    let mut checked = BTreeSet::new();
    let mut cells = 0;
    for (y, yname) in square.yp.names.iter().enumerate() {
        model.spend(1 + square.x.names.len())?;
        for (x, xname) in square.x.names.iter().enumerate() {
            if !square.in_future(y, x) {
                continue;
            }
            let points = square.comma(model, y, x)?;
            model.spend(1 + (points.len() + 1) * checked.len().max(1).ilog2() as usize)?;
            if checked.contains(&points) {
                continue;
            }
            let found = inspect(model, &square.xp, &points, prime, Request::First)?;
            if let Some(failure) = found.first() {
                let cycle = named_chain(&square.xp, &failure.cycle);
                return Ok(format!(
                    "shbcuniversal = fails at ({yname}, {xname}) in degree {} with cycle {cycle:?}",
                    failure.degree,
                ));
            }
            cells += 4 + points.len();
            sheaf::check_cells(cells)?;
            checked.insert(points);
        }
    }
    Ok("shbcuniversal = exact".to_owned())
}

pub(crate) fn ranks(
    model: &mut Model,
    square: &Square,
    y: usize,
    x: usize,
    prime: Prime,
) -> Result<String> {
    let points = square.comma(model, y, x)?;
    let out = inspect(model, &square.xp, &points, prime, Request::Ranks)?;
    let grades: Vec<_> = out.into_iter().map(|v| (v.degree, v.rank)).collect();
    Ok(format!("shbccomma = {grades:?}"))
}

pub(crate) fn witness(
    model: &mut Model,
    square: &Square,
    y: usize,
    x: usize,
    prime: Prime,
    degree: i64,
    index: usize,
) -> Result<String> {
    if degree < -1 {
        return Err(Error::plain("augmented comma homology starts in degree -1"));
    }
    let points = square.comma(model, y, x)?;
    let mut out = inspect(
        model,
        &square.xp,
        &points,
        prime,
        Request::Cycle(degree, index),
    )?;
    let actual = out
        .pop()
        .ok_or_else(|| Error::plain(format!("no comma homology generator in degree {degree}")))?;
    Ok(format!(
        "shbcwitness = {:?}",
        named_chain(&square.xp, &actual.cycle)
    ))
}

fn named_chain<'a>(site: &'a Poset, cycle: &Chain) -> Vec<(u32, Vec<&'a str>)> {
    cycle
        .iter()
        .map(|(c, chain)| {
            (
                *c,
                chain
                    .iter()
                    .map(|&point| site.names[point].as_str())
                    .collect(),
            )
        })
        .collect()
}

fn inspect(
    model: &mut Model,
    site: &Poset,
    points: &[usize],
    prime: Prime,
    request: Request,
) -> Result<Vec<Group>> {
    if points.is_empty() {
        let ranks = matches!(
            request,
            Request::Ranks | Request::First | Request::Cycle(-1, 0)
        );
        return Ok(if ranks {
            vec![Group {
                degree: -1,
                rank: 1,
                cycle: vec![(1, Vec::new())],
            }]
        } else {
            Vec::new()
        });
    }
    if matches!(request, Request::Cycle(-1, _)) || cone(model, site, points)? {
        return Ok(Vec::new());
    }
    let mut order = OrderComplex::induced(model, site, points)?;
    let mut previous = vec![Vec::new()];
    let mut current = order.vertices();
    let mut outgoing_rank = 1;
    let mut degree = 0i64;
    let mut output = Vec::new();
    while !current.is_empty() {
        if matches!(request,Request::Cycle(stop,_) if degree > stop) {
            break;
        }
        let next = order.extend(model, &current)?;
        let retained = layer_cells(&previous) + layer_cells(&current) + layer_cells(&next);
        sheaf::check_cells(retained)?;
        let mut incoming = Basis::default();
        boundaries(model, prime, &current, &next, |model, _, column| {
            incoming.insert(model, prime, column)?;
            sheaf::check_cells(retained + incoming.cells())
        })?;
        let rank = current
            .len()
            .checked_sub(outgoing_rank + incoming.rank())
            .ok_or_else(|| Error::plain("invalid augmented comma boundary ranks"))?;
        let requested_index = match request {
            Request::First if rank != 0 => Some(0),
            Request::Cycle(stop, i) if stop == degree && i < rank => Some(i),
            _ => None,
        };
        if let Some(index) = requested_index {
            let raw = cycle(
                model,
                prime,
                incoming,
                (&previous, &current),
                retained,
                index,
            )?;
            let chain = raw
                .into_iter()
                .map(|(column, c)| {
                    (
                        c,
                        current[column].iter().map(|&local| points[local]).collect(),
                    )
                })
                .collect();
            return Ok(vec![Group {
                degree,
                rank,
                cycle: chain,
            }]);
        }
        if rank != 0 && matches!(request, Request::Ranks) {
            output.push(Group {
                degree,
                rank,
                cycle: Vec::new(),
            });
        }
        outgoing_rank = incoming.rank();
        previous = current;
        current = next;
        degree += 1;
    }
    Ok(output)
}

fn cone(model: &mut Model, site: &Poset, points: &[usize]) -> Result<bool> {
    model.spend(2 * points.len() * points.len())?;
    Ok(points.iter().any(|&a| {
        points.iter().all(|&b| a == b || site.strict_less(a, b))
            || points.iter().all(|&b| a == b || site.strict_less(b, a))
    }))
}

fn layer_cells(layer: &[Simplex]) -> usize {
    layer.len() * (3 + layer.first().map_or(0, Vec::len))
}

fn boundaries(
    model: &mut Model,
    prime: Prime,
    source_faces: &[Simplex],
    chains: &[Simplex],
    mut visit: impl FnMut(&mut Model, usize, Vec<Coordinate>) -> Result<()>,
) -> Result<()> {
    complex::for_each_boundary(model, source_faces, chains, |model, column, faces| {
        let mut signed: Vec<_> = faces
            .into_iter()
            .enumerate()
            .map(|(face, row)| (row, if face % 2 == 0 { 1 } else { prime.negate(1) }))
            .collect();
        model.spend(signed.len() * (1 + signed.len().max(1).ilog2() as usize))?;
        signed.sort_unstable_by_key(|&(row, _)| row);
        visit(model, column, signed)
    })
}

fn cycle(
    model: &mut Model,
    prime: Prime,
    incoming: Basis,
    (previous, current): (&[Simplex], &[Simplex]),
    retained: usize,
    index: usize,
) -> Result<Vec<Coordinate>> {
    let mut outgoing = Vec::with_capacity(current.len());
    let mut cells = retained + incoming.cells();
    boundaries(model, prime, previous, current, |_, _, column| {
        cells += 3 + 2 * column.len();
        sheaf::check_cells(cells)?;
        outgoing.push(column);
        Ok(())
    })?;
    let group = QuotientBasis::build(model, prime, incoming, current.len(), |model, i| {
        model.spend(outgoing[i].len())?;
        Ok(outgoing[i].clone())
    })?;
    group.representative(model, prime, &[(index, 1)])
}
