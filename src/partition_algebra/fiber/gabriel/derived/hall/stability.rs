use std::cmp::Ordering;
use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::Model;
use crate::rv_qp_quad::Wide;

use super::super::auslander::tau::Fan;
use super::super::{Nakayama, Parameter, Path};
use super::{arithmetic, require_hall_length, Hall, Record, Shape, StabilityQuery};

#[derive(Clone, Copy, Eq, PartialEq)]
struct Slope(Rational);

type Factors = Vec<(Slope, Shape)>;

impl Ord for Slope {
    fn cmp(&self, other: &Self) -> Ordering {
        let a = Wide::<8>::product(self.0.numerator, other.0.denominator);
        let b = Wide::<8>::product(other.0.numerator, self.0.denominator);
        match (a.negative, b.negative) {
            (true, false) => Ordering::Less,
            (false, true) => Ordering::Greater,
            (true, true) => b.magnitude.cmp(a.magnitude),
            (false, false) => a.magnitude.cmp(b.magnitude),
        }
    }
}
impl PartialOrd for Slope {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Hall {
    pub(crate) fn stability(
        &self,
        model: &mut Model,
        query: StabilityQuery,
        theta: &[i128],
        value: &Record,
    ) -> Result<String> {
        if theta.len() != self.vertices {
            return Err(Error::plain(
                "partition Hall stability weight must have one entry per simple vertex",
            ));
        }
        let shape = value.module();
        let factors = self.hn(model, theta, shape)?;
        match query {
            StabilityQuery::Stability => {
                let semistable = factors.len() <= 1;
                let stable = semistable && shape.0.len() == 1 && shape.0[0].2 == 1 && {
                    let (i, l, _) = shape.0[0];
                    is_stable(theta, &self.path_heads[usize::from(i)][..usize::from(l)])?
                };
                let slope = (shape.length() != 0)
                    .then(|| slope_of_shape(self, theta, shape))
                    .transpose()?;
                let at_zero = slope.is_none_or(|s| s.0 == Rational::ZERO);
                let slope = slope.map_or_else(|| "none".into(), |s| s.0.to_string());
                Ok(format!("rvPartitionHallStability = slope={slope}, semistable={semistable}, stable={stable}, kingSemistable={}, kingStable={}",semistable&&at_zero,stable&&at_zero))
            }
            StabilityQuery::Hn => Ok(format!("rvPartitionHallHN = {}", display(&factors))),
            StabilityQuery::HnCheck => {
                require_hall_length(shape.length())?;
                let mut flags = Counter {
                    hall: self,
                    theta,
                    cache: BTreeMap::new(),
                };
                let actual = flags.flags(model, shape, None)?;
                if actual.number != 1 || actual.witness != factors {
                    return Err(Error::plain("partition Hall semistable convolution failed the unique Harder–Narasimhan flag identity"));
                }
                Ok(format!(
                    "rvPartitionHallHNCheck[q={}] = flags=1, factors={}",
                    self.field.order(),
                    display(&factors)
                ))
            }
        }
    }

    fn hn(&self, model: &mut Model, theta: &[i128], module: &Shape) -> Result<Factors> {
        let mut merged = BTreeMap::<Slope, BTreeMap<(u8, u8), u16>>::new();
        for &(i, l, m) in &module.0 {
            let heads = &self.path_heads[usize::from(i)][..usize::from(l)];
            let mut end = heads.len();
            let mut previous = None;
            while end != 0 {
                let mut best = None;
                for k in 0..end {
                    model.spend(end - k)?;
                    let slope = slope_of_layers(theta, &heads[k..end])?;
                    let candidate = (slope, end - k);
                    if best.is_none_or(|pair| pair < candidate) {
                        best = Some(candidate);
                    }
                }
                let (slope, length) = best.unwrap();
                end -= length;
                if previous.is_some_and(|last| last <= slope) {
                    return Err(Error::plain(
                        "partition Harder–Narasimhan slopes did not decrease",
                    ));
                }
                previous = Some(slope);
                let key = (heads[end] as u8, length as u8);
                *merged.entry(slope).or_default().entry(key).or_default() += m;
            }
        }
        Ok(merged
            .into_iter()
            .rev()
            .map(|(slope, parts)| {
                (
                    slope,
                    Shape(parts.into_iter().map(|((i, l), m)| (i, l, m)).collect()),
                )
            })
            .collect())
    }
}

fn is_stable(theta: &[i128], heads: &[usize]) -> Result<bool> {
    let all = slope_of_layers(theta, heads)?;
    for k in 1..heads.len() {
        if slope_of_layers(theta, &heads[k..])? >= all {
            return Ok(false);
        }
    }
    Ok(true)
}

fn slope_of_layers(theta: &[i128], layers: &[usize]) -> Result<Slope> {
    let weight = layers.iter().try_fold(0i128, |w, &i| {
        w.checked_add(theta[i]).ok_or_else(weight_overflow)
    })?;
    Ok(Slope(
        Rational::ratio(weight, layers.len() as i128).map_err(|_| weight_overflow())?,
    ))
}

fn slope_of_shape(hall: &Hall, theta: &[i128], shape: &Shape) -> Result<Slope> {
    let weight = hall
        .dim(shape)
        .into_iter()
        .zip(theta)
        .try_fold(0i128, |w, (d, &t)| {
            (d as i128)
                .checked_mul(t)
                .and_then(|x| w.checked_add(x))
                .ok_or_else(weight_overflow)
        })?;
    Ok(Slope(
        Rational::ratio(weight, shape.length() as i128).map_err(|_| weight_overflow())?,
    ))
}

fn weight_overflow() -> Error {
    Error::plain("partition Hall stability weight exceeds signed 128-bit arithmetic")
}

fn display(factors: &Factors) -> String {
    let values: Vec<_> = factors
        .iter()
        .map(|(s, x)| format!("({},{:?})", s.0, x.0))
        .collect();
    format!("[{}]", values.join(", "))
}

#[derive(Clone)]
struct Flags {
    number: u128,
    witness: Factors,
}

struct Counter<'a> {
    hall: &'a Hall,
    theta: &'a [i128],
    cache: BTreeMap<(Shape, Option<Slope>), Flags>,
}

impl Counter<'_> {
    fn flags(
        &mut self,
        model: &mut Model,
        module: &Shape,
        ceiling: Option<Slope>,
    ) -> Result<Flags> {
        if module.0.is_empty() {
            return Ok(Flags {
                number: 1,
                witness: Vec::new(),
            });
        }
        let key = (module.clone(), ceiling);
        if let Some(result) = self.cache.get(&key) {
            return Ok(result.clone());
        }
        let dims = self.hall.dim(module);
        let total = dims.iter().map(|&d| d + 1).product::<usize>();
        let mut result = Flags {
            number: 0,
            witness: Vec::new(),
        };
        for mut index in 1..total {
            model.spend(dims.len())?;
            let sub = dims
                .iter()
                .map(|&d| {
                    let v = index % (d + 1);
                    index /= d + 1;
                    v
                })
                .collect();
            for ((quotient, submodule), &coefficient) in
                self.hall.distribution(model, module, &sub)?.iter()
            {
                let stability = self.hall.hn(model, self.theta, submodule)?;
                let [(slope, _)] = stability.as_slice() else {
                    continue;
                };
                if ceiling.is_some_and(|c| c <= *slope) {
                    continue;
                }
                let tail = self.flags(model, quotient, Some(*slope))?;
                if tail.number == 0 {
                    continue;
                }
                if result.number == 0 {
                    result.witness.push((*slope, submodule.clone()));
                    result.witness.extend(tail.witness);
                }
                result.number = coefficient
                    .checked_mul(tail.number)
                    .and_then(|x| result.number.checked_add(x))
                    .ok_or_else(arithmetic)?;
            }
        }
        self.cache.insert(key, result.clone());
        Ok(result)
    }
}

pub(crate) fn chamber_report(algebra: &Nakayama, fan: &Fan, index: u128) -> Result<String> {
    let (theta, mut expected) = fan.stability_chamber(index)?;
    let mut actual = Vec::new();
    for state in algebra.resolutions.keys() {
        let mut positive = true;
        for length in 1..=state.length {
            let heads = layer_heads(algebra, [state.head as u128, length as u128]);
            if slope_of_layers(&theta, &heads)?.0.numerator <= 0 {
                positive = false;
                break;
            }
        }
        if positive {
            actual.push([state.head as u128, state.length as u128]);
        }
    }
    expected.sort_unstable();
    if actual != expected {
        return Err(Error::plain(
            "partition Hall chamber positive quotient test disagrees with Fac",
        ));
    }
    let Parameter { rank, n } = algebra.parameter;
    Ok(format!("rvPartitionHallChamber[{rank},{n},{index}] = theta={theta:?}, positiveQuotients={actual:?}"))
}

pub(crate) fn wall_report(algebra: &Nakayama, fan: &Fan, indices: [u128; 2]) -> Result<String> {
    let (theta, brick) = fan.stability_wall(indices)?;
    let heads = layer_heads(algebra, brick);
    if slope_of_layers(&theta, &heads)?.0 != Rational::ZERO || !is_stable(&theta, &heads)? {
        return Err(Error::plain(
            "partition Hall wall has no King-stable brick of the computed type",
        ));
    }
    let Parameter { rank, n } = algebra.parameter;
    let [i, j] = indices;
    Ok(format!(
        "rvPartitionHallWall[{rank},{n},{i},{j}] = theta={theta:?}, kingStableBrick={brick:?}"
    ))
}

fn layer_heads(algebra: &Nakayama, key: [u128; 2]) -> Vec<usize> {
    let [start, limit] = key.map(|x| x as usize);
    (0..limit)
        .map(|length| algebra.end(Path { start, length }))
        .collect()
}
