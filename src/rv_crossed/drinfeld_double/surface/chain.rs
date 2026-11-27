use std::collections::BTreeMap;

use crate::error::{Error, Result};

use super::{Double, Marking};

pub(super) type Chain<const N: usize> = BTreeMap<[usize; N], i128>;

#[derive(Clone, Copy)]
pub(crate) enum Move {
    T(usize, bool),
    S(usize, bool),
    Slide([usize; 2], bool),
    Swap([usize; 2]),
}

impl Move {
    pub(crate) fn check(self, genus: usize) -> Result<()> {
        let ([a, b], two) = match self {
            Self::T(i, _) | Self::S(i, _) => ([i, i], false),
            Self::Slide([i, j], _) | Self::Swap([i, j]) => ([i, j], true),
        };
        if a >= genus || b >= genus || (two && a == b) {
            return Err(Error::plain(
                "surface marking move requires valid distinct handle indices",
            ));
        }
        Ok(())
    }

    pub(super) fn apply(
        self,
        double: &Double,
        marks: &[usize],
    ) -> Result<(Marking, Chain<3>, usize)> {
        self.check(marks.len() / 2)?;
        let times = |a, b| double.group.mul(a, b);
        let inv = |a| double.group.inv(a);
        let back = matches!(
            self,
            Self::T(_, true) | Self::S(_, true) | Self::Slide(_, true)
        );
        let mut output = marks.to_vec();
        match self {
            Self::T(i, back) => {
                let a = marks[2 * i];
                output[2 * i + 1] = times(marks[2 * i + 1], if back { inv(a) } else { a });
            }
            Self::S(i, back) => {
                let (a, b) = (marks[2 * i], marks[2 * i + 1]);
                (output[2 * i], output[2 * i + 1]) = if back { (inv(b), a) } else { (b, inv(a)) };
            }
            Self::Slide([i, j], back) => {
                let (a, b, c, d) = (
                    marks[2 * i],
                    marks[2 * i + 1],
                    marks[2 * j],
                    marks[2 * j + 1],
                );
                (output[2 * i], output[2 * j + 1]) = if back {
                    (times(a, inv(c)), times(d, b))
                } else {
                    (times(a, c), times(d, inv(b)))
                };
            }
            Self::Swap([i, j]) => {
                output.swap(2 * i, 2 * j);
                output.swap(2 * i + 1, 2 * j + 1);
            }
        }
        let forward = if back { &output } else { marks };
        let sign = if back { -1 } else { 1 };
        let mut chain = Chain::new();
        match self {
            Self::T(i, _) => {
                let (a, b) = (forward[2 * i], forward[2 * i + 1]);
                accumulate(&mut chain, [a, b, a], sign);
            }
            Self::S(i, _) => {
                let (a, b) = (forward[2 * i], forward[2 * i + 1]);
                shuffle(&mut chain, [inv(a), a, b], -sign);
            }
            Self::Slide([i, j], _) => {
                let (a, b, c, d) = (
                    forward[2 * i],
                    forward[2 * i + 1],
                    forward[2 * j],
                    forward[2 * j + 1],
                );
                shuffle(&mut chain, [a, c, b], -sign);
                shuffle(&mut chain, [times(d, inv(b)), b, c], -sign);
            }
            Self::Swap(_) => {}
        }
        let mut expected = cycle(&output);
        for (p, s) in cycle(marks) {
            accumulate(&mut expected, p, -s);
        }
        let mut boundary = Chain::new();
        for (&[a, b, c], &s) in &chain {
            accumulate(&mut boundary, [b, c], s);
            accumulate(&mut boundary, [times(a, b), c], -s);
            accumulate(&mut boundary, [a, times(b, c)], s);
            accumulate(&mut boundary, [a, b], -s);
        }
        if boundary != expected {
            return Err(Error::plain(
                "marking move bar chain does not have the claimed integral boundary",
            ));
        }
        let sum = chain
            .iter()
            .map(|(&p, &s)| s * double.cocycle_phase(p) as i128)
            .sum::<i128>();
        let phase = sum.rem_euclid(double.conductor as i128) as usize;
        Ok((output, chain, phase))
    }
}

pub(super) fn cycle(marks: &[usize]) -> Chain<2> {
    let mut result = Chain::new();
    for pair in marks.chunks_exact(2) {
        accumulate(&mut result, [pair[0], pair[1]], 1);
        accumulate(&mut result, [pair[1], pair[0]], -1);
    }
    result
}

fn shuffle(chain: &mut Chain<3>, [a, b, c]: [usize; 3], s: i128) {
    accumulate(chain, [a, b, c], s);
    accumulate(chain, [a, c, b], -s);
    accumulate(chain, [c, a, b], s);
}

fn accumulate<const N: usize>(chain: &mut Chain<N>, term: [usize; N], delta: i128) {
    if term.contains(&0) || delta == 0 {
        return;
    }
    let n = chain.entry(term).or_default();
    *n += delta;
    if *n == 0 {
        chain.remove(&term);
    }
}
