use std::collections::BTreeMap;

use super::{
    append_crossing, arithmetic, projected_quantum_trace, signed_word, storage, Block, Context,
    Generator, Grade, MAX_WORD,
};
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_lens::{self, Rational};

struct Link {
    braid: Vec<Generator>,
    component_of: Vec<usize>,
    self_writhe: Vec<i128>,
    links: BTreeMap<(usize, usize), i128>,
}

impl Link {
    fn new(model: &mut Model, strands: usize, word: &[i128]) -> Result<Self> {
        let braid = signed_word(model, strands, word)?;
        storage(strands.saturating_mul(6).saturating_add(5 * braid.len()))?;
        model.spend(5 * strands + 2 * braid.len())?;
        let mut permutation: Vec<_> = (0..strands).collect();
        for &letter in braid.iter().rev() {
            let (index, _) = crossing(letter);
            permutation.swap(index - 1, index);
        }

        // The first unvisited original strand is the least index of its
        // permutation cycle and therefore determines the component order.
        let mut component_of = vec![usize::MAX; strands];
        let mut self_writhe = Vec::new();
        for first in 0..strands {
            if component_of[first] != usize::MAX {
                continue;
            }
            let component = self_writhe.len();
            self_writhe.push(0);
            let mut current = first;
            while component_of[current] == usize::MAX {
                component_of[current] = component;
                current = permutation[current];
            }
        }

        for (original, label) in permutation.iter_mut().enumerate() {
            *label = component_of[original];
        }
        let mut links = BTreeMap::<_, i128>::new();
        for &letter in braid.iter().rev() {
            let (index, positive) = crossing(letter);
            let (left, right) = (permutation[index - 1], permutation[index]);
            let sign = if positive { 1 } else { -1 };
            if left == right {
                self_writhe[left] += sign;
            } else {
                model.spend(1 + links.len().max(1).ilog2() as usize)?;
                *links.entry((left.min(right), left.max(right))).or_default() += sign;
            }
            permutation.swap(index - 1, index);
        }
        links.retain(|_, count| *count != 0);
        for count in links.values_mut() {
            debug_assert_eq!(count.rem_euclid(2), 0);
            *count /= 2;
        }

        Ok(Self {
            braid,
            component_of,
            self_writhe,
            links,
        })
    }

    fn components(&self) -> usize {
        self.self_writhe.len()
    }

    fn check_framings(&self, framings: &[i128]) -> Result<()> {
        if framings.len() != self.components() {
            return Err(Error::plain(format!(
                "framed braid closure has {} components and requires one framing per component",
                self.components()
            )));
        }
        Ok(())
    }

    fn matrix(&self, model: &mut Model, framings: &[i128]) -> Result<Vec<Vec<i128>>> {
        let components = self.components();
        let cells = components.saturating_mul(components + 1);
        storage(cells)?;
        model.spend(cells + self.links.len())?;
        let mut matrix = vec![vec![0; components]; components];
        for (i, &framing) in framings.iter().enumerate() {
            matrix[i][i] = framing;
        }
        for (&(i, j), &count) in &self.links {
            matrix[i][j] = count;
            matrix[j][i] = count;
        }
        Ok(matrix)
    }
}

fn crossing(letter: Generator) -> (usize, bool) {
    match letter {
        Generator::Braid(index, positive) => (index, positive),
        Generator::Tl(_) => unreachable!(),
    }
}

struct Cable {
    strands: usize,
    blocks: Vec<Block>,
    word: Vec<Generator>,
}

impl Cable {
    fn new(model: &mut Model, link: &Link, colors: &[usize]) -> Result<Self> {
        model.spend(2 * link.component_of.len())?;
        let mut strands = 0;
        let mut blocks = Vec::with_capacity(link.component_of.len());
        for &component in &link.component_of {
            let length = colors[component];
            blocks.push(Block {
                start: strands,
                length,
            });
            strands = strands.checked_add(length).ok_or_else(arithmetic)?;
        }
        let mut moving = blocks.clone();
        let mut word = Vec::new();
        for &letter in link.braid.iter().rev() {
            let (index, positive) = crossing(letter);
            let Block {
                start,
                length: left,
            } = moving[index - 1];
            let right = moving[index].length;
            let crossings = left.checked_mul(right).ok_or_else(arithmetic)?;
            if crossings > MAX_WORD - word.len() {
                return Err(Error::plain(
                    "cabled braid exceeds the limit of 4096 elementary crossings",
                ));
            }
            model.spend(1 + crossings)?;
            append_crossing(&mut word, start, left, right, positive);
            moving[index - 1].length = right;
            moving[index] = Block {
                start: start + right,
                length: left,
            };
        }
        word.reverse();
        debug_assert!(moving
            .iter()
            .zip(&blocks)
            .all(|(a, b)| a.start == b.start && a.length == b.length));
        Ok(Self {
            strands,
            blocks,
            word,
        })
    }
}

fn twist_exponent(grade: Grade, color: usize, power: i128) -> i128 {
    let order = grade.order as i128;
    let color = color as i128;
    let twist = (color * (color + 2) + color.rem_euclid(2) * (order / 2)).rem_euclid(order);
    (twist * power.rem_euclid(order)).rem_euclid(order)
}

fn trace(
    model: &mut Model,
    grade: Grade,
    link: &Link,
    colors: &[usize],
    framings: &[i128],
    context: &mut Context,
) -> Result<Rational> {
    let cable = Cable::new(model, link, colors)?;
    let sum = projected_quantum_trace(
        model,
        Grade {
            strands: cable.strands,
            ..grade
        },
        cable.blocks.iter().copied(),
        &cable.word,
        context,
    )?;
    let order = grade.order as i128;
    let exponent = colors
        .iter()
        .zip(framings.iter().zip(&link.self_writhe))
        .fold(0, |sum, (&color, (&framing, &self_writhe))| {
            let power = framing.rem_euclid(order) - self_writhe.rem_euclid(order);
            (sum + twist_exponent(grade, color, power)).rem_euclid(order)
        });
    rv_lens::shifted(model, grade.order, &sum, exponent)
}

pub(crate) fn framed(
    model: &mut Model,
    grade: Grade,
    word: &[i128],
    colors: &[usize],
    framings: &[i128],
) -> Result<Value> {
    let link = Link::new(model, grade.strands, word)?;
    link.check_framings(framings)?;
    if colors.len() != link.components() {
        return Err(Error::plain(format!(
            "framed braid closure has {} components and requires one color per component",
            link.components()
        )));
    }
    if colors.iter().any(|&color| color > grade.level) {
        return Err(Error::plain(
            "framed braid colors must be at most the level",
        ));
    }
    let value = trace(
        model,
        grade,
        &link,
        colors,
        framings,
        &mut Context::default(),
    )?;
    rv_lens::rational_value(model, value)
}

pub(crate) fn linking(
    model: &mut Model,
    strands: usize,
    word: &[i128],
    framings: &[i128],
) -> Result<Vec<Vec<i128>>> {
    let link = Link::new(model, strands, word)?;
    link.check_framings(framings)?;
    link.matrix(model, framings)
}
