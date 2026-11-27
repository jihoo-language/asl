use std::collections::BTreeSet;

use crate::complex::{Simplex, MAX_PIVOT_CELLS};
use crate::error::{Error, Result};
use crate::model::Model;

use crate::rv_crossed::Crossed;

type Word = Vec<i32>;

pub(in crate::rv_crossed) struct Quotient {
    products: Vec<Vec<usize>>,
    inverses: Vec<usize>,
}

impl Quotient {
    pub(in crate::rv_crossed) fn new(model: &mut Model, crossed: &Crossed) -> Result<Self> {
        let n = crossed.order();
        model.spend(n * n)?;
        let products: Vec<Vec<_>> = (0..n)
            .map(|a| (0..n).map(|b| crossed.quotient_mul(a, b)).collect())
            .collect();
        let inverses = products
            .iter()
            .map(|row| {
                row.iter()
                    .position(|&x| x == 0)
                    .expect("quotient group inverse")
            })
            .collect();
        Ok(Self { products, inverses })
    }

    pub(in crate::rv_crossed) fn order(&self) -> usize {
        self.inverses.len()
    }

    pub(in crate::rv_crossed) fn mul(&self, a: usize, b: usize) -> usize {
        self.products[a][b]
    }

    fn eval(&self, word: &[i32], values: &[usize]) -> usize {
        word.iter().fold(0, |product, &letter| {
            let input = values[letter.unsigned_abs() as usize];
            self.mul(
                product,
                if letter < 0 {
                    self.inverses[input]
                } else {
                    input
                },
            )
        })
    }
}

pub(super) struct Presentation {
    pub(super) edge_generators: Vec<usize>,
    pub(super) original_generators: usize,
    pub(super) generators: Vec<usize>,
    pub(super) relators: BTreeSet<Word>,
    eliminated: Vec<(usize, Word)>,
}

impl Presentation {
    pub(super) fn new(
        model: &mut Model,
        n: usize,
        tree: &[[usize; 2]],
        edges: &[Simplex],
        faces: &[Simplex],
    ) -> Result<Self> {
        let mut edge_variables = vec![vec![0; n]; n];
        let mut count = 0;
        let edge_generators = edges
            .iter()
            .map(|edge| {
                let &[x, y] = edge.as_slice() else {
                    unreachable!()
                };
                if tree.binary_search(&[x, y]).is_err() {
                    count += 1;
                    edge_variables[x][y] = count;
                }
                edge_variables[x][y]
            })
            .collect();
        model.spend(edges.len() + n * n)?;
        let mut relators = BTreeSet::new();
        for face in faces {
            let &[a, b, c] = face.as_slice() else {
                unreachable!()
            };
            let word = [
                edge_variables[a][b] as i32,
                edge_variables[b][c] as i32,
                -(edge_variables[a][c] as i32),
            ];
            insert_relator(
                model,
                &mut relators,
                freely_reduce(word.into_iter().filter(|&x| x != 0)),
            )?;
        }
        let mut eliminated = Vec::new();
        let mut eliminated_letters = 0;
        let mut active = vec![true; count + 1];
        while let Some((relator, index)) = eliminable(model, &relators)? {
            let letter = relator[index];
            let variable = letter.unsigned_abs() as usize;
            let mut replacement = [&relator[index + 1..], &relator[..index]].concat();
            if letter > 0 {
                replacement = inverse(&replacement).collect();
            }
            active[variable] = false;
            relators.remove(&relator);
            let inverse: Word = inverse(&replacement).collect();
            let mut next = BTreeSet::new();
            let mut letters = eliminated_letters + replacement.len();
            if letters > MAX_PIVOT_CELLS {
                return Err(presentation_limit());
            }
            for word in relators {
                model.spend(word.len() + 1)?;
                let size = word.len();
                let retained = if word
                    .iter()
                    .all(|x| x.unsigned_abs() != letter.unsigned_abs())
                {
                    next.insert(word).then_some(size)
                } else {
                    let mut substituted = Vec::with_capacity(size);
                    for x in word {
                        if x.unsigned_abs() == letter.unsigned_abs() {
                            model.spend(replacement.len())?;
                            let part = if x > 0 { &replacement } else { &inverse };
                            for &item in part {
                                append(&mut substituted, item);
                            }
                        } else {
                            append(&mut substituted, x);
                        }
                        if substituted.len() > MAX_PIVOT_CELLS {
                            return Err(presentation_limit());
                        }
                    }
                    let size = substituted.len();
                    insert_relator(model, &mut next, substituted)?.then_some(size)
                };
                letters += retained.unwrap_or(0);
                if letters > MAX_PIVOT_CELLS {
                    return Err(presentation_limit());
                }
            }
            eliminated_letters += replacement.len();
            eliminated.push((variable, replacement));
            relators = next;
        }
        let generators = (1..=count).filter(|&i| active[i]).collect();
        Ok(Self {
            edge_generators,
            original_generators: count,
            generators,
            relators,
            eliminated,
        })
    }

    pub(super) fn satisfies(
        &self,
        model: &mut Model,
        group: &Quotient,
        values: &[usize],
    ) -> Result<bool> {
        for word in &self.relators {
            model.spend(word.len())?;
            if group.eval(word, values) != 0 {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub(super) fn centralizer(
        &self,
        model: &mut Model,
        group: &Quotient,
        values: &[usize],
    ) -> Result<Option<usize>> {
        use std::cmp::Ordering;
        model.spend(group.order() * (1 + self.generators.len()))?;
        let mut order = 0;
        for g in 0..group.order() {
            let conjugates = self
                .generators
                .iter()
                .map(|&i| group.mul(group.mul(g, values[i]), group.inverses[g]));
            match conjugates.cmp(self.generators.iter().map(|&i| values[i])) {
                Ordering::Less => return Ok(None),
                Ordering::Equal => order += 1,
                Ordering::Greater => {}
            }
        }
        Ok(Some(order))
    }

    pub(super) fn recover(
        &self,
        model: &mut Model,
        group: &Quotient,
        values: &mut [usize],
    ) -> Result<()> {
        // Each eliminated generator is a word only in generators still
        // available when it was removed. Back substitution recovers the edge.
        for (i, word) in self.eliminated.iter().rev() {
            model.spend(1 + word.len())?;
            values[*i] = group.eval(word, values);
        }
        Ok(())
    }
}

fn eliminable(model: &mut Model, relators: &BTreeSet<Word>) -> Result<Option<(Word, usize)>> {
    let mut best = None;
    for word in relators {
        model.spend(1)?;
        if best
            .as_ref()
            .is_some_and(|&(size, _, _, _)| size < word.len())
        {
            continue;
        }
        model.spend(word.len().saturating_mul(word.len()))?;
        for (i, x) in word.iter().enumerate() {
            let variable = x.unsigned_abs();
            let option = (word.len(), variable, word, i);
            if best.is_none_or(|previous| option < previous)
                && !word[..i]
                    .iter()
                    .chain(&word[i + 1..])
                    .any(|a| a.unsigned_abs() == variable)
            {
                best = Some(option);
            }
        }
    }
    Ok(best.map(|(_, _, word, i)| (word.clone(), i)))
}

fn append(word: &mut Word, x: i32) {
    if word.last() == Some(&-x) {
        word.pop();
    } else {
        word.push(x);
    }
}

fn freely_reduce(word: impl IntoIterator<Item = i32>) -> Word {
    let mut reduced = Vec::new();
    for x in word {
        append(&mut reduced, x);
    }
    reduced
}

fn inverse(word: &[i32]) -> impl Iterator<Item = i32> + '_ {
    word.iter().rev().map(|x| -x)
}

fn insert_relator(model: &mut Model, set: &mut BTreeSet<Word>, word: Word) -> Result<bool> {
    let (mut start, mut end) = (0, word.len());
    while end > start + 1 && word[start] == -word[end - 1] {
        start += 1;
        end -= 1;
    }
    let word = &word[start..end];
    if word.is_empty() {
        return Ok(false);
    }
    let n = word.len();
    model.spend(n.saturating_mul(n).saturating_mul(2))?;
    let inverted: Word = inverse(word).collect();
    let canonical = [word, &inverted]
        .into_iter()
        .flat_map(|part| (0..n).map(move |i| part[i..].iter().chain(&part[..i])))
        .min_by(|a, b| a.clone().cmp(b.clone()))
        .expect("nonempty relator")
        .copied()
        .collect();
    Ok(set.insert(canonical))
}

fn presentation_limit() -> Error {
    Error::plain("homotopy Tietze presentation exceeds its word storage limit")
}
