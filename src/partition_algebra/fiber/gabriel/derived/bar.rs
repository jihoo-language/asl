use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::lee_linear::{Echelon, Rational};
use crate::model::Model;

use super::{Class, Grade, Nakayama, Parameter, Path};

type Word = Vec<Path>;
type Cochain = BTreeMap<Word, i128>;
const BAR_DEGREE: u128 = 6;

pub(crate) struct MasseyWitness {
    parameter: Parameter,
    arrows: [u128; 3],
    left: Cochain,
    right: Cochain,
    cycle: Cochain,
    class: Class,
    indeterminacy: usize,
}

impl MasseyWitness {
    pub(crate) fn report(self) -> String {
        let Parameter { rank, n } = self.parameter;
        let [a, b, c] = self.arrows;
        format!("rvPartitionMassey[{rank},{n},{a},{b},{c}] = left={}, right={}, cycle={}, class={}, indeterminacy={}", display(&self.left),display(&self.right),display(&self.cycle),self.class.render(),self.indeterminacy)
    }

    pub(crate) fn class(self) -> Result<Class> {
        if self.indeterminacy != 0 {
            return Err(Error::plain(
                "partition Massey product has nonzero indeterminacy",
            ));
        }
        Ok(self.class)
    }
}

fn display(value: &Cochain) -> String {
    let terms: Vec<_> = value
        .iter()
        .map(|(word, &coefficient)| {
            (
                word.iter().map(|p| (p.start, p.length)).collect::<Vec<_>>(),
                coefficient,
            )
        })
        .collect();
    format!("{terms:?}")
}

impl Nakayama {
    fn words(
        &self,
        model: &mut Model,
        start: usize,
        end: usize,
        length: usize,
    ) -> Result<Vec<Word>> {
        let mut words = Vec::new();
        self.extend_words(model, (start, end), length, &mut Vec::new(), &mut words)?;
        Ok(words)
    }

    fn extend_words(
        &self,
        model: &mut Model,
        endpoints: (usize, usize),
        until: usize,
        word: &mut Word,
        words: &mut Vec<Word>,
    ) -> Result<()> {
        let (start, end) = endpoints;
        if word.len() == until {
            if start == end {
                model.spend(1 + word.len())?;
                words.push(word.clone());
            }
            return Ok(());
        }
        for length in 1..self.paths[start].len() {
            model.spend(1)?;
            let path = Path { start, length };
            word.push(path);
            self.extend_words(model, (self.end(path), end), until, word, words)?;
            word.pop();
        }
        Ok(())
    }

    fn coboundary(
        &self,
        model: &mut Model,
        cochain: &Cochain,
        targets: &[Word],
    ) -> Result<Cochain> {
        let mut result = Cochain::new();
        for word in targets {
            let mut coefficient = 0;
            for k in 0..word.len().saturating_sub(1) {
                model.spend(word.len())?;
                if let Some(product) = self.then(word[k], word[k + 1]) {
                    let merged = word[..k]
                        .iter()
                        .copied()
                        .chain([product])
                        .chain(word[k + 2..].iter().copied())
                        .collect::<Word>();
                    let value = cochain.get(&merged).copied().unwrap_or(0);
                    coefficient += if k.is_multiple_of(2) { value } else { -value };
                }
            }
            if coefficient != 0 {
                result.insert(word.clone(), coefficient);
            }
        }
        Ok(result)
    }

    fn span_boundary(
        &self,
        model: &mut Model,
        inputs: &[Word],
        targets: &[Word],
    ) -> Result<Echelon> {
        let mut span = Echelon::default();
        for word in inputs {
            let boundary = self.coboundary(model, &atom(word.clone()), targets)?;
            insert(model, &mut span, targets, &boundary)?;
        }
        Ok(span)
    }

    pub(super) fn bar_dimension(&self, model: &mut Model, grade: Grade) -> Result<usize> {
        if grade.degree > BAR_DEGREE {
            return Err(Error::plain(
                "explicit partition relative bar degrees are limited to 6",
            ));
        }
        let i = self.vertex(grade.source)?;
        let j = self.vertex(grade.target)?;
        let degree = grade.degree as usize;
        let here = self.words(model, i, j, degree)?;
        let above = self.words(model, i, j, degree + 1)?;
        let incoming = if degree == 0 {
            0
        } else {
            let below = self.words(model, i, j, degree - 1)?;
            let mut span = Echelon::default();
            for word in below {
                let value = self.coboundary(model, &atom(word), &here)?;
                if !self.coboundary(model, &value, &above)?.is_empty() {
                    return Err(Error::plain(
                        "partition relative bar coboundary does not square to zero",
                    ));
                }
                insert(model, &mut span, &here, &value)?;
            }
            span.rank()
        };
        let outgoing = self.span_boundary(model, &here, &above)?.rank();
        here.len()
            .checked_sub(incoming + outgoing)
            .ok_or_else(|| Error::plain("partition relative bar ranks exceed their cochain space"))
    }

    pub(crate) fn massey(&self, model: &mut Model, arrows: [u128; 3]) -> Result<MasseyWitness> {
        let paths = arrows.map(|index| {
            usize::try_from(index)
                .ok()
                .and_then(|i| self.arrow_paths.get(i).copied())
        });
        let [Some(a), Some(b), Some(c)] = paths else {
            return Err(Error::plain("partition Massey arrow index is out of range"));
        };
        if self.end(a) != b.start || self.end(b) != c.start {
            return Err(Error::plain(
                "partition Massey arrows must be composable in traversal order",
            ));
        }
        let x = self.then(a, b).ok_or_else(adjacent)?;
        let y = self.then(b, c).ok_or_else(adjacent)?;
        let (source, target) = (a.start, self.end(c));
        let left = atom(vec![x]);
        let right = atom(vec![y]);
        let left_targets = self.words(model, a.start, self.end(b), 2)?;
        let right_targets = self.words(model, b.start, target, 2)?;
        if self.coboundary(model, &left, &left_targets)? != atom(vec![a, b])
            || self.coboundary(model, &right, &right_targets)? != atom(vec![b, c])
        {
            return Err(Error::plain(
                "partition Massey nullhomotopies failed the bar coboundary identity",
            ));
        }
        // With δ(dual(pq))(p,q)=1, the Leibniz rule gives the two terms the same sign.
        let cycle = Cochain::from([(vec![a, y], 1), (vec![x, c], 1)]);
        let words1 = self.words(model, source, target, 1)?;
        let words2 = self.words(model, source, target, 2)?;
        let words3 = self.words(model, source, target, 3)?;
        if !self.coboundary(model, &cycle, &words3)?.is_empty() {
            return Err(Error::plain(
                "partition Massey witness is not a relative bar cocycle",
            ));
        }
        let mut boundaries = self.span_boundary(model, &words1, &words2)?;
        let boundary_rank = boundaries.rank();
        let cycle_is_boundary = boundaries.spans(model, &linear(&words2, &cycle)?)?;
        if let Some(p) = self.one_arrow(b.start, target) {
            insert(model, &mut boundaries, &words2, &atom(vec![a, p]))?;
        }
        if let Some(p) = self.one_arrow(source, c.start) {
            insert(model, &mut boundaries, &words2, &atom(vec![p, c]))?;
        }
        let indeterminacy = boundaries.rank() - boundary_rank;
        let grade = Grade {
            parameter: self.parameter,
            source: source as u128,
            target: target as u128,
            degree: 2,
        };
        let dimension = self.dimension(grade)?;
        if self.bar_dimension(model, grade)? != dimension {
            return Err(Error::plain(
                "partition Massey comparison gives different bar and projective dimensions",
            ));
        }
        let coefficient = if dimension == 0 {
            0
        } else {
            let resolution = self.simple(source);
            let first = resolution.differential(1).unwrap();
            let second = resolution.differential(2).unwrap();
            if self.then(first, second).is_some() {
                return Err(Error::plain(
                    "partition Massey comparison paths fail to compose to zero",
                ));
            }
            // The relative bar differential of -[first,second] is second[first];
            // together with [first] in degree one this compares the minimal cover.
            -cycle.get(&vec![first, second]).copied().unwrap_or(0)
        };
        if (coefficient == 0) != cycle_is_boundary {
            return Err(Error::plain(
                "partition Massey comparison failed to detect the cohomology class",
            ));
        }
        Ok(MasseyWitness {
            parameter: self.parameter,
            arrows,
            left,
            right,
            cycle,
            class: Class {
                grade,
                coefficient: Rational::integer(coefficient),
            },
            indeterminacy,
        })
    }

    fn one_arrow(&self, source: usize, target: usize) -> Option<Path> {
        let path = Path {
            start: source,
            length: 1,
        };
        (self.paths[source].len() > 1 && self.end(path) == target).then_some(path)
    }
}

fn adjacent() -> Error {
    Error::plain("partition Massey product requires both adjacent Yoneda products to vanish")
}

fn atom(word: Word) -> Cochain {
    Cochain::from([(word, 1)])
}

fn linear(targets: &[Word], cochain: &Cochain) -> Result<crate::lee_linear::Vector> {
    let mut values = Vec::with_capacity(cochain.len());
    for (word, &coefficient) in cochain {
        let index = targets
            .binary_search(word)
            .map_err(|_| Error::plain("partition bar cochain has a noncomposable path word"))?;
        values.push((index, coefficient));
    }
    Ok(values)
}

fn insert(
    model: &mut Model,
    span: &mut Echelon,
    targets: &[Word],
    cochain: &Cochain,
) -> Result<()> {
    span.insert(model, linear(targets, cochain)?, Vec::new())?;
    Ok(())
}
