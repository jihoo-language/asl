use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::super::{Field, Hall, Key, KeySpec, Record};
use super::Matrix;

pub(crate) mod ambient;
mod path;
mod tangent;
pub(crate) use path::degree_bound;
use path::{accumulate, Terms, Word, WORD_LIMIT};

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Chain {
    pub(crate) frame: Rc<str>,
    pub(crate) degree: i128,
    terms: Terms,
}

impl Chain {
    pub(crate) fn zero(frame: Rc<str>, degree: i128) -> Self {
        Self {
            frame,
            degree,
            terms: Terms::new(),
        }
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(
            5 + self
                .terms
                .keys()
                .map(|w| 3 + w.letters.len())
                .sum::<usize>(),
        )?;
        Ok(Value::RVPartitionDG(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<std::cmp::Ordering> {
        model.spend(1 + (3 + WORD_LIMIT) * self.terms.len().min(rhs.terms.len()))?;
        Ok(self.cmp(rhs))
    }

    pub(crate) fn render(&self) -> String {
        format!(
            "RVPartitionDG({},{}){}",
            self.frame,
            self.degree,
            path::display(&self.terms)
        )
    }
}

pub(crate) enum Literal {
    Path(u128, Vec<u128>),
    Arrow(u128, Generator),
    Potential,
    Derivative(u128),
}

pub(crate) enum Generator {
    Original,
    Reversed,
    VertexLoop,
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Differential,
    Negate,
    Scale(u128),
}

pub(crate) enum Expression {
    Literal(Label, Literal),
    Unary(Unary, Box<ParsedExpr>),
    Product(Box<[ParsedExpr; 2]>),
}

pub(crate) enum Report {
    Declare(Label, KeySpec, bool),
    Info(Label),
    Trace(Box<ParsedExpr>),
    PotentialTrace(Label, Box<ParsedExpr>),
    Tangent(Label, TangentQuery, Box<ParsedExpr>),
}

pub(crate) enum TangentQuery {
    Summary,
    Differentials,
    Cohomology(u128),
}

pub(crate) struct Spec {
    name: Rc<str>,
    hall: Rc<Hall>,
    symbols: Vec<Symbol>,
    arrows: usize,
    gradients: Vec<Word>,
    potential: Option<Terms>,
    differential: Vec<Terms>,
    components: Vec<(Vec<usize>, usize, usize)>,
}

#[derive(Clone, Copy)]
struct Symbol {
    source: usize,
    target: usize,
    degree: i128,
}

impl Spec {
    pub(crate) fn new(
        model: &mut Model,
        name: Rc<str>,
        hall: Rc<Hall>,
        ordinary: bool,
    ) -> Result<Self> {
        let p = hall.field.characteristic();
        let arrows: Vec<_> = hall
            .successor
            .iter()
            .enumerate()
            .filter_map(|(source, &t)| {
                t.map(|target| Symbol {
                    source,
                    target,
                    degree: 0,
                })
            })
            .collect();
        let number = arrows.len();
        let mut indices = vec![None; hall.vertices];
        for (i, symbol) in arrows.iter().enumerate() {
            indices[symbol.source] = Some(i);
        }
        let mut potential = ordinary.then(Terms::new);
        let mut components = Vec::new();
        let mut gradients = vec![None; number];
        let mut used = BTreeSet::new();
        for start in 0..hall.vertices {
            if hall.successor[start].is_none() {
                if hall.predecessor[start].is_some() || hall.path_heads[start].len() != 1 {
                    return Err(Error::plain("partition Ginzburg construction requires cyclic or isolated basic components"));
                }
                continue;
            }
            if used.contains(&start) {
                continue;
            }
            let mut cycle = Vec::new();
            let mut at = start;
            loop {
                if !used.insert(at) {
                    if at != start {
                        return Err(Error::plain(
                            "partition Ginzburg arrows did not form oriented cycles",
                        ));
                    }
                    break;
                }
                cycle.push(at);
                at=hall.successor[at].ok_or_else(|| Error::plain("partition Ginzburg construction requires cyclic or isolated basic components"))?;
            }
            let length = hall.path_heads[start].len();
            if cycle.iter().any(|&i| hall.path_heads[i].len() != length) {
                return Err(Error::plain(
                    "partition Ginzburg cycle requires uniform zero-relation lengths",
                ));
            }
            let h = cycle.len();
            if !(length + 1).is_multiple_of(h) {
                return Err(Error::plain("partition Ginzburg relation length plus one must be divisible by the cycle length"));
            }
            let exponent = (length + 1) / h;
            let letters: Vec<_> = cycle
                .iter()
                .cycle()
                .take(length + 1)
                .map(|&i| indices[i].unwrap() as u8)
                .collect();
            if let Some(potential) = &mut potential {
                let multiplicity = (exponent as u32) % p;
                if multiplicity == 0 {
                    return Err(Error::plain(format!("partition Ginzburg potential denominator {exponent} vanishes in characteristic {p}")));
                }
                let coefficient = hall.field.inverse(model, multiplicity)?;
                accumulate(
                    &hall.field,
                    potential,
                    Word {
                        start,
                        letters: letters.clone(),
                    },
                    coefficient,
                )?;
            }
            for (position, &vertex) in cycle.iter().enumerate() {
                let complement = (1..=length)
                    .map(|k| letters[(position + k) % letters.len()])
                    .collect();
                gradients[indices[vertex].unwrap()] = Some(Word {
                    start: hall.successor[vertex].unwrap(),
                    letters: complement,
                });
            }
            components.push((cycle, length, exponent));
        }
        let gradients = gradients
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| Error::plain("partition Ginzburg gradient missed an arrow"))?;
        let mut symbols = arrows.clone();
        symbols.extend(arrows.iter().map(|a| Symbol {
            source: a.target,
            target: a.source,
            degree: -1,
        }));
        symbols.extend((0..hall.vertices).map(|i| Symbol {
            source: i,
            target: i,
            degree: -2,
        }));
        let mut differential = vec![Terms::new(); symbols.len()];
        for (i, word) in gradients.iter().enumerate() {
            differential[number + i].insert(word.clone(), 1);
        }
        for (i, a) in arrows.iter().enumerate() {
            let star = (number + i) as u8;
            accumulate(
                &hall.field,
                &mut differential[2 * number + a.source],
                Word {
                    start: a.source,
                    letters: vec![i as u8, star],
                },
                1,
            )?;
            accumulate(
                &hall.field,
                &mut differential[2 * number + a.target],
                Word {
                    start: a.target,
                    letters: vec![star, i as u8],
                },
                hall.field.sub(0, 1),
            )?;
        }
        let spec = Self {
            name,
            hall,
            symbols,
            arrows: number,
            gradients,
            potential,
            differential,
            components,
        };
        for i in 0..spec.symbols.len() {
            let basis = spec.symbol(i)?;
            let first = spec.differentiate(model, &basis)?;
            if !spec.differentiate(model, &first)?.terms.is_empty() {
                return Err(Error::plain(
                    "partition Ginzburg differential did not square to zero on a generator",
                ));
            }
        }
        if spec.potential.is_some() {
            for a in 0..number {
                let actual = spec.derivative(model, a)?;
                if actual.terms != spec.differential[number + a] {
                    return Err(Error::plain("partition Ginzburg cyclic derivative does not recover its defining relation"));
                }
            }
        }
        model.allocate(16 + 10 * spec.symbols.len())?;
        Ok(spec)
    }

    pub(crate) fn key(&self) -> &Key {
        &self.hall.key
    }

    fn chain(&self, degree: i128, terms: Terms) -> Chain {
        Chain {
            frame: self.name.clone(),
            degree,
            terms,
        }
    }

    fn arrow_index(&self, index: u128) -> Result<usize> {
        usize::try_from(index)
            .ok()
            .filter(|&i| i < self.arrows)
            .ok_or_else(|| Error::plain("partition Ginzburg arrow index is out of range"))
    }

    pub(crate) fn literal(&self, model: &mut Model, literal: Literal) -> Result<Chain> {
        match literal {
            Literal::Path(start, letters) => {
                let start = usize::try_from(start)
                    .ok()
                    .filter(|&i| i < self.hall.vertices)
                    .ok_or_else(|| Error::plain("partition Ginzburg vertex is out of range"))?;
                if letters.len() > WORD_LIMIT {
                    return Err(path::word_limit());
                }
                let letters = letters
                    .into_iter()
                    .map(|i| {
                        u8::try_from(i)
                            .ok()
                            .filter(|&j| usize::from(j) < self.symbols.len())
                            .ok_or_else(|| {
                                Error::plain("partition Ginzburg symbol is out of range")
                            })
                    })
                    .collect::<Result<Vec<_>>>()?;
                let word = Word { start, letters };
                let degree = self.validate(model, &word)?;
                Ok(self.chain(degree, BTreeMap::from([(word, 1)])))
            }
            Literal::Arrow(index, kind) => {
                let symbol = match kind {
                    Generator::Original => self.arrow_index(index)?,
                    Generator::Reversed => self.arrows + self.arrow_index(index)?,
                    Generator::VertexLoop => {
                        2 * self.arrows
                            + usize::try_from(index)
                                .ok()
                                .filter(|&i| i < self.hall.vertices)
                                .ok_or_else(|| {
                                    Error::plain("partition Ginzburg vertex is out of range")
                                })?
                    }
                };
                self.symbol(symbol)
            }
            Literal::Potential => {
                Ok(self.chain(0, self.potential.as_ref().ok_or_else(no_potential)?.clone()))
            }
            Literal::Derivative(i) => self.derivative(model, self.arrow_index(i)?),
        }
    }

    fn symbol(&self, i: usize) -> Result<Chain> {
        let a = self
            .symbols
            .get(i)
            .ok_or_else(|| Error::plain("partition Ginzburg symbol is out of range"))?;
        Ok(self.chain(
            a.degree,
            BTreeMap::from([(
                Word {
                    start: a.source,
                    letters: vec![i as u8],
                },
                1,
            )]),
        ))
    }

    fn derivative(&self, model: &mut Model, arrow: usize) -> Result<Chain> {
        let Some(potential) = &self.potential else {
            return Ok(self.chain(0, self.differential[self.arrows + arrow].clone()));
        };
        let mut answer = Terms::new();
        for (word, &coefficient) in potential {
            for (i, &symbol) in word.letters.iter().enumerate() {
                model.spend(1 + word.letters.len())?;
                if usize::from(symbol) == arrow {
                    let child = Word {
                        start: self.symbols[arrow].target,
                        letters: word.letters[i + 1..]
                            .iter()
                            .chain(&word.letters[..i])
                            .copied()
                            .collect(),
                    };
                    accumulate(&self.hall.field, &mut answer, child, coefficient)?;
                }
            }
        }
        Ok(self.chain(0, answer))
    }

    fn path_matrix<'a>(
        &self,
        model: &mut Model,
        dims: &[usize],
        word: &Word,
        arrow: impl Fn(usize) -> &'a Matrix,
    ) -> Result<Matrix> {
        let mut product = Matrix::identity(dims[word.start]);
        for &letter in &word.letters {
            product = arrow(usize::from(letter)).compose(model, &self.hall.field, &product)?;
        }
        Ok(product)
    }

    fn representation_trace<'a>(
        &self,
        model: &mut Model,
        dims: &[usize],
        arrow: impl Fn(usize) -> &'a Matrix,
    ) -> Result<u32> {
        let field = &self.hall.field;
        let mut trace = 0;
        for (word, &coefficient) in self.potential.as_ref().ok_or_else(no_potential)? {
            let matrix = self.path_matrix(model, dims, word, &arrow)?;
            trace = field.add(trace, field.mul(coefficient, matrix_trace(field, &matrix)));
        }
        Ok(trace)
    }

    pub(crate) fn info(&self) -> String {
        let symbols: Vec<_> = self
            .symbols
            .iter()
            .enumerate()
            .map(|(i, a)| (i, a.source, a.target, a.degree))
            .collect();
        format!("rvPartitionDGInfo[{}] = q={}, characteristic={}, symbols={symbols:?}, cyclicComponents={:?}, potential={}",self.name,self.hall.field.order(),self.hall.field.characteristic(),self.components,self.potential.is_some())
    }
}

fn no_potential() -> Error {
    Error::plain("partition Ginzburg declaration has integral differentiated data only")
}

fn matrix_trace(field: &Field, m: &Matrix) -> u32 {
    (0..m.rows).fold(0, |sum, i| field.add(sum, m.get(i, i)))
}
