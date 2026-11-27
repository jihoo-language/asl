use std::collections::btree_map::Entry;
use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::Model;

use super::{arithmetic, Chamber, Q, Q_MINUS_ONE};

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Word {
    pub(crate) length: usize,
    pub(crate) first: u8,
}

impl Word {
    pub(crate) const ID: Self = Self {
        length: 0,
        first: 0,
    };

    pub(crate) fn literal(panels: &[i128]) -> Result<Self> {
        if panels.iter().any(|&panel| !matches!(panel, 0 | 1))
            || panels.windows(2).any(|pair| pair[0] == pair[1])
        {
            return Err(Error::plain(
                "cubic affine Hecke words must alternate between panels 0 and 1",
            ));
        }
        Ok(Self {
            length: panels.len(),
            first: panels.first().copied().unwrap_or(0) as u8,
        })
    }

    pub(crate) fn panel(self, place: usize) -> u8 {
        (self.first + (place % 2) as u8) % 2
    }

    pub(crate) fn last(self) -> Option<u8> {
        self.length.checked_sub(1).map(|i| self.panel(i))
    }

    fn twist(self, power: i64) -> Self {
        if self.length == 0 || power.rem_euclid(2) == 0 {
            self
        } else {
            Self {
                length: self.length,
                first: 1 - self.first,
            }
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Basis {
    pub(crate) word: Word,
    pub(crate) omega: i64,
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Series<K>(BTreeMap<K, Rational>);

impl<K: Ord + Clone> Series<K> {
    pub(crate) fn zero() -> Self {
        Self(BTreeMap::new())
    }

    pub(crate) fn basis(key: K) -> Self {
        Self(BTreeMap::from([(key, Rational::integer(1))]))
    }

    pub(crate) fn coefficient(&self, key: &K) -> Rational {
        self.0.get(key).copied().unwrap_or(Rational::ZERO)
    }

    pub(crate) fn entries(&self) -> &BTreeMap<K, Rational> {
        &self.0
    }

    pub(crate) fn insert(&mut self, model: &mut Model, key: K, c: Rational) -> Result<()> {
        model.spend(1)?;
        if c.numerator == 0 {
            return Ok(());
        }
        match self.0.entry(key) {
            Entry::Vacant(slot) => {
                slot.insert(c);
            }
            Entry::Occupied(mut slot) => {
                let sum = slot.get().add(c).map_err(|_| arithmetic())?;
                if sum.numerator == 0 {
                    slot.remove();
                } else {
                    *slot.get_mut() = sum;
                }
            }
        }
        Ok(())
    }

    pub(crate) fn add_scaled(&mut self, model: &mut Model, rhs: &Self, c: Rational) -> Result<()> {
        if c.numerator == 0 {
            return Ok(());
        }
        for (key, &coefficient) in rhs.entries() {
            self.insert(
                model,
                key.clone(),
                coefficient.multiply(c).map_err(|_| arithmetic())?,
            )?;
        }
        Ok(())
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let mut out = self.clone();
        out.add_scaled(model, rhs, Rational::integer(1))?;
        Ok(out)
    }

    pub(crate) fn scale(&self, model: &mut Model, c: Rational) -> Result<Self> {
        let mut out = Self::zero();
        out.add_scaled(model, self, c)?;
        Ok(out)
    }
}

pub(crate) type Iwahori = Series<Basis>;
pub(crate) type Spherical = Series<usize>;
pub(crate) type Laurent = Series<i64>;
pub(crate) type Vector = Series<Chamber>;

fn multiply_panel(model: &mut Model, source: &Series<Word>, panel: u8) -> Result<Series<Word>> {
    let mut out = Series::zero();
    for (&word, &c) in source.entries() {
        if word.last() == Some(panel) {
            let former = if word.length == 1 {
                Word::ID
            } else {
                Word {
                    length: word.length - 1,
                    first: word.first,
                }
            };
            out.insert(
                model,
                word,
                c.multiply(Rational::integer(Q_MINUS_ONE))
                    .map_err(|_| arithmetic())?,
            )?;
            out.insert(
                model,
                former,
                c.multiply(Rational::integer(Q)).map_err(|_| arithmetic())?,
            )?;
        } else {
            let next = Word {
                length: word.length.checked_add(1).ok_or_else(arithmetic)?,
                first: if word.length == 0 { panel } else { word.first },
            };
            out.insert(model, next, c)?;
        }
    }
    Ok(out)
}

pub(crate) fn combine(model: &mut Model, lhs: &Iwahori, rhs: &Iwahori) -> Result<Iwahori> {
    let mut out = Iwahori::zero();
    for (&a, &x) in lhs.entries() {
        for (&b, &y) in rhs.entries() {
            let omega = a.omega.checked_add(b.omega).ok_or_else(arithmetic)?;
            let word_b = b.word.twist(a.omega);
            let c = x.multiply(y).map_err(|_| arithmetic())?;
            let mut words = Series::zero();
            words.insert(model, a.word, c)?;
            for i in 0..word_b.length {
                words = multiply_panel(model, &words, word_b.panel(i))?;
            }
            for (&word, &c) in words.entries() {
                out.insert(model, Basis { word, omega }, c)?;
            }
        }
    }
    Ok(out)
}

pub(crate) fn bernstein(model: &mut Model) -> Result<Iwahori> {
    let mut out = Iwahori::zero();
    for first in 0..=1 {
        out.insert(
            model,
            Basis {
                word: Word { length: 2, first },
                omega: 0,
            },
            Rational::integer(1),
        )?;
        out.insert(
            model,
            Basis {
                word: Word { length: 1, first },
                omega: 0,
            },
            Rational::integer(-Q_MINUS_ONE),
        )?;
    }
    Ok(out)
}

pub(crate) fn compression(model: &mut Model, input: &Spherical) -> Result<Iwahori> {
    let mut out = Iwahori::zero();
    for (&n, &c) in input.entries() {
        if n % 2 == 1 {
            return Err(Error::plain(
                "even-vertex Iwahori compression requires an even spherical support",
            ));
        }
        let weight = c
            .multiply(Rational::ratio(1, Q + 1).map_err(|_| arithmetic())?)
            .map_err(|_| arithmetic())?;
        if n == 0 {
            for length in [0, 1] {
                out.insert(
                    model,
                    Basis {
                        word: Word { first: 0, length },
                        omega: 0,
                    },
                    weight,
                )?;
            }
        } else {
            for (first, length) in [
                (1, n - 1),
                (0, n),
                (1, n),
                (0, n.checked_add(1).ok_or_else(arithmetic)?),
            ] {
                out.insert(
                    model,
                    Basis {
                        word: Word { first, length },
                        omega: 0,
                    },
                    weight,
                )?;
            }
        }
    }
    Ok(out)
}

fn radial_step(model: &mut Model, input: &Spherical) -> Result<Spherical> {
    let mut out = Spherical::zero();
    for (&n, &c) in input.entries() {
        let up = n.checked_add(1).ok_or_else(arithmetic)?;
        out.insert(model, up, c)?;
        if n > 0 {
            let factor = if n == 1 { Q + 1 } else { Q };
            out.insert(
                model,
                n - 1,
                c.multiply(Rational::integer(factor))
                    .map_err(|_| arithmetic())?,
            )?;
        }
    }
    Ok(out)
}

pub(crate) fn radial_product(
    model: &mut Model,
    lhs: &Spherical,
    rhs: &Spherical,
) -> Result<Spherical> {
    if lhs.entries().is_empty() || rhs.entries().is_empty() {
        return Ok(Spherical::zero());
    }
    let (source, factors) = if lhs.entries().last_key_value() > rhs.entries().last_key_value() {
        (lhs, rhs)
    } else {
        (rhs, lhs)
    };
    let (&max, _) = factors
        .entries()
        .last_key_value()
        .expect("nonempty spherical support");
    let mut result = Spherical::zero();
    result.add_scaled(model, source, factors.coefficient(&0))?;
    if max == 0 {
        return Ok(result);
    }
    let mut prior = source.clone();
    let mut current = radial_step(model, source)?;
    result.add_scaled(model, &current, factors.coefficient(&1))?;
    for n in 2..=max {
        model.spend(1)?;
        let mut next = radial_step(model, &current)?;
        let down = if n == 2 { Q + 1 } else { Q };
        next.add_scaled(model, &prior, Rational::integer(-down))?;
        result.add_scaled(model, &next, factors.coefficient(&n))?;
        prior = current;
        current = next;
    }
    Ok(result)
}

pub(crate) fn satake(model: &mut Model, source: &Spherical) -> Result<Laurent> {
    let mut output = Laurent::zero();
    for (&n, &weight) in source.entries() {
        let height = i64::try_from(n).map_err(|_| arithmetic())?;
        let power = Q
            .checked_pow(u32::try_from(n).map_err(|_| arithmetic())?)
            .ok_or_else(arithmetic)?;
        let mut add = |model: &mut Model, z, c| {
            output.insert(
                model,
                z,
                weight
                    .multiply(Rational::integer(c))
                    .map_err(|_| arithmetic())?,
            )
        };
        if n == 0 {
            add(model, 0, 1)?;
        } else {
            // A path can follow the reference end to level k, then leave
            // it through q-1 side branches. At the root there are q.
            add(model, -height, power)?;
            add(model, height, 1)?;
            let mut count = power / (Q * Q);
            for k in 1..height {
                add(model, 2 * k - height, Q_MINUS_ONE * count)?;
                count /= Q;
            }
        }
    }
    Ok(output)
}

pub(crate) fn laurent_product(model: &mut Model, lhs: &Laurent, rhs: &Laurent) -> Result<Laurent> {
    let mut output = Laurent::zero();
    for (&a, &x) in lhs.entries() {
        for (&b, &y) in rhs.entries() {
            output.insert(
                model,
                a.checked_add(b).ok_or_else(arithmetic)?,
                x.multiply(y).map_err(|_| arithmetic())?,
            )?;
        }
    }
    Ok(output)
}
