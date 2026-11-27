use std::cell::{Cell, RefCell};
use std::collections::{btree_map::Entry, BTreeMap};
use std::rc::Rc;

use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::{invert, Coordinate};

pub(crate) type Polynomial = Rc<[(i128, u32)]>;

pub(crate) struct Binomial {
    prime: Prime,
    factorials: Vec<u32>,
    inverse_factorials: Vec<u32>,
}

impl Binomial {
    pub(crate) fn build(model: &mut Model, prime: Prime, max: usize) -> Result<Self> {
        let mut factorials = vec![1; max + 1];
        let mut inverse_factorials = vec![1; max + 1];
        model.spend(2 * max)?;
        for r in 1..=max {
            factorials[r] = prime.mul(factorials[r - 1], r as u32);
        }
        inverse_factorials[max] = invert(model, prime, factorials[max])?;
        for r in (1..=max).rev() {
            inverse_factorials[r - 1] = prime.mul(inverse_factorials[r], r as u32);
        }
        Ok(Self {
            prime,
            factorials,
            inverse_factorials,
        })
    }

    pub(crate) fn cells(&self) -> usize {
        self.factorials.len() + self.inverse_factorials.len()
    }

    pub(crate) fn choose(&self, r: usize, j: usize) -> u32 {
        self.prime.mul(
            self.factorials[r],
            self.prime
                .mul(self.inverse_factorials[j], self.inverse_factorials[r - j]),
        )
    }

    pub(crate) fn translate(
        &self,
        model: &mut Model,
        r: usize,
        scalar: u32,
        mut coordinate: impl FnMut(usize) -> Result<usize>,
    ) -> Result<Vec<Coordinate>> {
        let mut power = 1;
        let mut column = Vec::with_capacity(r + 1);
        for j in (0..=r).rev() {
            let index = coordinate(j)?;
            let coefficient = self.prime.mul(self.choose(r, j), power);
            column.push((index, coefficient));
            power = self.prime.mul(power, scalar);
        }
        model.spend(column.len())?;
        column.reverse();
        Ok(column)
    }

    pub(crate) fn cartier_monomial(
        &self,
        model: &mut Model,
        b: i128,
        r: usize,
        mut power: impl FnMut(&mut Model, usize) -> Result<Polynomial>,
        retained: &Cell<usize>,
    ) -> Result<Vec<(i128, usize, u32)>> {
        let p = i128::from(self.prime.get());
        let mut result = Vec::new();
        // Expand y^r = (y^p − f)^r and select Laurent exponents ≡ −1 mod p.
        for j in 0..=r {
            let factor = power(model, r - j)?;
            let mut scalar = self.choose(r, j);
            if !(r - j).is_multiple_of(2) {
                scalar = self.prime.negate(scalar);
            }
            model.spend(1 + factor.len())?;
            for &(exponent, coefficient) in factor.iter() {
                let exponent = exponent.checked_add(b).ok_or_else(exponent_error)?;
                if exponent.rem_euclid(p) == p - 1 {
                    result.push((
                        exponent.div_euclid(p),
                        j,
                        self.prime.mul(coefficient, scalar),
                    ));
                    check_storage(retained.get() + 4 * result.len())?;
                }
            }
        }
        Ok(result)
    }
}

pub(crate) struct Powers(RefCell<Vec<Polynomial>>);

impl Default for Powers {
    fn default() -> Self {
        Self(RefCell::new(vec![Rc::from([(0, 1)])]))
    }
}

impl Powers {
    pub(crate) fn get(
        &self,
        model: &mut Model,
        prime: Prime,
        terms: impl ExactSizeIterator<Item = (i128, u32)> + Clone,
        exponent: usize,
        retained: &Cell<usize>,
    ) -> Result<Polynomial> {
        let mut powers = self.0.borrow_mut();
        while powers.len() <= exponent {
            let previous = powers.last().unwrap();
            model.spend(1usize.saturating_add(previous.len().saturating_mul(terms.len())))?;
            let mut next = BTreeMap::new();
            for &(i, a) in previous.iter() {
                for (j, b) in terms.clone() {
                    let value = prime.mul(a, b);
                    model.spend(1 + next.len().max(1).ilog2() as usize)?;
                    let degree = i.checked_add(j).ok_or_else(exponent_error)?;
                    match next.entry(degree) {
                        Entry::Occupied(mut entry) => {
                            let sum = prime.add(*entry.get(), value);
                            if sum == 0 {
                                entry.remove();
                            } else {
                                *entry.get_mut() = sum;
                            }
                        }
                        Entry::Vacant(entry) => {
                            entry.insert(value);
                            check_storage(retained.get() + 8 * next.len())?;
                        }
                    }
                }
            }
            retain(model, retained, 3 + 4 * next.len())?;
            powers.push(next.into_iter().collect());
        }
        Ok(powers[exponent].clone())
    }
}

pub(crate) fn retain(model: &mut Model, cells: &Cell<usize>, extra: usize) -> Result<()> {
    let total = cells.get() + extra;
    check_storage(total)?;
    model.allocate(extra)?;
    cells.set(total);
    Ok(())
}

fn exponent_error() -> Error {
    Error::plain("Artin–Schreier differential Laurent exponent is out of range")
}
