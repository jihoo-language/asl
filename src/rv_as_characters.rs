use std::cell::OnceCell;
use std::rc::Rc;

use crate::binary_linear::check_storage;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::{rv_cyclotomic as cyclo, rv_newton};

type Character = Rc<[Vec<i128>]>;

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Character(usize),
    CharacterNewton(usize),
    LPolynomial,
    Newton,
    Points(usize),
}

pub(crate) struct Characters {
    prime: Prime,
    factor_degree: usize,
    genus: usize,
    character: OnceCell<Character>,
    numerator: OnceCell<Rc<[i128]>>,
}

impl Characters {
    pub(crate) fn new(prime: Prime, factor_degree: usize, genus: usize) -> Self {
        Self {
            prime,
            factor_degree,
            genus,
            character: OnceCell::new(),
            numerator: OnceCell::new(),
        }
    }

    pub(crate) fn order(&self) -> usize {
        self.prime.get() as usize
    }

    fn checked_character(&self, character: usize) -> Result<()> {
        if !(1..self.order()).contains(&character) {
            return Err(Error::plain(
                "odd Artin–Schreier character index must lie in 1..p-1",
            ));
        }
        Ok(())
    }

    pub(crate) fn coefficient(
        &self,
        model: &mut Model,
        j: usize,
        k: usize,
        build: impl FnOnce(&mut Model) -> Result<Vec<Vec<i128>>>,
    ) -> Result<Value> {
        self.checked_character(j)?;
        if k > self.factor_degree {
            return Ok(cyclo::zero());
        }
        if k == 0 {
            return cyclo::persist(model, vec![1]);
        }
        let factor = self.primitive(model, build)?;
        let coefficient = self.conjugate(model, j, &factor[k])?;
        cyclo::persist(model, coefficient)
    }

    fn conjugate(
        &self,
        model: &mut Model,
        character: usize,
        coefficient: &[i128],
    ) -> Result<Vec<i128>> {
        if character == 1 || coefficient.len() <= 1 {
            model.spend(coefficient.len())?;
            return Ok(coefficient.to_vec());
        }
        cyclo::substitute_root(model, self.order(), coefficient, character)
    }

    pub(crate) fn report(
        &self,
        model: &mut Model,
        name: &str,
        operation: &str,
        query: Query,
        build: impl FnOnce(&mut Model) -> Result<Vec<Vec<i128>>>,
    ) -> Result<String> {
        let value = match query {
            Query::Character(j) => {
                self.checked_character(j)?;
                let factor = self.primitive(model, build)?;
                let mut conjugate = Vec::with_capacity(factor.len());
                for coefficient in factor.iter() {
                    let mut image = self.conjugate(model, j, coefficient)?;
                    model.spend(self.order())?;
                    image.resize(self.order() - 1, 0);
                    conjugate.push(image);
                }
                return Ok(format!("{operation}[{name}, {j}] = {conjugate:?}"));
            }
            Query::CharacterNewton(j) => {
                self.checked_character(j)?;
                // The valuation extends uniquely to Q(zeta_p), so all
                // cyclotomic conjugates have the same Newton polygon.
                let factor = self.primitive(model, build)?;
                let valuations = factor
                    .iter()
                    .map(|a| uniformizer_order(model, self.prime, a))
                    .collect::<Result<Vec<_>>>()?;
                let slopes = rv_newton::scaled_slopes(model, valuations, self.order() - 1)?;
                return Ok(format!(
                    "{operation}[{name}, {j}] = [{}]",
                    slopes.join(", ")
                ));
            }
            Query::LPolynomial => format!("{:?}", self.l_polynomial(model, build)?),
            Query::Newton => {
                let polynomial = self.l_polynomial(model, build)?;
                let valuations = polynomial
                    .iter()
                    .map(|&c| p_order(model, self.prime, c))
                    .collect::<Result<Vec<_>>>()?;
                let slopes = rv_newton::scaled_slopes(model, valuations, 1)?;
                format!("[{}]", slopes.join(", "))
            }
            Query::Points(m) => {
                let points = self.points(model, m, build)?;
                return Ok(format!("{operation}[{name}, {m}] = {points}"));
            }
        };
        Ok(format!("{operation}[{name}] = {value}"))
    }

    fn primitive(
        &self,
        model: &mut Model,
        build: impl FnOnce(&mut Model) -> Result<Vec<Vec<i128>>>,
    ) -> Result<Character> {
        if let Some(character) = self.character.get() {
            return Ok(character.clone());
        }
        let n = self.factor_degree;
        let factor = build(model)?;
        if factor.len() != n + 1 || factor[n].is_empty() {
            return Err(Error::plain(
                "odd Artin–Schreier character factor has the wrong degree",
            ));
        }
        model.allocate(factor.iter().map(|c| 3 + 2 * c.len()).sum())?;
        let factor: Character = factor.into();
        Ok(self.character.get_or_init(|| factor).clone())
    }

    fn l_polynomial(
        &self,
        model: &mut Model,
        build: impl FnOnce(&mut Model) -> Result<Vec<Vec<i128>>>,
    ) -> Result<Rc<[i128]>> {
        if let Some(polynomial) = self.numerator.get() {
            return Ok(polynomial.clone());
        }
        let p = self.order();
        let leading = u32::try_from(self.genus)
            .ok()
            .and_then(|g| (p as i128).checked_pow(g))
            .ok_or_else(coefficient_error)?;
        let mut polynomial = vec![vec![1]];
        if self.genus != 0 {
            let factor = self.primitive(model, build)?;
            for j in 1..p {
                let conjugate = factor
                    .iter()
                    .map(|coefficient| self.conjugate(model, j, coefficient))
                    .collect::<Result<Vec<_>>>()?;
                polynomial = polynomial_product(model, p, &polynomial, &conjugate)?;
            }
        }
        let polynomial = polynomial
            .into_iter()
            .map(|coefficient| match coefficient.as_slice() {
                [] => Ok(0),
                &[c] => Ok(c),
                _ => Err(Error::plain(
                    "odd Artin–Schreier character norm did not descend to integers",
                )),
            })
            .collect::<Result<Vec<_>>>()?;
        if polynomial.len() != 2 * self.genus + 1 || polynomial.last() != Some(&leading) {
            return Err(Error::plain(
                "odd Artin–Schreier numerator disagrees with the curve genus",
            ));
        }
        model.allocate(2 * polynomial.len())?;
        Ok(self.numerator.get_or_init(|| polynomial.into()).clone())
    }

    fn points(
        &self,
        model: &mut Model,
        extension: usize,
        build: impl FnOnce(&mut Model) -> Result<Vec<Vec<i128>>>,
    ) -> Result<u128> {
        if extension == 0 {
            return Err(Error::plain(
                "odd Artin–Schreier extension degree must be positive",
            ));
        }
        let p = self.order();
        let base = u32::try_from(extension)
            .ok()
            .and_then(|m| (p as u128).checked_pow(m))
            .and_then(|power| power.checked_add(1))
            .ok_or_else(points_error)?;
        if self.genus == 0 {
            return Ok(base);
        }
        let factor = self.primitive(model, build)?;
        let sum = power_trace(model, p, &factor, extension)?;
        // Q(zeta_p)/Q has trace p-1 on 1 and -1 on the other basis powers.
        let mut trace = sum
            .first()
            .copied()
            .unwrap_or(0)
            .checked_mul((p - 1) as i128)
            .ok_or_else(coefficient_error)?;
        for &c in sum.iter().skip(1) {
            trace = trace.checked_sub(c).ok_or_else(coefficient_error)?;
        }
        let result = if trace >= 0 {
            base.checked_add(trace as u128)
        } else {
            base.checked_sub(trace.unsigned_abs())
        };
        result.ok_or_else(points_error)
    }
}

fn power_trace(
    model: &mut Model,
    p: usize,
    factor: &[Vec<i128>],
    extension: usize,
) -> Result<Vec<i128>> {
    let degree = factor.len() - 1;
    check_storage((degree + 1).saturating_mul(3 + 2 * p))?;
    let mut window = vec![Vec::new(); degree];
    for power in 1..=extension {
        model.spend(1)?;
        let mut sum = Vec::new();
        if let Some(coefficient) = factor.get(power) {
            model.spend(coefficient.len())?;
            sum = coefficient
                .iter()
                .map(|&c| c.checked_mul(power as i128).ok_or_else(coefficient_error))
                .collect::<Result<Vec<_>>>()?;
        }
        // Negative Frobenius power trace S_k = k c_k - sum c_i S_(k-i).
        for (i, c) in factor
            .iter()
            .enumerate()
            .take(power.min(degree + 1))
            .skip(1)
        {
            let product = cyclo::product_coefficients(model, p, c, &window[(power - i) % degree])?;
            cyclo::accumulate(model, &mut sum, &product, true)?;
        }
        if power == extension {
            return Ok(sum);
        }
        window[power % degree] = sum;
    }
    unreachable!()
}

fn polynomial_product(
    model: &mut Model,
    p: usize,
    left: &[Vec<i128>],
    right: &[Vec<i128>],
) -> Result<Vec<Vec<i128>>> {
    let length = left.len() + right.len() - 1;
    check_storage((left.len() + right.len() + length).saturating_mul(3 + 2 * p))?;
    model.spend(left.len() * right.len())?;
    let mut product = vec![Vec::new(); length];
    for (i, x) in left.iter().enumerate() {
        for (j, y) in right.iter().enumerate() {
            let coefficient = cyclo::product_coefficients(model, p, x, y)?;
            cyclo::accumulate(model, &mut product[i + j], &coefficient, false)?;
        }
    }
    Ok(product)
}

fn p_order(model: &mut Model, prime: Prime, value: i128) -> Result<Option<i64>> {
    model.spend(1)?;
    if value == 0 {
        return Ok(None);
    }
    let p = i128::from(prime.get());
    let mut value = value;
    let mut order = 0;
    while value % p == 0 {
        model.spend(1)?;
        value /= p;
        order += 1;
    }
    Ok(Some(order))
}

fn uniformizer_order(model: &mut Model, prime: Prime, coefficient: &[i128]) -> Result<Option<i64>> {
    let mut common = None;
    for &value in coefficient {
        if let Some(order) = p_order(model, prime, value)? {
            common = Some(common.map_or(order, |old: i64| old.min(order)));
        }
    }
    let Some(common) = common else {
        return Ok(None);
    };
    let p = i128::from(prime.get());
    let divisor = p.checked_pow(common as u32).ok_or_else(coefficient_error)?;
    model.spend(coefficient.len())?;
    let mut residue: Vec<_> = coefficient
        .iter()
        .map(|&a| (a / divisor).rem_euclid(p) as u32)
        .collect();
    // Divide by x-1 over F_p until the remainder is nonzero. The degree
    // is less than p-1, so its first nonzero shifted coefficient gives
    // the residual uniformizer order without computing large binomials.
    let mut order = common * (i64::from(prime.get()) - 1);
    while !residue.is_empty() {
        model.spend(residue.len())?;
        for i in (1..residue.len()).rev() {
            residue[i - 1] = prime.add(residue[i - 1], residue[i]);
        }
        if residue[0] != 0 {
            return Ok(Some(order));
        }
        residue.rotate_left(1);
        residue.pop();
        order += 1;
    }
    Err(Error::plain(
        "odd Artin–Schreier normalized cyclotomic residue vanished",
    ))
}

pub(crate) fn coefficient_error() -> Error {
    Error::plain("odd Artin–Schreier cyclotomic coefficient is out of range")
}

fn points_error() -> Error {
    Error::plain("odd Artin–Schreier point count is out of range")
}
