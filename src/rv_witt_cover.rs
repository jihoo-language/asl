use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::binary_field::BinaryField;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_as_cover::{affine_value, integral_l_polynomial, point_count_field};
use crate::rv_newton;
use crate::rv_poly::{self, Polynomial};

type Gaussian = (i128, i128);

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Breaks,
    Genus,
    Points(usize),
    LPolynomial,
    Newton,
    CharacterL,
}

impl Query {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Breaks => "rvWittBreaks",
            Self::Genus => "rvWittGenus",
            Self::Points(_) => "rvWittPoints",
            Self::LPolynomial => "rvWittLPolynomial",
            Self::Newton => "rvWittNewton",
            Self::CharacterL => "rvWittCharacterL",
        }
    }
}

pub(crate) struct Cover {
    f: Polynomial,
    g: Polynomial,
    breaks: [usize; 2],
    genus: usize,
    traces: RefCell<BTreeMap<usize, [usize; 4]>>,
    l_polynomial: RefCell<Option<Rc<[i128]>>>,
    character_l: RefCell<Option<Rc<[Gaussian]>>>,
}

impl Cover {
    pub(crate) fn build(model: &mut Model, forms: Vec<Polynomial>) -> Result<Self> {
        let [f, g] = forms
            .try_into()
            .map_err(|_| Error::plain("rvwittcover requires exactly two polynomials"))?;
        let m = f.last().map_or(0, |term| term[0]);
        if m.is_multiple_of(2) {
            return Err(Error::plain(
                "rvwittcover first polynomial must have positive odd degree",
            ));
        }
        let n = g.last().map_or(0, |term| term[0]);
        if n != 0 && n.is_multiple_of(2) {
            return Err(Error::plain(
                "rvwittcover second polynomial must be zero, constant, or of odd degree",
            ));
        }
        let twice_m = m
            .checked_mul(2)
            .ok_or_else(|| Error::plain("Witt upper ramification break is out of range"))?;
        let upper = twice_m.max(n);
        let genus = (m as u128 + 2 * upper as u128 - 3) / 2;
        let genus =
            usize::try_from(genus).map_err(|_| Error::plain("Witt curve genus is out of range"))?;
        model.allocate(5 + rv_poly::cells(&f) + rv_poly::cells(&g))?;
        Ok(Self {
            f,
            g,
            breaks: [m, upper],
            genus,
            traces: RefCell::new(BTreeMap::new()),
            l_polynomial: RefCell::new(None),
            character_l: RefCell::new(None),
        })
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, query: Query) -> Result<String> {
        let value = match query {
            Query::Breaks => format!("{:?}", self.breaks),
            Query::Genus => self.genus.to_string(),
            Query::Points(extension) => {
                let points = self.points(model, extension)?;
                return Ok(format!("rvWittPoints[{name}, {extension}] = {points}"));
            }
            Query::LPolynomial => format!("{:?}", self.l_polynomial(model)?),
            Query::Newton => {
                let coefficients = self.l_polynomial(model)?;
                format!(
                    "[{}]",
                    rv_newton::integral_slopes(model, &coefficients)?.join(", ")
                )
            }
            Query::CharacterL => format!("{:?}", self.character_l(model)?),
        };
        Ok(format!("{}[{name}] = {value}", query.name()))
    }

    fn points(&self, model: &mut Model, extension: usize) -> Result<usize> {
        Ok(1 + 4 * self.trace_histogram(model, extension)?[0])
    }

    fn trace_histogram(&self, model: &mut Model, extension: usize) -> Result<[usize; 4]> {
        if let Some(&histogram) = self.traces.borrow().get(&extension) {
            return Ok(histogram);
        }
        let field = point_count_field(model, extension)?;
        let mut histogram = [0; 4];
        for parameter in 0..field.order() {
            let a = affine_value(model, field, &self.f, parameter)?;
            let b = affine_value(model, field, &self.g, parameter)?;
            let character_index = witt_trace_index(model, field, a, b)?;
            histogram[character_index] += 1;
        }
        model.allocate(5)?;
        self.traces.borrow_mut().insert(extension, histogram);
        Ok(histogram)
    }

    fn l_polynomial(&self, model: &mut Model) -> Result<Rc<[i128]>> {
        if let Some(coefficients) = self.l_polynomial.borrow().as_ref() {
            return Ok(coefficients.clone());
        }
        let coefficients = integral_l_polynomial(
            model,
            self.genus,
            "Witt curve L-polynomial",
            |model, extension| self.points(model, extension),
        )?;
        model.allocate(1 + 2 * coefficients.len())?;
        let coefficients: Rc<[_]> = coefficients.into();
        *self.l_polynomial.borrow_mut() = Some(coefficients.clone());
        Ok(coefficients)
    }

    fn character_l(&self, model: &mut Model) -> Result<Rc<[Gaussian]>> {
        if let Some(coefficients) = self.character_l.borrow().as_ref() {
            return Ok(coefficients.clone());
        }
        let degree = self.breaks[1] - 1;
        if degree >= usize::BITS as usize {
            return Err(Error::plain(
                "Witt character L-polynomial exceeds the point-count extension range",
            ));
        }
        let mut sums = vec![(0, 0); degree + 1];
        for (extension, sum) in sums.iter_mut().enumerate().skip(1) {
            let [one, i, minus_one, minus_i] = self.trace_histogram(model, extension)?;
            *sum = (one as i128 - minus_one as i128, i as i128 - minus_i as i128);
        }
        let mut coefficients = vec![(1, 0)];
        for n in 1..=degree {
            model.spend(n)?;
            let mut numerator = (0, 0);
            for j in 1..=n {
                numerator =
                    gaussian_add(numerator, gaussian_multiply(sums[j], coefficients[n - j])?)?;
            }
            let divisor = n as i128;
            if numerator.0 % divisor != 0 || numerator.1 % divisor != 0 {
                return Err(Error::plain(
                    "Witt character identity did not give a Gaussian integer coefficient",
                ));
            }
            coefficients.push((numerator.0 / divisor, numerator.1 / divisor));
        }
        model.allocate(1 + 4 * coefficients.len())?;
        let coefficients: Rc<[_]> = coefficients.into();
        *self.character_l.borrow_mut() = Some(coefficients.clone());
        Ok(coefficients)
    }
}

fn witt_trace_index(
    model: &mut Model,
    field: BinaryField,
    mut a: usize,
    mut b: usize,
) -> Result<usize> {
    let (mut first, mut second) = (a, b);
    model.spend(3 * field.degree() * (field.degree() - 1))?;
    for _ in 1..field.degree() {
        a = field.multiply(a, a);
        b = field.multiply(b, b);
        // Add the next Frobenius conjugate, using the original first
        // component in the Witt carry into the second component.
        second ^= b ^ field.multiply(first, a);
        first ^= a;
    }
    if first > 1 || second > 1 {
        return Err(Error::plain(
            "Witt trace did not land in the length-two Witt ring of F2",
        ));
    }
    Ok(first + 2 * second)
}

fn coefficient_limit() -> Error {
    Error::plain("Witt character L-polynomial coefficient is out of range")
}

fn gaussian_add((a, b): Gaussian, (c, d): Gaussian) -> Result<Gaussian> {
    Ok((
        a.checked_add(c).ok_or_else(coefficient_limit)?,
        b.checked_add(d).ok_or_else(coefficient_limit)?,
    ))
}

fn gaussian_multiply((a, b): Gaussian, (c, d): Gaussian) -> Result<Gaussian> {
    let ac = a.checked_mul(c).ok_or_else(coefficient_limit)?;
    let bd = b.checked_mul(d).ok_or_else(coefficient_limit)?;
    let ad = a.checked_mul(d).ok_or_else(coefficient_limit)?;
    let bc = b.checked_mul(c).ok_or_else(coefficient_limit)?;
    Ok((
        ac.checked_sub(bd).ok_or_else(coefficient_limit)?,
        ad.checked_add(bc).ok_or_else(coefficient_limit)?,
    ))
}
