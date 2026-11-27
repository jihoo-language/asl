use super::super::nakajima::{Factor, Mode, State, COLORS};
use super::{accumulate, checked, monomial_index, nonzero, ratio, Class, Form, MIDDLE, RANK};
use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::Model;
use crate::rv_lens::gcd;

const POINT: u8 = COLORS as u8 - 1;
const DELTA: usize = RANK - 1;
const DELTA_SQUARE: usize = MIDDLE - 1;

#[derive(Clone, Copy)]
enum BasisTerm {
    Pair(u8, u8),
    Double(u8),
}

impl BasisTerm {
    fn from_factors(factors: &[Factor]) -> Self {
        match factors {
            [Factor {
                mode: Mode { length: 2, color },
                exponent: 1,
            }] => Self::Double(*color),
            [Factor {
                mode: Mode { length: 1, color },
                exponent: 2,
            }] => Self::Pair(*color, *color),
            [Factor {
                mode:
                    Mode {
                        length: 1,
                        color: i,
                    },
                exponent: 1,
            }, Factor {
                mode:
                    Mode {
                        length: 1,
                        color: j,
                    },
                exponent: 1,
            }] => Self::Pair(*i, *j),
            _ => unreachable!("the Nakajima state has weight two"),
        }
    }

    fn factors(self) -> Vec<Factor> {
        let factor = |length, color| Factor {
            mode: Mode { length, color },
            exponent: 1,
        };
        match self {
            Self::Pair(i, j) => vec![factor(1, i), factor(1, j)],
            Self::Double(i) => vec![factor(2, i)],
        }
    }
}

fn basis_index(color: u8) -> usize {
    usize::from(color) - 1
}

fn basis_color(index: usize) -> u8 {
    (index + 1) as u8
}

pub(crate) fn forward(
    model: &mut Model,
    form: &Form,
    state: &State,
    mut result: Class,
) -> Result<Class> {
    use BasisTerm::{Double, Pair};
    let mut point_coefficient = Rational::ZERO;
    let one = Rational::integer(1);
    let two = Rational::integer(2);
    let half = ratio(1, 2)?;
    for (factors, coefficient) in state.terms() {
        model.spend(1)?;
        let scalar = Rational::integer(coefficient);
        let (index, multiplier) = match BasisTerm::from_factors(factors) {
            Pair(0, 0) => (0, two),
            Pair(0, POINT) => {
                point_coefficient = checked(point_coefficient.add(scalar))?;
                continue;
            }
            Pair(POINT, POINT) => (0, one),
            Pair(0, j) | Pair(j, POINT) => (basis_index(j), one),
            Pair(i, j) => {
                let (i, j) = (basis_index(i), basis_index(j));
                accumulate(
                    &mut point_coefficient,
                    Rational::integer(-form.gram[i][j]),
                    scalar,
                )?;
                (monomial_index(i, j), one)
            }
            Double(0) => (DELTA, two),
            Double(POINT) => (DELTA, half),
            Double(i) => (monomial_index(basis_index(i), DELTA), one),
        };
        accumulate(&mut result.coordinates[index], multiplier, scalar)?;
    }

    if point_coefficient != Rational::ZERO {
        // P = (G^vee + 2 delta^2)/20. The off-diagonal coefficients
        // of G^vee combine the two entries of the inverse Gram matrix.
        let inverse = form.inverse(model)?;
        model.spend(MIDDLE)?;
        for (index, &(i, j)) in form.monomials[..DELTA_SQUARE].iter().enumerate() {
            if j == DELTA {
                continue;
            }
            let weight = ratio(1, if i == j { 20 } else { 10 })?;
            let weight = checked(weight.multiply(*inverse.entry(i, j)))?;
            accumulate(&mut result.coordinates[index], point_coefficient, weight)?;
        }
        accumulate(
            &mut result.coordinates[DELTA_SQUARE],
            point_coefficient,
            ratio(1, 10)?,
        )?;
    }
    Ok(result)
}

fn inverse_terms(
    model: &mut Model,
    form: &Form,
    class: &Class,
) -> Result<Vec<(BasisTerm, Rational)>> {
    use BasisTerm::{Double, Pair};
    model.spend(1 + class.coordinates.len())?;
    let mut terms = Vec::with_capacity(class.coordinates.len());
    let two = Rational::integer(2);
    let half = ratio(1, 2)?;
    match class.degree {
        0 => terms.push((Pair(0, 0), checked(half.multiply(class.coordinates[0]))?)),
        2 => {
            for (i, value) in nonzero(&class.coordinates[..DELTA]) {
                terms.push((Pair(0, basis_color(i)), value));
            }
            terms.push((Double(0), checked(half.multiply(class.coordinates[DELTA]))?));
        }
        4 => {
            let delta_square = class.coordinates[DELTA_SQUARE];
            let inverse = if delta_square == Rational::ZERO {
                None
            } else {
                Some(form.inverse(model)?)
            };
            let mut point_coefficient = checked(delta_square.multiply(Rational::integer(-1)))?;
            for (&(i, j), &coefficient) in form.monomials[..DELTA_SQUARE]
                .iter()
                .zip(&class.coordinates)
            {
                if j == DELTA {
                    terms.push((Double(basis_color(i)), coefficient));
                    continue;
                }
                accumulate(
                    &mut point_coefficient,
                    Rational::integer(form.gram[i][j]),
                    coefficient,
                )?;
                let mut corrected = coefficient;
                if let Some(inverse) = inverse {
                    let weight = ratio(-1, if i == j { 2 } else { 1 })?;
                    let weight = checked(weight.multiply(*inverse.entry(i, j)))?;
                    accumulate(&mut corrected, delta_square, weight)?;
                }
                terms.push((Pair(basis_color(i), basis_color(j)), corrected));
            }
            terms.push((Pair(0, POINT), point_coefficient));
        }
        6 => {
            for (i, value) in nonzero(&class.coordinates[..DELTA]) {
                terms.push((Pair(basis_color(i), POINT), value));
            }
            terms.push((
                Double(POINT),
                checked(two.multiply(class.coordinates[DELTA]))?,
            ));
        }
        8 => terms.push((Pair(POINT, POINT), class.coordinates[0])),
        _ => {}
    }
    Ok(terms)
}

pub(crate) fn backward(
    model: &mut Model,
    form: &Form,
    class: &Class,
    scale: i128,
) -> Result<State> {
    if scale == 0 {
        return Ok(State::zero());
    }
    let mut monomials = Vec::new();
    for (basis, coefficient) in inverse_terms(model, form, class)? {
        if scale % coefficient.denominator != 0 {
            return Err(Error::plain(
                "rvHilbSquareNakajima scale leaves nonintegral Fock coefficients",
            ));
        }
        if coefficient.numerator != 0 {
            let integer = checked(coefficient.multiply(Rational::integer(scale)))?.numerator;
            monomials.push((basis.factors(), integer));
        }
    }
    State::from_monomial_coefficients(model, monomials)
}

pub(crate) fn denominator(model: &mut Model, form: &Form, class: &Class) -> Result<i128> {
    let mut denominator = 1i128;
    for (_, coefficient) in inverse_terms(model, form, class)? {
        let common = gcd(denominator as u128, coefficient.denominator as u128) as i128;
        denominator = (denominator / common)
            .checked_mul(coefficient.denominator)
            .ok_or_else(|| {
                Error::plain("Hilbert-square Fock denominator exceeds the signed 128-bit range")
            })?;
    }
    Ok(denominator)
}
