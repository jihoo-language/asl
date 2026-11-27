use super::{arithmetic, inconsistent, Matrix, Wide};
use crate::error::Result;
use crate::model::Model;
use crate::rv_qp_quad::{self, Element};
use crate::wide_integer::Natural;

struct Decomposition {
    lower: [[Element; 4]; 4],
    diagonal: [Element; 4],
}

pub(super) fn each_norm_two(
    model: &mut Model,
    gram: Matrix,
    mut visit: impl FnMut(&mut Model, [i128; 4]) -> Result<()>,
) -> Result<()> {
    let decomposition = Decomposition::new(model, gram)?;
    decomposition.search(model, &mut [0; 4], 4, number(2)?, &mut visit)
}

impl Decomposition {
    fn new(model: &mut Model, gram: Matrix) -> Result<Self> {
        let mut lower = [[Element::ZERO; 4]; 4];
        let mut diagonal = [Element::ZERO; 4];
        for j in 0..4 {
            lower[j][j] = Element::ONE;
            let mut pivot = number(gram[j][j])?;
            for (k, &d) in diagonal.iter().enumerate().take(j) {
                let square = multiply(model, lower[j][k], lower[j][k])?;
                let change = multiply(model, d, square)?;
                pivot = rv_qp_quad::affine_combination(model, pivot, change, true)?;
            }
            if pivot.a <= 0 {
                return Err(inconsistent());
            }
            diagonal[j] = pivot;
            let inverse = pivot.rational_inverse(model)?;
            for i in j + 1..4 {
                let mut entry = number(gram[i][j])?;
                for (k, &d) in diagonal.iter().enumerate().take(j) {
                    let product = multiply(model, lower[i][k], lower[j][k])?;
                    let change = multiply(model, d, product)?;
                    entry = rv_qp_quad::affine_combination(model, entry, change, true)?;
                }
                lower[i][j] = multiply(model, entry, inverse)?;
            }
        }
        Ok(Self { lower, diagonal })
    }

    fn search(
        &self,
        model: &mut Model,
        coefficients: &mut [i128; 4],
        variables: usize,
        remaining: Element,
        visit: &mut impl FnMut(&mut Model, [i128; 4]) -> Result<()>,
    ) -> Result<()> {
        let Some(j) = variables.checked_sub(1) else {
            if remaining.is_zero() {
                visit(model, *coefficients)?;
            }
            return Ok(());
        };
        let mut bias = Element::ZERO;
        for (k, &x) in coefficients.iter().enumerate().skip(j + 1) {
            let term = multiply(model, self.lower[k][j], number(x)?)?;
            bias = rv_qp_quad::affine_combination(model, bias, term, false)?;
        }
        let inverse = self.diagonal[j].rational_inverse(model)?;
        let radius_squared = multiply(model, remaining, inverse)?;
        let (first, last) = interval(model, bias, radius_squared)?;
        for x in first..=last {
            model.spend(1)?;
            coefficients[j] = x;
            let coordinate = rv_qp_quad::affine_combination(model, number(x)?, bias, false)?;
            let square = multiply(model, coordinate, coordinate)?;
            let contribution = multiply(model, self.diagonal[j], square)?;
            let room = rv_qp_quad::affine_combination(model, remaining, contribution, true)?;
            if room.a < 0 {
                return Err(inconsistent());
            }
            self.search(model, coefficients, j, room, visit)?;
        }
        Ok(())
    }
}

// For bias p/q and squared radius u/v, the integer n=q*x+p must satisfy
// |n| <= floor(sqrt(floor(u*q^2/v))). Convert these exact bounds back to x.
fn interval(model: &mut Model, bias: Element, squared: Element) -> Result<(i128, i128)> {
    if squared.a < 0 {
        return Err(inconsistent());
    }
    let numerator = Wide::<16>::product(bias.c, bias.c).multiply(squared.a)?;
    let (bound, _) = numerator.magnitude.divide(model, squared.c as u128)?;
    let radius = bound.integer_sqrt(model)?;
    let center = Wide::from_i128(bias.a).negate();
    let signed_radius = Wide {
        magnitude: radius,
        negative: false,
    };
    let lower = center.add(signed_radius.negate())?;
    let upper = center.add(signed_radius)?;
    Ok((
        divide_bound(model, lower, bias.c as u128, true)?,
        divide_bound(model, upper, bias.c as u128, false)?,
    ))
}

fn divide_bound(model: &mut Model, x: Wide<16>, d: u128, ceiling: bool) -> Result<i128> {
    let (mut quotient, remainder) = x.magnitude.divide(model, d)?;
    if remainder != 0 && ceiling != x.negative {
        quotient
            .add_multiple(Natural::from_small(1), 1)
            .ok_or_else(arithmetic)?;
    }
    quotient.signed(x.negative).ok_or_else(arithmetic)
}

fn number(x: i128) -> Result<Element> {
    Element::new(x, 0, 1)
}

fn multiply(model: &mut Model, a: Element, b: Element) -> Result<Element> {
    rv_qp_quad::product(model, 0, a, b)
}

#[cfg(test)]
mod tests {
    use super::{each_norm_two, interval, Element, Model};

    #[test]
    fn exact_shell_search_handles_a_skew_integral_basis() {
        let mut model = Model::new();
        let gram = [
            [10001, 100, 0, 0],
            [100, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ];
        let mut result = Vec::new();
        each_norm_two(&mut model, gram, |_, [a, b, c, d]| {
            result.push([a, a * 100 + b, c, d]);
            Ok(())
        })
        .unwrap_or_else(|err| panic!("{err}"));
        result.sort_unstable();
        result.dedup();
        assert_eq!(result.len(), 24);
        assert!(result
            .iter()
            .all(|x| x.iter().map(|a| a * a).sum::<i128>() == 2));
    }

    #[test]
    fn rational_search_interval_includes_exact_boundary() {
        let mut model = Model::new();
        let bias = Element::new(2, 0, 3).unwrap_or_else(|err| panic!("{err}"));
        let radius_squared = Element::new(1, 0, 9).unwrap_or_else(|err| panic!("{err}"));
        let bounds =
            interval(&mut model, bias, radius_squared).unwrap_or_else(|err| panic!("{err}"));
        assert_eq!(bounds, (-1, -1));
    }
}
