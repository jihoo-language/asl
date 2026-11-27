use super::{field, CantorFactor, Divisor};
use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::model::Model;

pub(crate) fn tate(
    model: &mut Model,
    curve: &Curve,
    order: u64,
    torsion: &Divisor,
    divisor: &Divisor,
) -> Result<Scalar> {
    if order < 2 {
        return Err(Error::plain("rvJacTate requires an order of at least 2"));
    }
    let field = field(curve);
    let p = field.ring.prime();
    // Long division of q−1, whose base-p digits are all p−1.
    let mut exponent = vec![0; field.ring.degree()];
    model.spend(exponent.len())?;
    let mut remainder = 0u128;
    for digit in &mut exponent {
        let dividend = remainder * u128::from(p) + u128::from(p - 1);
        *digit = (dividend / u128::from(order)) as u64;
        remainder = dividend % u128::from(order);
    }
    if remainder != 0 {
        return Err(Error::plain("rvJacTate order must divide q−1"));
    }
    let value = leading_value(model, curve, order, torsion, divisor)?
        .ok_or_else(|| torsion_error("first", order))?;
    let one = field.one();
    let mut result = one.clone();
    if value != one {
        for digit in exponent {
            if result != one {
                result = field.ring.frobenius(model, &result, 1)?;
            }
            if digit != 0 {
                let factor = field.ring.power(model, &value, digit)?;
                result = field.multiply(model, &result, &factor)?;
            }
        }
    }
    Ok(result)
}

pub(crate) fn weil(
    model: &mut Model,
    curve: &Curve,
    order: u64,
    left: &Divisor,
    right: &Divisor,
) -> Result<Scalar> {
    if order < 2 {
        return Err(Error::plain("rvJacWeil requires an order of at least 2"));
    }
    let field = field(curve);
    if order.is_multiple_of(field.ring.prime()) {
        return Err(Error::plain(
            "rvJacWeil order must be coprime to the residue characteristic",
        ));
    }
    let forward = leading_value(model, curve, order, left, right)?
        .ok_or_else(|| torsion_error("first", order))?;
    let reverse = leading_value(model, curve, order, right, left)?
        .ok_or_else(|| torsion_error("second", order))?;
    let value = divide(model, field, forward, &reverse)?;
    // Tame reciprocity corrects the ratio by (−1) raised to n times the
    // intersection number of the two degree-zero divisors.
    if !order.is_multiple_of(2) && odd_intersection(model, field, left, right)? {
        field.ring.negate(model, &value)
    } else {
        Ok(value)
    }
}

fn torsion_error(argument: &str, order: u64) -> Error {
    Error::plain(format!(
        "{argument} argument must have order dividing {order}"
    ))
}

fn odd_intersection(
    model: &mut Model,
    field: Field<'_>,
    left: &Divisor,
    right: &Divisor,
) -> Result<bool> {
    // Both representatives contain infinity with coefficient −deg(u).
    let mut odd = left.u.len().is_multiple_of(2) && right.u.len().is_multiple_of(2);
    let difference = poly::sum(model, field, &left.v, &right.v, true)?;
    model.spend(right.u.len() * field.cells())?;
    let (equal_ordinates, _) = partition(model, field, right.u.clone(), &difference)?;
    visit_leading_terms(
        model,
        field,
        equal_ordinates,
        &left.u,
        |_, order, places, _| {
            // The Hasse order is the multiplicity in the left divisor. The
            // degree of places weights residue degree by right multiplicity.
            odd ^= !order.is_multiple_of(2) && places.len().is_multiple_of(2);
            Ok(())
        },
    )?;
    Ok(odd)
}

/// Evaluate the leading unit of a function with divisor `order * torsion` on
/// the degree-zero divisor represented by `divisor`. The uniformizers are
/// x−α at ordinary points, y at branch points, and x^g/y at infinity.
/// Return `None` when `order` does not kill `torsion`.
pub(crate) fn leading_value(
    model: &mut Model,
    curve: &Curve,
    order: u64,
    torsion: &Divisor,
    divisor: &Divisor,
) -> Result<Option<Scalar>> {
    if order == 0 {
        return Err(Error::plain("Miller order must be positive"));
    }
    let field = field(curve);
    if torsion.is_zero() {
        return Ok(Some(field.one()));
    }
    if divisor.is_zero() {
        let multiple = torsion.scale(model, curve, i128::from(order))?;
        return Ok(multiple.is_zero().then(|| field.one()));
    }
    let support = Support::new(model, curve, divisor)?;
    model.spend(torsion.cells(curve))?;
    let mut multiple = torsion.clone();
    let mut value = field.one();
    for bit in (0..order.ilog2()).rev() {
        value = field.multiply(model, &value, &value)?;
        multiple = support.add(model, curve, &multiple, &multiple, &mut value)?;
        if (order >> bit) & 1 != 0 {
            multiple = support.add(model, curve, &multiple, torsion, &mut value)?;
        }
    }
    Ok(multiple.is_zero().then_some(value))
}

struct Support<'a> {
    divisor: &'a Divisor,
    field: Field<'a>,
    genus: usize,
    branch: Polynomial,
    ordinary: Polynomial,
    derivative: Polynomial,
}

impl<'a> Support<'a> {
    fn new(model: &mut Model, curve: &'a Curve, divisor: &'a Divisor) -> Result<Self> {
        let field = field(curve);
        model.spend(divisor.u.len() * field.cells())?;
        let (branch, ordinary) = partition(model, field, divisor.u.clone(), curve.polynomial())?;
        let derivative = if branch.len() > 1 {
            poly::derivative(model, field, curve.polynomial())?
        } else {
            Vec::new()
        };
        Ok(Self {
            divisor,
            field,
            genus: curve.grade.genus,
            branch,
            ordinary,
            derivative,
        })
    }

    fn add(
        &self,
        model: &mut Model,
        curve: &Curve,
        left: &Divisor,
        right: &Divisor,
        value: &mut Scalar,
    ) -> Result<Divisor> {
        left.add_with(model, curve, right, |model, factor| {
            let unit = match factor {
                CantorFactor::Composition(gcd) => self.polynomial_unit(model, gcd)?,
                CantorFactor::Reduction {
                    interpolation,
                    complement,
                    norm,
                } => {
                    let numerator = self.interpolation_unit(model, interpolation, norm)?;
                    let denominator = self.polynomial_unit(model, complement)?;
                    divide(model, self.field, numerator, &denominator)?
                }
            };
            *value = self.field.multiply(model, value, &unit)?;
            Ok(())
        })
    }

    fn polynomial_unit(&self, model: &mut Model, polynomial: &[Scalar]) -> Result<Scalar> {
        if polynomial.len() == 1 {
            return Ok(self.field.one());
        }
        // All polynomial factors emitted by Cantor are monic at infinity.
        model.spend((self.ordinary.len() + self.branch.len()) * self.field.cells())?;
        let ordinary =
            finite_polynomial_unit(model, self.field, self.ordinary.clone(), polynomial, None)?;
        let branch = finite_polynomial_unit(
            model,
            self.field,
            self.branch.clone(),
            polynomial,
            Some(&self.derivative),
        )?;
        self.field.multiply(model, &ordinary, &branch)
    }

    fn interpolation_unit(
        &self,
        model: &mut Model,
        interpolation: &[Scalar],
        norm: &[Scalar],
    ) -> Result<Scalar> {
        let field = self.field;
        let difference = poly::sum(model, field, &self.divisor.v, interpolation, true)?;
        model.spend(self.divisor.u.len() * field.cells())?;
        let (equal, unequal) = partition(model, field, self.divisor.u.clone(), &difference)?;
        let mut result = norm_on(model, field, unequal, &difference)?;
        if equal.len() > 1 {
            // Where y=s=0, the leading unit of y−s in parameter y is 1.
            let (_, ordinary) = partition(model, field, equal, interpolation)?;
            if ordinary.len() > 1 {
                // At y=s≠0 use (y−s)(y+s)=f−s² with unit y+s=2s.
                let double = poly::scale(model, field, interpolation, &field.signed(2))?;
                model.spend(ordinary.len() * field.cells())?;
                let denominator = norm_on(model, field, ordinary.clone(), &double)?;
                let numerator = finite_polynomial_unit(model, field, ordinary, norm, None)?;
                result = field.multiply(model, &result, &numerator)?;
                result = divide(model, field, result, &denominator)?;
            }
        }
        if interpolation.len() > self.genus + 1 {
            let leading = field.ring.negate(model, interpolation.last().unwrap())?;
            let infinity = field
                .ring
                .power(model, &leading, (self.divisor.u.len() - 1) as u64)?;
            result = divide(model, field, result, &infinity)?;
        }
        Ok(result)
    }
}

// Partition by complete primary factors, retaining the multiplicities of the
// original Mumford support. Raising value to at least deg(support) saturates
// its vanishing multiplicity at every common root.
fn partition(
    model: &mut Model,
    field: Field<'_>,
    support: Polynomial,
    value: &[Scalar],
) -> Result<(Polynomial, Polynomial)> {
    let one = || vec![field.one()];
    if value.is_empty() {
        return Ok((support, one()));
    }
    if support.len() == 1 || value.len() == 1 {
        return Ok((one(), support));
    }
    model.spend((support.len() + value.len()) * field.cells())?;
    let power = poly::power_mod(model, field, value, (support.len() - 1) as u64, &support)?;
    let mut vanishing = poly::gcd(model, field, support.clone(), power)?;
    if vanishing.len() == support.len() {
        return Ok((support, one()));
    }
    if vanishing.len() == 1 {
        return Ok((one(), support));
    }
    let leading = vanishing.last().expect("nonzero gcd with monic support");
    if leading != &field.one() {
        let inverse = field.inverse_unit(model, leading)?;
        vanishing = poly::scale(model, field, &vanishing, &inverse)?;
    }
    let remaining = poly::exact_divide(model, field, &support, &vanishing)?;
    Ok((vanishing, remaining))
}

fn norm_on(
    model: &mut Model,
    field: Field<'_>,
    support: Polynomial,
    polynomial: &[Scalar],
) -> Result<Scalar> {
    if support.len() == 1 {
        return Ok(field.one());
    }
    poly::resultant(model, field, support, polynomial)
}

fn finite_polynomial_unit(
    model: &mut Model,
    field: Field<'_>,
    pending: Polynomial,
    polynomial: &[Scalar],
    branch_derivative: Option<&[Scalar]>,
) -> Result<Scalar> {
    let mut result = field.one();
    visit_leading_terms(
        model,
        field,
        pending,
        polynomial,
        |model, order, places, coefficient| {
            let denominator = if let Some(derivative) = branch_derivative.filter(|_| order > 0) {
                // In parameter y, (x−α)/y² has leading unit 1/f′(α).
                model.spend(places.len() * field.cells())?;
                let norm = norm_on(model, field, places.clone(), derivative)?;
                field.ring.power(model, &norm, order as u64)?
            } else {
                field.one()
            };
            let unit = norm_on(model, field, places, coefficient)?;
            let unit = divide(model, field, unit, &denominator)?;
            result = field.multiply(model, &result, &unit)?;
            Ok(())
        },
    )?;
    Ok(result)
}

fn visit_leading_terms(
    model: &mut Model,
    field: Field<'_>,
    mut pending: Polynomial,
    polynomial: &[Scalar],
    mut visit: impl FnMut(&mut Model, usize, Polynomial, &[Scalar]) -> Result<()>,
) -> Result<()> {
    if pending.len() == 1 {
        return Ok(());
    }
    let mut binomial = Vec::new();
    for order in 0..polynomial.len() {
        let hasse;
        let coefficient = if order == 0 {
            polynomial
        } else {
            if order == 1 {
                model.spend(polynomial.len())?;
                binomial = vec![1; polynomial.len()];
            }
            // Update the Pascal column from binom(i,k−1) to binom(i,k).
            model.spend(polynomial.len() - order + 1)?;
            let mut sum = 0;
            for choice in &mut binomial[order - 1..] {
                let next = (sum + *choice) % field.ring.prime();
                *choice = sum;
                sum = next;
            }
            hasse = poly::trim(
                polynomial
                    .iter()
                    .zip(&binomial)
                    .skip(order)
                    .map(|(c, &choice)| field.ring.scale(model, c, choice))
                    .collect::<Result<_>>()?,
            );
            &hasse
        };
        let (vanishing, nonvanishing) = partition(model, field, pending, coefficient)?;
        if nonvanishing.len() > 1 {
            visit(model, order, nonvanishing, coefficient)?;
        }
        if vanishing.len() == 1 {
            return Ok(());
        }
        pending = vanishing;
    }
    Err(Error::plain("Miller factor has no nonzero leading unit"))
}

fn divide(model: &mut Model, field: Field<'_>, value: Scalar, by: &Scalar) -> Result<Scalar> {
    if by == &field.one() {
        return Ok(value);
    }
    let inverse = field.inverse_unit(model, by)?;
    field.multiply(model, &value, &inverse)
}
