use crate::error::{Error, Result};
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};

use super::eigenpacket::{Element, Packet};
use super::hecke::{integer, Operator};
use super::residual_hecke::power_mod;
use super::{decimal, multiply, Integer, Magnitude, INTEGER_WORDS};

pub(crate) enum Query {
    Data(Label),
    Fourier(Label, u128),
    Local(Label, u64),
    Integral(Box<ParsedExpr>),
}

pub(crate) struct Order {
    field_discriminant: u128,
    generator_discriminant: Magnitude,
    generator_index: Magnitude,
    generator_offset: Integer,
    conductor: Magnitude,
    order_discriminant: Magnitude,
    second_coefficient: [Integer; 2],
}

impl Order {
    pub(super) fn new(model: &mut Model, packet: &Packet) -> Result<Self> {
        if packet.degree() != 2 {
            return Err(Error::plain(
                "quadratic Hecke order requires a packet of degree two",
            ));
        }
        let generator_discriminant = discriminant(model, packet.polynomial())?;
        let root = generator_discriminant.integer_sqrt(model)?;
        if root.checked_product(root) == Some(generator_discriminant) {
            return Err(not_quadratic_field());
        }

        // The integral cusp basis has initial coefficients (1, *) and
        // (0, 1). Every a_n is thus an integral combination of 1 and a_2.
        // Factor the intrinsic a_2 discriminant even for a scaled generator.
        let t2 = Operator::basis(2)?;
        let second_polynomial = t2.characteristic(model, packet.weight, 2)?;
        let second_discriminant = discriminant(model, &second_polynomial)?;
        let field_discriminant = square_class(
            model,
            second_discriminant
                .small()
                .ok_or_else(factorization_limit)?,
        )?;
        let disc = Magnitude::from_small(field_discriminant);
        let square_index = generator_discriminant
            .divide_exact(model, disc)?
            .ok_or_else(inconsistent)?;
        let generator_index = square_index.integer_sqrt(model)?;
        if generator_index.is_zero()
            || generator_index.checked_product(generator_index) != Some(square_index)
        {
            return Err(inconsistent());
        }

        let trace = packet.polynomial()[1].negate();
        let offset_numerator = if field_discriminant.is_multiple_of(2) {
            trace
        } else {
            trace
                .add(positive(generator_index).negate())
                .map_err(|_| arithmetic())?
        };
        let generator_offset = exact_divide(model, offset_numerator, integer(2))?;
        let mut order = Self {
            field_discriminant,
            generator_discriminant,
            generator_index,
            generator_offset,
            conductor: Magnitude::ZERO,
            order_discriminant: second_discriminant,
            second_coefficient: [integer(0); 2],
        };
        let second = packet.image_element(model, &t2)?;
        order.second_coefficient = order.coordinates(model, &second)?;
        order.conductor = order.second_coefficient[1].magnitude;
        let square = order
            .conductor
            .checked_product(order.conductor)
            .ok_or_else(arithmetic)?;
        if order.conductor.is_zero() || square.checked_product(disc) != Some(second_discriminant) {
            return Err(inconsistent());
        }
        model.allocate(1 + 8 * INTEGER_WORDS)?;
        Ok(order)
    }

    pub(crate) fn data(&self, name: &str) -> String {
        let d = self.field_discriminant;
        let parity = d % 2;
        let offset = decimal(self.generator_offset);
        let u = decimal(positive(self.generator_index));
        let g = decimal(positive(self.conductor));
        let generator_disc = decimal(positive(self.generator_discriminant));
        let order_disc = decimal(positive(self.order_discriminant));
        format!(
            "rvHeckeQuadraticOrder[{name}] = D_K={d}, omega=({parity}+sqrt({d}))/2, theta=({offset}, {u}), theta_index={u}, theta_discriminant={generator_disc}, conductor={g}, order_discriminant={order_disc}, congruence_ideal=({g}*sqrt({d})), congruence_norm={order_disc}"
        )
    }

    pub(crate) fn fourier(
        &self,
        model: &mut Model,
        packet: &Packet,
        name: &str,
        n: u128,
    ) -> Result<String> {
        let [a, b] = match n {
            0 => [integer(0); 2],
            1 => [integer(1), integer(0)],
            2 => self.second_coefficient,
            _ => {
                let value = packet.image_element(model, &Operator::basis(n)?)?;
                let coordinates = self.coordinates(model, &value)?;
                if coordinates[1]
                    .magnitude
                    .divide_exact(model, self.conductor)?
                    .is_none()
                {
                    return Err(inconsistent());
                }
                coordinates
            }
        };
        Ok(format!(
            "rvHeckeQuadraticCoeff[{name}, {n}] = ({}, {})",
            decimal(a),
            decimal(b)
        ))
    }

    pub(crate) fn integral(&self, model: &mut Model, value: &Element) -> Result<String> {
        let [a, b] = self.coordinates(model, value)?;
        Ok(format!(
            "rvHeckeQuadraticIntegral = ({}, {})",
            decimal(a),
            decimal(b)
        ))
    }

    pub(crate) fn local(&self, model: &mut Model, name: &str, p: u64) -> Result<String> {
        check_local_prime(model, p)?;
        let d = self.field_discriminant;
        let conductor_valuation = valuation(model, self.conductor, p)?;
        let (splitting, e, f, count) = match kronecker(model, d, p)? {
            0 => ("ramified", 2, 1, 1),
            1 => ("split", 1, 1, 2),
            _ => ("inert", 1, 2, 1),
        };
        let discriminant_valuation = if e == 2 {
            valuation(model, Magnitude::from_small(d), p)?
        } else {
            0
        };
        let depth = e * conductor_valuation + discriminant_valuation;
        Ok(format!(
            "rvHeckeQuadraticLocal[{name}, {p}] = {splitting}, e={e}, f={f}, primes={count}, depth={depth}, first_sharp_index=2"
        ))
    }

    fn coordinates(&self, model: &mut Model, value: &Element) -> Result<[Integer; 2]> {
        let (numerator, denominator) = value.fraction_parts();
        let [a, b] = *numerator else {
            return Err(inconsistent());
        };
        let beta = multiply(model, b, positive(self.generator_index))?;
        let alpha = multiply(model, b, self.generator_offset)?
            .add(a)
            .map_err(|_| arithmetic())?;
        let denominator = positive(denominator);
        Ok([
            exact_divide(model, alpha, denominator)?,
            exact_divide(model, beta, denominator)?,
        ])
    }
}

fn discriminant(model: &mut Model, polynomial: &[Integer]) -> Result<Magnitude> {
    let trace = polynomial[1];
    let square = multiply(model, trace, trace)?;
    let four_norm = multiply(model, polynomial[2], integer(4))?;
    let difference = square.add(four_norm.negate()).map_err(|_| arithmetic())?;
    if difference.negative || difference.is_zero() {
        return Err(not_quadratic_field());
    }
    Ok(difference.magnitude)
}

fn square_class(model: &mut Model, mut n: u128) -> Result<u128> {
    let mut radicand = 1u128;
    let mut divisor = 2;
    while divisor <= n / divisor {
        model.spend(1).map_err(|_| factorization_limit())?;
        let mut odd = false;
        while n.is_multiple_of(divisor) {
            model.spend(1).map_err(|_| factorization_limit())?;
            odd = !odd;
            n /= divisor;
        }
        if odd {
            radicand *= divisor;
        }
        divisor += if divisor == 2 { 1 } else { 2 };
    }
    // The residual factor is 1 or prime by the completed trial bound.
    radicand *= n;
    if radicand == 1 {
        return Err(not_quadratic_field());
    }
    if radicand % 4 == 1 {
        Ok(radicand)
    } else {
        radicand.checked_mul(4).ok_or_else(factorization_limit)
    }
}

fn check_local_prime(model: &mut Model, p: u64) -> Result<()> {
    let not_prime = || Error::plain(format!("quadratic Hecke local index {p} is not prime"));
    if p == 2 {
        return Ok(());
    }
    if p < 2 || p.is_multiple_of(2) {
        return Err(not_prime());
    }
    let mut divisor = 3;
    while divisor <= p / divisor {
        model.spend(1)?;
        if p.is_multiple_of(divisor) {
            return Err(not_prime());
        }
        divisor += 2;
    }
    Ok(())
}

fn kronecker(model: &mut Model, d: u128, p: u64) -> Result<i8> {
    let residue = (d % u128::from(p)) as u64;
    if residue == 0 {
        return Ok(0);
    }
    if p == 2 {
        return Ok(if d % 8 == 1 { 1 } else { -1 });
    }
    let symbol = power_mod(model, residue, u128::from((p - 1) / 2), p)?;
    Ok(if symbol == 1 { 1 } else { -1 })
}

fn valuation(model: &mut Model, mut n: Magnitude, p: u64) -> Result<u32> {
    let mut value = 0;
    loop {
        model.spend(1)?;
        let (quotient, residue) = n.divide(model, u128::from(p))?;
        if residue != 0 {
            return Ok(value);
        }
        n = quotient;
        value += 1;
    }
}

fn positive(magnitude: Magnitude) -> Integer {
    Integer {
        magnitude,
        negative: false,
    }
}

fn exact_divide(model: &mut Model, numerator: Integer, divisor: Integer) -> Result<Integer> {
    numerator.divide_exact(model, divisor)?.ok_or_else(|| {
        Error::plain("Hecke eigenpacket element is nonintegral in the quadratic maximal order")
    })
}

fn not_quadratic_field() -> Error {
    Error::plain("quadratic Hecke order requires an irreducible real quadratic packet")
}

fn factorization_limit() -> Error {
    Error::plain("quadratic Hecke discriminant factorization exceeds the exact computation range")
}

fn arithmetic() -> Error {
    Error::plain("quadratic Hecke exact intermediate exceeds the 512-bit range")
}

fn inconsistent() -> Error {
    Error::plain("quadratic Hecke discriminant or coefficient order identities are inconsistent")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn square_class_and_dyadic_ramification_are_exact() -> std::result::Result<(), String> {
        let model = &mut Model::new();
        assert_eq!(
            square_class(model, 200).map_err(|error| error.to_string())?,
            8
        );
        assert_eq!(
            square_class(model, 48).map_err(|error| error.to_string())?,
            12
        );
        assert_eq!(
            square_class(model, 45).map_err(|error| error.to_string())?,
            5
        );
        assert!(square_class(model, 36).is_err());

        // Exercise the dyadic formula in Q(sqrt(2)), since the level-one
        // dimension-two weights in the integration tests have odd D_K.
        let order = Order {
            field_discriminant: 8,
            generator_discriminant: Magnitude::from_small(128),
            generator_index: Magnitude::from_small(4),
            generator_offset: integer(1),
            conductor: Magnitude::from_small(4),
            order_discriminant: Magnitude::from_small(128),
            second_coefficient: [integer(1), integer(4)],
        };
        assert_eq!(
            order.local(model, "Test", 2).map_err(|error| error.to_string())?,
            "rvHeckeQuadraticLocal[Test, 2] = ramified, e=2, f=1, primes=1, depth=7, first_sharp_index=2"
        );
        Ok(())
    }
}
