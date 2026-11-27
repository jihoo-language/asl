use std::fmt::Write;
use std::rc::Rc;

use crate::cm_class_poly::Polynomial;
use crate::ec_poly::{self as poly, Field};
use crate::elliptic::integer_power_mod;
use crate::elliptic_end_order::{canonical::Lift, Certificate};
use crate::error::{Error, Result};
use crate::finite_field_factor;
use crate::hw_algebra::CoefficientRing;
use crate::hyperelliptic_kummer_surface::cm_ring::{Class, Order, GROUP_LIMIT};
use crate::model::{Model, MAX_SHOWN_BYTES};
use crate::odd_cobar::Prime;
use crate::rv_ore::gcd;
use crate::rv_unram::{mul_mod, Ring};
use crate::syntax::{Label, ParsedExpr};

pub(crate) struct Address(pub(crate) Label, pub(crate) usize, pub(crate) bool);

impl Address {
    pub(crate) fn ideal(&self, model: &mut Model) -> Result<(Order, PrimeIdeal)> {
        let order = model.cm_ring_order(&self.0)?;
        let p = prime(model, self.1)?;
        Ok((order, PrimeIdeal::new(model, order, p, self.2)?))
    }
}

pub(crate) enum Statement {
    Prime(Address),
    Orbit(Address, Box<ParsedExpr>),
    Factors(Label, Label),
    Lift(Label, Label),
}

pub(crate) enum Expression {
    Class(Address),
    Translate(Address, i128, Box<ParsedExpr>),
}

pub(crate) struct PrimeIdeal {
    pub(crate) class: Class,
    form: [i128; 3],
}

pub(crate) fn prime(model: &mut Model, p: usize) -> Result<Prime> {
    if p == 2 {
        return Ok(Prime::validated(2));
    }
    if p < 2 || p.is_multiple_of(2) {
        return Err(Error::plain("CM Artin characteristic must be prime"));
    }
    Prime::new(model, p)
}

impl PrimeIdeal {
    pub(crate) fn generator_residue(&self, order: Order) -> i128 {
        (self.form[1] + order.discriminant().rem_euclid(2)) / 2
    }

    pub(crate) fn new(
        model: &mut Model,
        order: Order,
        prime: Prime,
        conjugate: bool,
    ) -> Result<Self> {
        let d = order.discriminant();
        let p = i128::from(prime.get());
        let residue = d.rem_euclid(p) as u64;
        if residue == 0 {
            return Err(Error::plain(
                "CM Artin prime must not divide the order discriminant",
            ));
        }
        let b = if p == 2 {
            if d.rem_euclid(8) != 1 {
                return Err(Error::plain(
                    "CM Artin prime must split in the quadratic order",
                ));
            }
            1
        } else {
            let r = split_root(model, residue, p as u64)?;
            let r = i128::from(r.min(p as u64 - r));
            if r.rem_euclid(2) == d.rem_euclid(2) {
                r
            } else {
                r - p
            }
        };
        let b = if conjugate { -b } else { b };
        let numerator = d.unsigned_abs() + b.unsigned_abs().pow(2);
        let denominator = (4 * p) as u128;
        if !numerator.is_multiple_of(denominator) {
            return Err(Error::plain(
                "CM Artin prime form has inconsistent discriminant",
            ));
        }
        let form = [p, b, (numerator / denominator) as i128];
        Ok(Self {
            class: order.class(model, form)?,
            form,
        })
    }

    pub(crate) fn report(
        &self,
        model: &mut Model,
        order: Order,
        name: &str,
        p: usize,
        conjugate: bool,
    ) -> Result<String> {
        let orientation = usize::from(conjugate);
        let period = order.class_order(model, self.class)?;
        Ok(format!(
            "rvCMArtinPrime[{name},{p},{orientation}] = form={:?}, class={:?}, order={period}",
            self.form, self.class.0
        ))
    }

    pub(crate) fn orbit(
        &self,
        model: &mut Model,
        order: Order,
        name: &str,
        p: usize,
        conjugate: bool,
        start: Class,
    ) -> Result<String> {
        let period = order.class_order(model, self.class)?;
        if period > GROUP_LIMIT {
            return Err(Error::plain(
                "CM Artin orbit exceeds the finite-value limit",
            ));
        }
        let orientation = usize::from(conjugate);
        let mut output =
            format!("rvCMArtinOrbit[{name},{p},{orientation}] = period={period}, classes=[");
        let mut value = start;
        for i in 0..period {
            if i != 0 {
                output.push_str(", ");
            }
            write!(output, "{:?}", value.0).expect("writing an Artin orbit to a string");
            if output.len() >= MAX_SHOWN_BYTES {
                return Err(Error::plain("CM Artin orbit exceeds the display limit"));
            }
            value = order.compose(model, self.class, value)?;
        }
        if value != start {
            return Err(Error::plain(
                "CM Artin orbit failed to return to its starting class",
            ));
        }
        output.push(']');
        Ok(output)
    }
}

pub(crate) fn factor_report(
    model: &mut Model,
    polynomial: &Polynomial,
    ring: &Ring,
    polynomial_name: &str,
    field_name: &str,
) -> Result<String> {
    if ring.precision() != 1 {
        return Err(Error::plain(
            "CM Artin factorization requires a precision-one RVUnram field",
        ));
    }
    let order = Order::new(polynomial.discriminant)?;
    let p = ring.prime();
    let action = PrimeIdeal::new(model, order, Prime::validated(p as u32), false)?;
    let period = order.class_order(model, action.class)?;
    let residue_degree = ring.degree();
    let factor_degree = period / gcd(period as u128, residue_degree as u128) as usize;
    let h = polynomial.degree();
    if !h.is_multiple_of(factor_degree) {
        return Err(Error::plain(
            "CM class polynomial degree differs from the Artin cycle prediction",
        ));
    }
    let factor_count = h / factor_degree;

    let field = Field { ring };
    field.check_polynomial(h + 1)?;
    let coefficients = polynomial.unramified_coefficients(model, ring)?;
    let derivative = poly::derivative(model, field, &coefficients)?;
    let gcd = poly::monic_gcd(model, field, coefficients.clone(), derivative)?;
    if gcd != [field.one()] {
        return Err(Error::plain(
            "CM Artin class polynomial reduction is not squarefree",
        ));
    }
    let mut factors = finite_field_factor::factor_squarefree(model, field, coefficients.clone())?;
    if factors.len() != factor_count
        || factors
            .iter()
            .any(|factor| factor.len() != factor_degree + 1)
    {
        return Err(Error::plain(
            "CM Artin factorization disagrees with the predicted cycle lengths",
        ));
    }
    factors.sort();
    let mut product = vec![field.one()];
    for factor in &factors {
        product = poly::multiply(model, field, &product, factor)?;
    }
    if product != coefficients {
        return Err(Error::plain(
            "CM Artin factors failed to reconstruct the polynomial",
        ));
    }

    let mut output = format!(
        "rvCMArtinFactors[{polynomial_name},{field_name}] = p={p}, residue_degree={residue_degree}, artin_order={period}, factor_degree={factor_degree}, factor_count={factor_count}, factors=["
    );
    for (i, factor) in factors.iter().enumerate() {
        if i != 0 {
            output.push_str(", ");
        }
        write!(output, "{factor:?}").expect("writing an Artin factor to a string");
        if output.len() >= MAX_SHOWN_BYTES {
            return Err(Error::plain(
                "CM Artin factorization exceeds the display limit",
            ));
        }
    }
    output.push(']');
    Ok(output)
}

pub(crate) fn lift_report(
    model: &mut Model,
    certificate: Rc<Certificate>,
    ring: Rc<Ring>,
    certificate_name: &str,
    ring_name: &str,
) -> Result<String> {
    let lift = Lift::new(model, Rc::clone(&certificate), Rc::clone(&ring))?;
    let order = Order::new(certificate.discriminant())?;
    let action = PrimeIdeal::new(model, order, Prime::validated(ring.prime() as u32), false)?;
    let period = order.class_order(model, action.class)?;
    if !ring.degree().is_multiple_of(period) {
        return Err(Error::plain(
            "CM Artin period does not divide the canonical residue degree",
        ));
    }
    let j = lift.j_invariant(model)?;
    for degree in 1..=period {
        if (ring.frobenius(model, &j, degree)? == j) != (degree == period) {
            return Err(Error::plain(
                "canonical j Frobenius period differs from the CM Artin period",
            ));
        }
    }
    Ok(format!(
        "rvCMArtinLift[{certificate_name},{ring_name}] = p={}, residue_degree={}, precision={}, artin_order={period}, exact_period={period}, j={j:?}",
        ring.prime(), ring.degree(), ring.precision(),
    ))
}

fn split_root(model: &mut Model, residue: u64, p: u64) -> Result<u64> {
    let midpoint = (p - 1) / 2;
    if integer_power_mod(model, residue, midpoint, p)? != 1 {
        return Err(Error::plain(
            "CM Artin prime must split in the quadratic order",
        ));
    }
    if p % 4 == 3 {
        return integer_power_mod(model, residue, (p + 1) / 4, p);
    }

    // Tonelli–Shanks in the 2-primary subgroup of F_p^×.
    let mut exponent = (p - 1).trailing_zeros();
    let odd_part = (p - 1) >> exponent;
    let mut nonsquare = 2;
    while integer_power_mod(model, nonsquare, midpoint, p)? == 1 {
        nonsquare += 1;
    }
    let mut generator = integer_power_mod(model, nonsquare, odd_part, p)?;
    let mut root = integer_power_mod(model, residue, odd_part.div_ceil(2), p)?;
    let mut remainder = integer_power_mod(model, residue, odd_part, p)?;
    while remainder != 1 {
        let mut step = 0;
        let mut check = remainder;
        while check != 1 && step < exponent {
            model.spend(1)?;
            check = mul_mod(check, check, p);
            step += 1;
        }
        if step == exponent {
            return Err(Error::plain(
                "CM Artin modular square root failed verification",
            ));
        }
        let correction = integer_power_mod(model, generator, 1 << (exponent - step - 1), p)?;
        root = mul_mod(root, correction, p);
        generator = mul_mod(correction, correction, p);
        remainder = mul_mod(remainder, generator, p);
        exponent = step;
    }
    if mul_mod(root, root, p) != residue {
        return Err(Error::plain(
            "CM Artin modular square root failed verification",
        ));
    }
    Ok(root)
}
