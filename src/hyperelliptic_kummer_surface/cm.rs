use super::rational_sparse::{decimal, Integer, Magnitude, WORK_LIMBS};
use super::unimodular::{pair_coordinates, Vector, RANK};
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_ore::gcd;
use crate::syntax::Label;

type Coordinates = [Integer; RANK];

pub(crate) enum Command {
    Form(Label),
    Periods(Label),
    Equivalent(Label, Label),
    ClassNumber(Label),
}

pub(crate) struct Reduced {
    basis: [Coordinates; 2],
    a: Magnitude,
    b: Integer,
    c: Magnitude,
    discriminant: Magnitude,
}

impl Reduced {
    pub(super) fn new(
        model: &mut Model,
        u: &Vector,
        v: &Vector,
        index: Magnitude,
        discriminant: Magnitude,
    ) -> Result<Self> {
        let mut basis = saturate(model, u, v, index)?;
        let [a, b, c] = gram_form(model, &basis)?;
        let (mut a, mut b, mut c) = (a.magnitude, b, c.magnitude);
        loop {
            model.spend(WORK_LIMBS * (3 * RANK + 8))?;
            // Shear v by q*u so that -a < b+2aq <= a, with the positive boundary.
            let two_a = a.multiply_u128(2).ok_or_else(arithmetic)?;
            let q = floor_quotient(model, positive(a).add(b.negate())?, two_a)?;
            if !q.is_zero() {
                let aq = q.multiply_wide(positive(a))?;
                c = require_positive(positive(c).add(q.multiply_wide(b.add(aq)?)?)?)?;
                b = b.add(aq.multiply(2)?)?;
                let [u, v] = &mut basis;
                for (v, &u) in v.iter_mut().zip(u.iter()) {
                    *v = v.add(q.multiply_wide(u)?)?;
                }
            }
            // The swap (u,v) -> (v,-u) has determinant one.
            if a.cmp(c).is_gt() {
                swap(&mut basis, &mut a, &mut b, &mut c);
            } else {
                if a == c && b.negative && !b.is_zero() {
                    swap(&mut basis, &mut a, &mut b, &mut c);
                }
                break;
            }
        }
        if b.is_zero() {
            b.negative = false;
        }
        let checked = gram_form(model, &basis)?;
        let ac4 = a
            .checked_product(c)
            .and_then(|ac| ac.multiply_u128(4))
            .ok_or_else(arithmetic)?;
        let b2 = b
            .magnitude
            .checked_product(b.magnitude)
            .ok_or_else(arithmetic)?;
        if !ac4.cmp(b2).is_gt()
            || ac4.difference(b2) != discriminant
            || checked[0].magnitude != a
            || checked[1].magnitude != b.magnitude
            || (!b.is_zero() && checked[1].negative != b.negative)
            || checked[2].magnitude != c
        {
            return Err(Error::plain(
                "reduced K3 CM basis has an inconsistent Gram form",
            ));
        }
        Ok(Self {
            basis,
            a,
            b,
            c,
            discriminant,
        })
    }

    pub(crate) fn basis(&self, index: usize) -> Result<Vector> {
        let mut w = [0; 6];
        let mut b = [0; 16];
        for (target, &entry) in w.iter_mut().chain(&mut b).zip(&self.basis[index]) {
            *target = entry.signed().ok_or_else(|| {
                Error::plain("reduced K3 CM basis coordinate exceeds the signed i128 range")
            })?;
        }
        Vector::new(w, b)
    }

    pub(crate) fn equivalent(&self, other: &Self) -> bool {
        self.a == other.a
            && self.b.magnitude == other.b.magnitude
            && self.b.negative == other.b.negative
            && self.c == other.c
    }

    fn order_data(&self, model: &mut Model) -> Result<(Magnitude, Magnitude)> {
        let content = self.a.gcd(model, self.b.magnitude)?.gcd(model, self.c)?;
        let square = content.checked_product(content).ok_or_else(arithmetic)?;
        let order_discriminant = self
            .discriminant
            .divide_exact(model, square)?
            .ok_or_else(arithmetic)?;
        Ok((content, order_discriminant))
    }

    pub(crate) fn form_report(&self, model: &mut Model, name: &str) -> Result<String> {
        let (content, order_discriminant) = self.order_data(model)?;
        Ok(format!(
            "rvK3CMForm[{name}] = form=[{}, {}, {}], content={}, discriminant=-{}, order_discriminant=-{}",
            decimal(self.a), signed_text(self.b), decimal(self.c), decimal(content),
            decimal(self.discriminant), decimal(order_discriminant),
        ))
    }

    pub(crate) fn periods_report(&self, name: &str) -> Result<String> {
        let d = decimal(self.discriminant);
        let two_a = self.a.multiply_u128(2).ok_or_else(arithmetic)?;
        Ok(format!(
            "rvK3CMPeriods[{name}] = tau1=({} + sqrt(-{d}))/{}, tau2=({} + sqrt(-{d}))/2, branch=positive_imaginary",
            signed_text(self.b.negate()), decimal(two_a), signed_text(self.b),
        ))
    }

    pub(crate) fn class_number(&self, model: &mut Model) -> Result<usize> {
        let (_, discriminant) = self.order_data(model)?;
        proper_class_number(model, discriminant)
    }
}

fn proper_class_number(model: &mut Model, discriminant: Magnitude) -> Result<usize> {
    let d = discriminant.small().ok_or_else(class_number_limit)?;
    // Reduction implies |d| >= 3a^2. There are a candidate b of the
    // required parity in the half-open interval -a < b <= a.
    let max_a = (d / 3).isqrt();
    let candidates = max_a * (max_a + 1) / 2;
    let cost = usize::try_from(candidates)
        .ok()
        .and_then(|n| n.checked_mul(24))
        .ok_or_else(class_number_limit)?;
    model.spend(cost).map_err(|_| class_number_limit())?;
    let mut count = 0;
    for a in 1..=max_a {
        let mut b = 1 - a as i128;
        if (b & 1) as u128 != d & 1 {
            b += 1;
        }
        while b <= a as i128 {
            let square = b.unsigned_abs().pow(2);
            let numerator = square.checked_add(d).ok_or_else(class_number_limit)?;
            let denominator = 4 * a;
            if numerator % denominator == 0 {
                let c = numerator / denominator;
                if (a < c || (a == c && b >= 0)) && gcd(gcd(a, b.unsigned_abs()), c) == 1 {
                    count += 1;
                }
            }
            b += 2;
        }
    }
    Ok(count)
}

fn saturate(
    model: &mut Model,
    u: &Vector,
    v: &Vector,
    index: Magnitude,
) -> Result<[Coordinates; 2]> {
    let u_lattice = u.wide_basis_coordinates(model)?;
    let mut content = Magnitude::default();
    for entry in u_lattice {
        content = content.gcd(model, entry.resize()?.magnitude)?;
    }
    let quotient = index.divide_exact(model, content)?.ok_or_else(arithmetic)?;
    let projection = if quotient.small() == Some(1) {
        Magnitude::default()
    } else {
        let v_lattice = v.wide_basis_coordinates(model)?;
        // After splitting off primitive u, the second saturation factor is
        // the content of the projected v. Its residue suffices to divide in Γ.
        projection_residue(model, u_lattice, v_lattice, content, quotient)?
    };
    let mut saturated = [[Integer::from_i128(0); RANK]; 2];
    let ((uw, ub), (vw, vb)) = (u.coordinates(), v.coordinates());
    for (i, (&u, &v)) in uw.iter().chain(ub).zip(vw.iter().chain(vb)).enumerate() {
        let u = exact_divide(model, Integer::from_i128(u), content)?;
        let v = Integer::from_i128(v).add(positive(projection).multiply_wide(u)?.negate())?;
        saturated[0][i] = u;
        saturated[1][i] = exact_divide(model, v, quotient)?;
    }
    Ok(saturated)
}

fn projection_residue(
    model: &mut Model,
    u: [super::unimodular::Integer; RANK],
    v: [super::unimodular::Integer; RANK],
    content: Magnitude,
    modulus: Magnitude,
) -> Result<Magnitude> {
    // After incorporating a coordinate, gcd is the gcd of the primitive
    // u coordinates seen so far. residue is the same Bézout combination
    // applied to v, modulo the second saturation factor.
    let mut gcd = Magnitude::default();
    let mut residue = Magnitude::default();
    for (u, v) in u.into_iter().zip(v) {
        if u.is_zero() {
            continue;
        }
        let coordinate = exact_divide(model, u.resize()?, content)?;
        let (mut a, mut b) = (gcd, coordinate.magnitude);
        let (mut left, mut right) = (
            residue,
            residue_mod(
                model,
                if coordinate.negative { v.negate() } else { v }.resize()?,
                modulus,
            )?,
        );
        while !b.is_zero() {
            let (q, remainder) = a.div_rem(model, b)?;
            let next = positive(left).add(positive(q).multiply_wide(positive(right))?.negate())?;
            (a, b) = (b, remainder);
            (left, right) = (right, residue_mod(model, next, modulus)?);
        }
        gcd = a;
        residue = left;
        if gcd.small() == Some(1) {
            return Ok(residue);
        }
    }
    Err(Error::plain(
        "K3 CM saturation requires a primitive first vector",
    ))
}

fn residue_mod(model: &mut Model, value: Integer, modulus: Magnitude) -> Result<Magnitude> {
    let (_, remainder) = value.magnitude.div_rem(model, modulus)?;
    Ok(if value.negative && !remainder.is_zero() {
        modulus.difference(remainder)
    } else {
        remainder
    })
}

fn floor_quotient(model: &mut Model, value: Integer, divisor: Magnitude) -> Result<Integer> {
    let (mut magnitude, remainder) = value.magnitude.div_rem(model, divisor)?;
    if value.negative && !remainder.is_zero() {
        magnitude
            .add_multiple(Magnitude::from_small(1), 1)
            .ok_or_else(arithmetic)?;
    }
    Ok(Integer {
        magnitude,
        negative: value.negative,
    })
}

fn exact_divide(model: &mut Model, value: Integer, divisor: Magnitude) -> Result<Integer> {
    value
        .divide_exact(model, positive(divisor))?
        .ok_or_else(arithmetic)
}

fn swap(basis: &mut [Coordinates; 2], a: &mut Magnitude, b: &mut Integer, c: &mut Magnitude) {
    basis.swap(0, 1);
    basis[1] = basis[1].map(Integer::negate);
    std::mem::swap(a, c);
    *b = b.negate();
}

fn gram_form(model: &mut Model, [u, v]: &[Coordinates; 2]) -> Result<[Integer; 3]> {
    model.spend(WORK_LIMBS * 4 * RANK * 3)?;
    let a = require_positive(pair_coordinates(|i, j| u[i].multiply_wide(u[j]))?)?;
    let b = pair_coordinates(|i, j| u[i].multiply_wide(v[j]))?;
    let c = require_positive(pair_coordinates(|i, j| v[i].multiply_wide(v[j]))?)?;
    let (a, remainder_a) = a.divide_u32(2);
    let (c, remainder_c) = c.divide_u32(2);
    if remainder_a != 0 || remainder_c != 0 {
        return Err(Error::plain("K3 CM basis must have even square"));
    }
    Ok([positive(a), b, positive(c)])
}

fn positive(magnitude: Magnitude) -> Integer {
    Integer {
        magnitude,
        negative: false,
    }
}

fn require_positive(value: Integer) -> Result<Magnitude> {
    if value.negative || value.is_zero() {
        return Err(Error::plain("K3 CM reduction requires a positive form"));
    }
    Ok(value.magnitude)
}

fn signed_text(value: Integer) -> String {
    let magnitude = decimal(value.magnitude);
    if value.negative && !value.is_zero() {
        format!("-{magnitude}")
    } else {
        magnitude
    }
}

fn arithmetic() -> Error {
    Error::plain("K3 CM calculation exceeds the supported exact arithmetic range")
}

fn class_number_limit() -> Error {
    Error::plain("K3 CM class number enumeration exceeds the computation limit")
}
