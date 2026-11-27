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

pub(super) struct Form {
    pub(super) a: Magnitude,
    pub(super) b: Integer,
    pub(super) c: Magnitude,
}

pub(super) enum BasisChange {
    Shear(Integer),
    Rotate,
}

impl Form {
    pub(super) fn reduce(
        &mut self,
        model: &mut Model,
        mut transform: impl FnMut(&mut Model, BasisChange) -> Result<()>,
    ) -> Result<()> {
        let Self { a, b, c } = self;
        loop {
            model.spend(WORK_LIMBS * 8)?;
            // Shear v by q*u so that -a < b+2aq <= a, with the positive boundary.
            let two_a = a.multiply_u128(2).ok_or_else(arithmetic)?;
            let q = floor_quotient(model, positive(*a).add(b.negate())?, two_a)?;
            if !q.is_zero() {
                let aq = q.multiply_wide(positive(*a))?;
                *c = require_positive(positive(*c).add(q.multiply_wide(b.add(aq)?)?)?)?;
                *b = b.add(aq.multiply(2)?)?;
                transform(model, BasisChange::Shear(q))?;
            }
            let rotate = a.cmp(*c).is_gt();
            // The rotation (u,v) -> (v,-u) has determinant one.
            if rotate || (*a == *c && b.negative && !b.is_zero()) {
                std::mem::swap(a, c);
                *b = b.negate();
                transform(model, BasisChange::Rotate)?;
            }
            if !rotate {
                break;
            }
        }
        if b.is_zero() {
            b.negative = false;
        }
        Ok(())
    }

    pub(super) fn discriminant(&self) -> Result<Option<Magnitude>> {
        let ac4 = self
            .a
            .checked_product(self.c)
            .and_then(|ac| ac.multiply_u128(4))
            .ok_or_else(arithmetic)?;
        let b2 = self
            .b
            .magnitude
            .checked_product(self.b.magnitude)
            .ok_or_else(arithmetic)?;
        Ok(ac4.cmp(b2).is_gt().then(|| ac4.difference(b2)))
    }

    pub(super) fn content(&self, model: &mut Model) -> Result<Magnitude> {
        self.a.gcd(model, self.b.magnitude)?.gcd(model, self.c)
    }

    fn equivalent(&self, other: &Self) -> bool {
        self.a == other.a
            && self.b.magnitude == other.b.magnitude
            && (self.b.is_zero() || self.b.negative == other.b.negative)
            && self.c == other.c
    }
}

pub(crate) fn reduce_binary(
    model: &mut Model,
    [a, b, c]: [i128; 3],
) -> Result<([i128; 3], [[i128; 2]; 2])> {
    let mut form = Form {
        a: Magnitude::from_small(u128::try_from(a).map_err(|_| arithmetic())?),
        b: Integer::from_i128(b),
        c: Magnitude::from_small(u128::try_from(c).map_err(|_| arithmetic())?),
    };
    if a <= 0 || c <= 0 || form.discriminant()?.is_none() {
        return Err(Error::plain(
            "Gauss reduction needs a positive binary quadratic form",
        ));
    }
    let zero = Integer::from_i128(0);
    let one = Integer::from_i128(1);
    let mut columns = [[one, zero], [zero, one]];
    form.reduce(model, |_, change| {
        match change {
            BasisChange::Shear(q) => {
                let [u, v] = &mut columns;
                for (v, &u) in v.iter_mut().zip(u.iter()) {
                    *v = v.add(q.multiply_wide(u)?)?;
                }
            }
            BasisChange::Rotate => {
                columns.swap(0, 1);
                columns[1] = columns[1].map(Integer::negate);
            }
        }
        Ok(())
    })?;
    let a = positive(form.a).signed().ok_or_else(arithmetic)?;
    let b = form.b.signed().ok_or_else(arithmetic)?;
    let c = positive(form.c).signed().ok_or_else(arithmetic)?;
    let mut basis = [[0; 2]; 2];
    for (row, coefficients) in basis.iter_mut().zip(columns) {
        for (x, c) in row.iter_mut().zip(coefficients) {
            *x = c.signed().ok_or_else(arithmetic)?;
        }
    }
    Ok(([a, b, c], basis))
}

pub(crate) struct Reduced {
    basis: [Coordinates; 2],
    form: Form,
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
        let mut form = gram_form(model, &basis)?;
        form.reduce(model, |model, operation| {
            match operation {
                BasisChange::Shear(q) => {
                    model.spend(WORK_LIMBS * 3 * RANK)?;
                    let [u, v] = &mut basis;
                    for (v, &u) in v.iter_mut().zip(u.iter()) {
                        *v = v.add(q.multiply_wide(u)?)?;
                    }
                }
                BasisChange::Rotate => {
                    basis.swap(0, 1);
                    basis[1] = basis[1].map(Integer::negate);
                }
            }
            Ok(())
        })?;
        if form.discriminant()? != Some(discriminant)
            || !form.equivalent(&gram_form(model, &basis)?)
        {
            return Err(Error::plain(
                "reduced K3 CM basis has an inconsistent Gram form",
            ));
        }
        Ok(Self {
            basis,
            form,
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
        self.form.equivalent(&other.form)
    }

    fn order_data(&self, model: &mut Model) -> Result<(Magnitude, Magnitude)> {
        let content = self.form.content(model)?;
        let square = content.checked_product(content).ok_or_else(arithmetic)?;
        let order_discriminant = self
            .discriminant
            .divide_exact(model, square)?
            .ok_or_else(arithmetic)?;
        Ok((content, order_discriminant))
    }

    pub(super) fn primitive_form(
        &self,
        model: &mut Model,
        expected_discriminant: Magnitude,
    ) -> Result<Form> {
        let (content, discriminant) = self.order_data(model)?;
        if discriminant != expected_discriminant {
            return Err(Error::plain(
                "K3 CM primitive form discriminant differs from the named order",
            ));
        }
        Ok(Form {
            a: self
                .form
                .a
                .divide_exact(model, content)?
                .ok_or_else(arithmetic)?,
            b: exact_divide(model, self.form.b, content)?,
            c: self
                .form
                .c
                .divide_exact(model, content)?
                .ok_or_else(arithmetic)?,
        })
    }

    pub(crate) fn form_report(&self, model: &mut Model, name: &str) -> Result<String> {
        let (content, order_discriminant) = self.order_data(model)?;
        Ok(format!(
            "rvK3CMForm[{name}] = form=[{}, {}, {}], content={}, discriminant=-{}, order_discriminant=-{}",
            decimal(self.form.a), signed_text(self.form.b), decimal(self.form.c), decimal(content),
            decimal(self.discriminant), decimal(order_discriminant),
        ))
    }

    pub(crate) fn periods_report(&self, name: &str) -> Result<String> {
        let d = decimal(self.discriminant);
        let two_a = self.form.a.multiply_u128(2).ok_or_else(arithmetic)?;
        Ok(format!(
            "rvK3CMPeriods[{name}] = tau1=({} + sqrt(-{d}))/{}, tau2=({} + sqrt(-{d}))/2, branch=positive_imaginary",
            signed_text(self.form.b.negate()), decimal(two_a), signed_text(self.form.b),
        ))
    }

    pub(crate) fn class_number(&self, model: &mut Model) -> Result<usize> {
        let (_, discriminant) = self.order_data(model)?;
        let d = discriminant.small().ok_or_else(class_number_limit)?;
        primitive_class_number(model, d)
    }
}

pub(super) fn primitive_class_number(model: &mut Model, d: u128) -> Result<usize> {
    let mut count = 0;
    each_primitive_reduced(model, d, |_| {
        count += 1;
        Ok(())
    })?;
    Ok(count)
}

pub(super) fn each_primitive_reduced(
    model: &mut Model,
    d: u128,
    mut visit: impl FnMut([i128; 3]) -> Result<()>,
) -> Result<()> {
    // Reduction implies |d| >= 3a^2. There are a candidate b of the
    // required parity in the half-open interval -a < b <= a.
    let max_a = (d / 3).isqrt();
    let candidates = max_a * (max_a + 1) / 2;
    let cost = usize::try_from(candidates)
        .ok()
        .and_then(|n| n.checked_mul(24))
        .ok_or_else(class_number_limit)?;
    model.spend(cost).map_err(|_| class_number_limit())?;
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
                    visit([a as i128, b, c as i128])?;
                }
            }
            b += 2;
        }
    }
    Ok(())
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
    let mut coordinates = Vec::with_capacity(RANK);
    for (u, v) in u.into_iter().zip(v) {
        coordinates.push((exact_divide(model, u.resize()?, content)?, v.resize()?));
    }
    bezout_projection(model, coordinates, modulus)?
        .ok_or_else(|| Error::plain("K3 CM saturation requires a primitive first vector"))
}

pub(super) fn bezout_projection(
    model: &mut Model,
    coordinates: impl IntoIterator<Item = (Integer, Integer)>,
    modulus: Magnitude,
) -> Result<Option<Magnitude>> {
    // Track the gcd of the first coordinates and apply the same Bézout
    // combination to the second coordinates modulo the supplied modulus.
    let mut gcd = Magnitude::default();
    let mut residue = Magnitude::default();
    for (coordinate, value) in coordinates {
        if coordinate.is_zero() {
            continue;
        }
        let (mut a, mut b) = (gcd, coordinate.magnitude);
        let (mut left, mut right) = (
            residue,
            residue_mod(
                model,
                if coordinate.negative {
                    value.negate()
                } else {
                    value
                },
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
            return Ok(Some(residue));
        }
    }
    Ok(None)
}

pub(super) fn residue_mod(
    model: &mut Model,
    value: Integer,
    modulus: Magnitude,
) -> Result<Magnitude> {
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

pub(super) fn exact_divide(
    model: &mut Model,
    value: Integer,
    divisor: Magnitude,
) -> Result<Integer> {
    value
        .divide_exact(model, positive(divisor))?
        .ok_or_else(arithmetic)
}

fn gram_form(model: &mut Model, [u, v]: &[Coordinates; 2]) -> Result<Form> {
    model.spend(WORK_LIMBS * 4 * RANK * 3)?;
    let a = require_positive(pair_coordinates(|i, j| u[i].multiply_wide(u[j]))?)?;
    let b = pair_coordinates(|i, j| u[i].multiply_wide(v[j]))?;
    let c = require_positive(pair_coordinates(|i, j| v[i].multiply_wide(v[j]))?)?;
    let (a, remainder_a) = a.divide_u32(2);
    let (c, remainder_c) = c.divide_u32(2);
    if remainder_a != 0 || remainder_c != 0 {
        return Err(Error::plain("K3 CM basis must have even square"));
    }
    Ok(Form { a, b, c })
}

pub(super) fn positive(magnitude: Magnitude) -> Integer {
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
