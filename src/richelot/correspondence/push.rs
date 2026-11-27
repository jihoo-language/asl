use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hyperelliptic_jacobian::{quadratic_trace, Divisor};
use crate::model::Model;

use super::Correspondence;

impl Correspondence {
    pub(crate) fn push(&self, model: &mut Model, input: &Divisor) -> Result<Divisor> {
        let field = Field {
            ring: &self.source.ring,
        };
        match input.u.as_slice() {
            [_] => Ok(Divisor::zero(&self.target)),
            [constant, _] => {
                let x = field.ring.negate(model, constant)?;
                let y = poly::evaluate(model, field, &input.v, &x)?;
                self.point(model, &x, &y)
            }
            [constant, linear, _] => self.push_quadratic(model, field, input, constant, linear),
            _ => unreachable!("a reduced genus-two divisor has degree at most two"),
        }
    }

    fn push_quadratic(
        &self,
        model: &mut Model,
        field: Field<'_>,
        input: &Divisor,
        constant: &Scalar,
        linear: &Scalar,
    ) -> Result<Divisor> {
        // The source support is (X-h)²-k, so over the quadratic algebra the
        // two points have coordinates X = h±ω with ω² = k.
        let minus_half = field.inverse_unit(model, &field.signed(-2))?;
        let h = field.multiply(model, linear, &minus_half)?;
        let h_squared = field.multiply(model, &h, &h)?;
        let k = field.subtract(model, &h_squared, constant)?;
        if field.is_zero(&k) {
            let y = poly::evaluate(model, field, &input.v, &h)?;
            let image = self.point(model, &h, &y)?;
            return image.add(model, &self.target, &image);
        }
        if input.v.is_empty() {
            return self.push_branch_support(model, field, &input.u);
        }

        let algebra = Quadratic { field, square: &k };
        let y = algebra.shifted_value(model, &input.v, &h)?;
        let norm = algebra.norm(model, y.each_ref())?;
        if field.is_zero(&norm) {
            return self.push_split(model, field, input, &h, y.each_ref());
        }

        let [i, j] = self.quadratic_indices();
        let g1 = algebra.shifted_value(model, &self.isogeny.source_factors[i], &h)?;
        let g2 = algebra.shifted_value(model, &self.isogeny.source_factors[j], &h)?;
        let inverse = algebra.inverse(model, y.each_ref(), &norm)?;
        let ratio = algebra.multiply(model, g1.each_ref(), inverse.each_ref())?;
        let x_ratio = algebra.multiply(model, [&h, &field.one()], ratio.each_ref())?;
        let (u_real, v_real) =
            self.fiber(model, field, [&g1[0], &g2[0]], [&x_ratio[0], &ratio[0]])?;
        let (u_imag, v_imag) =
            self.fiber(model, field, [&g1[1], &g2[1]], [&x_ratio[1], &ratio[1]])?;

        let width = u_real.len().max(u_imag.len());
        let index = width
            .checked_sub(1)
            .ok_or_else(|| Error::plain("Richelot fiber support vanished"))?;
        let zero = field.zero();
        let leading = [
            u_real.get(index).unwrap_or(&zero),
            u_imag.get(index).unwrap_or(&zero),
        ];
        let norm = algebra.norm(model, leading)?;
        if field.is_zero(&norm) {
            // One component loses a point at target infinity. The leading
            // zero divisor splits the algebra so each component can be reduced.
            return self.push_split(model, field, input, &h, leading);
        }
        let inverse = algebra.inverse(model, leading, &norm)?;
        let u = algebra.scale_support(model, [&u_real, &u_imag], inverse.each_ref())?;
        quadratic_trace::trace(model, &self.target, &k, u, [v_real, v_imag])
    }

    fn push_branch_support(
        &self,
        model: &mut Model,
        field: Field<'_>,
        support: &Polynomial,
    ) -> Result<Divisor> {
        let [i, j] = self.quadratic_indices();
        let mut image = Divisor::zero(&self.target);
        for (source, target) in [(i, j), (j, i)] {
            let overlap = poly::gcd(
                model,
                field,
                support.clone(),
                self.isogeny.source_factors[source].clone(),
            )?;
            if overlap.len().is_multiple_of(2) {
                let pair = self.branch_pair(model, field, target)?;
                image = image.add(model, &self.target, &pair)?;
            }
        }
        Ok(image)
    }

    fn push_split(
        &self,
        model: &mut Model,
        field: Field<'_>,
        input: &Divisor,
        h: &Scalar,
        [real, imaginary]: [&Scalar; 2],
    ) -> Result<Divisor> {
        // A nonzero zero divisor a+bω gives k=(a/b)². Its two evaluations
        // at ω=±a/b give the rational points in the source support.
        let inverse = field.inverse_unit(model, imaginary)?;
        let offset = field.multiply(model, real, &inverse)?;
        let x_left = field.add(model, h, &offset)?;
        let x_right = field.subtract(model, h, &offset)?;
        let y_left = poly::evaluate(model, field, &input.v, &x_left)?;
        let y_right = poly::evaluate(model, field, &input.v, &x_right)?;
        let left = self.point(model, &x_left, &y_left)?;
        let right = self.point(model, &x_right, &y_right)?;
        left.add(model, &self.target, &right)
    }
}

struct Quadratic<'a> {
    field: Field<'a>,
    square: &'a Scalar,
}

impl Quadratic<'_> {
    fn shifted_value(
        &self,
        model: &mut Model,
        polynomial: &[Scalar],
        h: &Scalar,
    ) -> Result<[Scalar; 2]> {
        let field = self.field;
        let mut real = poly::evaluate(model, field, polynomial, h)?;
        let mut imaginary = polynomial.get(1).cloned().unwrap_or_else(|| field.zero());
        if let Some(quadratic) = polynomial.get(2) {
            let shift = field.multiply(model, self.square, quadratic)?;
            real = field.add(model, &real, &shift)?;
            let term = field.multiply(model, h, quadratic)?;
            let twice = field.ring.scale(model, &term, 2)?;
            imaginary = field.add(model, &imaginary, &twice)?;
        }
        Ok([real, imaginary])
    }

    fn norm(&self, model: &mut Model, [real, imaginary]: [&Scalar; 2]) -> Result<Scalar> {
        let field = self.field;
        let a_squared = field.multiply(model, real, real)?;
        let b_squared = field.multiply(model, imaginary, imaginary)?;
        let term = field.multiply(model, self.square, &b_squared)?;
        field.subtract(model, &a_squared, &term)
    }

    fn inverse(
        &self,
        model: &mut Model,
        [real, imaginary]: [&Scalar; 2],
        norm: &Scalar,
    ) -> Result<[Scalar; 2]> {
        let field = self.field;
        let inverse = field.inverse_unit(model, norm)?;
        let real = field.multiply(model, &inverse, real)?;
        let imaginary = field.multiply(model, &inverse, imaginary)?;
        Ok([real, field.ring.negate(model, &imaginary)?])
    }

    fn multiply(
        &self,
        model: &mut Model,
        [a, b]: [&Scalar; 2],
        [c, d]: [&Scalar; 2],
    ) -> Result<[Scalar; 2]> {
        let field = self.field;
        let ac = field.multiply(model, a, c)?;
        let bd = field.multiply(model, b, d)?;
        let kbd = field.multiply(model, self.square, &bd)?;
        let ad = field.multiply(model, a, d)?;
        let bc = field.multiply(model, b, c)?;
        Ok([field.add(model, &ac, &kbd)?, field.add(model, &ad, &bc)?])
    }

    fn scale_support(
        &self,
        model: &mut Model,
        [real, imaginary]: [&[Scalar]; 2],
        scale: [&Scalar; 2],
    ) -> Result<[Polynomial; 2]> {
        let zero = self.field.zero();
        let mut output = [Vec::new(), Vec::new()];
        for i in 0..real.len().max(imaginary.len()) {
            let value = [
                real.get(i).unwrap_or(&zero),
                imaginary.get(i).unwrap_or(&zero),
            ];
            let [a, b] = self.multiply(model, value, scale)?;
            output[0].push(a);
            output[1].push(b);
        }
        Ok(output.map(poly::trim))
    }
}
