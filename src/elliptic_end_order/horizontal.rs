use std::rc::Rc;

use super::{Certificate, Evidence};
use crate::cm_artin::PrimeIdeal;
use crate::ec_poly::{self as poly, Field, Polynomial, Scalar};
use crate::elliptic::{self, Curve, DivisionPolynomials, Kernel};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::cm_ring::Order;
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::syntax::Label;

pub(crate) struct Address(pub(crate) Label, pub(crate) usize, pub(crate) bool);

pub(crate) enum Statement {
    Data(Address),
    Define(Label, Label, Address),
}

pub(crate) struct Action {
    pub(crate) certificate: Rc<Certificate>,
    pub(super) curve: Rc<Curve>,
    pub(super) ideal: PrimeIdeal,
    pub(crate) ell: u64,
    lambda: u64,
    other: u64,
    q: u64,
}

impl Action {
    pub(crate) fn new(model: &mut Model, address: &Address) -> Result<Self> {
        let certificate = model.ec_end_order(&address.0.text)?;
        let curve = model.ec_curve(certificate.curve_name())?;
        let ell = elliptic::validate_odd_ell(model, &curve, address.1)?;
        let modulus = i128::from(ell);
        let frobenius = certificate.frobenius;
        if frobenius.generator.rem_euclid(modulus) == 0 {
            return Err(Error::plain(
                "CM horizontal prime divides the Frobenius order index and has repeated eigenvalues",
            ));
        }
        let order = Order::new(certificate.discriminant())?;
        let ideal = PrimeIdeal::new(model, order, Prime::validated(ell as u32), address.2)?;
        let (q, trace) = curve.trace(model)?;
        let omega = ideal.generator_residue(order);
        let lambda = (frobenius.constant + frobenius.generator * omega).rem_euclid(modulus);
        let trace = i128::from(trace).rem_euclid(modulus);
        let other = (trace - lambda).rem_euclid(modulus);
        if lambda == 0
            || lambda == other
            || (lambda * lambda - trace * lambda + i128::from(q)).rem_euclid(modulus) != 0
        {
            return Err(Error::plain(
                "CM horizontal prime requires distinct nonzero Frobenius eigenvalues",
            ));
        }
        Ok(Self {
            certificate,
            curve,
            ideal,
            ell,
            lambda: lambda as u64,
            other: other as u64,
            q,
        })
    }

    pub(crate) fn report(&self, model: &mut Model, address: &Address) -> Result<String> {
        let (kernel, x_degree) = self.kernel_with_degree(model, &self.curve)?;
        Ok(format!(
            "rvECCMHorizontalData[{},{},{}] = curve={}, eigenvalue={}, other_eigenvalue={}, x_candidate_degree={x_degree}, kernel={:?}",
            address.0.text,
            self.ell,
            usize::from(address.2),
            self.certificate.curve_name(),
            self.lambda,
            self.other,
            kernel.polynomial,
        ))
    }

    pub(crate) fn kernel(&self, model: &mut Model) -> Result<Kernel> {
        self.kernel_on(model, &self.curve)
    }

    pub(super) fn kernel_on(&self, model: &mut Model, curve: &Curve) -> Result<Kernel> {
        self.kernel_with_degree(model, curve)
            .map(|(kernel, _)| kernel)
    }

    fn kernel_with_degree(&self, model: &mut Model, curve: &Curve) -> Result<(Kernel, usize)> {
        let field = curve.field();
        let ell = self.ell;
        let size = (ell - 1) as usize / 2;
        let psi = elliptic::division_polynomial(model, curve, ell)?;
        // Evaluate [lambda] with its smaller absolute representative. The
        // sign is invisible in x and must be restored in the y multiplier.
        let negative = self.lambda > ell / 2;
        let n = if negative {
            ell - self.lambda
        } else {
            self.lambda
        };
        let (x_num, x_den) =
            DivisionPolynomials::new(model, curve, Some(&psi))?.x_image(model, n)?;
        let x = field.polynomial(&[0, 1]);
        let x_frobenius = poly::power_mod(model, field, &x, self.q, &psi)?;
        let x_equation = equation(model, field, &psi, &x_frobenius, &x_num, &x_den)?;
        let candidate = poly::monic_gcd(model, field, psi, x_equation)?;
        let x_degree = candidate.len() - 1;
        if x_degree < size
            || !poly::coprime(model, field, candidate.clone(), x_den)?
            || !poly::coprime(model, field, candidate.clone(), curve.polynomial().to_vec())?
        {
            return Err(Error::plain(
                "CM horizontal x equation has inconsistent degree or poles",
            ));
        }
        let (mut y_num, y_den) =
            DivisionPolynomials::new(model, curve, Some(&candidate))?.y_multiplier(model, n)?;
        if negative {
            y_num = poly::scale(model, field, &y_num, &field.signed(-1))?;
        }
        if !poly::coprime(model, field, candidate.clone(), y_den.clone())? {
            return Err(Error::plain(
                "CM horizontal y multiplier has a pole on the torsion candidate",
            ));
        }
        // For a nonzero odd-torsion point y != 0 and
        // y^q / y = (x^3+A*x+B)^((q-1)/2).
        let y_frobenius = poly::power_mod(
            model,
            field,
            curve.polynomial(),
            (self.q - 1) / 2,
            &candidate,
        )?;
        let y_equation = equation(model, field, &candidate, &y_frobenius, &y_num, &y_den)?;
        let kernel = poly::monic_gcd(model, field, candidate, y_equation)?;
        if kernel.len() != size + 1 {
            return Err(Error::plain(
                "signed CM horizontal equations failed the expected kernel degree",
            ));
        }
        Ok((Kernel::new(model, curve, ell as usize, kernel)?, x_degree))
    }

    pub(crate) fn check_quotient(
        &self,
        model: &mut Model,
        source: &Curve,
        target: &Curve,
    ) -> Result<()> {
        if source.trace(model)? != target.trace(model)? {
            return Err(Error::plain(
                "CM horizontal quotient changed the Frobenius polynomial",
            ));
        }
        let polynomial = model
            .cm_class_polynomial_at(self.certificate.discriminant())
            .ok_or_else(|| {
                Error::plain("CM horizontal quotient requires the certified ring class polynomial")
            })?;
        if !Evidence::read(model, target, &polynomial)?.certified() {
            return Err(Error::plain(
                "CM horizontal quotient failed the original endomorphism-order certificate",
            ));
        }
        if matches!(self.ell, 3 | 5 | 7) {
            let modular = model.elliptic_modular_polynomial(Prime::validated(self.ell as u32))?;
            let x = source.j_invariant(model)?;
            let y = target.j_invariant(model)?;
            if !poly::zero(&modular.evaluate(model, &source.ring, &x, &y)?) {
                return Err(Error::plain(
                    "CM horizontal quotient failed the modular polynomial",
                ));
            }
        }
        Ok(())
    }
}

fn equation(
    model: &mut Model,
    field: Field<'_>,
    modulus: &[Scalar],
    frobenius: &[Scalar],
    numerator: &[Scalar],
    denominator: &[Scalar],
) -> Result<Polynomial> {
    let left = poly::multiply_mod(model, field, frobenius, denominator, modulus)?;
    poly::sum(model, field, &left, numerator, true)
}
