use std::rc::Rc;

use super::Certificate;
use crate::error::{Error, Result};
use crate::legendre_cm::Element;
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::rv_unram::{self, Ring};
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Expression {
    JInvariant(Label, Label),
    ModularResidual(Label, Label),
    Root(Label, Label, bool),
    Embed(Label, Label, bool, Box<ParsedExpr>),
}

pub(crate) struct Lift {
    certificate: Rc<Certificate>,
    ring: Rc<Ring>,
}

impl Lift {
    pub(crate) fn new(model: &Model, certificate: Rc<Certificate>, ring: Rc<Ring>) -> Result<Self> {
        let source = model.ec_curve(&certificate.curve)?;
        if !ring.reduces_to(&source.ring) {
            return Err(Error::plain(
                "canonical lift requires an unramified ring with the curve's residue presentation",
            ));
        }
        Ok(Self { certificate, ring })
    }

    pub(crate) fn j_invariant(&self, model: &mut Model) -> Result<Vec<u64>> {
        let certificate = &self.certificate;
        let source = model.ec_curve(&certificate.curve)?;
        let class = model
            .cm_class_polynomial_at(certificate.discriminant)
            .ok_or_else(|| {
                Error::plain("canonical lift requires the certified ring class polynomial")
            })?;
        let coefficients = class.unramified_coefficients(model, &self.ring)?;
        let j = simple_root(model, &self.ring, &certificate.j, &coefficients)?;
        let sigma = self.ring.frobenius(model, &j, 1)?;
        let expected = source.ring.frobenius(model, &certificate.j, 1)?;
        if self.ring.residue(&sigma) != expected
            || !zero(
                &self
                    .ring
                    .evaluate_polynomial(model, &coefficients, &sigma)?,
            )
        {
            return Err(Error::plain(
                "canonical j lift failed Frobenius equivariance",
            ));
        }
        Ok(j)
    }

    pub(crate) fn modular_residual(&self, model: &mut Model) -> Result<Vec<u64>> {
        let prime = Prime::validated(self.ring.prime() as u32);
        let polynomial = model.elliptic_modular_polynomial(prime)?;
        let j = self.j_invariant(model)?;
        let sigma = self.ring.frobenius(model, &j, 1)?;
        polynomial.evaluate(model, &self.ring, &j, &sigma)
    }

    pub(crate) fn root(&self, model: &mut Model, complementary: bool) -> Result<Vec<u64>> {
        let (alpha, beta) = self.roots(model)?;
        Ok(if complementary { beta } else { alpha })
    }

    fn roots(&self, model: &mut Model) -> Result<(Vec<u64>, Vec<u64>)> {
        let ring = &self.ring;
        let c = &self.certificate;
        let q = constant(ring, i128::from(c.q));
        let coefficients = vec![q.clone(), constant(ring, -c.trace), ring.scalar(1)];
        let seed = ring.scalar(c.trace.rem_euclid(i128::from(ring.prime())) as u64);
        let alpha = simple_root(model, ring, &seed, &coefficients)?;
        let inverse = ring.inverse(model, &alpha)?;
        let beta = ring.multiply(model, &q, &inverse)?;
        if ring.add(model, &alpha, &beta)? != constant(ring, c.trace)
            || ring.multiply(model, &alpha, &beta)? != q
        {
            return Err(Error::plain(
                "lifted ordinary Frobenius roots failed trace or norm",
            ));
        }
        Ok((alpha, beta))
    }

    pub(crate) fn embed(
        &self,
        model: &mut Model,
        value: Element,
        complementary: bool,
    ) -> Result<Vec<u64>> {
        let ring = &self.ring;
        let c = &self.certificate;
        let b = constant(ring, c.frobenius.generator);
        if zero(&ring.residue(&b)) {
            return Err(Error::plain(
                "Frobenius coefficient must be a unit for the ordinary CM embedding",
            ));
        }
        let root = self.root(model, complementary)?;
        let difference = ring.add(model, &root, &constant(ring, -c.frobenius.constant))?;
        let inverse = ring.inverse(model, &b)?;
        let omega = ring.multiply(model, &difference, &inverse)?;
        let relation = [
            constant(ring, c.generator_norm),
            constant(ring, -c.generator_trace),
            ring.scalar(1),
        ];
        if !zero(&ring.evaluate_polynomial(model, &relation, &omega)?) {
            return Err(Error::plain(
                "ordinary CM embedding failed the order relation",
            ));
        }
        let generator_part = ring.scale(
            model,
            &omega,
            rv_unram::reduce(value.generator, ring.modulus()),
        )?;
        ring.add(model, &constant(ring, value.constant), &generator_part)
    }
}

fn simple_root(
    model: &mut Model,
    ring: &Ring,
    seed: &[u64],
    coefficients: &[Vec<u64>],
) -> Result<Vec<u64>> {
    let mut root = seed.to_vec();
    let (mut value, mut derivative) = evaluate_pair(model, ring, coefficients, &root)?;
    if !zero(&ring.residue(&value)) {
        return Err(Error::plain("Hensel lift requires a residue root"));
    }
    if zero(&ring.residue(&derivative)) {
        return Err(Error::plain(
            "canonical Hensel lift requires a unit derivative",
        ));
    }
    let mut precision = 1;
    while precision < ring.precision() && !zero(&value) {
        let inverse = ring.inverse(model, &derivative)?;
        let correction = ring.multiply(model, &value, &inverse)?;
        let correction = ring.negate(model, &correction)?;
        root = ring.add(model, &root, &correction)?;
        (value, derivative) = evaluate_pair(model, ring, coefficients, &root)?;
        precision *= 2;
    }
    if !zero(&value) || ring.residue(&root) != seed {
        return Err(Error::plain(
            "canonical Hensel lift failed residual verification",
        ));
    }
    Ok(root)
}

fn evaluate_pair(
    model: &mut Model,
    ring: &Ring,
    coefficients: &[Vec<u64>],
    input: &[u64],
) -> Result<(Vec<u64>, Vec<u64>)> {
    let mut value = ring.scalar(0);
    let mut derivative = ring.scalar(0);
    for coefficient in coefficients.iter().rev() {
        let next = ring.multiply(model, &derivative, input)?;
        derivative = ring.add(model, &next, &value)?;
        let next = ring.multiply(model, &value, input)?;
        value = ring.add(model, &next, coefficient)?;
    }
    Ok((value, derivative))
}

fn constant(ring: &Ring, value: i128) -> Vec<u64> {
    ring.scalar(rv_unram::reduce(value, ring.modulus()))
}

fn zero(value: &[u64]) -> bool {
    value.iter().all(|&c| c == 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repeated_residue_root_is_rejected_even_when_the_initial_value_is_zero(
    ) -> std::result::Result<(), String> {
        let mut model = Model::new();
        let ring = Ring::new(&mut model, 7, 3, vec![0, 1]).map_err(|error| error.to_string())?;
        let coefficients = [ring.scalar(1), constant(&ring, -2), ring.scalar(1)];
        let result = simple_root(&mut model, &ring, &ring.scalar(1), &coefficients);
        assert!(matches!(result, Err(error) if error.message.contains("unit derivative")));
        Ok(())
    }
}
