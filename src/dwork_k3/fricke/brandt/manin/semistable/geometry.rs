use std::collections::BTreeMap;

use crate::ec_poly::{Field, Scalar};
use crate::error::Result;
use crate::finite_field_factor;
use crate::hw_algebra::CoefficientRing;
use crate::model::Model;

use super::{certificate, Basis};

pub(super) struct Edges {
    pub(super) parameters: Vec<Scalar>,
    pub(super) endpoints: Vec<[Scalar; 2]>,
    pub(super) widths: Vec<i128>,
    pub(super) frobenius: Vec<usize>,
    pub(super) w2: Vec<usize>,
    pub(super) quotient: Vec<usize>,
}

impl Edges {
    pub(super) fn new(model: &mut Model, basis: &Basis) -> Result<Self> {
        let field = Field { ring: &basis.ring };
        let mut oriented = BTreeMap::new();
        // (t+16)^3-j_2 t. Its roots represent isomorphism classes
        // of ordered cyclic 2-structures, independently of ramification.
        let base = field.polynomial(&[4096, 768, 48, 1]);
        for j in &basis.js {
            let mut cubic = base.clone();
            cubic[1] = field.subtract(model, &cubic[1], j)?;
            let factors = finite_field_factor::factor_with_multiplicity(model, field, &cubic)?;
            let mut degree = 0;
            for (d, multiplicity, f) in factors {
                if d != 1 {
                    return Err(certificate(
                        "ordered supersingular level-2 points split over Fp²",
                    ));
                }
                degree += multiplicity;
                let t = field.subtract(model, &field.zero(), &f[0])?;
                if field.is_zero(&t) || oriented.insert(t.clone(), j.clone()).is_some() {
                    return Err(certificate("distinct nonzero ordered level-2 parameters"));
                }
                let expected = if t == field.signed(-16) {
                    3
                } else if t == field.signed(8) {
                    2
                } else {
                    1
                };
                if multiplicity != expected {
                    return Err(certificate("coarse ramification at j=0 and j=1728"));
                }
            }
            if degree != 3 {
                return Err(certificate("complete degree-three j projection"));
            }
        }
        let mut endpoints = Vec::new();
        let mut widths = Vec::new();
        let mut frobenius = Vec::new();
        let mut w2 = Vec::new();
        let mut quotient = Vec::new();
        let (parameters, j2): (Vec<_>, Vec<_>) = oriented.into_iter().unzip();
        for (i, t) in parameters.iter().enumerate() {
            let inverse = field.inverse_unit(model, t)?;
            let opposite = basis.ring.scale(model, &inverse, 4096 % basis.p)?;
            let partner = parameters
                .binary_search(&opposite)
                .map_err(|_| certificate("W2 permutes ordered nodes"))?;
            let image = basis.ring.frobenius(model, t, 1)?;
            let image = parameters
                .binary_search(&image)
                .map_err(|_| certificate("Frobenius permutes ordered nodes"))?;
            let h = field.add(model, t, &opposite)?;
            let h = field.add(model, &h, &field.signed(24))?;
            let h = basis
                .vertices
                .binary_search(&h)
                .map_err(|_| certificate("ordered t projects to Fricke H"))?;
            let j = [j2[partner].clone(), j2[i].clone()];
            let mut unordered = j.clone();
            unordered.sort();
            if unordered != basis.endpoints[h] {
                return Err(certificate(
                    "ordered elliptic endpoints give the Fricke pair",
                ));
            }
            let width = if *t == field.signed(-64) {
                if j[1] != field.signed(1728) || partner != i {
                    return Err(certificate(
                        "the exceptional thickness-two supersingular node",
                    ));
                }
                2
            } else {
                1
            };
            widths.push(width);
            endpoints.push(j);
            frobenius.push(image);
            w2.push(partner);
            quotient.push(h);
        }
        for i in 0..parameters.len() {
            if frobenius[frobenius[i]] != i
                || w2[w2[i]] != i
                || frobenius[w2[i]] != w2[frobenius[i]]
                || widths[frobenius[i]] != widths[i]
                || widths[w2[i]] != widths[i]
                || quotient[w2[i]] != quotient[i]
            {
                return Err(certificate(
                    "the ordered involutions commute and preserve widths",
                ));
            }
        }
        if widths.iter().filter(|&&w| w == 2).count()
            != usize::from(basis.js.contains(&field.signed(1728)))
            || (0..basis.rank()).any(|i| !quotient.contains(&i))
        {
            return Err(certificate("complete exceptional and Fricke node support"));
        }
        let n = parameters.len();
        model.allocate(15 + 18 * n)?;
        Ok(Self {
            parameters,
            endpoints,
            widths,
            frobenius,
            w2,
            quotient,
        })
    }

    pub(super) fn oriented_matrix(
        &self,
        model: &mut Model,
        basis: &Basis,
        equation: &crate::elliptic_modular::Polynomial,
        ell: usize,
    ) -> Result<super::linear::Matrix> {
        let field = Field { ring: &basis.ring };
        let n = self.parameters.len();
        let factors = self
            .parameters
            .iter()
            .map(|t| Ok([field.subtract(model, &field.zero(), t)?, field.one()]))
            .collect::<Result<Vec<_>>>()?;
        let mut push = vec![vec![super::linear::ZERO; n]; n];
        for (i, t) in self.parameters.iter().enumerate() {
            let mut specialization = (0..=ell + 1)
                .map(|j| equation.specialization(model, &basis.ring, t, j))
                .collect::<Result<Vec<_>>>()?;
            if specialization.last() != Some(&field.one()) {
                return Err(certificate(
                    "monic specialization of oriented Hecke equation",
                ));
            }
            let mut sum = 0;
            for (j, factor) in factors.iter().enumerate() {
                let count =
                    finite_field_factor::remove_factor(model, field, &mut specialization, factor)?;
                sum += count;
                push[j][i] = crate::lee_linear::Rational::integer(count as i128);
            }
            if sum != ell + 1 || specialization != [field.one()] {
                return Err(certificate(
                    "ordered Hecke targets are complete with multiplicity",
                ));
            }
        }
        Ok(push)
    }
}
