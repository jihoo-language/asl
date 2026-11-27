use std::rc::Rc;

use super::{inconsistent, Hom};
use crate::ec_morphism::Morphism;
use crate::elliptic::{Curve, Point};
use crate::elliptic_cm_recognition as elliptic_maps;
use crate::elliptic_gluing::{Factor, Gluing, Side};
use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::hyperelliptic::Residues;
use crate::hyperelliptic_jacobian::Divisor;
use crate::hyperelliptic_kummer_surface::k3;
use crate::jacobian_level::brauer::{exterior_square, identity, Exterior};
use crate::jacobian_level::cm_cycles::{orthogonal_projector, symplectic, wedge, Cycle};
use crate::jacobian_level::linking::{
    add_vectors, exterior_change, image, intersection, product, scaled, verify_isometry, Fixed,
};
use crate::jacobian_level::Basis;
use crate::model::Model;
use crate::rv_unram::{add_mod, mul_mod, sub_mod};
use crate::syntax::Label;

pub(crate) struct Marking {
    pub(crate) hom_name: Label,
    base: Label,
    jacobian: Rc<Basis>,
    ring: Residues,
    vectors: [Fixed; 4],
    action: Exterior,
    projector: Option<Exterior>,
    map_columns: [[[u64; 2]; 2]; 2],
    polarization: [u64; 2],
}

impl Marking {
    pub(crate) fn new(
        model: &mut Model,
        hom_name: Label,
        hom: &Hom,
        base: Label,
        jacobian: Rc<Basis>,
        elliptic: [&Basis; 2],
    ) -> Result<Self> {
        if jacobian.len() != 4
            || jacobian.prime.get() == 2
            || elliptic
                .iter()
                .any(|b| b.len() != 2 || b.order != jacobian.order)
        {
            return Err(Error::plain("CM Hom marking needs two elliptic and one genus-two basis at the same odd prime-power level"));
        }
        let chart = Label {
            text: jacobian.curve_name.to_string(),
            span: base.span,
        };
        let (extended, maps) = Gluing::picard_maps(model, &chart)?;
        let original = model.elliptic_gluing(&base)?;
        let base_factors = sides(&original);
        let ext_factors = sides(&extended);
        if base_factors
            .iter()
            .zip(&hom.endpoints)
            .any(|(factor, endpoint)| {
                factor.curve.grade.ring != endpoint.ring_name
                    || factor.curve.polynomial() != endpoint.polynomial()
            })
            || ext_factors
                .iter()
                .zip(elliptic)
                .any(|(factor, basis)| factor.name != basis.curve_name)
        {
            return Err(Error::plain("CM Hom marking must use the actual product factors of the maps and gluing, in order"));
        }
        require_gluing_extension(model, &original, &extended)?;
        let [left, right] = ext_factors;
        let source = short_model(model, left)?;
        let target = short_model(model, right)?;
        let extended_models = [&source, &target];
        let names = [left.name.clone(), right.name.clone()];
        let old_models = hom.endpoints.each_ref().map(Rc::as_ref);
        let [phi, psi] = [
            elliptic_maps::extend_over(
                model,
                &hom.maps[0],
                names.clone(),
                old_models,
                extended_models,
            )?,
            elliptic_maps::extend_over(model, &hom.maps[1], names, old_models, extended_models)?,
        ];
        let ring = Residues::new(jacobian.prime, jacobian.height)?;
        let n = ring.modulus;
        let coordinates_of_image =
            |model: &mut Model, side: usize, divisor: &Divisor| -> Result<[u64; 4]> {
                let image = maps[side].pull(model, divisor)?;
                jacobian
                    .checked_coordinates(model, &image)?
                    .try_into()
                    .map_err(|_| inconsistent())
            };
        let [p, q] = elliptic[0].generators() else {
            unreachable!()
        };
        let [r, t] = elliptic[1].generators() else {
            unreachable!()
        };
        let p1 = coordinates_of_image(model, 0, p)?;
        let q1 = coordinates_of_image(model, 0, q)?;
        let p2 = coordinates_of_image(model, 1, r)?;
        let q2 = coordinates_of_image(model, 1, t)?;
        let polarization = [symplectic(p1, q1, n), symplectic(p2, q2, n)];
        for (i, &unit) in polarization.iter().enumerate() {
            if jacobian.polarization_scale(model, elliptic[i], 2)? != unit {
                return Err(inconsistent());
            }
        }
        if [p1, q1]
            .into_iter()
            .any(|x| [p2, q2].into_iter().any(|y| symplectic(x, y, n) != 0))
        {
            return Err(inconsistent());
        }
        let left_unit = ring.inverse_unit(model, &polarization[0])?;
        let right_unit = ring.inverse_unit(model, &polarization[1])?;
        let mut map_columns = [[[0; 2]; 2]; 2];
        let mut vectors = [
            scaled(wedge(p1, q1, n), left_unit, n),
            scaled(wedge(p2, q2, n), right_unit, n),
            [0; 6],
            [0; 6],
        ];
        for (i, morphism) in [phi, psi].iter().enumerate() {
            let p_image = apply(model, morphism, extended_models, p, right)?;
            let q_image = apply(model, morphism, extended_models, q, right)?;
            let matrix = &mut map_columns[i];
            for (entry, divisor) in matrix.iter_mut().zip([&p_image, &q_image]) {
                *entry = elliptic[1]
                    .checked_coordinates(model, divisor)?
                    .try_into()
                    .map_err(|_| inconsistent())?;
            }
            let [a, b] = *matrix;
            let det = sub_mod(mul_mod(a[0], b[1], n), mul_mod(a[1], b[0], n), n);
            if mul_mod(det, polarization[1], n)
                != mul_mod(morphism.degree() as u64, polarization[0], n)
            {
                return Err(inconsistent());
            }
            let graph_point = |v: [u64; 4], coefficients: [u64; 2]| -> [u64; 4] {
                std::array::from_fn(|j| {
                    add_mod(
                        v[j],
                        add_mod(
                            mul_mod(p2[j], coefficients[0], n),
                            mul_mod(q2[j], coefficients[1], n),
                            n,
                        ),
                        n,
                    )
                })
            };
            vectors[2 + i] = scaled(
                wedge(graph_point(p1, a), graph_point(q1, b), n),
                left_unit,
                n,
            );
        }

        // Transport the actual base q-Frobenius through the verified degree-four
        // product isogeny. This works when the odd-degree branch chart itself
        // lives only over an extension of the base field.
        let left_frob = elliptic[0].frobenius(model, &base_factors[0].curve, 1)?;
        let right_frob = elliptic[1].frobenius(model, &base_factors[1].curve, 1)?;
        let q_mod = hom.q % n;
        if left_frob.multiplier != q_mod || right_frob.multiplier != q_mod {
            return Err(inconsistent());
        }
        let mut entries = vec![0; 16];
        for (offset, f) in [&left_frob.matrix, &right_frob.matrix]
            .into_iter()
            .enumerate()
        {
            for i in 0..2 {
                for j in 0..2 {
                    entries[(2 * i + offset) * 4 + (2 * j + offset)] = *f.entry(i, j);
                }
            }
        }
        let product_frob = Matrix::from_entries(4, entries)?;
        let slots = [p1, p2, q1, q2];
        let s = Matrix::from_entries(
            4,
            (0..4)
                .flat_map(|i| slots.iter().map(move |col| col[i]))
                .collect(),
        )?;
        let inverse = hw_algebra::inverse(model, &ring, &s)?.ok_or_else(inconsistent)?;
        let first = hw_algebra::multiply(model, &ring, &s, &product_frob)?;
        let frobenius = hw_algebra::multiply(model, &ring, &first, &inverse)?;
        let twist = ring.inverse_unit(model, &q_mod)?;
        let action = exterior_square(&frobenius, n).map(|row| scaled(row, twist, n));
        verify_isometry(model, n, &action)?;
        for (i, &v) in vectors.iter().enumerate() {
            model.spend(36)?;
            if image(&action, v, n) != v {
                return Err(inconsistent());
            }
            for (j, &w) in vectors.iter().enumerate() {
                let pair = intersection(v, w, n);
                if add_mod(pair, pair, n) != hom.gram[i][j].rem_euclid(i128::from(n)) as u64 {
                    return Err(inconsistent());
                }
            }
        }
        let (q_wide, trace) = (i128::from(hom.q), hom.trace);
        let rank = k3::rational_abelian_rank(
            model,
            hom.q,
            &[
                1,
                -2 * trace,
                trace * trace + 2 * q_wide,
                -2 * q_wide * trace,
                q_wide.checked_mul(q_wide).ok_or_else(inconsistent)?,
            ],
            1,
        )?;
        if rank != 4 {
            return Err(inconsistent());
        }
        let projector = if hom
            .span_discriminant
            .rem_euclid(i128::from(ring.prime.get()))
            != 0
        {
            let operator = orthogonal_projector(model, ring, &vectors)?;
            model.spend(432)?;
            if product(&action, &operator, n) != operator
                || product(&operator, &action, n) != operator
            {
                return Err(inconsistent());
            }
            Some(operator)
        } else {
            None
        };
        model.allocate(128 + hom_name.text.len() + base.text.len())?;
        Ok(Self {
            hom_name,
            base,
            jacobian,
            ring,
            vectors,
            action,
            projector,
            map_columns,
            polarization,
        })
    }

    pub(crate) fn modulus(&self) -> u64 {
        self.ring.modulus
    }

    pub(crate) fn vector(&self, vector: [i128; 6]) -> Fixed {
        vector.map(|x| x.rem_euclid(i128::from(self.modulus())) as u64)
    }

    pub(crate) fn cycle(&self, model: &mut Model, vector: Cycle) -> Result<Fixed> {
        model.spend(24)?;
        let modulus = self.modulus();
        Ok(self
            .vectors
            .iter()
            .zip(vector)
            .fold([0; 6], |sum, (&v, c)| {
                add_vectors(
                    sum,
                    scaled(v, c.rem_euclid(i128::from(modulus)) as u64, modulus),
                    modulus,
                )
            }))
    }

    pub(crate) fn projector(&self, transcendental: bool) -> Result<Exterior> {
        let projector = self.projector.ok_or_else(|| {
            Error::plain(
                "CM Hom marked Gram is not a local unit; no orthogonal projector is certified",
            )
        })?;
        let n = self.modulus();
        Ok(if transcendental {
            std::array::from_fn(|i| {
                std::array::from_fn(|j| sub_mod(u64::from(i == j), projector[i][j], n))
            })
        } else {
            projector
        })
    }

    pub(crate) fn power(&self, model: &mut Model, mut power: usize) -> Result<Exterior> {
        let n = self.modulus();
        let mut operator = self.action;
        let mut result = identity();
        while power != 0 {
            if power & 1 == 1 {
                model.spend(216)?;
                result = product(&result, &operator, n);
            }
            power >>= 1;
            if power != 0 {
                model.spend(216)?;
                operator = product(&operator, &operator, n);
            }
        }
        Ok(result)
    }

    pub(crate) fn change(&self, model: &mut Model, source: &Self) -> Result<Exterior> {
        if self.base.text != source.base.text || self.ring != source.ring {
            return Err(Error::plain(
                "CM Hom marking change requires the same base gluing and level",
            ));
        }
        let matrix = exterior_change(model, &self.jacobian, &source.jacobian, self.ring)?;
        model.spend(432)?;
        let n = self.modulus();
        if product(&self.action, &matrix, n) != product(&matrix, &source.action, n) {
            return Err(inconsistent());
        }
        Ok(matrix)
    }

    pub(crate) fn report(&self, name: &str) -> String {
        let n = self.modulus();
        let mut output = format!("rvCMHomMarkData[{name}] = level={n}, polarization_units={:?}, map_columns={:?}, cycle_lifts={:?}, twisted_frobenius={:?}, abelian_tate_rank=4, locally_saturated={}",
            self.polarization,self.map_columns,self.vectors,self.action,self.projector.is_some());
        if let Some(projector) = self.projector {
            output.push_str(&format!(", algebraic_projector={projector:?}"));
        }
        output
    }
}

fn sides(gluing: &Gluing) -> [&Factor; 2] {
    [gluing.factor(Side::Left), gluing.factor(Side::Right)]
}

fn short_model(model: &mut Model, factor: &Factor) -> Result<Curve> {
    let f = factor.curve.polynomial();
    if f.len() != 4 || f[2] != factor.curve.ring.scalar(0) || f[3] != factor.curve.ring.scalar(1) {
        return Err(inconsistent());
    }
    Curve::new(
        model,
        factor.curve.grade.ring.clone(),
        factor.curve.ring.clone(),
        f[1].clone(),
        f[0].clone(),
    )
}

fn require_gluing_extension(model: &mut Model, base: &Gluing, extension: &Gluing) -> Result<()> {
    let coefficients = |g: &Gluing| {
        [
            g.pole.clone(),
            g.image_of_infinity.clone(),
            g.numerator.clone(),
            g.scale.clone(),
        ]
    };
    let embedded = base.curve.ring.embed_residue_coefficients(
        model,
        &extension.curve.ring,
        &coefficients(base),
    )?;
    if embedded != coefficients(extension) {
        return Err(Error::plain(
            "CM Hom marking needs the coefficient extension of the named base gluing",
        ));
    }
    Ok(())
}

fn apply(
    model: &mut Model,
    map: &Morphism,
    [source, target]: [&Curve; 2],
    divisor: &Divisor,
    range: &Factor,
) -> Result<Divisor> {
    let point = if divisor.u.len() == 1 {
        Point::Infinity
    } else {
        let x = source.ring.negate(model, &divisor.u[0])?;
        let y = divisor
            .v
            .first()
            .cloned()
            .unwrap_or_else(|| source.ring.scalar(0));
        source.point(model, x, y)?
    };
    match map.apply(model, source, target, &point)? {
        Point::Infinity => Ok(Divisor::zero(&range.curve)),
        Point::Affine(x, y) => Divisor::from_point(model, &range.curve, &x, y),
    }
}

pub(crate) fn apply_operator(
    model: &mut Model,
    operator: &Exterior,
    modulus: u64,
    vector: Fixed,
) -> Result<Fixed> {
    model.spend(36)?;
    Ok(image(operator, vector, modulus))
}

pub(crate) fn pair(x: Fixed, y: Fixed, n: u64) -> u64 {
    intersection(x, y, n)
}
