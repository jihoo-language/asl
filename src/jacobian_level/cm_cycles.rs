use std::rc::Rc;

use super::brauer::{identity, Exterior, PAIRS};
use super::linking::{column, image, intersection, product, Fixed, Pairing};
use super::Basis;
use crate::ec_poly::{self as poly, Scalar};
use crate::elliptic_gluing::{Gluing, Side};
use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::hyperelliptic_jacobian::{require_extension, Divisor};
use crate::hyperelliptic_kummer_surface::k3;
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::rv_unram::{add_mod, mul_mod, sub_mod};
use crate::syntax::{Label, ParsedExpr};

pub(crate) type Cycle = [i128; 4];

pub(crate) enum Statement {
    Define(Label, Label, Label, usize, Box<ParsedExpr>),
    Data(Label),
    Pair(Box<[ParsedExpr; 2]>),
    Hermitian(Box<ParsedExpr>),
    Defect(Label, usize),
}

pub(crate) enum Expression {
    Cycle(Label, Cycle),
    Lift(usize, Box<ParsedExpr>),
    Project(Label, bool, Box<ParsedExpr>),
}

pub(crate) struct Cycles {
    pub(crate) pairing_name: Label,
    pub(crate) pairing: Rc<Pairing>,
    index: i128,
    gram: [[i128; 4]; 4],
    generators: [Fixed; 4],
    projector: Option<Exterior>,
    gaussian_columns: [[u64; 2]; 2],
    weil_scale: u64,
    trace: i128,
    q: u64,
}

impl Cycles {
    pub(crate) fn new(
        model: &mut Model,
        pairing_name: Label,
        pairing: Rc<Pairing>,
        elliptic: &Basis,
        index: usize,
        gaussian: Scalar,
    ) -> Result<Self> {
        if index == 0 {
            return Err(Error::plain(
                "Kummer CM graph generator index must be positive",
            ));
        }
        let index = index as i128;
        let norm = index.checked_mul(index).ok_or_else(index_overflow)?;
        norm.checked_mul(4).ok_or_else(index_overflow)?;
        let jacobian = &pairing.basis;
        if elliptic.len() != 2 || elliptic.order != jacobian.order {
            return Err(Error::plain(
                "Kummer CM cycles require an elliptic basis at the same odd prime-power level",
            ));
        }

        let base_chart = model.jacobian_chart(&pairing.base)?;
        let ext_chart_name = Label {
            text: jacobian.curve_name.to_string(),
            span: pairing_name.span,
        };
        let extended_chart = model.jacobian_chart(&ext_chart_name)?;
        let base = model.elliptic_gluing(&Label {
            text: base_chart.source.to_string(),
            span: pairing_name.span,
        })?;
        let (extended, maps) = Gluing::picard_maps(model, &ext_chart_name)?;
        let e = &base.factor(Side::Left).curve;
        let e_left = extended.factor(Side::Left);
        let e_right = extended.factor(Side::Right);
        if elliptic.curve_name != e_left.name {
            return Err(Error::plain(
                "Kummer CM elliptic basis must mark the left factor of the extended gluing",
            ));
        }
        if e.polynomial() != base.factor(Side::Right).curve.polynomial()
            || e_left.curve.polynomial() != e_right.curve.polynomial()
        {
            return Err(Error::plain(
                "Kummer CM cycles require a glued square of the same elliptic equation",
            ));
        }
        require_extension(model, e, &e_left.curve)?;
        let ring = &e.ring;
        let ring_ext = &e_left.curve.ring;
        let [branch, pole, infinity, numerator, scale, gaussian]: [Scalar; 6] = ring
            .embed_residue_coefficients(
                model,
                ring_ext,
                &[
                    base_chart.branch.clone(),
                    base.pole.clone(),
                    base.image_of_infinity.clone(),
                    base.numerator.clone(),
                    base.scale.clone(),
                    gaussian,
                ],
            )?
            .try_into()
            .map_err(|_| inconsistent())?;
        if (branch, pole, infinity, numerator, scale)
            != (
                extended_chart.branch.clone(),
                extended.pole.clone(),
                extended.image_of_infinity.clone(),
                extended.numerator.clone(),
                extended.scale.clone(),
            )
        {
            return Err(Error::plain("Kummer CM cycles require the base gluing and branch chart to extend coefficientwise"));
        }
        let f = e.polynomial();
        if !poly::zero(&f[0]) || !poly::zero(&f[2]) {
            return Err(Error::plain(
                "Gaussian Kummer CM cycles require y^2=x^3+A*x",
            ));
        }
        if ring_ext.multiply(model, &gaussian, &gaussian)? != ring_ext.scalar(ring_ext.prime() - 1)
        {
            return Err(Error::plain(
                "Gaussian Kummer CM cycles require a base-field square root of -1",
            ));
        }
        let trace = -e.weil(model)?[1];
        if trace.rem_euclid(i128::from(ring.prime())) == 0 {
            return Err(Error::plain(
                "Kummer CM cycles require an ordinary elliptic square",
            ));
        }
        let q = ring.size_of_residue_field()?;
        // The certified (2,2) isogeny gives the square of the elliptic
        // characteristic polynomial, independent of the observed cycle vectors.
        let q_wide = i128::from(q);
        let rank = k3::rational_abelian_rank(
            model,
            q,
            &[
                1,
                -2 * trace,
                trace * trace + 2 * q_wide,
                -2 * q_wide * trace,
                q_wide.checked_mul(q_wide).ok_or_else(inconsistent)?,
            ],
            pairing.degree as u128,
        )?;
        if rank != 4 {
            return Err(inconsistent());
        }

        let [p, t] = elliptic.generators() else {
            unreachable!("the elliptic basis has rank two")
        };
        let [p1, t1, p2, t2] =
            [(0, p), (0, t), (1, p), (1, t)].map(|(side, point)| -> Result<[u64; 4]> {
                let pullback = maps[side].pull(model, point)?;
                jacobian
                    .checked_coordinates(model, &pullback)?
                    .try_into()
                    .map_err(|_| inconsistent())
            });
        let [p1, t1, p2, t2] = [p1?, t1?, p2?, t2?];
        let residue = pairing.ring;
        let n = residue.modulus;
        let unit = symplectic(p1, t1, n);
        if jacobian.polarization_scale(model, elliptic, 2)? != unit
            || symplectic(p2, t2, n) != unit
            || [p1, t1]
                .into_iter()
                .any(|x| [p2, t2].into_iter().any(|y| symplectic(x, y, n) != 0))
        {
            return Err(inconsistent());
        }
        let inverse_unit = residue.inverse_unit(model, &unit)?;

        let gaussian_action = elliptic.action_matrix(model, elliptic, |model, point| {
            let x = point.u[0].clone(); // x mapped to -x
            let y =
                ring_ext.multiply(model, &gaussian, &point.coefficient(&e_left.curve, true, 0))?;
            Divisor::from_point(model, &e_left.curve, &x, y)
        })?;
        let gaussian_columns: [[u64; 2]; 2] =
            std::array::from_fn(|i| [*gaussian_action.entry(0, i), *gaussian_action.entry(1, i)]);
        let [a, b] = gaussian_columns;
        if add_mod(a[0], b[1], n) != 0
            || sub_mod(mul_mod(a[0], b[1], n), mul_mod(a[1], b[0], n), n) != 1
        {
            return Err(inconsistent());
        }
        let graph = |col: [u64; 2]| -> [u64; 4] {
            let index_mod = index.rem_euclid(i128::from(n)) as u64;
            std::array::from_fn(|j| {
                mul_mod(
                    index_mod,
                    add_mod(mul_mod(col[0], p2[j], n), mul_mod(col[1], t2[j], n), n),
                    n,
                )
            })
        };
        let plus = |v: [u64; 4], w: [u64; 4]| std::array::from_fn(|j| add_mod(v[j], w[j], n));
        let generators = [
            wedge(p1, t1, n),
            wedge(p2, t2, n),
            wedge(plus(p1, p2), plus(t1, t2), n),
            wedge(plus(p1, graph(a)), plus(t1, graph(b)), n),
        ]
        .map(|v| super::linking::scaled(v, inverse_unit, n));
        let gram = [
            [0, 1, 1, norm],
            [1, 0, 1, 1],
            [1, 1, 0, norm + 1],
            [norm, 1, norm + 1, 0],
        ];
        for (j, left) in generators.iter().enumerate() {
            pairing.require_fixed(model, n, *left)?;
            for (&right, &expected) in generators.iter().zip(&gram[j]) {
                let actual = intersection(*left, right, n);
                if add_mod(actual, actual, n) != expected.rem_euclid(i128::from(n)) as u64 {
                    return Err(inconsistent());
                }
            }
        }
        let projector = if index.rem_euclid(i128::from(residue.prime.get())) == 0 {
            None
        } else {
            Some(projector(model, &pairing, &generators)?)
        };
        model.allocate(112 + pairing_name.text.len())?;
        Ok(Self {
            pairing_name,
            pairing,
            index,
            gram,
            generators,
            projector,
            gaussian_columns,
            weil_scale: unit,
            trace,
            q,
        })
    }

    pub(crate) fn pair(&self, model: &mut Model, left: Cycle, right: Cycle) -> Result<i128> {
        model.spend(16)?;
        let mut pairing = Wide::<16>::from_i128(0);
        for (a, row) in left.into_iter().zip(self.gram) {
            for (b, g) in right.into_iter().zip(row) {
                pairing = pairing.add(Wide::product(a, b).multiply(g)?)?;
            }
        }
        pairing.signed().ok_or_else(|| {
            Error::plain("Kummer CM product intersection exceeds signed 128-bit range")
        })
    }

    pub(crate) fn hermitian(
        &self,
        model: &mut Model,
        [left, right, diagonal, graph]: Cycle,
    ) -> Result<String> {
        model.spend(6)?;
        let narrow = |value: Wide<16>| {
            value.signed().ok_or_else(|| {
                Error::plain("Kummer CM Hermitian coefficient exceeds signed 128-bit range")
            })
        };
        let a = narrow(
            Wide::from_i128(right)
                .add(Wide::from_i128(diagonal))?
                .add(Wide::product(self.gram[0][3], graph))?,
        )?;
        let d = narrow(
            Wide::from_i128(left)
                .add(Wide::from_i128(diagonal))?
                .add(Wide::from_i128(graph))?,
        )?;
        let b = [narrow(Wide::product(-1, diagonal))?, graph];
        Ok(format!(
            "rvKummerCMCycleHermitian = a={a}, d={d}, b_in_1_omega={b:?}"
        ))
    }

    pub(crate) fn lift(&self, model: &mut Model, s: usize, cycle: Cycle) -> Result<Fixed> {
        let n = self.pairing.modulus(s)?;
        model.spend(24)?;
        let mut coordinates = [0; 6];
        for (coefficient, vector) in cycle.into_iter().zip(self.generators) {
            let scaled =
                super::linking::scaled(vector, coefficient.rem_euclid(i128::from(n)) as u64, n);
            coordinates = super::linking::add_vectors(coordinates, scaled, n);
        }
        self.pairing.require_fixed(model, n, coordinates)?;
        Ok(coordinates)
    }

    pub(crate) fn require_projector(&self) -> Result<&Exterior> {
        self.projector.as_ref().ok_or_else(|| Error::plain("Kummer CM graph Gram is not a local unit; the projectors require a saturated cycle lattice"))
    }

    pub(crate) fn project(
        &self,
        model: &mut Model,
        transcendental: bool,
        s: usize,
        value: Fixed,
    ) -> Result<Fixed> {
        let n = self.pairing.modulus(s)?;
        model.spend(42)?;
        let algebraic = image(self.require_projector()?, value, n);
        let result = if transcendental {
            std::array::from_fn(|i| sub_mod(value[i], algebraic[i], n))
        } else {
            algebraic
        };
        self.pairing.require_fixed(model, n, result)?;
        Ok(result)
    }

    pub(crate) fn defect(&self, model: &mut Model, name: &str, r: usize) -> Result<String> {
        if r == 0 {
            return Err(Error::plain("Kummer CM Frobenius degree must be positive"));
        }
        let a = self.trace / 2;
        let square = u128::try_from(i128::from(self.q) - a * a).map_err(|_| inconsistent())?;
        let b = square.isqrt();
        if self.trace % 2 != 0 || square == 0 || b * b != square {
            return Err(inconsistent());
        }
        let prime = u64::from(self.pairing.ring.prime.get());
        let mut modulus = prime;
        let mut valuation = 0;
        loop {
            let pi = (
                a.rem_euclid(i128::from(modulus)) as u64,
                (b % u128::from(modulus)) as u64,
            );
            let mut power = (1, 0);
            let mut factor = gaussian_product(pi, pi, modulus);
            let mut remaining = r;
            while remaining != 0 {
                model.spend(12)?;
                if remaining & 1 == 1 {
                    power = gaussian_product(power, factor, modulus);
                }
                remaining >>= 1;
                if remaining > 0 {
                    factor = gaussian_product(factor, factor, modulus);
                }
            }
            let norm = super::power_mod(model, self.q, r, modulus)?;
            if power != (norm, 0) {
                break;
            }
            valuation += 1;
            modulus = modulus.checked_mul(prime).ok_or_else(|| {
                Error::plain(
                    "Kummer CM Frobenius defect exceeds the independent local precision range",
                )
            })?;
        }
        Ok(format!("rvKummerCMDefect[{name}, {r}] = elliptic_frobenius=({a}, ±{b}), prime={prime}, transcendental_smith_exponents=[{valuation}, {valuation}]"))
    }

    pub(crate) fn report(&self, name: &str) -> String {
        let gram = self.gram;
        let modulus = self.pairing.ring.modulus;
        let generators = self.generators;
        let det = -4 * self.index * self.index;
        let saturated = self.projector.is_some();
        let (gaussian_columns, weil_scale) = (self.gaussian_columns, self.weil_scale);
        let mut out = format!("rvKummerCMCycleData[{name}] = level={modulus}, generator_index={}, product_gram={gram:?}, product_discriminant={det}, tate_rank=4, locally_saturated={saturated}, gaussian_columns={gaussian_columns:?}, weil_scale={weil_scale}, cycle_lifts={generators:?}", self.index);
        if let Some(projector) = self.projector {
            out.push_str(&format!(", algebraic_projector={projector:?}"));
        }
        out
    }
}

pub(crate) fn cycle_add(left: Cycle, right: Cycle) -> Result<Cycle> {
    cycle_coefficients(|i| left[i].checked_add(right[i]))
}

pub(crate) fn cycle_negate(cycle: Cycle) -> Result<Cycle> {
    cycle_coefficients(|i| cycle[i].checked_neg())
}

fn cycle_coefficients(mut entry: impl FnMut(usize) -> Option<i128>) -> Result<Cycle> {
    let mut result = [0; 4];
    for (i, x) in result.iter_mut().enumerate() {
        *x = entry(i).ok_or_else(|| {
            Error::plain("Kummer CM cycle coefficient exceeds signed 128-bit range")
        })?;
    }
    Ok(result)
}

fn gaussian_product((a, b): (u64, u64), (c, d): (u64, u64), n: u64) -> (u64, u64) {
    (
        sub_mod(mul_mod(a, c, n), mul_mod(b, d, n), n),
        add_mod(mul_mod(a, d, n), mul_mod(b, c, n), n),
    )
}

fn symplectic(x: [u64; 4], y: [u64; 4], modulus: u64) -> u64 {
    (0..2).fold(0, |value, i| {
        add_mod(
            value,
            sub_mod(
                mul_mod(x[i], y[i + 2], modulus),
                mul_mod(x[i + 2], y[i], modulus),
                modulus,
            ),
            modulus,
        )
    })
}

fn wedge(x: [u64; 4], y: [u64; 4], modulus: u64) -> Fixed {
    PAIRS.map(|(i, j)| {
        sub_mod(
            mul_mod(x[i], y[j], modulus),
            mul_mod(x[j], y[i], modulus),
            modulus,
        )
    })
}

fn projector(model: &mut Model, pairing: &Pairing, generators: &[Fixed; 4]) -> Result<Exterior> {
    let ring = pairing.ring;
    let n = ring.modulus;
    let cells = generators
        .iter()
        .flat_map(|&x| generators.iter().map(move |&y| intersection(x, y, n)))
        .collect();
    let gram = Matrix::from_entries(4, cells)?;
    let inverse = hw_algebra::inverse(model, &ring, &gram)?.ok_or_else(inconsistent)?;
    let id = identity();
    let mut result = [[0; 6]; 6];
    model.spend(1300)?;
    for j in 0..6 {
        let u = column(&id, j);
        let covector = generators.map(|v| intersection(v, u, n));
        let coefficients: [u64; 4] = std::array::from_fn(|i| {
            covector.iter().enumerate().fold(0, |sum, (k, &v)| {
                add_mod(sum, mul_mod(*inverse.entry(i, k), v, n), n)
            })
        });
        for (i, row) in result.iter_mut().enumerate() {
            row[j] = generators
                .iter()
                .zip(coefficients)
                .fold(0, |sum, (v, c)| add_mod(sum, mul_mod(v[i], c, n), n));
        }
    }
    if product(&result, &result, n) != result
        || product(&pairing.difference, &result, n) != [[0; 6]; 6]
        || product(&result, &pairing.difference, n) != [[0; 6]; 6]
    {
        return Err(inconsistent());
    }
    for i in 0..6 {
        for j in i..6 {
            if intersection(column(&result, i), column(&id, j), n)
                != intersection(column(&id, i), column(&result, j), n)
            {
                return Err(inconsistent());
            }
        }
    }
    Ok(result)
}

fn index_overflow() -> Error {
    Error::plain("Kummer CM graph index exceeds the intersection arithmetic range")
}

fn inconsistent() -> Error {
    Error::plain(
        "Kummer CM geometry disagrees with its marked intersection or Frobenius certificate",
    )
}
