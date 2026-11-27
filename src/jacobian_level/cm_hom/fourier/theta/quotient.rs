use std::cell::OnceCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use super::{
    alternating_smith, identity, lattice_matrix, multiply, omega, transpose, Kernel, Matrix,
    Subgroup, J,
};
use crate::ec_morphism::Morphism;
use crate::elliptic::Curve;
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::integer_smith::ColumnHermite;
use crate::model::Model;
use crate::rv_ore::gcd;
use crate::rv_qp_quad::{Element, Wide};
use crate::rv_unram::{add_mod, mul_mod};
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod order;

pub(crate) enum Statement {
    Define(Label, Box<ParsedExpr>),
    Data(Label),
    Order(Label),
    Rosati(Box<ParsedExpr>),
    Frobenius(Label),
    Profiles(Label),
    Pair(Box<[ParsedExpr; 2]>),
    Transition(Label, Label),
}

pub(crate) enum Expression {
    Literal(Lattice, Label, Vector),
    Map(Map, Label, Box<ParsedExpr>),
    Polarize(Box<ParsedExpr>),
    Endomorphism(Label, Matrix),
    EndomorphismBasis(Label, usize),
    EndomorphismIdentity(Label),
    EndomorphismFrobenius(Label),
    EndomorphismRosati(Box<ParsedExpr>),
    EndomorphismBinary(bool, Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Copy)]
pub(crate) enum Map {
    Include,
    Pullback,
    Transition,
    TransitionPullback,
    TorsionLift,
}

#[derive(Clone, Copy)]
pub(crate) enum Lattice {
    Base,
    BaseDual,
    Quotient,
    QuotientDual,
}

pub(crate) type Vector = [i128; 4];

pub(crate) struct Quotient {
    pub(crate) kernel_name: Label,
    kernel: Rc<Kernel>,
    subgroup: Rc<Subgroup>,
    numerator: Matrix,
    denominator: i128,
    inclusion: Matrix,
    form: Matrix,
    divisors: [i128; 2],
    omega_numerator: Matrix,
    conductor: u64,
    centralizer: OnceCell<[order::Endomorphism; 8]>,
}

struct Frobenius {
    matrix: Matrix,
    coefficients: [i128; 2],
    tangent_omega: Vec<u64>,
    field_size: u64,
    trace: i128,
}

impl Quotient {
    pub(crate) fn new(
        model: &mut Model,
        kernel_name: Label,
        kernel: Rc<Kernel>,
        subgroup: Rc<Subgroup>,
    ) -> Result<Self> {
        let denominator = i128::from(kernel.divisors[1]);
        let form = multiply(J, kernel.matrix)?;
        let mut hermite = ColumnHermite::from_diagonal(model, &[denominator; 4])?;
        let n = denominator as u64;
        let factors = kernel.moduli();
        for &point in subgroup.members() {
            let generator = torsion_lift(&kernel, point, n, factors);
            hermite.insert(model, generator.map(i128::from).to_vec())?;
        }
        hermite.canonicalize(model)?;
        let columns = hermite.columns();
        let numerator = std::array::from_fn(|i| std::array::from_fn(|j| columns[j][i]));
        let inclusion = triangular_inverse(numerator, denominator)?;
        if multiply(numerator, inclusion)? != diagonal(denominator)
            || determinant_triangular(inclusion)? != subgroup.members().len() as i128
        {
            return Err(inconsistent());
        }
        let form_numerator = multiply(multiply(transpose(numerator), form)?, numerator)?;
        let square = denominator
            .checked_mul(denominator)
            .ok_or_else(arithmetic)?;
        let descended = divide_matrix(form_numerator, square).ok_or_else(|| {
            Error::plain("CM polarized quotient requires an isotropic theta subgroup")
        })?;
        if multiply(multiply(transpose(inclusion), descended)?, inclusion)? != form {
            return Err(inconsistent());
        }
        let (elementary, _) = alternating_smith(model, descended)?;
        let divisors = [elementary[0][1], elementary[2][3]];
        let degree = subgroup.members().len() as i128;
        if divisors[0]
            .checked_mul(divisors[1])
            .and_then(|d| d.checked_mul(degree))
            != Some(kernel.dimension)
        {
            return Err(inconsistent());
        }
        let omega = omega_matrix(model, &kernel)?;
        let omega_numerator = multiply(multiply(inclusion, omega)?, numerator)?;
        let content = omega_numerator
            .iter()
            .flatten()
            .fold(denominator as u128, |d, c| gcd(d, c.unsigned_abs()));
        let conductor = n / content as u64;
        let omega_numerator =
            divide_matrix(omega_numerator, content as i128).ok_or_else(inconsistent)?;
        model.allocate(88 + kernel_name.text.len())?;
        Ok(Self {
            kernel_name,
            kernel,
            subgroup,
            numerator,
            denominator,
            inclusion,
            form: descended,
            divisors,
            omega_numerator,
            conductor,
            centralizer: OnceCell::new(),
        })
    }

    pub(crate) fn report(&self) -> String {
        format!(
            "rvCMQuotientData = basis_numerator={:?}, basis_denominator={}, inclusion={:?}, degree={}, riemann_form={:?}, type={:?}, principal={}, diagonal_cm_numerator={:?}, diagonal_cm_conductor={}",
            self.numerator, self.denominator, self.inclusion, self.subgroup.members().len(), self.form,
            self.divisors, self.divisors == [1, 1], self.omega_numerator, self.conductor
        )
    }

    pub(crate) fn include(&self, vector: Vector) -> Result<Vector> {
        apply_matrix(self.inclusion, vector)
    }

    pub(crate) fn pullback(&self, vector: Vector) -> Result<Vector> {
        apply_matrix(transpose(self.inclusion), vector)
    }

    pub(crate) fn pair(&self, x: Vector, y: Vector) -> Result<i128> {
        pairing(self.form, x, y)
    }

    pub(crate) fn polarize(&self, x: Vector) -> Result<Vector> {
        apply_matrix(self.form, x)
    }

    fn transition(&self, target: &Self) -> Result<Matrix> {
        if self.kernel_name.text != target.kernel_name.text {
            return Err(Error::plain(
                "CM quotient transition requires the same named theta kernel",
            ));
        }
        if self
            .subgroup
            .members()
            .iter()
            .any(|x| target.subgroup.members().binary_search(x).is_err())
        {
            return Err(Error::plain(
                "CM quotient transition requires a nested subgroup",
            ));
        }
        let matrix = multiply(target.inclusion, self.numerator)?;
        let matrix = divide_matrix(matrix, self.denominator).ok_or_else(inconsistent)?;
        if multiply(matrix, self.inclusion)? != target.inclusion
            || multiply(multiply(transpose(matrix), target.form)?, matrix)? != self.form
        {
            return Err(inconsistent());
        }
        Ok(matrix)
    }

    pub(crate) fn move_to(&self, target: &Self, x: Vector) -> Result<Vector> {
        apply_matrix(self.transition(target)?, x)
    }

    pub(crate) fn pullback_transition(&self, target: &Self, x: Vector) -> Result<Vector> {
        apply_matrix(transpose(self.transition(target)?), x)
    }

    pub(crate) fn lift(&self, point: super::Point) -> Result<Vector> {
        if self.subgroup.members().binary_search(&point).is_err() {
            return Err(Error::plain(
                "theta point does not belong to the quotient subgroup",
            ));
        }
        let n = self.denominator as u64;
        let vector = torsion_lift(&self.kernel, point, n, self.kernel.moduli());
        let scaled = apply_matrix(self.inclusion, vector.map(i128::from))?;
        if scaled.iter().any(|x| x % self.denominator != 0) {
            return Err(inconsistent());
        }
        Ok(scaled.map(|x| x / self.denominator))
    }

    pub(crate) fn transition_report(&self, target: &Self) -> Result<String> {
        let matrix = self.transition(target)?;
        let degree = target.subgroup.members().len() / self.subgroup.members().len();
        Ok(format!(
            "rvCMQuotientTransitionData = matrix={matrix:?}, dual={:?}, degree={degree}",
            transpose(matrix)
        ))
    }

    fn certified_frobenius(&self, model: &mut Model) -> Result<(Frobenius, usize, Matrix, i128)> {
        let p = self.kernel.characteristic;
        if (self.subgroup.members().len() as u64).is_multiple_of(p) {
            return Err(Error::plain(
                "CM quotient Frobenius requires a prime-to-characteristic subgroup",
            ));
        }
        let frobenius = Frobenius::new(model, &self.kernel)?;
        let degree = frobenius.orbit_degree(model, &self.kernel, &self.subgroup)?;
        let power = matrix_power(frobenius.matrix, degree)?;
        let image = multiply(multiply(self.inclusion, power)?, self.numerator)?;
        let descended = divide_matrix(image, self.denominator).ok_or_else(inconsistent)?;
        let q_power = i128::from(frobenius.field_size)
            .checked_pow(degree.try_into().map_err(|_| arithmetic())?)
            .ok_or_else(arithmetic)?;
        if multiply(multiply(transpose(descended), self.form)?, descended)?
            != multiply(diagonal(q_power), self.form)?
        {
            return Err(inconsistent());
        }
        Ok((frobenius, degree, descended, q_power))
    }

    pub(crate) fn frobenius_endomorphism(&self, model: &mut Model) -> Result<order::Endomorphism> {
        let (_, _, matrix, _) = self.certified_frobenius(model)?;
        self.endomorphism(matrix)
    }

    pub(crate) fn frobenius_report(&self, model: &mut Model) -> Result<String> {
        let (frobenius, degree, descended, q_power) = self.certified_frobenius(model)?;
        Ok(format!(
            "rvCMQuotientFrob = field_size={}, trace={}, frobenius_anchor={:?}, omega_differential={:?}, degree={degree}, base_defined={}, induced={descended:?}, multiplier={q_power}",
            frobenius.field_size, frobenius.trace, frobenius.coefficients, frobenius.tangent_omega, degree == 1
        ))
    }
}

pub(crate) fn profiles(model: &mut Model, name: Label, kernel: Rc<Kernel>) -> Result<String> {
    if !kernel.etale() {
        return Err(Error::plain(
            "CM quotient orbit profiles require a prime-to-characteristic polarization degree",
        ));
    }
    let frobenius = Frobenius::new(model, &kernel)?;
    let mut conductors = BTreeMap::new();
    let mut counts = BTreeMap::new();
    let mut assignments = Vec::new();
    for subgroup in kernel.lagrangians(model)? {
        let entry = Quotient::new(model, name.clone(), kernel.clone(), subgroup.clone())?;
        let degree = frobenius.orbit_degree(model, &kernel, subgroup)?;
        *conductors.entry(entry.conductor).or_insert(0_usize) += 1;
        *counts.entry(degree).or_insert(0_usize) += 1;
        assignments.push([entry.conductor, degree as u64]);
    }
    for (&degree, count) in &mut counts {
        if *count % degree != 0 {
            return Err(inconsistent());
        }
        *count /= degree;
    }
    Ok(format!(
        "rvCMQuotientProfiles = count={}, diagonal_cm_conductors={conductors:?}, frobenius_orbits={counts:?}, frobenius_anchor={:?}, conductor_and_degree={assignments:?}",
        assignments.len(), frobenius.coefficients
    ))
}

impl Frobenius {
    fn new(model: &mut Model, kernel: &Kernel) -> Result<Self> {
        let surface = &kernel.surface;
        let hom = &surface.frames[2];
        let q = hom.q;
        let trace = hom.trace;
        let fundamental = hom.fundamental;
        let field = hom.endpoints[0].field();
        let f_dual = differential(model, &surface.frames[1].maps[0], &hom.endpoints[1])?;
        let g = differential(model, &hom.maps[1], &hom.endpoints[0])?;
        let tangent_alpha = field.multiply(model, &f_dual, &g)?;
        let parity = fundamental.rem_euclid(2);
        let offset = (hom.form[1] - parity) / 2;
        let tangent_omega = field.subtract(model, &tangent_alpha, &field.signed(offset))?;
        let conductor = i128::from(hom.conductor);
        let mut candidates = Vec::with_capacity(2);
        for b in [conductor, -conductor] {
            let a = (trace - b * parity) / 2;
            let weighted = field.multiply(model, &field.signed(b), &tangent_omega)?;
            if field.add(model, &field.signed(a), &weighted)? == field.zero() {
                candidates.push([a, b]);
            }
        }
        let [[a, b]] = candidates.as_slice() else {
            return Err(Error::plain(
                "actual CM differential does not select a unique ordinary Frobenius orientation",
            ));
        };
        let omega = omega_matrix(model, kernel)?;
        let matrix = sum(diagonal(*a), multiply(diagonal(*b), omega)?)?;
        let residual = sum(
            sum(
                multiply(matrix, matrix)?,
                multiply(diagonal(-trace), matrix)?,
            )?,
            diagonal(i128::from(q)),
        )?;
        let form = multiply(J, kernel.matrix)?;
        if residual != [[0; 4]; 4]
            || multiply(matrix, kernel.matrix)? != multiply(kernel.matrix, matrix)?
            || multiply(multiply(transpose(matrix), form)?, matrix)?
                != multiply(diagonal(i128::from(q)), form)?
        {
            return Err(inconsistent());
        }
        Ok(Self {
            matrix,
            coefficients: [*a, *b],
            tangent_omega,
            field_size: q,
            trace,
        })
    }

    fn orbit_degree(
        &self,
        model: &mut Model,
        kernel: &Kernel,
        subgroup: &Subgroup,
    ) -> Result<usize> {
        let mut points = subgroup.members().to_vec();
        let n = kernel.divisors[1];
        let moduli = kernel.moduli();
        let matrix = multiply(
            multiply(kernel.to_symplectic, self.matrix)?,
            kernel.from_symplectic,
        )?;
        for degree in 1.. {
            model.spend(16 * points.len())?;
            for point in &mut points {
                let image = modular_vector(matrix, *point, n);
                *point = std::array::from_fn(|i| image[i] % moduli[i]);
            }
            points.sort_unstable();
            if points.windows(2).any(|p| p[0] == p[1]) {
                return Err(inconsistent());
            }
            if points == subgroup.members() {
                return Ok(degree);
            }
        }
        unreachable!()
    }
}

pub(crate) fn pair_base(kernel: &Kernel, x: Vector, y: Vector) -> Result<i128> {
    pairing(multiply(J, kernel.matrix)?, x, y)
}

pub(crate) fn polarize_base(kernel: &Kernel, x: Vector) -> Result<Vector> {
    apply_matrix(multiply(J, kernel.matrix)?, x)
}

fn omega_matrix(model: &mut Model, kernel: &Kernel) -> Result<Matrix> {
    let w = omega(&kernel.surface)?;
    lattice_matrix(
        model,
        &kernel.surface,
        [[w, Element::ZERO], [Element::ZERO, w]],
    )
}

fn differential(model: &mut Model, map: &Morphism, curve: &Curve) -> Result<Vec<u64>> {
    let field = curve.field();
    let Some(coordinates) = &map.coordinates else {
        return Ok(field.zero());
    };
    let x = &coordinates.x;
    if x.numerator.len() != x.denominator.len() + 1 {
        return Ok(field.zero());
    }
    let y = &coordinates.y_multiplier;
    // The invariant tangent parameter at infinity is x/y.
    let numerator = field.multiply(
        model,
        x.numerator.last().ok_or_else(inconsistent)?,
        y.denominator.last().ok_or_else(inconsistent)?,
    )?;
    let denominator = field.multiply(
        model,
        x.denominator.last().ok_or_else(inconsistent)?,
        y.numerator.last().ok_or_else(inconsistent)?,
    )?;
    let inverse = field.inverse_unit(model, &denominator)?;
    field.multiply(model, &numerator, &inverse)
}

fn triangular_inverse(matrix: Matrix, factor: i128) -> Result<Matrix> {
    let mut columns = [[0; 4]; 4];
    for (column, output) in columns.iter_mut().enumerate() {
        for (i, row) in matrix.iter().enumerate() {
            let rhs = if i == column { factor } else { 0 };
            let mut value = Wide::<8>::from_i128(rhs);
            for (&c, &x) in row[..i].iter().zip(&output[..i]) {
                value = value.add(Wide::product(c, x).negate())?;
            }
            let value = value.signed().ok_or_else(arithmetic)?;
            if row[i] <= 0 || value % row[i] != 0 {
                return Err(inconsistent());
            }
            output[i] = value / row[i];
        }
    }
    Ok(transpose(columns))
}

fn torsion_lift(
    kernel: &Kernel,
    point: super::Point,
    n: u64,
    factors: super::Point,
) -> super::Point {
    let scaled = std::array::from_fn(|i| (n / factors[i]) * point[i]);
    modular_vector(kernel.symplectic_basis, scaled, n)
}

fn modular_vector(matrix: Matrix, x: [u64; 4], n: u64) -> [u64; 4] {
    matrix.map(|row| {
        row.into_iter().zip(x).fold(0, |sum, (a, b)| {
            let a = a.rem_euclid(i128::from(n)) as u64;
            add_mod(sum, mul_mod(a, b, n), n)
        })
    })
}

fn apply_matrix(matrix: Matrix, vector: Vector) -> Result<Vector> {
    let mut result = [0; 4];
    for (entry, row) in result.iter_mut().zip(matrix) {
        *entry = super::dot(row, vector)?;
    }
    Ok(result)
}

fn pairing(form: Matrix, x: Vector, y: Vector) -> Result<i128> {
    let mut value = Wide::<16>::from_i128(0);
    for (row, coefficient) in form.into_iter().zip(x) {
        for (entry, right) in row.into_iter().zip(y) {
            value = value.add(Wide::product(coefficient, entry).multiply(right)?)?;
        }
    }
    value.signed().ok_or_else(arithmetic)
}

fn diagonal(n: i128) -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| if i == j { n } else { 0 }))
}

fn determinant_triangular(matrix: Matrix) -> Result<i128> {
    matrix
        .iter()
        .enumerate()
        .try_fold(1_i128, |result, (i, row)| {
            result.checked_mul(row[i]).ok_or_else(arithmetic)
        })
}

fn sum(mut a: Matrix, b: Matrix) -> Result<Matrix> {
    for (row, rhs) in a.iter_mut().zip(b) {
        for (entry, value) in row.iter_mut().zip(rhs) {
            *entry = entry.checked_add(value).ok_or_else(arithmetic)?;
        }
    }
    Ok(a)
}

fn divide_matrix(mut matrix: Matrix, divisor: i128) -> Option<Matrix> {
    for entry in matrix.iter_mut().flatten() {
        if *entry % divisor != 0 {
            return None;
        }
        *entry /= divisor;
    }
    Some(matrix)
}

fn matrix_power(mut matrix: Matrix, mut exponent: usize) -> Result<Matrix> {
    let mut result = identity();
    while exponent != 0 {
        if exponent % 2 == 1 {
            result = multiply(result, matrix)?;
        }
        exponent /= 2;
        if exponent != 0 {
            matrix = multiply(matrix, matrix)?;
        }
    }
    Ok(result)
}

fn arithmetic() -> Error {
    Error::plain("CM quotient lattice arithmetic exceeds the supported exact integer range")
}

fn inconsistent() -> Error {
    Error::plain("CM polarized quotient lattice certificate is inconsistent")
}

pub(crate) fn add(mut a: Vector, b: Vector) -> Result<Vector> {
    for (x, y) in a.iter_mut().zip(b) {
        *x = x.checked_add(y).ok_or_else(arithmetic)?;
    }
    Ok(a)
}

pub(crate) fn negate(mut value: Vector) -> Result<Vector> {
    for x in &mut value {
        *x = x.checked_neg().ok_or_else(arithmetic)?;
    }
    Ok(value)
}
