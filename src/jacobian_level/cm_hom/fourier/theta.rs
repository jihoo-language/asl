use std::cell::OnceCell;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use super::{Abelian, Divisor};
use crate::elliptic_cm_recognition as elliptic_maps;
use crate::elliptic_gluing::Factor;
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hyperelliptic::Residues;
use crate::integer_smith::{self, Track};
use crate::jacobian_level::cm_hom::marked;
use crate::model::Model;
use crate::rv_ore::gcd;
use crate::rv_qp_quad::{Element, Wide};
use crate::rv_unram::{add_mod, mul_mod, sub_mod};
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod metaplectic;
pub(crate) mod quotient;

pub(crate) enum Statement {
    Quotient(Box<quotient::Statement>),
    Metaplectic(metaplectic::Statement),
    Define(Label, Box<ParsedExpr>),
    Named(NamedReport, Label),
    Element(ElementReport, Box<ParsedExpr>),
    Pair(Box<[ParsedExpr; 2]>),
    Action([u64; 2], Box<ParsedExpr>),
    Subgroup(Box<ParsedExpr>),
    Torsion(Label, [Label; 2]),
}

pub(crate) enum NamedReport {
    Data,
    Lagrangians,
    Elements,
}

pub(crate) enum ElementReport {
    Lift,
    Order,
}

pub(crate) enum Expression {
    Quotient(Box<quotient::Expression>),
    Metaplectic(Box<metaplectic::Expression>),
    Element(Label, [i128; 4], bool),
    Operator(u64, Box<ParsedExpr>),
    Multiply(Box<[ParsedExpr; 2]>),
    Subgroup(Label, Vec<ParsedExpr>),
    Perpendicular(Box<ParsedExpr>),
    Lagrangian(Label, usize),
}

pub(crate) type Point = [u64; 4];
type Matrix = [[i128; 4]; 4];
const J: Matrix = [[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]];
const GROUP_LIMIT: usize = 4096;

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Operator {
    point: Point,
    phase: u64,
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Subgroup {
    members: Rc<[Point]>,
}

pub(crate) struct Kernel {
    surface: Rc<Abelian>,
    divisor: Divisor,
    matrix: Matrix,
    symplectic_basis: Matrix,
    to_symplectic: Matrix,
    from_symplectic: Matrix,
    divisors: [u64; 2],
    dimension: i128,
    characteristic: u64,
    lagrangians: OnceCell<Vec<Rc<Subgroup>>>,
}

impl Kernel {
    pub(crate) fn new(model: &mut Model, surface: Rc<Abelian>, divisor: Divisor) -> Result<Self> {
        let a = as_integer(divisor.a)?;
        as_integer(divisor.d)?;
        let coordinates = surface
            .field
            .coordinates(model, surface.bases[2], divisor.z)?;
        for coordinate in coordinates {
            as_integer(coordinate)?;
        }
        let dimension = as_integer(surface.determinant(model, divisor)?)?;
        if a <= 0 || dimension <= 0 {
            return Err(Error::plain(
                "CM theta kernel requires a positive definite divisor",
            ));
        }
        let hermitian = surface.matrix(model, divisor)?;
        let matrix = lattice_matrix(model, &surface, hermitian)?;
        let form = multiply(J, matrix)?;
        for (i, row) in form.iter().enumerate() {
            if row[i] != 0
                || row
                    .iter()
                    .enumerate()
                    .any(|(j, &x)| x.checked_neg() != Some(form[j][i]))
            {
                return Err(inconsistent());
            }
        }
        let pfaffian = Wide::<8>::product(form[0][1], form[2][3])
            .add(Wide::product(form[0][2], form[1][3]).negate())?
            .add(Wide::product(form[0][3], form[1][2]))?;
        if pfaffian.magnitude != Wide::<8>::from_i128(dimension).magnitude {
            return Err(inconsistent());
        }
        let (reduced, basis) = alternating_smith(model, form)?;
        let [a, b] = [reduced[0][1], reduced[2][3]];
        if a <= 0 || b <= 0 || b % a != 0 || a.checked_mul(b) != Some(dimension) {
            return Err(inconsistent());
        }
        let mut smith = integer_smith::Matrix::new(model, 4, 4)?;
        for (i, row) in matrix.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                smith.write(model, i, j, value)?;
            }
        }
        if integer_smith::reduce(model, smith, Track::None, Track::None)?.factors != [a, a, b, b] {
            return Err(inconsistent());
        }
        let [a64, b64] = [a, b].map(|x| {
            u64::try_from(x).map_err(|_| {
                Error::plain(
                    "CM theta kernel exponent exceeds the supported unsigned integer range",
                )
            })
        });
        let divisors = [a64?, b64?];
        // If M^(-1)x = T c, then T^t J x = diag(a J_2,b J_2)c.
        let negative_j = J.map(|row| row.map(|x| -x));
        let to_symplectic = multiply(multiply(negative_j, transpose(basis))?, J)?;
        let mut from_symplectic = multiply(matrix, basis)?;
        for row in &mut from_symplectic {
            for (value, divisor) in row.iter_mut().zip([a, a, b, b]) {
                if *value % divisor != 0 {
                    return Err(inconsistent());
                }
                *value /= divisor;
            }
        }
        if multiply(to_symplectic, from_symplectic)? != identity() {
            return Err(inconsistent());
        }
        model.allocate(80)?;
        let characteristic = surface.frames[2].endpoints[0].ring.prime();
        Ok(Self {
            surface,
            divisor,
            matrix,
            symplectic_basis: basis,
            to_symplectic,
            from_symplectic,
            divisors,
            dimension,
            characteristic,
            lagrangians: OnceCell::new(),
        })
    }

    pub(crate) fn cardinality_bound(&self) -> usize {
        usize::try_from(self.dimension)
            .unwrap_or(usize::MAX)
            .saturating_pow(2)
    }

    fn moduli(&self) -> Point {
        let [a, b] = self.divisors;
        [a, a, b, b]
    }

    fn etale(&self) -> bool {
        self.dimension % i128::from(self.characteristic) != 0
    }

    pub(crate) fn report(&self) -> String {
        let [d1, d2] = self.divisors;
        let size = Wide::<8>::product(self.dimension, self.dimension);
        // (d1*d2)^2 may exceed i128. Its factorization is an exact size in that case.
        let size = size
            .signed()
            .map_or_else(|| format!("{}^2", self.dimension), |n| n.to_string());
        format!(
            "rvCMThetaData = matrix={:?}, smith=[{d1}, {d1}, {d2}, {d2}], type={:?}, dimension={}, kernel_size={size}, characteristic={}, geometric_etale={}",
            self.matrix, self.divisors, self.dimension, self.characteristic, self.etale()
        )
    }

    pub(crate) fn torsion_report(&self, model: &mut Model, names: [Label; 2]) -> Result<String> {
        let elliptic = [
            model.jacobian_level(&names[0])?,
            model.jacobian_level(&names[1])?,
        ];
        let n = elliptic[0].order;
        if elliptic
            .iter()
            .any(|level| level.len() != 2 || level.order != n)
        {
            return Err(Error::plain("CM theta torsion check needs two full elliptic bases at the same prime-power level"));
        }
        let mut factors = Vec::with_capacity(2);
        for (basis, label) in elliptic.iter().zip(names) {
            let curve = model.hww_curve(&Label {
                text: basis.curve_name.to_string(),
                span: label.span,
            })?;
            factors.push(Factor {
                name: basis.curve_name.clone(),
                curve,
            });
        }
        let models = [
            marked::short_model(model, &factors[0])?,
            marked::short_model(model, &factors[1])?,
        ];
        let surface = &self.surface;
        let old_models = surface.frames[2].endpoints.each_ref().map(Rc::as_ref);
        let hermitian = surface.matrix(model, self.divisor)?;
        let endomorphism = surface.make_endomorphism(model, hermitian)?;
        let mut actual = [[0; 4]; 4];
        for i in 0..2 {
            for j in 0..2 {
                let map = surface.block(model, endomorphism, 2 * i + j)?;
                let map = elliptic_maps::extend_over(
                    model,
                    &map,
                    [factors[j].name.clone(), factors[i].name.clone()],
                    [old_models[j], old_models[i]],
                    [&models[j], &models[i]],
                )?;
                for (k, point) in elliptic[j].generators().iter().enumerate() {
                    let image =
                        marked::apply(model, &map, [&models[j], &models[i]], point, &factors[i])?;
                    let coordinates = elliptic[i].checked_coordinates(model, &image)?;
                    for (l, coordinate) in coordinates.into_iter().enumerate() {
                        actual[2 * i + l][2 * j + k] = coordinate;
                    }
                }
            }
        }
        // The two primitive Weil roots may differ. Express the right pairing
        // in units of the left root before inspecting the geometric alternating form.
        let residues = Residues::new(elliptic[0].prime, elliptic[0].height)?;
        let left_in_right = elliptic[1].polarization_scale(model, &elliptic[0], 1)?;
        let right_in_left = residues.inverse_unit(model, &left_in_right)?;
        let actual_form: [Point; 4] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let row = i ^ 1;
                let coefficient = if i < 2 { 1 } else { right_in_left };
                let entry = mul_mod(coefficient, actual[row][j], n);
                if i.is_multiple_of(2) {
                    entry
                } else {
                    sub_mod(0, entry, n)
                }
            })
        });
        for (i, row) in actual_form.iter().enumerate() {
            if row[i] != 0
                || row
                    .iter()
                    .enumerate()
                    .any(|(j, &x)| add_mod(x, actual_form[j][i], n) != 0)
            {
                return Err(inconsistent());
            }
        }
        // The cokernel of [actual | n I] is the cokernel of the level map,
        // with the same elementary divisors as its kernel over Z/n.
        let mut augmented = integer_smith::Matrix::new(model, 4, 8)?;
        for (i, row) in actual.iter().enumerate() {
            for (j, &entry) in row.iter().enumerate() {
                augmented.write(model, i, j, i128::from(entry))?;
            }
            augmented.write(model, i, i + 4, i128::from(n))?;
        }
        let smith = integer_smith::reduce(model, augmented, Track::None, Track::None)?.factors;
        let expected = self
            .moduli()
            .map(|d| gcd(u128::from(d), u128::from(n)) as i128);
        if smith != expected {
            return Err(inconsistent());
        }
        Ok(format!("rvCMThetaTorsion = level={n}, matrix={actual:?}, kernel_smith={smith:?}, right_pairing_scale={right_in_left}"))
    }

    pub(crate) fn point(&self, x: [i128; 4], lattice: bool) -> Point {
        let moduli = self.moduli();
        std::array::from_fn(|i| {
            let n = moduli[i];
            let residue = |v: i128| v.rem_euclid(i128::from(n)) as u64;
            if lattice {
                self.to_symplectic[i]
                    .into_iter()
                    .zip(x)
                    .fold(0, |sum, (a, b)| {
                        add_mod(sum, mul_mod(residue(a), residue(b), n), n)
                    })
            } else {
                residue(x[i])
            }
        })
    }

    pub(crate) fn lift(&self, x: Point) -> Result<[i128; 4]> {
        let mut answer = [0; 4];
        for (out, row) in answer.iter_mut().zip(self.from_symplectic) {
            *out = dot(row, x.map(i128::from))?;
        }
        Ok(answer)
    }

    pub(crate) fn add(&self, x: Point, y: Point) -> Point {
        let moduli = self.moduli();
        std::array::from_fn(|i| add_mod(x[i], y[i], moduli[i]))
    }

    pub(crate) fn negate(&self, x: Point) -> Point {
        let moduli = self.moduli();
        std::array::from_fn(|i| sub_mod(0, x[i], moduli[i]))
    }

    pub(crate) fn order(&self, x: Point) -> u64 {
        x.into_iter().zip(self.moduli()).fold(1, |order, (x, n)| {
            let m = n / gcd(u128::from(n), u128::from(x)) as u64;
            lcm(order, m)
        })
    }

    fn clock(&self, point: Point, position: [u64; 2]) -> u64 {
        let [a, b] = self.divisors;
        let low = (b / a) * mul_mod(point[1], position[0], a);
        let high = mul_mod(point[3], position[1], b);
        add_mod(low, high, b)
    }

    fn pairing(&self, x: Point, y: Point) -> u64 {
        sub_mod(
            self.clock(y, [x[0], x[2]]),
            self.clock(x, [y[0], y[2]]),
            self.divisors[1],
        )
    }

    pub(crate) fn pair(&self, x: Point, y: Point) -> String {
        let n = self.divisors[1];
        let numerator = self.pairing(x, y);
        let common = gcd(u128::from(numerator), u128::from(n)) as u64;
        format!("{}/{}", numerator / common, n / common)
    }

    pub(crate) fn operator(&self, point: Point, phase: u64) -> Operator {
        Operator {
            point,
            phase: phase % self.divisors[1],
        }
    }

    pub(crate) fn multiply(&self, a: Operator, b: Operator) -> Operator {
        let n = self.divisors[1];
        let shift = self.clock(a.point, [b.point[0], b.point[2]]);
        Operator {
            point: self.add(a.point, b.point),
            phase: sub_mod(add_mod(a.phase, b.phase, n), shift, n),
        }
    }

    pub(crate) fn action(&self, operator: Operator, [x, y]: [u64; 2]) -> Result<String> {
        let [a, b] = self.divisors;
        if x >= a || y >= b {
            return Err(Error::plain(
                "CM theta position is outside the Schrödinger basis",
            ));
        }
        let phase = sub_mod(operator.phase, self.clock(operator.point, [x, y]), b);
        let position = [
            add_mod(x, operator.point[0], a),
            add_mod(y, operator.point[2], b),
        ];
        Ok(format!(
            "rvCMThetaAction = position={position:?}, phase={phase}, root_order={b}"
        ))
    }

    pub(crate) fn elements(&self, model: &mut Model) -> Result<Vec<Point>> {
        let n = self.cardinality_bound();
        if n > GROUP_LIMIT {
            return Err(Error::plain(
                "CM theta kernel exceeds the 4096-element enumeration limit",
            ));
        }
        model.spend(n)?;
        let moduli = self.moduli();
        Ok((0..n)
            .map(|mut i| {
                moduli.map(|d| {
                    let x = i as u64 % d;
                    i /= d as usize;
                    x
                })
            })
            .collect())
    }

    fn adjoin(&self, model: &mut Model, subgroup: &[Point], x: Point) -> Result<Vec<Point>> {
        let mut members = subgroup.to_vec();
        let mut t = x;
        while subgroup.binary_search(&t).is_err() {
            model.spend(subgroup.len())?;
            members.extend(subgroup.iter().map(|&y| self.add(t, y)));
            t = self.add(t, x);
        }
        members.sort_unstable();
        Ok(members)
    }

    pub(crate) fn subgroup(&self, model: &mut Model, generators: &[Point]) -> Result<Subgroup> {
        // Each generated subgroup stays in the same bounded finite ambient group.
        if self.cardinality_bound() > GROUP_LIMIT {
            return Err(Error::plain(
                "CM theta kernel exceeds the 4096-element subgroup limit",
            ));
        }
        let mut members = vec![[0; 4]];
        for &x in generators {
            model.spend(1)?;
            if members.binary_search(&x).is_err() {
                members = self.adjoin(model, &members, x)?;
            }
        }
        model.allocate(4 * members.len())?;
        Ok(Subgroup {
            members: members.into(),
        })
    }

    pub(crate) fn perpendicular(&self, model: &mut Model, subgroup: &Subgroup) -> Result<Subgroup> {
        let mut members = Vec::new();
        for x in self.elements(model)? {
            if self.orthogonal(model, &subgroup.members, x)? {
                members.push(x);
            }
        }
        members.sort_unstable();
        model.allocate(4 * members.len())?;
        Ok(Subgroup {
            members: members.into(),
        })
    }

    fn orthogonal(&self, model: &mut Model, points: &[Point], x: Point) -> Result<bool> {
        for &y in points {
            model.spend(1)?;
            if self.pairing(x, y) != 0 {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn exponent(&self, subgroup: &Subgroup) -> u64 {
        subgroup
            .members
            .iter()
            .fold(1, |order, &x| lcm(order, self.order(x)))
    }

    pub(crate) fn subgroup_report(&self, model: &mut Model, subgroup: &Subgroup) -> Result<String> {
        let perp = self.perpendicular(model, subgroup)?;
        let isotropic = subgroup
            .members
            .iter()
            .all(|x| perp.members.binary_search(x).is_ok());
        let maximal = perp == *subgroup;
        let size = subgroup.members.len();
        let descent = if isotropic && self.etale() {
            size.to_string()
        } else {
            "unavailable".to_string()
        };
        Ok(format!(
            "rvCMThetaSubgroupData = size={size}, exponent={}, isotropic={isotropic}, perpendicular_size={}, maximal_isotropic={maximal}, geometric_descent_degree={descent}",
            self.exponent(subgroup), perp.members.len()
        ))
    }

    fn lagrangians(&self, model: &mut Model) -> Result<&[Rc<Subgroup>]> {
        if let Some(known) = self.lagrangians.get() {
            return Ok(known);
        }
        let group = self.elements(model)?;
        let zero: Rc<[Point]> = vec![[0; 4]].into();
        let mut seen = BTreeSet::from([zero.clone()]);
        let mut pending = vec![zero];
        let mut lagrangians = Vec::new();
        while let Some(points) = pending.pop() {
            if points.len() as i128 == self.dimension {
                let subgroup = Subgroup { members: points };
                if self.perpendicular(model, &subgroup)? != subgroup {
                    return Err(inconsistent());
                }
                lagrangians.push(Rc::new(subgroup));
                continue;
            }
            for &x in &group {
                model.spend(1)?;
                if points.binary_search(&x).is_ok() || !self.orthogonal(model, &points, x)? {
                    continue;
                }
                let larger: Rc<[Point]> = self.adjoin(model, &points, x)?.into();
                if seen.insert(larger.clone()) {
                    model.allocate(4 * larger.len())?;
                    pending.push(larger);
                }
            }
        }
        lagrangians.sort_unstable();
        Ok(self.lagrangians.get_or_init(|| lagrangians))
    }

    pub(crate) fn lagrangian(&self, model: &mut Model, i: usize) -> Result<Rc<Subgroup>> {
        self.lagrangians(model)?
            .get(i)
            .cloned()
            .ok_or_else(|| Error::plain("CM theta Lagrangian index is out of range"))
    }

    pub(crate) fn lagrangian_report(&self, model: &mut Model) -> Result<String> {
        let lagrangians = self.lagrangians(model)?;
        let mut exponents = BTreeMap::new();
        for subgroup in lagrangians {
            *exponents.entry(self.exponent(subgroup)).or_insert(0) += 1;
        }
        Ok(format!(
            "rvCMThetaLagrangians = count={}, by_exponent={exponents:?}",
            lagrangians.len()
        ))
    }
}

impl Operator {
    pub(crate) fn coordinates(self) -> (Point, u64) {
        (self.point, self.phase)
    }
}

impl Subgroup {
    pub(crate) fn members(&self) -> &[Point] {
        &self.members
    }
}

fn omega(surface: &Abelian) -> Result<Element> {
    if surface.frames[2].fundamental % 4 == 0 {
        Element::new(0, 1, 1)
    } else {
        Element::new(1, 1, 2)
    }
}

fn lattice_matrix(model: &mut Model, surface: &Abelian, matrix: super::Matrix) -> Result<Matrix> {
    let bases = [[Element::ONE, omega(surface)?], surface.bases[2]];
    let mut output_matrix = [[0; 4]; 4];
    for (i, row) in matrix.into_iter().enumerate() {
        for (j, entry) in row.into_iter().enumerate() {
            for (k, &input) in bases[j].iter().enumerate() {
                let value = surface.field.mul(model, entry, input)?;
                let output = surface.field.coordinates(model, bases[i], value)?;
                for (l, coordinate) in output.into_iter().enumerate() {
                    output_matrix[2 * i + l][2 * j + k] = as_integer(coordinate)?;
                }
            }
        }
    }
    Ok(output_matrix)
}

fn lcm(a: u64, b: u64) -> u64 {
    // Here both numbers divide the exponent, so their least common multiple fits too.
    (a / gcd(u128::from(a), u128::from(b)) as u64) * b
}

fn as_integer(value: Element) -> Result<i128> {
    if value.b == 0 && value.c == 1 {
        Ok(value.a)
    } else {
        Err(Error::plain(
            "CM theta kernel requires an integral divisor in the certified Hom lattice",
        ))
    }
}

fn identity() -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| i128::from(i == j)))
}

fn transpose(matrix: Matrix) -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| matrix[j][i]))
}

fn dot(a: [i128; 4], b: [i128; 4]) -> Result<i128> {
    let mut result = Wide::<8>::from_i128(0);
    for (x, y) in a.into_iter().zip(b) {
        result = result.add(Wide::product(x, y))?;
    }
    result.signed().ok_or_else(arithmetic)
}

fn multiply(a: Matrix, b: Matrix) -> Result<Matrix> {
    let mut result = [[0; 4]; 4];
    let columns = transpose(b);
    for (out, row) in result.iter_mut().zip(a) {
        for (entry, column) in out.iter_mut().zip(columns) {
            *entry = dot(row, column)?;
        }
    }
    Ok(result)
}

fn swap(form: &mut Matrix, basis: &mut Matrix, i: usize, j: usize) {
    form.swap(i, j);
    for row in form.iter_mut().chain(basis) {
        row.swap(i, j);
    }
}

// Add coefficient times column source to column target, with the matching row operation on the form.
fn shear(
    form: &mut Matrix,
    basis: &mut Matrix,
    target: usize,
    source: usize,
    coefficient: i128,
) -> Result<()> {
    let linear = |x, y| {
        Wide::<8>::from_i128(x)
            .add(Wide::product(y, coefficient))?
            .signed()
            .ok_or_else(arithmetic)
    };
    for row in basis {
        row[target] = linear(row[target], row[source])?;
    }
    let mut target_row = form[target];
    for (i, row) in form.iter_mut().enumerate() {
        if i != target {
            row[target] = linear(row[target], row[source])?;
            target_row[i] = row[target].checked_neg().ok_or_else(arithmetic)?;
        }
    }
    form[target] = target_row;
    Ok(())
}

fn alternating_smith(model: &mut Model, original: Matrix) -> Result<(Matrix, Matrix)> {
    let mut form = original;
    let mut basis = identity();
    for i in [0, 2] {
        let (r, s) = (i..4)
            .find_map(|r| ((r + 1)..4).find(|&s| form[r][s] != 0).map(|s| (r, s)))
            .ok_or_else(inconsistent)?;
        swap(&mut form, &mut basis, i, r);
        swap(&mut form, &mut basis, i + 1, s);
        loop {
            model.spend(16)?;
            if form[i][i + 1] < 0 {
                swap(&mut form, &mut basis, i, i + 1);
            }
            let pivot = form[i][i + 1];
            if let Some(j) = (i + 2..4).find(|&j| form[i][j] != 0) {
                let coefficient = form[i][j]
                    .div_euclid(pivot)
                    .checked_neg()
                    .ok_or_else(arithmetic)?;
                shear(&mut form, &mut basis, j, i + 1, coefficient)?;
                if form[i][j] != 0 {
                    swap(&mut form, &mut basis, j, i + 1);
                }
            } else if let Some(j) = (i + 2..4).find(|&j| form[i + 1][j] != 0) {
                let coefficient = form[i + 1][j].div_euclid(pivot);
                shear(&mut form, &mut basis, j, i, coefficient)?;
                if form[i + 1][j] != 0 {
                    swap(&mut form, &mut basis, j, i);
                }
            } else if i == 0 && form[2][3] % pivot != 0 {
                shear(&mut form, &mut basis, i, 2, 1)?;
            } else {
                break;
            }
        }
    }
    if multiply(multiply(transpose(basis), original)?, basis)? != form {
        return Err(inconsistent());
    }
    Ok((form, basis))
}

fn arithmetic() -> Error {
    Error::plain("CM theta lattice arithmetic exceeds the supported exact integer range")
}

fn inconsistent() -> Error {
    Error::plain("CM theta lattice certificate is inconsistent")
}

#[cfg(test)]
mod tests {
    use super::{alternating_smith, Matrix, Model, J};
    use crate::rv_ore::gcd;

    #[test]
    fn alternating_smith_recovers_gcd_and_pfaffian_across_integral_forms() {
        let mut model = Model::new();
        let mut seed = 1_u64;
        let mut cases = vec![[6, 0, 0, 0, 0, 10]];
        for _ in 0..256 {
            cases.push(std::array::from_fn(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                i128::from((seed >> 48) % 47) - 23
            }));
        }
        for entries @ [a, b, c, d, e, f] in cases {
            let pfaffian = (a * f - b * e + c * d).abs();
            if pfaffian == 0 {
                continue;
            }
            let first = entries
                .into_iter()
                .fold(0, |g, entry| gcd(g, entry.unsigned_abs())) as i128;
            let second = pfaffian / first;
            let form: Matrix = [[0, a, b, c], [-a, 0, d, e], [-b, -d, 0, f], [-c, -e, -f, 0]];
            let (actual, _) =
                alternating_smith(&mut model, form).unwrap_or_else(|err| panic!("{err}"));
            let mut expected = J;
            for (i, row) in expected.iter_mut().enumerate() {
                for entry in row {
                    *entry *= if i < 2 { first } else { second };
                }
            }
            assert_eq!(actual, expected);
        }
    }
}
