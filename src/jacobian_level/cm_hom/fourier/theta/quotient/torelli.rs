use std::collections::BTreeMap;
use std::rc::Rc;

use super::order::{alternating_inverse, hermite_basis, Endomorphism};
use super::{divide_matrix, identity, multiply, sum, transpose, Kernel, Matrix, Quotient, Wide};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::cm;
use crate::integer_smith::{self, Track};
use crate::model::Model;
use crate::rv_ore::gcd;
use crate::syntax::{Label, ParsedExpr};

mod shell;

pub(crate) enum Statement {
    Data(Label),
    Transcendental(Label),
    Torelli(Label),
    Frobenius(Label),
    Profiles(Label),
    Pair(Box<[ParsedExpr; 2]>),
}

pub(crate) enum Expression {
    Wedge(Label, Wedge),
    Basis(bool, Label, usize),
    Rosati(Box<ParsedExpr>),
    Projector(Label, usize),
}

pub(crate) type Wedge = [i128; 6];
const PAIRS: [[usize; 2]; 6] = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
const CUP_SIGNS: Wedge = [1, -1, 1, 1, -1, 1];

pub(super) struct Hodge {
    orientation: i128,
    algebraic: [Wedge; 4],
    algebraic_gram: [[i128; 4]; 4],
    transcendental: [Wedge; 2],
    transcendental_gram: [[i128; 2]; 2],
    determinant: i128,
    reduced: [i128; 3],
    reduced_basis: [Wedge; 2],
    proper_change: [[i128; 2]; 2],
}

pub(super) struct Splitting {
    gram: [[i128; 4]; 4],
    shell_size: usize,
    product: Option<Product>,
}

struct Product {
    projectors: [Endomorphism; 2],
    basis: Matrix,
    inverse: Matrix,
    conductors: [u64; 2],
}

impl Quotient {
    fn hodge_data(&self, model: &mut Model) -> Result<&Hodge> {
        if let Some(data) = self.hodge.get() {
            return Ok(data);
        }
        let e = self.form;
        let (_, pfaffian) = alternating_inverse(e)?;
        let orientation = pfaffian.signum();
        if orientation == 0 {
            return Err(inconsistent());
        }
        let omega = self.omega_numerator;
        let w_trace =
            self.kernel.surface.frames[2].fundamental.rem_euclid(2) * i128::from(self.conductor);
        let mut equations = integer_smith::Matrix::new(model, 6, 6)?;
        for j in 0..6 {
            let mut generator = [0; 6];
            generator[j] = 1;
            let s = alternating(generator)?;
            let residual = cm_relation(omega, w_trace, s)?;
            for (i, coefficient) in from_alternating(residual)?.into_iter().enumerate() {
                equations.write(model, i, j, coefficient)?;
            }
        }
        let algebraic: [Wedge; 4] = kernel_basis(model, equations)?;
        for &form in &algebraic {
            if cm_relation(omega, w_trace, alternating(form)?)? != [[0; 4]; 4] {
                return Err(inconsistent());
            }
        }
        if cm_relation(omega, w_trace, e)? != [[0; 4]; 4] {
            return Err(inconsistent());
        }
        let algebraic_gram = gram(orientation, algebraic)?;
        let algebraic_discriminant = determinant(algebraic_gram)?;
        let mut orthogonality = integer_smith::Matrix::new(model, 4, 6)?;
        for (i, form) in algebraic.iter().enumerate() {
            for (j, (&x, &sign)) in form.iter().rev().zip(&CUP_SIGNS).enumerate() {
                let value = x.checked_mul(sign * orientation).ok_or_else(arithmetic)?;
                orthogonality.write(model, i, j, value)?;
            }
        }
        let mut transcendental: [Wedge; 2] = kernel_basis(model, orthogonality)?;
        let complex = sum(
            multiply(super::diagonal(2), omega)?,
            super::diagonal(-w_trace),
        )?;
        let numerator = multiply(transpose(complex), alternating(transcendental[0])?)?;
        let oriented = cup(orientation, from_alternating(numerator)?, transcendental[1])?;
        // Applying I to the first argument sends Re(sigma) to -Im(sigma).
        if oriented == 0 {
            return Err(inconsistent());
        }
        if oriented > 0 {
            transcendental[1] = negate_wedge(transcendental[1])?;
        }
        for &a in &algebraic {
            for &t in &transcendental {
                if cup(orientation, a, t)? != 0 {
                    return Err(inconsistent());
                }
            }
        }
        let transcendental_gram = gram(orientation, transcendental)?;
        let [[u, v], [_, w]] = transcendental_gram;
        let discriminant = determinant(transcendental_gram)?;
        if u <= 0
            || w <= 0
            || u % 2 != 0
            || w % 2 != 0
            || discriminant <= 0
            || algebraic_discriminant.checked_neg() != Some(discriminant)
        {
            return Err(inconsistent());
        }
        let (reduced, proper_change) = cm::reduce_binary(model, [u / 2, v, w / 2])?;
        if determinant(proper_change)? != 1 {
            return Err(inconsistent());
        }
        let mut reduced_basis = [[0; 6]; 2];
        for (entry, coefficients) in reduced_basis.iter_mut().zip(proper_change) {
            *entry = combine(transcendental, coefficients)?;
        }
        let g = gram(orientation, reduced_basis)?;
        if [g[0][0] / 2, g[0][1], g[1][1] / 2] != reduced {
            return Err(inconsistent());
        }
        let data = Hodge {
            orientation,
            algebraic,
            algebraic_gram,
            transcendental,
            transcendental_gram,
            determinant: discriminant,
            reduced,
            reduced_basis,
            proper_change,
        };
        model.allocate(90)?;
        Ok(self.hodge.get_or_init(|| data))
    }

    pub(crate) fn ns_report(&self, model: &mut Model) -> Result<String> {
        let data = self.hodge_data(model)?;
        Ok(format!(
            "rvCMQuotientNS = rank=4, signature=(1,3), determinant={}, cup_orientation={}, basis={:?}, gram={:?}",
            -data.determinant, data.orientation, data.algebraic, data.algebraic_gram,
        ))
    }

    pub(crate) fn transcendental_report(&self, model: &mut Model) -> Result<String> {
        let data = self.hodge_data(model)?;
        let [a, b, c] = data.reduced;
        Ok(format!(
            "rvCMQuotientTrans = rank=2, signature=(2,0), determinant={}, basis={:?}, gram={:?}, oriented_form={:?}, unoriented_form={:?}, proper_change_columns={:?}, reduced_basis={:?}, modulus=(-({b})+sqrt(-{}))/(2*{a})",
            data.determinant, data.transcendental, data.transcendental_gram, data.reduced,
            [a, b.abs(), c], data.proper_change, data.reduced_basis, data.determinant,
        ))
    }

    pub(crate) fn wedge_basis(
        &self,
        model: &mut Model,
        transcendental: bool,
        i: usize,
    ) -> Result<Wedge> {
        let hodge = self.hodge_data(model)?;
        let basis: &[_] = if transcendental {
            &hodge.reduced_basis
        } else {
            &hodge.algebraic
        };
        basis
            .get(i)
            .copied()
            .ok_or_else(|| Error::plain("CM quotient wedge basis index is out of range"))
    }

    pub(crate) fn wedge_pair(&self, x: Wedge, y: Wedge) -> Result<i128> {
        let (_, pfaffian) = alternating_inverse(self.form)?;
        cup(pfaffian.signum(), x, y)
    }

    pub(crate) fn wedge_rosati(&self, wedge: Wedge) -> Result<Endomorphism> {
        self.need_principal()?;
        let s = alternating(wedge)?;
        let trace =
            self.kernel.surface.frames[2].fundamental.rem_euclid(2) * i128::from(self.conductor);
        if cm_relation(self.omega_numerator, trace, s)? != [[0; 4]; 4] {
            return Err(Error::plain(
                "CM quotient Rosati conversion requires an algebraic wedge form",
            ));
        }
        let (adjugate, pfaffian) = alternating_inverse(self.form)?;
        let x = divide_matrix(multiply(adjugate, s)?, pfaffian).ok_or_else(inconsistent)?;
        if multiply(self.form, x)? != s || multiply(transpose(x), self.form)? != s {
            return Err(inconsistent());
        }
        self.endomorphism(x)
    }

    fn need_principal(&self) -> Result<()> {
        if self.divisors != [1, 1] {
            return Err(Error::plain(
                "CM quotient Torelli requires a principal polarization",
            ));
        }
        Ok(())
    }

    fn splitting(&self, model: &mut Model) -> Result<&Splitting> {
        if let Some(data) = self.splitting.get() {
            return Ok(data);
        }
        self.need_principal()?;
        let data = self.hodge_data(model)?;
        let [a, b, c, d] = data.algebraic.map(|x| self.wedge_rosati(x));
        let basis = [a?.0, b?.0, c?.0, d?.0];
        let mut gram = [[0; 4]; 4];
        for (i, row) in gram.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value = trace_product(basis[i], basis[j])?;
            }
        }
        let mut shell_size = 0;
        let mut projectors = Vec::new();
        shell::each_norm_two(model, gram, |model, coefficients| {
            shell_size += 1;
            let mut x = [[0; 4]; 4];
            for (i, row) in x.iter_mut().enumerate() {
                for (j, entry) in row.iter_mut().enumerate() {
                    let column = basis.map(|a| a[i][j]);
                    *entry = super::super::dot(column, coefficients)?;
                }
            }
            if trace_product(x, x)? != 2 {
                return Err(inconsistent());
            }
            if trace_product(x, identity())? != 2 {
                return Ok(());
            }
            let map = self.endomorphism(x)?;
            if multiply(x, x)? != x
                || self.rosati_integral(map)? != map
                || matrix_rank(model, x)? != 2
            {
                return Err(inconsistent());
            }
            projectors.push(map);
            Ok(())
        })?;
        let product = match projectors.as_slice() {
            [] => None,
            &[x, y] => Some(Product::new(model, self, x, y)?),
            _ => return Err(inconsistent()),
        };
        let result = Splitting {
            gram,
            shell_size,
            product,
        };
        model.allocate(80)?;
        Ok(self.splitting.get_or_init(|| result))
    }

    pub(crate) fn torelli_report(&self, model: &mut Model) -> Result<String> {
        let splitting = self.splitting(model)?;
        let prefix = format!(
            "rvCMQuotientTorelli = trace_gram={:?}, norm_two_vectors={}",
            splitting.gram, splitting.shell_size
        );
        let Some(product) = &splitting.product else {
            return Ok(format!(
                "{prefix}, splitting=indecomposable, geometric_complex_type=genus_two_jacobian"
            ));
        };
        Ok(format!(
            "{prefix}, splitting=elliptic_product, projectors={:?}, factor_inclusion={:?}, inverse={:?}, factor_polarization_degrees=[1, 1], factor_cm_conductors={:?}",
            product.projectors.map(|p| p.0), product.basis, product.inverse, product.conductors,
        ))
    }

    pub(crate) fn factor_frobenius_report(&self, model: &mut Model) -> Result<String> {
        let product = self.splitting(model)?.product.as_ref().ok_or_else(|| {
            Error::plain("CM quotient factor Frobenius requires a polarized elliptic product")
        })?;
        let (_, degree, frobenius, multiplier) = self.certified_frobenius(model)?;
        let action = multiply(multiply(product.inverse, frobenius)?, product.basis)?;
        let mut blocks = [[[0; 2]; 2]; 2];
        for (i, row) in action.iter().enumerate() {
            for (j, &entry) in row.iter().enumerate() {
                if i / 2 == j / 2 {
                    blocks[i / 2][i % 2][j % 2] = entry;
                } else if entry != 0 {
                    return Err(inconsistent());
                }
            }
        }
        for &block in &blocks {
            if determinant(block)? != multiplier {
                return Err(inconsistent());
            }
        }
        Ok(format!(
            "rvCMQuotientFactorFrob = quotient_definition_degree={degree}, fixed_factors=[0, 1], matrices={blocks:?}, multiplier={multiplier}"
        ))
    }

    pub(crate) fn torelli_projector(&self, model: &mut Model, i: usize) -> Result<Endomorphism> {
        self.splitting(model)?
            .product
            .as_ref()
            .and_then(|p| p.projectors.get(i))
            .copied()
            .ok_or_else(|| {
                Error::plain("CM quotient Torelli splitting projector does not exist at this index")
            })
    }
}

impl Product {
    fn new(model: &mut Model, q: &Quotient, a: Endomorphism, b: Endomorphism) -> Result<Self> {
        if sum(a.0, b.0)? != identity()
            || multiply(a.0, b.0)? != [[0; 4]; 4]
            || multiply(b.0, a.0)? != [[0; 4]; 4]
        {
            return Err(inconsistent());
        }
        let mut columns = image_basis(model, a.0)?;
        let second = image_basis(model, b.0)?;
        columns[2] = second[0];
        columns[3] = second[1];
        for offset in [0, 2] {
            let form = super::pairing(q.form, columns[offset], columns[offset + 1])?;
            if form == -1 {
                columns[offset + 1] = super::negate(columns[offset + 1])?;
            } else if form != 1 {
                return Err(inconsistent());
            }
        }
        let basis = transpose(columns);
        let inverse = unimodular_inverse(basis)?;
        if multiply(inverse, basis)? != identity()
            || multiply(multiply(transpose(basis), q.form)?, basis)? != super::J
        {
            return Err(inconsistent());
        }
        let actions = multiply(multiply(inverse, q.omega_numerator)?, basis)?;
        let expected = multiply(multiply(inverse, a.0)?, basis)?;
        if expected != [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]] {
            return Err(inconsistent());
        }
        let mut conductors = [q.conductor; 2];
        for (i, row) in actions.iter().enumerate() {
            for (j, &entry) in row.iter().enumerate() {
                if i / 2 == j / 2 {
                    conductors[i / 2] =
                        gcd(u128::from(conductors[i / 2]), entry.unsigned_abs()) as u64;
                } else if entry != 0 {
                    return Err(inconsistent());
                }
            }
        }
        let conductors = conductors.map(|common| q.conductor / common);
        Ok(Self {
            projectors: [a, b],
            basis,
            inverse,
            conductors,
        })
    }
}

pub(crate) fn profiles(model: &mut Model, name: Label, kernel: Rc<Kernel>) -> Result<String> {
    let mut profile = BTreeMap::new();
    for lagrangian in kernel.lagrangians(model)? {
        let q = Quotient::new(model, name.clone(), kernel.clone(), lagrangian.clone())?;
        let determinant = q.hodge_data(model)?.determinant;
        let data = q.splitting(model)?;
        let mut conductors = data.product.as_ref().map(|p| p.conductors);
        if let Some(conductors) = &mut conductors {
            conductors.sort_unstable();
        }
        *profile
            .entry((q.conductor, determinant, conductors))
            .or_insert(0) += 1;
    }
    Ok(format!(
        "rvCMQuotientTorelliProfiles = (diagonal_conductor,transcendental_determinant,elliptic_factor_conductors): {profile:?}"
    ))
}

fn kernel_basis<const R: usize>(
    model: &mut Model,
    equations: integer_smith::Matrix,
) -> Result<[Wedge; R]> {
    let reduction = integer_smith::reduce(model, equations, Track::None, Track::Basis)?;
    let rank = reduction.factors.len();
    if rank + R != 6 {
        return Err(inconsistent());
    }
    let mut basis = std::array::from_fn(|i| {
        let coordinates = reduction.columns.basis(rank + i);
        std::array::from_fn(|j| coordinates.get(&j).copied().unwrap_or(0))
    });
    hermite_basis(model, &mut basis)?;
    Ok(basis)
}

fn cm_relation(omega: Matrix, trace: i128, s: Matrix) -> Result<Matrix> {
    sum(
        sum(multiply(transpose(omega), s)?, multiply(s, omega)?)?,
        multiply(super::diagonal(-trace), s)?,
    )
}

fn alternating(coordinates: Wedge) -> Result<Matrix> {
    let mut result = [[0; 4]; 4];
    for ([i, j], value) in PAIRS.into_iter().zip(coordinates) {
        result[i][j] = value;
        result[j][i] = value.checked_neg().ok_or_else(arithmetic)?;
    }
    Ok(result)
}

fn from_alternating(matrix: Matrix) -> Result<Wedge> {
    for (i, row) in matrix.iter().enumerate() {
        for (j, &entry) in row.iter().enumerate() {
            if entry.checked_neg() != Some(matrix[j][i]) {
                return Err(inconsistent());
            }
        }
    }
    Ok(PAIRS.map(|[i, j]| matrix[i][j]))
}

pub(crate) fn negate_wedge(mut x: Wedge) -> Result<Wedge> {
    for entry in &mut x {
        *entry = entry.checked_neg().ok_or_else(arithmetic)?;
    }
    Ok(x)
}

fn cup(orientation: i128, x: Wedge, y: Wedge) -> Result<i128> {
    let mut answer = Wide::<16>::from_i128(0);
    for ((a, b), sign) in x.into_iter().zip(y.into_iter().rev()).zip(CUP_SIGNS) {
        answer = answer.add(Wide::product(a, b).multiply(sign * orientation)?)?;
    }
    answer.signed().ok_or_else(arithmetic)
}

fn gram<const N: usize>(orientation: i128, rows: [Wedge; N]) -> Result<[[i128; N]; N]> {
    let mut result = [[0; N]; N];
    for (out, x) in result.iter_mut().zip(rows) {
        for (value, y) in out.iter_mut().zip(rows) {
            *value = cup(orientation, x, y)?;
        }
    }
    Ok(result)
}

fn combine(basis: [Wedge; 2], [a, b]: [i128; 2]) -> Result<Wedge> {
    let mut vector = [0; 6];
    for (i, entry) in vector.iter_mut().enumerate() {
        *entry = Wide::<8>::product(a, basis[0][i])
            .add(Wide::product(b, basis[1][i]))?
            .signed()
            .ok_or_else(arithmetic)?;
    }
    Ok(vector)
}

fn trace_product(x: Matrix, y: Matrix) -> Result<i128> {
    let mut trace = Wide::<16>::from_i128(0);
    for (row, column) in x.into_iter().zip(transpose(y)) {
        for (a, b) in row.into_iter().zip(column) {
            trace = trace.add(Wide::product(a, b))?;
        }
    }
    trace.signed().ok_or_else(arithmetic)
}

fn smith(model: &mut Model, matrix: Matrix, rows: Track) -> Result<integer_smith::Smith> {
    let mut input = integer_smith::Matrix::new(model, 4, 4)?;
    for (i, row) in matrix.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            input.write(model, i, j, value)?;
        }
    }
    integer_smith::reduce(model, input, rows, Track::None)
}

fn matrix_rank(model: &mut Model, matrix: Matrix) -> Result<usize> {
    Ok(smith(model, matrix, Track::None)?.factors.len())
}

fn image_basis(model: &mut Model, matrix: Matrix) -> Result<Matrix> {
    let reduction = smith(model, matrix, Track::Basis)?;
    if reduction.factors.as_slice() != [1, 1] {
        return Err(inconsistent());
    }
    let mut actual = [[0; 4]; 4];
    for (i, column) in actual.iter_mut().enumerate().take(2) {
        let basis = reduction.rows.basis(i);
        *column = std::array::from_fn(|j| basis.get(&j).copied().unwrap_or(0));
        if super::apply_matrix(matrix, *column)? != *column {
            return Err(inconsistent());
        }
    }
    Ok(actual)
}

fn determinant<const N: usize>(matrix: [[i128; N]; N]) -> Result<i128> {
    minor(&matrix, (1 << N) - 1, (1 << N) - 1)?
        .signed()
        .ok_or_else(arithmetic)
}

fn minor<const N: usize>(matrix: &[[i128; N]; N], rows: usize, columns: usize) -> Result<Wide<16>> {
    if rows == 0 {
        return Ok(Wide::from_i128(1));
    }
    let i = rows.trailing_zeros() as usize;
    let mut result = Wide::<16>::from_i128(0);
    let mut offset = 0_u32;
    let mut remaining = columns;
    while remaining != 0 {
        let j = remaining.trailing_zeros() as usize;
        let contribution =
            minor(matrix, rows ^ (1 << i), columns ^ (1 << j))?.multiply(matrix[i][j])?;
        result = result.add(if offset.is_multiple_of(2) {
            contribution
        } else {
            contribution.negate()
        })?;
        offset += 1;
        remaining ^= 1 << j;
    }
    Ok(result)
}

fn unimodular_inverse(matrix: Matrix) -> Result<Matrix> {
    let determinant = determinant(matrix)?;
    if determinant.unsigned_abs() != 1 {
        return Err(inconsistent());
    }
    let mut inverse = [[0; 4]; 4];
    for (i, row) in inverse.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            let sign = if (i + j).is_multiple_of(2) {
                determinant
            } else {
                -determinant
            };
            *value = minor(&matrix, 15 ^ (1 << j), 15 ^ (1 << i))?
                .multiply(sign)?
                .signed()
                .ok_or_else(arithmetic)?;
        }
    }
    Ok(inverse)
}

fn arithmetic() -> Error {
    Error::plain("CM quotient Hodge lattice exceeds the supported exact integer range")
}

fn inconsistent() -> Error {
    Error::plain("CM quotient Hodge or Torelli lattice certificate is inconsistent")
}

pub(crate) fn add_wedge(mut a: Wedge, b: Wedge) -> Result<Wedge> {
    for (x, y) in a.iter_mut().zip(b) {
        *x = x.checked_add(y).ok_or_else(arithmetic)?;
    }
    Ok(a)
}
