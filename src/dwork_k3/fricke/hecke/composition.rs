use crate::ec_poly::{self as poly, Field, Polynomial as FieldPolynomial, Scalar};
use crate::elliptic_modular::{Integer, Polynomial};
use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::model::Model;
use crate::rv_unram::Ring;

use super::{factors, prime_level, require_field};

pub(crate) fn pair(
    model: &mut Model,
    levels: [usize; 2],
    ring: &Ring,
    x: &[u64],
    y: &[u64],
) -> Result<Scalar> {
    let [a, b] = levels;
    prime_level(a)?;
    prime_level(b)?;
    let first = model.dwork_hecke_equation(a)?;
    let second = model.dwork_hecke_equation(b)?;
    let f = specialize(model, &first, ring, x, a + 1)?;
    let g = specialize(model, &second, ring, y, b + 1)?;
    resultant(model, ring, &f, &g)
}

pub(crate) fn local(
    model: &mut Model,
    levels: [usize; 2],
    ring: &Ring,
    x: &[u64],
) -> Result<FieldPolynomial> {
    let [a, b] = levels;
    prime_level(a)?;
    prime_level(b)?;
    require_field(ring)?;
    let full = local_resultant(model, levels, ring, x)?;
    let field = Field { ring };
    if full.len() != (a + 1) * (b + 1) + 1 || full.last() != Some(&field.one()) {
        return Err(certificate());
    }
    if a != b {
        let reverse = local_resultant(model, [b, a], ring, x)?;
        if full != reverse {
            return Err(certificate());
        }
        return Ok(full);
    }

    let linear = vec![ring.negate(model, x)?, field.one()];
    let mut diagonal = vec![field.one()];
    for _ in 0..=a {
        diagonal = poly::multiply(model, field, &diagonal, &linear)?;
    }
    let primitive = poly::exact_divide(model, field, &full, &diagonal)?;
    if a == 3 {
        let global = model.dwork_hecke_equation(9)?;
        if primitive != specialize(model, &global, ring, x, 12)? {
            return Err(certificate());
        }
    }
    Ok(primitive)
}

pub(crate) fn report(
    model: &mut Model,
    levels: [usize; 2],
    ring: &Ring,
    x: &[u64],
) -> Result<String> {
    let [a, b] = levels;
    let polynomial = local(model, levels, ring, x)?;
    let data = factors(model, ring, &polynomial)?;
    let degree = polynomial.len() - 1;
    let resultant_degree = (a + 1) * (b + 1);
    let characteristic = if [a, b].contains(&(ring.prime() as usize)) {
        ", level characteristic with schematic multiplicities"
    } else {
        ""
    };
    let detail = if a == b {
        let cyclic_diagonal = data
            .roots
            .iter()
            .find(|(root, _)| root == x)
            .map_or(0, |(_, multiplicity)| *multiplicity);
        format!(
            "backtracking divisor degree={}, cyclic diagonal multiplicity={cyclic_diagonal}",
            a + 1
        )
    } else {
        "coprime composition agrees in both orders".into()
    };
    Ok(format!(
        "rvDworkHeckeTower[{a}, {b}] = H={x:?}, resultant degree={resultant_degree}, cyclic degree={degree}, {detail}{characteristic}, {}",
        data.describe()
    ))
}

fn local_resultant(
    model: &mut Model,
    [a, b]: [usize; 2],
    ring: &Ring,
    x: &[u64],
) -> Result<FieldPolynomial> {
    let first = model.dwork_hecke_equation(a)?;
    let second = model.dwork_hecke_equation(b)?;
    let f = specialize(model, &first, ring, x, a + 1)?
        .into_iter()
        .map(|c| poly::trim(vec![c]))
        .collect::<Vec<_>>();
    let mut g = Vec::new();
    for i in 0..=b + 1 {
        let row = (0..=b + 1)
            .map(|j| second.field_coefficient(model, ring, i, j))
            .collect::<Result<_>>()?;
        g.push(poly::trim(row));
    }
    resultant(model, &FieldPolynomials(Field { ring }), &f, &g)
}

pub(super) fn certify_nine(model: &mut Model, nine: &Polynomial) -> Result<()> {
    let three = model.dwork_hecke_equation(3)?;
    // Each side has degree <=16 in each variable. Both are symmetric:
    // swapping two monic quartics in the resultant has sign (-1)^16.
    let specializations = (0..=16)
        .map(|x| {
            (0..=4)
                .map(|j| integral_specialization(&three, x, j, 4))
                .map(|value| value.map(Integral))
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    for (x, first) in specializations.iter().enumerate() {
        let row = (0..=12)
            .map(|j| integral_specialization(nine, x as i128, j, 12))
            .collect::<Result<Vec<_>>>()?;
        for (y, second) in specializations.iter().enumerate().skip(x) {
            let actual = resultant(model, &Integers, first, second)?.0;
            let expected =
                integer_value(&row, y as i128)?.multiply((x as i128 - y as i128).pow(4))?;
            if actual.magnitude != expected.magnitude
                || (!actual.is_zero() && actual.negative != expected.negative)
            {
                return Err(certificate());
            }
        }
    }
    Ok(())
}

fn integral_specialization(p: &Polynomial, x: i128, j: usize, degree: usize) -> Result<Integer> {
    (0..=degree)
        .rev()
        .try_fold(Integer::from_i128(0), |sum, i| {
            sum.multiply(x)?.add(p.wide_coefficient::<40>(i, j)?)
        })
}

fn integer_value(coefficients: &[Integer], x: i128) -> Result<Integer> {
    coefficients
        .iter()
        .rev()
        .try_fold(Integer::from_i128(0), |sum, &c| sum.multiply(x)?.add(c))
}

pub(super) fn specialize(
    model: &mut Model,
    p: &Polynomial,
    ring: &Ring,
    x: &[u64],
    degree: usize,
) -> Result<FieldPolynomial> {
    (0..=degree)
        .map(|j| p.specialization(model, ring, x, j))
        .collect()
}

// For monic f, Res(f,g) is the determinant of multiplication by g
// on R[Z]/f. All arithmetic, including the determinant, is division-free.
fn resultant<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    f: &[R::Element],
    g: &[R::Element],
) -> Result<R::Element> {
    let Some((leading, lower)) = f.split_last() else {
        return Err(certificate());
    };
    if *leading != ring.one() || lower.is_empty() || !lower.len().is_multiple_of(2) {
        return Err(certificate());
    }
    let n = lower.len();
    let mut column = vec![ring.zero(); n];
    for coefficient in g.iter().rev() {
        multiply_by_z(model, ring, lower, &mut column)?;
        column[0] = ring.add(model, &column[0], coefficient)?;
    }
    let mut matrix = vec![ring.zero(); n * n];
    for j in 0..n {
        for (i, coefficient) in column.iter().enumerate() {
            matrix[i * n + j] = coefficient.clone();
        }
        if j + 1 < n {
            multiply_by_z(model, ring, lower, &mut column)?;
        }
    }
    let polynomial = hw_algebra::characteristic(model, ring, &Matrix::from_entries(n, matrix)?)?;
    Ok(polynomial.into_iter().last().unwrap())
}

fn multiply_by_z<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    lower: &[R::Element],
    column: &mut [R::Element],
) -> Result<()> {
    let leading = column.last().unwrap().clone();
    for i in (0..column.len()).rev() {
        let term = ring.multiply(model, &leading, &lower[i])?;
        let shifted = if i == 0 {
            ring.zero()
        } else {
            column[i - 1].clone()
        };
        column[i] = ring.subtract(model, &shifted, &term)?;
    }
    Ok(())
}

fn certificate() -> Error {
    Error::plain("Dwork Fricke Hecke resultant failed its exact composition certificate")
}

#[derive(Clone, Copy)]
struct Integral(Integer);

impl PartialEq for Integral {
    fn eq(&self, other: &Self) -> bool {
        self.0.magnitude == other.0.magnitude
            && (self.0.is_zero() || self.0.negative == other.0.negative)
    }
}

impl Eq for Integral {}

struct Integers;

impl CoefficientRing for Integers {
    type Element = Integral;

    fn cells(&self) -> usize {
        41
    }
    fn zero(&self) -> Integral {
        Integral(Integer::from_i128(0))
    }
    fn one(&self) -> Integral {
        Integral(Integer::from_i128(1))
    }
    fn is_zero(&self, value: &Integral) -> bool {
        value.0.is_zero()
    }
    fn is_unit(&self, value: &Integral) -> bool {
        matches!(value.0.signed(), Some(1 | -1))
    }

    fn add(&self, model: &mut Model, a: &Integral, b: &Integral) -> Result<Integral> {
        model.spend(2)?;
        Ok(Integral(a.0.add(b.0)?))
    }

    fn subtract(&self, model: &mut Model, a: &Integral, b: &Integral) -> Result<Integral> {
        self.add(model, a, &Integral(b.0.negate()))
    }

    fn multiply(&self, model: &mut Model, a: &Integral, b: &Integral) -> Result<Integral> {
        model.spend(8)?;
        Ok(Integral(a.0.multiply_wide(b.0)?))
    }

    fn inverse_unit(&self, _: &mut Model, value: &Integral) -> Result<Integral> {
        self.is_unit(value)
            .then_some(*value)
            .ok_or_else(certificate)
    }
}

struct FieldPolynomials<'a>(Field<'a>);

impl CoefficientRing for FieldPolynomials<'_> {
    type Element = FieldPolynomial;

    fn cells(&self) -> usize {
        self.0.cells()
    }
    fn zero(&self) -> FieldPolynomial {
        Vec::new()
    }
    fn one(&self) -> FieldPolynomial {
        vec![self.0.one()]
    }
    fn is_zero(&self, value: &FieldPolynomial) -> bool {
        value.is_empty()
    }
    fn is_unit(&self, value: &FieldPolynomial) -> bool {
        value.len() == 1
    }

    fn add(
        &self,
        model: &mut Model,
        a: &FieldPolynomial,
        b: &FieldPolynomial,
    ) -> Result<FieldPolynomial> {
        poly::sum(model, self.0, a, b, false)
    }

    fn subtract(
        &self,
        model: &mut Model,
        a: &FieldPolynomial,
        b: &FieldPolynomial,
    ) -> Result<FieldPolynomial> {
        poly::sum(model, self.0, a, b, true)
    }

    fn multiply(
        &self,
        model: &mut Model,
        a: &FieldPolynomial,
        b: &FieldPolynomial,
    ) -> Result<FieldPolynomial> {
        poly::multiply(model, self.0, a, b)
    }

    fn inverse_unit(&self, model: &mut Model, value: &FieldPolynomial) -> Result<FieldPolynomial> {
        if let [constant] = value.as_slice() {
            return Ok(vec![self.0.inverse_unit(model, constant)?]);
        }
        Err(certificate())
    }
}
