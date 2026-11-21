use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::binary_field::BinaryField;
use crate::binary_linear::{check_storage, Entry};
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::rv_hom::Grade;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Term {
    pub(crate) coordinate: usize,
    pub(crate) power: i64,
    pub(crate) coefficient: usize,
}

pub(crate) type Support = Rc<[Term]>;
pub(crate) type Polynomial = Vec<(i64, usize)>;
pub(crate) type PolynomialMatrix = Vec<Vec<Polynomial>>;

pub(crate) fn grade(model: &Model, ty: TypeId, name: &str) -> Result<Grade> {
    match model.kind(ty) {
        Kind::RVLaurent(grade) => Ok(grade.clone()),
        _ => Err(Error::plain(format!("{name} requires an RVLaurent map"))),
    }
}

pub(crate) fn dimensions(model: &mut Model, ty: TypeId) -> Result<(usize, usize)> {
    let grade = grade(model, ty, "Laurent matrix")?;
    let hom = model.rv_hom_type(grade)?;
    Ok(model.rv_hom(hom)?.dimensions())
}

pub(crate) fn persist(model: &mut Model, terms: Vec<Term>) -> Result<Value> {
    check_storage(3 * terms.len())?;
    model.allocate(3 * terms.len())?;
    Ok(Value::RVLaurent(terms.into()))
}

pub(crate) fn zero() -> Value {
    Value::RVLaurent([].into())
}

pub(crate) fn monomial(model: &mut Model, terms: &[Entry], power: i64) -> Result<Value> {
    model.spend(terms.len())?;
    check_storage(3 * terms.len())?;
    persist(
        model,
        terms
            .iter()
            .map(|&(coordinate, coefficient)| Term {
                coordinate,
                power,
                coefficient,
            })
            .collect(),
    )
}

pub(crate) fn compare(model: &mut Model, lhs: &[Term], rhs: &[Term]) -> Result<Ordering> {
    model.spend(1 + lhs.len().min(rhs.len()))?;
    Ok(lhs.iter().cmp(rhs.iter()))
}

fn key(term: &Term) -> (usize, i64) {
    (term.coordinate, term.power)
}

pub(crate) fn add(model: &mut Model, lhs: &[Term], rhs: &[Term]) -> Result<Value> {
    let sum = add_support(model, lhs, rhs)?;
    persist(model, sum)
}

pub(crate) fn add_support(
    model: &mut Model,
    mut lhs: &[Term],
    mut rhs: &[Term],
) -> Result<Vec<Term>> {
    model.spend(lhs.len() + rhs.len())?;
    let mut result = Vec::new();
    while let (Some((a, rest_left)), Some((b, rest_right))) = (lhs.split_first(), rhs.split_first())
    {
        match key(a).cmp(&key(b)) {
            Ordering::Less => {
                result.push(*a);
                lhs = rest_left;
            }
            Ordering::Equal => {
                let coefficient = a.coefficient ^ b.coefficient;
                if coefficient != 0 {
                    result.push(Term { coefficient, ..*a });
                }
                lhs = rest_left;
                rhs = rest_right;
            }
            Ordering::Greater => {
                result.push(*b);
                rhs = rest_right;
            }
        }
    }
    check_storage(3 * (result.len() + lhs.len() + rhs.len()))?;
    result.extend_from_slice(lhs);
    result.extend_from_slice(rhs);
    Ok(result)
}

pub(crate) fn scale(
    model: &mut Model,
    field: BinaryField,
    terms: &[Term],
    scalar: usize,
) -> Result<Value> {
    if scalar == 0 {
        return Ok(zero());
    }
    model.spend(field.degree() * terms.len())?;
    persist(
        model,
        terms
            .iter()
            .map(|&term| Term {
                coefficient: field.multiply(scalar, term.coefficient),
                ..term
            })
            .collect(),
    )
}

pub(crate) fn frobenius(
    model: &mut Model,
    field: BinaryField,
    terms: &[Term],
    iterations: usize,
) -> Result<Vec<Term>> {
    model.spend(terms.len())?;
    terms
        .iter()
        .map(|&term| {
            Ok(Term {
                coefficient: field.frobenius(model, term.coefficient, iterations)?,
                ..term
            })
        })
        .collect()
}

pub(crate) fn frobenius_type(
    model: &mut Model,
    mut grade: Grade,
    iterations: usize,
) -> Result<TypeId> {
    model.spend(grade.point.len())?;
    let point = grade
        .point
        .iter()
        .map(|&x| grade.field.frobenius(model, x, iterations))
        .collect::<Result<Vec<_>>>()?;
    model.allocate(1 + point.len())?;
    grade.point = point.into();
    model.rv_laurent_type(grade)
}

pub(crate) fn composition_type(model: &mut Model, outer: TypeId, inner: TypeId) -> Result<TypeId> {
    let large = grade(model, outer, "ecomp")?;
    let small = grade(model, inner, "ecomp")?;
    if small.target != large.source {
        return Err(Error::plain(
            "ecomp requires composable Laurent Margolis endpoints",
        ));
    }
    if small.field != large.field || small.selected != large.selected || small.point != large.point
    {
        return Err(Error::plain(
            "ecomp requires the same affine Margolis field point",
        ));
    }
    model.rv_laurent_type(Grade {
        source: small.source,
        ..large
    })
}

pub(crate) fn compose(
    model: &mut Model,
    outer_ty: TypeId,
    outer: &[Term],
    inner_ty: TypeId,
    inner: &[Term],
) -> Result<Value> {
    if outer.is_empty() || inner.is_empty() {
        return Ok(zero());
    }
    let field = grade(model, outer_ty, "ecomp")?.field;
    let (_, middle) = dimensions(model, inner_ty)?;
    let (_, output) = dimensions(model, outer_ty)?;
    let value = product(model, field, output, middle, outer, inner)?;
    persist(model, value)
}

pub(crate) fn add_powers(left: i64, right: i64) -> Result<i64> {
    left.checked_add(right)
        .ok_or_else(|| Error::plain("Laurent exponent is out of the signed 64-bit range"))
}

pub(crate) fn accumulate<K: Ord>(
    model: &mut Model,
    terms: &mut BTreeMap<K, usize>,
    key: K,
    coefficient: usize,
) -> Result<()> {
    model.spend(1 + terms.len().max(1).ilog2() as usize)?;
    if let Some(old) = terms.get_mut(&key) {
        *old ^= coefficient;
        if *old == 0 {
            terms.remove(&key);
        }
    } else if coefficient != 0 {
        check_storage(3 * (1 + terms.len()))?;
        terms.insert(key, coefficient);
    }
    Ok(())
}

pub(crate) fn product(
    model: &mut Model,
    field: BinaryField,
    output: usize,
    middle: usize,
    outer: &[Term],
    inner: &[Term],
) -> Result<Vec<Term>> {
    let mut terms = BTreeMap::new();
    if outer.is_empty() || inner.is_empty() {
        return Ok(Vec::new());
    }
    let search = 1 + outer.len().ilog2() as usize;
    for &a in inner {
        let (source, matched) = (a.coordinate / middle, a.coordinate % middle);
        model.spend(search)?;
        let start = outer.partition_point(|entry| entry.coordinate / output < matched);
        for &b in &outer[start..] {
            model.spend(1)?;
            if b.coordinate / output != matched {
                break;
            }
            model.spend(field.degree())?;
            let coordinate = source
                .checked_mul(output)
                .and_then(|offset| offset.checked_add(b.coordinate % output))
                .ok_or_else(|| Error::plain("Laurent matrix coordinate is out of range"))?;
            let power = add_powers(a.power, b.power)?;
            accumulate(
                model,
                &mut terms,
                (coordinate, power),
                field.multiply(a.coefficient, b.coefficient),
            )?;
        }
    }
    Ok(terms
        .into_iter()
        .map(|((coordinate, power), coefficient)| Term {
            coordinate,
            power,
            coefficient,
        })
        .collect())
}

pub(crate) fn add_polynomial_product(
    model: &mut Model,
    field: BinaryField,
    sum: &mut BTreeMap<i64, usize>,
    lhs: &Polynomial,
    rhs: &Polynomial,
) -> Result<()> {
    for &(a, x) in lhs {
        model.spend(1)?;
        for &(b, y) in rhs {
            let power = add_powers(a, b)?;
            model.spend(field.degree())?;
            accumulate(model, sum, power, field.multiply(x, y))?;
        }
    }
    Ok(())
}

pub(crate) fn polynomial_matrix(
    model: &mut Model,
    size: usize,
    terms: &[Term],
) -> Result<(PolynomialMatrix, usize)> {
    let cells = size
        .checked_mul(size)
        .and_then(|cells| cells.checked_add(2 * terms.len()))
        .ok_or_else(|| Error::plain("Laurent polynomial matrix exceeds the storage limit"))?;
    check_storage(cells)?;
    model.spend(cells)?;
    let mut matrix = vec![vec![Polynomial::new(); size]; size];
    for &Term {
        coordinate,
        power,
        coefficient,
    } in terms
    {
        matrix[coordinate % size][coordinate / size].push((power, coefficient));
    }

    Ok((matrix, cells))
}
