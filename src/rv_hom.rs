use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::binary_field::BinaryField;
use crate::binary_linear::{add_multiple, check_storage, Basis, Entry, Quotient, Vector};
use crate::comodule::Generator;
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::q_koszul::{Indices, Term};
use crate::rank_variety::Endpoint;

pub(crate) type Support = Rc<[Entry]>;

#[derive(Clone, PartialEq, Eq, Hash)]
pub(crate) struct Grade {
    pub(crate) source: TypeId,
    pub(crate) target: TypeId,
    pub(crate) selected: Indices,
    pub(crate) field: BinaryField,
    pub(crate) point: Rc<[usize]>,
}

impl Grade {
    pub(crate) fn new(
        model: &mut Model,
        source: TypeId,
        target: TypeId,
        selected: Indices,
        modulus: usize,
        point: Vec<usize>,
    ) -> Result<Self> {
        let field = BinaryField::new(model, modulus)?;
        field.check_point(&point, selected.len(), "Margolis point")?;
        model.check_q_endpoints(source, target, &selected)?;
        model.allocate(1 + point.len())?;
        Ok(Self {
            source,
            target,
            selected,
            field,
            point: point.into(),
        })
    }
}

struct Margolis {
    offsets: BTreeMap<usize, (usize, usize)>,
    differential: Vec<Vector>,
    quotient: Quotient,
}

impl Margolis {
    fn build(model: &mut Model, grade: &Grade, base: TypeId) -> Result<Self> {
        let (active, point): (Vec<_>, Vec<_>) = grade
            .selected
            .iter()
            .copied()
            .zip(grade.point.iter().copied())
            .filter(|&(_, value)| value != 0)
            .unzip();
        let endpoint = Endpoint::build(model, base, &active)?;
        let size = endpoint.columns.len();
        let field = grade.field;
        let mut boundaries = Basis::default();
        let mut differential = Vec::with_capacity(size);
        let mut cells = size;
        for column in 0..size {
            let image = endpoint.weighted_column(model, column, &point)?;
            cells += 2 * image.len();
            check_storage(cells)?;
            boundaries.insert(model, field, image.clone())?;
            differential.push(image);
        }
        let quotient = Quotient::build(model, field, boundaries, size, |model, column| {
            let image = &differential[column];
            model.spend(image.len())?;
            Ok(image.clone())
        })?;
        let offsets = endpoint.into_groups();
        model.allocate(cells + 3 * offsets.len())?;
        Ok(Self {
            offsets,
            differential,
            quotient,
        })
    }

    fn rank(&self) -> usize {
        self.quotient.rank()
    }

    fn index(&self, generator: Generator) -> Result<usize> {
        self.offsets
            .get(&generator.degree)
            .and_then(|&(start, rank)| (generator.index < rank).then(|| start + generator.index))
            .ok_or_else(|| Error::plain("coefficient generator is outside the Margolis endpoint"))
    }
}

pub(crate) struct Hom {
    source: Rc<Margolis>,
    target: Rc<Margolis>,
    rank: usize,
}

impl Hom {
    pub(crate) fn build(model: &mut Model, grade: &Grade) -> Result<Self> {
        let source = Rc::new(Margolis::build(model, grade, grade.source)?);
        let target = if grade.target == grade.source {
            source.clone()
        } else {
            Rc::new(Margolis::build(model, grade, grade.target)?)
        };
        let rank = source
            .rank()
            .checked_mul(target.rank())
            .ok_or_else(|| Error::plain("Margolis Hom dimension is out of range"))?;
        Ok(Self {
            source,
            target,
            rank,
        })
    }

    pub(crate) fn rank(&self) -> usize {
        self.rank
    }

    pub(crate) fn dimensions(&self) -> (usize, usize) {
        (self.source.rank(), self.target.rank())
    }
}

pub(crate) fn zero() -> Value {
    Value::RVHom([].into())
}

pub(crate) fn persist(model: &mut Model, coordinates: Vector) -> Result<Value> {
    check_storage(2 * coordinates.len())?;
    model.allocate(2 * coordinates.len())?;
    Ok(Value::RVHom(coordinates.into()))
}

pub(crate) fn identity(model: &mut Model, ty: TypeId) -> Result<Value> {
    let hom = model.rv_hom(ty)?;
    let size = hom.source.rank();
    model.spend(size)?;
    check_storage(2 * size)?;
    persist(
        model,
        (0..size).map(|index| (index * size + index, 1)).collect(),
    )
}

pub(crate) fn add(
    model: &mut Model,
    field: BinaryField,
    left: &[Entry],
    right: &[Entry],
) -> Result<Value> {
    let value = add_multiple(model, field, left, right, 1)?;
    persist(model, value)
}

pub(crate) fn scale(
    model: &mut Model,
    field: BinaryField,
    value: &[Entry],
    scalar: usize,
) -> Result<Value> {
    if scalar == 0 || value.is_empty() {
        return Ok(zero());
    }
    model.spend(field.degree() * value.len())?;
    persist(
        model,
        value
            .iter()
            .map(|&(index, coefficient)| (index, field.multiply(scalar, coefficient)))
            .collect(),
    )
}

pub(crate) fn compare(model: &mut Model, left: &[Entry], right: &[Entry]) -> Result<Ordering> {
    model.spend(1 + left.len().min(right.len()))?;
    Ok(left.iter().rev().cmp(right.iter().rev()))
}

pub(crate) fn elements(
    model: &mut Model,
    ty: TypeId,
    field: BinaryField,
    limit: usize,
) -> Result<Vec<Value>> {
    let rank = model.rv_hom(ty)?.rank();
    let q = field.order();
    let mut size = 1;
    for _ in 0..rank {
        if size > limit / q {
            return Err(Error::plain(format!(
                "type exceeds the enumeration limit of {limit} elements"
            )));
        }
        size *= q;
    }
    model.allocate(size)?;
    model.spend(size * (1 + rank))?;
    let mut values = Vec::with_capacity(size);
    for mut index in 0..size {
        let mut coefficients = Vec::new();
        for generator in 0..rank {
            let digit = index % q;
            if digit != 0 {
                coefficients.push((generator, digit));
            }
            index /= q;
        }
        values.push(persist(model, coefficients)?);
    }
    Ok(values)
}

pub(crate) fn composition_type(model: &mut Model, outer: TypeId, inner: TypeId) -> Result<TypeId> {
    let (Kind::RVHom(g), Kind::RVHom(f)) = (model.kind(outer), model.kind(inner)) else {
        return Err(Error::plain("ecomp requires two RVHom maps"));
    };
    if f.target != g.source {
        return Err(Error::plain("ecomp requires composable Margolis endpoints"));
    }
    if f.selected != g.selected || f.field != g.field || f.point != g.point {
        return Err(Error::plain(
            "ecomp requires the same affine Margolis field point",
        ));
    }
    model.rv_hom_type(Grade {
        source: f.source,
        ..g.clone()
    })
}

pub(crate) fn compose(
    model: &mut Model,
    outer_ty: TypeId,
    outer: &[Entry],
    inner_ty: TypeId,
    inner: &[Entry],
) -> Result<Value> {
    if outer.is_empty() || inner.is_empty() {
        return Ok(zero());
    }
    let Kind::RVHom(grade) = model.kind(outer_ty) else {
        unreachable!()
    };
    let field = grade.field;
    let middle = model.rv_hom(inner_ty)?.target.rank();
    let output = model.rv_hom(outer_ty)?.target.rank();
    let mut matrix = BTreeMap::new();
    let search = 1 + outer.len().max(1).ilog2() as usize;
    for &(index, coefficient) in inner {
        let (source, matched) = (index / middle, index % middle);
        model.spend(search)?;
        let start = outer.partition_point(|&(entry, _)| entry / output < matched);
        for &(entry, scalar) in &outer[start..] {
            model.spend(1)?;
            if entry / output != matched {
                break;
            }
            model.spend(field.degree())?;
            let target = source
                .checked_mul(output)
                .and_then(|offset| offset.checked_add(entry % output))
                .ok_or_else(|| Error::plain("Margolis Hom coordinate is out of range"))?;
            accumulate(
                model,
                &mut matrix,
                target,
                field.multiply(coefficient, scalar),
            )?;
        }
    }
    persist(model, matrix.into_iter().collect())
}

pub(crate) fn specialize_type(
    model: &mut Model,
    input: TypeId,
    modulus: usize,
    point: Vec<usize>,
) -> Result<TypeId> {
    let (source, target, selected, localized) = match model.kind(input) {
        Kind::QKoszul(g) | Kind::QExt(g) => (g.source, g.target, g.selected.clone(), false),
        Kind::VKoszul(g) | Kind::VExt(g) => (g.source, g.target, g.selected.clone(), true),
        _ => {
            return Err(Error::plain(
                "rvSpec requires a polynomial or localized Koszul cochain or class",
            ))
        }
    };
    let grade = Grade::new(model, source, target, selected, modulus, point)?;
    if localized && grade.point.last() == Some(&0) {
        return Err(Error::plain(
            "rvSpec requires a nonzero top coefficient for localized Koszul values",
        ));
    }
    model.rv_hom_type(grade)
}

pub(crate) fn specialize(
    model: &mut Model,
    input: TypeId,
    value: Value,
    output: TypeId,
) -> Result<Value> {
    let Kind::RVHom(grade) = model.kind(output) else {
        unreachable!()
    };
    let grade = grade.clone();
    let localized = matches!(model.kind(input), Kind::VKoszul(..) | Kind::VExt(..));
    let field = grade.field;
    let (terms, raw) = match value {
        Value::QKoszul(terms) => (terms, true),
        Value::Class(class) => {
            if class.is_empty() {
                return Ok(zero());
            }
            let terms = model.q_ext(input)?.representative(model, &class)?;
            (terms.into(), false)
        }
        _ => unreachable!(),
    };
    if terms.is_empty() {
        return Ok(zero());
    }
    let hom = model.rv_hom(output)?;
    let evaluator = Evaluator::new(model, &grade, localized)?;
    let evaluated = matrix(model, &hom, &evaluator, &terms)?;
    if evaluated.is_empty() {
        return Ok(zero());
    }
    let mut columns = vec![Vec::new(); hom.source.differential.len()];
    for ((source, target), coefficient) in evaluated {
        columns[source].push((target, coefficient));
    }
    if raw {
        verify_chain_map(model, field, &hom, &columns)?;
    }
    if hom.rank() == 0 {
        return Ok(zero());
    }
    let width = hom.target.rank();
    let mut result = Vec::new();
    for generator in 0..hom.source.rank() {
        model.spend(1)?;
        let mut image = Vec::new();
        for &(column, scalar) in hom.source.quotient.representative(generator) {
            image = add_multiple(model, field, &image, &columns[column], scalar)?;
        }
        let image = hom.target.quotient.project(model, field, image)?;
        check_storage(2 * (result.len() + image.len()))?;
        result.extend(
            image
                .into_iter()
                .map(|(index, coefficient)| (generator * width + index, coefficient)),
        );
    }
    persist(model, result)
}

fn verify_chain_map(
    model: &mut Model,
    field: BinaryField,
    hom: &Hom,
    columns: &[Vector],
) -> Result<()> {
    for (column, source_edge) in columns.iter().zip(&hom.source.differential) {
        model.spend(1)?;
        let mut image = Vec::new();
        for &(index, scalar) in column {
            image = add_multiple(
                model,
                field,
                &image,
                &hom.target.differential[index],
                scalar,
            )?;
        }
        for &(index, scalar) in source_edge {
            image = add_multiple(model, field, &image, &columns[index], scalar)?;
        }
        if !image.is_empty() {
            return Err(Error::plain(
                "rvSpec requires a chain map at the given Margolis point",
            ));
        }
    }
    Ok(())
}

struct Evaluator<'a> {
    grade: &'a Grade,
    inverse: Option<usize>,
}

impl<'a> Evaluator<'a> {
    fn new(model: &mut Model, grade: &'a Grade, localized: bool) -> Result<Self> {
        let inverse = if localized {
            Some(grade.field.inverse(model, *grade.point.last().unwrap())?)
        } else {
            None
        };
        Ok(Self { grade, inverse })
    }

    fn coefficient(&self, model: &mut Model, powers: &[i128]) -> Result<usize> {
        let field = self.grade.field;
        let mut scalar = 1;
        for (&power, &value) in powers.iter().zip(self.grade.point.iter()) {
            if power == 0 {
                continue;
            }
            let base = if power < 0 {
                self.inverse
                    .ok_or_else(|| Error::plain("negative Koszul power requires localization"))?
            } else {
                value
            };
            let coefficient = field.power(model, base, power.unsigned_abs())?;
            if coefficient == 0 {
                return Ok(0);
            }
            model.spend(field.degree())?;
            scalar = field.multiply(scalar, coefficient);
        }
        Ok(scalar)
    }
}

fn matrix(
    model: &mut Model,
    hom: &Hom,
    evaluator: &Evaluator,
    terms: &[Term],
) -> Result<BTreeMap<(usize, usize), usize>> {
    let mut values = BTreeMap::new();
    for term in terms {
        model.spend(term.cells())?;
        let coefficient = evaluator.coefficient(model, &term.powers)?;
        if coefficient != 0 {
            let source = hom.source.index(term.source)?;
            let target = hom.target.index(term.target)?;
            accumulate(model, &mut values, (source, target), coefficient)?;
        }
    }
    Ok(values)
}

fn accumulate<K: Ord>(
    model: &mut Model,
    values: &mut BTreeMap<K, usize>,
    index: K,
    coefficient: usize,
) -> Result<()> {
    model.spend(1 + values.len().max(1).ilog2() as usize)?;
    if let Some(previous) = values.get_mut(&index) {
        *previous ^= coefficient;
        if *previous == 0 {
            values.remove(&index);
        }
    } else {
        check_storage(3 * (1 + values.len()))?;
        values.insert(index, coefficient);
    }
    Ok(())
}
