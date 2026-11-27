use std::collections::BTreeMap;

use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;

use super::super::graded_basis::{complement_monomial, graded_basis};
use super::super::{
    grade_range, Accumulator, Grade, Integer, Monomial, State, COLORS, LIMBS, POINT,
};
use super::{cancel, factorial, normalize, positive, pure, Class, Denominator};

mod sparse;
use sparse::{Column, Echelon};

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Inverse,
    Lower,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Inverse => "rvHilbCohLefschetzInverse",
            Self::Lower => "rvHilbCohLower",
        }
    }

    pub(crate) fn output_degree(self, grade: &Grade) -> Result<i128> {
        match self {
            Self::Lower => grade.degree.checked_sub(2).ok_or_else(grade_range),
            Self::Inverse => {
                let degree = Integer::product(4, grade.length)
                    .add(Integer::from_i128(grade.degree).negate())?
                    .signed()
                    .ok_or_else(grade_range)?;
                if degree > grade.degree {
                    return Err(Error::plain(
                        "rvHilbCohLefschetzInverse requires a middle or upper cohomological degree",
                    ));
                }
                Ok(degree)
            }
        }
    }

    pub(crate) fn apply(
        self,
        model: &mut Model,
        h: &Class,
        x: &Class,
        grade: &Grade,
    ) -> Result<Class> {
        if !grade.supports_states() {
            return Ok(Class::zero());
        }
        if matches!(self, Self::Lower) {
            if grade.degree == 0 || is_zero(x) {
                return Ok(Class::zero());
            }
            if grade.degree == 2 {
                return lower_divisor(model, h, x, grade.length);
            }
        }
        let dimension = grade.length.checked_mul(2).ok_or_else(grade_range)?;
        let mut action = Multiplication {
            h,
            length: grade.length,
            squares: Vec::new(),
            cells: 0,
        };
        match self {
            Self::Inverse => action.inverse(model, x, grade, self.output_degree(grade)?, dimension),
            Self::Lower => action.lower(model, x, grade, dimension),
        }
    }
}

struct Multiplication<'a> {
    h: &'a Class,
    length: i128,
    squares: Vec<Class>,
    cells: usize,
}

impl Multiplication<'_> {
    fn power(&mut self, model: &mut Model, mut exponent: i128) -> Result<Class> {
        if exponent == 0 {
            return Class::unit(model, self.length as u128);
        }
        if is_zero(self.h) {
            return Ok(Class::zero());
        }
        if self.squares.is_empty() {
            self.squares.push(self.h.scale(model, 1, 1)?);
            self.cells = self.squares[0].cells();
        }
        let mut index = 0;
        let mut result: Option<Class> = None;
        while exponent > 0 {
            if index == self.squares.len() {
                let previous = self.squares.last().unwrap();
                let square = previous.cup(model, previous, self.length)?;
                self.cells += square.cells();
                guard_cells(self.cells)?;
                self.squares.push(square);
            }
            if exponent & 1 == 1 {
                let next = &self.squares[index];
                result = Some(if let Some(prefix) = result {
                    prefix.cup(model, next, self.length)?
                } else {
                    next.scale(model, 1, 1)?
                });
            }
            exponent /= 2;
            index += 1;
        }
        Ok(result.expect("positive exponent"))
    }

    fn raise(&mut self, model: &mut Model, x: &Class, exponent: i128) -> Result<Class> {
        if is_zero(x) {
            return Ok(Class::zero());
        }
        if exponent == 0 {
            return x.scale(model, 1, 1);
        }
        self.power(model, exponent)?.cup(model, x, self.length)
    }

    fn inverse(
        &mut self,
        model: &mut Model,
        x: &Class,
        grade: &Grade,
        low: i128,
        dimension: i128,
    ) -> Result<Class> {
        let exponent = dimension - low;
        if exponent == 0 {
            return x.scale(model, 1, 1);
        }
        if low == 0 && is_zero(x) {
            // Fujiki certifies the scalar map without constructing a high cup power.
            return if divisor_pair(model, self.h, self.h, self.length)?.is_zero() {
                Err(singular())
            } else {
                Ok(Class::zero())
            };
        }
        let operator = self.power(model, exponent)?;
        if is_zero(&operator) {
            return Err(singular());
        }
        if low == 0 {
            return invert_top(model, x, &operator, self.length as u128);
        }
        let source_grade = Grade {
            degree: low,
            ..grade.clone()
        };
        solve_slice(model, &source_grade, &operator, x)
    }

    fn lower(
        &mut self,
        model: &mut Model,
        x: &Class,
        grade: &Grade,
        dimension: i128,
    ) -> Result<Class> {
        if is_zero(self.h) {
            return Err(singular());
        }
        let first = (grade.degree - dimension).max(1);
        let mut remainder = x.scale(model, 1, 1)?;
        let mut output = Class::zero();
        for j in (first..=grade.degree / 2).rev() {
            if is_zero(&remainder) {
                break;
            }
            model.spend(1)?;
            let d = grade.degree - 2 * j;
            // This power kills every remaining string with a smaller j.
            let raised = self.raise(model, &remainder, dimension - grade.degree + j)?;
            let primitive = self.inverse(model, &raised, grade, d, dimension)?;
            if is_zero(&primitive) {
                continue;
            }
            let lowered = self.raise(model, &primitive, j - 1)?;
            let weight = Integer::product(j, dimension - d - j + 1);
            let term = Class::scaled_wide(
                model,
                &lowered.numerator,
                weight,
                Denominator::from_small(lowered.denominator),
            )?;
            output = output.add(model, &term)?;
            if j > first {
                let component = self.h.cup(model, &lowered, self.length)?;
                remainder = remainder.add_signed(model, &component, true)?;
            }
        }
        Ok(output)
    }
}

fn solve_slice(model: &mut Model, source: &Grade, operator: &Class, x: &Class) -> Result<Class> {
    let basis = graded_basis(model, source)?;
    let size = basis.len();
    let basis_cells = basis
        .iter()
        .map(|monomial| 12 + 10 * monomial.len())
        .sum::<usize>();
    guard_cells(basis_cells.saturating_add(size.saturating_mul(2 * (LIMBS + 4))))?;
    model
        .spend(size.saturating_mul(size))
        .map_err(|_| work_limit())?;
    let mut target = BTreeMap::new();
    for (index, monomial) in basis.iter().enumerate() {
        let mut dual = monomial.clone();
        complement_monomial(&mut dual);
        target.insert(dual, index);
    }
    let mut echelon = Echelon::default();
    for (index, monomial) in basis.iter().enumerate() {
        let state = State {
            terms: BTreeMap::from([(monomial.clone(), 1)]),
        };
        let image = operator
            .numerator
            .accumulate_cup(model, &state, source.length)?;
        let column = Column {
            image: numbered(model, &target, image.terms)?,
            witness: vec![(index, positive(1))],
        };
        if !echelon.insert(model, column)? {
            return Err(singular());
        }
    }
    if is_zero(x) {
        return Ok(Class::zero());
    }
    let image = numbered(
        model,
        &target,
        x.numerator
            .terms
            .iter()
            .map(|(key, &value)| (key, Integer::from_i128(value))),
    )?;
    let relation = echelon.solve(
        model,
        Column {
            image,
            witness: vec![(size, positive(1))],
        },
    )?;
    if !relation.image.is_empty() {
        return Err(singular());
    }
    let mut witness = relation.witness;
    let (_, mut last) = witness
        .pop()
        .expect("the right-hand side has a nonzero witness");
    let mut op = Denominator::from_small(operator.denominator);
    let mut rhs = Denominator::from_small(x.denominator);
    cancel(model, &mut op, &mut rhs)?;
    cancel(model, &mut op, &mut last.magnitude)?;
    let scalar = Integer {
        magnitude: op,
        negative: !last.negative,
    };
    let denominator = last.magnitude.checked_product(rhs).ok_or_else(arithmetic)?;
    let mut output = Accumulator::default();
    for (index, coefficient) in witness {
        model.spend(1 + basis[index].len())?;
        output.add(
            basis[index].clone(),
            coefficient
                .multiply_wide(scalar)
                .map_err(|_| arithmetic())?,
        )?;
    }
    normalize(model, output, denominator)
}

fn numbered<K: std::borrow::Borrow<Monomial>>(
    model: &mut Model,
    target: &BTreeMap<Monomial, usize>,
    terms: impl IntoIterator<Item = (K, Integer)>,
) -> Result<sparse::Vector> {
    let mut result = Vec::new();
    for (key, coefficient) in terms {
        if coefficient.is_zero() {
            continue;
        }
        let key = key.borrow();
        model.spend(1 + key.len())?;
        let column = *target.get(key).ok_or_else(|| {
            Error::plain("Lefschetz cup image lies outside the complementary Fock basis")
        })?;
        result.push((column, coefficient));
    }
    result.sort_unstable_by_key(|&(coordinate, _)| coordinate);
    Ok(result)
}

fn invert_top(model: &mut Model, x: &Class, operator: &Class, n: u128) -> Result<Class> {
    let image = operator
        .numerator
        .terms
        .values()
        .next()
        .copied()
        .ok_or_else(singular)?;
    let Some(&target) = x.numerator.terms.values().next() else {
        return Ok(Class::zero());
    };
    let numerator = [Integer::from_i128(target), positive(operator.denominator)];
    let denominator = [
        Denominator::from_small(image.unsigned_abs()),
        Denominator::from_small(x.denominator),
        factorial(model, n)?,
    ];
    scalar_unit(model, n, numerator, denominator, image < 0)
}

fn lower_divisor(model: &mut Model, h: &Class, x: &Class, length: i128) -> Result<Class> {
    let square = divisor_pair(model, h, h, length)?;
    if square.is_zero() {
        return Err(singular());
    }
    let pairing = divisor_pair(model, h, x, length)?;
    if pairing.is_zero() {
        return Ok(Class::zero());
    }
    // Polarized Fujiki identifies the degree-two primitive space with h-perp.
    // Thus Λ_h(x) = (2n) q(h,x)/q(h) in H^0.
    let numerator = [
        Integer::product(2, length),
        pairing,
        positive(h.denominator),
    ];
    let n = length as u128;
    let denominator = [
        square.magnitude,
        Denominator::from_small(x.denominator),
        factorial(model, n)?,
    ];
    scalar_unit(model, n, numerator, denominator, square.negative)
}

pub(super) fn divisor_pair(
    model: &mut Model,
    h: &Class,
    x: &Class,
    length: i128,
) -> Result<Integer> {
    let coordinates = |class: &Class| {
        let mut coords = [0i128; COLORS];
        for (monomial, &value) in &class.numerator.terms {
            let mode = monomial.last().expect("degree-two mode").mode;
            coords[usize::from(if mode.length == 2 { POINT } else { mode.color })] = value;
        }
        coords
    };
    let left = coordinates(h);
    let right = coordinates(x);
    let gram = model.nakajima_gram()?;
    let mut surface = Integer::from_i128(0);
    model.spend(COLORS * COLORS * LIMBS)?;
    for a in 1..POINT {
        for b in 1..POINT {
            let product = Integer::product(left[usize::from(a)], right[usize::from(b)])
                .multiply(gram.entry(a, b))
                .map_err(|_| arithmetic())?;
            surface = surface.add(product).map_err(|_| arithmetic())?;
        }
    }
    if length == 1 {
        return Ok(surface);
    }
    // In Fock coordinates μ(a) has denominator (n−1)! and δ has denominator
    // 2(n−2)!, with q(δ)=−2(n−1). The common factor cancels in the ratio.
    let boundary = Integer::product(left[usize::from(POINT)], right[usize::from(POINT)])
        .multiply(-8)
        .map_err(|_| arithmetic())?;
    surface
        .multiply(length - 1)
        .and_then(|scaled| scaled.add(boundary))
        .map_err(|_| arithmetic())
}

fn scalar_unit<const A: usize, const B: usize>(
    model: &mut Model,
    n: u128,
    mut numerator: [Integer; A],
    mut denominator: [Denominator; B],
    negative: bool,
) -> Result<Class> {
    for value in &mut numerator {
        for divisor in &mut denominator {
            cancel(model, &mut value.magnitude, divisor)?;
        }
    }
    let mut scale = positive(1);
    scale.negative = negative;
    for factor in numerator {
        scale = scale.multiply_wide(factor).map_err(|_| arithmetic())?;
    }
    let mut divisor = Denominator::from_small(1);
    for factor in denominator {
        divisor = divisor.checked_product(factor).ok_or_else(arithmetic)?;
    }
    Class::scaled_wide(model, &pure(n, 0), scale, divisor)
}

fn is_zero(class: &Class) -> bool {
    class.numerator.terms.is_empty()
}

fn guard_cells(cells: usize) -> Result<()> {
    if cells > MAX_PIVOT_CELLS {
        Err(work_limit())
    } else {
        Ok(())
    }
}

fn work_limit() -> Error {
    Error::plain("Hilbert Lefschetz computation exceeds the graded matrix work or storage limit")
}

fn arithmetic() -> Error {
    Error::plain("Hilbert Lefschetz computation exceeds the wide exact arithmetic range")
}

fn singular() -> Error {
    Error::plain("Hilbert Lefschetz cup map is singular on the requested slice")
}
