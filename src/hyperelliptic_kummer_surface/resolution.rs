use std::cmp::Ordering;
use std::rc::Rc;

use super::{field, Point};
use crate::ec_poly::{self as poly, Field, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::model::{Model, Value};

pub(crate) struct ResolvedPoint {
    base: Rc<Point>,
    direction: Option<[Scalar; 4]>,
}

impl ResolvedPoint {
    pub(crate) fn smooth(model: &mut Model, curve: &Curve, base: Rc<Point>) -> Result<Self> {
        if base.is_node(model, curve)? {
            return Err(Error::plain(
                "rvKummerSmooth requires a nonsingular Kummer point",
            ));
        }
        Ok(Self {
            base,
            direction: None,
        })
    }

    pub(crate) fn exceptional(
        model: &mut Model,
        curve: &Curve,
        base: Rc<Point>,
        mut direction: [Scalar; 4],
    ) -> Result<Self> {
        if !base.is_node(model, curve)? {
            return Err(Error::plain("rvKummerExceptional requires a Kummer node"));
        }
        let field = field(curve);
        let pivot = pivot(&base);
        let shift = direction[pivot].clone();
        for (coordinate, base_coordinate) in direction.iter_mut().zip(base.coordinates()) {
            let correction = field.multiply(model, &shift, base_coordinate)?;
            *coordinate = field.subtract(model, coordinate, &correction)?;
        }
        let first = direction
            .iter()
            .find(|coordinate| !poly::zero(coordinate))
            .ok_or_else(|| {
                Error::plain("rvKummerExceptional requires a nonzero projective tangent direction")
            })?;
        if *first != field.one() {
            let inverse = field.inverse_unit(model, first)?;
            for coordinate in &mut direction {
                *coordinate = field.multiply(model, coordinate, &inverse)?;
            }
        }
        let quartic = quartic(model, curve)?;
        let cone = TangentCone::new(model, field, &base, &quartic)?;
        if !cone.contains(model, field, cone.axes.map(|axis| &direction[axis]))? {
            return Err(Error::plain(
                "rvKummerExceptional direction does not lie on the quadratic tangent cone",
            ));
        }
        Ok(Self {
            base,
            direction: Some(direction),
        })
    }

    pub(crate) fn base(&self) -> &Rc<Point> {
        &self.base
    }

    pub(crate) fn direction(&self) -> Option<&[Scalar; 4]> {
        self.direction.as_ref()
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        let order = self.base.compare(model, &other.base)?;
        if !order.is_eq() {
            return Ok(order);
        }
        if self.direction.is_some() && other.direction.is_some() {
            model.spend(4 * (1 + self.base.coordinates()[0].len()))?;
        }
        Ok(self.direction.cmp(&other.direction))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        let cells = self.direction.as_ref().map_or(0, |h| 4 * (3 + h[0].len()));
        model.allocate(2 + cells)?;
        Ok(Value::KummerK3(Rc::new(self)))
    }
}

fn pivot(base: &Point) -> usize {
    base.coordinates()
        .iter()
        .position(|coordinate| !poly::zero(coordinate))
        .expect("normalized projective point")
}

type QuarticTerm = ([usize; 4], Scalar);

const W: usize = 0;
const S: usize = 1;
const P: usize = 2;
const T: usize = 3;

fn quartic(model: &mut Model, curve: &Curve) -> Result<Vec<QuarticTerm>> {
    let field = field(curve);
    let f = curve.polynomial();
    let mut terms = vec![
        ([S, S, T, T], field.one()),
        ([W, P, T, T], field.signed(-4)),
    ];
    // The remaining terms containing T are −2Mt, homogenized to degree four.
    let middle = [
        ([W, W, W, T], -4),
        ([W, W, S, T], -2),
        ([W, W, P, T], -4),
        ([W, S, P, T], -2),
        ([W, P, P, T], -4),
        ([S, P, P, T], -2),
    ];
    for ((slots, multiplier), coefficient) in middle.into_iter().zip(f) {
        let coefficient = field.multiply(model, &field.signed(multiplier), coefficient)?;
        terms.push((slots, coefficient));
    }
    // These are the homogenized coefficients of (M²−4N)/(s²−4p), obtained
    // by expanding the five power sums of the two roots of x²−sx+p.
    let constant = [
        ([W, W, W, W], -4, 0, 2),
        ([W, W, W, W], 1, 1, 1),
        ([W, W, W, P], -2, 1, 3),
        ([W, W, P, P], 2, 1, 5),
        ([W, W, P, P], -4, 2, 4),
        ([W, W, P, P], 1, 3, 3),
        ([W, P, P, P], -2, 3, 5),
        ([P, P, P, P], 1, 5, 5),
        ([W, W, W, S], -4, 0, 3),
        ([W, W, S, P], 4, 0, 5),
        ([W, W, S, P], -4, 1, 4),
        ([W, S, P, P], -4, 2, 5),
        ([W, W, S, S], -4, 0, 4),
        ([W, S, S, P], -4, 1, 5),
        ([W, S, S, S], -4, 0, 5),
    ];
    for (slots, multiplier, left, right) in constant {
        let coefficient = field.multiply(model, &f[left], &f[right])?;
        let coefficient = field.multiply(model, &field.signed(multiplier), &coefficient)?;
        terms.push((slots, coefficient));
    }
    terms.retain(|(_, coefficient)| !poly::zero(coefficient));
    Ok(terms)
}

struct TangentCone {
    axes: [usize; 3],
    coefficients: [[Scalar; 3]; 3],
}

impl TangentCone {
    fn new(
        model: &mut Model,
        field: Field<'_>,
        base: &Point,
        quartic: &[QuarticTerm],
    ) -> Result<Self> {
        let pivot = pivot(base);
        let axes = std::array::from_fn(|i| i + usize::from(i >= pivot));
        let mut axis_index = [None; 4];
        for (index, &axis) in axes.iter().enumerate() {
            axis_index[axis] = Some(index);
        }
        let mut coefficients = std::array::from_fn(|_| std::array::from_fn(|_| field.zero()));
        let point = base.coordinates();
        for (slots, coefficient) in quartic {
            // Choose the two factors contributing direction coordinates. This
            // computes the second Hasse coefficient in every odd characteristic.
            model.spend(6)?;
            for first in 0..4 {
                for second in first + 1..4 {
                    let (Some(i), Some(j)) = (axis_index[slots[first]], axis_index[slots[second]])
                    else {
                        continue;
                    };
                    let mut term = coefficient.clone();
                    for (index, &axis) in slots.iter().enumerate() {
                        if index != first && index != second {
                            term = field.multiply(model, &term, &point[axis])?;
                        }
                    }
                    coefficients[i][j] = field.add(model, &coefficients[i][j], &term)?;
                }
            }
        }
        Ok(Self { axes, coefficients })
    }

    fn contains(
        &self,
        model: &mut Model,
        field: Field<'_>,
        direction: [&Scalar; 3],
    ) -> Result<bool> {
        let mut value = field.zero();
        model.spend(6)?;
        for (i, row) in self.coefficients.iter().enumerate() {
            for (coefficient, second) in row[i..].iter().zip(&direction[i..]) {
                if !poly::zero(coefficient) {
                    let product = field.multiply(model, direction[i], second)?;
                    let term = field.multiply(model, coefficient, &product)?;
                    value = field.add(model, &value, &term)?;
                }
            }
        }
        Ok(field.is_zero(&value))
    }
}

pub(crate) fn elements(model: &mut Model, curve: &Curve, limit: usize) -> Result<Vec<Value>> {
    let points = super::elements(model, curve, limit)?;
    let field = field(curve);
    let q = usize::try_from(field.ring.size_of_residue_field()?)
        .map_err(|_| Error::plain("Kummer resolution exceeds the finite element limit"))?;
    model.spend(q * field.cells())?;
    let scalars: Vec<_> = (0..q).map(|i| field.ring.coordinate_row(i)).collect();
    let quartic = quartic(model, curve)?;
    let zero = field.zero();
    let one = field.one();
    let mut output = Vec::new();
    let mut insert = |model: &mut Model, point: ResolvedPoint| -> Result<()> {
        if output.len() == limit {
            return Err(Error::plain(
                "Kummer resolution exceeds the finite element limit",
            ));
        }
        output.push(point.persist(model)?);
        Ok(())
    };
    for point in points {
        let Value::Kummer(base) = point else {
            unreachable!()
        };
        if !base.is_node(model, curve)? {
            insert(
                model,
                ResolvedPoint {
                    base,
                    direction: None,
                },
            )?;
            continue;
        }
        let cone = TangentCone::new(model, field, &base, &quartic)?;
        let mut visit = |model: &mut Model, coordinates: [&Scalar; 3]| -> Result<()> {
            if cone.contains(model, field, coordinates)? {
                let mut direction = std::array::from_fn(|_| zero.clone());
                for (axis, value) in cone.axes.into_iter().zip(coordinates) {
                    direction[axis] = value.clone();
                }
                insert(
                    model,
                    ResolvedPoint {
                        base: base.clone(),
                        direction: Some(direction),
                    },
                )?;
            }
            Ok(())
        };
        // Visit the normalized projective plane in lexicographic order.
        visit(model, [&zero, &zero, &one])?;
        for last in &scalars {
            visit(model, [&zero, &one, last])?;
        }
        for middle in &scalars {
            for last in &scalars {
                visit(model, [&one, middle, last])?;
            }
        }
    }
    Ok(output)
}
