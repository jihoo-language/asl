use std::rc::Rc;

use crate::ec_poly::{Field, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::hw_closed::Curve;
use crate::jacobian_level::Basis;
use crate::model::Model;

use super::{require_curve, Point};

const WIDTH: usize = 4;
const NODE_COUNT: usize = 16;

pub(crate) enum Query {
    Entry(usize, usize, usize),
    Factor(usize, usize),
    Trope(usize, usize),
}

pub(crate) struct Group {
    pub(crate) curve: Rc<Curve>,
    matrices: Box<[Matrix<Scalar>]>,
    tropes: Box<[[Scalar; WIDTH]]>,
}

impl Group {
    pub(crate) fn build(model: &mut Model, basis: &Basis, curve: Rc<Curve>) -> Result<Self> {
        if curve.grade.genus != 2 || basis.prime.get() != 2 {
            return Err(Error::plain(
                "Kummer theta operations require a certified genus-two level-2 basis",
            ));
        }
        require_curve(&curve)?;
        let field = Field { ring: &curve.ring };
        let nodes = (0..NODE_COUNT)
            .map(|mask| {
                let divisor = basis.binary_combination(model, mask as u8)?;
                Point::project(model, &curve, &divisor)
            })
            .collect::<Result<Vec<_>>>()?;
        let anchors = anchors(field, &nodes)?;
        let source = columns(&nodes, anchors)?;
        let source_inverse = inverse(model, field, &source)?;
        let source_coordinates = nodes
            .iter()
            .map(|node| apply(model, field, &source_inverse, node.coordinates()))
            .collect::<Result<Vec<_>>>()?;

        let mut matrices = Vec::with_capacity(NODE_COUNT);
        let mut tropes = Vec::with_capacity(NODE_COUNT);
        for shift in 0..NODE_COUNT {
            let target = columns(&nodes, anchors.map(|index| index ^ shift))?;
            let target_inverse = inverse(model, field, &target)?;
            let target_coordinates = (0..NODE_COUNT)
                .map(|mask| {
                    apply(
                        model,
                        field,
                        &target_inverse,
                        nodes[mask ^ shift].coordinates(),
                    )
                })
                .collect::<Result<Vec<_>>>()?;
            let scales = solve_scales(model, field, &source_coordinates, &target_coordinates)?;
            let diagonal = (0..WIDTH * WIDTH)
                .map(|index| {
                    if index / WIDTH == index % WIDTH {
                        scales[index % WIDTH].clone()
                    } else {
                        field.zero()
                    }
                })
                .collect();
            // The two frames are invertible and every diagonal entry is nonzero.
            let diagonal = Matrix::from_entries(WIDTH, diagonal)?;
            let change = hw_algebra::multiply(model, &field, &diagonal, &source_inverse)?;
            let matrix = hw_algebra::multiply(model, &field, &target, &change)?;
            let mut entries = matrix.entries().to_vec();
            normalize(model, field, &mut entries)?;
            let matrix = Matrix::from_entries(WIDTH, entries)?;
            for (mask, node) in nodes.iter().enumerate() {
                let mut image = apply(model, field, &matrix, node.coordinates())?;
                normalize(model, field, &mut image)?;
                if image != *nodes[mask ^ shift].coordinates() {
                    return Err(Error::plain(
                        "Kummer theta matrix failed to translate a node",
                    ));
                }
            }
            // The first row is the pullback of w=0. Translation is its own
            // projective inverse, so this is also the translated trope.
            let mut trope = std::array::from_fn(|i| matrix.entries()[i].clone());
            normalize(model, field, &mut trope)?;
            tropes.push(trope);
            matrices.push(matrix);
        }
        model.allocate(NODE_COUNT * (5 + WIDTH * (WIDTH + 1) * field.cells()))?;
        Ok(Self {
            curve,
            matrices: matrices.into_boxed_slice(),
            tropes: tropes.into_boxed_slice(),
        })
    }

    pub(crate) fn matrix(&self, mask: usize) -> &Matrix<Scalar> {
        &self.matrices[mask]
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        mask: usize,
        coordinates: &[Scalar; WIDTH],
    ) -> Result<[Scalar; WIDTH]> {
        let field = Field {
            ring: &self.curve.ring,
        };
        apply(model, field, self.matrix(mask), coordinates)
    }

    pub(crate) fn translate(&self, model: &mut Model, mask: usize, point: &Point) -> Result<Point> {
        let coordinates = self.apply(model, mask, point.coordinates())?;
        Point::new(model, &self.curve, coordinates)
    }

    pub(crate) fn query(&self, model: &mut Model, query: Query) -> Result<Scalar> {
        match query {
            Query::Entry(mask, row, column) => {
                Ok(self.matrix(mask).entries()[row * WIDTH + column].clone())
            }
            Query::Trope(mask, index) => Ok(self.tropes[mask][index].clone()),
            Query::Factor(left, right) => self.factor(model, left, right),
        }
    }

    fn factor(&self, model: &mut Model, left: usize, right: usize) -> Result<Scalar> {
        let field = Field {
            ring: &self.curve.ring,
        };
        let product = hw_algebra::multiply(model, &field, self.matrix(left), self.matrix(right))?;
        let normalized = self.matrix(left ^ right).entries();
        let first = normalized
            .iter()
            .position(|c| !field.is_zero(c))
            .expect("an invertible matrix is nonzero");
        let factor = &product.entries()[first];
        if field.is_zero(factor) {
            return Err(Error::plain("Kummer theta product has zero central factor"));
        }
        for (actual, expected) in product.entries().iter().zip(normalized) {
            if *actual != field.multiply(model, factor, expected)? {
                return Err(Error::plain(
                    "Kummer theta product failed its projective relation",
                ));
            }
        }
        Ok(factor.clone())
    }
}

fn anchors(field: Field<'_>, nodes: &[Point]) -> Result<[usize; WIDTH]> {
    // The identity and two finite branch points span the plane w=0. Any
    // node with degree-two support completes the linear frame.
    let mut branch = nodes.iter().enumerate().filter(|(_, node)| {
        let [w, x, ..] = node.coordinates();
        field.is_zero(w) && !field.is_zero(x)
    });
    let interior = nodes
        .iter()
        .position(|node| !field.is_zero(&node.coordinates()[0]));
    match (branch.next(), branch.next(), interior) {
        (Some((first, _)), Some((second, _)), Some(interior)) => Ok([0, first, second, interior]),
        _ => Err(Error::plain("Kummer nodes do not contain a linear frame")),
    }
}

fn columns(nodes: &[Point], indices: [usize; WIDTH]) -> Result<Matrix<Scalar>> {
    Matrix::from_entries(
        WIDTH,
        (0..WIDTH)
            .flat_map(|row| indices.map(|index| nodes[index].coordinates()[row].clone()))
            .collect(),
    )
}

fn inverse(model: &mut Model, field: Field<'_>, matrix: &Matrix<Scalar>) -> Result<Matrix<Scalar>> {
    hw_algebra::inverse(model, &field, matrix)?
        .ok_or_else(|| Error::plain("Kummer theta node frame is singular"))
}

fn solve_scales(
    model: &mut Model,
    field: Field<'_>,
    source: &[[Scalar; WIDTH]],
    target: &[[Scalar; WIDTH]],
) -> Result<[Scalar; WIDTH]> {
    for (x, y) in source.iter().zip(target) {
        if x.iter()
            .zip(y)
            .any(|(a, b)| field.is_zero(a) != field.is_zero(b))
        {
            return Err(Error::plain(
                "Kummer theta node coordinates have incompatible supports",
            ));
        }
    }
    // In these frames the unknown map is diagonal. The equations
    // d_i x_i y_j = d_j x_j y_i fix the ratios on each support. Their
    // kernel is one-dimensional precisely when propagation reaches all
    // four entries. The full node test checks consistency of all ratios.
    let mut scales = [Some(field.one()), None, None, None];
    for _ in 1..WIDTH {
        for (x, y) in source.iter().zip(target) {
            let Some(pivot) = (0..WIDTH).find(|&i| scales[i].is_some() && !field.is_zero(&x[i]))
            else {
                continue;
            };
            for i in 0..WIDTH {
                if scales[i].is_some() || field.is_zero(&x[i]) {
                    continue;
                }
                let numerator =
                    field.multiply(model, scales[pivot].as_ref().unwrap(), &x[pivot])?;
                let numerator = field.multiply(model, &numerator, &y[i])?;
                let denominator = field.multiply(model, &y[pivot], &x[i])?;
                let inverse = field.inverse_unit(model, &denominator)?;
                scales[i] = Some(field.multiply(model, &numerator, &inverse)?);
            }
        }
    }
    let [Some(a), Some(b), Some(c), Some(d)] = scales else {
        return Err(Error::plain(
            "Kummer theta node equations must have a one-dimensional kernel",
        ));
    };
    Ok([a, b, c, d])
}

fn apply(
    model: &mut Model,
    field: Field<'_>,
    matrix: &Matrix<Scalar>,
    vector: &[Scalar; WIDTH],
) -> Result<[Scalar; WIDTH]> {
    let mut image = std::array::from_fn(|_| field.zero());
    for (row, coordinate) in matrix.entries().chunks_exact(WIDTH).zip(&mut image) {
        for (a, b) in row.iter().zip(vector) {
            let term = field.multiply(model, a, b)?;
            *coordinate = field.add(model, coordinate, &term)?;
        }
    }
    Ok(image)
}

fn normalize(model: &mut Model, field: Field<'_>, values: &mut [Scalar]) -> Result<()> {
    let scale = values
        .iter()
        .find(|c| !field.is_zero(c))
        .ok_or_else(|| Error::plain("Kummer theta projective coordinates cannot all vanish"))?;
    if scale != &field.one() {
        let inverse = field.inverse_unit(model, scale)?;
        for value in values {
            *value = field.multiply(model, value, &inverse)?;
        }
    }
    Ok(())
}
