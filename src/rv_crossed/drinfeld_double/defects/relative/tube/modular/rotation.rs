use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_lens::Rational;

use super::super::super::linear::{self, GradedMap};
use super::super::super::tensor::Layout;
use super::super::super::tree::{self, Tree};
use super::super::super::{decompose, Ambient, Balanced, Category, Module};
use super::super::build::{self, Builder};
use super::super::{Double, Tube};
use super::{Field, Matrix, Torus};

const MAX_ROTATION_FIBERS: usize = 4096;
const MAX_ROTATION_CHANNELS: usize = 256;

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Rotation {
    simple: usize,
    dual: usize,
    power: usize,
    matrix: Matrix,
    trace: Rational,
}

impl Rotation {
    pub(crate) fn data(&self) -> String {
        format!("rvDoubleTubeRotationData = simple={}, dual={}, power={}, dimension={}, trace={:?}/{}, verified=true",
            self.simple,self.dual,self.power,self.matrix.rows,self.trace.numerator,self.trace.denominator)
    }

    pub(crate) fn entry(&self, key: [usize; 2]) -> Result<Rational> {
        self.matrix.entry(key)
    }

    pub(crate) fn trace(&self) -> Rational {
        self.trace.clone()
    }
}

impl Tube {
    pub(crate) fn rotate(
        &self,
        model: &mut Model,
        double: &Double,
        category: &Category,
        torus: &Torus,
        [i, n]: [usize; 2],
    ) -> Result<Value> {
        if !(2..=4).contains(&n) {
            return Err(Error::plain(
                "direct pivotal rotation supports indicator degrees 2 through 4",
            ));
        }
        let expected = self.indicator(model, torus, [i, n])?;
        let mut work = Builder::new(model, double, category)?;
        let unit_index = self.blocks[0].slots[0][0];
        let x = work.objects[i].clone();
        let unit = work.objects[unit_index].clone();
        let h = category.boundaries[0].members()?.len();
        let dimension = x.dimension() / h;
        if dimension
            .checked_pow(n as u32)
            .and_then(|d| d.checked_mul(h))
            .is_none_or(|d| d > MAX_ROTATION_FIBERS)
        {
            return Err(Error::plain(
                "direct pivotal rotation exceeds its graded tensor bound",
            ));
        }
        let mut dual = None;
        for j in 0..work.objects.len() {
            if !work.channels([j, i, unit_index])?.is_empty() && dual.replace(j).is_some() {
                return Err(Error::plain("boundary simple has more than one dual"));
            }
        }
        let dual = dual.ok_or_else(|| Error::plain("boundary simple has no dual unit channel"))?;
        let d = work.objects[dual].clone();
        let right = work.channels([dual, i, unit_index])?;
        let left = work.channels([i, dual, unit_index])?;
        let ([right], [left]) = (right.as_ref(), left.as_ref()) else {
            return Err(Error::plain(
                "pivotal dual requires one evaluation channel in each orientation",
            ));
        };
        let coev_r = &right.inclusion;
        let mut ev_l = right.projection.clone();
        let dimension_scalar = Rational::integer(dimension as i128);
        ev_l.scale(&mut work.field, &dimension_scalar)?;
        let mut coev_l = left.inclusion.clone();
        let mut ev_r = left.projection.clone();
        ev_r.scale(&mut work.field, &dimension_scalar)?;
        let scalar = snake(
            &mut work.field,
            &work.ambient,
            [&x, &d, &unit],
            [&coev_l, &ev_l],
            0,
        )?;
        let scalar = scalar_endomorphism(&mut work.field, &x, &scalar)?;
        let inv = work.field.inverse(&scalar)?;
        coev_l.scale(&mut work.field, &inv)?;
        ev_r.scale(&mut work.field, &scalar)?;
        let ids = [GradedMap::identity(&x)?, GradedMap::identity(&d)?];
        for (objects, duals, id) in [
            (
                [x.as_ref(), d.as_ref(), unit.as_ref()],
                [&coev_l, &ev_l, coev_r, &ev_r],
                &ids[0],
            ),
            (
                [d.as_ref(), x.as_ref(), unit.as_ref()],
                [coev_r, &ev_r, &coev_l, &ev_l],
                &ids[1],
            ),
        ] {
            for side in 0..2 {
                let start = 2 * side;
                let check = snake(
                    &mut work.field,
                    &work.ambient,
                    objects,
                    [duals[start], duals[start + 1]],
                    side,
                )?;
                if check.grades != id.grades {
                    return Err(Error::plain(
                        "pivotal dual cups and caps do not satisfy both snake identities",
                    ));
                }
            }
        }
        let mut tail: Vec<Balanced> = Vec::new();
        for _ in 2..n {
            let y = tail.last().map_or(x.as_ref(), |p| &p.module);
            tail.push(Balanced::new(work.field.model, &work.ambient, &x, y)?);
        }
        let y = tail.last().map_or(x.as_ref(), |p| &p.module);
        let range = Balanced::new(work.field.model, &work.ambient, &x, y)?;
        let moved = Balanced::new(work.field.model, &work.ambient, y, &x)?;
        let homs = decompose(work.field.model, &work.ambient, category, &range.module)?;
        if homs[unit_index].orbits.len() > MAX_ROTATION_CHANNELS {
            return Err(Error::plain(
                "direct pivotal rotation exceeds its invariant Hom bound",
            ));
        }
        let channels = linear::hom_basis(
            &mut work.field,
            &work.ambient,
            &unit,
            &range.module,
            category.simple(unit_index)?.simple.flux,
            &homs[unit_index],
        )?;
        let alignment = align(&mut work.field, &work.ambient, &x, n, [&moved, &range])?;
        let kernel = Mate::new(
            &mut work.field,
            &work.ambient,
            [&x, &d, y, &unit],
            [&range, &moved],
        )?;
        let mut first_to_last = Matrix::zero(channels.len(), channels.len())?;
        for (j, channel) in channels.iter().enumerate() {
            let vector = kernel.rotate(
                &mut work.field,
                &work.ambient,
                &channel.inclusion,
                coev_r,
                &ev_l,
            )?;
            let vector = alignment.compose(&mut work.field, &vector)?;
            let mut check = GradedMap::zero(&unit, &range.module)?;
            for (r, basis) in channels.iter().enumerate() {
                let a = basis.inclusion.inner(&mut work.field, &vector)?;
                let c = work.field.multiply(&a, &basis.inverse_norm)?;
                first_to_last.set(r, j, c.clone());
                let mut term = basis.inclusion.clone();
                term.scale(&mut work.field, &c)?;
                check.add_assign(&mut work.field, &term)?;
            }
            if check.grades != vector.grades {
                return Err(Error::plain(
                    "pivotal rotation escaped the invariant Hom space",
                ));
            }
        }
        // The positive twist convention is paired with the last-to-first
        // rotation. The two directions differ in the odd carry test.
        let matrix = first_to_last.inverse(&mut work.field)?;
        let identity = Matrix::identity(matrix.rows)?;
        let mut trace = Rational::zero();
        for i in 0..matrix.rows {
            trace = work.field.add(&trace, matrix.get(i, i))?;
        }
        let mut power = identity.clone();
        for _ in 0..n {
            power = matrix.multiply(&mut work.field, &power)?;
        }
        if power != identity || trace != expected {
            return Err(Error::plain(format!(
                "pivotal rotation fails its order or indicator trace check at ({i},{n})"
            )));
        }
        work.field
            .model
            .allocate(matrix.cells() + 2 * trace.numerator.len() + 8)?;
        Ok(Value::DoubleTubeRotation(Rc::new(Rotation {
            simple: i,
            dual,
            power: n,
            matrix,
            trace,
        })))
    }
}

fn scalar_endomorphism(
    field: &mut Field<'_>,
    module: &Module,
    value: &GradedMap,
) -> Result<Rational> {
    let identity = GradedMap::identity(module)?;
    let trace = identity.inner(field, value)?;
    let size = Rational::integer(module.dimension() as i128);
    let inverse = field.inverse(&size)?;
    let scalar = field.multiply(&trace, &inverse)?;
    let mut target = identity;
    target.scale(field, &scalar)?;
    if target.grades != value.grades {
        return Err(Error::plain("pivotal snake is not a scalar endomorphism"));
    }
    Ok(scalar)
}

fn snake(
    field: &mut Field<'_>,
    ambient: &Ambient<'_>,
    [x, d, u]: [&Module; 3],
    [born, cap]: [&GradedMap; 2],
    side: usize,
) -> Result<GradedMap> {
    let inputs = [x, d, x];
    let [a, b] = tree::TRIPLE
        .each_ref()
        .map(|shape| Tree::build(field.model, ambient, &inputs, shape));
    let [left, right] = [a?, b?];
    let routes = tree::rotation(field.model, ambient, &inputs, &left, &right, [0, 1, 2, 3])?;
    let forward = GradedMap::monomial(field, &routes)?;
    let plans = [&left, &right];
    let [from, to] = [plans[side].layout(), plans[1 - side].layout()];
    let id = GradedMap::identity(x)?;
    let first_factors = if side == 0 { [u, x] } else { [x, u] };
    let last_factors = if side == 0 { [x, u] } else { [u, x] };
    let first = Balanced::new(field.model, ambient, first_factors[0], first_factors[1])?;
    let last = Balanced::new(field.model, ambient, last_factors[0], last_factors[1])?;
    let first = Layout {
        product: &first,
        factors: first_factors,
    };
    let last = Layout {
        product: &last,
        factors: last_factors,
    };
    let born_factors = if side == 0 { [born, &id] } else { [&id, born] };
    let cap_factors = if side == 0 { [&id, cap] } else { [cap, &id] };
    let inc = linear::tensor(field, ambient, first, from, born_factors)?;
    let end = linear::tensor(field, ambient, to, last, cap_factors)?;
    let left_unit = unitor(field, ambient, first, x, side)?.adjoint(field)?;
    let right_unit = unitor(field, ambient, last, x, 1 - side)?;
    let assoc = if side == 0 {
        forward
    } else {
        forward.adjoint(field)?
    };
    let m = inc.compose(field, &left_unit)?;
    let m = assoc.compose(field, &m)?;
    let m = end.compose(field, &m)?;
    right_unit.compose(field, &m)
}

struct Mate<'a> {
    objects: [&'a Module; 4],
    pair: [&'a Balanced; 2],
    dx: Balanced,
    du: Balanced,
    uy: Balanced,
    first: Balanced,
    second: Balanced,
    associator: GradedMap,
    rho_inv: GradedMap,
    lambda: GradedMap,
}

impl<'a> Mate<'a> {
    fn new(
        field: &mut Field<'_>,
        ambient: &Ambient<'_>,
        objects: [&'a Module; 4],
        pair: [&'a Balanced; 2],
    ) -> Result<Self> {
        let [x, d, y, u] = objects;
        let leaves = [d, x, y];
        let [l, r] = tree::TRIPLE
            .each_ref()
            .map(|shape| Tree::build(field.model, ambient, &leaves, shape));
        let [l, r] = [l?, r?];
        let rotation = tree::rotation(field.model, ambient, &leaves, &l, &r, [0, 1, 2, 3])?;
        let associator = GradedMap::monomial(field, &rotation)?.adjoint(field)?;
        let dx = Balanced::new(field.model, ambient, d, x)?;
        let du = Balanced::new(field.model, ambient, d, u)?;
        let uy = Balanced::new(field.model, ambient, u, y)?;
        let first = Balanced::new(field.model, ambient, d, &pair[0].module)?;
        let second = Balanced::new(field.model, ambient, &dx.module, y)?;
        if first.bases != r.layout().product.bases || second.bases != l.layout().product.bases {
            return Err(Error::plain(
                "pivotal rotation disagrees with its balanced associator bases",
            ));
        }
        let rho_inv = unitor(
            field,
            ambient,
            Layout {
                product: &du,
                factors: [d, u],
            },
            d,
            1,
        )?
        .adjoint(field)?;
        let lambda = unitor(
            field,
            ambient,
            Layout {
                product: &uy,
                factors: [u, y],
            },
            y,
            0,
        )?;
        Ok(Self {
            objects,
            pair,
            dx,
            du,
            uy,
            first,
            second,
            associator,
            rho_inv,
            lambda,
        })
    }

    fn rotate(
        &self,
        field: &mut Field<'_>,
        ambient: &Ambient<'_>,
        vector: &GradedMap,
        born: &GradedMap,
        cap: &GradedMap,
    ) -> Result<GradedMap> {
        let [x, d, y, u] = self.objects;
        let d0 = GradedMap::identity(d)?;
        let y0 = GradedMap::identity(y)?;
        let x0 = GradedMap::identity(x)?;
        let inc = linear::tensor(
            field,
            ambient,
            Layout {
                product: &self.du,
                factors: [d, u],
            },
            Layout {
                product: &self.first,
                factors: [d, &self.pair[0].module],
            },
            [&d0, vector],
        )?;
        let eval = linear::tensor(
            field,
            ambient,
            Layout {
                product: &self.second,
                factors: [&self.dx.module, y],
            },
            Layout {
                product: &self.uy,
                factors: [u, y],
            },
            [cap, &y0],
        )?;
        let mate = inc.compose(field, &self.rho_inv)?;
        let mate = self.associator.compose(field, &mate)?;
        let mate = eval.compose(field, &mate)?;
        let mate = self.lambda.compose(field, &mate)?;
        let map = linear::tensor(
            field,
            ambient,
            Layout {
                product: &self.dx,
                factors: [d, x],
            },
            Layout {
                product: self.pair[1],
                factors: [y, x],
            },
            [&mate, &x0],
        )?;
        map.compose(field, born)
    }
}

fn unitor(
    field: &mut Field<'_>,
    ambient: &Ambient<'_>,
    layout: Layout<'_>,
    module: &Module,
    side: usize,
) -> Result<GradedMap> {
    let rows = build::unitor(ambient, layout, module, side);
    for (row, &dim) in rows.iter().zip(&module.fibers) {
        if row.len() != dim {
            return Err(Error::plain(
                "pivotal unitors do not preserve graded dimensions",
            ));
        }
        let mut seen = vec![false; dim];
        for image in row {
            if image.target >= dim || std::mem::replace(&mut seen[image.target], true) {
                return Err(Error::plain("pivotal unitor is not a graded bijection"));
            }
        }
    }
    GradedMap::monomial(field, &rows)
}

fn align(
    field: &mut Field<'_>,
    ambient: &Ambient<'_>,
    x: &Module,
    n: usize,
    [left, right]: [&Balanced; 2],
) -> Result<GradedMap> {
    if n == 2 {
        if left.bases != right.bases {
            return Err(Error::plain(
                "twofold pivotal rotation has different pair bases",
            ));
        }
        return GradedMap::identity(&left.module);
    }
    let leaves = vec![x; n];
    let (shapes, moves) = if n == 3 {
        (vec![&tree::TRIPLE[0], &tree::TRIPLE[1]], vec![[0, 1, 2, 3]])
    } else {
        (
            vec![&tree::FOUR[1], &tree::FOUR[2], &tree::FOUR[3]],
            vec![[0, 1, 3, 4], [1, 2, 3, 4]],
        )
    };
    let trees = shapes
        .into_iter()
        .map(|shape| Tree::build(field.model, ambient, &leaves, shape))
        .collect::<Result<Vec<_>>>()?;
    if left.bases != trees[0].layout().product.bases
        || right.bases != trees.last().unwrap().layout().product.bases
    {
        return Err(Error::plain(
            "pivotal rotation does not return to the chosen tensor basis",
        ));
    }
    let mut map = GradedMap::identity(&left.module)?;
    for (edge, phase) in trees.windows(2).zip(moves) {
        let rotation = tree::rotation(field.model, ambient, &leaves, &edge[0], &edge[1], phase)?;
        map = GradedMap::monomial(field, &rotation)?.compose(field, &map)?;
    }
    Ok(map)
}
