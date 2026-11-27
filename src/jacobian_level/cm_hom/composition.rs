use std::rc::Rc;

use super::Hom;
use crate::ec_morphism::{self, Morphism};
use crate::elliptic::Curve;
use crate::elliptic_cm_recognition as maps;
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hyperelliptic_kummer_surface::cm_ring::Order;
use crate::model::{Model, Value};
use crate::rv_ore::gcd;
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Statement {
    Dual(Label, Label, Option<Box<[ParsedExpr; 2]>>),
    Compose(Label, Label, Label, Label),
    Report(Label),
    Coordinates(Box<ParsedExpr>),
}

pub(crate) enum Expression {
    Coordinates(Label, Box<ParsedExpr>),
    Literal(Label, [i128; 2]),
    Map(Box<ParsedExpr>),
    Compose(Label, Box<[ParsedExpr; 2]>),
    Inverse(Label, Box<ParsedExpr>),
    Product(Label, usize, bool),
}

pub(crate) struct Coordinates {
    pub(crate) numerator: [i128; 2],
    pub(crate) denominator: i128,
    pub(crate) map: Rc<Morphism>,
}

impl Coordinates {
    pub(crate) fn read(model: &mut Model, frame: &Hom, map: Rc<Morphism>) -> Result<Self> {
        let [first, second] = &frame.maps;
        if (map.source.as_ref(), map.target.as_ref())
            != (first.source.as_ref(), first.target.as_ref())
        {
            return Err(Error::plain(
                "CM Hom coordinate map has different named endpoints from the frame",
            ));
        }
        let [domain, range] = frame.endpoints.each_ref().map(Rc::as_ref);
        let [a, b, c] = frame.form;
        let d = map.degree() as i128;
        let u = maps::sum(model, domain, range, &map, first)?.degree() as i128 - d - a;
        let v = maps::sum(model, domain, range, &map, second)?.degree() as i128 - d - c;
        let mut numerator = [det([u, b], [v, 2 * c])?, det([2 * a, u], [b, v])?];
        // The symmetric degree-polarization Gram is [[2a,b],[b,2c]].
        let mut denominator = -frame.span_discriminant;
        let common = gcd(
            gcd(denominator.unsigned_abs(), numerator[0].unsigned_abs()),
            numerator[1].unsigned_abs(),
        ) as i128;
        denominator /= common;
        numerator.iter_mut().for_each(|k| *k /= common);
        let recovered = frame.combination(model, numerator)?;
        let scaled = scale_map(model, [domain, range], &map, denominator)?;
        if !recovered.compare(model, &scaled)?.is_eq() {
            return Err(Error::plain(
                "CM Hom coordinates failed to reconstruct the exact elliptic morphism",
            ));
        }
        model.allocate(5)?;
        Ok(Self {
            numerator,
            denominator,
            map,
        })
    }

    pub(crate) fn literal(model: &mut Model, frame: &Hom, numerator: [i128; 2]) -> Result<Self> {
        let map = frame.combination(model, numerator)?;
        Ok(Self {
            numerator,
            denominator: 1,
            map: persist(model, map)?,
        })
    }

    pub(crate) fn add(model: &mut Model, frame: &Hom, left: &Self, right: &Self) -> Result<Self> {
        let [source, target] = frame.endpoints.each_ref().map(Rc::as_ref);
        let map = maps::sum(model, source, target, &left.map, &right.map)?;
        let map = persist(model, map)?;
        Self::read(model, frame, map)
    }

    pub(crate) fn negate(model: &mut Model, frame: &Hom, input: &Self) -> Result<Self> {
        let [source, target] = frame.endpoints.each_ref().map(Rc::as_ref);
        let map = maps::negate(model, source, target, &input.map)?;
        let map = persist(model, map)?;
        Self::read(model, frame, map)
    }

    pub(crate) fn compare(&self, other: &Self) -> std::cmp::Ordering {
        (self.numerator, self.denominator).cmp(&(other.numerator, other.denominator))
    }

    pub(crate) fn report(&self) -> String {
        format!(
            "rvCMHomCoordinates = numerator={:?}, denominator={}, in_integral_span={}, degree={}",
            self.numerator,
            self.denominator,
            self.denominator == 1,
            self.map.degree()
        )
    }
}

pub(crate) fn dual(
    model: &mut Model,
    frame: &Hom,
    supplied: Option<[Rc<Morphism>; 2]>,
) -> Result<Hom> {
    let [domain, range] = frame.endpoints.each_ref().map(Rc::as_ref);
    let supplied = if let Some(duals) = supplied {
        duals
    } else {
        let [a, b] = &frame.maps;
        let a = a.dual(model, domain, range)?;
        let b = b.dual(model, domain, range)?;
        [persist(model, a)?, persist(model, b)?]
    };
    for (map, dual) in frame.maps.iter().zip(&supplied) {
        verify_dual(model, map, dual, [domain, range])?;
    }
    let mut child = Hom::new(model, supplied, None, false)?;
    if child.form != frame.form || child.fundamental != frame.fundamental {
        return Err(invalid());
    }
    child.maximal = frame.maximal;
    Ok(child)
}

pub(super) fn verify_dual(
    model: &mut Model,
    map: &Morphism,
    dual: &Morphism,
    [domain, range]: [&Curve; 2],
) -> Result<()> {
    if dual.source != map.target || dual.target != map.source || dual.degree() != map.degree() {
        return Err(Error::plain(
            "CM Hom dual requires reversed named endpoints and the same degree",
        ));
    }
    dual.certify(model, range, domain)?;
    let a = dual.compose(model, domain, map)?;
    let b = map.compose(model, range, dual)?;
    let degree = map.degree() as i128;
    let expected_a = Morphism::multiple(model, map.source.clone(), domain, degree)?;
    let expected_b = Morphism::multiple(model, map.target.clone(), range, degree)?;
    if !a.compare(model, &expected_a)?.is_eq() || !b.compare(model, &expected_b)?.is_eq() {
        return Err(Error::plain(
            "CM Hom dual fails an exact Rosati sandwich identity",
        ));
    }
    Ok(())
}

pub(crate) fn inverse(model: &mut Model, target: &Hom, map: &Morphism) -> Result<Coordinates> {
    if map.degree() != 1 {
        return Err(Error::plain(
            "only a degree-one elliptic morphism has an integral geometric inverse",
        ));
    }
    let [range, domain] = target.endpoints.each_ref().map(Rc::as_ref);
    let first = &target.maps[0];
    if map.source != first.target || map.target != first.source {
        return Err(Error::plain(
            "CM Hom inverse frame must have reversed named endpoints",
        ));
    }
    let coords = map
        .coordinates
        .as_ref()
        .expect("a degree-one map is nonconstant");
    let Some(x_scale) = coords.x.numerator.get(1) else {
        return Err(invalid());
    };
    let Some(y_scale) = coords.y_multiplier.numerator.first() else {
        return Err(invalid());
    };
    let field = domain.field();
    // For a short-Weierstrass isomorphism these are u² and u³.
    let inverse_x = field.inverse_unit(model, x_scale)?;
    let u = field.multiply(model, y_scale, &inverse_x)?;
    let u_inverse = field.inverse_unit(model, &u)?;
    let result = Morphism::chart(
        model,
        map.target.clone(),
        map.source.clone(),
        range,
        domain,
        &u_inverse,
    )?;
    verify_dual(model, map, &result, [domain, range])?;
    let result = persist(model, result)?;
    Coordinates::read(model, target, result)
}

pub(crate) fn compose(
    model: &mut Model,
    target: &Hom,
    outer: &Morphism,
    inner: &Morphism,
) -> Result<Coordinates> {
    endpoints(target, outer, inner)?;
    let result = outer.compose(model, &target.endpoints[0], inner)?;
    let result = persist(model, result)?;
    Coordinates::read(model, target, result)
}

fn endpoints(target: &Hom, outer: &Morphism, inner: &Morphism) -> Result<()> {
    let witness = &target.maps[0];
    if inner.target != outer.source
        || inner.source != witness.source
        || outer.target != witness.target
    {
        return Err(Error::plain(
            "CM Hom composition requires composable named endpoints and the correct target frame",
        ));
    }
    Ok(())
}

type HermiteBasis = [[i128; 2]; 2];
type HermiteWitnesses = [[i128; 4]; 2];

pub(crate) struct Product {
    pub(crate) target_name: Label,
    target: Rc<Hom>,
    columns: [Rc<Coordinates>; 4],
    denominator: i128,
    hnf: HermiteBasis,
    witnesses: HermiteWitnesses,
    smith: [u128; 2],
    intersection_indices: [u128; 2],
    orientation: [i8; 2],
    ideal_product: bool,
}

impl Product {
    pub(crate) fn new(
        model: &mut Model,
        outer: &Hom,
        inner: &Hom,
        target_name: Label,
        target: Rc<Hom>,
    ) -> Result<Self> {
        endpoints(&target, &outer.maps[0], &inner.maps[0])?;
        let [a, b] = &outer.maps;
        let [c, d] = &inner.maps;
        let columns = [
            compose(model, &target, a, c)?,
            compose(model, &target, a, d)?,
            compose(model, &target, b, c)?,
            compose(model, &target, b, d)?,
        ]
        .map(Rc::new);
        let denominator = columns.iter().try_fold(1, |d, x| lcm(d, x.denominator))?;
        let mut integers = [[0; 2]; 4];
        for (row, coordinate) in integers.iter_mut().zip(&columns) {
            for (entry, value) in row.iter_mut().zip(coordinate.numerator) {
                *entry = mul(value, denominator / coordinate.denominator)?;
            }
        }
        let mut minors = [[0; 4]; 4];
        let mut determinant = 0;
        for (i, &v) in integers.iter().enumerate() {
            for (j, &w) in integers.iter().enumerate().skip(i + 1) {
                model.spend(2)?;
                minors[i][j] = det(v, w)?;
                determinant = gcd(determinant, minors[i][j].unsigned_abs());
            }
        }
        let d1 = integers
            .into_iter()
            .flatten()
            .fold(0, |g, x| gcd(g, x.unsigned_abs()));
        if d1 == 0 || determinant == 0 {
            return Err(invalid());
        }
        let smith = [d1, determinant / d1];
        let (hnf, witnesses) = hermite(model, integers)?;
        if (det(hnf[0], hnf[1])? as u128) != determinant || !smith[1].is_multiple_of(smith[0]) {
            return Err(invalid());
        }
        let orientation = [
            orientation_pair(
                [minors[0][2], minors[1][3]],
                inner.maps.each_ref().map(|m| m.degree() as i128),
                denominator,
                outer.span_discriminant,
                target.span_discriminant,
            )?,
            orientation_pair(
                [minors[0][1], minors[2][3]],
                outer.maps.each_ref().map(|m| m.degree() as i128),
                denominator,
                inner.span_discriminant,
                target.span_discriminant,
            )?,
        ];
        let mut intersection_indices = [1u128; 2];
        let den = denominator as u128;
        for value in smith {
            let common = gcd(value, den);
            intersection_indices[0] = intersection_indices[0]
                .checked_mul(value / common)
                .ok_or_else(overflow)?;
            intersection_indices[1] = intersection_indices[1]
                .checked_mul(den / common)
                .ok_or_else(overflow)?;
        }
        let ideal_product = [outer, inner, target.as_ref()]
            .into_iter()
            .all(|h| h.maximal && h.index == 1);
        if ideal_product {
            if denominator != 1
                || determinant != 1
                || outer.fundamental != target.fundamental
                || inner.fundamental != target.fundamental
            {
                return Err(invalid());
            }
            let order = Order::new(target.fundamental)?;
            let [outer_class, inner_class] =
                orient_classes(model, order, [outer, inner], orientation)?;
            if order.compose(model, outer_class, inner_class)? != target.class(model, order)? {
                return Err(Error::plain(
                    "CM Hom product disagrees with oriented ideal class composition",
                ));
            }
        }
        model.allocate(62 + target_name.text.len())?;
        Ok(Self {
            target_name,
            target,
            columns,
            denominator,
            hnf,
            witnesses,
            smith,
            intersection_indices,
            orientation,
            ideal_product,
        })
    }

    pub(crate) fn coordinate(
        &self,
        model: &mut Model,
        index: usize,
        basis: bool,
    ) -> Result<Rc<Coordinates>> {
        if !basis {
            return self
                .columns
                .get(index)
                .cloned()
                .ok_or_else(|| Error::plain("CM Hom product column index is out of range"));
        }
        let witness = self
            .witnesses
            .get(index)
            .ok_or_else(|| Error::plain("CM Hom product basis index is out of range"))?;
        let endpoints = self.target.endpoints.each_ref().map(Rc::as_ref);
        let first = &self.target.maps[0];
        let mut result = Morphism::zero(first.source.clone(), first.target.clone());
        for (&coefficient, column) in witness.iter().zip(&self.columns) {
            if coefficient != 0 {
                let term = scale_map(model, endpoints, &column.map, coefficient)?;
                result = maps::sum(model, endpoints[0], endpoints[1], &result, &term)?;
            }
        }
        let result = persist(model, result)?;
        let coordinates = Coordinates::read(model, &self.target, result)?;
        for (a, b) in coordinates.numerator.into_iter().zip(self.hnf[index]) {
            if mul(a, self.denominator)? != mul(b, coordinates.denominator)? {
                return Err(invalid());
            }
        }
        Ok(Rc::new(coordinates))
    }

    pub(crate) fn report(&self, name: &str) -> String {
        let columns = self
            .columns
            .each_ref()
            .map(|x| (x.numerator, x.denominator));
        let degrees = self.columns.each_ref().map(|x| x.map.degree());
        let [target_index, product_index] = self.intersection_indices;
        let index = if self.denominator == 1 {
            target_index.to_string()
        } else {
            "nonintegral".into()
        };
        format!("rvCMHomProductData[{name}] = rational_rank=2, product_degrees={degrees:?}, coordinates={columns:?}, common_denominator={}, hermite_numerators={:?}, smith_numerators={:?}, integral_index={index}, target_intersection_index={target_index}, product_intersection_index={product_index}, orientation_outer_inner={:?}, full_ideal_product_certified={}",
            self.denominator,self.hnf,self.smith,self.orientation,self.ideal_product)
    }
}

fn orientation_pair(
    determinants: [i128; 2],
    degrees: [i128; 2],
    denominator: i128,
    source_discriminant: i128,
    target_discriminant: i128,
) -> Result<i8> {
    let sign = determinants[0].signum();
    if sign == 0 || determinants[1].signum() != sign {
        return Err(invalid());
    }
    for (determinant, degree) in determinants.into_iter().zip(degrees) {
        let left = Wide::<16>::product(determinant, determinant).multiply(-target_discriminant)?;
        let right = Wide::<16>::product(degree, degree)
            .multiply(-source_discriminant)?
            .multiply(denominator)?
            .multiply(denominator)?
            .multiply(denominator)?
            .multiply(denominator)?;
        if left.magnitude != right.magnitude {
            return Err(invalid());
        }
    }
    Ok(sign as i8)
}

fn orient_classes(
    model: &mut Model,
    order: Order,
    frames: [&Hom; 2],
    signs: [i8; 2],
) -> Result<[crate::hyperelliptic_kummer_surface::cm_ring::Class; 2]> {
    let [a, b] = std::array::from_fn(|i| -> Result<_> {
        let class = frames[i].class(model, order)?;
        if signs[i] < 0 {
            order.inverse(model, class)
        } else {
            Ok(class)
        }
    });
    Ok([a?, b?])
}

pub(super) fn scale_map(
    model: &mut Model,
    [domain, range]: [&Curve; 2],
    map: &Morphism,
    coefficient: i128,
) -> Result<Morphism> {
    if coefficient == 1 {
        return Ok(map.clone());
    }
    if coefficient == 0 || map.degree() == 0 {
        return Ok(Morphism::zero(map.source.clone(), map.target.clone()));
    }
    if coefficient == -1 {
        return maps::negate(model, domain, range, map);
    }
    let multiple = Morphism::multiple(model, map.target.clone(), range, coefficient)?;
    multiple.compose(model, domain, map)
}

fn persist(model: &mut Model, map: Morphism) -> Result<Rc<Morphism>> {
    let Value::ECMap(map) = ec_morphism::persist(model, map)? else {
        unreachable!()
    };
    Ok(map)
}

fn det([a, b]: [i128; 2], [c, d]: [i128; 2]) -> Result<i128> {
    Wide::<16>::product(a, d)
        .add(Wide::product(b, c).negate())?
        .signed()
        .ok_or_else(overflow)
}

fn mul(a: i128, b: i128) -> Result<i128> {
    a.checked_mul(b).ok_or_else(overflow)
}

fn lcm(a: i128, b: i128) -> Result<i128> {
    mul(a / (gcd(a.unsigned_abs(), b.unsigned_abs()) as i128), b)
}

fn hermite(model: &mut Model, columns: [[i128; 2]; 4]) -> Result<(HermiteBasis, HermiteWitnesses)> {
    let ordinates = columns.map(|x| x[1]);
    let (g, mut vertical) = bezout_four(model, ordinates)?;
    if g == 0 {
        return Err(invalid());
    }
    let abscissas = columns.map(|x| x[0]);
    let head = dot(abscissas, vertical)?;
    let mut horizontal_rows = [0; 4];
    for (i, entry) in horizontal_rows.iter_mut().enumerate() {
        *entry = abscissas[i]
            .checked_sub(mul(ordinates[i] / g, head)?)
            .ok_or_else(overflow)?;
    }
    let (h, mut horizontal) = bezout_four(model, horizontal_rows)?;
    if h == 0 {
        return Err(invalid());
    }
    let vertical_coefficient = dot(ordinates.map(|y| y / g), horizontal)?;
    let quotient = head.div_euclid(h);
    for (horizontal, vertical) in horizontal.iter_mut().zip(&mut vertical) {
        *horizontal = horizontal
            .checked_sub(mul(vertical_coefficient, *vertical)?)
            .ok_or_else(overflow)?;
        *vertical = vertical
            .checked_sub(mul(quotient, *horizontal)?)
            .ok_or_else(overflow)?;
    }
    let remainder = head.rem_euclid(h);
    let hnf = [[h, 0], [remainder, g]];
    for (basis, coefficients) in hnf.into_iter().zip([horizontal, vertical]) {
        if basis != [dot(abscissas, coefficients)?, dot(ordinates, coefficients)?] {
            return Err(invalid());
        }
    }
    Ok((hnf, [horizontal, vertical]))
}

fn dot(left: [i128; 4], right: [i128; 4]) -> Result<i128> {
    left.into_iter()
        .zip(right)
        .try_fold(Wide::<16>::from_i128(0), |sum, (a, b)| {
            sum.add(Wide::product(a, b))
        })?
        .signed()
        .ok_or_else(overflow)
}

fn bezout_four(model: &mut Model, values: [i128; 4]) -> Result<(i128, [i128; 4])> {
    let mut gcd = 0;
    let mut coefficients = [0i128; 4];
    for (i, value) in values.into_iter().enumerate() {
        let (mut a, mut b) = (gcd, value.checked_abs().ok_or_else(overflow)?);
        let (mut u, mut v) = ([1i128, 0], [0i128, value.signum()]);
        while b != 0 {
            model.spend(1)?;
            let quotient = a / b;
            (a, b) = (b, a % b);
            for (x, y) in u.iter_mut().zip(&mut v) {
                (*x, *y) = (*y, x.checked_sub(mul(quotient, *y)?).ok_or_else(overflow)?);
            }
        }
        for c in &mut coefficients[..i] {
            *c = mul(*c, u[0])?;
        }
        coefficients[i] = u[1];
        gcd = a;
    }
    Ok((gcd, coefficients))
}

fn invalid() -> Error {
    Error::plain("CM Hom exact composition, norm orientation or lattice certificate disagrees")
}

fn overflow() -> Error {
    Error::plain("CM Hom composition arithmetic exceeds the signed integer range")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hermite_retains_noncoordinate_generators_and_integral_witnesses() {
        let columns = [[6, 0], [3, 4], [0, -8], [-9, -12]];
        let (basis, witnesses) = hermite(&mut Model::new(), columns)
            .ok()
            .expect("rank-two column lattice");
        assert_eq!(basis, [[6, 0], [3, 4]]);
        assert_eq!(det(basis[0], basis[1]).ok(), Some(24));
        for (vector, coefficients) in basis.into_iter().zip(witnesses) {
            assert_eq!(
                dot(columns.map(|v| v[0]), coefficients).ok(),
                Some(vector[0])
            );
            assert_eq!(
                dot(columns.map(|v| v[1]), coefficients).ok(),
                Some(vector[1])
            );
        }
    }
}
