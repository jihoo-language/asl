use std::cmp::Ordering;
use std::fmt::Write;
use std::rc::Rc;

use super::super::{raw_coordinates, K3Vector, COLUMN_LIMBS};
use super::{add_scaled, column_integer, Class, ColumnInteger, Lattice, Vector, AMBIENT_RANK};
use crate::error::{Error, Result};
use crate::model::{Model, MAX_SHOWN_BYTES};
use crate::rv_lens::gcd;
use crate::wide_integer::Natural;

type Magnitude = Natural<COLUMN_LIMBS>;

#[derive(Clone, Copy, PartialEq, Eq)]
enum WallKind {
    Divisorial,
    Small5,
    Small3,
    Regular,
}

impl WallKind {
    fn classify(square: Option<i128>, divisibility: i128) -> Self {
        match (square, divisibility) {
            (Some(-2), 1) | (Some(-6), 3) => Self::Divisorial,
            (Some(-4), 1) => Self::Small5,
            (Some(-24), 3) => Self::Small3,
            _ => Self::Regular,
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Divisorial => "divisorial",
            Self::Small5 => "small codimension 5",
            Self::Small3 => "small codimension 3",
            Self::Regular => "regular",
        }
    }
}

pub(crate) fn report(model: &mut Model, lattice: &Lattice, class: &Class) -> Result<String> {
    let (primitive, content) = primitive(model, class)?;
    let square = primitive.wide_pair(model, lattice, &primitive)?;
    let divisibility = primitive.divisibility(model, lattice)?;
    let kind = WallKind::classify(square.signed(), divisibility);
    let sign = if square.negative && !square.is_zero() {
        "-"
    } else {
        ""
    };
    Ok(format!(
        "rvOG10Wall = content={content}, square={sign}{}, div={divisibility}, kind={}",
        decimal(model, square.magnitude)?,
        kind.label()
    ))
}

pub(crate) fn reflect(model: &mut Model, lattice: &Lattice, d: &Class, x: &Class) -> Result<Class> {
    let (d, _) = primitive(model, d)?;
    let square = d.wide_pair(model, lattice, &d)?.signed();
    let divisibility = d.divisibility(model, lattice)?;
    if WallKind::classify(square, divisibility) != WallKind::Divisorial {
        return Err(Error::plain("rvOG10Reflect requires a divisorial wall ray"));
    }
    let pairing = x.wide_pair(model, lattice, &d)?;
    let shift = pairing
        .divide_exact(model, column_integer(divisibility))?
        .ok_or_else(|| Error::plain("OG10 wall reflection is not integral"))?;
    linear_combination(model, x, column_integer(1), &d, shift, column_integer(1))?
        .ok_or_else(|| Error::plain("OG10 wall reflection is not integral"))
}

pub(crate) struct Segment {
    crossings: Vec<Crossing>,
    scales: [u128; 2],
    reversed: bool,
}

struct Crossing {
    u: usize,
    v: Magnitude,
    square: i128,
    divisibility: i128,
    kind: WallKind,
    class: Class,
}

impl Segment {
    pub(crate) fn new(
        model: &mut Model,
        lattice: &Lattice,
        h0: &Class,
        h1: &Class,
    ) -> Result<Self> {
        let (mut h0, scale0) = primitive(model, h0)?;
        let (mut h1, scale1) = primitive(model, h1)?;
        let mut a = h0.wide_pair(model, lattice, &h0)?;
        let b = h0.wide_pair(model, lattice, &h1)?;
        let mut c = h1.wide_pair(model, lattice, &h1)?;
        if [a, b, c].iter().any(|x| x.negative || x.is_zero()) {
            return Err(Error::plain(
                "OG10 wall search requires positive endpoints in the same cone component",
            ));
        }
        let delta = b.multiply_wide(b)?.add(a.multiply_wide(c)?.negate())?;
        if delta.negative || delta.is_zero() {
            return Err(Error::plain("OG10 wall search requires independent endpoints spanning a plane of signature (1,1)"));
        }

        // A larger C gives a shorter scan. Restore the input orientation
        // when reconstructing walls and reporting crossing parameters.
        let reversed = a.magnitude.cmp(c.magnitude).is_gt();
        let mut scales = [scale0, scale1];
        if reversed {
            std::mem::swap(&mut h0, &mut h1);
            std::mem::swap(&mut a, &mut c);
            scales.swap(0, 1);
        }
        let orientation = if reversed { -1 } else { 1 };
        let mut crossings = Vec::new();
        for s in [2i128, 4, 6, 24] {
            let (bound, _) = delta.multiply(s)?.magnitude.div_rem(model, c.magnitude)?;
            let max_u = bound
                .integer_sqrt(model)?
                .small()
                .and_then(|u| usize::try_from(u).ok())
                .ok_or_else(work_limit)?;
            let work = max_u
                .checked_add(1)
                .and_then(|n| n.checked_mul(8))
                .ok_or_else(work_limit)?;
            model.spend(work).map_err(|_| work_limit())?;

            let as_term = a.multiply(s)?;
            for u in 0..=max_u {
                let u_int = u as i128;
                let radicand =
                    delta.multiply_wide(ColumnInteger::product(u_int, u_int).add(as_term)?)?;
                let root = radicand.magnitude.integer_sqrt(model)?;
                if root.checked_product(root) != Some(radicand.magnitude) {
                    continue;
                }
                let bu = b.multiply(u_int)?;
                let numerator = ColumnInteger {
                    magnitude: root,
                    negative: false,
                }
                .add(bu.negate())?;
                let Some(v) = numerator.divide_exact(model, a)? else {
                    continue;
                };
                if v.negative && !v.is_zero() {
                    continue;
                }
                let first = c.multiply(u_int)?.add(b.multiply_wide(v)?)?.negate();
                let second = bu.add(a.multiply_wide(v)?)?;
                let Some(class) = linear_combination(
                    model,
                    &h0,
                    first.multiply(orientation)?,
                    &h1,
                    second.multiply(orientation)?,
                    delta,
                )?
                else {
                    continue;
                };
                if primitive_content(model, &class)? != 1 {
                    continue;
                }
                let divisibility = class.divisibility(model, lattice)?;
                let kind = WallKind::classify(Some(-s), divisibility);
                if kind != WallKind::Regular {
                    crate::hw_algebra::storage(
                        crossings.len() + 1,
                        2 * AMBIENT_RANK + COLUMN_LIMBS,
                    )?;
                    crossings.push(Crossing {
                        u,
                        v: v.magnitude,
                        square: -s,
                        divisibility,
                        kind,
                        class,
                    });
                }
            }
        }
        model.spend(crossings.len() * (1 + crossings.len().max(1).ilog2() as usize))?;
        crossings.sort_by(compare_crossings);
        if reversed {
            crossings.reverse();
        }
        Ok(Self {
            crossings,
            scales,
            reversed,
        })
    }

    pub(crate) fn wall_at(self, index: usize) -> Result<Class> {
        self.crossings
            .into_iter()
            .nth(index)
            .map(|crossing| crossing.class)
            .ok_or_else(|| Error::plain(format!("rvOG10WallAt index {index} is out of range")))
    }

    pub(crate) fn report(&self, model: &mut Model) -> Result<String> {
        let mut out = format!("rvOG10Walls = count={}, crossings=[", self.crossings.len());
        for (index, wall) in self.crossings.iter().enumerate() {
            if index != 0 {
                out.push_str(", ");
            }
            let t = wall.parameter(model, self.scales, self.reversed)?;
            write!(
                out,
                "{{t={t}, square={}, div={}, kind={}}}",
                wall.square,
                wall.divisibility,
                wall.kind.label()
            )
            .expect("writing to String");
            if out.len() >= MAX_SHOWN_BYTES {
                return Err(Error::plain("OG10 wall report exceeds the output limit"));
            }
        }
        out.push(']');
        Ok(out)
    }
}

impl Crossing {
    fn parameter(
        &self,
        model: &mut Model,
        [scale0, scale1]: [u128; 2],
        reversed: bool,
    ) -> Result<String> {
        if self.u == 0 {
            return Ok(if reversed { "1" } else { "0" }.into());
        }
        if self.v.is_zero() {
            return Ok(if reversed { "0" } else { "1" }.into());
        }
        // Restore the original endpoint lengths and undo the optional swap.
        let mut u = Magnitude::from_small(self.u as u128)
            .multiply_u128(scale0)
            .ok_or_else(work_limit)?;
        let mut v = self.v.multiply_u128(scale1).ok_or_else(work_limit)?;
        if reversed {
            std::mem::swap(&mut u, &mut v);
        }
        let gcd = u.gcd(model, v)?;
        let (numerator, _) = u.div_rem(model, gcd)?;
        let (mut denominator, _) = v.div_rem(model, gcd)?;
        denominator
            .add_multiple(numerator, 1)
            .ok_or_else(work_limit)?;
        Ok(format!(
            "{}/{}",
            decimal(model, numerator)?,
            decimal(model, denominator)?
        ))
    }
}

fn compare_crossings(left: &Crossing, right: &Crossing) -> Ordering {
    // Positive rescaling of endpoints preserves crossing order. Use the
    // primitive parameters and extra limbs for multiplying by a usize.
    let cross = |u: usize, v: Magnitude| {
        v.resize::<{ COLUMN_LIMBS + 2 }>()
            .and_then(|v| v.multiply_u128(u as u128))
            .expect("two extra limbs hold a usize multiplication")
    };
    cross(left.u, right.v).cmp(cross(right.u, left.v))
}

fn primitive(model: &mut Model, class: &Class) -> Result<(Class, u128)> {
    let divisor = primitive_content(model, class)?;
    let scale = ColumnInteger {
        magnitude: Magnitude::from_small(divisor),
        negative: false,
    };
    let reduced = linear_combination(
        model,
        class,
        column_integer(1),
        class,
        column_integer(0),
        scale,
    )?
    .expect("primitive content divides the class");
    Ok((reduced, divisor))
}

fn primitive_content(model: &mut Model, class: &Class) -> Result<u128> {
    model.spend(AMBIENT_RANK + 1)?;
    let coordinates = raw_coordinates(&class.perpendicular);
    let content = coordinates
        .iter()
        .fold(class.coefficient.unsigned_abs(), |content, entry| {
            gcd(content, entry.unsigned_abs())
        });
    if content == 0 {
        return Err(Error::plain("OG10 wall rays require a nonzero class"));
    }
    // Integer divisibility must also preserve the marked affine gluing
    // condition. Its only possible additional obstruction is parity.
    let reduced = coordinates.map(|x| signed_divide(x, content));
    Ok(if vector_from_raw(reduced).is_some() {
        content
    } else {
        content / 2
    })
}

fn signed_divide(value: i128, divisor: u128) -> i128 {
    if divisor > i128::MAX as u128 {
        if value == i128::MIN {
            -1
        } else {
            0
        }
    } else {
        value / divisor as i128
    }
}

fn linear_combination(
    model: &mut Model,
    left: &Class,
    left_scale: ColumnInteger,
    right: &Class,
    right_scale: ColumnInteger,
    denominator: ColumnInteger,
) -> Result<Option<Class>> {
    let mut coordinates = [column_integer(0); AMBIENT_RANK];
    add_scaled(model, &mut coordinates, &left.perpendicular, left_scale)?;
    add_scaled(model, &mut coordinates, &right.perpendicular, right_scale)?;
    let coefficient = left_scale
        .multiply(left.coefficient)?
        .add(right_scale.multiply(right.coefficient)?)?;
    let mut quotient = [column_integer(0); AMBIENT_RANK + 1];
    for (out, value) in quotient
        .iter_mut()
        .zip(coordinates.into_iter().chain([coefficient]))
    {
        let Some(value) = value.divide_exact(model, denominator)? else {
            return Ok(None);
        };
        *out = value;
    }
    let parity = std::array::from_fn(|i| i128::from(quotient[i].magnitude.residue(2)));
    if vector_from_raw(parity).is_none() {
        return Ok(None);
    }
    let mut reduced = [0; AMBIENT_RANK + 1];
    for (out, value) in reduced.iter_mut().zip(quotient) {
        *out = value.signed().ok_or_else(super::coordinate_range)?;
    }
    let raw = reduced[..AMBIENT_RANK]
        .try_into()
        .expect("Mukai coordinate count");
    Ok(vector_from_raw(raw).map(|perpendicular| Class {
        perpendicular: Rc::new(perpendicular),
        coefficient: reduced[AMBIENT_RANK],
    }))
}

fn vector_from_raw(coordinates: [i128; AMBIENT_RANK]) -> Option<Vector> {
    let (wedge, exceptional) = coordinates[2..].split_at(6);
    Some(Vector {
        r: coordinates[0],
        s: coordinates[1],
        c: Rc::new(
            K3Vector::new(
                wedge.try_into().expect("six wedge coordinates"),
                exceptional
                    .try_into()
                    .expect("sixteen exceptional coordinates"),
            )
            .ok()?,
        ),
    })
}

fn decimal(model: &mut Model, mut magnitude: Magnitude) -> Result<String> {
    let mut tail = Vec::new();
    while magnitude.small().is_none() {
        model.spend(COLUMN_LIMBS)?;
        let (quotient, remainder) = magnitude.divide_u32(1_000_000_000);
        magnitude = quotient;
        tail.push(remainder);
    }
    let mut out = magnitude.small().expect("final decimal digits").to_string();
    for digit in tail.into_iter().rev() {
        write!(out, "{digit:09}").expect("writing to String");
    }
    Ok(out)
}

fn work_limit() -> Error {
    Error::plain("OG10 wall segment search exceeds the work limit")
}
