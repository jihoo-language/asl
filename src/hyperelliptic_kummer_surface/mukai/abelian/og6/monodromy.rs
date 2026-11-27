use std::rc::Rc;

use super::{abelian_range, coordinate_range, Class, Integer, Lattice, Vector, RANK};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::rational_sparse::decimal;
use crate::model::Model;
use crate::rv_lens::gcd;
use crate::wide_integer::Natural;

#[derive(Clone, Copy, PartialEq, Eq)]
struct Residue([u8; 2]);

impl Residue {
    fn for_primitive(lattice: &Lattice, class: &Class, divisibility: i128) -> Result<Self> {
        if divisibility < 2 {
            return Ok(Self([0, 0]));
        }
        let defining = &lattice.vector.coordinates;
        let first_odd = defining
            .iter()
            .position(|entry| entry & 1 != 0)
            .ok_or_else(residue_error)?;
        let b = class.perpendicular.coordinates[first_odd] & 1;
        let k = class.exceptional & 1;
        // In U^4/2, the annihilator of w^perp is precisely the span
        // of w. The corresponding glue generator is (2z−w)/2,
        // where (z,w)=1, so its reduction mod 2 is read from w.
        if divisibility != 2
            || b == 0 && k == 0
            || !class
                .perpendicular
                .coordinates
                .iter()
                .zip(defining)
                .all(|(&entry, &w)| (entry & 1) == b * (w & 1))
        {
            return Err(residue_error());
        }
        Ok(Self([b as u8, k as u8]))
    }

    fn full_invariant(self) -> u8 {
        self.0[0] + self.0[1]
    }

    fn quadratic_value(self) -> &'static str {
        match self.full_invariant() {
            0 => "0",
            1 => "-1/2",
            _ => "-1",
        }
    }
}

struct Orbit {
    primitive: Class,
    content: u128,
    square: Integer,
    divisibility: i128,
    residue: Residue,
}

impl Orbit {
    fn new(model: &mut Model, lattice: &Lattice, class: &Class) -> Result<Self> {
        model.spend(RANK + 1)?;
        // The perpendicular of an integral vector is saturated in U^4.
        let content = class
            .perpendicular
            .coordinates
            .iter()
            .fold(class.exceptional.unsigned_abs(), |content, &entry| {
                gcd(content, entry.unsigned_abs())
            });
        let divide = |entry: i128| {
            Natural::<13>::from_small(entry.unsigned_abs() / content.max(1))
                .signed(entry < 0)
                .expect("dividing a signed coordinate by a positive integer")
        };
        let primitive = Class {
            perpendicular: Rc::new(Vector {
                coordinates: class.perpendicular.coordinates.map(divide),
            }),
            exceptional: divide(class.exceptional),
        };
        let square = primitive.wide_pair(model, &primitive)?;
        let divisibility = primitive.divisibility(model, lattice)?;
        let residue = Residue::for_primitive(lattice, &primitive, divisibility)?;
        Ok(Self {
            primitive,
            content,
            square,
            divisibility,
            residue,
        })
    }

    fn same_common_invariants(&self, other: &Self) -> bool {
        self.content == other.content
            && self.square.magnitude == other.square.magnitude
            && (self.square.is_zero() || self.square.negative == other.square.negative)
            && self.divisibility == other.divisibility
    }

    fn wall_kind(&self) -> WallKind {
        match (self.square.signed(), self.divisibility) {
            (Some(-2 | -4), 2) => WallKind::Divisorial,
            (Some(-2), 1) => WallKind::Small,
            _ => WallKind::Regular,
        }
    }

    fn report(&self) -> String {
        let [b, k] = self.residue.0;
        let sign = if self.square.negative && !self.square.is_zero() {
            "-"
        } else {
            ""
        };
        format!(
            "content={}, square={sign}{}, div={}, residue=[{b},{k}], discq={}",
            self.content,
            decimal(self.square.magnitude),
            self.divisibility,
            self.residue.quadratic_value()
        )
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum WallKind {
    Divisorial,
    Small,
    Regular,
}

impl WallKind {
    fn label(self) -> &'static str {
        match self {
            Self::Divisorial => "divisorial",
            Self::Small => "small",
            Self::Regular => "regular",
        }
    }
}

pub(crate) fn orbit_report(model: &mut Model, lattice: &Lattice, class: &Class) -> Result<String> {
    let orbit = Orbit::new(model, lattice, class)?;
    Ok(format!("rvOG6Orbit = {}", orbit.report()))
}

pub(crate) fn compare_orbits(
    model: &mut Model,
    lattice: &Lattice,
    left: &Class,
    right: &Class,
) -> Result<String> {
    let left = Orbit::new(model, lattice, left)?;
    let right = Orbit::new(model, lattice, right)?;
    let common = left.same_common_invariants(&right);
    let stable = common && left.residue == right.residue;
    let full = common && left.residue.full_invariant() == right.residue.full_invariant();
    Ok(format!("rvOG6SameOrbit = stable={stable}, full={full}"))
}

pub(crate) fn wall_report(model: &mut Model, lattice: &Lattice, class: &Class) -> Result<String> {
    let orbit = Orbit::new(model, lattice, class)?;
    if orbit.content == 0 {
        return Err(Error::plain("rvOG6Wall requires a nonzero ray"));
    }
    Ok(format!(
        "rvOG6Wall = {}, kind={}",
        orbit.report(),
        orbit.wall_kind().label()
    ))
}

pub(crate) fn reflect(
    model: &mut Model,
    lattice: &Lattice,
    ray: &Class,
    class: &Class,
) -> Result<Class> {
    let orbit = Orbit::new(model, lattice, ray)?;
    if orbit.wall_kind() != WallKind::Divisorial {
        return Err(Error::plain(
            "rvOG6Reflect requires a divisorial OG6 wall ray",
        ));
    }
    let d = orbit.primitive;
    let divisor = if orbit.square.signed() == Some(-4) {
        2
    } else {
        1
    };
    let shift = class
        .wide_pair(model, &d)?
        .divide_exact(model, Integer::from_i128(divisor))?
        .ok_or_else(|| Error::plain("OG6 reflection does not preserve the integral lattice"))?;
    let combine = |existing, axis| shift.multiply(axis)?.add(Integer::from_i128(existing));
    let mut perpendicular = Vector::zero();
    for ((out, &existing), &axis) in perpendicular
        .coordinates
        .iter_mut()
        .zip(&class.perpendicular.coordinates)
        .zip(&d.perpendicular.coordinates)
    {
        *out = combine(existing, axis)?
            .signed()
            .ok_or_else(abelian_range)?;
    }
    let exceptional = combine(class.exceptional, d.exceptional)?
        .signed()
        .ok_or_else(coordinate_range)?;
    Ok(Class {
        perpendicular: Rc::new(perpendicular),
        exceptional,
    })
}

fn residue_error() -> Error {
    Error::plain("primitive OG6 vector has an inconsistent discriminant residue")
}
