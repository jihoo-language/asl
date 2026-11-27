use std::cmp::Ordering;
use std::rc::Rc;

use super::super::moduli::{
    binomial_with_square, column_integer, eliminate, ColumnInteger, COLUMN_CELLS,
};
use super::{coordinate_range as abelian_range, Integer, Vector, RANK};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::polarized;
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod monodromy;

const N: usize = 3;
type Column = [ColumnInteger; RANK];

pub(crate) enum Expression {
    Build(Label, i128, Box<ParsedExpr>),
    Exceptional(Label, bool),
    Basis(Label, usize),
    Reflect(Box<[ParsedExpr; 2]>),
}

pub(crate) enum Command {
    Define(Label, Box<ParsedExpr>),
    Lattice(Label),
    Pair(Box<[ParsedExpr; 2]>),
    Divisibility(Box<ParsedExpr>),
    Fujiki(usize, Box<[ParsedExpr; 2]>),
    RiemannRoch(Box<ParsedExpr>),
    Orbit(Box<ParsedExpr>),
    SameOrbit(Box<[ParsedExpr; 2]>),
    Wall(Box<ParsedExpr>),
}

pub(crate) struct Lattice {
    pub(crate) marking: Rc<str>,
    vector: Rc<Vector>,
    columns: Vec<Column>,
}

impl Lattice {
    pub(crate) fn new(model: &mut Model, marking: Rc<str>, vector: Rc<Vector>) -> Result<Self> {
        if vector.wide_pair(model, &vector)?.signed() != Some(2) {
            return Err(Error::plain(
                "rvog6 requires a primitive abelian Mukai vector of square two",
            ));
        }
        let row = ambient_pairings(model, &vector)?.to_vec();
        let mut columns = vec![[column_integer(0); RANK]; RANK];
        for (index, column) in columns.iter_mut().enumerate() {
            column[index] = column_integer(1);
        }
        // Square two already forces primitivity in an even lattice.
        eliminate(
            model,
            &mut columns,
            row,
            "rvog6 requires a primitive abelian Mukai vector",
        )?;
        Ok(Self {
            marking,
            vector,
            columns,
        })
    }

    pub(crate) fn cells(&self) -> usize {
        2 + self.columns.len() * RANK * COLUMN_CELLS
    }

    fn require_perpendicular(&self, model: &mut Model, vector: &Vector) -> Result<()> {
        if !vector.wide_pair(model, &self.vector)?.is_zero() {
            return Err(Error::plain(
                "rvOG6 requires an abelian Mukai vector perpendicular to the defining vector",
            ));
        }
        Ok(())
    }
}

fn ambient_pairings(model: &mut Model, vector: &Vector) -> Result<Column> {
    let mut row = [column_integer(0); RANK];
    for (i, entry) in row.iter_mut().enumerate() {
        *entry = vector.wide_pair(model, &Vector::basis(i)?)?.resize()?;
    }
    Ok(row)
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Class {
    perpendicular: Rc<Vector>,
    exceptional: i128,
}

impl Class {
    pub(crate) fn new(
        model: &mut Model,
        lattice: &Lattice,
        perpendicular: Rc<Vector>,
        exceptional: i128,
    ) -> Result<Self> {
        lattice.require_perpendicular(model, &perpendicular)?;
        Ok(Self {
            perpendicular,
            exceptional,
        })
    }

    pub(crate) fn zero() -> Self {
        Self {
            perpendicular: Rc::new(Vector::zero()),
            exceptional: 0,
        }
    }

    pub(crate) fn exceptional(half: bool) -> Self {
        Self {
            exceptional: if half { 1 } else { 2 },
            ..Self::zero()
        }
    }

    pub(crate) fn basis(model: &mut Model, lattice: &Lattice, index: usize) -> Result<Self> {
        match index {
            7 => Ok(Self::exceptional(true)),
            0..=6 => {
                model.spend(RANK * COLUMN_CELLS)?;
                let column = &lattice.columns[index];
                let mut perpendicular = Vector::zero();
                for (out, coordinate) in perpendicular.coordinates.iter_mut().zip(column) {
                    *out = coordinate.signed().ok_or_else(abelian_range)?;
                }
                Ok(Self {
                    perpendicular: Rc::new(perpendicular),
                    exceptional: 0,
                })
            }
            _ => Err(Error::plain("rvOG6Basis index must be between 0 and 7")),
        }
    }

    pub(crate) fn coordinates(&self) -> (&Vector, i128) {
        (&self.perpendicular, self.exceptional)
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(2 * (RANK + 1) + 1)?;
        Ok(Value::OG6(Rc::new(self)))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        model.spend(RANK + 1)?;
        Ok(self.cmp(other))
    }

    pub(crate) fn add(&self, model: &mut Model, other: &Self) -> Result<Self> {
        let exceptional = self
            .exceptional
            .checked_add(other.exceptional)
            .ok_or_else(coordinate_range)?;
        let perpendicular = Rc::new(self.perpendicular.add(model, &other.perpendicular)?);
        Ok(Self {
            perpendicular,
            exceptional,
        })
    }

    pub(crate) fn negate(&self, model: &mut Model) -> Result<Self> {
        let exceptional = self
            .exceptional
            .checked_neg()
            .ok_or_else(coordinate_range)?;
        let perpendicular = Rc::new(self.perpendicular.negate(model)?);
        Ok(Self {
            perpendicular,
            exceptional,
        })
    }

    fn wide_pair(&self, model: &mut Model, other: &Self) -> Result<Integer> {
        let correction = Integer::product(self.exceptional, other.exceptional).multiply(-2)?;
        self.perpendicular
            .wide_pair(model, &other.perpendicular)?
            .add(correction)
    }

    pub(crate) fn pair(&self, model: &mut Model, other: &Self) -> Result<i128> {
        self.wide_pair(model, other)?
            .signed()
            .ok_or_else(|| Error::plain("OG6 pairing exceeds the signed 128-bit range"))
    }

    pub(crate) fn divisibility(&self, model: &mut Model, lattice: &Lattice) -> Result<i128> {
        let row = ambient_pairings(model, &self.perpendicular)?;
        let mut gcd = column_integer(self.exceptional).multiply(-2)?.magnitude;
        for column in &lattice.columns {
            if gcd.small() == Some(1) {
                return Ok(1);
            }
            model.spend(RANK * COLUMN_CELLS)?;
            let mut pairing = column_integer(0);
            for (&coordinate, &entry) in column.iter().zip(&row) {
                if !coordinate.is_zero() && !entry.is_zero() {
                    pairing = pairing.add(coordinate.multiply_wide(entry)?)?;
                }
            }
            gcd = gcd.gcd(model, pairing.magnitude)?;
        }
        gcd.signed(false)
            .ok_or_else(|| Error::plain("OG6 divisibility exceeds the signed 128-bit range"))
    }

    pub(crate) fn fujiki(&self, model: &mut Model, other: &Self, j: usize) -> Result<i128> {
        let left = polarized::fujiki_left_exponent(N, j)?;
        let xx = self.wide_pair(model, self)?;
        let xy = self.wide_pair(model, other)?;
        let yy = other.wide_pair(model, other)?;
        polarized::fujiki_pairings(model, left, j, xx, xy, yy)?
            .checked_mul(4)
            .ok_or_else(|| Error::plain("OG6 Fujiki intersection exceeds the signed 128-bit range"))
    }

    pub(crate) fn riemann_roch(&self, model: &mut Model) -> Result<i128> {
        let q = self.wide_pair(model, self)?;
        binomial_with_square(model, q, N as i128, 3)?
            .checked_mul(4)
            .ok_or_else(|| Error::plain("OG6 Riemann–Roch value exceeds the signed 128-bit range"))
    }
}

pub(crate) fn report(name: &str) -> String {
    format!(
        "rvOG6Lattice[{name}] = n=3, dimension=6, rank=8, signature=(3,5), determinant=-4, discriminant=(Z/2)^2, q=diag(-1/2,-1/2) mod 2"
    )
}

fn coordinate_range() -> Error {
    Error::plain("OG6 half-exceptional coefficient exceeds the signed 128-bit range")
}
