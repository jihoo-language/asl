use std::cmp::Ordering;
use std::rc::Rc;

use super::{
    add_scaled, ambient_pairings, bezout_pairings, column_integer, column_pair, realize,
    riemann_roch_polynomial, ColumnInteger, Lattice, PairingRow, Vector, AMBIENT_RANK,
};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::polarized;
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod walls;

const N: usize = 5;

#[derive(Clone, Copy)]
pub(crate) enum WallOperation {
    Reflect,
    At(usize),
}

impl WallOperation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Reflect => "rvOG10Reflect",
            Self::At(_) => "rvOG10WallAt",
        }
    }
}

pub(crate) enum Expression {
    Build(Label, i128, Box<ParsedExpr>),
    Exceptional(Label),
    Basis(Label, usize),
    Reflect(Box<[ParsedExpr; 2]>),
    WallAt(usize, Box<[ParsedExpr; 2]>),
}

pub(crate) enum Command {
    Lattice(Label),
    Pair(Box<[ParsedExpr; 2]>),
    Divisibility(Box<ParsedExpr>),
    Fujiki(usize, Box<[ParsedExpr; 2]>),
    RiemannRoch(Box<ParsedExpr>),
    Wall(Box<ParsedExpr>),
    Walls(Box<[ParsedExpr; 2]>),
}

pub(super) struct Form {
    bezout_pairings: PairingRow,
    generator_square: ColumnInteger,
}

impl Lattice {
    pub(crate) fn require_og10(&self) -> Result<()> {
        if self.square != 2 {
            return Err(Error::plain(
                "RVOG10 requires a Mukai moduli lattice of square two",
            ));
        }
        Ok(())
    }

    fn og10_form(&self, model: &mut Model) -> Result<&Form> {
        if let Some(form) = self.og10.get() {
            return Ok(form);
        }
        let bezout_pairings = bezout_pairings(model, &self.ambient, &self.bezout)?;
        let bezout_square = column_pair(model, &self.bezout, &bezout_pairings)?;
        // g = 2b − w, so g² = 4b² − 2 and G² = b² − 2.
        let generator_square = bezout_square.add(column_integer(-2))?;
        model.allocate((AMBIENT_RANK + 1) * super::COLUMN_CELLS)?;
        Ok(self.og10.get_or_init(|| {
            Box::new(Form {
                bezout_pairings,
                generator_square,
            })
        }))
    }
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Class {
    perpendicular: Rc<Vector>,
    coefficient: i128,
}

impl Class {
    pub(crate) fn new(
        model: &mut Model,
        lattice: &Lattice,
        perpendicular: Rc<Vector>,
        coefficient: i128,
    ) -> Result<Self> {
        lattice.require_og10()?;
        lattice.require_perpendicular(model, &perpendicular, "rvOG10")?;
        Ok(Self {
            perpendicular,
            coefficient,
        })
    }

    pub(crate) fn zero() -> Self {
        Self {
            perpendicular: Rc::new(Vector::zero()),
            coefficient: 0,
        }
    }

    pub(crate) fn exceptional(model: &mut Model, lattice: &Lattice) -> Result<Self> {
        lattice.require_og10()?;
        let mut coordinates = lattice.raw_column(model, &lattice.bezout, -2)?;
        add_scaled(model, &mut coordinates, &lattice.vector, column_integer(1))?;
        Ok(Self {
            perpendicular: Rc::new(realize(coordinates)?),
            coefficient: 2,
        })
    }

    pub(crate) fn basis(model: &mut Model, lattice: &Lattice, index: usize) -> Result<Self> {
        lattice.require_og10()?;
        match index {
            0..=22 => Ok(Self {
                perpendicular: Rc::new(lattice.basis(model, index)?),
                coefficient: 0,
            }),
            23 => Ok(Self {
                coefficient: 1,
                ..Self::zero()
            }),
            _ => Err(Error::plain("rvOG10Basis index must be between 0 and 23")),
        }
    }

    pub(crate) fn coordinates(&self) -> (&Vector, i128) {
        (&self.perpendicular, self.coefficient)
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        model.spend(AMBIENT_RANK + 1)?;
        Ok(self.cmp(rhs))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(2 * (AMBIENT_RANK + 1) + 1)?;
        Ok(Value::OG10(Rc::new(self)))
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        let coefficient = self
            .coefficient
            .checked_add(rhs.coefficient)
            .ok_or_else(coordinate_range)?;
        let perpendicular = Rc::new(self.perpendicular.add(model, &rhs.perpendicular)?);
        Ok(Self {
            perpendicular,
            coefficient,
        })
    }

    pub(crate) fn negate(&self, model: &mut Model) -> Result<Self> {
        let coefficient = self
            .coefficient
            .checked_neg()
            .ok_or_else(coordinate_range)?;
        let perpendicular = Rc::new(self.perpendicular.negate(model)?);
        Ok(Self {
            perpendicular,
            coefficient,
        })
    }

    fn wide_pair(&self, model: &mut Model, lattice: &Lattice, rhs: &Self) -> Result<ColumnInteger> {
        let mut pairing: ColumnInteger = self
            .perpendicular
            .wide_pair(model, &rhs.perpendicular)?
            .resize()?;
        // On W, (g, u)/2 = (b, u). Use the integral pairings directly.
        if self.coefficient != 0 {
            let row = ambient_pairings(model, &lattice.ambient, &rhs.perpendicular)?;
            let cross = column_pair(model, &lattice.bezout, &row)?.multiply(self.coefficient)?;
            pairing = pairing.add(cross)?;
        }
        if rhs.coefficient != 0 {
            let row = ambient_pairings(model, &lattice.ambient, &self.perpendicular)?;
            let cross = column_pair(model, &lattice.bezout, &row)?.multiply(rhs.coefficient)?;
            pairing = pairing.add(cross)?;
        }
        if self.coefficient != 0 && rhs.coefficient != 0 {
            let mixed = lattice
                .og10_form(model)?
                .generator_square
                .multiply(self.coefficient)?
                .multiply(rhs.coefficient)?;
            pairing = pairing.add(mixed)?;
        }
        Ok(pairing)
    }

    pub(crate) fn pair(&self, model: &mut Model, lattice: &Lattice, rhs: &Self) -> Result<i128> {
        self.wide_pair(model, lattice, rhs)?
            .signed()
            .ok_or_else(|| Error::plain("OG10 pairing exceeds the signed 128-bit range"))
    }

    pub(crate) fn divisibility(&self, model: &mut Model, lattice: &Lattice) -> Result<i128> {
        let form = lattice.og10_form(model)?;
        let mut row = ambient_pairings(model, &lattice.ambient, &self.perpendicular)?;
        let g_pair = column_pair(model, &lattice.bezout, &row)?
            .add(form.generator_square.multiply(self.coefficient)?)?;
        if self.coefficient != 0 {
            for (entry, &correction) in row.iter_mut().zip(&form.bezout_pairings) {
                *entry = entry.add(correction.multiply(self.coefficient)?)?;
            }
        }
        let mut divisor = g_pair.magnitude;
        for column in &lattice.columns {
            if divisor.small() == Some(1) {
                return Ok(1);
            }
            let pairing = column_pair(model, column, &row)?;
            divisor = divisor.gcd(model, pairing.magnitude)?;
        }
        divisor
            .signed(false)
            .ok_or_else(|| Error::plain("OG10 divisibility exceeds the signed 128-bit range"))
    }

    pub(crate) fn fujiki(
        &self,
        model: &mut Model,
        lattice: &Lattice,
        rhs: &Self,
        j: usize,
    ) -> Result<i128> {
        let left = polarized::fujiki_left_exponent(N, j)?;
        let xx = self.wide_pair(model, lattice, self)?;
        let xy = self.wide_pair(model, lattice, rhs)?;
        let yy = rhs.wide_pair(model, lattice, rhs)?;
        polarized::fujiki_pairings(model, left, j, xx, xy, yy)
    }

    pub(crate) fn riemann_roch(&self, model: &mut Model, lattice: &Lattice) -> Result<i128> {
        let square = self.wide_pair(model, lattice, self)?;
        let square = square.resize().map_err(|_| {
            Error::plain("OG10 Riemann–Roch value exceeds the signed 128-bit range")
        })?;
        riemann_roch_polynomial(model, square, N as i128)
    }
}

pub(crate) fn report(name: &str) -> String {
    format!(
        "rvOG10Lattice[{name}] = n=5, dimension=10, rank=24, signature=(3,21), determinant=-3, discriminant=Z/3, q=-2/3 mod 2, index=2"
    )
}

fn coordinate_range() -> Error {
    Error::plain("OG10 gluing coefficient exceeds the signed 128-bit range")
}
