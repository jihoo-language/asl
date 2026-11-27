use std::cell::OnceCell;
use std::rc::Rc;

use super::{polarized, unimodular, Integer, K3Vector, Vector};
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod og10;

const AMBIENT_RANK: usize = unimodular::RANK + 2;
const RANK: usize = AMBIENT_RANK - 1;
// Each ambient pairing has at most 132 bits. A kernel column uses at
// most 23 successive Bézout multipliers of that size, with one more
// multiplication and a 24-term sum to pair or realize the column.
// The isotropic case also checks for overflow in its second elimination.
const COLUMN_LIMBS: usize = 104;
const COLUMN_CELLS: usize = COLUMN_LIMBS + 1;
const PAIRING_CELLS: usize = std::mem::size_of::<Integer>().div_ceil(std::mem::size_of::<u32>());
type ColumnInteger = Wide<COLUMN_LIMBS>;
type Column = [ColumnInteger; AMBIENT_RANK];
type PairingRow = [ColumnInteger; AMBIENT_RANK];

pub(crate) enum Expression {
    Basis(Label, usize),
    Glue(Label),
    Dual(Label),
    Quotient(Label, Box<ParsedExpr>),
    Lift(Box<ParsedExpr>),
}

pub(crate) enum Command {
    Define(Label, Box<ParsedExpr>),
    Lattice(Label),
    Divisibility(Label, Box<ParsedExpr>),
    Fujiki(Label, usize, Box<[ParsedExpr; 2]>),
    RiemannRoch(Label, Box<ParsedExpr>),
    QuotientPair(Box<[ParsedExpr; 2]>),
}

pub(crate) struct Lattice {
    pub(crate) marking: Rc<str>,
    vector: Rc<Vector>,
    square: i128,
    ambient: Vec<Vector>,
    columns: Vec<Column>,
    bezout: Column,
    dual_correction: Option<ColumnInteger>,
    og10: OnceCell<Box<og10::Form>>,
}

impl Lattice {
    pub(crate) fn new(model: &mut Model, marking: Rc<str>, vector: Rc<Vector>) -> Result<Self> {
        let square = vector.pair(model, &vector)?;
        if square < 0 || square % 2 != 0 {
            return Err(Error::plain(
                "rvmukaimoduli requires a Mukai vector of positive even square or a primitive nonzero isotropic vector",
            ));
        }
        let zero = Rc::new(K3Vector::zero());
        let mut ambient = vec![
            Vector {
                r: 1,
                c: zero.clone(),
                s: 0,
            },
            Vector {
                r: 0,
                c: zero,
                s: 1,
            },
        ];
        ambient.extend(
            unimodular::integral_basis(model)?
                .into_iter()
                .map(|c| Vector {
                    r: 0,
                    c: Rc::new(c),
                    s: 0,
                }),
        );
        let row = ambient_pairings(model, &ambient, &vector)?;
        let mut columns = vec![[column_integer(0); AMBIENT_RANK]; AMBIENT_RANK];
        for (i, column) in columns.iter_mut().enumerate() {
            column[i] = column_integer(1);
        }
        let bezout = eliminate(
            model,
            &mut columns,
            row.to_vec(),
            "rvmukaimoduli requires a primitive Mukai vector",
        )?;
        let dual_correction = if square == 0 {
            let row = bezout_pairings(model, &ambient, &bezout)?;
            let square = column_pair(model, &bezout, &row)?;
            let (half, remainder) = square.magnitude.divide_u32(2);
            if remainder != 0 {
                return Err(Error::plain("Mukai lattice must have even squares"));
            }
            let kernel_row = columns
                .iter()
                .map(|column| column_pair(model, column, &row))
                .collect::<Result<_>>()?;
            // On v-perp, (w, x) = (w0, x). This covector has value 1
            // on v, so its integral kernel is the desired complement.
            eliminate(
                model,
                &mut columns,
                kernel_row,
                "Mukai isotropic kernel covector is not primitive",
            )?;
            Some(ColumnInteger {
                magnitude: half,
                negative: square.negative,
            })
        } else {
            None
        };
        Ok(Self {
            marking,
            vector,
            square,
            ambient,
            columns,
            bezout,
            dual_correction,
            og10: OnceCell::new(),
        })
    }

    pub(crate) fn cells(&self) -> usize {
        4 + AMBIENT_RANK * (AMBIENT_RANK * COLUMN_CELLS + AMBIENT_RANK + 1) + COLUMN_CELLS
    }

    pub(crate) fn basis(&self, model: &mut Model, index: usize) -> Result<Vector> {
        let column = self.columns.get(index).ok_or_else(|| {
            Error::plain(format!(
                "rvMukaiModuliBasis index must be between 0 and {}",
                self.columns.len() - 1
            ))
        })?;
        realize(self.raw_column(model, column, 1)?)
    }

    pub(crate) fn glue(&self, model: &mut Model) -> Result<Vector> {
        if self.square == 0 {
            return Ok(Vector::zero());
        }
        let mut coordinates = self.raw_column(model, &self.bezout, self.square)?;
        add_scaled(model, &mut coordinates, &self.vector, column_integer(-1))?;
        realize(coordinates)
    }

    pub(crate) fn require_isotropic(&self, command: &str) -> Result<()> {
        self.isotropic_correction(command).map(|_| ())
    }

    fn isotropic_correction(&self, command: &str) -> Result<ColumnInteger> {
        self.dual_correction.ok_or_else(|| {
            Error::plain(format!(
                "{command} requires an isotropic Mukai moduli lattice"
            ))
        })
    }

    pub(crate) fn dual(&self, model: &mut Model) -> Result<Vector> {
        let correction = self.isotropic_correction("rvMukaiModuliDual")?;
        let mut coordinates = self.raw_column(model, &self.bezout, 1)?;
        add_scaled(model, &mut coordinates, &self.vector, correction.negate())?;
        realize(coordinates)
    }

    pub(crate) fn project(&self, model: &mut Model, x: &Vector) -> Result<Vector> {
        self.require_isotropic("rvMukaiQuotient")?;
        self.require_perpendicular(model, x, "rvMukaiQuotient")?;
        let row = ambient_pairings(model, &self.ambient, x)?;
        let shift = column_pair(model, &self.bezout, &row)?;
        let mut coordinates = raw_coordinates(x).map(column_integer);
        add_scaled(model, &mut coordinates, &self.vector, shift.negate())?;
        realize(coordinates)
    }

    pub(crate) fn report(&self, name: &str) -> String {
        let n = self.n();
        let square = self.square;
        if square == 0 {
            return format!(
                "rvMukaiModuliLattice[{name}] = n=1, dimension=2, rank=22, signature=(3,19), determinant=-1, discriminant=trivial"
            );
        }
        format!(
            "rvMukaiModuliLattice[{name}] = n={n}, dimension={}, rank={RANK}, signature=(3,20), determinant={square}, discriminant=Z/{square}, q=-1/{square} mod 2",
            2 * n as u128
        )
    }

    pub(crate) fn divisibility(&self, model: &mut Model, x: &Vector) -> Result<i128> {
        self.require_perpendicular(model, x, "rvMukaiModuliDiv")?;
        let pairings = ambient_pairings(model, &self.ambient, x)?;
        let mut divisor = column_integer(0).magnitude;
        for column in &self.columns {
            let pairing = column_pair(model, column, &pairings)?;
            divisor = divisor.gcd(model, pairing.magnitude)?;
            if divisor.small() == Some(1) {
                return Ok(1);
            }
        }
        divisor.signed(false).ok_or_else(|| {
            Error::plain("Mukai moduli divisibility exceeds the signed 128-bit range")
        })
    }

    pub(crate) fn fujiki(
        &self,
        model: &mut Model,
        x: &Vector,
        y: &Vector,
        j: usize,
    ) -> Result<i128> {
        self.require_perpendicular(model, x, "rvMukaiModuliFujiki")?;
        self.require_perpendicular(model, y, "rvMukaiModuliFujiki")?;
        let n = usize::try_from(self.n())
            .map_err(|_| Error::plain("Mukai moduli Fujiki dimension exceeds the work limit"))?;
        let left = polarized::fujiki_left_exponent(n, j)?;
        // For isotropic v, replacing either perpendicular input by
        // its canonical representative does not change any pairing.
        let xx = x.wide_pair(model, x)?;
        let xy = x.wide_pair(model, y)?;
        let yy = y.wide_pair(model, y)?;
        polarized::fujiki_pairings(model, left, j, xx, xy, yy)
    }

    pub(crate) fn riemann_roch(&self, model: &mut Model, x: &Vector) -> Result<i128> {
        self.require_perpendicular(model, x, "rvMukaiModuliRR")?;
        let square = x.wide_pair(model, x)?;
        riemann_roch_polynomial(model, square, self.n())
    }

    fn n(&self) -> i128 {
        self.square / 2 + 1
    }

    fn require_perpendicular(&self, model: &mut Model, x: &Vector, command: &str) -> Result<()> {
        if !x.wide_pair(model, &self.vector)?.is_zero() {
            return Err(Error::plain(format!(
                "{command} requires a Mukai vector perpendicular to the moduli vector"
            )));
        }
        Ok(())
    }

    fn raw_column(&self, model: &mut Model, column: &Column, scale: i128) -> Result<Column> {
        let mut result = [column_integer(0); AMBIENT_RANK];
        for (&coefficient, vector) in column.iter().zip(&self.ambient) {
            if !coefficient.is_zero() {
                add_scaled(model, &mut result, vector, coefficient.multiply(scale)?)?;
            }
        }
        Ok(result)
    }
}

fn add_scaled(
    model: &mut Model,
    coordinates: &mut Column,
    vector: &Vector,
    scale: ColumnInteger,
) -> Result<()> {
    model.spend(AMBIENT_RANK * COLUMN_CELLS)?;
    for (target, value) in coordinates.iter_mut().zip(raw_coordinates(vector)) {
        if value != 0 {
            *target = target.add(scale.multiply(value)?)?;
        }
    }
    Ok(())
}

fn realize(coordinates: Column) -> Result<Vector> {
    let mut output = [0; AMBIENT_RANK];
    for (out, value) in output.iter_mut().zip(coordinates) {
        *out = value.signed().ok_or_else(super::coordinate_range)?;
    }
    let (wedge, exceptional) = output[2..].split_at(6);
    Ok(Vector {
        r: output[0],
        c: Rc::new(K3Vector::new(
            wedge.try_into().unwrap(),
            exceptional.try_into().unwrap(),
        )?),
        s: output[1],
    })
}

fn column_integer(value: i128) -> ColumnInteger {
    ColumnInteger::from_i128(value)
}

fn raw_coordinates(vector: &Vector) -> [i128; AMBIENT_RANK] {
    let mut output = [0; AMBIENT_RANK];
    output[0] = vector.r;
    output[1] = vector.s;
    let (wedge, exceptional) = vector.c.coordinates();
    output[2..8].copy_from_slice(wedge);
    output[8..].copy_from_slice(exceptional);
    output
}

fn ambient_pairings(model: &mut Model, ambient: &[Vector], vector: &Vector) -> Result<PairingRow> {
    let mut row = [column_integer(0); AMBIENT_RANK];
    for (entry, basis) in row.iter_mut().zip(ambient) {
        *entry = vector.wide_pair(model, basis)?.resize()?;
    }
    Ok(row)
}

fn bezout_pairings(model: &mut Model, ambient: &[Vector], bezout: &Column) -> Result<PairingRow> {
    let mut row = [column_integer(0); AMBIENT_RANK];
    for (entry, value) in row.iter_mut().zip(ambient) {
        let pairings = ambient_pairings(model, ambient, value)?;
        *entry = column_pair(model, bezout, &pairings)?;
    }
    Ok(row)
}

fn column_pair(model: &mut Model, column: &Column, row: &[ColumnInteger]) -> Result<ColumnInteger> {
    debug_assert_eq!(row.len(), AMBIENT_RANK);
    model.spend(AMBIENT_RANK * COLUMN_CELLS)?;
    let mut result = column_integer(0);
    for (&coordinate, &entry) in column.iter().zip(row) {
        if !coordinate.is_zero() && !entry.is_zero() {
            result = result.add(coordinate.multiply_wide(entry)?)?;
        }
    }
    Ok(result)
}

// Split a primitive covector of a given integral basis as (1, 0, …).
fn eliminate(
    model: &mut Model,
    columns: &mut Vec<Column>,
    mut row: Vec<ColumnInteger>,
    nonprimitive: &str,
) -> Result<Column> {
    let pivot = (0..row.len())
        .filter(|&i| !row[i].is_zero())
        .min_by(|&i, &j| row[i].magnitude.cmp(row[j].magnitude))
        .ok_or_else(|| Error::plain(nonprimitive))?;
    columns.swap(0, pivot);
    row.swap(0, pivot);
    let mut a = row[0];
    if a.negative {
        a = a.negate();
        for entry in &mut columns[0] {
            *entry = entry.negate();
        }
    }

    for (i, &b) in row.iter().enumerate().skip(1) {
        if b.is_zero() {
            continue;
        }
        let (g, p, q) = bezout(model, a, b)?;
        let x = a.divide_exact(model, g)?.expect("gcd divides left entry");
        let y = b
            .divide_exact(model, g)?
            .expect("gcd divides right entry")
            .negate();
        let old = columns[0];
        model.spend(AMBIENT_RANK * 4 * COLUMN_CELLS)?;
        for (j, entry) in columns[i].into_iter().enumerate() {
            columns[0][j] = old[j].multiply_wide(p)?.add(entry.multiply_wide(q)?)?;
            columns[i][j] = old[j].multiply_wide(y)?.add(entry.multiply_wide(x)?)?;
        }
        // This two-column operation has determinant (pa + qb)/g = 1.
        a = g;
    }
    if a.magnitude.small() != Some(1) {
        return Err(Error::plain(nonprimitive));
    }
    Ok(columns.remove(0))
}

fn bezout(
    model: &mut Model,
    a: ColumnInteger,
    b: ColumnInteger,
) -> Result<(ColumnInteger, ColumnInteger, ColumnInteger)> {
    let [zero, one] = [column_integer(0), column_integer(1)];
    if a.magnitude == one.magnitude {
        return Ok((one, one, zero));
    }
    let mut left = a.magnitude;
    let mut right = b.magnitude;
    let mut left_coefficients = [one, zero];
    let mut right_coefficients = [zero, one];
    while !right.is_zero() {
        let (quotient, remainder) = left.div_rem(model, right)?;
        let quotient = ColumnInteger {
            magnitude: quotient,
            negative: false,
        };
        (left, right) = (right, remainder);
        model.spend(4 * COLUMN_CELLS)?;
        for (x, y) in left_coefficients.iter_mut().zip(&mut right_coefficients) {
            (*x, *y) = (*y, x.add(quotient.multiply_wide(*y)?.negate())?);
        }
    }
    let [p, q] = left_coefficients;
    let q = if b.negative { q.negate() } else { q };
    Ok((
        ColumnInteger {
            magnitude: left,
            negative: false,
        },
        p,
        q,
    ))
}

fn riemann_roch_polynomial(model: &mut Model, q: Integer, n: i128) -> Result<i128> {
    let (half, remainder) = q.magnitude.divide_u32(2);
    if remainder != 0 {
        return Err(Error::plain("Mukai lattice must have even squares"));
    }
    let upper = Integer {
        magnitude: half,
        negative: q.negative,
    }
    .add(Integer::from_i128(n + 1))?;
    generalized_binomial(model, upper, n)
}

fn generalized_binomial(model: &mut Model, mut upper: Integer, n: i128) -> Result<i128> {
    let negative = upper.negative && !upper.is_zero();
    let negative_output = negative && n % 2 != 0;
    if negative {
        // binom(m,n) = (−1)^n binom(n−m−1,n) for m < 0.
        upper = upper.negate().add(Integer::from_i128(n - 1))?;
    }
    let n_wide = Integer::from_i128(n);
    if upper.magnitude.cmp(n_wide.magnitude).is_lt() {
        return Ok(0);
    }
    let complement = upper.magnitude.difference(n_wide.magnitude);
    let k = complement.small().unwrap_or(n as u128).min(n as u128);
    let start = upper.add(Integer::from_i128(-(k as i128)))?;
    let mut value = one();
    for i in 1..=k {
        model.spend(PAIRING_CELLS * 4)?;
        let numerator = start.add(Integer::from_i128(i as i128))?;
        value = value
            .multiply_wide(numerator)?
            .divide_exact(model, Integer::from_i128(i as i128))?
            .expect("binomial recurrence is integral");
        if value.magnitude.signed(negative_output).is_none() {
            return Err(rr_range());
        }
    }
    value.magnitude.signed(negative_output).ok_or_else(rr_range)
}

fn rr_range() -> Error {
    Error::plain("Riemann–Roch value exceeds the signed 128-bit range")
}

fn one() -> Integer {
    Integer::from_i128(1)
}
