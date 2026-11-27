use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::syntax::Label;

use super::{decimal, inverse, overflow, storage, Integer};

// The initial path frame can have entries of up to 384 bits. Each wall
// reflection preserves the target determinant and never increases its trace.
// Relative to the original period this bounds the basis magnitudes by 256
// bits, so transforming the original Fourier index takes fewer than 650 bits.
type Coordinate = Wide<24>;
type Basis = [[i128; 2]; 2];
const IDENTITY: Basis = [[1, 0], [0, 1]];
const MAX_WALLS: usize = 256;

pub(crate) struct Period {
    original: [i128; 3],
    reduced: [i128; 3],
    basis: Basis,
}

pub(crate) enum Command {
    Declare(Label, [i128; 3]),
    Data(Label),
    Coefficient(Label, [i128; 3]),
    Path(Label, Label, [i128; 3]),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<Option<String>> {
        let output = match self {
            Self::Declare(name, original) => {
                model.register_siegel_chamber(name, original)?;
                return Ok(None);
            }
            Self::Data(name) => {
                let period = model.siegel_chamber(&name)?;
                format!(
                    "rvSiegelChamberData[{}] = reduced={:?}, basis={:?}",
                    name.text, period.reduced, period.basis
                )
            }
            Self::Coefficient(name, index @ [m, n, r]) => {
                let period = model.siegel_chamber(&name)?;
                let transformed =
                    Symmetric::fourier(index)?.congruent(dual_basis(period.basis)?)?;
                let value = match transformed.indices()? {
                    Some((m, n, r)) => inverse::coefficient(
                        model,
                        m,
                        n,
                        r.resize().map_err(|_| overflow())?,
                        inverse::Chamber::Inside,
                    )?,
                    None => Integer::from_i128(0),
                };
                format!(
                    "rvIgusaChamberCoeff[{}, {m}, {n}, {r}] = {}",
                    name.text,
                    decimal(value)
                )
            }
            Self::Path(start, target, index @ [m, n, r]) => {
                let first = model.siegel_chamber(&start)?;
                let last = model.siegel_chamber(&target)?;
                let (walls, jump) = path(model, &first, &last, index)?;
                format!(
                    "rvIgusaChamberPath[{}, {}, {m}, {n}, {r}] = walls={walls}, jump={}",
                    start.text,
                    target.text,
                    decimal(jump)
                )
            }
        };
        Ok(Some(output))
    }
}

impl Period {
    pub(crate) fn new(model: &mut Model, original: [i128; 3]) -> Result<Self> {
        let [mut a, mut b, mut c] = original;
        if a <= 0
            || c <= 0
            || !Coordinate::product(a, c)
                .magnitude
                .cmp(Coordinate::product(b, b).magnitude)
                .is_gt()
        {
            return Err(Error::plain(
                "Siegel chamber period must be positive definite",
            ));
        }
        let mut basis = IDENTITY;
        loop {
            model.spend(1)?;
            if a > c {
                std::mem::swap(&mut a, &mut c);
                for row in &mut basis {
                    row.swap(0, 1);
                }
            }
            // Subtract the nearest multiple of the shorter first column.
            // The new diagonal is c − k(b + b′); exact cancellation can
            // be large even though this diagonal never increases.
            let remainder = b.rem_euclid(a);
            let round_up = remainder > a - remainder;
            let k = b.div_euclid(a) + i128::from(round_up);
            let new_b = if round_up { remainder - a } else { remainder };
            if k != 0 {
                c = subtract_products(c, k, b, new_b)?;
                for [x, y] in &mut basis {
                    *y = subtract_products(*y, k, *x, 0)?;
                }
                b = new_b;
            }
            if a <= c {
                break;
            }
        }
        if b == 0 {
            return Err(boundary());
        }
        if b < 0 {
            b = -b;
            for row in &mut basis {
                row[1] = checked_negate(row[1])?;
            }
        }
        Ok(Self {
            original,
            reduced: [a, b, c],
            basis,
        })
    }
}

fn subtract_products(value: i128, k: i128, first: i128, second: i128) -> Result<i128> {
    Coordinate::from_i128(value)
        .add(Coordinate::product(k, first).negate())?
        .add(Coordinate::product(k, second).negate())?
        .signed()
        .ok_or_else(geometry_overflow)
}

fn checked_negate(value: i128) -> Result<i128> {
    value.checked_neg().ok_or_else(geometry_overflow)
}

fn dual_basis([[a, b], [c, d]]: Basis) -> Result<Basis> {
    // U⁻ᵀ equals this adjugate transpose up to an overall sign, which
    // disappears in the symmetric congruence for the Fourier index.
    Ok([[d, checked_negate(c)?], [checked_negate(b)?, a]])
}

#[derive(Clone, Copy)]
struct Symmetric([Coordinate; 3]);

impl Symmetric {
    fn fourier([m, n, r]: [i128; 3]) -> Result<Self> {
        let entries = [
            Coordinate::from_i128(m).multiply(2)?,
            Coordinate::from_i128(r),
            Coordinate::from_i128(n).multiply(2)?,
        ];
        Ok(Self(entries))
    }

    fn congruent(self, [[a, b], [c, d]]: Basis) -> Result<Self> {
        let (left, right) = ([a, c], [b, d]);
        Ok(Self([
            self.pair(left, left)?,
            self.pair(left, right)?,
            self.pair(right, right)?,
        ]))
    }

    fn pair(self, [a, b]: [i128; 2], [c, d]: [i128; 2]) -> Result<Coordinate> {
        let [first, cross, last] = self.0;
        let x = first.multiply(a)?.add(cross.multiply(b)?)?;
        let y = cross.multiply(a)?.add(last.multiply(b)?)?;
        x.multiply(c)?.add(y.multiply(d)?)
    }

    fn indices(self) -> Result<Option<(i128, i128, Coordinate)>> {
        let [a, r, c] = self.0;
        // The diagonal entries of 2T stay even under every integral frame.
        let degree = |value: Coordinate| {
            let (magnitude, remainder) = value.magnitude.divide_u32(2);
            debug_assert_eq!(remainder, 0);
            Coordinate { magnitude, ..value }
        };
        let (m, n) = (degree(a), degree(c));
        if [m, n]
            .iter()
            .any(|degree| degree.negative && degree.magnitude.small().is_none_or(|value| value > 1))
        {
            return Ok(None);
        }
        Ok(Some((
            m.signed().ok_or_else(storage)?,
            n.signed().ok_or_else(storage)?,
            r,
        )))
    }

    fn next_wall(self) -> Result<Option<Wall>> {
        let [a, b, c] = self.0;
        if b.is_zero() {
            return Err(boundary());
        }
        if b.negative {
            return Ok(Some(Wall::Sign));
        }
        if b.magnitude == a.magnitude || b.magnitude == c.magnitude {
            return Err(boundary());
        }
        if b.magnitude.cmp(a.magnitude).is_gt() {
            Ok(Some(Wall::First))
        } else if b.magnitude.cmp(c.magnitude).is_gt() {
            Ok(Some(Wall::Second))
        } else {
            Ok(None)
        }
    }
}

#[derive(Clone, Copy)]
enum Wall {
    Sign,
    First,
    Second,
}

impl Wall {
    fn reflection(self) -> Basis {
        match self {
            Self::Sign => [[1, 0], [0, -1]],
            Self::First => [[1, 2], [0, -1]],
            Self::Second => [[-1, 0], [2, 1]],
        }
    }

    fn frame(self) -> Basis {
        match self {
            Self::Sign => IDENTITY,
            Self::First => [[1, 1], [0, -1]],
            Self::Second => [[0, -1], [1, 1]],
        }
    }
}

fn path(
    model: &mut Model,
    start: &Period,
    target: &Period,
    index: [i128; 3],
) -> Result<(usize, Integer)> {
    let mut relative =
        Symmetric(target.original.map(Coordinate::from_i128)).congruent(start.basis)?;
    let mut walls = Vec::new();
    while let Some(wall) = relative.next_wall()? {
        if walls.len() == MAX_WALLS {
            return Err(Error::plain(format!(
                "Siegel chamber path exceeds the {MAX_WALLS}-wall limit"
            )));
        }
        model.spend(1)?;
        relative = relative.congruent(wall.reflection())?;
        walls.push(wall);
    }
    // In every wall frame, the current chamber is on the positive side.
    // Add its inside-minus-outside jump, then change the current basis.
    let mut fourier = Symmetric::fourier(index)?.congruent(dual_basis(start.basis)?)?;
    let mut total = Integer::from_i128(0);
    for wall in &walls {
        let local = fourier.congruent(dual_basis(wall.frame())?)?;
        if !local.0[1].is_zero() {
            if let Some((m, n, r)) = local.indices()? {
                let jump = inverse::wall_jump(model, m, n, r.resize().map_err(|_| overflow())?)?;
                total = total.add(jump).map_err(|_| overflow())?;
            }
        }
        fourier = fourier.congruent(dual_basis(wall.reflection())?)?;
    }
    Ok((walls.len(), total))
}

fn boundary() -> Error {
    Error::plain("Siegel chamber period lies on an integral separating wall")
}

fn geometry_overflow() -> Error {
    Error::plain("Siegel chamber reduction exceeds the signed 128-bit coordinate range")
}
