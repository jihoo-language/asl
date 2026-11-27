use std::cell::OnceCell;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::jacobian_level::Basis;
use crate::model::Model;
use crate::rv_qp_quad::Wide;
use crate::syntax::{Label, ParsedExpr};

use super::nikulin::{self, AFFINE};
use super::polarized::{self, Divisor, COEFFICIENTS};

mod points;
mod spectrum;

const RANK: usize = COEFFICIENTS;
const GLUE_DIMENSION: usize = 1 + AFFINE.len();
type Coordinates = [i128; RANK];
type Integer = Wide<8>;

#[derive(Clone, Copy)]
pub(crate) enum Generator {
    Translate(u8),
    Switch(u8),
    Project(u8),
}

pub(crate) enum Command {
    Define(Label, Label, Vec<Generator>),
    Entropy(Label),
    Cycles(Label),
}

pub(crate) struct Expression(pub(crate) Label, pub(crate) Box<ParsedExpr>);
pub(crate) struct PointExpression(
    pub(crate) Label,
    pub(crate) i128,
    pub(crate) Box<ParsedExpr>,
);

pub(crate) struct Word {
    pub(crate) basis_name: Rc<str>,
    pub(crate) generators: Box<[Generator]>,
    basis: Rc<Basis>,
    origin: u16,
    point_geometry: OnceCell<points::Geometry>,
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<Option<String>> {
        match self {
            Self::Define(name, basis_name, generators) => {
                let basis = polarized::basis(model, &basis_name)?;
                let origin = polarized::trope_origin(model, &basis)?;
                let word = Word {
                    basis_name: basis_name.text.into(),
                    generators: generators.into_boxed_slice(),
                    basis,
                    origin,
                    point_geometry: OnceCell::new(),
                };
                model.define_k3_automorphism(name, word)?;
                Ok(None)
            }
            Self::Entropy(name) => {
                let word = model.k3_automorphism(&name)?;
                Ok(Some(word.entropy(model, &name.text)?))
            }
            Self::Cycles(name) => {
                let word = model.k3_automorphism(&name)?;
                Ok(Some(word.cycles(model, &name.text)?))
            }
        }
    }
}

impl Word {
    pub(crate) fn apply(&self, model: &mut Model, divisor: &Divisor) -> Result<Divisor> {
        let image = self.action(model, *divisor.coefficients())?;
        Divisor::literal(model, &self.basis, image)
    }

    fn action(&self, model: &mut Model, mut coordinates: Coordinates) -> Result<Coordinates> {
        for &generator in &self.generators {
            coordinates = generator.apply(model, self.origin, &coordinates)?;
        }
        Ok(coordinates)
    }

    fn entropy(&self, model: &mut Model, name: &str) -> Result<String> {
        let lattice = IntegralBasis::new(model, self.origin)?;
        let mut matrix = [[0; RANK]; RANK];
        for (column, vector) in lattice.vectors().enumerate() {
            let image = self.action(model, vector)?;
            let coordinates = lattice.coordinates(model, &image)?;
            for (row, coefficient) in matrix.iter_mut().zip(coordinates) {
                row[column] = coefficient;
            }
        }
        spectrum::entropy(model, name, matrix)
    }
}

impl Generator {
    fn apply(self, model: &mut Model, origin: u16, x: &Coordinates) -> Result<Coordinates> {
        let mut result = *x;
        match self {
            Self::Translate(mask) => {
                model.spend(RANK)?;
                for (node, &coefficient) in x[1..].iter().enumerate() {
                    result[(node ^ usize::from(mask)) + 1] = coefficient;
                }
            }
            Self::Switch(mask) => {
                model.spend(RANK * 8 * 8)?;
                let total = x[1..].iter().try_fold(Integer::from_i128(0), |sum, &c| {
                    sum.add(Integer::from_i128(c))
                })?;
                result[0] = signed(Integer::product(3, x[0]).add(half(total)?)?)?;
                let hyperplane = Integer::from_i128(x[0]);
                for (node, coefficient) in result[1..].iter_mut().enumerate() {
                    // t_mask follows the Gauss switch. The membership test is
                    // for this output node in the translated trope of each E_v.
                    let translated = node ^ usize::from(mask);
                    let sum = x[1..]
                        .iter()
                        .enumerate()
                        .filter(|(v, _)| origin & (1 << (translated ^ v)) != 0)
                        .try_fold(Integer::from_i128(0), |sum, (_, &c)| {
                            sum.add(Integer::from_i128(c))
                        })?;
                    *coefficient = signed(hyperplane.add(half(sum)?)?.negate())?;
                }
            }
            Self::Project(mask) => {
                model.spend(80)?;
                let index = usize::from(mask) + 1;
                result[0] = signed(Integer::product(3, x[0]).add(Integer::product(2, x[index]))?)?;
                result[index] =
                    signed(Integer::product(-4, x[0]).add(Integer::product(-3, x[index]))?)?;
            }
        }
        Ok(result)
    }
}

fn half(value: Integer) -> Result<Integer> {
    let (magnitude, parity) = value.magnitude.divide_u32(2);
    if parity != 0 {
        return Err(Error::plain(
            "K3 automorphism must preserve integral marked lattice coordinates",
        ));
    }
    Ok(Integer {
        magnitude,
        negative: value.negative,
    })
}

fn signed(value: Integer) -> Result<i128> {
    value
        .signed()
        .ok_or_else(|| Error::plain("K3 automorphism coefficient exceeds signed 128-bit range"))
}

struct IntegralBasis {
    rows: [u32; GLUE_DIMENSION],
    pivots: Vec<usize>,
    nonpivots: Vec<usize>,
}

impl IntegralBasis {
    fn new(model: &mut Model, origin: u16) -> Result<Self> {
        model.spend(RANK * GLUE_DIMENSION * 2)?;
        let mut rows = [0; GLUE_DIMENSION];
        rows[0] = 1 | (u32::from(origin) << 1);
        for (row, affine) in rows[1..].iter_mut().zip(AFFINE) {
            *row = u32::from(affine) << 1;
        }
        let pivots = nikulin::binary_rref(&mut rows, RANK);
        if pivots.len() != GLUE_DIMENSION {
            return Err(Error::plain(
                "marked Kummer divisor gluing code has the wrong rank",
            ));
        }
        let nonpivots = (0..RANK)
            .filter(|column| !pivots.contains(column))
            .collect();
        Ok(Self {
            rows,
            pivots,
            nonpivots,
        })
    }

    fn vectors(&self) -> impl Iterator<Item = Coordinates> + '_ {
        self.rows
            .iter()
            .map(|row| std::array::from_fn(|i| i128::from((row >> i) & 1)))
            .chain(self.nonpivots.iter().map(|&i| {
                let mut vector = [0; RANK];
                vector[i] = 2;
                vector
            }))
    }

    fn coordinates(&self, model: &mut Model, vector: &Coordinates) -> Result<Coordinates> {
        model.spend(8 * RANK * GLUE_DIMENSION)?;
        let mut result = [0; RANK];
        for (i, &pivot) in self.pivots.iter().enumerate() {
            result[i] = vector[pivot];
        }
        for (i, &column) in self.nonpivots.iter().enumerate() {
            let mut remainder = Integer::from_i128(vector[column]);
            for (&row, &pivot) in self.rows.iter().zip(&self.pivots) {
                if row & (1 << column) != 0 {
                    remainder = remainder.add(Integer::from_i128(vector[pivot]).negate())?;
                }
            }
            result[GLUE_DIMENSION + i] = signed(half(remainder)?)?;
        }
        Ok(result)
    }
}
