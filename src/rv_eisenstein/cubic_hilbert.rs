use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_unram::Ring;
use crate::syntax::{Label, ParsedExpr};

use super::kummer::{divide_by_lambda, inverse_base_unit, Extension as Kummer, Generator};
use super::{Element as BaseElement, Extension as Base};

pub(crate) const CLASS_PRECISION: usize = 4;

pub(crate) struct LocalSpec(pub(crate) Label, pub(crate) usize);

pub(crate) enum Expression {
    Local(LocalSpec, i128, Box<ParsedExpr>),
    FromIntegral(usize, Box<ParsedExpr>),
    LocalMultiply(Box<[ParsedExpr; 2]>),
    LocalInverse(Box<ParsedExpr>),
    OneMinus(usize, Box<ParsedExpr>),
    Class(Box<ParsedExpr>),
    ClassLiteral(Label, [i128; 4]),
    Pair(Pair, Box<[ParsedExpr; 2]>),
}

pub(crate) enum Pair {
    Symbol,
    IsNorm,
    Splits,
}

pub(crate) enum Report {
    Certificate(Label),
    NormSubgroup(Box<ParsedExpr>),
    Invariant(Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) field: Box<str>,
    pub(crate) precision: usize,
}

impl Grade {
    pub(crate) fn resolve(
        model: &Model,
        LocalSpec(name, precision): LocalSpec,
    ) -> Result<(Self, Rc<Base>)> {
        let field = model.eisenstein(&name)?;
        require_field(&field).map_err(|error| error.at(name.span))?;
        field.algebra(precision)?;
        Ok((
            Self {
                field: name.text.into(),
                precision,
            },
            field,
        ))
    }
}

pub(crate) fn require_field(field: &Base) -> Result<()> {
    if field.prime() != 3 || field.base.degree() != 1 || field.cyclotomic_order != Some(3) {
        return Err(Error::plain(
            "wild cubic Hilbert pairing requires a named Q_3(zeta_3) cyclotomic base",
        ));
    }
    Ok(())
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Local {
    valuation: i128,
    unit: BaseElement,
}

impl Local {
    pub(crate) fn new(
        model: &mut Model,
        field: &Base,
        precision: usize,
        valuation: i128,
        unit: BaseElement,
    ) -> Result<Self> {
        let algebra = field.algebra(precision)?;
        let unit = algebra.reduce(unit);
        if algebra.valuation(model, &unit)? != 0 {
            return Err(Error::plain(
                "the normalized cubic local factor must be a unit",
            ));
        }
        Ok(Self { valuation, unit })
    }

    pub(crate) fn from_integral(
        model: &mut Model,
        field: &Base,
        source: usize,
        precision: usize,
        value: BaseElement,
    ) -> Result<Self> {
        let algebra = field.algebra(source)?;
        let valuation = algebra.valuation(model, &value)?;
        if valuation == source {
            return Err(Error::plain(
                "cannot certify a nonzero cubic local element from zero residue",
            ));
        }
        if source - valuation < precision {
            return Err(Error::plain(
                "insufficient relative precision after removing the local valuation",
            ));
        }
        let unit = divide_by_lambda(model, field, source, value, valuation)?;
        Self::new(model, field, precision, valuation as i128, unit)
    }

    pub(crate) fn multiply(
        &self,
        model: &mut Model,
        field: &Base,
        precision: usize,
        right: &Self,
    ) -> Result<Self> {
        let valuation = self
            .valuation
            .checked_add(right.valuation)
            .ok_or_else(valuation_range)?;
        let unit = field
            .algebra(precision)?
            .multiply(model, &self.unit, &right.unit)?;
        Ok(Self { valuation, unit })
    }

    pub(crate) fn inverse(
        &self,
        model: &mut Model,
        field: &Base,
        precision: usize,
    ) -> Result<Self> {
        let valuation = self.valuation.checked_neg().ok_or_else(valuation_range)?;
        let unit = inverse_base_unit(model, &field.algebra(precision)?, &self.unit)?;
        Ok(Self { valuation, unit })
    }

    pub(crate) fn one_minus(
        &self,
        model: &mut Model,
        field: &Base,
        source: usize,
        target: usize,
    ) -> Result<Self> {
        if self.valuation == 0 {
            let algebra = field.algebra(source)?;
            let negative = algebra.negate(model, &self.unit)?;
            let value = algebra.add(model, &algebra.scalar(1), &negative)?;
            return Self::from_integral(model, field, source, target, value);
        }

        let algebra = field.algebra(target)?;
        let power = self.valuation.unsigned_abs();
        let shift = if power >= target as u128 {
            algebra.zero()
        } else {
            algebra.power(model, algebra.pi(), power as u64)?
        };
        let unit = algebra.reduce(self.unit.clone());
        let (valuation, unit) = if self.valuation > 0 {
            let small = algebra.multiply(model, &unit, &shift)?;
            let negative = algebra.negate(model, &small)?;
            (0, algebra.add(model, &algebra.scalar(1), &negative)?)
        } else {
            let negative = algebra.negate(model, &unit)?;
            (self.valuation, algebra.add(model, &shift, &negative)?)
        };
        Self::new(model, field, target, valuation, unit)
    }

    pub(crate) fn class(&self, model: &mut Model, field: &Base) -> Result<Class> {
        let algebra = field.algebra(CLASS_PRECISION)?;
        let one = algebra.scalar(1);
        let negative_one = algebra.negate(model, &one)?;
        let mut unit = algebra.reduce(self.unit.clone());
        if unit.0[0][0] % 3 == 2 {
            unit = algebra.negate(model, &unit)?;
        }
        let mut coordinates = [0; 4];
        coordinates[0] = self.valuation.rem_euclid(3) as u8;
        for (level, coordinate) in coordinates.iter_mut().enumerate().skip(1) {
            let discrepancy = algebra.add(model, &unit, &negative_one)?;
            let digit =
                divide_by_lambda(model, field, CLASS_PRECISION, discrepancy, level)?.0[0][0] % 3;
            *coordinate = digit as u8;
            if digit != 0 {
                let power = algebra.power(model, algebra.pi(), level as u64)?;
                let generator = algebra.add(model, &one, &power)?;
                let inverse = inverse_base_unit(model, &algebra, &generator)?;
                let peel = algebra.power(model, inverse, digit)?;
                unit = algebra.multiply(model, &unit, &peel)?;
            }
        }
        if unit != one {
            return Err(inconsistent());
        }
        Ok(Class(coordinates))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(7)?;
        Ok(Value::CubicLocal(Rc::new(self)))
    }

    pub(crate) fn write(&self) -> String {
        format!(
            "valuation={}, unit={:?}",
            self.valuation,
            [self.unit.0[0][0], self.unit.0[1][0]]
        )
    }
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Class(pub(crate) [u8; 4]);

impl Class {
    pub(crate) fn literal(values: [i128; 4]) -> Self {
        Self(values.map(|value| value.rem_euclid(3) as u8))
    }

    pub(crate) fn add(self, other: Self) -> Self {
        Self(std::array::from_fn(|i| (self.0[i] + other.0[i]) % 3))
    }

    pub(crate) fn negate(self) -> Self {
        Self(self.0.map(|a| (3 - a) % 3))
    }
}

pub(crate) struct Certificate {
    pub(crate) gram: [[u8; 4]; 4],
    boundary_action: u8,
    norms_of_uniformizers: [BaseElement; 2],
    unit_norm_mod27: u64,
    cyclotomic_action: u64,
}

impl Certificate {
    pub(crate) fn derive(model: &mut Model) -> Result<Self> {
        let coefficients = Rc::new(Ring::new(model, 3, 4, vec![0, 1])?);
        let base = Rc::new(Base::new(model, coefficients.clone(), Vec::new(), Some(1))?);

        let boundary = Kummer::new(model, "K".into(), base.clone(), 3, vec![vec![1]])?;
        let boundary = boundary.algebra(model, 1)?;
        let theta = boundary.generator(model, Generator::Theta)?;
        let square = boundary.multiply(model, &theta, &theta)?;
        let frobenius = boundary.multiply(model, &square, &theta)?;
        let one = boundary.embed(base.algebra(1)?.scalar(1));
        let expected = boundary.add(model, &theta, &one)?;
        if frobenius != expected {
            return Err(inconsistent());
        }
        let mut boundary_action = None;
        for c in 0..3 {
            if boundary.sigma(model, &theta, c)? == frobenius {
                boundary_action = Some(c as u8);
                break;
            }
        }
        let boundary_action = boundary_action.ok_or_else(inconsistent)?;

        let mut norms_of_uniformizers = [base.algebra(5)?.zero(), base.algebra(5)?.zero()];
        for (index, norm) in norms_of_uniformizers.iter_mut().enumerate() {
            let tower = Kummer::new(model, "K".into(), base.clone(), index + 1, vec![vec![1]])?;
            let algebra = tower.algebra(model, 13)?;
            let pi = algebra.generator(model, Generator::Uniformizer)?;
            *norm = algebra.conjugate_invariant(model, &pi, true)?;
            if *norm != base.algebra(5)?.pi() {
                return Err(inconsistent());
            }
        }

        let algebra = base.algebra(6)?;
        let lambda = algebra.pi();
        let square = algebra.multiply(model, &lambda, &lambda)?;
        let u2 = algebra.add(model, &algebra.scalar(1), &square)?;
        let conjugate = algebra.automorphism(model, &u2, 2)?;
        let norm = algebra.multiply(model, &u2, &conjugate)?;
        let unit_norm_mod27 = norm.0[0][0];
        if norm.0[1][0] != 0 || unit_norm_mod27 % 3 == 0 {
            return Err(inconsistent());
        }
        // Arithmetic reciprocity on p-power roots of unity uses the inverse
        // unit. Norm compatibility reduces the K-Artin action to Q_3.
        let cyclotomic_action = (1..9)
            .find(|&x| (unit_norm_mod27 * x) % 9 == 1)
            .ok_or_else(inconsistent)?;
        let ninth = Base::new(model, coefficients, Vec::new(), Some(2))?;
        let algebra = ninth.algebra(12)?;
        let root = algebra.add(model, &algebra.scalar(1), &algebra.pi())?;
        let third = algebra.power(model, root.clone(), 3)?;
        let image = algebra.automorphism(model, &root, cyclotomic_action)?;
        let mut cyclotomic_exponent = None;
        for c in 0..3 {
            let multiplier = algebra.power(model, third.clone(), c)?;
            if algebra.multiply(model, &root, &multiplier)? == image {
                cyclotomic_exponent = Some(c as u8);
                break;
            }
        }
        let cyclotomic_exponent = cyclotomic_exponent.ok_or_else(inconsistent)?;
        if boundary_action == 0 || cyclotomic_exponent == 0 {
            return Err(inconsistent());
        }

        // The two checked relative norms kill the entries with lambda in
        // those rows. Unramified reciprocity kills all unit entries in
        // the boundary row. Alternation determines the transpose entries.
        let mut gram = [[0; 4]; 4];
        gram[3][0] = boundary_action;
        gram[0][3] = 3 - boundary_action;
        gram[1][2] = cyclotomic_exponent;
        gram[2][1] = 3 - cyclotomic_exponent;
        Ok(Self {
            gram,
            boundary_action,
            norms_of_uniformizers,
            unit_norm_mod27,
            cyclotomic_action,
        })
    }

    pub(crate) fn exponent(&self, left: Class, right: Class) -> u8 {
        let row = self.norm_functional(left);
        row.iter()
            .zip(right.0)
            .map(|(&a, b)| u16::from(a) * u16::from(b))
            .sum::<u16>()
            .rem_euclid(3) as u8
    }

    fn norm_functional(&self, left: Class) -> [u8; 4] {
        std::array::from_fn(|j| {
            left.0
                .iter()
                .zip(&self.gram)
                .map(|(&a, row)| a * row[j])
                .sum::<u8>()
                % 3
        })
    }

    pub(crate) fn norm_subgroup(&self, left: Class) -> String {
        let row = self.norm_functional(left);
        let pivot = row.iter().position(|&a| a != 0);
        let mut basis = Vec::new();
        for i in 0..4 {
            if pivot == Some(i) {
                continue;
            }
            let mut vector = [0; 4];
            vector[i] = 1;
            if let Some(j) = pivot {
                // Nonzero elements of F3 are their own inverses.
                vector[j] = (3 - (row[i] * row[j]) % 3) % 3;
            }
            basis.push(vector);
        }
        format!(
            "rvCubicNormSubgroup = functional={row:?}, basis={basis:?}, order={}",
            3_usize.pow(basis.len() as u32)
        )
    }

    pub(crate) fn write(&self, name: &str) -> String {
        let norms = self
            .norms_of_uniformizers
            .each_ref()
            .map(|x| [x.0[0][0], x.0[1][0]]);
        let determinant = (self.gram[3][0] * self.gram[1][2]).pow(2) % 3;
        format!("rvCubicHilbertCertificate[{name}] = basis=[lambda, zeta, U2, U3], gram={:?}, determinant_mod3={determinant}, boundary_frobenius={}, uniformizer_norms_mod_lambda5={norms:?}, norm_U2_mod27={}, zeta9_action={}", self.gram, self.boundary_action, self.unit_norm_mod27, self.cyclotomic_action)
    }
}

#[derive(Clone, Copy)]
pub(crate) enum LocalUnary {
    Inverse,
    OneMinus { source: usize },
}

fn valuation_range() -> Error {
    Error::plain("signed cubic local valuation exceeds the i128 range")
}

fn inconsistent() -> Error {
    Error::plain("inconsistent wild cubic Hilbert certificate or unit filtration")
}
