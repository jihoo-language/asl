use std::cell::OnceCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use super::nilpotent::Limit;
use super::{
    add_scaled, column_integer, dot, eliminate, pairing_row, polarized_basis, raw, realize, Class,
    Column, Cusp, Lattice, COLUMN_CELLS, RAW_RANK,
};
use crate::error::{Error, Result};
use crate::hyperelliptic_kummer_surface::rational_sparse::{
    decimal, Coefficients, Domain, Integer, Magnitude, WORK_LIMBS,
};
use crate::integer_smith::{self, Matrix, Track};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod exterior;
pub(crate) mod projector;

const RANK: usize = 7;
const DIMENSION: usize = 1 << RANK;
const EVEN_RANK: usize = DIMENSION / 2;
type Terms = BTreeMap<u32, Integer>;

pub(crate) enum Expression {
    Basis(Label, u32),
    PolarizedBasis(Label, usize),
    Vector(Label, Box<ParsedExpr>),
    Multiply(Box<[ParsedExpr; 2]>),
    Scale(i128, Box<ParsedExpr>),
    Log(Label),
    Exp(Label, i128),
}

pub(crate) enum Command {
    Data(Label),
    Coefficient(u32, Box<ParsedExpr>),
}

pub(crate) enum Og6Clifford {}

impl Domain for Og6Clifford {
    const NAME: &'static str = "OG6 integral Clifford";
}

pub(crate) type Element = Coefficients<Og6Clifford>;

impl Element {
    pub(crate) fn blade(mask: u32) -> Result<Self> {
        check_mask(mask)?;
        Ok(Self::singleton(mask, 1))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        self.charge_storage(model)?;
        Ok(Value::OG6Clifford(Rc::new(self)))
    }

    fn integral(model: &mut Model, terms: &Terms) -> Result<Self> {
        Self::normalized(
            model,
            terms.iter().map(|(&mask, value)| (mask, value)),
            Magnitude::from_small(1),
        )
    }
}

pub(crate) fn check_mask(mask: u32) -> Result<()> {
    if (mask as usize) >= DIMENSION {
        return Err(Error::plain(
            "OG6 Clifford blade mask must be between 0 and 127",
        ));
    }
    Ok(())
}

pub(crate) struct Algebra {
    basis: Vec<Column>,
    dual: Vec<Column>,
    gram: [[Integer; RANK]; RANK],
    projector: OnceCell<Box<projector::Marking>>,
}

type RegularColumns = Vec<Vec<(usize, Integer)>>;

pub(super) struct Action {
    columns: RegularColumns,
    factors: OnceCell<Vec<Magnitude>>,
}

impl Cusp {
    pub(crate) fn spin_algebra(&self, model: &mut Model, lattice: &Lattice) -> Result<&Algebra> {
        if let Some(algebra) = self.spin.get() {
            return Ok(algebra);
        }
        let algebra = Algebra::new(model, lattice, &self.polarization)?;
        model.allocate(2 * RANK * RAW_RANK * COLUMN_CELLS + RANK * RANK * (WORK_LIMBS + 1))?;
        Ok(self.spin.get_or_init(|| Box::new(algebra)))
    }
}

impl Algebra {
    fn new(model: &mut Model, lattice: &Lattice, h: &Class) -> Result<Self> {
        let basis = polarized_basis(model, lattice, &raw(h))?;
        let mut forms = vec![[column_integer(0); RAW_RANK]; RAW_RANK];
        for (i, form) in forms.iter_mut().enumerate() {
            form[i] = column_integer(1);
        }
        let mut dual = Vec::with_capacity(RANK);
        for column in &basis {
            let row = forms
                .iter()
                .map(|form| dot(model, form, column))
                .collect::<Result<_>>()?;
            let selected = eliminate(
                model,
                &mut forms,
                row,
                "polarized OG6 lattice is not saturated",
            )?;
            for prior in &mut dual {
                let coefficient = dot(model, prior, column)?.negate();
                add_scaled(model, prior, &selected, coefficient)?;
            }
            dual.push(selected);
        }
        let zero = Integer::from_i128(0);
        let mut gram = [[zero; RANK]; RANK];
        for (i, x) in basis.iter().enumerate() {
            let row = pairing_row(x)?;
            for (j, y) in basis[..=i].iter().enumerate() {
                let value: Integer = dot(model, &row, y)?.resize()?;
                let value = if i == j {
                    // Diagonal entries in the rewrite rule are q/2.
                    value
                        .divide_exact(model, Integer::from_i128(2))?
                        .ok_or_else(|| Error::plain("polarized OG6 Clifford form must be even"))?
                } else {
                    value
                };
                gram[i][j] = value;
                gram[j][i] = value;
            }
        }
        Ok(Self {
            basis,
            dual,
            gram,
            projector: OnceCell::new(),
        })
    }

    pub(crate) fn basis(&self, index: usize) -> Result<Class> {
        let column = self
            .basis
            .get(index)
            .ok_or_else(|| Error::plain("rvOG6PolarizedBasis index must be between 0 and 6"))?;
        realize(column)
    }

    pub(crate) fn vector(&self, model: &mut Model, cusp: &Cusp, x: &Class) -> Result<Element> {
        if !x.wide_pair(model, &cusp.polarization)?.is_zero() {
            return Err(Error::plain(
                "OG6 Clifford vector must be perpendicular to the polarization",
            ));
        }
        let terms = self.vector_terms(model, x)?;
        Element::integral(model, &terms)
    }

    fn vector_terms(&self, model: &mut Model, x: &Class) -> Result<Terms> {
        let x = raw(x);
        let mut terms = Terms::new();
        for (i, dual) in self.dual.iter().enumerate() {
            let value = dot(model, dual, &x)?.resize()?;
            Element::accumulate(&mut terms, 1 << i, value)?;
        }
        Ok(terms)
    }

    pub(crate) fn multiply(
        &self,
        model: &mut Model,
        left: &Element,
        right: &Element,
    ) -> Result<Element> {
        let left = widen(left)?;
        let right = widen(right)?;
        let result = self.product(model, &left, &right)?;
        Element::integral(model, &result)
    }

    fn product(&self, model: &mut Model, left: &Terms, right: &Terms) -> Result<Terms> {
        let mut result = Terms::new();
        for (&mask, &scalar) in left {
            let mut terms = right.clone();
            let mut remaining = mask;
            // In a sorted PBW word, the highest generator acts first.
            while remaining != 0 && !terms.is_empty() {
                let index = remaining.ilog2() as usize;
                remaining ^= 1 << index;
                terms = self.generator_left(model, index, &terms)?;
            }
            model.spend((1 + terms.len()) * WORK_LIMBS)?;
            for (key, value) in terms {
                Element::accumulate(&mut result, key, scalar.multiply_wide(value)?)?;
            }
        }
        Ok(result)
    }

    fn generator_left(&self, model: &mut Model, index: usize, terms: &Terms) -> Result<Terms> {
        model.spend(terms.len() * (RANK + 1) * WORK_LIMBS)?;
        let bit = 1 << index;
        let mut result = Terms::new();
        for (&mask, &scalar) in terms {
            let mut coefficient = scalar;
            let mut smaller = mask & (bit - 1);
            while smaller != 0 {
                let j = smaller.trailing_zeros() as usize;
                let gram = self.gram[index][j];
                if !gram.is_zero() {
                    let contraction = coefficient.multiply_wide(gram)?;
                    Element::accumulate(&mut result, mask ^ (1 << j), contraction)?;
                }
                coefficient = coefficient.negate();
                smaller &= smaller - 1;
            }
            if mask & bit != 0 {
                coefficient = coefficient.multiply_wide(self.gram[index][index])?;
            }
            Element::accumulate(&mut result, mask ^ bit, coefficient)?;
        }
        Ok(result)
    }

    fn lift(&self, model: &mut Model, limit: &Limit) -> Result<Terms> {
        let (left, right, coefficient) = limit.spin_factors();
        let left = self.vector_terms(model, left)?;
        let right = self.vector_terms(model, right)?;
        let mut lift = self.product(model, &left, &right)?;
        for value in lift.values_mut() {
            *value = value.multiply(coefficient)?;
        }
        Ok(lift)
    }

    pub(crate) fn spin(
        &self,
        model: &mut Model,
        limit: &Limit,
        exponent: Option<i128>,
    ) -> Result<Element> {
        let mut lift = if exponent == Some(0) {
            Terms::new()
        } else {
            self.lift(model, limit)?
        };
        if let Some(m) = exponent {
            for value in lift.values_mut() {
                *value = value.multiply(m)?;
            }
            Element::accumulate(&mut lift, 0, Integer::from_i128(1))?;
        }
        Element::integral(model, &lift)
    }

    pub(crate) fn report(&self, model: &mut Model, name: &str, limit: &Limit) -> Result<String> {
        let factors = self.action(model, limit)?.smith(model)?;
        let rank = factors.len();
        let middle = EVEN_RANK
            .checked_sub(2 * rank)
            .ok_or_else(|| Error::plain("OG6 spin logarithm has incompatible square-zero rank"))?;
        let factors = factors
            .iter()
            .copied()
            .map(decimal)
            .collect::<Vec<_>>()
            .join(", ");
        Ok(format!(
            "rvOG6SpinData[{name}] = rank=64, abelian_dimension=32, log_rank={rank}, graded=[{rank}, {middle}, {rank}], smith=[{factors}]"
        ))
    }

    fn action<'a>(&self, model: &mut Model, limit: &'a Limit) -> Result<&'a Action> {
        if let Some(action) = limit.spin_action.get() {
            return Ok(action);
        }
        let lift = self.lift(model, limit)?;
        if !self.product(model, &lift, &lift)?.is_empty() {
            return Err(Error::plain("OG6 spin logarithm must square to zero"));
        }
        let columns = self.regular_columns(model, &lift)?;
        let entries: usize = columns.iter().map(Vec::len).sum();
        model.allocate(3 * EVEN_RANK + entries * (WORK_LIMBS + 3))?;
        Ok(limit.spin_action.get_or_init(|| {
            Box::new(Action {
                columns,
                factors: OnceCell::new(),
            })
        }))
    }

    fn regular_columns(&self, model: &mut Model, left: &Terms) -> Result<RegularColumns> {
        let mut columns = Vec::with_capacity(EVEN_RANK);
        for index in 0..EVEN_RANK {
            let right = Terms::from([(even_mask(index), Integer::from_i128(1))]);
            columns.push(
                self.product(model, left, &right)?
                    .into_iter()
                    .map(|(mask, coefficient)| (even_index(mask), coefficient))
                    .collect(),
            );
        }
        Ok(columns)
    }
}

impl Action {
    fn smith(&self, model: &mut Model) -> Result<&[Magnitude]> {
        if let Some(factors) = self.factors.get() {
            return Ok(factors);
        }
        let factors = self.reduce(model)?;
        model.allocate(1 + factors.len() * WORK_LIMBS)?;
        Ok(self.factors.get_or_init(|| factors))
    }

    fn reduce(&self, model: &mut Model) -> Result<Vec<Magnitude>> {
        let mut content = Magnitude::from_small(0);
        for (_, entry) in self.columns.iter().flatten() {
            if content.small() == Some(1) {
                break;
            }
            content = content.gcd(model, entry.magnitude)?;
        }
        if content.is_zero() {
            return Ok(Vec::new());
        }
        // Extract common content in wider arithmetic before invoking the
        // existing exact Smith routine on the actual ordered PBW matrix.
        let mut matrix = Matrix::new(model, EVEN_RANK, EVEN_RANK)?;
        for (j, column) in self.columns.iter().enumerate() {
            for &(i, entry) in column {
                let reduced = entry
                    .magnitude
                    .divide_exact(model, content)?
                    .expect("matrix content divides each entry");
                let coefficient = reduced.signed(entry.negative).ok_or_else(|| {
                    Error::plain(
                        "normalized OG6 spin matrix exceeds the signed 128-bit coefficient range",
                    )
                })?;
                matrix.write(model, i, j, coefficient)?;
            }
        }
        let reduced = integer_smith::reduce(model, matrix, Track::None, Track::None)?;
        reduced
            .factors
            .into_iter()
            .map(|factor| {
                content
                    .multiply_u128(factor.unsigned_abs())
                    .ok_or_else(Element::arithmetic)
            })
            .collect()
    }
}

// In each consecutive pair of PBW masks exactly one has even parity.
fn even_mask(index: usize) -> u32 {
    ((index as u32) << 1) | (index.count_ones() & 1)
}

fn even_index(mask: u32) -> usize {
    debug_assert!(mask.count_ones().is_multiple_of(2));
    (mask >> 1) as usize
}

fn widen(element: &Element) -> Result<Terms> {
    element
        .terms()
        .iter()
        .map(|&(mask, scalar)| Ok((mask, scalar.resize()?)))
        .collect()
}
