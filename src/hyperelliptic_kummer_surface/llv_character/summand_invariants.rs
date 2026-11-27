use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt::Write;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, MAX_SHOWN_BYTES};

use super::multiplication::{Factors, Group};
use super::orthogonal::DIMENSION;
use super::{
    arithmetic, decimal, guard, report_limit, validate, validate_length, weight_exceeds, weyl,
    work_limit, Magnitude, Oscillators, WEYL_LIMBS,
};

mod littlewood;

#[derive(Clone, Copy)]
pub(crate) enum HodgeFactors {
    Generic,
    Real(u128),
    Complex(u128),
}

#[derive(Clone, Copy)]
struct Blocks {
    rho: usize,
    number: usize,
    rank: usize,
    group: Group,
}

impl Blocks {
    fn new(rho: u128, kind: HodgeFactors) -> Result<Self> {
        let (degree, group) = match kind {
            HodgeFactors::Generic => {
                if !(1..=(DIMENSION - 2) as u128).contains(&rho) {
                    return Err(Error::plain("Hodge polarization rank must be in 1..=21"));
                }
                return Ok(Self {
                    rho: rho as usize,
                    number: 1,
                    rank: DIMENSION - rho as usize,
                    group: Group::Orthogonal,
                });
            }
            HodgeFactors::Real(degree) => (degree, Group::Orthogonal),
            HodgeFactors::Complex(degree) => (degree, Group::Unitary),
        };
        let factors = Factors::new(group, rho, degree)?;
        Ok(Self {
            rho: factors.rho,
            number: factors.degree,
            rank: factors.dimension,
            group,
        })
    }
}

pub(super) fn report(
    model: &mut Model,
    n: u128,
    d: i128,
    rho: u128,
    factors: HodgeFactors,
    shape: Option<Vec<u128>>,
) -> Result<String> {
    validate_length(n)?;
    let blocks = Blocks::new(rho, factors)?;
    let header = match factors {
        HodgeFactors::Generic => format!("{n}, {d}, {rho}"),
        HodgeFactors::Real(e) => format!("{n}, {d}, {rho}, real, {e}"),
        HodgeFactors::Complex(e) => format!("{n}, {d}, {rho}, cm, {e}"),
    };
    let charge = degree_charge(n, d);
    if let Some(shape) = shape {
        let shape = validate(shape)?;
        let irrep = if let Some(a) = charge {
            Character::new(blocks, a).dimension(model, &shape)?
        } else {
            Magnitude::default()
        };
        let mut contribution = Magnitude::default();
        if charge.is_some_and(|a| a <= n) && !irrep.is_zero() && !weight_exceeds(&shape, n) {
            let copies = if is_verbitsky(n, &shape) {
                1
            } else {
                Oscillators::new(model, n)?.multiplicity(model, &shape)?
            };
            contribution = irrep.multiply_u128(copies as u128).ok_or_else(arithmetic)?;
        }
        return Ok(format!(
            "rvHilbLLVHodgeInvariant[{header}]{shape:?} = irrep={} contribution={}",
            decimal(model, irrep)?,
            decimal(model, contribution)?
        ));
    }

    let mut result = format!("rvHilbLLVHodgeInvariantDecomposition[{header}] = components=[");
    let mut total = Magnitude::default();
    let mut verbitsky = Magnitude::default();
    let mut divisors = Magnitude::default();
    if let Some(a) = charge.filter(|&a| a <= n) {
        divisors = symmetric_dimension(model, blocks.rho, n - a)?;
        if a == n {
            let shape: Vec<_> = (n != 0).then_some(n).into_iter().collect();
            write!(result, "({shape:?}, 1, 1, 1)").unwrap();
            verbitsky = Magnitude::from_small(1);
            total = verbitsky;
        } else {
            let mut character = Character::new(blocks, a);
            let oscillators = Oscillators::new(model, n)?;
            let mut separator = "";
            oscillators.visit_summands(model, |model, shape, copies| {
                let irrep = character.dimension(model, shape)?;
                if irrep.is_zero() {
                    return Ok(());
                }
                let contribution = irrep.multiply_u128(copies as u128).ok_or_else(arithmetic)?;
                if is_verbitsky(n, shape) {
                    verbitsky = irrep;
                }
                total.add_multiple(contribution, 1).ok_or_else(arithmetic)?;
                write!(
                    result,
                    "{separator}({shape:?}, {copies}, {}, {})",
                    decimal(model, irrep)?,
                    decimal(model, contribution)?,
                )
                .unwrap();
                separator = ", ";
                check_report(&result)
            })?;
        }
    }
    let quadratic = subtract(verbitsky, divisors)?;
    let other = subtract(total, verbitsky)?;
    write!(
        result,
        "] total={} divisors={} quadraticContractions={} nonVerbitsky={}",
        decimal(model, total)?,
        decimal(model, divisors)?,
        decimal(model, quadratic)?,
        decimal(model, other)?,
    )
    .unwrap();
    check_report(&result)?;
    Ok(result)
}

fn check_report(result: &str) -> Result<()> {
    if result.len() > MAX_SHOWN_BYTES {
        Err(report_limit())
    } else {
        Ok(())
    }
}

fn is_verbitsky(n: u128, shape: &[u128]) -> bool {
    shape.len() <= 1 && shape.first().copied().unwrap_or(0) == n
}

fn degree_charge(n: u128, d: i128) -> Option<u128> {
    if d % 2 != 0 {
        return None;
    }
    let half = d / 2;
    if half < 0 {
        n.checked_add(half.unsigned_abs())
    } else {
        Some(n.abs_diff(half as u128))
    }
}

fn subtract(value: Magnitude, minus: Magnitude) -> Result<Magnitude> {
    if value.cmp(minus) == Ordering::Less {
        Err(Error::plain(
            "LLV Hodge invariant subtraction gives a negative dimension",
        ))
    } else {
        Ok(value.difference(minus))
    }
}

fn symmetric_dimension(model: &mut Model, variables: usize, degree: u128) -> Result<Magnitude> {
    let mut value = Magnitude::from_small(1);
    for i in 1..variables {
        let mut factor = Magnitude::from_small(degree);
        factor
            .add_multiple(Magnitude::from_small(i as u128), 1)
            .ok_or_else(arithmetic)?;
        value = value.checked_product(factor).ok_or_else(arithmetic)?;
        let (quotient, remainder) = value.divide(model, i as u128)?;
        if remainder != 0 {
            return Err(arithmetic());
        }
        value = quotient;
    }
    Ok(value)
}

// The trace term cancels the denominator from the two opposite charged
// lines of the standard representation. Each neutral block contributes
// one invariant quadratic form in symmetric powers of that block.
fn harmonic_symmetric(model: &mut Model, blocks: Blocks, excess: u128) -> Result<Magnitude> {
    if excess < 2 {
        return symmetric_dimension(model, blocks.rho, excess);
    }
    let terms = usize::try_from(excess / 2 + 1).map_err(|_| work_limit())?;
    model
        .spend(terms.saturating_mul((blocks.rho + blocks.number) * WEYL_LIMBS))
        .map_err(|_| work_limit())?;
    let mut total = Magnitude::default();
    for j in 0..terms {
        let products = symmetric_dimension(model, blocks.number, j as u128)?;
        let divisors = symmetric_dimension(model, blocks.rho, excess - 2 * j as u128)?;
        total
            .add_multiple(
                products.checked_product(divisors).ok_or_else(arithmetic)?,
                1,
            )
            .ok_or_else(arithmetic)?;
    }
    Ok(total)
}

type Shape = Vec<usize>;
type DiagramSizes = Vec<Vec<Shape>>;
type Memo = HashMap<Shape, Magnitude>;

struct Character {
    blocks: Blocks,
    charge: u128,
    diagrams: HashMap<Shape, Rc<DiagramSizes>>,
    irreps: Memo,
    neutral: Vec<Memo>,
    unitary: HashMap<Shape, u128>,
    cells: usize,
}

impl Character {
    fn new(blocks: Blocks, charge: u128) -> Self {
        Self {
            blocks,
            charge,
            diagrams: HashMap::new(),
            irreps: HashMap::new(),
            neutral: (0..=blocks.number).map(|_| HashMap::new()).collect(),
            unitary: HashMap::new(),
            cells: 0,
        }
    }

    fn dimension(&mut self, model: &mut Model, shape: &[u128]) -> Result<Magnitude> {
        let first = shape.first().copied().unwrap_or(0);
        let Some(excess) = first.checked_sub(self.charge) else {
            return Ok(Magnitude::default());
        };
        if shape.len() <= 1 {
            return harmonic_symmetric(model, self.blocks, excess);
        }
        let shape = shape
            .iter()
            .map(|&part| usize::try_from(part).map_err(|_| work_limit()))
            .collect::<Result<Vec<_>>>()?;
        self.orthogonal_irrep(model, &shape)
    }

    fn reserve(&mut self, cells: usize) -> Result<()> {
        self.cells = self.cells.saturating_add(cells);
        guard(self.cells)
    }

    fn subdiagrams(&mut self, model: &mut Model, shape: &[usize]) -> Result<Rc<DiagramSizes>> {
        if let Some(value) = self.diagrams.get(shape) {
            return Ok(Rc::clone(value));
        }
        let boxes = shape
            .iter()
            .try_fold(0usize, |n, &x| n.checked_add(x))
            .ok_or_else(work_limit)?;
        self.reserve(boxes.saturating_add(1).saturating_mul(8))?;
        let mut diagrams = vec![Vec::new(); boxes + 1];
        collect_subdiagrams(model, shape, &mut Vec::new(), 0, &mut |shape, size| {
            self.reserve(4 * shape.len() + 8)?;
            diagrams[size].push(shape.to_vec());
            Ok(())
        })?;
        let diagrams = Rc::new(diagrams);
        self.diagrams.insert(shape.to_vec(), Rc::clone(&diagrams));
        Ok(diagrams)
    }

    fn orthogonal_irrep(&mut self, model: &mut Model, shape: &[usize]) -> Result<Magnitude> {
        if shape.len() <= 1 {
            let shape = shape.iter().map(|&part| part as u128).collect::<Vec<_>>();
            return self.dimension(model, &shape);
        }
        if (shape[0] as u128) < self.charge {
            return Ok(Magnitude::default());
        }
        if let Some(&value) = self.irreps.get(shape) {
            return Ok(value);
        }
        let gl_value = self.charged_schur(model, shape)?;
        let mut correction = Magnitude::default();
        let diagrams = self.subdiagrams(model, shape)?;
        let boxes = diagrams.len() - 1;
        // Triangular Littlewood restriction from GL(25) to SO(25).
        for even_size in (2..=boxes).step_by(2) {
            for even in diagrams[even_size]
                .iter()
                .filter(|x| x.iter().all(|v| v.is_multiple_of(2)))
            {
                for inner in &diagrams[boxes - even_size] {
                    let lr = littlewood::coefficient(model, shape, inner, even)?;
                    if lr != 0 {
                        correction
                            .add_multiple(self.orthogonal_irrep(model, inner)?, lr)
                            .ok_or_else(arithmetic)?;
                    }
                }
            }
        }
        let count = subtract(gl_value, correction)?;
        self.reserve(4 * shape.len() + WEYL_LIMBS + 8)?;
        self.irreps.insert(shape.to_vec(), count);
        Ok(count)
    }

    fn charged_schur(&mut self, model: &mut Model, shape: &[usize]) -> Result<Magnitude> {
        let mut sum = Magnitude::default();
        let diagrams = self.subdiagrams(model, shape)?;
        let boxes = diagrams.len() - 1;
        for (size, candidates) in diagrams.iter().enumerate() {
            let plus = (boxes - size) as u128;
            let Some(minus) = plus.checked_sub(self.charge).filter(|&m| m <= size as u128) else {
                continue;
            };
            for first in candidates
                .iter()
                .filter(|inner| littlewood::horizontal(shape, inner))
            {
                model.spend(1)?;
                let rest = self.subdiagrams(model, first)?;
                for neutral in rest[size - minus as usize]
                    .iter()
                    .filter(|inner| littlewood::horizontal(first, inner))
                {
                    model.spend(WEYL_LIMBS)?;
                    let coeff = self.neutral_schur(model, neutral, self.blocks.number)?;
                    sum.add_multiple(coeff, 1).ok_or_else(arithmetic)?;
                }
            }
        }
        Ok(sum)
    }

    fn neutral_schur(
        &mut self,
        model: &mut Model,
        shape: &[usize],
        factors: usize,
    ) -> Result<Magnitude> {
        if shape.is_empty() {
            return Ok(Magnitude::from_small(1));
        }
        if let Some(&value) = self.neutral[factors].get(shape) {
            return Ok(value);
        }
        let value = if factors == 0 {
            weyl::gl_dimension(model, shape, self.blocks.rho)?
        } else {
            let mut sum = Magnitude::default();
            let diagrams = self.subdiagrams(model, shape)?;
            let boxes = diagrams.len() - 1;
            for (size, right) in diagrams.iter().enumerate() {
                for block in right {
                    let coeff = self.block_schur(model, block)?;
                    if coeff == 0 {
                        continue;
                    }
                    for inner in &diagrams[boxes - size] {
                        let lr = littlewood::coefficient(model, shape, inner, block)?;
                        if lr != 0 {
                            let value = self.neutral_schur(model, inner, factors - 1)?;
                            let product = value.multiply_u128(coeff).ok_or_else(arithmetic)?;
                            sum.add_multiple(product, lr).ok_or_else(arithmetic)?;
                        }
                    }
                }
            }
            sum
        };
        self.reserve(4 * shape.len() + WEYL_LIMBS + 8)?;
        self.neutral[factors].insert(shape.to_vec(), value);
        Ok(value)
    }

    fn block_schur(&mut self, model: &mut Model, shape: &[usize]) -> Result<u128> {
        model.spend(1)?;
        match self.blocks.group {
            Group::Orthogonal => Ok(u128::from(
                shape.len() <= self.blocks.rank && shape.iter().all(|x| x.is_multiple_of(2))
                    || shape.len() == self.blocks.rank
                        && shape.iter().all(|x| !x.is_multiple_of(2)),
            )),
            Group::Unitary => {
                if let Some(&value) = self.unitary.get(shape) {
                    return Ok(value);
                }
                let boxes = shape.iter().sum::<usize>();
                if !boxes.is_multiple_of(2) {
                    return Ok(0);
                }
                let diagrams = self.subdiagrams(model, shape)?;
                let mut value = 0u128;
                for half in diagrams[boxes / 2]
                    .iter()
                    .filter(|x| x.len() <= self.blocks.rank)
                {
                    value = value
                        .checked_add(u128::from(littlewood::coefficient(
                            model, shape, half, half,
                        )?))
                        .ok_or_else(arithmetic)?;
                }
                self.reserve(4 * shape.len() + 12)?;
                self.unitary.insert(shape.to_vec(), value);
                Ok(value)
            }
        }
    }
}

fn collect_subdiagrams(
    model: &mut Model,
    outer: &[usize],
    shape: &mut Shape,
    size: usize,
    emit: &mut impl FnMut(&[usize], usize) -> Result<()>,
) -> Result<()> {
    model.spend(1)?;
    emit(shape, size)?;
    if let Some(&row) = outer.get(shape.len()) {
        let bound = shape.last().copied().unwrap_or(row).min(row);
        for part in 1..=bound {
            shape.push(part);
            collect_subdiagrams(model, outer, shape, size + part, emit)?;
            shape.pop();
        }
    }
    Ok(())
}
