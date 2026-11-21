use crate::bp::{self, Shape};
use crate::bp_polynomial::{Monomial, Polynomial};
use crate::bp_rel;
use crate::bp_smith::{self, Matrix};
use crate::complex::MAX_PIVOT_CELLS;
use crate::error::{Error, Result};
use crate::model::Model;

const MAX_BASIS: usize = 32_768;

pub(crate) struct Group {
    pub(crate) free: usize,
    pub(crate) torsion: Vec<i128>,
}

#[derive(Clone, Copy)]
pub(crate) struct Domain {
    pub(crate) shape: Shape,
    lower: Option<usize>,
}

impl Domain {
    pub(crate) fn absolute(shape: Shape) -> Self {
        Self { shape, lower: None }
    }

    pub(crate) fn relative(shape: Shape, lower: usize) -> Self {
        Self {
            shape,
            lower: Some(lower),
        }
    }

    fn at_arrows(self, arrows: usize) -> Self {
        Self {
            shape: Shape {
                arrows,
                ..self.shape
            },
            ..self
        }
    }

    pub(crate) fn previous(self) -> Result<Self> {
        Ok(Self {
            shape: previous_shape(self.shape)?,
            ..self
        })
    }

    fn retains(self, monomial: &Monomial) -> bool {
        self.lower
            .is_none_or(|lower| bp_rel::retains(monomial, lower, self.shape.height))
    }

    pub(crate) fn project(self, model: &mut Model, input: &Polynomial) -> Result<Polynomial> {
        if self.lower.is_none() {
            model.spend(input.cells())?;
            return Ok(input.clone());
        }
        model.spend(input.cells())?;
        Polynomial::from_terms(
            model,
            input.terms().filter(|(monomial, _)| self.retains(monomial)),
        )
    }

    pub(crate) fn differential(self, model: &mut Model, input: &Polynomial) -> Result<Polynomial> {
        let image = bp::differential(model, self.shape, input)?;
        if self.lower.is_some() {
            self.project(model, &image)
        } else {
            Ok(image)
        }
    }

    pub(crate) fn group(self, model: &mut Model, t: usize) -> Result<Group> {
        let shape = self.shape;
        if t == 0 || shape.arrows == 0 || !self.has_monomials(t) {
            return Ok(Group {
                free: usize::from(self.lower.is_none() && t == 0 && shape.arrows == 0),
                torsion: Vec::new(),
            });
        }

        // Over Z[1/p] the right Hazewinkel log is monic triangular in
        // the arrow coordinates. Inserting the origin contracts the pair
        // groupoid in positive degree and preserves every lower height.
        // In degree zero the height inclusions are identities. Therefore
        // absolute and relative torsion are both p-primary.
        let previous = self.previous()?;
        let (source, rank) = previous.incoming_basis(model, t)?;
        if rank == 0 {
            return Ok(Group {
                free: 0,
                torsion: Vec::new(),
            });
        }
        let target = self.basis(model, t)?;
        let incoming = previous.matrix(model, &source, &target)?;
        let torsion = bp_smith::primary_factors(model, &incoming, shape.prime, rank)?
            .into_iter()
            .filter(|&factor| factor > 1)
            .collect();
        Ok(Group { free: 0, torsion })
    }

    pub(crate) fn incoming_basis(
        self,
        model: &mut Model,
        t: usize,
    ) -> Result<(Vec<Monomial>, usize)> {
        let mut rank = 0usize;
        for arrows in 0..self.shape.arrows {
            let dimension = self.at_arrows(arrows).basis(model, t)?.len();
            rank = exact_rank(dimension, rank)?;
        }
        let source = self.basis(model, t)?;
        rank = exact_rank(source.len(), rank)?;
        Ok((source, rank))
    }
}

pub(crate) fn group(model: &mut Model, shape: Shape, t: usize) -> Result<Group> {
    Domain::absolute(shape).group(model, t)
}

fn exact_rank(dimension: usize, incoming: usize) -> Result<usize> {
    dimension
        .checked_sub(incoming)
        .ok_or_else(|| Error::plain("BP rational ranks are incompatible with exactness"))
}

pub(crate) fn previous_shape(shape: Shape) -> Result<Shape> {
    let arrows = shape
        .arrows
        .checked_sub(1)
        .ok_or_else(|| Error::plain("bpNull requires positive BP cobar arity"))?;
    Ok(Shape { arrows, ..shape })
}

pub(crate) fn order(
    model: &mut Model,
    shape: Shape,
    t: usize,
    cycle: &Polynomial,
) -> Result<Option<i128>> {
    shape.require_integral("bporder")?;
    Domain::absolute(shape).order(model, t, cycle)
}

pub(crate) fn null(
    model: &mut Model,
    shape: Shape,
    t: usize,
    cycle: &Polynomial,
) -> Result<Polynomial> {
    Domain::absolute(shape).null(model, t, cycle)
}

impl Domain {
    pub(crate) fn order(
        self,
        model: &mut Model,
        t: usize,
        cycle: &Polynomial,
    ) -> Result<Option<i128>> {
        if cycle.is_zero() {
            return Ok(Some(1));
        }
        self.validate_cycle(model, t, cycle)?;
        if self.shape.arrows == 0 {
            return Ok(None);
        }
        let ClosedSystem {
            incoming,
            target,
            rank,
            ..
        } = self.closed_system(model, t, cycle)?;
        finite_order(model, &incoming, &target, self.shape.prime, rank).map(Some)
    }

    pub(crate) fn null(
        self,
        model: &mut Model,
        t: usize,
        cycle: &Polynomial,
    ) -> Result<Polynomial> {
        // Without an order restriction every torsion class has a witness.
        Ok(self.integral_witness(model, t, cycle, false)?.unwrap())
    }

    pub(crate) fn primitive(
        self,
        model: &mut Model,
        t: usize,
        cycle: &Polynomial,
    ) -> Result<Option<Polynomial>> {
        self.integral_witness(model, t, cycle, true)
    }

    fn integral_witness(
        self,
        model: &mut Model,
        t: usize,
        cycle: &Polynomial,
        require_boundary: bool,
    ) -> Result<Option<Polynomial>> {
        let previous = self.previous()?;
        if cycle.is_zero() {
            return Ok(Some(Polynomial::default()));
        }
        self.validate_cycle(model, t, cycle)?;
        let ClosedSystem {
            basis,
            incoming,
            target,
            rank,
        } = self.closed_system(model, t, cycle)?;
        let order = finite_order(model, &incoming, &target, self.shape.prime, rank)?;
        if require_boundary && order != 1 {
            return Ok(None);
        }
        drop(incoming);
        let mut multiple = cycle.clone();
        multiple.scale(model, order)?;
        let (mut numerator, mut power) = bp::contracting_numerator(model, self.shape, &multiple)?;
        if self.lower.is_some() {
            numerator = self.project(model, &numerator)?;
            let prime = u128::from(self.shape.prime);
            while power > 0 && numerator.divide_if_exact(model, prime)? {
                power -= 1;
            }
        }
        let witness = if power == 0 {
            numerator
        } else {
            previous.integral_lift(model, t, &basis, &numerator, power)?
        };
        let image = previous.differential(model, &witness)?;
        if image.compare(model, &multiple)? != std::cmp::Ordering::Equal {
            return Err(Error::plain(
                "BP integral witness does not have the required differential",
            ));
        }
        Ok(Some(witness))
    }

    fn integral_lift(
        self,
        model: &mut Model,
        t: usize,
        basis: &[Monomial],
        numerator: &Polynomial,
        power: u32,
    ) -> Result<Polynomial> {
        // For the rational lift U/M choose a closed Z congruent to U mod M.
        // The exponent E of the previous cohomology kills every class, so
        // Z can be found as d(y)/E by solving d(y) = E U modulo E M.
        let preceding = self.previous()?;
        let (source, rank) = preceding.incoming_basis(model, t)?;
        let matrix = preceding.matrix(model, &source, basis)?;
        let exponent = bp_smith::primary_factors(model, &matrix, self.shape.prime, rank)?
            .last()
            .copied()
            .unwrap_or(1) as u64;
        let denominator = u64::from(self.shape.prime)
            .checked_pow(power)
            .ok_or_else(precision_overflow)?;
        let modulus = denominator
            .checked_mul(exponent)
            .ok_or_else(precision_overflow)?;
        let mut target = coordinates(model, numerator, basis)?;
        let correction = bp_smith::modular_preimage(
            model,
            &matrix,
            self.shape.prime,
            modulus,
            &target,
            exponent,
            rank,
        )?;
        let image = matrix.apply(model, &correction)?;
        let exponent = i128::from(exponent);
        let denominator = i128::from(denominator);
        model.spend(target.len())?;
        for (coefficient, image) in target.iter_mut().zip(image) {
            if image % exponent != 0 {
                return Err(lift_divisibility());
            }
            let reduced = coefficient
                .checked_sub(image / exponent)
                .ok_or_else(bp_smith::witness_overflow)?;
            if reduced % denominator != 0 {
                return Err(lift_divisibility());
            }
            *coefficient = reduced / denominator;
        }
        Polynomial::from_terms(model, basis.iter().zip(target))
    }
}

fn precision_overflow() -> Error {
    Error::plain("BP integral witness precision exceeds the u64 modulus limit")
}

fn lift_divisibility() -> Error {
    Error::plain("BP integral witness correction is not exactly divisible")
}

struct ClosedSystem {
    basis: Vec<Monomial>,
    incoming: Matrix,
    target: Vec<i128>,
    rank: usize,
}

impl Domain {
    fn closed_system(
        self,
        model: &mut Model,
        t: usize,
        cycle: &Polynomial,
    ) -> Result<ClosedSystem> {
        let current = self.basis(model, t)?;
        let previous = self.previous()?;
        let (preceding, rank) = previous.incoming_basis(model, t)?;
        let target = coordinates(model, cycle, &current)?;
        let incoming = previous.matrix(model, &preceding, &current)?;
        Ok(ClosedSystem {
            basis: preceding,
            incoming,
            target,
            rank,
        })
    }
}

fn finite_order(
    model: &mut Model,
    incoming: &Matrix,
    target: &[i128],
    prime: u32,
    rank: usize,
) -> Result<i128> {
    let factors = bp_smith::primary_factors(model, incoming, prime, rank)?;
    let enlarged = incoming.augmented(model, target)?;
    let augmented = bp_smith::primary_factors(model, &enlarged, prime, rank)?;
    // Both lattices have the same rational span. Their Smith indices in
    // that saturation have ratio equal to the order of the adjoined class.
    let prime = i128::from(prime);
    model.spend(factors.len() + augmented.len())?;
    let index = |factors: Vec<i128>| {
        factors
            .into_iter()
            .map(|factor| factor.ilog(prime))
            .sum::<u32>()
    };
    let exponent = index(factors)
        .checked_sub(index(augmented))
        .ok_or_else(|| {
            Error::plain("BP Smith indices disagree with the incoming lattice inclusion")
        })?;
    prime.checked_pow(exponent).ok_or_else(order_overflow)
}

fn order_overflow() -> Error {
    Error::plain("BP integral class order or its multiple exceeds the signed i128 range")
}

pub(crate) fn check_normalized(
    model: &mut Model,
    shape: Shape,
    polynomial: &Polynomial,
) -> Result<()> {
    for (monomial, _) in polynomial.terms() {
        model.spend(1 + monomial.len())?;
        normalized_monomial(shape, monomial)?;
    }
    Ok(())
}

fn normalized_monomial(shape: Shape, monomial: &Monomial) -> Result<()> {
    let mut slots = 0u32;
    for &(variable, _) in monomial {
        slots |= 1 << (variable / shape.height);
    }
    require_arrow_slots(shape, slots)
}

fn require_arrow_slots(shape: Shape, slots: u32) -> Result<()> {
    let arrows = (1u32 << (shape.arrows + 1)) - 2;
    if slots & arrows != arrows {
        return Err(Error::plain(
            "BP Ext requires a normalized cochain using every arrow slot",
        ));
    }
    Ok(())
}

pub(crate) fn validate_homogeneous(
    model: &mut Model,
    shape: Shape,
    t: usize,
    cycle: &Polynomial,
) -> Result<()> {
    let weight = unit_weight(shape.prime, t).ok_or_else(wrong_weight)?;
    let weights = weights(shape, weight);
    for (monomial, _) in cycle.terms() {
        model.spend(1 + monomial.len())?;
        let mut remaining = weight;
        let mut slots = 0u32;
        for &(variable, exponent) in monomial {
            slots |= 1 << (variable / shape.height);
            let degree = weights
                .get(variable % shape.height)
                .and_then(|degree| degree.checked_mul(exponent))
                .ok_or_else(wrong_weight)?;
            remaining = remaining.checked_sub(degree).ok_or_else(wrong_weight)?;
        }
        require_arrow_slots(shape, slots)?;
        if remaining != 0 {
            return Err(wrong_weight());
        }
    }
    Ok(())
}

impl Domain {
    pub(crate) fn validate_cycle(
        self,
        model: &mut Model,
        t: usize,
        cycle: &Polynomial,
    ) -> Result<()> {
        let shape = self.shape;
        validate_homogeneous(model, shape, t, cycle)?;
        if self.at_arrows(shape.arrows + 1).has_monomials(t)
            && !self.differential(model, cycle)?.is_zero()
        {
            return Err(Error::plain("BP Ext requires a closed cochain"));
        }
        Ok(())
    }

    pub(crate) fn has_monomials(self, t: usize) -> bool {
        if t == 0 {
            return self.lower.is_none() && self.shape.arrows == 0;
        }
        if self.shape.height == 0 {
            return false;
        }
        let Some(weight) = unit_weight(self.shape.prime, t) else {
            return false;
        };
        let minimum = if let Some(lower) = self.lower {
            let Some(&high) = weights(self.shape, weight).get(lower) else {
                return false;
            };
            high.saturating_add(self.shape.arrows.saturating_sub(1))
        } else {
            self.shape.arrows
        };
        weight >= minimum
    }
}

fn wrong_weight() -> Error {
    Error::plain("BP Ext cochain does not have the requested homogeneous internal degree")
}

fn unit_weight(prime: u32, t: usize) -> Option<usize> {
    let first_degree = 2 * (u128::from(prime) - 1);
    let t = t as u128;
    t.is_multiple_of(first_degree)
        .then(|| (t / first_degree) as usize)
}

fn weights(shape: Shape, weight: usize) -> Vec<usize> {
    let mut result = Vec::new();
    let mut next = 1u128;
    for _ in 0..shape.height {
        if next > weight as u128 {
            break;
        }
        result.push(next as usize);
        next = next * u128::from(shape.prime) + 1;
    }
    result
}

pub(crate) fn basis(model: &mut Model, shape: Shape, t: usize) -> Result<Vec<Monomial>> {
    Domain::absolute(shape).basis(model, t)
}

impl Domain {
    pub(crate) fn basis(self, model: &mut Model, t: usize) -> Result<Vec<Monomial>> {
        let shape = self.shape;
        if !self.has_monomials(t) {
            return Ok(Vec::new());
        }
        if t == 0 {
            return Ok(vec![Box::new([])]);
        }
        let weight = unit_weight(shape.prime, t).unwrap();
        let mut search = Monomials {
            shape,
            weights: weights(shape, weight),
            lower: self.lower,
            high_terms: 0,
            current: Vec::new(),
            collected: Vec::new(),
            cells: 0,
        };
        search.enumerate(model, 0, search.weights.len() - 1, weight, false)?;
        model.spend(search.cells * (1 + search.collected.len().max(1).ilog2() as usize))?;
        search.collected.sort_unstable();
        Ok(search.collected)
    }
}

struct Monomials {
    shape: Shape,
    weights: Vec<usize>,
    lower: Option<usize>,
    high_terms: usize,
    current: Vec<(usize, usize)>,
    collected: Vec<Monomial>,
    cells: usize,
}

impl Monomials {
    fn enumerate(
        &mut self,
        model: &mut Model,
        slot: usize,
        variable: usize,
        remaining: usize,
        used: bool,
    ) -> Result<()> {
        model.spend(1)?;
        if let Some(lower) = self.lower {
            if self.high_terms == 0
                && (remaining < self.weights[lower]
                    || (slot == self.shape.arrows && variable < lower))
            {
                return Ok(());
            }
        }
        let reserved = self.shape.arrows - slot;
        let Some(available) = remaining.checked_sub(reserved) else {
            return Ok(());
        };
        let degree = self.weights[variable];
        let minimum = usize::from(slot != 0 && variable == 0 && !used);
        let maximum = available / degree;
        // The last variable is v1 or t1 and has unit weight, so its
        // exponent is forced in the final slot.
        if reserved == 0 && variable == 0 {
            if maximum >= minimum {
                self.with_power(slot, variable, maximum);
                self.collect(model)?;
                self.without_power(variable, maximum);
            }
            return Ok(());
        }
        for exponent in minimum..=maximum {
            self.with_power(slot, variable, exponent);
            let remaining = remaining - exponent * degree;
            if let Some(previous) = variable.checked_sub(1) {
                self.enumerate(model, slot, previous, remaining, used || exponent != 0)?;
            } else {
                self.enumerate(model, slot + 1, self.weights.len() - 1, remaining, false)?;
            }
            self.without_power(variable, exponent);
        }
        Ok(())
    }

    fn with_power(&mut self, slot: usize, variable: usize, exponent: usize) {
        if exponent != 0 {
            self.current
                .push((slot * self.shape.height + variable, exponent));
            self.high_terms += usize::from(self.lower.is_some_and(|lower| variable >= lower));
        }
    }

    fn without_power(&mut self, variable: usize, exponent: usize) {
        if exponent != 0 {
            self.high_terms -= usize::from(self.lower.is_some_and(|lower| variable >= lower));
            self.current.pop();
        }
    }

    fn collect(&mut self, model: &mut Model) -> Result<()> {
        if self.lower.is_some() && self.high_terms == 0 {
            return Ok(());
        }
        self.cells += 1 + 2 * self.current.len();
        if self.collected.len() == MAX_BASIS || self.cells > MAX_PIVOT_CELLS {
            return Err(Error::plain(
                "BP normalized basis exceeds the storage limit",
            ));
        }
        model.spend(1 + self.current.len() * (1 + self.current.len().max(1).ilog2() as usize))?;
        let mut monomial = self.current.clone().into_boxed_slice();
        monomial.sort_unstable();
        self.collected.push(monomial);
        Ok(())
    }
}

pub(crate) fn matrix(
    model: &mut Model,
    shape: Shape,
    source: &[Monomial],
    target: &[Monomial],
) -> Result<Matrix> {
    Domain::absolute(shape).matrix(model, source, target)
}

impl Domain {
    pub(crate) fn matrix(
        self,
        model: &mut Model,
        source: &[Monomial],
        target: &[Monomial],
    ) -> Result<Matrix> {
        restricted_matrix(
            model,
            self.shape,
            source,
            target,
            self.lower.map(|height| (height, true)),
        )
    }
}

pub(crate) fn low_matrix(
    model: &mut Model,
    shape: Shape,
    source: &[Monomial],
    target: &[Monomial],
    cutoff: usize,
) -> Result<Matrix> {
    restricted_matrix(model, shape, source, target, Some((cutoff, false)))
}

fn restricted_matrix(
    model: &mut Model,
    shape: Shape,
    source: &[Monomial],
    target: &[Monomial],
    filter: Option<(usize, bool)>,
) -> Result<Matrix> {
    let mut matrix = Matrix::new(target.len(), source.len())?;
    if target.is_empty() {
        return Ok(matrix);
    }
    model.spend(target.len() * source.len())?;
    for (column, monomial) in source.iter().enumerate() {
        let cochain = Polynomial::from_terms(model, [(monomial, 1)])?;
        let image = bp::differential(model, shape, &cochain)?;
        if filter.is_some() {
            model.spend(image.cells())?;
        }
        for (monomial, coefficient) in image.terms() {
            if filter.is_some_and(|(height, keep_high)| {
                bp_rel::retains(monomial, height, shape.height) != keep_high
            }) {
                continue;
            }
            let row = index_of(model, target, monomial)?;
            matrix.set(row, column, coefficient);
        }
    }
    Ok(matrix)
}

pub(crate) fn coordinates(
    model: &mut Model,
    cochain: &Polynomial,
    basis: &[Monomial],
) -> Result<Vec<i128>> {
    model.spend(basis.len())?;
    let mut result = vec![0; basis.len()];
    for (monomial, coefficient) in cochain.terms() {
        result[index_of(model, basis, monomial)?] = coefficient;
    }
    Ok(result)
}

pub(crate) fn index_of(
    model: &mut Model,
    basis: &[Monomial],
    monomial: &Monomial,
) -> Result<usize> {
    model.spend((1 + 2 * monomial.len()) * (1 + basis.len().max(1).ilog2() as usize))?;
    basis
        .binary_search(monomial)
        .map_err(|_| Error::plain("BP term is outside its normalized cobar basis"))
}
