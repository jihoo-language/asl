use std::fmt;
use std::rc::Rc;

use crate::bp_polynomial::Polynomial;
use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;

const MAX_HEIGHT: usize = 8;
const MAX_ARROWS: usize = 16;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Shape {
    pub(crate) prime: u32,
    pub(crate) height: usize,
    pub(crate) arrows: usize,
    pub(crate) depth: usize,
}

impl fmt::Display for Shape {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self {
            prime,
            height,
            arrows,
            depth,
        } = *self;
        if depth == 0 {
            write!(out, "BP({prime}, {height}, {arrows})")
        } else {
            write!(out, "BPI({prime}, {height}, {arrows}, {depth})")
        }
    }
}

impl Shape {
    pub(crate) fn new(model: &mut Model, p: usize, height: usize, arrows: usize) -> Result<Self> {
        let prime = match p {
            2 => 2,
            n if n < 2 || n % 2 == 0 => return Err(Error::plain("BP requires a prime")),
            n => Prime::new(model, n)?.get(),
        };
        Self {
            prime,
            height: 0,
            arrows: 0,
            depth: 0,
        }
        .raise(height)?
        .with_arrows(arrows)
    }

    pub(crate) fn raise(self, height: usize) -> Result<Self> {
        if height < self.height {
            return Err(Error::plain("bpRaise cannot decrease the BP height"));
        }
        if height > MAX_HEIGHT {
            return Err(Error::plain(format!(
                "BP height exceeds the limit of {MAX_HEIGHT}"
            )));
        }
        Ok(Self { height, ..self })
    }

    pub(crate) fn quotient(self, depth: usize) -> Result<Self> {
        if depth == 0 || depth > self.height + 1 {
            return Err(Error::plain(format!(
                "BPI ideal depth must be in 1..={}",
                self.height + 1
            )));
        }
        if depth < self.depth {
            return Err(Error::plain("bpReduce cannot decrease the ideal depth"));
        }
        Ok(Self { depth, ..self })
    }

    fn coefficient_prime(self) -> Option<u32> {
        (self.depth != 0).then_some(self.prime)
    }

    pub(crate) fn normalize(self, model: &mut Model, polynomial: &mut Polynomial) -> Result<()> {
        if let Some(prime) = self.coefficient_prime() {
            polynomial.reduce(model, prime, self.depth - 1)?;
        }
        Ok(())
    }

    pub(crate) fn require_integral(self, operation: &str) -> Result<()> {
        if self.depth != 0 {
            return Err(Error::plain(format!(
                "{operation} requires an integral BP cochain"
            )));
        }
        Ok(())
    }

    fn require_quotient(self, operation: &str) -> Result<()> {
        if self.depth == 0 {
            return Err(Error::plain(format!("{operation} requires a BPI cochain")));
        }
        Ok(())
    }

    pub(crate) fn next(self) -> Result<Self> {
        self.with_arrows(self.arrows + 1)
    }

    fn with_arrows(self, arrows: usize) -> Result<Self> {
        if arrows > MAX_ARROWS {
            return Err(Error::plain(format!(
                "BP arrow arity exceeds the limit of {MAX_ARROWS}"
            )));
        }
        Ok(Self { arrows, ..self })
    }

    fn same_ring(self, other: Self) -> bool {
        self.prime == other.prime && self.height == other.height && self.depth == other.depth
    }

    pub(crate) fn cup(self, other: Self) -> Result<Self> {
        if !self.same_ring(other) {
            return Err(Error::plain(
                "bpCup requires the same prime, height, and ideal depth",
            ));
        }
        self.with_arrows(self.arrows + other.arrows)
    }

    pub(crate) fn variable(self, slot: usize, index: usize) -> Result<Polynomial> {
        if slot > self.arrows {
            return Err(Error::plain(format!(
                "BP arrow slot must be in 1..={}",
                self.arrows
            )));
        }
        if index == 0 || index > self.height {
            return Err(Error::plain(format!(
                "BP generator index must be in 1..={}",
                self.height
            )));
        }
        Ok(Polynomial::variable(slot * self.height + index - 1))
    }

    fn variables(self, slot: usize, count: usize) -> Vec<Polynomial> {
        (0..count)
            .map(|i| Polynomial::variable(slot * self.height + i))
            .collect()
    }
}

pub(crate) fn persist(model: &mut Model, polynomial: Polynomial) -> Result<Value> {
    model.allocate(polynomial.cells())?;
    Ok(Value::BP(Rc::new(polynomial)))
}

pub(crate) fn zero() -> Value {
    Value::BP(Rc::new(Polynomial::default()))
}

// Widen or narrow the offsets of the source and arrow slots. At lower
// width this retains the component involving only lower-height variables.
pub(crate) fn reindex_height(
    model: &mut Model,
    source_width: usize,
    target_width: usize,
    input: &Polynomial,
) -> Result<Polynomial> {
    if source_width == 0 || source_width == target_width {
        model.spend(input.cells())?;
        return Ok(input.clone());
    }
    input.reindex(model, |index| {
        let level = index % source_width;
        (level < target_width).then(|| (index / source_width) * target_width + level)
    })
}

// A Hazewinkel log polynomial has an integral numerator and a denominator
// that is a power of p. Logs are temporary. Every public coordinate formula
// is converted back by exact division before it is used in a coface.
#[derive(Clone)]
struct Log {
    numerator: Polynomial,
    denominator: u32,
}

impl Log {
    fn one() -> Self {
        Self {
            numerator: Polynomial::integer(1),
            denominator: 0,
        }
    }

    fn reduce(&mut self, model: &mut Model, prime: u32) -> Result<()> {
        while self.denominator != 0 && self.numerator.divide_if_exact(model, u128::from(prime))? {
            self.denominator -= 1;
        }
        Ok(())
    }

    fn add(
        &mut self,
        model: &mut Model,
        prime: u32,
        right: &Log,
        factor: &Polynomial,
        subtract: bool,
    ) -> Result<()> {
        if right.numerator.is_zero() || factor.is_zero() {
            return Ok(());
        }
        let mut term = right.numerator.multiply(model, factor)?;
        if self.numerator.is_zero() {
            self.denominator = right.denominator;
        }
        if self.denominator < right.denominator {
            let scale = denominator(prime, right.denominator - self.denominator)?;
            self.numerator.scale(model, scale)?;
            self.denominator = right.denominator;
        }
        let scale = denominator(prime, self.denominator - right.denominator)?;
        term.scale(model, scale)?;
        self.numerator.add_assign(model, &term, subtract)?;
        self.reduce(model, prime)
    }

    fn times_prime(&mut self, model: &mut Model, prime: u32) -> Result<()> {
        if self.denominator == 0 {
            self.numerator.scale(model, i128::from(prime))
        } else {
            self.denominator -= 1;
            Ok(())
        }
    }

    fn integral(mut self, model: &mut Model, prime: u32) -> Result<Polynomial> {
        if self.denominator != 0 {
            self.numerator
                .divide_exact(model, denominator(prime, self.denominator)?)?;
        }
        Ok(self.numerator)
    }
}

fn denominator(prime: u32, exponent: u32) -> Result<i128> {
    i128::from(prime)
        .checked_pow(exponent)
        .ok_or_else(|| Error::plain("BP log denominator exceeds the signed i128 range"))
}

fn prime_power(prime: u32, index: usize) -> Result<usize> {
    (prime as usize)
        .checked_pow(index as u32)
        .ok_or_else(|| Error::plain("BP prime-power exponent is out of range"))
}

fn logs(model: &mut Model, prime: u32, generators: &[Polynomial]) -> Result<Vec<Log>> {
    let mut result = vec![Log::one()];
    for n in 1..=generators.len() {
        model.spend(generators[n - 1].cells())?;
        let mut sum = Log {
            numerator: generators[n - 1].clone(),
            denominator: 0,
        };
        for (i, previous) in result.iter().enumerate().skip(1) {
            let factor = generators[n - i - 1].power(model, prime_power(prime, i)?)?;
            sum.add(model, prime, previous, &factor, false)?;
        }
        sum.denominator += 1;
        sum.reduce(model, prime)?;
        result.push(sum);
    }
    Ok(result)
}

fn right_logs(
    model: &mut Model,
    prime: u32,
    logs: &[Log],
    arrow: &[Polynomial],
) -> Result<Vec<Log>> {
    let mut result = vec![Log::one()];
    for (n, log) in logs.iter().enumerate().skip(1) {
        model.spend(log.numerator.cells())?;
        let mut sum = log.clone();
        for (i, previous) in logs[..n].iter().enumerate() {
            let factor = arrow[n - i - 1].power(model, prime_power(prime, i)?)?;
            sum.add(model, prime, previous, &factor, false)?;
        }
        result.push(sum);
    }
    Ok(result)
}

fn right_unit(
    model: &mut Model,
    prime: u32,
    base: &[Polynomial],
    arrow: &[Polynomial],
) -> Result<Vec<Polynomial>> {
    let logs = logs(model, prime, base)?;
    let right = right_logs(model, prime, &logs, arrow)?;
    let mut result: Vec<Polynomial> = Vec::new();
    // The right-unit images satisfy the same recursion as the original v's.
    for (n, log) in right.iter().enumerate().skip(1) {
        model.spend(log.numerator.cells())?;
        let mut sum = log.clone();
        sum.times_prime(model, prime)?;
        for (i, previous) in right[..n].iter().enumerate().skip(1) {
            let factor = result[n - i - 1].power(model, prime_power(prime, i)?)?;
            sum.add(model, prime, previous, &factor, true)?;
        }
        result.push(sum.integral(model, prime)?);
    }
    Ok(result)
}

fn composition(
    model: &mut Model,
    prime: u32,
    base: &[Polynomial],
    first: &[Polynomial],
    second: &[Polynomial],
) -> Result<Vec<Polynomial>> {
    let logs = logs(model, prime, base)?;
    let first_right = right_logs(model, prime, &logs, first)?;
    let second_right = right_logs(model, prime, &first_right, second)?;
    let mut result: Vec<Polynomial> = Vec::new();
    // Solve for the composite coordinates in increasing index order.
    for (n, log) in second_right.into_iter().enumerate().skip(1) {
        let mut sum = log;
        for (i, previous) in logs[..=n].iter().enumerate().skip(1) {
            let factor = if i == n {
                Polynomial::integer(1)
            } else {
                result[n - i - 1].power(model, prime_power(prime, i)?)?
            };
            sum.add(model, prime, previous, &factor, true)?;
        }
        result.push(sum.integral(model, prime)?);
    }
    Ok(result)
}

pub(crate) fn contracting_numerator(
    model: &mut Model,
    shape: Shape,
    input: &Polynomial,
) -> Result<(Polynomial, u32)> {
    // Insert the origin before the first vertex of the simplex. The arrow
    // from the origin to v has coordinates t_i = L_i(v). Pullback along
    // this insertion contracts the normalized pair-groupoid cochains.
    let height = input.max_used_index(model, shape.height, 1)?;
    let logarithms = logs(model, shape.prime, &shape.variables(0, height))?;
    let mut output = Log {
        numerator: Polynomial::default(),
        denominator: 0,
    };
    for (monomial, coefficient) in input.terms() {
        model.spend(1 + monomial.len())?;
        if monomial
            .first()
            .is_some_and(|&(variable, _)| variable < shape.height)
        {
            continue;
        }
        let mut term = Log {
            numerator: Polynomial::integer(coefficient),
            denominator: 0,
        };
        for &(variable, exponent) in monomial {
            let shifted = variable - shape.height;
            if shifted < shape.height {
                let logarithm = &logarithms[shifted + 1];
                let power = logarithm.numerator.power(model, exponent)?;
                term.numerator = term.numerator.multiply(model, &power)?;
                if logarithm.denominator != 0 {
                    let exponent = u32::try_from(exponent).map_err(|_| contraction_overflow())?;
                    let denominator = logarithm
                        .denominator
                        .checked_mul(exponent)
                        .ok_or_else(contraction_overflow)?;
                    term.denominator = term
                        .denominator
                        .checked_add(denominator)
                        .ok_or_else(contraction_overflow)?;
                }
            } else {
                let power = Polynomial::variable(shifted).power(model, exponent)?;
                term.numerator = term.numerator.multiply(model, &power)?;
            }
        }
        output.add(model, shape.prime, &term, &Polynomial::integer(1), false)?;
    }
    Ok((output.numerator, output.denominator))
}

fn contraction_overflow() -> Error {
    Error::plain("BP rational contraction denominator exponent is out of range")
}

pub(crate) fn check_face(shape: Shape, face: usize) -> Result<()> {
    if face > shape.arrows + 1 {
        return Err(Error::plain(format!(
            "BP coface index must be in 0..={}",
            shape.arrows + 1
        )));
    }
    Ok(())
}

pub(crate) fn coface(
    model: &mut Model,
    shape: Shape,
    face: usize,
    input: &Polynomial,
) -> Result<Polynomial> {
    check_face(shape, face)?;
    let target = shape.next()?;
    if face == shape.arrows + 1 {
        model.spend(input.cells())?;
        return Ok(input.clone());
    }
    let used_height = input.max_used_index(model, shape.height, face)?;
    let mut images = Vec::with_capacity((shape.arrows + 1) * shape.height);
    for slot in 0..=shape.arrows {
        let mut variables = target.variables(slot + usize::from(slot > face), shape.height);
        if slot == face && used_height != 0 && (face != 0 || used_height > target.depth) {
            let arrow = target.variables(slot + 1, used_height);
            let transported = transport_base(model, target, face.saturating_sub(1), used_height)?;
            let coordinates = if face == 0 {
                right_unit(model, shape.prime, &transported, &arrow)?
            } else {
                let first = target.variables(slot, used_height);
                composition(model, shape.prime, &transported, &first, &arrow)?
            };
            variables.splice(..used_height, coordinates);
        }
        images.extend(variables);
    }
    for image in &mut images {
        target.normalize(model, image)?;
    }
    let mut output = input.substitute(model, &images, target.coefficient_prime())?;
    target.normalize(model, &mut output)?;
    Ok(output)
}

fn transport_base(
    model: &mut Model,
    target: Shape,
    arrows: usize,
    height: usize,
) -> Result<Vec<Polynomial>> {
    let mut base = target.variables(0, height);
    for variable in base.iter_mut().take(target.depth.saturating_sub(1)) {
        *variable = Polynomial::default();
    }
    // The first depth generators are invariant in this quotient.
    if height > target.depth {
        for slot in 1..=arrows {
            let arrow = target.variables(slot, height);
            base = right_unit(model, target.prime, &base, &arrow)?;
            for coordinate in &mut base {
                target.normalize(model, coordinate)?;
            }
        }
    }
    Ok(base)
}

pub(crate) fn differential(
    model: &mut Model,
    shape: Shape,
    input: &Polynomial,
) -> Result<Polynomial> {
    shape.next()?;
    let mut result = Polynomial::default();
    for face in 0..=shape.arrows + 1 {
        let image = coface(model, shape, face, input)?;
        result.add_assign(model, &image, face % 2 != 0)?;
    }
    shape.normalize(model, &mut result)?;
    Ok(result)
}

pub(crate) fn cup(
    model: &mut Model,
    left: Shape,
    right: Shape,
    f: &Polynomial,
    g: &Polynomial,
) -> Result<Polynomial> {
    let target = left.cup(right)?;
    if f.is_zero() || g.is_zero() {
        return Ok(Polynomial::default());
    }
    let used_height = g.max_used_index(model, target.height, 0)?;
    let base = transport_base(model, target, left.arrows, used_height)?;
    let mut images = target.variables(0, target.height);
    images.splice(..used_height, base);
    for slot in 1..=right.arrows {
        images.extend(target.variables(left.arrows + slot, target.height));
    }
    let mut g = g.substitute(model, &images, target.coefficient_prime())?;
    target.normalize(model, &mut g)?;
    let mut output = f.multiply(model, &g)?;
    target.normalize(model, &mut output)?;
    Ok(output)
}

fn degeneracy_target(shape: Shape, index: usize) -> Result<Shape> {
    if shape.arrows == 0 {
        return Err(Error::plain("bpDegeneracy requires at least one arrow"));
    }
    if index >= shape.arrows {
        return Err(Error::plain(format!(
            "BP degeneracy index must be in 0..={}",
            shape.arrows - 1
        )));
    }
    Ok(Shape {
        arrows: shape.arrows - 1,
        ..shape
    })
}

fn degeneracy(
    model: &mut Model,
    shape: Shape,
    index: usize,
    input: &Polynomial,
) -> Result<Polynomial> {
    let start = (index + 1) * shape.height;
    let end = start + shape.height;
    input.reindex(model, |variable| {
        if variable >= end {
            Some(variable - shape.height)
        } else if variable >= start {
            None
        } else {
            Some(variable)
        }
    })
}

fn normalized(model: &mut Model, shape: Shape, input: &Polynomial) -> Result<Polynomial> {
    model.spend(input.cells())?;
    let mut output = input.clone();
    for index in 0..shape.arrows {
        let lower = degeneracy(model, shape, index, &output)?;
        if lower.is_zero() {
            continue;
        }
        // This coface is a section of the degeneracy and preserves the earlier kernels.
        let correction = coface(model, degeneracy_target(shape, index)?, index + 1, &lower)?;
        output.add_assign(model, &correction, true)?;
        shape.normalize(model, &mut output)?;
    }
    Ok(output)
}

pub(crate) struct Massey {
    shapes: [Shape; 5],
    pub(crate) target: Shape,
}

impl Massey {
    pub(crate) fn new(shapes: [Shape; 5]) -> Result<Self> {
        let [a, b, c, u, v] = shapes;
        if shapes[1..].iter().any(|&other| !a.same_ring(other)) {
            return Err(Error::plain(
                "bpMassey requires the same prime, height, and ideal depth",
            ));
        }
        if u.arrows + 1 != a.arrows + b.arrows {
            return Err(Error::plain("bpMassey witness U must have arity a+b-1"));
        }
        if v.arrows + 1 != b.arrows + c.arrows {
            return Err(Error::plain("bpMassey witness V must have arity b+c-1"));
        }
        for shape in shapes {
            shape.next()?;
        }
        let target = a.cup(v)?;
        Ok(Self { shapes, target })
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        [a, b, c, u, v]: [&Polynomial; 5],
    ) -> Result<Polynomial> {
        let [a_shape, b_shape, c_shape, u_shape, v_shape] = self.shapes;
        for (shape, value, name) in [(a_shape, a, "A"), (b_shape, b, "B"), (c_shape, c, "C")] {
            if !differential(model, shape, value)?.is_zero() {
                return Err(Error::plain(format!(
                    "bpMassey input {name} must be closed"
                )));
            }
        }
        for (shape, witness, left_shape, left, right_shape, right, equation) in [
            (u_shape, u, a_shape, a, b_shape, b, "d(U) = bpCup(A, B)"),
            (v_shape, v, b_shape, b, c_shape, c, "d(V) = bpCup(B, C)"),
        ] {
            let actual = differential(model, shape, witness)?;
            let expected = cup(model, left_shape, right_shape, left, right)?;
            if !actual.compare(model, &expected)?.is_eq() {
                return Err(Error::plain(format!("bpMassey requires {equation}")));
            }
        }
        let mut result = cup(model, a_shape, v_shape, a, v)?;
        let other = cup(model, u_shape, c_shape, u, c)?;
        result.add_assign(model, &other, a_shape.arrows % 2 == 0)?;
        self.target.normalize(model, &mut result)?;
        Ok(result)
    }
}

fn connecting(model: &mut Model, shape: Shape, input: &Polynomial) -> Result<Polynomial> {
    // The canonical residues can also be read as a lift into the preceding
    // quotient. Its differential lies in the kernel precisely for cycles.
    let previous = Shape {
        depth: shape.depth - 1,
        ..shape
    };
    let mut numerator = differential(model, previous, input)?;
    let noncycle = || Error::plain("bpDelta requires a closed BPI cochain");
    if previous.depth == 0 {
        if !numerator.divide_if_exact(model, u128::from(shape.prime))? {
            return Err(noncycle());
        }
        Ok(numerator)
    } else {
        numerator
            .divide_variable(model, previous.depth - 1)?
            .ok_or_else(noncycle)
    }
}

fn p_bockstein(
    model: &mut Model,
    shape: Shape,
    input: &Polynomial,
    exponent: usize,
) -> Result<Polynomial> {
    let mut numerator = differential(model, shape, input)?;
    if numerator.is_zero() {
        return Ok(numerator);
    }
    let divisor = u32::try_from(exponent)
        .ok()
        .and_then(|power| u128::from(shape.prime).checked_pow(power));
    if let Some(divisor) = divisor {
        if numerator.divide_if_exact(model, divisor)? {
            return Ok(numerator);
        }
    }
    Err(Error::plain(format!(
        "bpBockP requires the differential to be divisible by p^{exponent}"
    )))
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Negate,
    Power(usize),
    Coface(usize),
    Differential,
    Reduce(usize),
    Lift,
    Delta,
    BockP(usize),
    Degeneracy(usize),
    Normalize,
    Null(usize),
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Negate => "bpNeg",
            Self::Power(_) => "bpPow",
            Self::Coface(_) => "bpCoface",
            Self::Differential => "bpD",
            Self::Reduce(_) => "bpReduce",
            Self::Lift => "bpLift",
            Self::Delta => "bpDelta",
            Self::BockP(_) => "bpBockP",
            Self::Degeneracy(_) => "bpDegeneracy",
            Self::Normalize => "bpNormalize",
            Self::Null(_) => "bpNull",
        }
    }

    pub(crate) fn target(self, shape: Shape) -> Result<Shape> {
        match self {
            Self::Coface(face) => {
                check_face(shape, face)?;
                shape.next()
            }
            Self::Differential => shape.next(),
            Self::Null(_) => {
                shape.require_integral(self.name())?;
                crate::bp_ext::previous_shape(shape)
            }
            Self::Negate | Self::Power(_) | Self::Normalize => Ok(shape),
            Self::Degeneracy(index) => degeneracy_target(shape, index),
            Self::Reduce(depth) => shape.quotient(depth),
            Self::Lift => {
                shape.require_quotient(self.name())?;
                Ok(Shape { depth: 0, ..shape })
            }
            Self::Delta => {
                shape.require_quotient(self.name())?;
                Ok(Shape {
                    depth: shape.depth - 1,
                    ..shape.next()?
                })
            }
            Self::BockP(_) => {
                shape.require_integral(self.name())?;
                shape.next()
            }
        }
    }

    pub(crate) fn apply(
        self,
        model: &mut Model,
        shape: Shape,
        input: &Polynomial,
    ) -> Result<Polynomial> {
        let target = self.target(shape)?;
        match self {
            Self::Negate => {
                model.spend(input.cells())?;
                let mut output = input.clone();
                output.scale(model, -1)?;
                target.normalize(model, &mut output)?;
                Ok(output)
            }
            Self::Power(power) => input.power_mod(model, power, shape.coefficient_prime()),
            Self::Coface(face) => coface(model, shape, face, input),
            Self::Differential => differential(model, shape, input),
            Self::Reduce(_) | Self::Lift => {
                model.spend(input.cells())?;
                let mut output = input.clone();
                target.normalize(model, &mut output)?;
                Ok(output)
            }
            Self::Delta => connecting(model, shape, input),
            Self::BockP(exponent) => p_bockstein(model, shape, input, exponent),
            Self::Degeneracy(index) => degeneracy(model, shape, index, input),
            Self::Normalize => normalized(model, shape, input),
            Self::Null(t) => crate::bp_ext::null(model, shape, t, input),
        }
    }
}
