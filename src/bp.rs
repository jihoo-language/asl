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
        if height > MAX_HEIGHT {
            return Err(Error::plain(format!(
                "BP height exceeds the limit of {MAX_HEIGHT}"
            )));
        }
        if arrows > MAX_ARROWS {
            return Err(Error::plain(format!(
                "BP arrow arity exceeds the limit of {MAX_ARROWS}"
            )));
        }
        Ok(Self {
            prime,
            height,
            arrows,
            depth: 0,
        })
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

    fn require_quotient(self, operation: &str) -> Result<()> {
        if self.depth == 0 {
            return Err(Error::plain(format!("{operation} requires a BPI cochain")));
        }
        Ok(())
    }

    pub(crate) fn next(self) -> Result<Self> {
        if self.arrows == MAX_ARROWS {
            return Err(Error::plain(format!(
                "BP arrow arity exceeds the limit of {MAX_ARROWS}"
            )));
        }
        Ok(Self {
            arrows: self.arrows + 1,
            ..self
        })
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
    let mut transported = target.variables(0, used_height);
    for variable in transported.iter_mut().take(target.depth.saturating_sub(1)) {
        *variable = Polynomial::default();
    }
    // The first depth generators are invariant in this quotient.
    if used_height > target.depth {
        for slot in 1..face {
            let arrow = target.variables(slot, used_height);
            transported = right_unit(model, shape.prime, &transported, &arrow)?;
            for coordinate in &mut transported {
                target.normalize(model, coordinate)?;
            }
        }
    }
    let mut images = Vec::with_capacity((shape.arrows + 1) * shape.height);
    for slot in 0..=shape.arrows {
        let mut variables = target.variables(slot + usize::from(slot > face), shape.height);
        if slot == face && used_height != 0 && (face != 0 || used_height > target.depth) {
            let arrow = target.variables(slot + 1, used_height);
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
        }
    }

    pub(crate) fn target(self, shape: Shape) -> Result<Shape> {
        match self {
            Self::Coface(face) => {
                check_face(shape, face)?;
                shape.next()
            }
            Self::Differential => shape.next(),
            Self::Negate | Self::Power(_) => Ok(shape),
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
                if shape.depth != 0 {
                    return Err(Error::plain("bpBockP requires an integral BP cochain"));
                }
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
        }
    }
}
