use std::rc::Rc;

use crate::error::{Error, Result};
use crate::legendre_closed;
use crate::model::{Model, Value};
use crate::rv_qp_quad::{self, Element as QuadElement, Field, Wide};
use crate::rv_unram::{self, Ring};

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) family: Rc<str>,
    pub(crate) exponent: usize,
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Element {
    pub(crate) constant: i128,
    pub(crate) generator: i128,
}

impl Element {
    pub(crate) const ZERO: Self = Self::new(0, 0);
    const ONE: Self = Self::new(1, 0);
    const GENERATOR: Self = Self::new(0, 1);

    pub(crate) const fn new(constant: i128, generator: i128) -> Self {
        Self {
            constant,
            generator,
        }
    }

    pub(crate) fn add(self, model: &mut Model, other: Self) -> Result<Self> {
        model.spend(1)?;
        Ok(Self::new(
            add(self.constant, other.constant)?,
            add(self.generator, other.generator)?,
        ))
    }

    pub(crate) fn negate(self, model: &mut Model) -> Result<Self> {
        self.scale(model, -1)
    }

    fn scale(self, model: &mut Model, factor: i128) -> Result<Self> {
        model.spend(1)?;
        Ok(Self::new(
            mul(self.constant, factor)?,
            mul(self.generator, factor)?,
        ))
    }

    pub(crate) fn product(
        self,
        model: &mut Model,
        other: Self,
        trace: i128,
        norm: i128,
    ) -> Result<Element> {
        model.spend(80)?;
        type Integer = Wide<13>;
        let product = Integer::product(self.generator, other.generator);
        let constant = Integer::product(self.constant, other.constant)
            .add(product.multiply(norm)?.negate())?;
        let generator = Integer::product(self.constant, other.generator)
            .add(Integer::product(self.generator, other.constant))?
            .add(product.multiply(trace)?)?;
        Ok(Self::new(
            constant.signed().ok_or_else(arithmetic_limit)?,
            generator.signed().ok_or_else(arithmetic_limit)?,
        ))
    }
}

#[derive(Clone, Copy, Eq, PartialEq)]
pub(crate) struct Family {
    pub(crate) prime: usize,
    pub(crate) trace: i128,
    norm: i128,
    pub(crate) discriminant: i128,
}

impl Family {
    pub(crate) fn new(model: &mut Model, ring: &Ring, lambda: &[u64]) -> Result<Self> {
        legendre_closed::validate(ring, lambda)?;
        let prime = ring.prime() as usize;
        let norm = i128::from(ring.size_of_residue_field()?);
        let trace = i128::from(ring.legendre_trace(model, lambda)?);
        Self::from_frobenius(prime, trace, norm)
    }

    pub(crate) fn from_frobenius(prime: usize, trace: i128, norm: i128) -> Result<Self> {
        if trace % prime as i128 == 0 {
            return Err(Error::plain("Legendre CM order requires an ordinary curve"));
        }
        let discriminant = sub(mul(trace, trace)?, mul(4, norm)?)?;
        if discriminant >= 0 {
            return Err(Error::plain(
                "ordinary Legendre Frobenius discriminant must be negative",
            ));
        }
        Ok(Self {
            prime,
            trace,
            norm,
            discriminant,
        })
    }

    // Work in the fixed basis (1, pi) with relation pi^2 = trace*pi - norm.
    // This computes ( -q*u_(m-1), u_m ) without ever computing q^m.
    fn power(self, model: &mut Model, mut exponent: usize) -> Result<Element> {
        let mut output = Element::ONE;
        let mut factor = Element::GENERATOR;
        while exponent != 0 {
            if exponent & 1 != 0 {
                output = output.product(model, factor, self.trace, self.norm)?;
            }
            exponent >>= 1;
            if exponent != 0 {
                factor = factor.product(model, factor, self.trace, self.norm)?;
            }
        }
        Ok(output)
    }

    fn power_trace(self, power: Element) -> Result<i128> {
        add(mul(2, power.constant)?, mul(self.trace, power.generator)?)
    }

    fn power_norm(self, model: &mut Model, mut exponent: usize) -> Result<i128> {
        let mut output = 1;
        let mut factor = self.norm;
        while exponent != 0 {
            model.spend(1)?;
            if exponent & 1 != 0 {
                output = mul(output, factor)?;
            }
            exponent >>= 1;
            if exponent != 0 {
                factor = mul(factor, factor)?;
            }
        }
        Ok(output)
    }

    fn field_element(self, value: Element) -> Result<QuadElement> {
        let Element {
            constant,
            generator,
        } = value;
        if generator & 1 == 0 {
            let half = generator / 2;
            QuadElement::new(add(constant, mul(self.trace, half)?)?, half, 1)
        } else {
            let a = add(mul(2, constant)?, mul(self.trace, generator)?)?;
            QuadElement::new(a, generator, 2)
        }
    }

    pub(crate) fn pi(self, model: &mut Model, exponent: usize) -> Result<Value> {
        let power = self.power(model, exponent)?;
        rv_qp_quad::persist(model, self.field_element(power)?)
    }

    pub(crate) fn lift(
        self,
        model: &mut Model,
        exponent: usize,
        coefficient: Element,
    ) -> Result<Value> {
        let value = if coefficient.generator == 0 {
            coefficient
        } else {
            self.power(model, exponent)?
                .scale(model, coefficient.generator)?
                .add(model, Element::new(coefficient.constant, 0))?
        };
        rv_qp_quad::persist(model, self.field_element(value)?)
    }

    pub(crate) fn multiply(
        self,
        model: &mut Model,
        exponent: usize,
        left: Element,
        right: Element,
    ) -> Result<Element> {
        if left.generator == 0 {
            return right.scale(model, left.constant);
        }
        if right.generator == 0 {
            return left.scale(model, right.constant);
        }
        let trace = self.power_trace(self.power(model, exponent)?)?;
        let norm = self.power_norm(model, exponent)?;
        left.product(model, right, trace, norm)
    }

    fn relative_generator(
        self,
        model: &mut Model,
        target: usize,
        source: usize,
    ) -> Result<(Element, Element, i128)> {
        let target = self.power(model, target)?;
        let source = self.power(model, source)?;
        let generator = source
            .generator
            .checked_div(target.generator)
            .ok_or_else(arithmetic_limit)?;
        if generator == 0 || mul(generator, target.generator)? != source.generator {
            return Err(Error::plain(
                "CM order index is not a nonzero integral quotient",
            ));
        }
        Ok((target, source, generator))
    }

    pub(crate) fn include(
        self,
        model: &mut Model,
        target: usize,
        source: usize,
        value: Element,
    ) -> Result<Element> {
        require_inclusion(target, source)?;
        if target == source || value.generator == 0 {
            return Ok(value);
        }
        let (target_power, source_power, coefficient) =
            self.relative_generator(model, target, source)?;
        let constant = sub(
            source_power.constant,
            mul(coefficient, target_power.constant)?,
        )?;
        Element::new(constant, coefficient)
            .scale(model, value.generator)?
            .add(model, Element::new(value.constant, 0))
    }

    pub(crate) fn report(
        self,
        model: &mut Model,
        family: &str,
        exponent: usize,
        query: Query,
    ) -> Result<String> {
        require_positive(exponent)?;
        match query {
            Query::Order => {
                let power = self.power(model, exponent)?;
                let trace = self.power_trace(power)?;
                let index = power.generator.unsigned_abs();
                let discriminant = mul(mul(self.discriminant, power.generator)?, power.generator)?;
                Ok(format!("rvCMOrder[{family}, {exponent}] = trace={trace}, index={index}, discriminant={discriminant}"))
            }
            Query::Trace => {
                let trace = self.power_trace(self.power(model, exponent)?)?;
                Ok(format!("rvCMTrace[{family}, {exponent}] = {trace}"))
            }
            Query::Index(source) => {
                require_inclusion(exponent, source)?;
                let index = if exponent == source {
                    1
                } else {
                    self.relative_generator(model, exponent, source)?
                        .2
                        .unsigned_abs()
                };
                Ok(format!(
                    "rvCMIndex[{family}, {exponent}, {source}] = {index}"
                ))
            }
            Query::Norm => {
                let norm = self.power_norm(model, exponent)?;
                Ok(format!("rvCMNorm[{family}, {exponent}] = {norm}"))
            }
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Order,
    Trace,
    Index(usize),
    Norm,
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Lift,
    Include(usize),
    Negate,
}

impl Unary {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Lift => "rvCMLift",
            Self::Include(_) => "rvCMInclude",
            Self::Negate => "neg",
        }
    }
}

pub(crate) fn require_positive(exponent: usize) -> Result<()> {
    if exponent == 0 {
        return Err(Error::plain("CM order exponent must be positive"));
    }
    Ok(())
}

pub(crate) fn require_inclusion(target: usize, source: usize) -> Result<()> {
    require_positive(target)?;
    require_positive(source)?;
    if !source.is_multiple_of(target) {
        return Err(Error::plain(
            "CM order inclusion requires the target exponent to divide the source exponent",
        ));
    }
    Ok(())
}

pub(crate) fn require_same_prime(field: Field, ring: &Ring) -> Result<()> {
    if u64::from(field.prime().get()) != ring.prime() {
        return Err(Error::plain(
            "quadratic specialization requires the same residue characteristic",
        ));
    }
    Ok(())
}

pub(crate) fn specialize(
    model: &mut Model,
    field: Field,
    ring: &Ring,
    value: QuadElement,
) -> Result<Value> {
    let residue = field.residue(model, value, ring.precision())?;
    let coefficient = u64::try_from(residue)
        .map_err(|_| Error::plain("quadratic residue exceeds unramified coefficient range"))?;
    rv_unram::persist(model, ring.scalar(coefficient))
}

pub(crate) fn persist(model: &mut Model, value: Element) -> Result<Value> {
    model.allocate(4)?;
    Ok(Value::CMOrder(value))
}

fn arithmetic_limit() -> Error {
    Error::plain("CM order arithmetic exceeds signed i128 range")
}

fn mul(left: i128, right: i128) -> Result<i128> {
    left.checked_mul(right).ok_or_else(arithmetic_limit)
}

fn add(left: i128, right: i128) -> Result<i128> {
    left.checked_add(right).ok_or_else(arithmetic_limit)
}

fn sub(left: i128, right: i128) -> Result<i128> {
    left.checked_sub(right).ok_or_else(arithmetic_limit)
}
