use super::{accumulate, checked, monomial_index, nonzero, ratio, Class, Form, MIDDLE, RANK};
use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::Model;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Lower,
    Lambda,
    Primitive,
}

impl Operation {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Lower => "rvHilbSquareLower",
            Self::Lambda => "rvHilbSquareLambda",
            Self::Primitive => "rvHilbSquarePrimitive",
        }
    }

    pub(crate) fn degree(self, degree: i128) -> Result<i128> {
        match self {
            Self::Lower | Self::Lambda => degree.checked_sub(2).ok_or_else(|| {
                Error::plain("Hilbert-square lowering degree exceeds the signed 128-bit range")
            }),
            Self::Primitive if matches!(degree, 0 | 2 | 4) => Ok(degree),
            Self::Primitive => Err(Error::plain(
                "rvHilbSquarePrimitive requires cohomological degree zero, two, or four",
            )),
        }
    }

    pub(crate) fn apply(
        self,
        model: &mut Model,
        form: &Form,
        generator: &Class,
        class: &Class,
    ) -> Result<Class> {
        if matches!(self, Self::Lower) {
            return lower(model, form, generator, class);
        }
        let contractions = form.contract(model, &generator.coordinates)?;
        let square = dot(model, &contractions, &generator.coordinates)?;
        if square == Rational::ZERO {
            return Err(Error::plain(format!(
                "{} requires a nonisotropic generator",
                self.name()
            )));
        }
        let inverse = ratio(square.denominator, square.numerator)?;
        match self {
            Self::Lambda => lower(model, form, generator, class)?.scale(model, inverse),
            Self::Primitive => primitive(model, form, generator, class, &contractions, inverse),
            Self::Lower => unreachable!(),
        }
    }
}

fn lower(model: &mut Model, form: &Form, generator: &Class, class: &Class) -> Result<Class> {
    let degree = Operation::Lower.degree(class.degree)?;
    let two = Rational::integer(2);
    let four = Rational::integer(4);
    match class.degree {
        2 => {
            let contractions = form.contract(model, &generator.coordinates)?;
            let pairing = dot(model, &contractions, &class.coordinates)?;
            Ok(Class {
                degree,
                coordinates: vec![checked(pairing.multiply(four))?],
            })
        }
        4 => {
            let mut image = generator.cup(model, form, class)?;
            image.degree = degree;
            image.scale(model, two)
        }
        6 => {
            // H⁶ and H² use the same Beauville coordinates.
            let vector = class.with_vector_degree(model, 2)?;
            generator.cup(model, form, &vector)?.scale(model, two)
        }
        8 => {
            let scalar = checked(four.multiply(class.coordinates[0]))?;
            let mut image = generator.scale(model, scalar)?;
            image.degree = degree;
            Ok(image)
        }
        _ => Ok(Class::zero(degree)),
    }
}

fn primitive(
    model: &mut Model,
    form: &Form,
    generator: &Class,
    class: &Class,
    contractions: &[Rational; RANK],
    inverse: Rational,
) -> Result<Class> {
    let negative_inverse = checked(inverse.multiply(Rational::integer(-1)))?;
    match class.degree {
        0 => class.scale(model, Rational::integer(1)),
        2 => {
            let pairing = dot(model, contractions, &class.coordinates)?;
            let coefficient = checked(pairing.multiply(negative_inverse))?;
            let correction = generator.scale(model, coefficient)?;
            class.add(model, &correction)
        }
        4 => {
            let mut w = generator.cup(model, form, class)?;
            w.degree = 2;
            let pairing = dot(model, contractions, &w.coordinates)?;
            let aw = generator
                .cup(model, form, &w)?
                .scale(model, negative_inverse)?;
            let coefficient = checked(pairing.multiply(inverse))?;
            let coefficient = checked(coefficient.multiply(inverse))?;
            let coefficient = checked(coefficient.multiply(ratio(2, 3)?))?;
            let square = generator
                .cup(model, form, generator)?
                .scale(model, coefficient)?;
            class.add(model, &aw)?.add(model, &square)
        }
        _ => unreachable!("primitive degree checked by the interpreter"),
    }
}

pub(crate) fn rotate(
    model: &mut Model,
    form: &Form,
    a: &Class,
    b: &Class,
    class: &Class,
) -> Result<Class> {
    let mut output = Class::zero(class.degree);
    if !matches!(class.degree, 2 | 4 | 6) {
        return Ok(output);
    }
    let covector_a = form.contract(model, &a.coordinates)?;
    let covector_b = form.contract(model, &b.coordinates)?;
    if class.degree != 4 {
        let q_b = dot(model, &covector_b, &class.coordinates)?;
        let q_a = dot(model, &covector_a, &class.coordinates)?;
        let q_a = checked(q_a.multiply(Rational::integer(-1)))?;
        for (coordinate, (&ai, &bi)) in output
            .coordinates
            .iter_mut()
            .zip(a.coordinates.iter().zip(&b.coordinates))
        {
            accumulate(coordinate, q_b, ai)?;
            accumulate(coordinate, q_a, bi)?;
        }
        return Ok(output);
    }

    model.spend(RANK * RANK + MIDDLE * (2 * RANK + 1))?;
    // Each column is R_ab(e_j) = q(b,e_j)a − q(a,e_j)b. Apply the
    // same matrix to both factors of every symmetric monomial.
    let mut columns = [[Rational::ZERO; RANK]; RANK];
    for (j, column) in columns.iter_mut().enumerate() {
        let q_a = checked(covector_a[j].multiply(Rational::integer(-1)))?;
        for (coordinate, (&ai, &bi)) in column
            .iter_mut()
            .zip(a.coordinates.iter().zip(&b.coordinates))
        {
            accumulate(coordinate, covector_b[j], ai)?;
            accumulate(coordinate, q_a, bi)?;
        }
    }
    for (index, coefficient) in nonzero(&class.coordinates) {
        let (i, j) = form.monomials[index];
        for (axis, (&left, &right)) in columns[i].iter().zip(&columns[j]).enumerate() {
            let first = monomial_index(axis.min(j), axis.max(j));
            let second = monomial_index(i.min(axis), i.max(axis));
            accumulate(&mut output.coordinates[first], coefficient, left)?;
            accumulate(&mut output.coordinates[second], coefficient, right)?;
        }
    }
    Ok(output)
}

fn dot(model: &mut Model, left: &[Rational], right: &[Rational]) -> Result<Rational> {
    model.spend(RANK)?;
    let mut result = Rational::ZERO;
    for (&a, &b) in left.iter().zip(right) {
        accumulate(&mut result, a, b)?;
    }
    Ok(result)
}

impl Class {
    fn with_vector_degree(&self, model: &mut Model, degree: i128) -> Result<Self> {
        debug_assert!(matches!(self.degree, 2 | 6) && matches!(degree, 2 | 6));
        model.spend(self.coordinates.len())?;
        Ok(Self {
            degree,
            coordinates: self.coordinates.clone(),
        })
    }
}
