use super::mixed::{connes_terms, differential_terms};
use super::{
    grading_error, linear, storage, Algebra, Chain, Complex, Error, Flavor, Form, Grade, Key,
    Model, Operation, Rational, Result,
};

impl Operation {
    pub(crate) fn transition(self, grade: &Grade, class: bool) -> Result<(Grade, bool)> {
        use Flavor::{Cyclic, Hochschild, Negative, Periodic};
        let result = match self {
            Self::Diff if !class => (grade.with_offset(-1)?, false),
            Self::Connes if grade.flavor == Hochschild => (grade.with_offset(1)?, class),
            Self::Class if !class => (grade.clone(), true),
            Self::Representative if class => (grade.clone(), false),
            Self::Primitive if !class => (grade.with_offset(1)?, false),
            Self::S if grade.flavor == Cyclic => (grade.with_offset(-2)?, class),
            Self::I if grade.flavor == Hochschild => (grade.with_flavor(Cyclic), class),
            Self::Boundary if grade.flavor == Cyclic && class => {
                (grade.with_flavor(Hochschild).with_offset(1)?, true)
            }
            Self::U if grade.flavor == Negative => (grade.with_offset(-2)?, class),
            Self::Eval if grade.flavor == Negative => (grade.with_flavor(Hochschild), class),
            Self::Periodic if grade.flavor == Negative => (grade.with_flavor(Periodic), class),
            Self::Lift if grade.flavor == Hochschild && class => {
                (grade.with_flavor(Negative), true)
            }
            _ => {
                return Err(Error::plain(format!(
                    "{} does not accept this {} {}",
                    self.name(),
                    grade.flavor.name(),
                    if class { "homology class" } else { "chain" }
                )));
            }
        };
        Ok(result)
    }
}

fn raw_map(model: &mut Model, operation: Operation, input: &Chain) -> Result<Chain> {
    let mut result = Chain::default();
    for (key, &coefficient) in &input.entries {
        match operation {
            Operation::Connes | Operation::Boundary => {
                if key.power == 0 {
                    connes_terms(model, key, 0, |model, key, sign| {
                        let value = coefficient.multiply(model, sign)?;
                        result.insert(model, key, value)
                    })?;
                }
            }
            Operation::S | Operation::U => {
                if matches!(operation, Operation::U) || key.power < 0 {
                    let power = key.power.checked_add(1).ok_or_else(grading_error)?;
                    result.insert(
                        model,
                        Key {
                            power,
                            word: key.word.clone(),
                        },
                        coefficient,
                    )?;
                }
            }
            Operation::Eval if key.power != 0 => (),
            Operation::I | Operation::Periodic | Operation::Eval => {
                result.insert(model, key.clone(), coefficient)?;
            }
            _ => unreachable!(),
        }
    }
    Ok(result)
}

fn boundary(model: &mut Model, algebra: &Algebra, grade: &Grade, input: &Chain) -> Result<Chain> {
    let mut result = Chain::default();
    for (key, &coefficient) in &input.entries {
        differential_terms(model, algebra, grade.flavor, key, |model, key, value| {
            let value = coefficient.multiply(model, value)?;
            result.insert(model, key, value)
        })?;
    }
    Ok(result)
}

impl Form {
    pub(crate) fn cyclic_operation(
        &self,
        model: &mut Model,
        grade: &Grade,
        operation: Operation,
        class: bool,
        value: &Chain,
    ) -> Result<Chain> {
        if value.empty() {
            return Ok(Chain::default());
        }
        let (target, target_class) = operation.transition(grade, class)?;
        let algebra = self.cyclic_algebra(model, grade)?;
        match operation {
            Operation::Diff => boundary(model, &algebra, grade, value),
            Operation::Representative => unreachable!(),
            Operation::Class => Complex::new(model, &algebra, grade)?.project(model, value),
            Operation::Primitive => {
                let primitive = Complex::new(model, &algebra, grade)?
                    .primitive(model, &algebra, grade, value)?;
                let verified = boundary(model, &algebra, &target, &primitive)?;
                if !verified.compare(model, value)?.is_eq() {
                    return Err(Error::plain(
                        "Verbitsky cyclic primitive failed its differential equation",
                    ));
                }
                Ok(primitive)
            }
            Operation::Lift => lift(model, &algebra, grade, &target, value),
            _ => {
                let result = raw_map(model, operation, value)?;
                if target_class && !result.empty() {
                    Complex::new(model, &algebra, &target)?.project(model, &result)
                } else {
                    Ok(result)
                }
            }
        }
    }

    pub(crate) fn cyclic_sbi(&self, model: &mut Model, grade: &Grade) -> Result<String> {
        use Flavor::{Cyclic, Hochschild};
        let grades = [
            grade.with_flavor(Cyclic).with_offset(-1)?,
            grade.with_flavor(Hochschild),
            grade.with_flavor(Cyclic),
            grade.with_flavor(Cyclic).with_offset(-2)?,
            grade.with_flavor(Hochschild).with_offset(-1)?,
            grade.with_flavor(Cyclic).with_offset(-1)?,
        ];
        let operators = [
            Operation::Boundary,
            Operation::I,
            Operation::S,
            Operation::Boundary,
            Operation::I,
        ];
        let algebra = self.cyclic_algebra(model, grade)?;
        let mut source = Complex::new(model, &algebra, &grades[0])?;
        let mut previous_image = Vec::new();
        let mut ranks = [0; 5];
        for (i, &operation) in operators.iter().enumerate() {
            let target = Complex::new(model, &algebra, &grades[i + 1])?;
            for chain in &previous_image {
                let composition = raw_map(model, operation, chain)?;
                if !target.project(model, &composition)?.empty() {
                    return Err(Error::plain(
                        "Verbitsky SBI composition is nonzero on homology",
                    ));
                }
            }
            let (image, rank) = image(model, &source, &target, operation)?;
            ranks[i] = rank;
            if i > 0 && ranks[i - 1] + ranks[i] != source.dimensions()[3] {
                return Err(Error::plain(
                    "Verbitsky SBI image and kernel dimensions disagree",
                ));
            }
            previous_image = image;
            source = target;
        }
        Ok(format!(
            "B_into_HH={} I={} S={} B={} I_next={} exact=true",
            ranks[0], ranks[1], ranks[2], ranks[3], ranks[4]
        ))
    }
}

fn image(
    model: &mut Model,
    source: &Complex,
    target: &Complex,
    operation: Operation,
) -> Result<(Vec<Chain>, usize)> {
    let representatives = source.representatives(model)?;
    let mut matrix = linear::matrix(representatives.len(), target.axis.len())?;
    let mut images = Vec::new();
    let mut cells = 0;
    for (representative, row) in representatives.iter().zip(&mut matrix) {
        let raw = raw_map(model, operation, representative)?;
        let value = target.project(model, &raw)?;
        *row = target.axis.coordinates(&value)?;
        cells += value.cells + 8;
        storage(cells)?;
        images.push(value);
    }
    let rank = linear::rref(model, matrix, target.axis.len())?.pivots.len();
    Ok((images, rank))
}

fn lift(
    model: &mut Model,
    algebra: &Algebra,
    grade: &Grade,
    target: &Grade,
    value: &Chain,
) -> Result<Chain> {
    let negative = Complex::new(model, algebra, target)?;
    let hochschild = Complex::new(model, algebra, grade)?;
    let representatives = negative.representatives(model)?;
    let mut matrix = linear::matrix(hochschild.axis.len(), representatives.len())?;
    for (j, representative) in representatives.iter().enumerate() {
        let evaluated = raw_map(model, Operation::Eval, representative)?;
        let evaluated = hochschild.project(model, &evaluated)?;
        let column = hochschild.axis.coordinates(&evaluated)?;
        for (row, coordinate) in matrix.iter_mut().zip(column) {
            row[j] = coordinate;
        }
    }
    let rhs = hochschild.axis.coordinates(value)?;
    let solution = linear::solve(model, matrix, representatives.len(), &rhs)?.ok_or_else(|| {
        Error::plain("Verbitsky Hochschild class has no full negative cyclic lift")
    })?;
    let result = combination(model, &representatives, &solution)?;
    let evaluated = raw_map(model, Operation::Eval, &result)?;
    if !hochschild
        .project(model, &evaluated)?
        .compare(model, value)?
        .is_eq()
    {
        return Err(Error::plain(
            "Verbitsky negative cyclic lift failed verification",
        ));
    }
    Ok(result)
}

fn combination(model: &mut Model, values: &[Chain], coefficients: &[Rational]) -> Result<Chain> {
    let mut result = Chain::default();
    for (value, &coefficient) in values.iter().zip(coefficients) {
        for (key, &entry) in &value.entries {
            let product = entry.multiply(model, coefficient)?;
            result.insert(model, key.clone(), product)?;
        }
    }
    Ok(result)
}
