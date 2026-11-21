use crate::binary_field::{BinaryField, Embedding};
use crate::binary_linear::{add_multiple, check_storage, Basis, Entry, Vector};
use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::rv_descent::{self, Matrix};
use crate::rv_hom::Grade;

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Norm,
    Companion,
    Frame,
    Factors,
}

impl Operation {
    fn name(self) -> &'static str {
        match self {
            Self::Norm => "rvFrobNorm",
            Self::Companion => "rvShintani",
            Self::Frame => "rvShintaniFrame",
            Self::Factors => "rvInvariantFactors",
        }
    }
}

pub(crate) fn plan(
    model: &mut Model,
    ty: TypeId,
    modulus: usize,
    root: usize,
    operation: Operation,
) -> Result<(TypeId, Embedding)> {
    let (grade, embedding, point) =
        rv_descent::endomorphism_point(model, ty, modulus, root, operation.name())?;
    let output = if matches!(operation, Operation::Companion) {
        model.allocate(1 + point.len())?;
        model.rv_hom_type(Grade {
            field: embedding.small(),
            point: point.into(),
            ..grade
        })?
    } else {
        ty
    };
    Ok((output, embedding))
}

fn norm(
    model: &mut Model,
    ty: TypeId,
    embedding: &Embedding,
    support: &[Entry],
    operation: Operation,
) -> Result<Matrix> {
    let field = embedding.large();
    let step = embedding.small().degree();
    let degree = field.degree() / step;
    let size = model.rv_hom(ty)?.dimensions().0;
    let input = Matrix::read(model, size, size, support, |_, coefficient| Ok(coefficient))?;
    let mut output = Matrix::new(size, size)?;
    let mut basis = Basis::default();
    for column in 0..size {
        let mut image = vec![(column, 1)];
        for _ in 0..degree {
            image = input.semilinear_step(model, field, step, &image)?;
        }
        model.spend(image.len())?;
        basis.insert(model, field, image.clone())?;
        if basis.rank() != column + 1 {
            return Err(Error::plain(format!(
                "{} requires an invertible Margolis endomorphism",
                operation.name()
            )));
        }
        output.push(image)?;
    }
    Ok(output)
}

pub(crate) fn evaluate(
    model: &mut Model,
    ty: TypeId,
    embedding: &Embedding,
    support: &[Entry],
    operation: Operation,
) -> Result<Value> {
    let norm = norm(model, ty, embedding, support, operation)?;
    match operation {
        Operation::Norm => norm.value(model),
        Operation::Companion => {
            let canonical = Smith::compute(model, embedding, &norm, false)?;
            canonical.companion(model)?.value(model)
        }
        Operation::Frame => {
            let canonical = Smith::compute(model, embedding, &norm, true)?;
            canonical
                .frame(model, embedding.large(), &norm)?
                .value(model)
        }
        Operation::Factors => unreachable!(),
    }
}

pub(crate) fn invariant_factors(
    model: &mut Model,
    ty: TypeId,
    embedding: &Embedding,
    support: &[Entry],
) -> Result<Vec<Polynomial>> {
    let norm = norm(model, ty, embedding, support, Operation::Factors)?;
    Ok(Smith::compute(model, embedding, &norm, false)?.factors)
}

// Polynomials have coefficients in ascending order, with trailing zeros removed.
type Polynomial = Vec<usize>;

fn trim(polynomial: &mut Polynomial) {
    while polynomial.last() == Some(&0) {
        polynomial.pop();
    }
}

fn divide(
    model: &mut Model,
    field: BinaryField,
    numerator: &[usize],
    denominator: &[usize],
) -> Result<(Polynomial, Polynomial)> {
    model.spend(1 + numerator.len())?;
    if numerator.len() < denominator.len() {
        return Ok((Vec::new(), numerator.to_vec()));
    }
    let quotient_len = numerator.len() - denominator.len() + 1;
    check_storage(quotient_len + numerator.len())?;
    let leading = *denominator.last().unwrap();
    let inverse = if leading == 1 {
        1
    } else {
        field.inverse(model, leading)?
    };
    let mut quotient = vec![0; quotient_len];
    let mut remainder = numerator.to_vec();
    while remainder.len() >= denominator.len() {
        let offset = remainder.len() - denominator.len();
        model.spend(1 + field.degree() * (1 + denominator.len()))?;
        let coefficient = field.multiply(*remainder.last().unwrap(), inverse);
        quotient[offset] = coefficient;
        for (entry, &value) in remainder[offset..].iter_mut().zip(denominator) {
            *entry ^= field.multiply(coefficient, value);
        }
        trim(&mut remainder);
    }
    Ok((quotient, remainder))
}

fn add_product(
    model: &mut Model,
    field: BinaryField,
    target: &[usize],
    factor: &[usize],
    source: &[usize],
) -> Result<Polynomial> {
    if source.is_empty() || factor.is_empty() {
        model.spend(1 + target.len())?;
        return Ok(target.to_vec());
    }
    let size = target.len().max(factor.len() + source.len() - 1);
    check_storage(size)?;
    model.spend(1 + target.len() + factor.len())?;
    let mut result = target.to_vec();
    result.resize(size, 0);
    for (offset, &coefficient) in factor.iter().enumerate() {
        if coefficient == 0 {
            continue;
        }
        model.spend(field.degree() * source.len())?;
        for (entry, &value) in result[offset..].iter_mut().zip(source) {
            *entry ^= field.multiply(coefficient, value);
        }
    }
    trim(&mut result);
    Ok(result)
}

struct Canonical {
    size: usize,
    factors: Vec<Polynomial>,
    cyclic_vectors: Option<Vec<Vector>>,
}

impl Canonical {
    fn companion(&self, model: &mut Model) -> Result<Matrix> {
        let mut output = Matrix::new(self.size, self.size)?;
        let mut offset = 0;
        for factor in &self.factors {
            let degree = factor.len() - 1;
            model.spend(2 * degree)?;
            for row in 1..degree {
                output.push(vec![(offset + row, 1)])?;
            }
            let last = factor[..degree]
                .iter()
                .enumerate()
                .filter_map(|(row, &coefficient)| {
                    (coefficient != 0).then_some((offset + row, coefficient))
                })
                .collect();
            output.push(last)?;
            offset += degree;
        }
        Ok(output)
    }

    fn frame(self, model: &mut Model, field: BinaryField, norm: &Matrix) -> Result<Matrix> {
        let mut output = Matrix::new(self.size, self.size)?;
        let vectors = self.cyclic_vectors.unwrap();
        for (factor, mut vector) in self.factors.into_iter().zip(vectors) {
            let degree = factor.len() - 1;
            for _ in 1..degree {
                let next = norm.apply(model, field, &vector)?;
                output.push(vector)?;
                vector = next;
            }
            output.push(vector)?;
        }
        Ok(output)
    }
}

struct Smith<'a> {
    field: BinaryField,
    norm: &'a Matrix,
    table: Vec<Vec<Polynomial>>,
    cells: usize,
    generators: Option<Vec<Vector>>,
    generator_cells: usize,
}

impl<'a> Smith<'a> {
    fn new(model: &mut Model, field: BinaryField, norm: &'a Matrix, track: bool) -> Result<Self> {
        let size = norm.columns().len();
        let cells = size
            .checked_mul(size)
            .ok_or_else(|| Error::plain("Shintani polynomial matrix exceeds the storage limit"))?;
        let generator_cells = if track { 3 * size } else { 0 };
        check_storage(cells + generator_cells)?;
        model.spend(cells + generator_cells)?;
        let mut smith = Self {
            field,
            norm,
            table: vec![vec![Vec::new(); size]; size],
            cells,
            generators: track.then(|| (0..size).map(|index| vec![(index, 1)]).collect()),
            generator_cells,
        };
        for (column, entries) in norm.columns().iter().enumerate() {
            model.spend(1 + entries.len())?;
            for &(row, coefficient) in entries {
                smith.write(row, column, vec![coefficient])?;
            }
            let constant = smith.table[column][column].first().copied().unwrap_or(0);
            smith.write(column, column, vec![constant, 1])?;
        }
        Ok(smith)
    }

    fn write(&mut self, row: usize, column: usize, value: Polynomial) -> Result<()> {
        self.cells += value.len();
        self.cells -= self.table[row][column].len();
        check_storage(self.cells + self.generator_cells)?;
        self.table[row][column] = value;
        Ok(())
    }

    fn swap_rows(&mut self, first: usize, second: usize) {
        self.table.swap(first, second);
        if let Some(generators) = &mut self.generators {
            generators.swap(first, second);
        }
    }

    fn swap_columns(&mut self, model: &mut Model, first: usize, second: usize) -> Result<()> {
        if first != second {
            model.spend(self.table.len())?;
            for row in &mut self.table {
                row.swap(first, second);
            }
        }
        Ok(())
    }

    fn row_add(
        &mut self,
        model: &mut Model,
        target: usize,
        source: usize,
        factor: &[usize],
    ) -> Result<()> {
        for column in 0..self.table.len() {
            let value = add_product(
                model,
                self.field,
                &self.table[target][column],
                factor,
                &self.table[source][column],
            )?;
            self.write(target, column, value)?;
        }
        if let Some(generators) = &mut self.generators {
            // Row_target += factor(T) Row_source acts inversely on the
            // generator of coordinate source. The characteristic is two.
            let mut image = Vec::new();
            for &coefficient in factor.iter().rev() {
                model.spend(1)?;
                image = self.norm.apply(model, self.field, &image)?;
                image = add_multiple(model, self.field, &image, &generators[target], coefficient)?;
            }
            let updated = add_multiple(model, self.field, &generators[source], &image, 1)?;
            self.generator_cells += 2 * updated.len();
            self.generator_cells -= 2 * generators[source].len();
            check_storage(self.cells + self.generator_cells)?;
            generators[source] = updated;
        }
        Ok(())
    }

    fn column_add(
        &mut self,
        model: &mut Model,
        target: usize,
        source: usize,
        factor: &[usize],
    ) -> Result<()> {
        for row in 0..self.table.len() {
            let value = add_product(
                model,
                self.field,
                &self.table[row][target],
                factor,
                &self.table[row][source],
            )?;
            self.write(row, target, value)?;
        }
        Ok(())
    }

    fn pivot(&mut self, model: &mut Model, index: usize) -> Result<()> {
        let size = self.table.len();
        model.spend((size - index) * (size - index))?;
        let mut next = None;
        for (row, entries) in self.table.iter().enumerate().skip(index) {
            for (column, polynomial) in entries.iter().enumerate().skip(index) {
                let len = polynomial.len();
                if len != 0 {
                    let candidate = (len, row, column);
                    if next.is_none_or(|previous| candidate < previous) {
                        next = Some(candidate);
                    }
                }
            }
        }
        let (_, row, column) =
            next.ok_or_else(|| Error::plain("Shintani polynomial matrix lost rank"))?;
        self.swap_rows(index, row);
        self.swap_columns(model, index, column)?;

        'reduce: loop {
            model.spend(1 + 2 * (size - index))?;
            for row in index + 1..size {
                if self.table[row][index].is_empty() {
                    continue;
                }
                let (quotient, remainder) = divide(
                    model,
                    self.field,
                    &self.table[row][index],
                    &self.table[index][index],
                )?;
                if !quotient.is_empty() {
                    self.row_add(model, row, index, &quotient)?;
                }
                if !remainder.is_empty() {
                    self.swap_rows(index, row);
                    continue 'reduce;
                }
            }
            for column in index + 1..size {
                if self.table[index][column].is_empty() {
                    continue;
                }
                let (quotient, remainder) = divide(
                    model,
                    self.field,
                    &self.table[index][column],
                    &self.table[index][index],
                )?;
                if !quotient.is_empty() {
                    self.column_add(model, column, index, &quotient)?;
                }
                if !remainder.is_empty() {
                    self.swap_columns(model, index, column)?;
                    continue 'reduce;
                }
            }

            // Once the pivot row and column are clear, force the pivot to
            // divide every entry of the remaining block before proceeding.
            // Adding such a row exposes a smaller Euclidean remainder.
            if self.table[index][index].len() > 1 {
                for row in index + 1..size {
                    for column in index + 1..size {
                        let (_, remainder) = divide(
                            model,
                            self.field,
                            &self.table[row][column],
                            &self.table[index][index],
                        )?;
                        if !remainder.is_empty() {
                            self.row_add(model, index, row, &[1])?;
                            continue 'reduce;
                        }
                    }
                }
            }
            break;
        }

        let diagonal = &mut self.table[index][index];
        let leading = *diagonal.last().unwrap();
        if leading != 1 {
            let inverse = self.field.inverse(model, leading)?;
            model.spend(self.field.degree() * diagonal.len())?;
            for coefficient in diagonal {
                *coefficient = self.field.multiply(*coefficient, inverse);
            }
            if let Some(generators) = &mut self.generators {
                let generator = &mut generators[index];
                model.spend(self.field.degree() * generator.len())?;
                for (_, coefficient) in generator {
                    *coefficient = self.field.multiply(*coefficient, leading);
                }
            }
        }
        Ok(())
    }

    fn compute(
        model: &mut Model,
        embedding: &Embedding,
        norm: &'a Matrix,
        track: bool,
    ) -> Result<Canonical> {
        let size = norm.columns().len();
        let mut smith = Self::new(model, embedding.large(), norm, track)?;
        for index in 0..size {
            smith.pivot(model, index)?;
        }
        let mut factors = Vec::new();
        let mut cyclic_vectors = track.then(Vec::new);
        let mut total_degree = 0;
        for (index, row) in smith.table.into_iter().enumerate() {
            let factor = row.into_iter().nth(index).unwrap();
            if factor.len() < 2 {
                continue;
            }
            total_degree += factor.len() - 1;
            let factor = factor
                .into_iter()
                .map(|coefficient| {
                    embedding.preimage(model, coefficient)?.ok_or_else(|| {
                        Error::plain(
                            "Shintani invariant factor does not lie in the specified subfield",
                        )
                    })
                })
                .collect::<Result<Polynomial>>()?;
            factors.push(factor);
            if let (Some(generators), Some(vectors)) = (&mut smith.generators, &mut cyclic_vectors)
            {
                vectors.push(std::mem::take(&mut generators[index]));
            }
        }
        if total_degree != size {
            return Err(Error::plain(
                "Shintani invariant degrees do not match the matrix dimension",
            ));
        }
        Ok(Canonical {
            size,
            factors,
            cyclic_vectors,
        })
    }
}
