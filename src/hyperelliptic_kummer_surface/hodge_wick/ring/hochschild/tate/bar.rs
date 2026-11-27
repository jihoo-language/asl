use super::super::{linear, Algebra, Bilinear, Cochain, Form, Grade as BarGrade, Key};
use super::{grading_range, storage, Definition, Error, Grade, Model, Rational, Result, State};

type Coordinates = Vec<(usize, Rational)>;

struct Comparison {
    powers: Vec<usize>,
    strides: Vec<usize>,
    exponents: Vec<Vec<usize>>,
    forward: Vec<Coordinates>,
    inverse: Vec<Coordinates>,
}

impl Comparison {
    fn new(model: &mut Model, algebra: &Algebra, powers: &[usize]) -> Result<Self> {
        let mut size = 1usize;
        let mut strides = Vec::new();
        for &power in powers {
            strides.push(size);
            size = size.checked_mul(power).ok_or_else(incompatible)?;
        }
        let variables = algebra.weight_range(1);
        if size != algebra.len() || variables.len() != powers.len() {
            return Err(incompatible());
        }
        let mut forward = Vec::with_capacity(size);
        let mut exponents = Vec::with_capacity(size);
        let mut matrix = linear::matrix(size, 2 * size)?;
        forward.push(vec![(0, Rational::integer(1))]);
        exponents.push(vec![0; powers.len()]);
        for index in 1..size {
            let coordinates = strides
                .iter()
                .zip(powers)
                .map(|(&stride, &power)| index / stride % power)
                .collect::<Vec<_>>();
            let variable = coordinates
                .iter()
                .position(|&exponent| exponent > 0)
                .unwrap();
            let product = multiply_basis(
                model,
                algebra,
                &forward[index - strides[variable]],
                variables.start + variable,
            )?;
            exponents.push(coordinates);
            forward.push(product);
        }
        for (i, (&power, &stride)) in powers.iter().zip(&strides).enumerate() {
            let before_relation = (power - 1) * stride;
            if !multiply_basis(
                model,
                algebra,
                &forward[before_relation],
                variables.start + i,
            )?
            .is_empty()
            {
                return Err(incompatible());
            }
        }
        for (column, values) in forward.iter().enumerate() {
            matrix[column][size + column] = Rational::integer(1);
            for &(row, coefficient) in values {
                matrix[row][column] = coefficient;
            }
        }
        let matrix = linear::rref(model, matrix, 2 * size)?;
        if matrix.pivots.iter().copied().ne(0..size) {
            return Err(incompatible());
        }
        let mut inverse = vec![Vec::new(); size];
        for (i, row) in matrix.rows.into_iter().enumerate() {
            for (j, coefficient) in row.into_iter().skip(size).enumerate() {
                if !coefficient.is_zero() {
                    inverse[j].push((i, coefficient));
                }
            }
        }
        Ok(Self {
            powers: powers.to_vec(),
            strides,
            exponents,
            forward,
            inverse,
        })
    }

    fn monomial_index(&self, factors: &[super::Factor]) -> usize {
        factors
            .iter()
            .zip(&self.strides)
            .map(|(factor, stride)| factor.exponent as usize * stride)
            .sum()
    }

    fn constant(
        &self,
        model: &mut Model,
        grade: &BarGrade,
        index: usize,
        shift: i128,
    ) -> Result<Representative> {
        let mut value = Cochain::default();
        for &(output, coefficient) in &self.forward[index] {
            value.insert(
                model,
                Key {
                    word: Vec::new(),
                    output,
                },
                coefficient,
            )?;
        }
        let grade = BarGrade {
            arity: 0,
            shift,
            ..grade.clone()
        };
        Ok(Representative { grade, value })
    }

    fn euler(
        &self,
        model: &mut Model,
        grade: &BarGrade,
        variable: usize,
    ) -> Result<Representative> {
        let mut value = Cochain::default();
        for (input, inverse) in self.inverse.iter().enumerate().skip(1) {
            for &(monomial, coefficient) in inverse {
                let exponent = self.exponents[monomial][variable];
                if exponent == 0 {
                    continue;
                }
                let coefficient =
                    coefficient.multiply(model, Rational::integer(exponent as i128))?;
                for &(output, scalar) in &self.forward[monomial] {
                    let entry = coefficient.multiply(model, scalar)?;
                    value.insert(
                        model,
                        Key {
                            word: vec![input],
                            output,
                        },
                        entry,
                    )?;
                }
            }
        }
        let grade = BarGrade {
            arity: 1,
            shift: 0,
            ..grade.clone()
        };
        Ok(Representative { grade, value })
    }

    fn carry(
        &self,
        model: &mut Model,
        grade: &BarGrade,
        variable: usize,
    ) -> Result<Representative> {
        let mut value = Cochain::default();
        let d = self.powers[variable];
        for (i, left) in self.inverse.iter().enumerate().skip(1) {
            for (j, right) in self.inverse.iter().enumerate().skip(1) {
                model.spend(1)?;
                for &(a, f) in left {
                    for &(b, g) in right {
                        model.spend(self.powers.len())?;
                        if self.exponents[a][variable] + self.exponents[b][variable] < d {
                            continue;
                        }
                        let vanishes = self.exponents[a]
                            .iter()
                            .zip(&self.exponents[b])
                            .zip(&self.powers)
                            .enumerate()
                            .any(|(index, ((&a, &b), &power))| index != variable && a + b >= power);
                        if vanishes {
                            continue;
                        }
                        let monomial = a + b - d * self.strides[variable];
                        let coefficient = f.multiply(model, g)?;
                        for &(output, scalar) in &self.forward[monomial] {
                            let entry = coefficient.multiply(model, scalar)?;
                            value.insert(
                                model,
                                Key {
                                    word: vec![i, j],
                                    output,
                                },
                                entry,
                            )?;
                        }
                    }
                }
            }
        }
        let grade = BarGrade {
            arity: 2,
            shift: -(d as i128),
            ..grade.clone()
        };
        Ok(Representative { grade, value })
    }
}

fn multiply_basis(
    model: &mut Model,
    algebra: &Algebra,
    value: &Coordinates,
    basis: usize,
) -> Result<Coordinates> {
    let mut sums = vec![Rational::zero(); algebra.len()];
    for &(index, coefficient) in value {
        for &(k, scalar) in algebra.product(index, basis) {
            let term = coefficient.multiply(model, scalar)?;
            sums[k] = sums[k].add(model, term)?;
        }
    }
    Ok(sums
        .into_iter()
        .enumerate()
        .filter(|(_, value)| !value.is_zero())
        .collect())
}

fn incompatible() -> Error {
    Error::plain("the marked divisor axes do not present the named monomial complete intersection")
}

struct Representative {
    grade: BarGrade,
    value: Cochain,
}

impl Representative {
    fn cup(&self, model: &mut Model, form: &Form, right: &Self) -> Result<Self> {
        let grade = self.grade.product_grade(&right.grade, Bilinear::Cup)?;
        let value =
            form.hochschild_bilinear(model, &grade, Bilinear::Cup, &self.value, &right.value)?;
        Ok(Self { grade, value })
    }

    fn power(
        &self,
        model: &mut Model,
        form: &Form,
        comparison: &Comparison,
        mut exponent: u128,
    ) -> Result<Self> {
        let mut result = comparison.constant(model, &self.grade, 0, 0)?;
        let mut squared = None;
        while exponent > 0 {
            let factor = squared.as_ref().unwrap_or(self);
            if !exponent.is_multiple_of(2) {
                result = result.cup(model, form, factor)?;
            }
            exponent /= 2;
            if exponent > 0 {
                squared = Some(factor.cup(model, form, factor)?);
            }
        }
        Ok(result)
    }
}

impl Definition {
    pub(crate) fn bar_class(
        &self,
        model: &mut Model,
        form: &Form,
        grade: &Grade,
        target: &BarGrade,
        input: &State,
    ) -> Result<Cochain> {
        grade.bar_grade(target.basis.clone(), target.length)?;
        if !input.terms.is_empty() {
            let arity = usize::try_from(target.arity).map_err(|_| grading_range())?;
            storage(arity.saturating_add(1).saturating_mul(4))?;
        }
        let algebra = form.hochschild_algebra(model, target)?;
        let comparison = Comparison::new(model, &algebra, &self.powers)?;
        let mut eulers = Vec::new();
        let mut carries = Vec::new();
        for i in 0..self.powers.len() {
            let need_euler = input.terms.keys().any(|key| key[i].exterior);
            let need_carry = input.terms.keys().any(|key| key[i].periodic > 0);
            eulers.push(
                need_euler
                    .then(|| comparison.euler(model, target, i))
                    .transpose()?,
            );
            carries.push(
                need_carry
                    .then(|| comparison.carry(model, target, i))
                    .transpose()?,
            );
        }
        let mut sum = Cochain::default();
        for (key, &coefficient) in &input.terms {
            let index = comparison.monomial_index(key);
            let weight = key.iter().map(|factor| factor.exponent as i128).sum();
            let mut current = comparison.constant(model, target, index, weight)?;
            for (i, factor) in key.iter().enumerate() {
                if factor.exterior {
                    current = current.cup(model, form, eulers[i].as_ref().unwrap())?;
                }
                if factor.periodic > 0 {
                    let power = carries[i].as_ref().unwrap().power(
                        model,
                        form,
                        &comparison,
                        factor.periodic as u128,
                    )?;
                    current = current.cup(model, form, &power)?;
                }
            }
            if current.grade != *target {
                return Err(Error::plain(
                    "normalized bar comparison produced the wrong grade",
                ));
            }
            for (key, value) in current.value.entries {
                let value = value.multiply(model, coefficient)?;
                sum.insert(model, key, value)?;
            }
        }
        form.hochschild_class(model, target, &sum)
    }
}
