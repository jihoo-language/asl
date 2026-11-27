use super::mixed::differential_terms;
use super::{
    basis_error, linear, storage, work_limit, Algebra, Chain, Error, Flavor, Grade, Key, Model,
    Rational, Result,
};

#[derive(Default)]
pub(super) struct Axis(Vec<Key>);

impl Axis {
    fn new(model: &mut Model, algebra: &Algebra, grade: &Grade) -> Result<Self> {
        if !grade.possible() {
            return Ok(Self::default());
        }
        let mut min = match grade.flavor {
            Flavor::Hochschild => grade.degree,
            Flavor::Negative => 0.max(grade.degree),
            _ => 0,
        };
        let max = match grade.flavor {
            Flavor::Hochschild => grade.degree,
            Flavor::Cyclic => grade.degree.min(grade.weight),
            _ => grade.weight,
        };
        let top = algebra.top() as i128;
        // The leading factor can have weight zero, and no factor exceeds top.
        let min_by_weight = grade.weight / top - i128::from(grade.weight % top == 0);
        min = min.max(min_by_weight);
        if min > max {
            return Ok(Self::default());
        }
        let start = usize::try_from(min).map_err(|_| work_limit())?;
        let end = usize::try_from(max).map_err(|_| work_limit())?;
        let count = (end - start).checked_add(1).ok_or_else(work_limit)?;
        storage(count)?;
        model.spend(count)?;
        let mut cells = 0;
        let mut keys = Vec::new();
        for arity in start..=end {
            let Some(power) = grade.power(arity) else {
                continue;
            };
            let factors = arity.checked_add(1).ok_or_else(work_limit)?;
            storage(factors.saturating_mul(4))?;
            append_words(
                model,
                algebra,
                factors,
                grade.weight,
                power,
                &mut keys,
                &mut cells,
            )?;
        }
        Ok(Self(keys))
    }

    pub(super) fn len(&self) -> usize {
        self.0.len()
    }

    fn locate(&self, key: &Key) -> Result<usize> {
        self.0.binary_search(key).map_err(|_| {
            Error::plain("Verbitsky chain coefficient lies outside its graded complex")
        })
    }

    pub(super) fn coordinates(&self, value: &Chain) -> Result<Vec<Rational>> {
        let mut row = vec![Rational::zero(); self.len()];
        for (key, &coefficient) in &value.entries {
            row[self.locate(key)?] = coefficient;
        }
        Ok(row)
    }

    fn chain(&self, model: &mut Model, values: &[Rational]) -> Result<Chain> {
        let mut result = Chain::default();
        for (key, &coefficient) in self.0.iter().zip(values) {
            if !coefficient.is_zero() {
                result.insert(model, key.clone(), coefficient)?;
            }
        }
        Ok(result)
    }
}

fn append_words(
    model: &mut Model,
    algebra: &Algebra,
    factors: usize,
    weight: i128,
    power: i128,
    keys: &mut Vec<Key>,
    cells: &mut usize,
) -> Result<()> {
    let mut word = vec![0; factors];
    let mut prefix = vec![0i128; factors + 1];
    let mut next = vec![1; factors];
    next[0] = 0;
    let mut depth = 0;
    loop {
        if depth == factors {
            *cells = cells.saturating_add(factors + 10);
            storage(*cells)?;
            model.spend(factors)?;
            keys.push(Key {
                power,
                word: word.clone(),
            });
            depth -= 1;
            continue;
        }
        let rem = (factors - depth - 1) as i128;
        let mut found = false;
        for i in next[depth]..algebra.len() {
            model.spend(1)?;
            let current = prefix[depth] + algebra.weights[i] as i128;
            if current + rem * (algebra.top() as i128) < weight {
                continue;
            }
            if current + rem > weight {
                break;
            }
            prefix[depth + 1] = current;
            word[depth] = i;
            next[depth] = i + 1;
            depth += 1;
            if depth < factors {
                next[depth] = 1;
            }
            found = true;
            break;
        }
        if !found {
            if depth == 0 {
                break;
            }
            depth -= 1;
        }
    }
    Ok(())
}

fn differential(
    model: &mut Model,
    algebra: &Algebra,
    flavor: Flavor,
    from: &Axis,
    to: &Axis,
    mut emit: impl FnMut(&mut Model, usize, usize, Rational) -> Result<()>,
) -> Result<()> {
    for (column, key) in from.0.iter().enumerate() {
        differential_terms(model, algebra, flavor, key, |model, term, coefficient| {
            emit(model, to.locate(&term)?, column, coefficient)
        })?;
    }
    Ok(())
}

pub(super) struct Complex {
    pub(super) axis: Axis,
    cycles: linear::Echelon,
    boundaries: linear::Echelon,
}

impl Complex {
    pub(super) fn new(model: &mut Model, algebra: &Algebra, grade: &Grade) -> Result<Self> {
        let axis = Axis::new(model, algebra, grade)?;
        if axis.len() == 0 {
            return Ok(Self {
                axis,
                cycles: linear::rref(model, Vec::new(), 0)?,
                boundaries: linear::rref(model, Vec::new(), 0)?,
            });
        }
        let next = Axis::new(model, algebra, &grade.with_offset(-1)?)?;
        let prev = Axis::new(model, algebra, &grade.with_offset(1)?)?;
        let mut outgoing = linear::matrix(next.len(), axis.len())?;
        differential(
            model,
            algebra,
            grade.flavor,
            &axis,
            &next,
            |model, i, j, coefficient| {
                outgoing[i][j] = outgoing[i][j].add(model, coefficient)?;
                Ok(())
            },
        )?;
        let cycles = linear::rref(model, outgoing, axis.len())?;
        let mut incoming = linear::matrix(prev.len(), axis.len())?;
        differential(
            model,
            algebra,
            grade.flavor,
            &prev,
            &axis,
            |model, i, j, coefficient| {
                incoming[j][i] = incoming[j][i].add(model, coefficient)?;
                Ok(())
            },
        )?;
        let boundaries = linear::rref(model, incoming, axis.len())?;
        for row in &boundaries.rows {
            if !super::super::complex::annihilates(model, &cycles, row)? {
                return Err(Error::plain(
                    "Verbitsky mixed complex failed its square-zero check",
                ));
            }
        }
        Ok(Self {
            axis,
            cycles,
            boundaries,
        })
    }

    pub(super) fn dimensions(&self) -> [usize; 4] {
        let chains = self.axis.len();
        let cycles = chains - self.cycles.rows.len();
        let boundaries = self.boundaries.rows.len();
        [chains, cycles, boundaries, cycles - boundaries]
    }

    fn reduce(&self, model: &mut Model, coordinates: &mut [Rational]) -> Result<()> {
        for (&pivot, row) in self.boundaries.pivots.iter().zip(&self.boundaries.rows) {
            linear::add_multiple(model, coordinates, row, coordinates[pivot].negate())?;
        }
        Ok(())
    }

    pub(super) fn project(&self, model: &mut Model, value: &Chain) -> Result<Chain> {
        let mut coordinates = self.axis.coordinates(value)?;
        if !super::super::complex::annihilates(model, &self.cycles, &coordinates)? {
            return Err(Error::plain("Verbitsky cyclic chain is not a cycle"));
        }
        self.reduce(model, &mut coordinates)?;
        self.axis.chain(model, &coordinates)
    }

    pub(super) fn basis(&self, model: &mut Model, index: usize) -> Result<Chain> {
        if index >= self.dimensions()[3] {
            return Err(basis_error());
        }
        let basis = self.canonical_rows(model)?;
        self.axis.chain(model, &basis[index])
    }

    fn canonical_rows(&self, model: &mut Model) -> Result<linear::Matrix> {
        let mut cycles = self.cycles.kernel(self.axis.len())?;
        for row in &mut cycles {
            self.reduce(model, row)?;
        }
        let basis = linear::rref(model, cycles, self.axis.len())?.rows;
        if basis.len() != self.dimensions()[3] {
            return Err(Error::plain(
                "Verbitsky cyclic quotient has unexpected rank",
            ));
        }
        Ok(basis)
    }

    pub(super) fn representatives(&self, model: &mut Model) -> Result<Vec<Chain>> {
        let rows = self.canonical_rows(model)?;
        let mut cells = 0;
        rows.iter()
            .map(|row| {
                let chain = self.axis.chain(model, row)?;
                cells += chain.cells + 8;
                storage(cells)?;
                Ok(chain)
            })
            .collect()
    }

    pub(super) fn primitive(
        &self,
        model: &mut Model,
        algebra: &Algebra,
        grade: &Grade,
        boundary: &Chain,
    ) -> Result<Chain> {
        let from = Axis::new(model, algebra, &grade.with_offset(1)?)?;
        let mut matrix = linear::matrix(self.axis.len(), from.len())?;
        differential(
            model,
            algebra,
            grade.flavor,
            &from,
            &self.axis,
            |model, i, j, coefficient| {
                matrix[i][j] = matrix[i][j].add(model, coefficient)?;
                Ok(())
            },
        )?;
        let rhs = self.axis.coordinates(boundary)?;
        let solution = linear::solve(model, matrix, from.len(), &rhs)?
            .ok_or_else(|| Error::plain("Verbitsky cyclic chain has no boundary primitive"))?;
        from.chain(model, &solution)
    }
}
