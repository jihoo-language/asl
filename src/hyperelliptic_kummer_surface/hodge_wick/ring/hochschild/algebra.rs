use crate::error::Error;

use super::{
    linear, moments, storage, work_limit, Form, Grade, Key, Model, Monomial, Rational, Rc, RefCell,
    Result, Tensor, LIMBS,
};

pub(crate) struct Algebra {
    pub(super) weights: Vec<usize>,
    pub(super) bases: Vec<Monomial>,
    starts: Vec<usize>,
    product: Vec<Vec<(usize, Rational)>>,
    splits: Vec<Vec<(usize, usize, Rational)>>,
    trace: Rational,
    dual: RefCell<Option<Rc<linear::Matrix>>>,
}

impl Algebra {
    pub(super) fn len(&self) -> usize {
        self.weights.len()
    }

    pub(super) fn top(&self) -> usize {
        self.starts.len() - 2
    }

    pub(super) fn weight_range(&self, weight: usize) -> std::ops::Range<usize> {
        self.starts[weight]..self.starts[weight + 1]
    }

    pub(super) fn product(&self, left: usize, right: usize) -> &[(usize, Rational)] {
        &self.product[left * self.len() + right]
    }

    pub(super) fn coefficient_indices(
        &self,
        weight: usize,
        tensor: &Tensor,
    ) -> Result<Vec<(usize, Rational)>> {
        let indices = self.weight_range(weight);
        tensor
            .terms
            .iter()
            .map(|(powers, &coefficient)| {
                let k = self.bases[indices.clone()]
                    .binary_search(powers)
                    .map_err(|_| {
                        Error::plain("Verbitsky product is outside its canonical basis")
                    })?;
                Ok((indices.start + k, coefficient))
            })
            .collect()
    }

    pub(super) fn splits(&self, output: usize) -> &[(usize, usize, Rational)] {
        &self.splits[output]
    }

    pub(super) fn validate_key(&self, grade: &Grade, key: &Key) -> Result<()> {
        if grade.arity != key.word.len() as i128 {
            return Err(Error::plain("Hochschild coordinate has the wrong arity"));
        }
        if key.output >= self.len() || key.word.iter().any(|&a| a >= self.len()) {
            return Err(Error::plain(
                "Hochschild coordinate index is outside the ring basis",
            ));
        }
        if key.word.contains(&0) {
            return Err(Error::plain(
                "normalized Hochschild inputs must have positive weight",
            ));
        }
        let sum: i128 = key.word.iter().map(|&a| self.weights[a] as i128).sum();
        if sum.checked_add(grade.shift) != Some(self.weights[key.output] as i128) {
            return Err(Error::plain(
                "Hochschild coordinate has the wrong internal shift",
            ));
        }
        Ok(())
    }

    pub(super) fn pairing_inverse(&self, model: &mut Model) -> Result<Rc<linear::Matrix>> {
        if let Some(inverse) = self.dual.borrow().as_ref() {
            return Ok(Rc::clone(inverse));
        }
        let size = self.len();
        let mut gram = linear::matrix(size, 2 * size)?;
        for (i, row) in gram.iter_mut().enumerate() {
            row[size + i] = Rational::integer(1);
            for j in self.weight_range(self.top() - self.weights[i]) {
                let value = self
                    .product(i, j)
                    .first()
                    .map_or_else(Rational::zero, |p| p.1);
                row[j] = value.multiply(model, self.trace)?;
            }
        }
        let echelon = linear::rref(model, gram, 2 * size)?;
        if echelon.pivots.iter().copied().ne(0..size) {
            return Err(Error::plain(
                "Verbitsky Hochschild Frobenius pairing is singular",
            ));
        }
        let inverse = Rc::new(
            echelon
                .rows
                .into_iter()
                .map(|mut row| row.split_off(size))
                .collect(),
        );
        model.allocate(size * size * (2 * LIMBS + 1))?;
        *self.dual.borrow_mut() = Some(Rc::clone(&inverse));
        Ok(inverse)
    }

    pub(super) fn top_index(&self) -> usize {
        self.len() - 1
    }

    pub(super) fn trace_factor(&self) -> Rational {
        self.trace
    }
}

impl Form {
    pub(super) fn hochschild_algebra(
        &self,
        model: &mut Model,
        grade: &Grade,
    ) -> Result<Rc<Algebra>> {
        self.hochschild_validate(grade)?;
        if let Some(value) = self.hochschild.borrow().get(&grade.length) {
            return Ok(Rc::clone(value));
        }
        let top = grade
            .length
            .checked_mul(2)
            .and_then(|n| usize::try_from(n).ok())
            .ok_or_else(work_limit)?;
        storage(top.saturating_add(1))?;
        model.spend(top + 1)?;
        let mut starts = Vec::with_capacity(top + 2);
        let mut dim = 0usize;
        for i in 0..=top {
            starts.push(dim);
            dim = dim
                .checked_add(self.rank(&grade.ring(i as i128))?)
                .ok_or_else(work_limit)?;
            storage(dim.saturating_mul(dim).saturating_mul(4))?;
        }
        starts.push(dim);
        let mut bases = Vec::with_capacity(dim);
        let mut weights = Vec::with_capacity(dim);
        let mut cells = starts.len();
        for i in 0..=top {
            let slice = self.slice(model, &grade.ring(i as i128))?;
            for j in slice.pivots() {
                let powers = &slice.generators[j];
                cells += 8 + 5 * powers.len();
                storage(cells)?;
                bases.push(powers.clone());
                weights.push(i);
            }
        }
        let trace = moments::integral(model, self, &Tensor::monomial(bases[dim - 1].clone()))?;
        cells += 4 * dim * dim;
        storage(cells)?;
        let mut algebra = Algebra {
            weights,
            bases,
            starts,
            trace,
            dual: RefCell::default(),
            product: vec![Vec::new(); dim * dim],
            splits: vec![Vec::new(); dim],
        };
        model.spend(dim * dim)?;
        for i in 0..dim {
            for j in i..dim {
                let weight = algebra.weights[i] + algebra.weights[j];
                if weight > top {
                    continue;
                }
                let left = Tensor::monomial(algebra.bases[i].clone());
                let right = Tensor::monomial(algebra.bases[j].clone());
                let raw = left.cup(model, &right)?;
                let reduced = self.project(model, &grade.ring(weight as i128), &raw)?;
                let entries = algebra.coefficient_indices(weight, &reduced)?;
                cells += entries.len() * if i == j { 1 } else { 2 } * (2 * LIMBS + 2);
                storage(cells)?;
                if i > 0 {
                    for &(k, coefficient) in &entries {
                        cells += 2 * LIMBS + 3;
                        storage(cells)?;
                        algebra.splits[k].push((i, j, coefficient));
                        if i != j {
                            cells += 2 * LIMBS + 3;
                            storage(cells)?;
                            algebra.splits[k].push((j, i, coefficient));
                        }
                    }
                }
                if i != j {
                    algebra.product[j * dim + i] = entries.clone();
                }
                algebra.product[i * dim + j] = entries;
            }
        }
        let algebra = Rc::new(algebra);
        model.allocate(cells)?;
        self.hochschild
            .borrow_mut()
            .insert(grade.length, Rc::clone(&algebra));
        Ok(algebra)
    }
}
