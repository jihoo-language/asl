use super::{
    add_exponents, check_report, storage, BTreeMap, Error, Hom, Model, Morphism, Object, Rational,
    Result, Write, LIMBS,
};

#[derive(Clone, Default)]
struct Polynomial(BTreeMap<u128, Rational>);

impl Polynomial {
    fn term(exponent: u128, value: i128) -> Self {
        Self(BTreeMap::from([(exponent, Rational::integer(value))]))
    }

    fn empty(&self) -> bool {
        self.0.is_empty()
    }

    fn add(&mut self, model: &mut Model, exponent: u128, value: Rational) -> Result<()> {
        if value.is_zero() {
            return Ok(());
        }
        model.spend(1)?;
        let prior = self
            .0
            .get(&exponent)
            .copied()
            .unwrap_or_else(Rational::zero);
        let sum = prior.add(model, value)?;
        if sum.is_zero() {
            self.0.remove(&exponent);
        } else {
            storage((self.0.len() + usize::from(prior.is_zero())).saturating_mul(2 * LIMBS + 8))?;
            self.0.insert(exponent, sum);
        }
        Ok(())
    }

    fn add_product(&mut self, model: &mut Model, left: &Self, right: &Self) -> Result<()> {
        for (&a, &x) in &left.0 {
            for (&b, &y) in &right.0 {
                model.spend(1)?;
                let value = x.multiply(model, y)?;
                self.add(model, add_exponents(a, b)?, value)?;
            }
        }
        Ok(())
    }

    fn negate(&self) -> Self {
        Self(self.0.iter().map(|(&a, &v)| (a, v.negate())).collect())
    }

    fn scale(&mut self, model: &mut Model, factor: Rational) -> Result<()> {
        for value in self.0.values_mut() {
            *value = value.multiply(model, factor)?;
        }
        Ok(())
    }

    fn divide(&self, model: &mut Model, divisor: &Self) -> Result<(Self, Self)> {
        let (&degree, &leading) = divisor.0.last_key_value().ok_or_else(|| {
            Error::plain("division by zero in matrix factorization Smith calculation")
        })?;
        let inverse = leading.inverse()?;
        let mut quotient = Self::default();
        let mut remainder = self.clone();
        while let Some((&exponent, &value)) = remainder.0.last_key_value() {
            let Some(power) = exponent.checked_sub(degree) else {
                break;
            };
            let coefficient = value.multiply(model, inverse)?;
            quotient.add(model, power, coefficient)?;
            for (&e, &v) in &divisor.0 {
                let correction = v.multiply(model, coefficient)?.negate();
                remainder.add(model, add_exponents(e, power)?, correction)?;
            }
        }
        Ok((quotient, remainder))
    }

    fn monomial(&self, exponent: u128) -> bool {
        self.0.len() == 1 && self.0.get(&exponent) == Some(&Rational::integer(1))
    }
}

struct Matrix(Vec<Vec<Polynomial>>);

impl Matrix {
    fn new(dimension: usize) -> Result<Self> {
        storage(dimension.saturating_mul(dimension).saturating_mul(5))?;
        Ok(Self(vec![
            vec![Polynomial::default(); dimension];
            dimension
        ]))
    }

    fn cells(&self) -> usize {
        self.0
            .iter()
            .flat_map(|r| r.iter())
            .map(|p| 5 + p.0.len() * (2 * LIMBS + 8))
            .sum()
    }

    fn swap_columns(&mut self, a: usize, b: usize) {
        for row in &mut self.0 {
            row.swap(a, b);
        }
    }

    fn add_rows(
        &mut self,
        model: &mut Model,
        target: usize,
        source: usize,
        factor: &Polynomial,
    ) -> Result<()> {
        let (target, source) = if target < source {
            let (left, right) = self.0.split_at_mut(source);
            (&mut left[target], &right[0])
        } else {
            let (left, right) = self.0.split_at_mut(target);
            (&mut right[0], &left[source])
        };
        for (left, right) in target.iter_mut().zip(source) {
            model.spend(1)?;
            left.add_product(model, factor, right)?;
        }
        Ok(())
    }

    fn add_columns(
        &mut self,
        model: &mut Model,
        target: usize,
        source: usize,
        factor: &Polynomial,
    ) -> Result<()> {
        for row in &mut self.0 {
            let (target, source) = if target < source {
                let (left, right) = row.split_at_mut(source);
                (&mut left[target], &right[0])
            } else {
                let (left, right) = row.split_at_mut(target);
                (&mut right[0], &left[source])
            };
            model.spend(1)?;
            target.add_product(model, source, factor)?;
        }
        Ok(())
    }

    fn verify_product(&self, model: &mut Model, right: &Self, exponent: u128) -> Result<()> {
        for (i, row) in self.0.iter().enumerate() {
            for j in 0..self.0.len() {
                let mut product = Polynomial::default();
                for (k, left) in row.iter().enumerate() {
                    model.spend(1)?;
                    product.add_product(model, left, &right.0[k][j])?;
                }
                if (i == j && !product.monomial(exponent)) || (i != j && !product.empty()) {
                    return Err(Error::plain(
                        "matrix factorization cone does not square to the potential",
                    ));
                }
            }
        }
        Ok(())
    }

    fn text(&self, model: &mut Model) -> Result<String> {
        let mut text = "[".to_owned();
        let mut sep = "";
        for (i, row) in self.0.iter().enumerate() {
            for (j, poly) in row.iter().enumerate() {
                for (&exp, &value) in &poly.0 {
                    write!(text, "{sep}({i}, {j}, {exp}, {})", value.text(model)?).unwrap();
                    check_report(&text)?;
                    sep = ", ";
                }
            }
        }
        text.push(']');
        Ok(text)
    }
}

struct Cone {
    forward: Matrix,
    reverse: Matrix,
}

impl Cone {
    fn check_storage(&self) -> Result<()> {
        storage(self.forward.cells().saturating_add(self.reverse.cells()))
    }

    fn swap_rows(&mut self, a: usize, b: usize) {
        self.forward.0.swap(a, b);
        self.reverse.swap_columns(a, b);
    }

    fn swap_columns(&mut self, a: usize, b: usize) {
        self.forward.swap_columns(a, b);
        self.reverse.0.swap(a, b);
    }

    fn add_rows(
        &mut self,
        model: &mut Model,
        target: usize,
        source: usize,
        factor: &Polynomial,
    ) -> Result<()> {
        self.forward.add_rows(model, target, source, factor)?;
        self.reverse
            .add_columns(model, source, target, &factor.negate())?;
        self.check_storage()
    }

    fn add_columns(
        &mut self,
        model: &mut Model,
        target: usize,
        source: usize,
        factor: &Polynomial,
    ) -> Result<()> {
        self.forward.add_columns(model, target, source, factor)?;
        self.reverse
            .add_rows(model, source, target, &factor.negate())?;
        self.check_storage()
    }

    fn pivot(&mut self, model: &mut Model, index: usize) -> Result<()> {
        let size = self.forward.0.len();
        let mut candidate = None;
        for i in index..size {
            for j in index..size {
                model.spend(1)?;
                if let Some((&degree, _)) = self.forward.0[i][j].0.last_key_value() {
                    if candidate.is_none_or(|(prior, _, _)| degree < prior) {
                        candidate = Some((degree, i, j));
                    }
                }
            }
        }
        let (_, i, j) = candidate
            .ok_or_else(|| Error::plain("matrix factorization cone has a zero Smith pivot"))?;
        self.swap_rows(i, index);
        self.swap_columns(j, index);
        'euclidean: loop {
            for i in index + 1..size {
                if self.forward.0[i][index].empty() {
                    continue;
                }
                let (q, r) =
                    self.forward.0[i][index].divide(model, &self.forward.0[index][index])?;
                self.add_rows(model, i, index, &q.negate())?;
                if !r.empty() {
                    self.swap_rows(i, index);
                    continue 'euclidean;
                }
            }
            for j in index + 1..size {
                if self.forward.0[index][j].empty() {
                    continue;
                }
                let (q, r) =
                    self.forward.0[index][j].divide(model, &self.forward.0[index][index])?;
                self.add_columns(model, j, index, &q.negate())?;
                if !r.empty() {
                    self.swap_columns(j, index);
                    continue 'euclidean;
                }
            }
            for i in index + 1..size {
                for j in index + 1..size {
                    if self.forward.0[i][j].empty() {
                        continue;
                    }
                    if !self.forward.0[i][j]
                        .divide(model, &self.forward.0[index][index])?
                        .1
                        .empty()
                    {
                        self.add_rows(model, index, i, &Polynomial::term(0, 1))?;
                        continue 'euclidean;
                    }
                }
            }
            break;
        }
        let &leading = self.forward.0[index][index].0.last_key_value().unwrap().1;
        let inverse = leading.inverse()?;
        for poly in &mut self.forward.0[index] {
            poly.scale(model, inverse)?;
        }
        for row in &mut self.reverse.0 {
            row[index].scale(model, leading)?;
        }
        Ok(())
    }

    fn smith(&mut self, model: &mut Model, d: usize) -> Result<Vec<usize>> {
        let size = self.forward.0.len();
        for index in 0..size {
            self.pivot(model, index)?;
        }
        let mut exponents = Vec::with_capacity(size);
        for i in 0..size {
            for j in 0..size {
                if i != j && (!self.forward.0[i][j].empty() || !self.reverse.0[i][j].empty()) {
                    return Err(Error::plain(
                        "matrix factorization Smith result has an off-diagonal entry",
                    ));
                }
            }
            let exponent = self.forward.0[i][i]
                .0
                .first_key_value()
                .map(|(&power, _)| power)
                .filter(|&power| power <= d as u128)
                .ok_or_else(|| {
                    Error::plain("matrix factorization Smith divisor lies outside the potential")
                })?;
            if !self.forward.0[i][i].monomial(exponent)
                || !self.reverse.0[i][i].monomial(d as u128 - exponent)
            {
                return Err(Error::plain(
                    "matrix factorization Smith divisors do not multiply to the potential",
                ));
            }
            exponents.push(exponent as usize);
        }
        if exponents.windows(2).any(|pair| pair[0] > pair[1]) {
            return Err(Error::plain(
                "matrix factorization Smith divisors are not ordered by divisibility",
            ));
        }
        Ok(exponents)
    }
}

impl Hom<'_> {
    pub(crate) fn cone(
        &self,
        model: &mut Model,
        map: &Morphism,
        report: bool,
    ) -> Result<(Object, Option<String>)> {
        if self.odd {
            return Err(Error::plain(
                "matrix factorization cone requires an even morphism",
            ));
        }
        self.require_cycle(model, map)?;
        let d = self.source.d;
        let split = self.target.factors.len();
        let size = split + self.source.factors.len();
        let mut cone = Cone {
            forward: Matrix::new(size)?,
            reverse: Matrix::new(size)?,
        };
        for (i, &j) in self.target.factors.iter().enumerate() {
            cone.forward.0[i][i] = Polynomial::term(j as u128, 1);
            cone.reverse.0[i][i] = Polynomial::term((d - j) as u128, 1);
        }
        for (i, &j) in self.source.factors.iter().enumerate() {
            let i = i + split;
            cone.forward.0[i][i] = Polynomial::term((d - j) as u128, -1);
            cone.reverse.0[i][i] = Polynomial::term(j as u128, -1);
        }
        for (&key, &coefficient) in &map.terms {
            let matrix = if key.slot {
                &mut cone.forward
            } else {
                &mut cone.reverse
            };
            matrix.0[key.row][key.column + split].add(model, key.exponent, coefficient)?;
        }
        cone.check_storage()?;
        cone.forward
            .verify_product(model, &cone.reverse, d as u128)?;
        cone.reverse
            .verify_product(model, &cone.forward, d as u128)?;
        let original = if report {
            Some((cone.forward.text(model)?, cone.reverse.text(model)?))
        } else {
            None
        };
        let smith = cone.smith(model, d)?;
        let factors = smith
            .iter()
            .copied()
            .filter(|&i| i > 0 && i < d)
            .collect::<Vec<_>>();
        let description=original.map(|(forward,reverse)| format!(
            "even_to_odd={forward} odd_to_even={reverse} Smith={smith:?} factors={factors:?} contractible={}", smith.len()-factors.len()
        ));
        if let Some(text) = description.as_ref() {
            check_report(text)?;
        }
        let object = Object {
            mark: self.source.mark.clone(),
            d,
            factors,
        };
        Ok((object, description))
    }
}
