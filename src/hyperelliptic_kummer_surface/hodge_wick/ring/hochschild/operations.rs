use super::super::resolution::{deformation, Vector};
use super::complex::differential_terms;
use super::{storage, Algebra, Bilinear, Cochain, Form, Grade, Key, Model, Result, Tensor};

impl Form {
    pub(crate) fn hochschild_differential(
        &self,
        model: &mut Model,
        grade: &Grade,
        value: &Cochain,
    ) -> Result<Cochain> {
        let mut result = Cochain::default();
        if value.empty() || grade.arity <= 0 {
            return Ok(result);
        }
        let algebra = self.hochschild_algebra(model, grade)?;
        for (key, &value) in &value.entries {
            differential_terms(model, &algebra, key, |model, key, coefficient| {
                let coefficient = coefficient.multiply(model, value)?;
                result.insert(model, key, coefficient)
            })?;
        }
        Ok(result)
    }

    pub(crate) fn hochschild_bv(
        &self,
        model: &mut Model,
        grade: &Grade,
        value: &Cochain,
    ) -> Result<Cochain> {
        let mut result = Cochain::default();
        if value.empty() || grade.arity <= 0 {
            return Ok(result);
        }
        let algebra = self.hochschild_algebra(model, grade)?;
        let dual = algebra.pairing_inverse(model)?;
        for (key, &coefficient) in &value.entries {
            if key.output != algebra.top_index() {
                continue;
            }
            let trace = coefficient.multiply(model, algebra.trace_factor())?;
            let arity = key.word.len();
            for (i, &test) in key.word.iter().enumerate() {
                model.spend(arity)?;
                let sign = arity.is_multiple_of(2) && !i.is_multiple_of(2);
                let signed = if sign { trace.negate() } else { trace };
                // B rotates (test,a_1,...). A word with the marked test
                // at position i appears with the same parity after rotating
                // backwards because arity*(arity-1) is even.
                for output in algebra.weight_range(algebra.top() - algebra.weights[test]) {
                    let entry = signed.multiply(model, dual[output][test])?;
                    if !entry.is_zero() {
                        let mut word = Vec::with_capacity(arity - 1);
                        word.extend_from_slice(&key.word[i + 1..]);
                        word.extend_from_slice(&key.word[..i]);
                        result.insert(model, Key { word, output }, entry)?;
                    }
                }
            }
        }
        Ok(result)
    }

    pub(crate) fn hochschild_bilinear(
        &self,
        model: &mut Model,
        grade: &Grade,
        op: Bilinear,
        left: &Cochain,
        right: &Cochain,
    ) -> Result<Cochain> {
        let mut result = Cochain::default();
        if left.empty() || right.empty() || !grade.possible() {
            return Ok(result);
        }
        match op {
            Bilinear::Cup => {
                let algebra = self.hochschild_algebra(model, grade)?;
                cup(model, &algebra, &mut result, left, right)?;
            }
            Bilinear::Insert => insertions(model, &mut result, left, right, false)?,
            Bilinear::Bracket => {
                let p = left.entries.keys().next().unwrap().word.len();
                let q = right.entries.keys().next().unwrap().word.len();
                insertions(model, &mut result, left, right, false)?;
                let negative = !(p.is_multiple_of(2) && q.is_multiple_of(2));
                insertions(model, &mut result, right, left, negative)?;
            }
        }
        Ok(result)
    }

    pub(crate) fn hochschild_from_tangent(
        &self,
        model: &mut Model,
        grade: &Grade,
        value: &Vector,
    ) -> Result<Cochain> {
        let mut result = Cochain::default();
        if !grade.possible() || value.is_empty() {
            return Ok(result);
        }
        let algebra = self.hochschild_algebra(model, grade)?;
        let tangent_grade = deformation::Grade {
            basis: grade.basis.clone(),
            length: grade.length,
            shift: grade.shift,
        };
        let mut requests = Vec::new();
        let mut words = Vec::new();
        let mut cells = 0;
        for a in 1..algebra.len() {
            for b in a..algebra.len() {
                model.spend(1)?;
                let weight = algebra.weights[a] + algebra.weights[b];
                let Some(slot) = tangent_grade.slot(weight) else {
                    continue;
                };
                let left = Tensor::monomial(algebra.bases[a].clone());
                let right = Tensor::monomial(algebra.bases[b].clone());
                let mut raw = left.cup(model, &right)?;
                for &(index, coefficient) in algebra.product(a, b) {
                    raw.insert(model, algebra.bases[index].clone(), coefficient.negate())?;
                }
                if !raw.terms.is_empty() {
                    cells += raw.cells + 12;
                    storage(cells)?;
                    requests.push((weight, raw));
                    words.push((a, b, slot.degree as usize / 2));
                }
            }
        }
        let corrections =
            self.deformation_product_corrections(model, &tangent_grade, value, &requests)?;
        for ((a, b, weight), correction) in words.into_iter().zip(corrections) {
            for (output, coefficient) in algebra.coefficient_indices(weight, &correction)? {
                result.insert(
                    model,
                    Key {
                        word: vec![a, b],
                        output,
                    },
                    coefficient,
                )?;
                if a != b {
                    result.insert(
                        model,
                        Key {
                            word: vec![b, a],
                            output,
                        },
                        coefficient,
                    )?;
                }
            }
        }
        self.hochschild_class(model, grade, &result)
    }
}

fn cup(
    model: &mut Model,
    algebra: &Algebra,
    result: &mut Cochain,
    left: &Cochain,
    right: &Cochain,
) -> Result<()> {
    for (a, &f) in &left.entries {
        for (b, &g) in &right.entries {
            model.spend(1)?;
            let product = algebra.product(a.output, b.output);
            if product.is_empty() {
                continue;
            }
            let scalar = f.multiply(model, g)?;
            for &(output, coefficient) in product {
                let mut word = a.word.clone();
                word.extend_from_slice(&b.word);
                let entry = scalar.multiply(model, coefficient)?;
                result.insert(model, Key { word, output }, entry)?;
            }
        }
    }
    Ok(())
}

pub(super) fn insertions(
    model: &mut Model,
    result: &mut Cochain,
    outer: &Cochain,
    inner: &Cochain,
    negative: bool,
) -> Result<()> {
    for (a, &f) in &outer.entries {
        for (b, &g) in &inner.entries {
            model.spend(1 + a.word.len())?;
            // A normalized cochain vanishes when a unit is substituted;
            // this also handles the insertion of a scalar zero-cochain.
            if b.output == 0 {
                continue;
            }
            for (i, &slot) in a.word.iter().enumerate() {
                if slot != b.output {
                    continue;
                }
                let mut word = Vec::with_capacity(a.word.len() + b.word.len() - 1);
                word.extend_from_slice(&a.word[..i]);
                word.extend_from_slice(&b.word);
                word.extend_from_slice(&a.word[i + 1..]);
                let value = f.multiply(model, g)?;
                let sign = b.word.len().is_multiple_of(2) && !i.is_multiple_of(2);
                let value = if sign ^ negative {
                    value.negate()
                } else {
                    value
                };
                result.insert(
                    model,
                    Key {
                        word,
                        output: a.output,
                    },
                    value,
                )?;
            }
        }
    }
    Ok(())
}
