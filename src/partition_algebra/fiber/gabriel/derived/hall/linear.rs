use super::Field;
use crate::error::Result;
use crate::model::Model;

pub(super) type Row = Vec<u32>;

pub(super) struct Span {
    pub(super) rows: Vec<Option<Row>>,
    pub(super) rank: usize,
}

impl Span {
    pub(super) fn new(width: usize) -> Self {
        Self {
            rows: vec![None; width],
            rank: 0,
        }
    }

    pub(super) fn from_rows(width: usize, rows: impl Iterator<Item = (usize, Row)>) -> Self {
        let mut span = Self::new(width);
        for (pivot, row) in rows {
            span.rows[pivot] = Some(row);
            span.rank += 1;
        }
        span
    }

    pub(super) fn reduce(&self, model: &mut Model, field: &Field, mut row: Row) -> Result<Row> {
        for (pivot, old) in self.rows.iter().enumerate() {
            if let Some(old) = old {
                let scalar = row[pivot];
                if scalar == 0 {
                    continue;
                }
                model.spend(row.len() - pivot)?;
                for (x, &y) in row[pivot..].iter_mut().zip(&old[pivot..]) {
                    *x = field.sub(*x, field.mul(scalar, y));
                }
            }
        }
        Ok(row)
    }

    pub(super) fn insert(&mut self, model: &mut Model, field: &Field, row: Row) -> Result<()> {
        let mut row = self.reduce(model, field, row)?;
        if let Some(pivot) = row.iter().position(|&x| x != 0) {
            let scale = field.inverse(model, row[pivot])?;
            model.spend(row.len() - pivot)?;
            for term in &mut row[pivot..] {
                *term = field.mul(*term, scale);
            }
            self.rows[pivot] = Some(row);
            self.rank += 1;
        }
        Ok(())
    }
}

impl Span {
    pub(super) fn kernel(&self, model: &mut Model, field: &Field) -> Result<Vec<Row>> {
        let mut kernel = Vec::new();
        for (free, test) in self.rows.iter().enumerate() {
            if test.is_some() {
                continue;
            }
            let mut vector = vec![0; self.rows.len()];
            vector[free] = 1;
            for (pivot, row) in self.rows.iter().enumerate().rev() {
                if let Some(row) = row {
                    model.spend(row.len() - pivot)?;
                    let sum = row[pivot + 1..]
                        .iter()
                        .zip(&vector[pivot + 1..])
                        .fold(0, |sum, (&a, &b)| field.add(sum, field.mul(a, b)));
                    vector[pivot] = field.sub(0, sum);
                }
            }
            kernel.push(vector);
        }
        Ok(kernel)
    }

    pub(super) fn vanishes_on(
        &self,
        model: &mut Model,
        field: &Field,
        vector: &[u32],
    ) -> Result<bool> {
        for row in self.rows.iter().flatten() {
            model.spend(row.len())?;
            if row
                .iter()
                .zip(vector)
                .fold(0, |sum, (&a, &b)| field.add(sum, field.mul(a, b)))
                != 0
            {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub(super) fn coordinates(
        &self,
        model: &mut Model,
        field: &Field,
        mut value: Row,
    ) -> Result<Option<Row>> {
        let mut coefficients = Vec::new();
        for (pivot, basis) in self.rows.iter().enumerate() {
            if let Some(basis) = basis {
                let scalar = value[pivot];
                coefficients.push(scalar);
                if scalar == 0 {
                    continue;
                }
                model.spend(value.len() - pivot)?;
                for (a, &b) in value[pivot..].iter_mut().zip(&basis[pivot..]) {
                    *a = field.sub(*a, field.mul(scalar, b));
                }
            }
        }
        Ok(value.iter().all(|&x| x == 0).then_some(coefficients))
    }
}
