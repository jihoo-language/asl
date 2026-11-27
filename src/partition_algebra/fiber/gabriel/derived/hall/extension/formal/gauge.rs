use crate::error::{Error, Result};
use crate::model::Model;

use super::{homogeneous, Coefficients, Field, Frame, Jet, Matrix, Monomial, Series, Terms};

impl Frame {
    pub(crate) fn gauge(&self, model: &mut Model, jet: &Jet, terms: Terms) -> Result<Jet> {
        let terms = self.read_terms(terms, 1..=jet.order, self.group.setup.gauges)?;
        self.apply_gauge(model, jet, &terms)
    }

    pub(crate) fn slice(&self, model: &mut Model, jet: &Jet) -> Result<Jet> {
        let field = &self.group.hall.field;
        let mut jet = self.jet(jet.order, jet.coefficients.clone());
        for degree in 1..=jet.order {
            let mut corrections = Coefficients::new();
            let mut expected = Coefficients::new();
            for power in homogeneous(self.axes, degree) {
                let Some(vector) = jet.coefficients.get(&power) else {
                    continue;
                };
                let (h, reduced) = self.boundary.solve(model, field, vector)?;
                if h.iter().any(|&x| x != 0) {
                    corrections.insert(power, h);
                }
                if reduced.iter().any(|&x| x != 0) {
                    expected.insert(power, reduced);
                }
            }
            if !corrections.is_empty() {
                jet = self.apply_gauge(model, &jet, &corrections)?;
            }
            let verified: Coefficients = jet
                .coefficients
                .iter()
                .filter(|(p, _)| p.degree() == degree)
                .map(|(&p, v)| (p, v.clone()))
                .collect();
            if verified != expected {
                return Err(Error::plain(
                    "partition formal slice did not remove the degreewise coboundary pivots",
                ));
            }
        }
        Ok(jet)
    }

    fn apply_gauge(&self, model: &mut Model, jet: &Jet, terms: &Coefficients) -> Result<Jet> {
        let field = &self.group.hall.field;
        let mut gauges = Vec::new();
        for (vertex, &dimension) in self.group.setup.m.dims.iter().enumerate() {
            let mut h = Series::zero(dimension, dimension);
            for (&monomial, value) in terms {
                let entries = value[self.group.setup.splittings[vertex].clone()].to_vec();
                h.add_term(
                    model,
                    field,
                    monomial,
                    &Matrix::new(dimension, dimension, entries),
                    1,
                )?;
            }
            gauges.push(unipotent(model, field, h, jet.order)?);
        }
        let mut accumulated = Coefficients::new();
        let arrows = self.arrows(model, jet)?;
        for (source, arrow) in arrows.into_iter().enumerate() {
            let Some(arrow) = arrow else {
                continue;
            };
            let target = self.group.hall.successor[source].unwrap();
            let right = arrow.compose(model, field, &gauges[source].1, jet.order)?;
            let after = gauges[target].0.compose(model, field, &right, jet.order)?;
            let base = self.group.setup.m.arrows[source].as_ref().unwrap();
            let constant = after.terms.get(&Monomial::ONE);
            if !constant.map_or_else(|| base.is_zero(), |c| c == base) {
                return Err(Error::plain(
                    "partition formal gauge changed the constant representation",
                ));
            }
            let range = self.group.setup.arrows[source].as_ref().unwrap();
            for (power, matrix) in after.terms {
                if power == Monomial::ONE {
                    continue;
                }
                let block = accumulated
                    .entry(power)
                    .or_insert_with(|| vec![0; self.group.setup.variables]);
                for i in 0..matrix.rows {
                    for j in 0..matrix.cols {
                        block[range.start + i * matrix.cols + j] = matrix.get(i, j);
                    }
                }
            }
        }
        let transformed = self.jet(jet.order, accumulated);
        self.verify(model, &transformed)?;
        for axis in 0..self.axes as u128 {
            if self.first(model, jet, axis)? != self.first(model, &transformed, axis)? {
                return Err(Error::plain(
                    "partition formal gauge changed the first-order extension class",
                ));
            }
        }
        Ok(transformed)
    }
}

fn unipotent(model: &mut Model, field: &Field, h: Series, order: u8) -> Result<(Series, Series)> {
    let identity = Series::identity(h.rows);
    let mut g = identity.clone();
    g.add(model, field, &h, 1)?;
    let mut negative = Series::zero(h.rows, h.cols);
    negative.add(model, field, &h, field.sub(0, 1))?;
    let mut inverse = identity.clone();
    let mut term = identity.clone();
    for _ in 0..order {
        term = negative.compose(model, field, &term, order)?;
        if term.terms.is_empty() {
            break;
        }
        inverse.add(model, field, &term, 1)?;
    }
    if g.compose(model, field, &inverse, order)? != identity
        || inverse.compose(model, field, &g, order)? != identity
    {
        return Err(Error::plain(
            "partition formal vertex gauge failed invertibility",
        ));
    }
    Ok((g, inverse))
}
