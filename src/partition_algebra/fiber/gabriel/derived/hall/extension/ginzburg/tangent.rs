use crate::error::{Error, Result};
use crate::model::Model;

use super::super::{require_hall_length, Matrix, Representation, Setup, Span};
use super::{Field, Record, Spec, TangentQuery, Word};

impl Spec {
    pub(crate) fn potential_trace(&self, model: &mut Model, module: &Record) -> Result<String> {
        let potential = self.potential.as_ref().ok_or_else(super::no_potential)?;
        require_hall_length(module.module().length())?;
        let representation = Representation::new(&self.hall, module.module());
        let field = &self.hall.field;
        let mut answer = 0;
        for (word, &coefficient) in potential {
            let path = self.evaluate(model, &representation, word)?;
            answer = field.add(answer, field.mul(coefficient, matrix_trace(field, &path)));
        }
        Ok(format!(
            "rvPartitionDGPotentialTrace[{}] = {answer}",
            self.name
        ))
    }

    fn evaluate(
        &self,
        model: &mut Model,
        representation: &Representation,
        word: &Word,
    ) -> Result<Matrix> {
        let mut product = Matrix::identity(representation.dims[word.start]);
        for &letter in &word.letters {
            let arrow = self.symbols[usize::from(letter)];
            product = representation.arrows[arrow.source]
                .as_ref()
                .unwrap()
                .compose(model, &self.hall.field, &product)?;
        }
        Ok(product)
    }

    fn relation_derivative(
        &self,
        model: &mut Model,
        setup: &Setup,
        coordinate: usize,
    ) -> Result<Vec<u32>> {
        let middle = setup.middle(&unit(setup.variables, coordinate));
        let mut values = Vec::with_capacity(setup.variables);
        for (a, relation) in self.symbols[..self.arrows].iter().zip(&self.gradients) {
            let relation_value = self.evaluate(model, &middle, relation)?;
            // C² is dual to C¹ via Tr(B_a X_a), so store reversed-arrow entries transposed.
            for row in 0..setup.m.dims[a.target] {
                for col in 0..setup.m.dims[a.source] {
                    values.push(relation_value.get(col, setup.m.dims[a.target] + row));
                }
            }
        }
        Ok(values)
    }

    fn dual_gauge(&self, model: &mut Model, setup: &Setup, coordinate: usize) -> Result<Vec<u32>> {
        let field = &self.hall.field;
        let mut moment: Vec<_> = setup.m.dims.iter().map(|&d| Matrix::zero(d, d)).collect();
        for arrow in &self.symbols[..self.arrows] {
            let (s, t) = (arrow.source, arrow.target);
            let range = setup.arrows[s].as_ref().unwrap();
            if !range.contains(&coordinate) {
                continue;
            }
            let (ds, dt) = (setup.m.dims[s], setup.m.dims[t]);
            let offset = coordinate - range.start;
            let dual = Matrix::new(ds, dt, unit(ds * dt, (offset % ds) * dt + offset / ds));
            let a = setup.m.arrows[s].as_ref().unwrap();
            let target = a.compose(model, field, &dual)?;
            let source = dual.compose(model, field, a)?;
            moment[t].add_scaled(field, &target, 1);
            moment[s].add_scaled(field, &source, field.sub(0, 1));
        }
        Ok(moment
            .iter()
            .flat_map(|m| (0..m.cols).flat_map(move |j| (0..m.rows).map(move |i| m.get(i, j))))
            .collect())
    }

    fn potential_hessian(
        &self,
        model: &mut Model,
        setup: &Setup,
        directions: &[Direction],
    ) -> Result<Matrix> {
        let field = &self.hall.field;
        let mut hessian = vec![0; directions.len() * directions.len()];
        for (word, &coefficient) in self.potential.as_ref().unwrap() {
            for (i, x) in directions.iter().enumerate() {
                for (j, y) in directions.iter().enumerate().skip(i) {
                    let entry = &mut hessian[i * directions.len() + j];
                    for (k, &letter_x) in word.letters.iter().enumerate() {
                        if letter_x != x.arrow {
                            continue;
                        }
                        for (l, &letter_y) in word.letters.iter().enumerate() {
                            if l != k && letter_y == y.arrow {
                                let trace = self.mixed_trace(
                                    model,
                                    &setup.m,
                                    word,
                                    [(k, &x.matrix), (l, &y.matrix)],
                                )?;
                                *entry = field.add(*entry, field.mul(coefficient, trace));
                            }
                        }
                    }
                    hessian[j * directions.len() + i] = hessian[i * directions.len() + j];
                }
            }
        }
        Ok(Matrix::new(directions.len(), directions.len(), hessian))
    }

    fn mixed_trace(
        &self,
        model: &mut Model,
        representation: &Representation,
        word: &Word,
        variations: [(usize, &Matrix); 2],
    ) -> Result<u32> {
        let field = &self.hall.field;
        let mut product = Matrix::identity(representation.dims[word.start]);
        for (position, &letter) in word.letters.iter().enumerate() {
            let source = self.symbols[usize::from(letter)].source;
            let matrix = if let Some((_, matrix)) = variations.iter().find(|(i, _)| *i == position)
            {
                matrix
            } else {
                representation.arrows[source].as_ref().unwrap()
            };
            product = matrix.compose(model, field, &product)?;
        }
        Ok(matrix_trace(field, &product))
    }

    pub(crate) fn tangent(
        &self,
        model: &mut Model,
        query: TangentQuery,
        module: &Record,
    ) -> Result<String> {
        if matches!(query, TangentQuery::Cohomology(i) if i > 3) {
            return Err(Error::plain(
                "partition Ginzburg cohomology degree must be in 0..=3",
            ));
        }
        require_hall_length(module.module().length())?;
        let field = &self.hall.field;
        let shape = module.module();
        let setup = Setup::new(&self.hall, shape, shape);
        let (g, v) = (setup.gauges, setup.variables);
        model.spend(1 + 3 * (g + v).pow(2))?;
        let d0 = from_columns(v, g, |j| {
            let vector = setup.boundary(model, &self.hall, &unit(g, j))?;
            Ok(vector.into_iter().map(|x| field.sub(0, x)).collect())
        })?;
        let d1 = from_columns(v, v, |j| self.relation_derivative(model, &setup, j))?;
        let d2 = from_columns(g, v, |j| self.dual_gauge(model, &setup, j))?;
        if !is_transpose(&d1, &d1) || !is_transpose(&d0, &d2) {
            return Err(Error::plain(
                "partition Ginzburg Hessian or gauge adjoint certificate failed",
            ));
        }
        if !d1.compose(model, field, &d0)?.is_zero() || !d2.compose(model, field, &d1)?.is_zero() {
            return Err(Error::plain(
                "partition Ginzburg tangent differentials did not compose to zero",
            ));
        }
        if self.potential.is_some() {
            let mut directions = Vec::with_capacity(v);
            for (i, arrow) in self.symbols[..self.arrows].iter().enumerate() {
                let (ds, dt) = (setup.m.dims[arrow.source], setup.m.dims[arrow.target]);
                for j in 0..ds * dt {
                    directions.push(Direction {
                        arrow: i as u8,
                        matrix: Matrix::new(dt, ds, unit(ds * dt, j)),
                    });
                }
            }
            if d1 != self.potential_hessian(model, &setup, &directions)? {
                return Err(Error::plain(
                    "partition Ginzburg derivative differs from the potential trace Hessian",
                ));
            }
        }

        let dimensions = [g, v, v, g];
        let differentials = [d0, d1, d2];
        if matches!(query, TangentQuery::Differentials) {
            let matrices = differentials
                .each_ref()
                .map(|d| (0..d.rows).map(|i| row(d, i)).collect::<Vec<_>>());
            return Ok(format!(
                "rvPartitionDGComplex[{}] = complex={dimensions:?}, differentials={matrices:?}",
                self.name
            ));
        }
        let ranks = differentials
            .iter()
            .map(|d| d.rank(model, field))
            .collect::<Result<Vec<_>>>()?;
        let basis = (0usize..4)
            .map(|i| {
                let incoming = i.checked_sub(1).and_then(|j| differentials.get(j));
                cohomology(model, field, dimensions[i], incoming, differentials.get(i))
            })
            .collect::<Result<Vec<_>>>()?;
        let ext: Vec<_> = basis.iter().map(Vec::len).collect();
        if ext[1] != self.hall.projective_ext_dimension(shape, shape) {
            return Err(Error::plain(
                "partition Ginzburg Ext¹ disagrees with the partition module extension dimension",
            ));
        }
        let mut grams = Vec::with_capacity(4);
        let mut pairing = Vec::with_capacity(4);
        for i in 0..4 {
            let gram = gram(field, &basis[i], &basis[3 - i]);
            let rank = gram.rank(model, field)?;
            if rank != ext[i] || rank != ext[3 - i] {
                return Err(Error::plain(
                    "partition Ginzburg cohomology trace pairing is degenerate",
                ));
            }
            pairing.push(rank);
            grams.push(gram);
        }
        if let TangentQuery::Cohomology(i) = query {
            let i = i as usize;
            let coordinates: Vec<Vec<_>> = (0..grams[i].rows).map(|j| row(&grams[i], j)).collect();
            return Ok(format!("rvPartitionDGCohomology[{},{}] = basis={:?}, dualBasis={:?}, tracePairing={coordinates:?}", self.name, i, basis[i], basis[3-i]));
        }
        Ok(format!("rvPartitionDGTangent[{}] = complex={dimensions:?}, differentialRanks={ranks:?}, ext={ext:?}, tracePairingRanks={pairing:?}", self.name))
    }
}

struct Direction {
    arrow: u8,
    matrix: Matrix,
}

fn unit(width: usize, i: usize) -> Vec<u32> {
    let mut vector = vec![0; width];
    vector[i] = 1;
    vector
}

fn from_columns(
    rows: usize,
    cols: usize,
    mut f: impl FnMut(usize) -> Result<Vec<u32>>,
) -> Result<Matrix> {
    let mut entries = vec![0; rows * cols];
    for j in 0..cols {
        let column = f(j)?;
        debug_assert_eq!(column.len(), rows);
        for (i, value) in column.into_iter().enumerate() {
            entries[i * cols + j] = value;
        }
    }
    Ok(Matrix::new(rows, cols, entries))
}

fn row(m: &Matrix, i: usize) -> Vec<u32> {
    (0..m.cols).map(|j| m.get(i, j)).collect()
}

fn is_transpose(a: &Matrix, b: &Matrix) -> bool {
    a.rows == b.cols
        && a.cols == b.rows
        && (0..a.rows).all(|i| (0..a.cols).all(|j| a.get(i, j) == b.get(j, i)))
}

fn matrix_trace(field: &Field, m: &Matrix) -> u32 {
    (0..m.rows).fold(0, |sum, i| field.add(sum, m.get(i, i)))
}

fn cohomology(
    model: &mut Model,
    field: &Field,
    width: usize,
    incoming: Option<&Matrix>,
    outgoing: Option<&Matrix>,
) -> Result<Vec<Vec<u32>>> {
    let mut boundaries = Span::new(width);
    if let Some(d) = incoming {
        for j in 0..d.cols {
            boundaries.insert(model, field, (0..width).map(|i| d.get(i, j)).collect())?;
        }
    }
    let mut equations = Span::new(width);
    if let Some(d) = outgoing {
        for i in 0..d.rows {
            equations.insert(model, field, row(d, i))?;
        }
    }
    let mut quotient = Span::new(width);
    for cycle in equations.kernel(model, field)? {
        let representative = boundaries.reduce(model, field, cycle)?;
        quotient.insert(model, field, representative)?;
    }
    Ok(quotient.rows.into_iter().flatten().collect())
}

fn gram(field: &Field, a: &[Vec<u32>], b: &[Vec<u32>]) -> Matrix {
    Matrix::new(
        a.len(),
        b.len(),
        a.iter()
            .flat_map(|a| {
                b.iter().map(move |b| {
                    a.iter()
                        .zip(b)
                        .fold(0, |sum, (&x, &y)| field.add(sum, field.mul(x, y)))
                })
            })
            .collect(),
    )
}
