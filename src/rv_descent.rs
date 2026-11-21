use crate::binary_field::{BinaryField, Embedding};
use crate::binary_linear::{add_multiple, check_storage, Basis, Entry, Vector};
use crate::error::{Error, Result};
use crate::model::{Kind, Model, TypeId, Value};
use crate::rv_hom::{self, Grade};

#[derive(Clone, Copy)]
pub(crate) enum Twist {
    Descend,
    Extend,
}

impl Twist {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Descend => "rvTwistedDescend",
            Self::Extend => "rvTwistedExtend",
        }
    }
}

fn grade(model: &Model, ty: TypeId, name: &str) -> Result<Grade> {
    match model.kind(ty) {
        Kind::RVHom(grade) => Ok(grade.clone()),
        _ => Err(Error::plain(format!("{name} requires RVHom maps"))),
    }
}

fn endomorphism(grade: &Grade, name: &str) -> Result<()> {
    if grade.source != grade.target {
        return Err(Error::plain(format!(
            "{name} requires RVHom endomorphism data"
        )));
    }
    Ok(())
}

pub(crate) fn rational_point(
    model: &mut Model,
    grade: &Grade,
    modulus: usize,
    root: usize,
    name: &str,
) -> Result<(Embedding, Vec<usize>)> {
    let small = BinaryField::new(model, modulus)?;
    let embedding = Embedding::new(model, small, grade.field, root)?;
    model.spend(grade.point.len())?;
    let point = grade
        .point
        .iter()
        .map(|&value| {
            embedding.preimage(model, value)?.ok_or_else(|| {
                Error::plain(format!(
                    "{name} requires an affine point in the specified subfield"
                ))
            })
        })
        .collect::<Result<Vec<_>>>()?;
    Ok((embedding, point))
}

pub(crate) fn inverse_type(model: &mut Model, ty: TypeId) -> Result<TypeId> {
    let grade = grade(model, ty, "rvInverse")?;
    model.rv_hom_type(Grade {
        source: grade.target,
        target: grade.source,
        ..grade
    })
}

pub(crate) fn hilbert_embedding(
    model: &mut Model,
    ty: TypeId,
    modulus: usize,
    root: usize,
) -> Result<Embedding> {
    Ok(endomorphism_point(model, ty, modulus, root, "rvHilbert90")?.1)
}

pub(crate) fn endomorphism_point(
    model: &mut Model,
    ty: TypeId,
    modulus: usize,
    root: usize,
    name: &str,
) -> Result<(Grade, Embedding, Vec<usize>)> {
    let grade = grade(model, ty, name)?;
    endomorphism(&grade, name)?;
    let (embedding, point) = rational_point(model, &grade, modulus, root, name)?;
    Ok((grade, embedding, point))
}

pub(crate) fn twist_plan(
    model: &mut Model,
    operation: Twist,
    types: [TypeId; 3],
    modulus: usize,
    root: usize,
) -> Result<(TypeId, Embedding)> {
    let name = operation.name();
    let [left, input, right] = types;
    let left = grade(model, left, name)?;
    let input = grade(model, input, name)?;
    let right = grade(model, right, name)?;
    endomorphism(&left, name)?;
    endomorphism(&right, name)?;
    if left.selected != right.selected || left.point != right.point || left.field != right.field {
        return Err(Error::plain(format!(
            "{name} requires both data at the same affine Margolis field point"
        )));
    }
    let (embedding, point) = rational_point(model, &left, modulus, root, name)?;
    model.allocate(1 + point.len())?;
    let large = Grade {
        source: right.source,
        target: left.target,
        ..left
    };
    let small = Grade {
        field: embedding.small(),
        point: point.into(),
        ..large.clone()
    };
    let (expected, output) = match operation {
        Twist::Descend => (&large, small),
        Twist::Extend => (&small, large),
    };
    if input != *expected {
        return Err(Error::plain(format!(
            "{name} requires a map with the matching endpoints and affine field point"
        )));
    }
    Ok((model.rv_hom_type(output)?, embedding))
}

pub(crate) struct Matrix {
    rows: usize,
    columns: Vec<Vector>,
    cells: usize,
}

impl Matrix {
    pub(crate) fn new(rows: usize, columns: usize) -> Result<Self> {
        check_storage(columns)?;
        Ok(Self {
            rows,
            columns: Vec::with_capacity(columns),
            cells: columns,
        })
    }

    pub(crate) fn read(
        model: &mut Model,
        rows: usize,
        columns: usize,
        support: &[Entry],
        mut coefficient: impl FnMut(&mut Model, usize) -> Result<usize>,
    ) -> Result<Self> {
        model.spend(columns + support.len())?;
        let mut matrix = Self::new(rows, columns)?;
        matrix.cells += 2 * support.len();
        check_storage(matrix.cells)?;
        matrix.columns.resize_with(columns, Vec::new);
        for &(index, value) in support {
            matrix.columns[index / rows].push((index % rows, coefficient(model, value)?));
        }
        Ok(matrix)
    }

    pub(crate) fn push(&mut self, column: Vector) -> Result<()> {
        self.cells += 2 * column.len();
        check_storage(self.cells)?;
        self.columns.push(column);
        Ok(())
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        field: BinaryField,
        vector: &[Entry],
    ) -> Result<Vector> {
        model.spend(vector.len())?;
        let mut image = Vec::new();
        for &(column, coefficient) in vector {
            image = add_multiple(model, field, &image, &self.columns[column], coefficient)?;
        }
        Ok(image)
    }

    pub(crate) fn columns(&self) -> &[Vector] {
        &self.columns
    }

    pub(crate) fn semilinear_step(
        &self,
        model: &mut Model,
        field: BinaryField,
        step: usize,
        vector: &[Entry],
    ) -> Result<Vector> {
        model.spend(1 + vector.len())?;
        let mut conjugate = Vec::with_capacity(vector.len());
        for &(row, coefficient) in vector {
            conjugate.push((row, field.frobenius(model, coefficient, step)?));
        }
        self.apply(model, field, &conjugate)
    }

    pub(crate) fn value(self, model: &mut Model) -> Result<Value> {
        model.spend(self.cells)?;
        let mut support = Vec::new();
        for (column, entries) in self.columns.into_iter().enumerate() {
            for (row, coefficient) in entries {
                let coordinate = column
                    .checked_mul(self.rows)
                    .and_then(|offset| offset.checked_add(row))
                    .ok_or_else(|| Error::plain("Margolis Hom coordinate is out of range"))?;
                support.push((coordinate, coefficient));
            }
        }
        rv_hom::persist(model, support)
    }
}

struct Frame {
    matrix: Matrix,
    basis: Basis,
}

impl Frame {
    fn new(size: usize) -> Result<Self> {
        Ok(Self {
            matrix: Matrix::new(size, size)?,
            basis: Basis::default(),
        })
    }

    fn try_column(
        &mut self,
        model: &mut Model,
        field: BinaryField,
        column: Vector,
    ) -> Result<bool> {
        model.spend(column.len())?;
        let coordinate = self.basis.rank();
        let dependent =
            self.basis
                .insert_with_witness(model, field, column.clone(), vec![(coordinate, 1)])?;
        if dependent.is_none() {
            self.matrix.push(column)?;
            check_storage(self.matrix.cells + self.basis.cells())?;
            return Ok(true);
        }
        Ok(false)
    }

    fn solve(&self, model: &mut Model, field: BinaryField, vector: Vector) -> Result<Vector> {
        let mut coordinates = Vec::new();
        let remainder = self
            .basis
            .reduce(model, field, vector, Some(&mut coordinates))?;
        debug_assert!(remainder.is_empty());
        Ok(coordinates)
    }
}

pub(crate) fn inverse(model: &mut Model, ty: TypeId, support: &[Entry]) -> Result<Value> {
    let field = grade(model, ty, "rvInverse")?.field;
    let (source, target) = model.rv_hom(ty)?.dimensions();
    if source != target {
        return Err(Error::plain("rvInverse requires equal Margolis dimensions"));
    }
    let matrix = Matrix::read(model, target, source, support, |_, coefficient| {
        Ok(coefficient)
    })?;
    let mut frame = Frame::new(source)?;
    for column in matrix.columns {
        if !frame.try_column(model, field, column)? {
            return Err(Error::plain(
                "rvInverse requires an invertible Margolis map",
            ));
        }
    }
    let mut inverse = Matrix::new(source, target)?;
    for column in 0..target {
        inverse.push(frame.solve(model, field, vec![(column, 1)])?)?;
    }
    inverse.value(model)
}

fn fixed_frame(
    model: &mut Model,
    embedding: &Embedding,
    size: usize,
    support: &[Entry],
) -> Result<Frame> {
    let field = embedding.large();
    let step = embedding.small().degree();
    let extension = field.degree() / step;
    let matrix = Matrix::read(model, size, size, support, |_, coefficient| Ok(coefficient))?;
    let mut frame = Frame::new(size)?;
    if size == 0 {
        return Ok(frame);
    }

    // The polynomial generator of L also generates L over K. Its first
    // [L:K] powers have the indicated binary encodings in L.
    let mut scalars = Vec::with_capacity(extension);
    for power in 0..extension {
        let mut alpha = 1 << power;
        let mut conjugates = Vec::with_capacity(extension);
        for _ in 0..extension {
            conjugates.push(alpha);
            alpha = field.frobenius(model, alpha, step)?;
        }
        scalars.push(conjugates);
    }
    for column in 0..size {
        let unit = vec![(column, 1)];
        let mut orbit = Vec::with_capacity(extension);
        let mut cells = extension;
        let mut image = unit.clone();
        for _ in 0..extension {
            cells += 2 * image.len();
            check_storage(cells)?;
            let next = matrix.semilinear_step(model, field, step, &image)?;
            orbit.push(image);
            image = next;
        }
        // (A sigma)^[L:K] is linear and equals the cocycle norm.
        if image != unit {
            return Err(Error::plain(
                "Hilbert 90 requires the cocycle norm to equal the identity",
            ));
        }
        for conjugates in &scalars {
            if frame.basis.rank() == size {
                break;
            }
            // The trace sums of alpha^k e_j span the fixed frame because the
            // conjugate matrix of a K-basis is invertible over L.
            let mut fixed = Vec::new();
            for (image, &scalar) in orbit.iter().zip(conjugates) {
                fixed = add_multiple(model, field, &fixed, image, scalar)?;
            }
            frame.try_column(model, field, fixed)?;
        }
    }
    if frame.basis.rank() != size {
        return Err(Error::plain(
            "Hilbert 90 could not construct a full fixed frame",
        ));
    }
    Ok(frame)
}

pub(crate) fn hilbert(
    model: &mut Model,
    embedding: &Embedding,
    ty: TypeId,
    support: &[Entry],
) -> Result<Value> {
    let size = model.rv_hom(ty)?.dimensions().0;
    fixed_frame(model, embedding, size, support)?
        .matrix
        .value(model)
}

pub(crate) fn twisted(
    model: &mut Model,
    operation: Twist,
    embedding: &Embedding,
    types: [TypeId; 3],
    supports: [&[Entry]; 3],
) -> Result<Value> {
    let field = embedding.large();
    let source = model.rv_hom(types[2])?.dimensions().0;
    let target = model.rv_hom(types[0])?.dimensions().0;
    let [left, input, right] = supports;
    let left = fixed_frame(model, embedding, target, left)?;
    let right = fixed_frame(model, embedding, source, right)?;
    let input = Matrix::read(
        model,
        target,
        source,
        input,
        |model, coefficient| match operation {
            Twist::Descend => Ok(coefficient),
            Twist::Extend => embedding.image(model, coefficient),
        },
    )?;
    let mut output = Matrix::new(target, source)?;
    for column in 0..source {
        let image = match operation {
            Twist::Descend => {
                let image = input.apply(model, field, &right.matrix.columns[column])?;
                let mut coordinates = left.solve(model, field, image)?;
                // In fixed frames the intertwining equation says precisely
                // that every coefficient is fixed by relative Frobenius.
                for (_, coefficient) in &mut coordinates {
                    *coefficient = embedding.preimage(model, *coefficient)?.ok_or_else(|| {
                        Error::plain(
                            "rvTwistedDescend requires a map that intertwines the two descent data",
                        )
                    })?;
                }
                coordinates
            }
            Twist::Extend => {
                let coordinates = right.solve(model, field, vec![(column, 1)])?;
                let image = input.apply(model, field, &coordinates)?;
                left.matrix.apply(model, field, &image)?
            }
        };
        output.push(image)?;
    }
    output.value(model)
}
