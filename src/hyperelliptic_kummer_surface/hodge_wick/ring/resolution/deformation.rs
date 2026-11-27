use std::fmt::Write;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};

use super::super::super::rational::LIMBS;
use super::super::linear::{self, Echelon, Matrix};
use super::builder::{Axis, Catalog};
use super::{
    check_report, koszul, multiply_powers, storage, Entry, Form, Monomial, Rational, Resolution,
    RingGrade, Tensor, Vector,
};

pub(crate) enum Command {
    Dims(Label, i128, i128),
    Table(Label, i128),
    Obstruction(Box<ParsedExpr>),
}

pub(crate) enum Expression {
    Basis(bool, Label, i128, i128, usize),
    Literal(Label, i128, i128, Vec<(usize, ParsedExpr)>),
    Cast(bool, Box<ParsedExpr>),
    Scale(Entry, Box<ParsedExpr>),
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub(crate) struct Grade {
    pub(crate) basis: Rc<str>,
    pub(crate) length: i128,
    pub(crate) shift: i128,
}

impl Grade {
    fn resolution(&self) -> super::Grade {
        super::Grade {
            basis: self.basis.clone(),
            length: self.length,
            homological: 0,
            weight: 0,
        }
    }

    fn ring(&self, weight: i128) -> RingGrade {
        RingGrade {
            basis: self.basis.clone(),
            length: self.length,
            degree: 2 * weight,
        }
    }

    // Extreme signed shifts have no supported coordinates. Test the
    // bounds before multiplying a weight by two or allocating a frame.
    pub(crate) fn slot(&self, degree: usize) -> Option<RingGrade> {
        let weight = self.shift.checked_add(degree as i128)?;
        (weight >= 0 && weight <= self.length.saturating_mul(2) && weight <= i128::MAX / 2)
            .then(|| self.ring(weight))
    }
}

#[derive(Default)]
struct Frame(Vec<(usize, Monomial)>);

impl Frame {
    fn new(model: &mut Model, form: &Form, grade: &Grade, shifts: &[usize]) -> Result<Self> {
        let mut frame = Vec::new();
        let mut cells = 0;
        for (i, &degree) in shifts.iter().enumerate() {
            if let Some(slot) = grade.slot(degree) {
                let slice = form.slice(model, &slot)?;
                model.spend(slice.rank())?;
                for pivot in slice.pivots() {
                    let powers = &slice.generators[pivot];
                    cells += 10 + 5 * powers.len();
                    storage(cells)?;
                    frame.push((i, powers.clone()));
                }
            }
        }
        Ok(Self(frame))
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn fill(&self, i: usize, value: &Tensor, mut put: impl FnMut(usize, Rational)) -> Result<()> {
        for (powers, &coefficient) in &value.terms {
            let j = self
                .0
                .binary_search_by(|(row, monomial)| row.cmp(&i).then_with(|| monomial.cmp(powers)))
                .map_err(|_| {
                    Error::plain(
                        "Verbitsky deformation coordinate lies outside its graded Hom module",
                    )
                })?;
            put(j, coefficient);
        }
        Ok(())
    }

    fn coordinates(&self, value: &Vector) -> Result<Vec<Rational>> {
        let mut row = vec![Rational::zero(); self.len()];
        for (&i, polynomial) in &value.entries {
            self.fill(i, polynomial, |j, coefficient| row[j] = coefficient)?;
        }
        Ok(row)
    }

    fn vector(&self, model: &mut Model, row: &[Rational]) -> Result<Vector> {
        let mut value = Vector::default();
        for ((i, powers), &coefficient) in self.0.iter().zip(row) {
            value.insert(model, *i, powers.clone(), coefficient)?;
        }
        Ok(value)
    }
}

struct Hom<'a> {
    form: &'a Form,
    grade: &'a Grade,
    resolution: &'a Resolution,
    shifts: Vec<&'a [usize]>,
}

impl<'a> Hom<'a> {
    fn new(form: &'a Form, grade: &'a Grade, resolution: &'a Resolution) -> Self {
        Self {
            form,
            grade,
            resolution,
            shifts: resolution.shifts(),
        }
    }

    fn shifts(&self, p: usize) -> &[usize] {
        self.shifts.get(p).copied().unwrap_or_default()
    }

    fn project(&self, model: &mut Model, p: usize, j: usize, value: &Tensor) -> Result<Tensor> {
        let degree = self.resolution.shift(p, j)?;
        self.grade.slot(degree).map_or_else(
            || Ok(Tensor::default()),
            |slot| self.form.project(model, &slot, value),
        )
    }

    fn differential(
        &self,
        model: &mut Model,
        p: usize,
        from: &Frame,
        to: &Frame,
    ) -> Result<Matrix> {
        let mut matrix = linear::matrix(to.len(), from.len())?;
        for (i, (index, powers)) in from.0.iter().enumerate() {
            for (j, &shift) in self.shifts(p + 1).iter().enumerate() {
                if self.grade.slot(shift).is_none() {
                    continue;
                }
                if let Some(polynomial) = self.resolution.entry(p + 1, *index, j)? {
                    let product = monomial_times(model, powers, polynomial)?;
                    let reduced = self.project(model, p + 1, j, &product)?;
                    to.fill(j, &reduced, |row, coefficient| matrix[row][i] = coefficient)?;
                }
            }
        }
        Ok(matrix)
    }

    fn derivations(&self, model: &mut Model, targets: &Frame) -> Result<Matrix> {
        let weights = (0..self.form.variable_count())
            .map(|i| if i < self.form.gram.size() { 1 } else { 2 })
            .collect::<Vec<_>>();
        let sources = Frame::new(model, self.form, self.grade, &weights)?;
        let mut rows = linear::matrix(sources.len(), targets.len())?;
        for ((variable, powers), row) in sources.0.iter().zip(&mut rows) {
            for (j, &shift) in self.shifts(1).iter().enumerate() {
                if self.grade.slot(shift).is_none() {
                    continue;
                }
                if let Some(polynomial) = self.resolution.entry(1, 0, j)? {
                    let derivative = differentiate(model, polynomial, *variable)?;
                    let product = monomial_times(model, powers, &derivative)?;
                    let reduced = self.project(model, 1, j, &product)?;
                    targets.fill(j, &reduced, |index, coefficient| row[index] = coefficient)?;
                }
            }
        }
        Ok(rows)
    }
}

fn monomial_times(model: &mut Model, powers: &Monomial, polynomial: &Tensor) -> Result<Tensor> {
    let mut output = Tensor::default();
    for (monomial, &coefficient) in &polynomial.terms {
        model.spend(1)?;
        output.insert(model, multiply_powers(powers, monomial)?, coefficient)?;
    }
    Ok(output)
}

fn differentiate(model: &mut Model, polynomial: &Tensor, variable: usize) -> Result<Tensor> {
    let mut derivative = Tensor::default();
    for (powers, &coefficient) in &polynomial.terms {
        if let Some(i) = powers.iter().position(|p| usize::from(p.index) == variable) {
            let mut monomial = powers.clone();
            let exponent = monomial[i].exponent;
            if exponent == 1 {
                monomial.remove(i);
            } else {
                monomial[i].exponent -= 1;
            }
            let exponent = i128::try_from(exponent).map_err(|_| koszul::grading_range())?;
            let value = coefficient.multiply(model, Rational::integer(exponent))?;
            derivative.insert(model, monomial, value)?;
        }
    }
    Ok(derivative)
}

fn certify_cycle(model: &mut Model, equations: &Echelon, value: &[Rational]) -> Result<bool> {
    for row in &equations.rows {
        let mut sum = Rational::zero();
        for (&left, &right) in row.iter().zip(value) {
            let product = left.multiply(model, right)?;
            sum = sum.add(model, product)?;
        }
        if !sum.is_zero() {
            return Ok(false);
        }
    }
    Ok(true)
}

fn reduce(model: &mut Model, boundaries: &Echelon, value: &mut [Rational]) -> Result<()> {
    for (&pivot, row) in boundaries.pivots.iter().zip(&boundaries.rows) {
        linear::add_multiple(model, value, row, value[pivot].negate())?;
    }
    Ok(())
}

struct Quotient {
    constraints: Echelon,
    boundaries: Echelon,
    dimension: usize,
}

impl Quotient {
    fn new(
        model: &mut Model,
        constraints: Matrix,
        boundaries: Matrix,
        width: usize,
    ) -> Result<Self> {
        let constraints = linear::rref(model, constraints, width)?;
        let boundaries = linear::rref(model, boundaries, width)?;
        for row in &boundaries.rows {
            if !certify_cycle(model, &constraints, row)? {
                return Err(Error::plain(
                    "Verbitsky deformation complex failed the cycle check",
                ));
            }
        }
        let dimension = width - constraints.pivots.len() - boundaries.pivots.len();
        Ok(Self {
            constraints,
            boundaries,
            dimension,
        })
    }

    fn basis(&self, model: &mut Model, width: usize) -> Result<Matrix> {
        let mut cycles = self.constraints.kernel(width)?;
        for row in &mut cycles {
            reduce(model, &self.boundaries, row)?;
        }
        let basis = linear::rref(model, cycles, width)?.rows;
        if basis.len() != self.dimension {
            return Err(Error::plain(
                "Verbitsky deformation quotient has unexpected dimension",
            ));
        }
        Ok(basis)
    }

    fn canonical(&self, model: &mut Model, value: &mut [Rational]) -> Result<()> {
        reduce(model, &self.boundaries, value)
    }
}

struct Complex {
    generators: Frame,
    syzygies: Frame,
    differential: Matrix,
    abstract_space: Quotient,
    obstruction_space: Quotient,
}

impl Complex {
    fn new(model: &mut Model, hom: &Hom<'_>) -> Result<Self> {
        let generators = Frame::new(model, hom.form, hom.grade, hom.shifts(1))?;
        let syzygies = Frame::new(model, hom.form, hom.grade, hom.shifts(2))?;
        let identities = Frame::new(model, hom.form, hom.grade, hom.shifts(3))?;
        let differential = hom.differential(model, 1, &generators, &syzygies)?;
        let next = hom.differential(model, 2, &syzygies, &identities)?;
        let derivations = hom.derivations(model, &generators)?;
        let mut image = linear::matrix(generators.len(), syzygies.len())?;
        for (j, row) in differential.iter().enumerate() {
            for (i, &entry) in row.iter().enumerate() {
                image[i][j] = entry;
            }
        }
        let abstract_space =
            Quotient::new(model, differential.clone(), derivations, generators.len())?;
        let obstruction_space = Quotient::new(model, next, image, syzygies.len())?;
        Ok(Self {
            generators,
            syzygies,
            differential,
            abstract_space,
            obstruction_space,
        })
    }

    fn dims(&self) -> [usize; 3] {
        [
            self.generators.len() - self.abstract_space.constraints.pivots.len(),
            self.abstract_space.dimension,
            self.obstruction_space.dimension,
        ]
    }
}

impl Form {
    pub(crate) fn validate_deformation(&self, grade: &Grade) -> Result<()> {
        self.validate_free(&grade.resolution())
    }

    pub(crate) fn deformation_slot(
        &self,
        model: &mut Model,
        grade: &Grade,
        index: usize,
    ) -> Result<RingGrade> {
        let resolution = self.resolution(model, &grade.resolution())?;
        grade.slot(resolution.shift(1, index)?).ok_or_else(|| {
            Error::plain("Verbitsky tangent generator has no ring component at this shift")
        })
    }

    pub(crate) fn deformation_dims(&self, model: &mut Model, grade: &Grade) -> Result<[usize; 3]> {
        self.validate_deformation(grade)?;
        let top = grade.length.saturating_mul(2);
        let last_shift = top.saturating_add((self.variable_count() + self.field.degree()) as i128);
        if grade.shift < -last_shift || grade.shift > top {
            return Ok([0; 3]);
        }
        let resolution = self.resolution(model, &grade.resolution())?;
        Ok(Complex::new(model, &Hom::new(self, grade, &resolution))?.dims())
    }

    pub(crate) fn deformation_table(&self, model: &mut Model, grade: &Grade) -> Result<String> {
        let resolution = self.resolution(model, &grade.resolution())?;
        let shifts = resolution.shifts();
        let min = -(*shifts.iter().copied().flatten().max().unwrap_or(&0) as i128);
        let max = 2 * grade.length;
        let mut totals = [0; 3];
        let mut text = String::from("[(shift, tangent, T1, Ext1)");
        for shift in min..=max {
            model.spend(1)?;
            let grade = Grade {
                shift,
                ..grade.clone()
            };
            let dims = Complex::new(model, &Hom::new(self, &grade, &resolution))?.dims();
            if dims.iter().any(|&dim| dim != 0) {
                write!(text, ", ({shift}, {}, {}, {})", dims[0], dims[1], dims[2]).unwrap();
                check_report(&text)?;
                for (total, dim) in totals.iter_mut().zip(dims) {
                    *total += dim;
                }
            }
        }
        write!(
            text,
            "] totals = ({}, {}, {})",
            totals[0], totals[1], totals[2]
        )
        .unwrap();
        check_report(&text)?;
        Ok(text)
    }

    pub(crate) fn deformation_basis(
        &self,
        model: &mut Model,
        grade: &Grade,
        abstract_class: bool,
        index: usize,
    ) -> Result<Vector> {
        let resolution = self.resolution(model, &grade.resolution())?;
        let hom = Hom::new(self, grade, &resolution);
        let frame = Frame::new(model, self, grade, hom.shifts(1))?;
        let syzygies = Frame::new(model, self, grade, hom.shifts(2))?;
        let differential = hom.differential(model, 1, &frame, &syzygies)?;
        let basis = if abstract_class {
            let derivations = hom.derivations(model, &frame)?;
            Quotient::new(model, differential, derivations, frame.len())?
                .basis(model, frame.len())?
        } else {
            linear::rref(model, differential, frame.len())?.kernel(frame.len())?
        };
        let row = basis
            .get(index)
            .ok_or_else(|| Error::plain("Verbitsky deformation basis index out of range"))?;
        frame.vector(model, row)
    }

    pub(crate) fn tangent_literal<'a>(
        &self,
        model: &mut Model,
        grade: &Grade,
        entries: impl IntoIterator<Item = (usize, &'a Tensor)>,
    ) -> Result<Vector> {
        let mut value = Vector::default();
        for (i, tensor) in entries {
            value.add_tensor(model, i, tensor)?;
        }
        if value.entries.is_empty() {
            return Ok(value);
        }
        let resolution = self.resolution(model, &grade.resolution())?;
        let hom = Hom::new(self, grade, &resolution);
        let frame = Frame::new(model, self, grade, hom.shifts(1))?;
        let syzygies = Frame::new(model, self, grade, hom.shifts(2))?;
        let differential = hom.differential(model, 1, &frame, &syzygies)?;
        let equations = linear::rref(model, differential, frame.len())?;
        if !certify_cycle(model, &equations, &frame.coordinates(&value)?)? {
            return Err(Error::plain(
                "Verbitsky tangent coordinates fail the syzygy equations",
            ));
        }
        Ok(value)
    }

    pub(crate) fn deformation_class(
        &self,
        model: &mut Model,
        grade: &Grade,
        value: &Vector,
    ) -> Result<Vector> {
        if value.entries.is_empty() {
            return Ok(Vector::default());
        }
        let resolution = self.resolution(model, &grade.resolution())?;
        let hom = Hom::new(self, grade, &resolution);
        let frame = Frame::new(model, self, grade, hom.shifts(1))?;
        let mut row = frame.coordinates(value)?;
        let derivations = hom.derivations(model, &frame)?;
        let boundaries = linear::rref(model, derivations, frame.len())?;
        reduce(model, &boundaries, &mut row)?;
        frame.vector(model, &row)
    }
}

fn solve(
    model: &mut Model,
    mut matrix: Matrix,
    width: usize,
    rhs: &[Rational],
) -> Result<Option<Vec<Rational>>> {
    let mut solution = vec![Rational::zero(); width];
    if rhs.iter().all(|entry| entry.is_zero()) {
        return Ok(Some(solution));
    }
    storage(
        matrix
            .len()
            .saturating_mul(width.saturating_add(1))
            .saturating_mul(2 * LIMBS + 1),
    )?;
    for (row, &entry) in matrix.iter_mut().zip(rhs) {
        row.push(entry);
    }
    let echelon = linear::rref(model, matrix, width + 1)?;
    if echelon.pivots.last() == Some(&width) {
        return Ok(None);
    }
    for (&i, row) in echelon.pivots.iter().zip(echelon.rows) {
        solution[i] = row[width];
    }
    Ok(Some(solution))
}

fn polynomial_dot(model: &mut Model, row: &Vector, column: &Vector) -> Result<Tensor> {
    let mut sum = Tensor::default();
    for (i, value) in &column.entries {
        if let Some(multiplier) = row.entries.get(i) {
            let term = value.cup(model, multiplier)?;
            for (powers, coefficient) in term.terms {
                sum.insert(model, powers, coefficient)?;
            }
        }
    }
    Ok(sum)
}

fn lift_relation(
    model: &mut Model,
    hom: &Hom<'_>,
    catalog: &mut Catalog,
    weight: i128,
    residual: &Tensor,
) -> Result<Vector> {
    if residual.terms.is_empty() {
        return Ok(Vector::default());
    }
    let weight = usize::try_from(weight).map_err(|_| koszul::grading_range())?;
    let source = Axis::new(model, hom.form, catalog, hom.shifts(1), weight)?;
    let target = Axis::new(model, hom.form, catalog, &[0], weight)?;
    let mut matrix = linear::matrix(target.len(), source.len())?;
    for (i, (index, monomial)) in source.entries().iter().enumerate() {
        if let Some(polynomial) = hom.resolution.entry(1, 0, *index)? {
            for (powers, &coefficient) in &polynomial.terms {
                model.spend(1)?;
                let product = multiply_powers(powers, monomial)?;
                matrix[target.locate(0, &product)?][i] = coefficient;
            }
        }
    }
    let mut rhs = vec![Rational::zero(); target.len()];
    for (powers, &coefficient) in &residual.terms {
        rhs[target.locate(0, powers)?] = coefficient.negate();
    }
    let solution = solve(model, matrix, source.len(), &rhs)?.ok_or_else(|| {
        Error::plain("Verbitsky syzygy correction does not lift to the polynomial ideal")
    })?;
    let column = source.vector(model, &solution)?;
    let mut checked = residual.clone();
    let image = hom.resolution.differential(model, 1, &column)?;
    if let Some(value) = image.entries.get(&0) {
        for (powers, &coefficient) in &value.terms {
            checked.insert(model, powers.clone(), coefficient)?;
        }
    }
    if !checked.terms.is_empty() {
        return Err(Error::plain(
            "Verbitsky polynomial syzygy correction failed verification",
        ));
    }
    Ok(column)
}

fn format_columns(model: &mut Model, form: &Form, columns: &[Vector]) -> Result<String> {
    let mut text = "[".to_owned();
    let mut separator = "";
    for (j, column) in columns.iter().enumerate() {
        if !column.entries.is_empty() {
            write!(text, "{separator}({j}, {})", column.text(model, form)?).unwrap();
            check_report(&text)?;
            separator = ", ";
        }
    }
    text.push(']');
    Ok(text)
}

fn lift_columns(model: &mut Model, hom: &Hom<'_>, value: &Vector) -> Result<Vec<Vector>> {
    let mut catalog = Catalog::default();
    let mut columns = Vec::new();
    let mut cells = 0;
    for (j, &degree) in hom.shifts(2).iter().enumerate() {
        let residual = polynomial_dot(model, value, hom.resolution.column(2, j)?)?;
        let weight = hom
            .grade
            .shift
            .checked_add(degree as i128)
            .ok_or_else(koszul::grading_range)?;
        let correction = lift_relation(model, hom, &mut catalog, weight, &residual)?;
        cells += correction.cells + 8;
        storage(cells)?;
        columns.push(correction);
    }
    Ok(columns)
}

impl Form {
    pub(crate) fn deformation_obstruction(
        &self,
        model: &mut Model,
        grade: &Grade,
        g: &Vector,
    ) -> Result<String> {
        if g.entries.is_empty() {
            return Ok(
                "residual = []; Ext1 class = []; lift verified; h = []; B = []; C = []".to_owned(),
            );
        }
        let resolution = self.resolution(model, &grade.resolution())?;
        let hom = Hom::new(self, grade, &resolution);
        let b = lift_columns(model, &hom, g)?;
        let doubled = Grade {
            shift: grade
                .shift
                .checked_mul(2)
                .ok_or_else(koszul::grading_range)?,
            ..grade.clone()
        };
        let hom2 = Hom::new(self, &doubled, &resolution);
        let complex = Complex::new(model, &hom2)?;

        let mut products = Vec::new();
        let mut cells = 0;
        let mut residual = Vector::default();
        for (j, correction) in b.iter().enumerate() {
            let raw = polynomial_dot(model, g, correction)?;
            let reduced = hom2.project(model, 2, j, &raw)?;
            residual.add_tensor(model, j, &reduced)?;
            cells += raw.cells + 8;
            storage(cells)?;
            // Keep the raw polynomial residual until C is verified over S.
            products.push(raw);
        }
        let coords = complex.syzygies.coordinates(&residual)?;
        if !certify_cycle(model, &complex.obstruction_space.constraints, &coords)? {
            return Err(Error::plain(
                "Verbitsky quadratic obstruction failed its syzygy cocycle check",
            ));
        }
        let mut class = coords.clone();
        complex.obstruction_space.canonical(model, &mut class)?;
        let representative = complex.syzygies.vector(model, &class)?;
        let mut text = format!(
            "residual = {}; Ext1 class = {}",
            residual.text(model, self)?,
            representative.text(model, self)?,
        );
        if !representative.entries.is_empty() {
            write!(
                text,
                "; obstructed; B = {}",
                format_columns(model, self, &b)?
            )
            .unwrap();
        } else {
            let rhs = coords.iter().map(|c| c.negate()).collect::<Vec<_>>();
            let h = solve(model, complex.differential, complex.generators.len(), &rhs)?
                .ok_or_else(|| {
                    Error::plain("Verbitsky quadratic residual reduction failed verification")
                })?;
            let h = complex.generators.vector(model, &h)?;
            let mut catalog = Catalog::default();
            let mut c = Vec::new();
            let mut cells = 0;
            for (j, mut raw) in products.into_iter().enumerate() {
                let ha = polynomial_dot(model, &h, resolution.column(2, j)?)?;
                for (powers, coefficient) in ha.terms {
                    raw.insert(model, powers, coefficient)?;
                }
                let weight = doubled
                    .shift
                    .checked_add(resolution.shift(2, j)? as i128)
                    .ok_or_else(koszul::grading_range)?;
                let column = lift_relation(model, &hom2, &mut catalog, weight, &raw)?;
                cells += column.cells + 8;
                storage(cells)?;
                c.push(column);
            }
            write!(
                text,
                "; lift verified; h = {}; B = {}; C = {}",
                h.text(model, self)?,
                format_columns(model, self, &b)?,
                format_columns(model, self, &c)?,
            )
            .unwrap();
        }
        check_report(&text)?;
        Ok(text)
    }
}

impl Form {
    pub(in crate::hyperelliptic_kummer_surface::hodge_wick::ring) fn deformation_product_corrections(
        &self,
        model: &mut Model,
        grade: &Grade,
        tangent: &Vector,
        requests: &[(usize, Tensor)],
    ) -> Result<Vec<Tensor>> {
        let resolution = self.resolution(model, &grade.resolution())?;
        let hom = Hom::new(self, grade, &resolution);
        let mut catalog = Catalog::default();
        let mut products = Vec::with_capacity(requests.len());
        let mut cells = 0;
        for (weight, residual) in requests {
            let coefficients = lift_relation(model, &hom, &mut catalog, *weight as i128, residual)?;
            let correction = polynomial_dot(model, tangent, &coefficients)?;
            let slot = grade.slot(*weight).ok_or_else(koszul::grading_range)?;
            let reduced = self.project(model, &slot, &correction)?;
            cells += reduced.cells + 8;
            storage(cells)?;
            products.push(reduced);
        }
        Ok(products)
    }
}
