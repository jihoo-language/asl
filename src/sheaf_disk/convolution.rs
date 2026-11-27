use std::collections::HashSet;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, TypeId, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::{self, Basis, Coordinate, QuotientBasis};
use crate::sheaf::{self, Matrix};
use crate::syntax::{Label, ParsedExpr};

use super::braid::{Sequence, SequenceMap};
use super::{linear, parabolic, Disk, Origin};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub(crate) enum Form {
    Raw,
    Quotient,
}

impl Form {
    pub(crate) fn detect(text: &str) -> Option<Self> {
        match text {
            "ShConvRaw" => Some(Self::Raw),
            "ShConvVec" => Some(Self::Quotient),
            _ => None,
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Raw => "ShConvRaw",
            Self::Quotient => "ShConvVec",
        }
    }
}

struct Operators {
    source: Rc<Sequence>,
    lambda: u32,
    inverse_lambda: u32,
    width: usize,
}

impl Operators {
    fn action(
        &self,
        model: &mut Model,
        generator: Generator,
        raw: &[Coordinate],
    ) -> Result<Vec<Coordinate>> {
        if raw.is_empty() {
            return Ok(Vec::new());
        }
        let n = self.source.nearby;
        let p = self.source.prime;
        let i = generator.index;
        let origin = i * n;
        let local = &self.source.terms[i];
        model.spend(raw.len() + 1)?;
        let start = raw.partition_point(|&(j, _)| j < origin);
        let end = start + raw[start..].partition_point(|&(j, _)| j < origin + n);
        let diagonal = reindex(&raw[start..end], origin);
        let mut correction = Vec::new();
        for group in raw.chunk_by(|x, y| x.0 / n == y.0 / n) {
            let j = group[0].0 / n;
            if i == j {
                continue;
            }
            let raw_local = reindex(group, j * n);
            let transformed = self.source.terms[j].t.apply(model, p, &raw_local)?;
            let difference =
                prime_linear::add_scaled(model, p, &transformed, &raw_local, p.negate(1))?;
            let scale = if j < i { self.lambda } else { 1 };
            correction = prime_linear::add_scaled(model, p, &correction, &difference, scale)?;
        }
        let image = if generator.positive {
            let transformed = local.t.apply(model, p, &diagonal)?;
            prime_linear::add_scaled(model, p, &correction, &transformed, self.lambda)?
        } else {
            let difference =
                prime_linear::add_scaled(model, p, &diagonal, &correction, p.negate(1))?;
            let mut image = local.t_inverse.apply(model, p, &difference)?;
            model.spend(image.len())?;
            for (_, value) in &mut image {
                *value = p.mul(*value, self.inverse_lambda);
            }
            image
        };
        let count = raw.len() - (end - start) + image.len();
        sheaf::check_cells(2 * count)?;
        let mut result = Vec::with_capacity(count);
        result.extend_from_slice(&raw[..start]);
        result.extend(image.into_iter().map(|(j, value)| (origin + j, value)));
        result.extend_from_slice(&raw[end..]);
        Ok(result)
    }
}

fn reindex(group: &[Coordinate], origin: usize) -> Vec<Coordinate> {
    group
        .iter()
        .map(|&(j, value)| (j - origin, value))
        .collect()
}

#[derive(Clone, Copy)]
pub(crate) struct Generator {
    index: usize,
    positive: bool,
}

pub(crate) struct Convolution {
    pub(crate) sequence: Rc<Sequence>,
    source_name: Rc<str>,
    operators: Operators,
    quotient: QuotientBasis,
}

impl Convolution {
    pub(crate) fn build(
        model: &mut Model,
        source_name: Rc<str>,
        source: Rc<Sequence>,
        scalar: i64,
    ) -> Result<Self> {
        let p = source.prime;
        let n = source.nearby;
        let r = source.terms.len();
        let lambda = sheaf::residue(p, scalar);
        if lambda == 0 || lambda == 1 {
            return Err(Error::plain(
                "middle convolution scalar must be different from zero and one in the prime field",
            ));
        }
        let inverse_lambda = prime_linear::invert(model, p, lambda)?;
        let width = n.saturating_mul(r);
        sheaf::check_cells(width)?;
        let operators = Operators {
            source,
            lambda,
            inverse_lambda,
            width,
        };
        let source = &operators.source;
        let mut relations = Basis::default();
        let mut sum_ranks = 0;
        for (i, disk) in source.terms.iter().enumerate() {
            let shift = i * n;
            let difference = linear::one_minus(model, p, n, &disk.t)?;
            let fixed = kernels(model, p, n, |_, j| Ok(difference.column(j).to_vec()))?;
            sum_ranks += n - fixed.len();
            for vector in fixed {
                relations.insert(
                    model,
                    p,
                    vector.into_iter().map(|(j, x)| (shift + j, x)).collect(),
                )?;
            }
        }
        let fixed = kernels(model, p, n, |model, j| {
            prime_linear::add_scaled(
                model,
                p,
                &[(j, p.negate(1))],
                source.total.column(j),
                lambda,
            )
        })?;
        let product_rank = n - fixed.len();
        for mut vector in fixed {
            // Adjacent fixed-row equations give x_i = A_(i+1) x_(i+1).
            // The first remaining equation becomes (lambda A_1 ... A_r - I) x_r = 0.
            let mut suffix = Vec::with_capacity(r);
            let mut cells = r;
            for i in (0..r).rev() {
                model.spend(vector.len() + 1)?;
                cells += 2 * vector.len();
                sheaf::check_cells(cells)?;
                let shift = i * n;
                suffix.push(
                    vector
                        .iter()
                        .map(|&(j, x)| (shift + j, x))
                        .collect::<Vec<_>>(),
                );
                if i != 0 {
                    vector = source.terms[i].t.apply(model, p, &vector)?;
                }
            }
            let lifted = suffix.into_iter().rev().flatten().collect::<Vec<_>>();
            let independent = relations.insert_with_witness(model, p, lifted, Vec::new())?;
            if independent.is_some() {
                return Err(Error::plain(
                    "middle convolution fixed spaces K and L intersect",
                ));
            }
        }
        let expected = (sum_ranks + product_rank).checked_sub(n);
        let quotient = QuotientBasis::build(model, p, relations, width, |_, _| Ok(Vec::new()))?;
        let rank = quotient.rank();
        if expected != Some(rank) {
            return Err(Error::plain(
                "middle convolution quotient has inconsistent rank",
            ));
        }
        let mut terms = Vec::with_capacity(r);
        for index in 0..r {
            let generator = Generator {
                index,
                positive: true,
            };
            let monodromy = Matrix::from_action(model, rank, |model, j| {
                let lift = quotient.representative(model, p, &[(j, 1)])?;
                let raw = operators.action(model, generator, &lift)?;
                coordinate_projection(&quotient, model, p, raw)
            })?;
            let difference = linear::one_minus(model, p, rank, &monodromy)?;
            let (u, v, vanishing) = linear::image_factorization(model, p, rank, &difference)?;
            let disk = Disk::new(model, p, rank, vanishing, u, v, Origin::Declared)?;
            model.spend(monodromy.cells().min(disk.t.cells()))?;
            if monodromy != *disk.t {
                return Err(Error::plain(
                    "middle convolution quiver changed the quotient monodromy",
                ));
            }
            terms.push(Rc::new(disk));
        }
        let sequence = Rc::new(Sequence::build(model, terms)?);
        model.allocate(5)?;
        Ok(Self {
            sequence,
            source_name,
            operators,
            quotient,
        })
    }

    pub(crate) fn dimension(&self, form: Form) -> usize {
        match form {
            Form::Raw => self.operators.width,
            Form::Quotient => self.quotient.rank(),
        }
    }

    pub(crate) fn generator(&self, index: i64) -> Result<Generator> {
        let i = usize::try_from(index.unsigned_abs())
            .ok()
            .and_then(|i| i.checked_sub(1))
            .filter(|&i| i < self.sequence.terms.len())
            .ok_or_else(|| {
                Error::plain(format!(
                    "middle convolution generator {index} is out of range"
                ))
            })?;
        Ok(Generator {
            index: i,
            positive: index > 0,
        })
    }

    pub(crate) fn literal(&self, model: &mut Model, rows: Vec<Vec<i64>>) -> Result<Value> {
        let source = &self.operators.source;
        if rows.len() != source.terms.len() {
            return Err(Error::plain(format!(
                "middle convolution raw vector requires {} blocks, found {}",
                source.terms.len(),
                rows.len(),
            )));
        }
        let mut entries = Vec::new();
        for (i, row) in rows.into_iter().enumerate() {
            let block = parabolic::signed_vector(model, source.prime, source.nearby, row)?;
            entries.extend(block.into_iter().map(|(j, x)| (i * source.nearby + j, x)));
        }
        crate::odd_ext::persist(model, entries)
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        operation: Operation,
        form: Form,
        raw: &[Coordinate],
    ) -> Result<Value> {
        let p = self.sequence.prime;
        let value = match operation {
            Operation::Project => coordinate_projection(&self.quotient, model, p, raw.to_vec())?,
            Operation::Section => self.quotient.representative(model, p, raw)?,
            Operation::Action(generator) if form == Form::Raw => {
                self.operators.action(model, generator, raw)?
            }
            Operation::Action(generator) => {
                let disk = &self.sequence.terms[generator.index];
                let matrix = if generator.positive {
                    &disk.t
                } else {
                    &disk.t_inverse
                };
                matrix.apply(model, p, raw)?
            }
        };
        crate::odd_ext::persist(model, value)
    }

    pub(crate) fn report(
        &self,
        model: &mut Model,
        form: Form,
        raw: &[Coordinate],
    ) -> Result<String> {
        let width = self.operators.source.nearby;
        let strands = self.sequence.terms.len();
        let n = self.dimension(form);
        sheaf::check_cells(n + strands)?;
        model.spend(n + raw.len() + strands)?;
        let entries = parabolic::dense(n, raw);
        Ok(match form {
            Form::Quotient => format!("{entries:?}"),
            Form::Raw if width == 0 => format!("{:?}", vec![Vec::<u32>::new(); strands]),
            Form::Raw => format!("{:?}", entries.chunks_exact(width).collect::<Vec<_>>()),
        })
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Operation {
    Project,
    Section,
    Action(Generator),
}

fn kernels(
    model: &mut Model,
    p: Prime,
    width: usize,
    mut column: impl FnMut(&mut Model, usize) -> Result<Vec<Coordinate>>,
) -> Result<Vec<Vec<Coordinate>>> {
    let mut images = Basis::default();
    let mut kernels = Vec::new();
    let mut cells = 0;
    for j in 0..width {
        let image = column(model, j)?;
        if let Some(kernel) = images.insert_with_witness(model, p, image, vec![(j, 1)])? {
            cells += 1 + 2 * kernel.len();
            sheaf::check_cells(cells + images.cells())?;
            kernels.push(kernel);
        }
    }
    Ok(kernels)
}

fn coordinate_projection(
    space: &QuotientBasis,
    model: &mut Model,
    p: Prime,
    representative: Vec<Coordinate>,
) -> Result<Vec<Coordinate>> {
    space.project(model, p, representative)?.ok_or_else(|| {
        Error::plain("middle convolution projection failed to span the ambient free space")
    })
}

fn diagonal_map(
    model: &mut Model,
    p: Prime,
    matrix: &Matrix,
    input: usize,
    output: usize,
    raw: &[Coordinate],
) -> Result<Vec<Coordinate>> {
    model.spend(raw.len() + 1)?;
    let mut result = Vec::new();
    for group in raw.chunk_by(|x, y| x.0 / input == y.0 / input) {
        let i = group[0].0 / input;
        let image = matrix.apply(model, p, &reindex(group, i * input))?;
        sheaf::check_cells(2 * (result.len() + image.len()))?;
        result.extend(image.into_iter().map(|(j, x)| (i * output + j, x)));
    }
    Ok(result)
}

pub(crate) struct ConvolutionMap {
    pub(crate) descended: Rc<SequenceMap>,
    source: Rc<Convolution>,
    target: Rc<Convolution>,
    original: Rc<SequenceMap>,
}

impl ConvolutionMap {
    pub(crate) fn build(
        model: &mut Model,
        names: (Rc<str>, Rc<str>),
        source: Rc<Convolution>,
        target: Rc<Convolution>,
        original: Rc<SequenceMap>,
    ) -> Result<Self> {
        if original.source != source.source_name
            || original.target != target.source_name
            || source.operators.lambda != target.operators.lambda
            || source.sequence.prime != target.sequence.prime
        {
            return Err(Error::plain(
                "middle convolution map needs the given original endpoints and equal field scalars",
            ));
        }
        let p = source.sequence.prime;
        let n = source.operators.source.nearby;
        let m = target.operators.source.nearby;
        let matrix = Matrix::from_action(model, source.quotient.rank(), |model, j| {
            let lift = source.quotient.representative(model, p, &[(j, 1)])?;
            let mapped = diagonal_map(model, p, &original.a, n, m, &lift)?;
            coordinate_projection(&target.quotient, model, p, mapped)
        })?;
        let descended = Rc::new(SequenceMap::from_nearby_middle(
            model,
            names,
            &source.sequence,
            &target.sequence,
            matrix,
        )?);
        model.allocate(3)?;
        Ok(Self {
            descended,
            source,
            target,
            original,
        })
    }

    pub(crate) fn apply(&self, model: &mut Model, form: Form, raw: &[Coordinate]) -> Result<Value> {
        let p = self.source.sequence.prime;
        let support = match form {
            Form::Raw => diagonal_map(
                model,
                p,
                &self.original.a,
                self.source.operators.source.nearby,
                self.target.operators.source.nearby,
                raw,
            )?,
            Form::Quotient => self.descended.a.apply(model, p, raw)?,
        };
        crate::odd_ext::persist(model, support)
    }
}

pub(crate) fn infinity(model: &mut Model, sequence: &Sequence) -> Result<String> {
    let n = sequence.nearby;
    let monodromy = linear::inverse(model, sequence.prime, n, &sequence.total)?;
    linear::report(model, &monodromy, n, n)
}

pub(crate) fn rigidity(model: &mut Model, sequence: &Sequence) -> Result<i128> {
    let n = sequence.nearby;
    let cells = n.saturating_mul(n);
    sheaf::check_cells(cells)?;
    let p = sequence.prime;
    let matrices = sequence
        .terms
        .iter()
        .map(|disk| disk.t.as_ref())
        .chain(std::iter::once(&sequence.total));
    let mut total = 0i128;
    for matrix in matrices {
        let rows = linear::transpose(model, matrix)?;
        let mut constraints = Basis::default();
        for j in 0..n {
            for i in 0..n {
                model.spend(1 + rows.column(j).len() + matrix.column(i).len())?;
                let left = rows
                    .column(j)
                    .iter()
                    .map(|&(k, x)| (k * n + i, x))
                    .collect::<Vec<_>>();
                let right = matrix
                    .column(i)
                    .iter()
                    .map(|&(k, x)| (j * n + k, x))
                    .collect::<Vec<_>>();
                let commutator = prime_linear::add_scaled(model, p, &left, &right, p.negate(1))?;
                constraints.insert(model, p, commutator)?;
            }
        }
        // The centralizers of the total monodromy and its inverse are equal.
        total += (cells - constraints.rank()) as i128;
    }
    Ok(total + (1 - sequence.terms.len() as i128) * cells as i128)
}

pub(crate) fn group_order(model: &mut Model, sequence: &Sequence) -> Result<usize> {
    let identity = Rc::new(Matrix::identity(model, sequence.nearby)?);
    let mut cells = identity.cells() + 2;
    let mut queue = vec![identity.clone()];
    let mut known = HashSet::from([identity]);
    let mut position = 0;
    while let Some(element) = queue.get(position).cloned() {
        for disk in &sequence.terms {
            let next = Rc::new(Matrix::compose(model, sequence.prime, &disk.t, &element)?);
            model.spend(next.cells() + 1)?;
            if known.insert(next.clone()) {
                cells += 2 + next.cells();
                sheaf::check_cells(cells)?;
                queue.push(next);
            }
        }
        position += 1;
    }
    Ok(queue.len())
}

pub(crate) enum Command {
    Convolve(Label, Label, i64),
    Map(Label, [Label; 3]),
    Rank(Label),
    Sequence(Label, Query),
    Report(Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Query {
    Infinity,
    Rigidity,
    GroupOrder,
}

pub(crate) enum Expression {
    Raw(Label, Vec<Vec<i64>>),
    Generator(Label, usize),
    Project(Box<ParsedExpr>),
    Section(Box<ParsedExpr>),
    Action(Form, i64, Box<ParsedExpr>),
    Map(Label, Box<ParsedExpr>),
}

pub(crate) struct Transfer {
    map: Rc<ConvolutionMap>,
    form: Form,
}

impl Transfer {
    pub(crate) fn plan(
        model: &mut Model,
        label: &Label,
        form: Form,
        source: Rc<str>,
    ) -> Result<(TypeId, Self)> {
        let map = model.convolution_map(&label.text)?;
        if source != map.descended.source {
            return Err(Error::plain(
                "middle convolution map requires its declared source quotient or raw space",
            ));
        }
        let (target, _) = model.convolution_type(map.descended.target.clone(), form)?;
        Ok((target, Self { map, form }))
    }

    pub(crate) fn apply(&self, model: &mut Model, value: &[Coordinate]) -> Result<Value> {
        self.map.apply(model, self.form, value)
    }
}
