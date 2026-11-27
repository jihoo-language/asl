use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::{Basis, Coordinate, QuotientBasis};
use crate::sheaf::{self, Matrix};
use crate::syntax::{Label, ParsedExpr};

use super::braid::{Sequence, SequenceMap};
use super::linear;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Form {
    Nearby,
    Raw,
    Class,
}

impl Form {
    pub(crate) fn detect(text: &str) -> Option<Self> {
        Some(match text {
            "ShParVec" => Self::Nearby,
            "ShParRaw" => Self::Raw,
            "ShPar" => Self::Class,
            _ => return None,
        })
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Nearby => "ShParVec",
            Self::Raw => "ShParRaw",
            Self::Class => "ShPar",
        }
    }
}

pub(crate) struct Parabolic {
    pub(crate) sequence: Rc<Sequence>,
    charts: Vec<linear::Image>,
    starts: Vec<usize>,
    boundary: Matrix,
    sum: Matrix,
    quotient: QuotientBasis,
}

impl Parabolic {
    pub(crate) fn build(model: &mut Model, sequence: Rc<Sequence>) -> Result<Self> {
        let p = sequence.prime;
        let width = sequence.nearby;
        let mut prefix = Matrix::identity(model, width)?;
        model.spend(prefix.cells().min(sequence.total.cells()))?;
        if sequence.total != prefix {
            return Err(Error::plain(
                "parabolic cohomology requires identity total monodromy",
            ));
        }
        let mut charts = Vec::with_capacity(sequence.terms.len());
        let mut starts = Vec::with_capacity(sequence.terms.len() + 1);
        let mut sum_columns = Vec::new();
        let mut size = 0;
        let mut cells = 4;
        for disk in &sequence.terms {
            let local = linear::one_minus(model, p, width, &disk.t)?;
            let chart = linear::image(model, p, width, &local)?;
            let transported = Matrix::compose(model, p, &prefix, &chart.inclusion)?;
            cells += transported.cells()
                + chart.inclusion.cells()
                + chart.projection.cells()
                + chart.basis.cells()
                + 2;
            sheaf::check_cells(cells)?;
            sum_columns.extend(transported.columns().map(|(j, col)| (size + j, col.into())));
            starts.push(size);
            size += chart.basis.rank();
            sheaf::check_cells(size)?;
            charts.push(chart);
            prefix = Matrix::compose(model, p, &prefix, &disk.t)?;
        }
        starts.push(size);
        let sum = Matrix::from_columns(sum_columns)?;
        let boundary = Matrix::from_action(model, width, |model, index| {
            model.spend(charts.len())?;
            let mut result = Vec::new();
            for (start, chart) in starts.iter().zip(&charts) {
                let column = chart.projection.column(index);
                model.spend(column.len())?;
                result.extend(column.iter().map(|&(j, x)| (start + j, x)));
            }
            Ok(result)
        })?;
        cells += boundary.cells();
        sheaf::check_cells(cells)?;
        if !Matrix::compose(model, p, &sum, &boundary)?.is_zero() {
            return Err(Error::plain(
                "parabolic boundary failed the total cocycle constraint",
            ));
        }
        let mut boundaries = Basis::default();
        for (index, column) in boundary.columns() {
            let _ = boundaries.insert_with_witness(model, p, column.to_vec(), vec![(index, 1)])?;
            sheaf::check_cells(cells + boundaries.cells())?;
        }
        let quotient = QuotientBasis::build(model, p, boundaries, size, |_, index| {
            Ok(sum.column(index).to_vec())
        })?;
        model.allocate(cells)?;
        Ok(Self {
            sequence,
            charts,
            starts,
            boundary,
            sum,
            quotient,
        })
    }

    fn prime(&self) -> Prime {
        self.sequence.prime
    }

    pub(crate) fn dimension(&self, form: Form) -> usize {
        match form {
            Form::Nearby => self.sequence.nearby,
            Form::Raw => *self.starts.last().expect("final offset"),
            Form::Class => self.quotient.rank(),
        }
    }

    pub(crate) fn literal_vec(&self, model: &mut Model, values: Vec<i64>) -> Result<Value> {
        let entries = signed_vector(model, self.prime(), self.sequence.nearby, values)?;
        crate::odd_ext::persist(model, entries)
    }

    pub(crate) fn literal_rows(&self, model: &mut Model, rows: Vec<Vec<i64>>) -> Result<Value> {
        if rows.len() != self.charts.len() {
            return Err(Error::plain(format!(
                "parabolic raw tuple requires {} strands, found {}",
                self.charts.len(),
                rows.len()
            )));
        }
        let values = rows
            .into_iter()
            .map(|row| signed_vector(model, self.prime(), self.sequence.nearby, row))
            .collect::<Result<Vec<_>>>()?;
        let raw = self.encode_rows(model, values)?;
        crate::odd_ext::persist(model, raw)
    }

    fn encode_rows(
        &self,
        model: &mut Model,
        rows: Vec<Vec<Coordinate>>,
    ) -> Result<Vec<Coordinate>> {
        let mut result = Vec::new();
        for (i, (chart, row)) in self.charts.iter().zip(rows).enumerate() {
            let (remainder, coordinates) =
                chart.basis.reduce_with_witness(model, self.prime(), row)?;
            if !remainder.is_empty() {
                return Err(Error::plain(format!(
                    "parabolic tuple at strand {i} is outside im(I - T_i)"
                )));
            }
            model.spend(coordinates.len())?;
            sheaf::check_cells(2 * (result.len() + coordinates.len()))?;
            result.extend(
                coordinates
                    .into_iter()
                    .map(|(j, x)| (self.starts[i] + j, x)),
            );
        }
        Ok(result)
    }

    fn decode_rows(&self, model: &mut Model, raw: &[Coordinate]) -> Result<Vec<Vec<Coordinate>>> {
        let mut result = Vec::with_capacity(self.charts.len());
        let mut offset = 0;
        let mut cells = 0;
        for (i, chart) in self.charts.iter().enumerate() {
            let limit = self.starts[i + 1];
            let used = raw[offset..].partition_point(|&(index, _)| index < limit);
            let part: Vec<_> = raw[offset..offset + used]
                .iter()
                .map(|&(j, x)| (j - self.starts[i], x))
                .collect();
            model.spend(1 + part.len())?;
            let value = chart.inclusion.apply(model, self.prime(), &part)?;
            cells += 1 + 2 * value.len();
            sheaf::check_cells(cells)?;
            result.push(value);
            offset += used;
        }
        debug_assert_eq!(offset, raw.len());
        Ok(result)
    }

    pub(crate) fn report_rows(&self, model: &mut Model, raw: &[Coordinate]) -> Result<String> {
        let rows = self.decode_rows(model, raw)?;
        let width = self.sequence.nearby;
        sheaf::check_cells((width + 1).saturating_mul(rows.len()))?;
        model.spend((width + 1) * rows.len())?;
        let rows: Vec<_> = rows.iter().map(|row| dense(width, row)).collect();
        Ok(format!("{rows:?}"))
    }

    pub(crate) fn report_vec(&self, model: &mut Model, raw: &[Coordinate]) -> Result<String> {
        let width = self.sequence.nearby;
        model.spend(width)?;
        Ok(format!("{:?}", dense(width, raw)))
    }

    fn project(&self, model: &mut Model, raw: Vec<Coordinate>) -> Result<Vec<Coordinate>> {
        self.quotient
            .project(model, self.prime(), raw)?
            .ok_or_else(|| {
                Error::plain("parabolic tuple is not a cocycle: its transported sum is nonzero")
            })
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        unary: Unary,
        raw: &[Coordinate],
    ) -> Result<Value> {
        let image = match unary {
            Unary::Boundary => self.boundary.apply(model, self.prime(), raw)?,
            Unary::Sum => self.sum.apply(model, self.prime(), raw)?,
            Unary::Class => self.project(model, raw.to_vec())?,
            Unary::Rep => self.quotient.representative(model, self.prime(), raw)?,
            Unary::Primitive => self
                .quotient
                .boundary_witness(model, self.prime(), raw.to_vec())?
                .ok_or_else(|| Error::plain("parabolic tuple is not a coupled coboundary"))?,
        };
        crate::odd_ext::persist(model, image)
    }
}

fn signed_vector(
    model: &mut Model,
    prime: Prime,
    expected: usize,
    values: Vec<i64>,
) -> Result<Vec<Coordinate>> {
    if values.len() != expected {
        return Err(Error::plain(format!(
            "nearby vector requires {expected} coordinates, found {}",
            values.len()
        )));
    }
    model.spend(expected)?;
    Ok(values
        .into_iter()
        .enumerate()
        .filter_map(|(i, x)| {
            let value = sheaf::residue(prime, x);
            (value != 0).then_some((i, value))
        })
        .collect())
}

fn dense(width: usize, support: &[Coordinate]) -> Vec<u32> {
    let mut out = vec![0; width];
    for &(index, value) in support {
        out[index] = value;
    }
    out
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Boundary,
    Sum,
    Class,
    Rep,
    Primitive,
}

impl Unary {
    pub(crate) fn detect(text: &str) -> Option<Self> {
        Some(match text {
            "shParBoundary" => Self::Boundary,
            "shParSum" => Self::Sum,
            "shParClass" => Self::Class,
            "shParRep" => Self::Rep,
            "shParPrimitive" => Self::Primitive,
            _ => return None,
        })
    }

    pub(crate) fn forms(self) -> (Form, Form) {
        match self {
            Self::Boundary => (Form::Nearby, Form::Raw),
            Self::Sum | Self::Primitive => (Form::Raw, Form::Nearby),
            Self::Class => (Form::Raw, Form::Class),
            Self::Rep => (Form::Class, Form::Raw),
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Transport {
    Braid,
    Map,
    Identify,
}

impl Transport {
    pub(crate) fn detect(text: &str) -> Option<Self> {
        Some(match text {
            "shParBraid" => Self::Braid,
            "shParMap" => Self::Map,
            "shParIdentify" => Self::Identify,
            _ => return None,
        })
    }
}

pub(crate) struct Transfer {
    source: Rc<Parabolic>,
    target: Rc<Parabolic>,
    action: TransferKind,
}

enum TransferKind {
    Braid,
    Map(Rc<SequenceMap>),
    Identify,
}

impl Transport {
    pub(crate) fn plan(
        self,
        model: &mut Model,
        label: &Label,
        source: Rc<str>,
        form: Form,
    ) -> Result<(crate::model::TypeId, Transfer)> {
        let (target, action) = match self {
            Self::Braid => {
                let sequence = model.disk_sequence(&label.text)?;
                let braid = sequence
                    .braid
                    .as_ref()
                    .ok_or_else(|| Error::plain("shParBraid requires a declared Hurwitz image"))?;
                if braid.source != source {
                    return Err(Error::plain(
                        "shParBraid input must belong to the declared source sequence",
                    ));
                }
                (label.text.as_str().into(), TransferKind::Braid)
            }
            Self::Map => {
                let map = model.disk_sequence_map(&label.text)?;
                if map.source != source {
                    return Err(Error::plain(
                        "shParMap input must belong to the disk map's source sequence",
                    ));
                }
                (map.target.clone(), TransferKind::Map(map))
            }
            Self::Identify => {
                let lhs = model.disk_sequence(&source)?;
                let rhs = model.disk_sequence(&label.text)?;
                if !lhs.same_data(model, &rhs, false)? {
                    return Err(Error::plain(
                        "shParIdentify needs equal ordered based monodromies over the same field",
                    ));
                }
                (label.text.as_str().into(), TransferKind::Identify)
            }
        };
        let (_, source) = model.parabolic_type(source, form)?;
        let (ty, target) = model.parabolic_type(target, form)?;
        Ok((
            ty,
            Transfer {
                source,
                target,
                action,
            },
        ))
    }
}

impl Transfer {
    pub(crate) fn apply(&self, model: &mut Model, form: Form, raw: &[Coordinate]) -> Result<Value> {
        let p = self.source.prime();
        if form == Form::Nearby {
            let output = match &self.action {
                TransferKind::Map(map) => map.a.apply(model, p, raw)?,
                _ => raw.to_vec(),
            };
            return crate::odd_ext::persist(model, output);
        }
        let input = if form == Form::Class {
            self.source.quotient.representative(model, p, raw)?
        } else {
            raw.to_vec()
        };
        let mut rows = self.source.decode_rows(model, &input)?;
        match &self.action {
            TransferKind::Map(map) => {
                let mut cells = rows.len();
                for row in &mut rows {
                    *row = map.a.apply(model, p, row)?;
                    cells += 2 * row.len();
                    sheaf::check_cells(cells)?;
                }
            }
            TransferKind::Braid => {
                self.target
                    .sequence
                    .braid
                    .as_ref()
                    .expect("validated Hurwitz image")
                    .action(model, &mut rows)?;
            }
            TransferKind::Identify => {}
        }
        let output = self.target.encode_rows(model, rows)?;
        let output = if form == Form::Class {
            self.target.project(model, output)?
        } else {
            output
        };
        crate::odd_ext::persist(model, output)
    }
}

pub(crate) enum Command {
    Sequence(Label, Vec<Label>),
    Hurwitz(Label, Label, Vec<i64>),
    Point(Label, Label, usize),
    Product(Label),
    Same([Label; 2], bool),
    Map(Label, [Label; 2], Vec<ParsedExpr>),
    BraidMap(Label, [Label; 3]),
    MapAt(Label, usize),
    Rank(Label),
    Report(Box<ParsedExpr>, bool),
}

pub(crate) enum Expression {
    Nearby(Label, Vec<i64>),
    Raw(Label, Vec<Vec<i64>>),
    Generator(Label, usize),
    Unary(Unary, Box<ParsedExpr>),
    Transport(Transport, Label, Box<ParsedExpr>),
}
