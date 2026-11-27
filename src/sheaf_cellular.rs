use std::collections::{BTreeMap, HashMap, HashSet};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Poset, TypeId, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::{self, Basis, Coordinate};
use crate::sheaf::{self, Cochains, Cohomology, Generator, Layer, Place, Sink};
use crate::sheaf_hyper::Complex;
use crate::sheaf_page::{FilteredCochains, Grade, Page, Pages};
use crate::syntax::Label;

pub(crate) struct Presentation {
    pub(crate) cells: Vec<(Label, usize, Label)>,
    pub(crate) boundaries: Vec<(Label, Vec<(Label, i64)>)>,
    pub(crate) augmentation: Vec<(Label, i64)>,
}

struct Cell {
    label: String,
    degree: usize,
    point: usize,
}

pub(crate) struct Resolution {
    base: TypeId,
    prime: Prime,
    cells: Vec<Cell>,
    lookup: HashMap<String, usize>,
    degrees: BTreeMap<usize, Vec<usize>>,
    coboundary: Vec<Vec<Coordinate>>,
}

impl Resolution {
    pub(crate) fn build(
        model: &mut Model,
        base: TypeId,
        prime: Prime,
        objects: Presentation,
    ) -> Result<Self> {
        let poset = model.poset(base).expect("resolved finite poset");
        let mut cells = Vec::new();
        let mut lookup = HashMap::new();
        let mut degrees: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        let mut stored = 0;
        for (label, degree, point) in objects.cells {
            model.spend(1)?;
            stored += 10 + label.text.len();
            sheaf::check_cells(stored)?;
            if lookup.insert(label.text.clone(), cells.len()).is_some() {
                return Err(collision(&label, "generator"));
            }
            degrees.entry(degree).or_default().push(cells.len());
            cells.push(Cell {
                label: label.text,
                degree,
                point: model.point(base, &point)?,
            });
        }

        let mut boundaries = vec![Vec::new(); cells.len()];
        let mut seen = HashSet::new();
        for (label, entries) in objects.boundaries {
            model.spend(1)?;
            let index = locate(&lookup, &label)?;
            if !seen.insert(index) {
                return Err(collision(&label, "boundary"));
            }
            let source = &cells[index];
            let Some(previous) = source.degree.checked_sub(1) else {
                return Err(Error::new(
                    label.span,
                    "a degree-zero generator has no boundary",
                ));
            };
            let mut column = BTreeMap::new();
            for (term, value) in entries {
                model.spend(1 + column.len().max(1).ilog2() as usize)?;
                let row = locate(&lookup, &term)?;
                let target = &cells[row];
                if target.degree != previous || !at_or_below(&poset, target.point, source.point) {
                    return Err(Error::new(term.span, format!(
                        "boundary term '{}' must have degree {previous} and lie at or below the point of '{}'",
                        term.text, source.label
                    )));
                }
                let scalar = sheaf::residue(prime, value);
                let coefficient = prime.add(column.get(&row).copied().unwrap_or(0), scalar);
                if coefficient == 0 {
                    column.remove(&row);
                } else {
                    column.insert(row, coefficient);
                    sheaf::check_cells(stored + 4 * column.len())?;
                }
            }
            stored += 4 * column.len();
            boundaries[index] = column.into_iter().collect();
        }
        seen.clear();
        let mut augmentation = vec![0; cells.len()];
        for (label, value) in objects.augmentation {
            model.spend(1)?;
            let index = locate(&lookup, &label)?;
            if !seen.insert(index) {
                return Err(collision(&label, "augmentation entry"));
            }
            if cells[index].degree != 0 {
                return Err(Error::new(
                    label.span,
                    "the augmentation only accepts degree-zero generators",
                ));
            }
            augmentation[index] = sheaf::residue(prime, value);
        }

        for (index, column) in boundaries.iter().enumerate() {
            model.spend(1 + column.len())?;
            if cells[index].degree == 1 {
                let sum = column.iter().fold(0, |sum, &(row, factor)| {
                    prime.add(sum, prime.mul(factor, augmentation[row]))
                });
                if sum != 0 {
                    return Err(Error::plain(format!(
                        "sheaf resolution augmentation of the boundary of '{}' is nonzero",
                        cells[index].label
                    )));
                }
            } else if cells[index].degree > 1 {
                let mut square = Vec::new();
                for &(row, factor) in column {
                    square =
                        prime_linear::add_scaled(model, prime, &square, &boundaries[row], factor)?;
                }
                if !square.is_empty() {
                    return Err(Error::plain(format!(
                        "sheaf resolution boundary squared is nonzero at '{}'",
                        cells[index].label
                    )));
                }
            }
        }

        for (point, name) in poset.names.iter().enumerate() {
            let mut ranks = BTreeMap::new();
            for (&degree, group) in &degrees {
                let mut rank = Basis::default();
                let mut dimension = 0;
                for &index in group {
                    model.spend(1)?;
                    if !at_or_below(&poset, cells[index].point, point) {
                        continue;
                    }
                    dimension += 1;
                    if degree == 0 {
                        if augmentation[index] != 0 && rank.rank() == 0 {
                            rank.insert(model, prime, vec![(0, augmentation[index])])?;
                        }
                    } else {
                        rank.insert(model, prime, boundaries[index].clone())?;
                    }
                }
                if dimension != 0 {
                    ranks.insert(degree, (dimension, rank.rank()));
                }
            }
            if ranks.get(&0).is_none_or(|&(_, rank)| rank == 0) {
                return Err(Error::plain(format!(
                    "sheaf resolution augmentation is not surjective at '{name}'"
                )));
            }
            for (&degree, &(dimension, outgoing)) in &ranks {
                let incoming = degree.checked_add(1).and_then(|next| ranks.get(&next));
                let boundaries = incoming.map_or(0, |&(_, rank)| rank);
                if dimension != outgoing + boundaries {
                    return Err(Error::plain(format!(
                        "sheaf resolution is not exact in degree {degree} at '{name}'"
                    )));
                }
            }
        }

        let mut coboundary = vec![Vec::new(); cells.len()];
        for (target, column) in boundaries.iter().enumerate() {
            model.spend(1 + column.len())?;
            for &(source, coefficient) in column {
                coboundary[source].push((target, coefficient));
            }
        }
        model.allocate(stored)?;
        Ok(Self {
            base,
            prime,
            cells,
            lookup,
            degrees,
            coboundary,
        })
    }
}

pub(crate) struct Cellular {
    resolution: Rc<Resolution>,
    coefficient: Rc<Complex>,
    height: usize,
    cohomology: Cohomology,
    pages: Pages,
}

impl Cellular {
    pub(crate) fn new(
        model: &mut Model,
        resolution: Rc<Resolution>,
        coefficient: Rc<Complex>,
    ) -> Result<Self> {
        let first = &coefficient.components()[0].sheaf;
        if first.base() != resolution.base || first.prime != resolution.prime {
            return Err(Error::plain("shcellular requires the same named base and prime field for the resolution and complex"));
        }
        let mut height = 0;
        for (q, part) in coefficient.components().iter().enumerate() {
            model.spend(resolution.cells.len())?;
            for cell in &resolution.cells {
                if part.sheaf.dims()[cell.point] != 0 {
                    let top = cell
                        .degree
                        .checked_add(q)
                        .and_then(|d| d.checked_add(1))
                        .ok_or_else(|| {
                            Error::plain("cellular sheaf total degree is out of range")
                        })?;
                    height = height.max(top);
                }
            }
        }
        model.allocate(6)?;
        Ok(Self {
            resolution,
            coefficient,
            height,
            cohomology: Cohomology::default(),
            pages: Pages::default(),
        })
    }

    pub(crate) fn term(
        &self,
        model: &mut Model,
        vertical: usize,
        cell: &Label,
        index: usize,
        scalar: i64,
    ) -> Result<(usize, Value)> {
        let ordinal = locate(&self.resolution.lookup, cell)?;
        let data = &self.resolution.cells[ordinal];
        let dimensions = self.coefficient.component(vertical)?.sheaf.dims();
        if index >= dimensions[data.point] {
            return Err(Error::new(
                cell.span,
                format!(
                    "stalk coordinate {index} is out of range for cell '{}' of dimension {}",
                    cell.text, dimensions[data.point]
                ),
            ));
        }
        let total = vertical + data.degree;
        let coefficient = sheaf::residue(self.prime(), scalar);
        if coefficient == 0 {
            return Ok((total, Value::SheafC([].into())));
        }
        let generator = Generator {
            vertical,
            place: Place::Cell {
                degree: data.degree,
                ordinal,
            },
            index,
        };
        Ok((
            total,
            sheaf::persist(model, vec![(generator, coefficient)])?,
        ))
    }
}

impl Cochains for Cellular {
    fn prime(&self) -> Prime {
        self.resolution.prime
    }

    fn cofaces(
        &self,
        model: &mut Model,
        generator: &Generator,
        coefficient: u32,
        emit: &mut Sink<'_>,
    ) -> Result<()> {
        let Place::Cell { degree, ordinal } = generator.place else {
            unreachable!()
        };
        let q = generator.vertical;
        let point = self.resolution.cells[ordinal].point;
        let part = &self.coefficient.components()[q].sheaf;
        let prime = self.prime();
        for &(target, scalar) in &self.resolution.coboundary[ordinal] {
            let next = &self.resolution.cells[target];
            let identity = [(generator.index, 1)];
            let rows = if point == next.point {
                &identity[..]
            } else {
                part.restriction(point, next.point).column(generator.index)
            };
            model.spend(1 + rows.len())?;
            let factor = prime.mul(coefficient, scalar);
            for &(index, scalar) in rows {
                let term = Generator {
                    vertical: q,
                    place: Place::Cell {
                        degree: next.degree,
                        ordinal: target,
                    },
                    index,
                };
                emit(model, term, prime.mul(factor, scalar))?;
            }
        }
        if let Some(map) = self.coefficient.maps().get(q) {
            let rows = map.matrix(point).column(generator.index);
            model.spend(1 + rows.len())?;
            let factor = if degree.is_multiple_of(2) {
                coefficient
            } else {
                prime.negate(coefficient)
            };
            for &(index, scalar) in rows {
                let term = Generator {
                    vertical: q + 1,
                    place: generator.place.clone(),
                    index,
                };
                emit(model, term, prime.mul(factor, scalar))?;
            }
        }
        Ok(())
    }

    fn layer(&self, model: &mut Model, degree: usize, limit: usize) -> Result<Layer> {
        self.layer_window(model, degree, 0, degree, limit)
    }
}

impl FilteredCochains for Cellular {
    fn has_leading(&self, p: usize, q: usize) -> bool {
        q < self.coefficient.components().len() && self.resolution.degrees.contains_key(&p)
    }

    fn page(&self, model: &mut Model, grade: Grade) -> Result<Rc<Page>> {
        self.pages.get(model, self, grade)
    }

    fn layer_window(
        &self,
        model: &mut Model,
        degree: usize,
        lo: usize,
        hi: usize,
        limit: usize,
    ) -> Result<Layer> {
        let mut output = Layer::default();
        if degree >= self.height || lo > hi {
            return Ok(output);
        }
        let first = degree.saturating_sub(hi);
        let Some(last) = degree
            .checked_sub(lo)
            .map(|q| q.min(self.coefficient.components().len() - 1))
        else {
            return Ok(output);
        };
        for q in first..=last {
            model.spend(1)?;
            let Some(cells) = self.resolution.degrees.get(&(degree - q)) else {
                continue;
            };
            let dims = self.coefficient.components()[q].sheaf.dims();
            let layer = Layer::cellular(
                model,
                degree - q,
                cells
                    .iter()
                    .map(|&index| (index, dims[self.resolution.cells[index].point])),
                limit,
            )?;
            output.append_layer(model, layer, q, limit)?;
        }
        Ok(output)
    }
}

#[derive(Clone)]
pub(crate) enum Total {
    Roos(Rc<Complex>),
    Cellular(Rc<Cellular>),
}

impl Total {
    pub(crate) fn filtered(&self) -> &dyn FilteredCochains {
        match self {
            Self::Roos(complex) => complex.as_ref(),
            Self::Cellular(complex) => complex.as_ref(),
        }
    }

    pub(crate) fn cohomology(&self) -> &Cohomology {
        match self {
            Self::Roos(complex) => &complex.cohomology,
            Self::Cellular(complex) => &complex.cohomology,
        }
    }

    pub(crate) fn height(&self) -> usize {
        match self {
            Self::Roos(complex) => complex.height,
            Self::Cellular(complex) => complex.height,
        }
    }

    pub(crate) fn poset(&self) -> &Rc<Poset> {
        match self {
            Self::Roos(complex) => complex.poset(),
            Self::Cellular(complex) => complex.coefficient.poset(),
        }
    }

    pub(crate) fn cell_label(&self, ordinal: usize) -> &str {
        match self {
            Self::Cellular(complex) => &complex.resolution.cells[ordinal].label,
            Self::Roos(_) => unreachable!("expected a cellular total complex"),
        }
    }

    pub(crate) fn rank(&self, model: &mut Model, degree: usize) -> Result<usize> {
        if degree >= self.height() {
            return Ok(0);
        }
        self.cohomology().rank(model, self.filtered(), degree)
    }
}

fn at_or_below(poset: &Poset, point: usize, above: usize) -> bool {
    point == above || poset.strict_less(point, above)
}

fn locate(lookup: &HashMap<String, usize>, label: &Label) -> Result<usize> {
    lookup.get(&label.text).copied().ok_or_else(|| {
        Error::new(
            label.span,
            format!("unknown sheaf resolution generator '{}'", label.text),
        )
    })
}

fn collision(label: &Label, what: &str) -> Error {
    Error::new(
        label.span,
        format!("duplicate sheaf resolution {what} '{}'", label.text),
    )
}
