use std::cell::OnceCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::ec_poly::{Field, Scalar};
use crate::error::{Error, Result};
use crate::hw_algebra::CoefficientRing;
use crate::hw_closed::Curve;
use crate::jacobian_level::Basis;
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::super::{field, nikulin, polarized, theta, Point};
use super::{jets, quartic, ResolvedPoint};

pub(crate) mod mordell;

pub(crate) enum Command {
    Define(Label, Label, u8, u8),
    Fiber(Label, Box<ParsedExpr>),
    Census(Label),
}

pub(crate) struct Expression(pub(crate) Label, pub(crate) Operation);

pub(crate) enum Operation {
    Divisor,
    Param(Box<[ParsedExpr; 2]>),
    Map(Box<ParsedExpr>),
    Section(u8, Box<ParsedExpr>),
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Parameter([Scalar; 2]);

impl Parameter {
    pub(crate) fn coordinates(&self) -> &[Scalar; 2] {
        &self.0
    }

    fn new(model: &mut Model, field: Field<'_>, mut coordinates: [Scalar; 2]) -> Result<Self> {
        let Some(first) = coordinates.iter().find(|entry| !field.is_zero(entry)) else {
            return Err(Error::plain(
                "K3 pencil parameter must be a nonzero projective pair",
            ));
        };
        if *first != field.one() {
            let inverse = field.inverse_unit(model, first)?;
            for coordinate in &mut coordinates {
                *coordinate = field.multiply(model, coordinate, &inverse)?;
            }
        }
        Ok(Self(coordinates))
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<std::cmp::Ordering> {
        model.spend(2 * (1 + self.0[0].len()))?;
        Ok(self.cmp(rhs))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(1 + 2 * (3 + self.0[0].len()))?;
        Ok(Value::K3PencilBase(Rc::new(self)))
    }

    fn index(&self, field: Field<'_>) -> usize {
        // This index is used only when the whole base fits the element limit.
        if field.is_zero(&self.0[0]) {
            return 0;
        }
        1 + self.0[1].iter().fold(0, |index, &c| {
            index * field.ring.prime() as usize + c as usize
        })
    }
}

type Forms = [[Scalar; 4]; 2];
type Section = [[Scalar; 2]; 2];

pub(crate) struct Pencil {
    pub(crate) basis_name: Label,
    pub(crate) basis: Rc<Basis>,
    pub(crate) curve: Rc<Curve>,
    masks: [u8; 2],
    nodes: [Point; 2],
    forms: Forms,
    origin: u16,
    section_masks: u16,
    sections: [OnceCell<Section>; 16],
    counts: OnceCell<Box<[usize]>>,
}

impl Pencil {
    pub(crate) fn new(model: &mut Model, basis_name: Label, a: u8, b: u8) -> Result<Self> {
        if a == b {
            return Err(Error::plain("K3 pencil requires two distinct nodes"));
        }
        let basis = nikulin::basis_for(
            model,
            &basis_name,
            "K3 pencils require a certified genus-two level-2 basis",
        )?;
        let curve = model.kummer_curve(&basis.curve_name)?;
        let mut masks = [a, b];
        masks.sort_unstable();
        let left = basis.binary_combination(model, masks[0])?;
        let right = basis.binary_combination(model, masks[1])?;
        let left = Point::project(model, &curve, &left)?;
        let right = Point::project(model, &curve, &right)?;
        let coordinates = [left.coordinates(), right.coordinates()];
        let field = field(&curve);
        let mut midpoint = std::array::from_fn(|_| field.zero());
        for (i, entry) in midpoint.iter_mut().enumerate() {
            *entry = field.add(model, &coordinates[0][i], &coordinates[1][i])?;
        }
        // Both endpoints are double points, hence the restricted quartic
        // is a scalar times s²t². Its value at s=t=1 is that scalar.
        let terms = quartic(model, &curve)?;
        if field.is_zero(&jets::evaluate(model, field, &midpoint, &terms)?) {
            return Err(Error::plain(
                "K3 pencil secant is contained in the Kummer quartic",
            ));
        }
        let forms = annihilator(
            model,
            field,
            [coordinates[0].clone(), coordinates[1].clone()],
        )?;
        let origin = polarized::trope_origin(model, &basis)?;
        let section_masks = (0..16).fold(0u16, |set, mask| {
            let contains = |node: u8| origin & (1u16 << (node ^ mask)) != 0;
            set | (u16::from(contains(masks[0]) ^ contains(masks[1])) << mask)
        });
        if section_masks.count_ones() != 8 {
            return Err(Error::plain("K3 pencil must have eight trope sections"));
        }
        model.allocate(36 + basis_name.text.len() + 16 * field.cells())?;
        Ok(Self {
            basis_name,
            basis,
            curve,
            masks,
            nodes: [left, right],
            forms,
            origin,
            section_masks,
            sections: std::array::from_fn(|_| OnceCell::new()),
            counts: OnceCell::new(),
        })
    }

    pub(crate) fn count(&self) -> usize {
        self.curve
            .ring
            .size_of_residue_field()
            .ok()
            .and_then(|q| usize::try_from(q).ok())
            .map_or(usize::MAX, |q| q.saturating_add(1))
    }

    pub(crate) fn elements(&self, model: &mut Model, limit: usize) -> Result<Vec<Value>> {
        let count = self.count();
        if count > limit {
            return Err(Error::plain(
                "K3 pencil base exceeds the finite element limit",
            ));
        }
        let field = field(&self.curve);
        model.spend(count * field.cells())?;
        let mut elements = Vec::with_capacity(count);
        elements.push(Parameter([field.zero(), field.one()]).persist(model)?);
        for i in 0..count - 1 {
            let parameter = Parameter([field.one(), field.ring.coordinate_row(i)]);
            elements.push(parameter.persist(model)?);
        }
        Ok(elements)
    }

    pub(crate) fn parameter(&self, model: &mut Model, u: Scalar, v: Scalar) -> Result<Parameter> {
        Parameter::new(model, field(&self.curve), [u, v])
    }

    pub(crate) fn divisor(&self, model: &mut Model) -> Result<Value> {
        self.fiber_class(model)?.persist(model)
    }

    fn fiber_class(&self, model: &mut Model) -> Result<polarized::Divisor> {
        let mut coefficients = [0; polarized::COEFFICIENTS];
        coefficients[0] = 2;
        for mask in self.masks {
            coefficients[1 + usize::from(mask)] = -2;
        }
        polarized::Divisor::literal(model, &self.basis, coefficients)
    }

    pub(crate) fn map(&self, model: &mut Model, point: &ResolvedPoint) -> Result<Parameter> {
        let base = point.base().coordinates();
        model.spend(8 * field(&self.curve).cells())?;
        let coordinates = if self.nodes.iter().any(|node| node.coordinates() == base) {
            point
                .direction()
                .ok_or_else(|| Error::plain("K3 pencil base node lacks an exceptional direction"))?
        } else {
            base
        };
        let field = field(&self.curve);
        let first = jets::dot(model, field, self.forms[0].iter().zip(coordinates))?;
        let second = jets::dot(model, field, self.forms[1].iter().zip(coordinates))?;
        Parameter::new(model, field, [first, second])
    }

    pub(crate) fn check_section(&self, mask: u8) -> Result<()> {
        if self.section_masks & (1 << mask) == 0 {
            return Err(Error::plain(
                "K3 pencil section trope must contain exactly one of the two base nodes",
            ));
        }
        Ok(())
    }

    pub(crate) fn section(
        &self,
        model: &mut Model,
        mask: u8,
        parameter: &Parameter,
    ) -> Result<Value> {
        self.check_section(mask)?;
        let group = self.basis.theta(model)?;
        let field = field(&self.curve);
        let cache = &self.sections[usize::from(mask)];
        let columns = if let Some(columns) = cache.get() {
            columns
        } else {
            let zero = field.zero();
            let one = field.one();
            let mut image = |u, v| -> Result<_> {
                let point = trope_parameter(model, &group, mask, u, v)?;
                Ok(self.map(model, &point)?.0)
            };
            let left = image(&one, &zero)?;
            let right = image(&zero, &one)?;
            let diagonal = image(&one, &one)?;
            let det = determinant(model, field, &left, &right)?;
            let alpha = determinant(model, field, &diagonal, &right)?;
            let beta = determinant(model, field, &left, &diagonal)?;
            if [&det, &alpha, &beta].into_iter().any(|d| field.is_zero(d)) {
                return Err(Error::plain(
                    "K3 pencil trope does not give a degree-one section",
                ));
            }
            let mut columns = [left, right];
            for (column, scale) in columns.iter_mut().zip([alpha, beta]) {
                for entry in column {
                    *entry = field.multiply(model, entry, &scale)?;
                }
            }
            model.allocate(4 * field.cells())?;
            cache.get_or_init(|| columns)
        };
        // Cramer's rule gives the inverse projective pair without any
        // division, also at infinity and above exceptional directions.
        let u = determinant(model, field, &parameter.0, &columns[1])?;
        let v = determinant(model, field, &columns[0], &parameter.0)?;
        let point = trope_parameter(model, &group, mask, &u, &v)?;
        if &self.map(model, &point)? != parameter {
            return Err(Error::plain(
                "K3 pencil section failed to recover the input parameter",
            ));
        }
        point.persist(model)
    }

    pub(crate) fn fiber(&self, model: &mut Model, parameter: &Parameter) -> Result<usize> {
        Ok(self.fiber_counts(model)?[parameter.index(field(&self.curve))])
    }

    pub(crate) fn census(&self, model: &mut Model, name: &str) -> Result<String> {
        let counts = self.fiber_counts(model)?;
        model.spend(counts.len() * (1 + counts.len().ilog2() as usize))?;
        let mut histogram = BTreeMap::<usize, usize>::new();
        for &count in counts {
            *histogram.entry(count).or_default() += 1;
        }
        let histogram: Vec<_> = histogram.into_iter().collect();
        let total: usize = counts.iter().sum();
        Ok(format!(
            "rvK3PencilCensus[{name}] = base points={}, total={total}, fiber sizes={histogram:?}",
            counts.len()
        ))
    }

    fn fiber_counts(&self, model: &mut Model) -> Result<&[usize]> {
        if let Some(counts) = self.counts.get() {
            return Ok(counts);
        }
        let ty = model.kummer_k3_type(self.basis.curve_name.clone())?;
        let points = model.elements(ty)?;
        // A rational section meets every rational fiber. This also bounds
        // the storage and indexing by the already enumerated surface.
        let base_count = self.count();
        if base_count > points.len() {
            return Err(Error::plain(
                "K3 pencil has fewer resolved points than its base",
            ));
        }
        let field = field(&self.curve);
        model.spend(points.len() * field.cells())?;
        let mut counts = vec![0; base_count];
        for point in points.iter() {
            let Value::KummerK3(point) = point else {
                unreachable!()
            };
            counts[self.map(model, point)?.index(field)] += 1;
        }
        if counts.contains(&0) {
            return Err(Error::plain("K3 pencil rational section missed a fiber"));
        }
        model.allocate(base_count)?;
        Ok(self.counts.get_or_init(|| counts.into_boxed_slice()))
    }
}

fn trope_parameter(
    model: &mut Model,
    group: &theta::Group,
    mask: u8,
    u: &Scalar,
    v: &Scalar,
) -> Result<ResolvedPoint> {
    let point = ResolvedPoint::base_trope_parameter(model, &group.curve, u, v)?;
    point.projective_image(model, &group.curve, |model, coordinates| {
        group.apply(model, usize::from(mask), coordinates)
    })
}

fn determinant(
    model: &mut Model,
    field: Field<'_>,
    a: &[Scalar; 2],
    b: &[Scalar; 2],
) -> Result<Scalar> {
    let left = field.multiply(model, &a[0], &b[1])?;
    let right = field.multiply(model, &a[1], &b[0])?;
    field.subtract(model, &left, &right)
}

fn annihilator(model: &mut Model, field: Field<'_>, mut rows: Forms) -> Result<Forms> {
    let pivots = rref(model, field, &mut rows)?;
    let mut forms = std::array::from_fn(|_| std::array::from_fn(|_| field.zero()));
    for (form, free) in forms.iter_mut().zip((0..4).filter(|i| !pivots.contains(i))) {
        form[free] = field.one();
        for (row, pivot) in rows.iter().zip(pivots) {
            form[pivot] = field.ring.negate(model, &row[free])?;
        }
    }
    rref(model, field, &mut forms)?;
    Ok(forms)
}

fn rref(model: &mut Model, field: Field<'_>, rows: &mut Forms) -> Result<[usize; 2]> {
    let mut pivots = [0; 2];
    let mut rank = 0;
    for column in 0..4 {
        let Some(row) = (rank..2).find(|&row| !field.is_zero(&rows[row][column])) else {
            continue;
        };
        rows.swap(rank, row);
        let inverse = field.inverse_unit(model, &rows[rank][column])?;
        for coordinate in &mut rows[rank][column..] {
            *coordinate = field.multiply(model, coordinate, &inverse)?;
        }
        let other = 1 - rank;
        let factor = rows[other][column].clone();
        for j in column..4 {
            let change = field.multiply(model, &factor, &rows[rank][j])?;
            rows[other][j] = field.subtract(model, &rows[other][j], &change)?;
        }
        pivots[rank] = column;
        rank += 1;
        if rank == 2 {
            return Ok(pivots);
        }
    }
    Err(Error::plain("K3 pencil linear data must have rank two"))
}
