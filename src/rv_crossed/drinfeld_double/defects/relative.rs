use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_cyclotomic;
use crate::rv_lens::{self, Rational};
use crate::syntax::{Label, ParsedExpr};

use super::super::representations::{self, IntertwinerOrbit};
use super::{Category, Double, Image};

mod associator;
pub(crate) use associator::associate;
mod graded;
mod tensor;
use graded::{Ambient, Module};
use tensor::{Balanced, Pair};

pub(crate) enum Command {
    Data(Box<ParsedExpr>),
    Basis([usize; 2], Box<ParsedExpr>),
    InducedAction([usize; 4], Box<ParsedExpr>),
    Action([usize; 4], Box<ParsedExpr>),
    Hom(usize, Box<ParsedExpr>),
    Inclusion([usize; 4], Box<ParsedExpr>),
    Projector([usize; 3], Box<ParsedExpr>),
    Associate(Label, Option<[usize; 2]>, Box<[ParsedExpr; 3]>),
}

pub(crate) enum Expression {
    Compose(Label, Box<[ParsedExpr; 2]>),
}

pub(crate) enum Operand {
    Simple(usize),
    Graded(Rc<Relative>),
}

impl Operand {
    fn conductor(&self) -> Option<usize> {
        match self {
            Self::Simple(_) => None,
            Self::Graded(value) => Some(value.conductor),
        }
    }

    fn module(self, model: &mut Model, ambient: &Ambient<'_>, cat: &Category) -> Result<Module> {
        match self {
            Self::Simple(index) => ambient.induce(model, cat, index),
            Self::Graded(value) => {
                let scale = ambient.conductor / value.conductor;
                let ends = ambient.ends(cat);
                let old = &value.module;
                Module::build(model, ambient, ends, old.fibers.clone(), |side, x, g, i| {
                    let k = ambient.boundaries[ends[side]].local(x);
                    let image = old.actions[side][k][g][i];
                    Ok(Image {
                        target: image.target,
                        phase: image.phase * scale,
                    })
                })
            }
        }
    }
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Relative {
    conductor: usize,
    module: Module,
    bases: Vec<Vec<Pair>>,
    homs: Vec<Hom>,
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(super) struct Hom {
    rows: usize,
    cols: usize,
    orbits: Vec<IntertwinerOrbit>,
}

pub(crate) fn check_categories(sources: &[&Category], target: &Category) -> Result<()> {
    let first = sources[0];
    let last = sources.last().unwrap();
    if sources
        .windows(2)
        .any(|pair| !Rc::ptr_eq(&pair[0].boundaries[1], &pair[1].boundaries[0]))
    {
        return Err(Error::plain(
            "relative tensor requires the identical middle boundary",
        ));
    }
    if !Rc::ptr_eq(&first.boundaries[0], &target.boundaries[0])
        || !Rc::ptr_eq(&last.boundaries[1], &target.boundaries[1])
    {
        return Err(Error::plain(
            "relative tensor requires the ordered outer boundary category",
        ));
    }
    Ok(())
}

impl Relative {
    pub(crate) fn new(
        model: &mut Model,
        double: &Double,
        target: &Category,
        inputs: [(&Category, Operand); 2],
    ) -> Result<Self> {
        let [(left, a), (right, b)] = inputs;
        check_categories(&[left, right], target)?;
        let ambient = Ambient::new(
            model,
            double,
            &[left, right, target],
            [a.conductor(), b.conductor()].into_iter().flatten(),
        )?;
        let m = a.module(model, &ambient, left)?;
        let n = b.module(model, &ambient, right)?;
        let product = Balanced::new(model, &ambient, &m, &n)?;
        let homs = decompose(model, &ambient, target, &product.module)?;
        let mut module = product.module;
        module.ends = [0, 1];
        Ok(Self {
            conductor: ambient.conductor,
            module,
            bases: product.bases,
            homs,
        })
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        let cells = self.module.cells()
            + 10
            + self
                .bases
                .iter()
                .map(|basis| 1 + 4 * basis.len())
                .sum::<usize>()
            + self
                .homs
                .iter()
                .map(|h| 7 + h.orbits.iter().map(|o| 3 + 2 * o.len()).sum::<usize>())
                .sum::<usize>();
        model.allocate(cells)?;
        Ok(Value::DoubleRelative(Rc::new(self)))
    }

    pub(crate) fn data(&self) -> String {
        let grades = self
            .module
            .fibers
            .iter()
            .enumerate()
            .filter_map(|(g, &dim)| (dim != 0).then_some((g, dim)))
            .collect::<Vec<_>>();
        let sectors = sector_counts(&self.homs);
        format!("rvDoubleRelativeData = dimension={}, conductor={}, grades={grades:?}, sectors={sectors:?}",self.module.dimension(),self.conductor)
    }

    pub(crate) fn basis(&self, [g, i]: [usize; 2]) -> Result<String> {
        let pair = self
            .bases
            .get(g)
            .and_then(|basis| basis.get(i))
            .ok_or_else(|| Error::plain("relative grade or basis column is out of range"))?;
        Ok(format!(
            "rvDoubleRelativeBasis[{g},{i}] = grades={:?}, columns={:?}",
            pair.grades, pair.columns
        ))
    }

    pub(crate) fn action(
        &self,
        double: &Double,
        target: &Category,
        key: [usize; 4],
    ) -> Result<String> {
        report_action(
            double,
            target,
            &self.module,
            self.conductor,
            key,
            "rvDoubleRelativeAction",
        )
    }

    fn hom(&self, c: usize) -> Result<&Hom> {
        self.homs
            .get(c)
            .ok_or_else(|| Error::plain("relative output simple index is out of range"))
    }

    pub(crate) fn hom_report(&self, c: usize) -> Result<String> {
        let h = self.hom(c)?;
        Ok(format!(
            "rvDoubleRelativeHom[{c}] = multiplicity={}, rows={}, columns={}",
            h.orbits.len(),
            h.rows,
            h.cols
        ))
    }

    pub(crate) fn inclusion(&self, model: &mut Model, [c, b, r, s]: [usize; 4]) -> Result<String> {
        let h = self.hom(c)?;
        if r >= h.rows || s >= h.cols {
            return Err(Error::plain(
                "relative inclusion matrix coordinate is out of range",
            ));
        }
        let orbit = h
            .orbits
            .get(b)
            .ok_or_else(|| Error::plain("relative inclusion index is out of range"))?;
        model.spend(orbit.len())?;
        let coordinate = r * h.cols + s;
        let value = orbit
            .iter()
            .find_map(|&(i, p)| (i == coordinate).then_some(p));
        let entry = phase_entry(model, self.conductor, value, 1)?;
        Ok(format!(
            "rvDoubleRelativeInclusion[{c},{b},{r},{s}] = {entry}"
        ))
    }

    pub(crate) fn projector(&self, model: &mut Model, [c, r, s]: [usize; 3]) -> Result<String> {
        let h = self.hom(c)?;
        if r >= h.rows * h.cols || s >= h.rows * h.cols {
            return Err(Error::plain(
                "relative projector coordinate is out of range",
            ));
        }
        let mut value = None;
        let mut denominator = 1;
        for orbit in &h.orbits {
            model.spend(orbit.len())?;
            let find = |j| orbit.iter().find_map(|&(i, p)| (i == j).then_some(p));
            if let (Some(a), Some(b)) = (find(r), find(s)) {
                value = Some((a + self.conductor - b) % self.conductor);
                denominator = orbit.len();
                break;
            }
        }
        let entry = phase_entry(model, self.conductor, value, denominator)?;
        Ok(format!("rvDoubleRelativeProjector[{c},{r},{s}] = {entry}"))
    }
}

pub(crate) fn induced_action(
    model: &mut Model,
    double: &Double,
    source: &Category,
    simple: usize,
    key: [usize; 4],
) -> Result<String> {
    let ambient = Ambient::new(model, double, &[source], [])?;
    let module = ambient.induce(model, source, simple)?;
    report_action(
        double,
        source,
        &module,
        ambient.conductor,
        key,
        "rvDoubleRelativeInducedAction",
    )
}

fn report_action(
    double: &Double,
    category: &Category,
    module: &Module,
    conductor: usize,
    [side, x, g, i]: [usize; 4],
    op: &str,
) -> Result<String> {
    let boundary = category
        .boundaries
        .get(side)
        .ok_or_else(|| Error::plain("relative action side must be 0 (left) or 1 (right)"))?;
    let k = boundary
        .members()?
        .binary_search(&x)
        .map_err(|_| Error::plain("relative action element is outside its boundary"))?;
    let image = module.actions[side][k]
        .get(g)
        .and_then(|v| v.get(i))
        .ok_or_else(|| Error::plain("relative action grade or basis column is out of range"))?;
    Ok(format!(
        "{op}[{side},{x},{g},{i}] = grade={}, target={}, phase={}, conductor={conductor}",
        double.group.mul(x, g),
        image.target,
        image.phase
    ))
}

fn phase_entry(
    model: &mut Model,
    conductor: usize,
    phase: Option<usize>,
    denominator: usize,
) -> Result<String> {
    let val = if let Some(phase) = phase {
        let mut raw = vec![0; phase + 1];
        raw[phase] = 1;
        let numerator = rv_cyclotomic::reduced(model, conductor, raw)?;
        rv_lens::fraction(model, &numerator, denominator as u128)?
    } else {
        Rational::zero()
    };
    Ok(format!(
        "RVCycloQ({conductor}){{{:?} / {}}}",
        val.numerator, val.denominator
    ))
}

fn sector_counts(homs: &[Hom]) -> Vec<(usize, usize)> {
    homs.iter()
        .enumerate()
        .filter_map(|(i, h)| (!h.orbits.is_empty()).then_some((i, h.orbits.len())))
        .collect()
}

fn decompose(
    model: &mut Model,
    ambient: &Ambient<'_>,
    target: &Category,
    module: &Module,
) -> Result<Vec<Hom>> {
    let mut result = Vec::with_capacity(target.count());
    let scale = ambient.conductor / target.conductor;
    let mut dimension = 0;
    for block in &target.cosets {
        let g = block.support[0];
        let rows = module.fibers[g];
        if block.support.iter().any(|&x| module.fibers[x] != rows) {
            return Err(Error::plain(
                "relative module has inconsistent support dimensions in a coset",
            ));
        }
        let mut width = 0;
        let action = module.transport(model, ambient, target, g)?;
        let theta = block
            .theta
            .iter()
            .map(|row| row.iter().map(|p| p * scale).collect())
            .collect::<Vec<_>>();
        representations::verify_action(model, &target.group, &theta, &action, ambient.conductor)?;
        for i in block.first..block.first + block.count {
            let defect = target.simple(i)?;
            let cols = defect.simple.dimension;
            let rep = defect
                .representation
                .action
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|x| Image {
                            target: x.target,
                            phase: x.phase * scale,
                        })
                        .collect()
                })
                .collect::<Vec<_>>();
            let orbits =
                representations::intertwiner_orbits(model, &action, &rep, ambient.conductor)?;
            width += orbits.len() * cols;
            result.push(Hom { rows, cols, orbits });
        }
        if width != rows {
            return Err(Error::plain(
                "relative multiplicity spaces do not span the balanced fiber",
            ));
        }
        dimension += width * block.support.len();
    }
    if dimension != module.dimension() {
        return Err(Error::plain(
            "relative decomposition has the wrong dimension",
        ));
    }
    Ok(result)
}
