use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_lens;
use crate::syntax::{Label, ParsedExpr};

use super::super::super::matrix::Field;
use super::{check_categories, Ambient, Category, Double, Image, Operand, Pair};

mod branch;
mod coherence;
mod movement;
use branch::Bulk;
use movement::{Central, Map};

pub(crate) use branch::fusion as validate_branching;
pub(crate) use branch::self_boundary;

pub(crate) enum Command {
    Validate(Label, Option<[usize; 2]>),
    Unit(Label, usize),
    Whisker([Label; 4], usize, Box<[ParsedExpr; 2]>),
    Data(Box<ParsedExpr>),
    Image([usize; 2], Box<ParsedExpr>),
    Basis([usize; 3], Box<ParsedExpr>),
}

pub(crate) enum Expression {
    Induce(Label, usize),
    Move([Label; 2], usize, Option<usize>, Box<ParsedExpr>),
    Entry([usize; 3], bool, Box<ParsedExpr>),
}

pub(crate) fn induction(
    model: &mut Model,
    double: &Double,
    category: &Category,
    simple: usize,
) -> Result<Value> {
    let ambient = Ambient::new(model, double, &[category], [])?;
    Bulk::new(model, &ambient, category, simple)?.value(model, &ambient, category)
}

pub(crate) fn unit(
    model: &mut Model,
    double: &Double,
    category: &Category,
    simple: usize,
) -> Result<String> {
    let ambient = Ambient::new(model, double, &[category], [])?;
    coherence::unit(model, &ambient, category, simple)
}

pub(crate) fn check_path(boundaries: &[&Category], defects: &[&Category]) -> Result<()> {
    debug_assert_eq!(boundaries.len(), defects.len() + 1);
    for &cat in boundaries {
        self_boundary(cat)?;
    }
    for ([left, right], &actual) in boundaries
        .windows(2)
        .map(|cats| [cats[0], cats[1]])
        .zip(defects)
    {
        if !Rc::ptr_eq(&actual.boundaries[0], &left.boundaries[0])
            || !Rc::ptr_eq(&actual.boundaries[1], &right.boundaries[0])
        {
            return Err(Error::plain(
                "central map requires the exact ordered defect boundaries",
            ));
        }
    }
    Ok(())
}

pub(crate) fn move_conductor(
    model: &mut Model,
    double: &Double,
    boundaries: [&Category; 2],
    defect: &Category,
) -> Result<usize> {
    check_path(&boundaries, &[defect])?;
    Ok(Ambient::new(model, double, &[boundaries[0], boundaries[1], defect], [])?.conductor)
}

pub(crate) fn movement(
    model: &mut Model,
    double: &Double,
    boundaries: [&Category; 2],
    defect: (&Category, Operand),
    color: [usize; 2],
    names: [Rc<str>; 2],
) -> Result<Value> {
    check_path(&boundaries, &[defect.0])?;
    let [bulk, conductor] = color;
    if defect.1.conductor().is_some_and(|n| conductor % n != 0) {
        return Err(Error::plain(
            "central map conductor does not contain the graded operand's phases",
        ));
    }
    let ambient = Ambient::new(
        model,
        double,
        &[boundaries[0], boundaries[1], defect.0],
        [conductor],
    )?;
    let bh = Bulk::new(model, &ambient, boundaries[0], bulk)?;
    let bk = Bulk::new(model, &ambient, boundaries[1], bulk)?;
    let index = match &defect.1 {
        Operand::Simple(i) => Some(*i),
        Operand::Graded(_) => None,
    };
    let simple = defect.1.module(model, &ambient, defect.0)?;
    let crossing = Central::new(model, &ambient, [&bh, &bk], &simple)?;
    let bases = crossing.products.map(|p| p.bases);
    let value = CentralMap {
        conductor: ambient.conductor,
        bulk,
        defect: index,
        names,
        maps: crossing.maps,
        bases,
    };
    value.persist(model)
}

pub(crate) fn whisker(
    model: &mut Model,
    double: &Double,
    boundaries: [&Category; 3],
    output: &Category,
    defects: [(&Category, Operand); 2],
    bulk: usize,
) -> Result<String> {
    let [(m, i), (n, j)] = defects;
    check_path(&boundaries, &[m, n])?;
    check_categories(&[m, n], output)?;
    let ambient = Ambient::new(
        model,
        double,
        &[boundaries[0], boundaries[1], boundaries[2], m, n, output],
        [i.conductor(), j.conductor()].into_iter().flatten(),
    )?;
    let [h, k, l] = boundaries;
    let bh = Bulk::new(model, &ambient, h, bulk)?;
    let bk = Bulk::new(model, &ambient, k, bulk)?;
    let bl = Bulk::new(model, &ambient, l, bulk)?;
    let m = i.module(model, &ambient, m)?;
    let n = j.module(model, &ambient, n)?;
    coherence::whisker(model, &ambient, [&bh, &bk, &bl], [&m, &n])
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct CentralMap {
    conductor: usize,
    bulk: usize,
    defect: Option<usize>,
    names: [Rc<str>; 2],
    maps: [Map; 2],
    bases: [Vec<Vec<Pair>>; 2],
}

impl CentralMap {
    fn persist(self, model: &mut Model) -> Result<Value> {
        let cells = self
            .maps
            .iter()
            .map(|m| m.0.iter().map(|r| 1 + r.len() * 2).sum::<usize>())
            .sum::<usize>()
            + self
                .bases
                .iter()
                .map(|v| v.iter().map(|r| 1 + r.len() * 4).sum::<usize>())
                .sum::<usize>();
        model.allocate(cells + 12)?;
        Ok(Value::DoubleCentralMap(Rc::new(self)))
    }

    pub(crate) fn data(&self) -> String {
        let nonempty = self.maps[0]
            .0
            .iter()
            .enumerate()
            .filter_map(|(g, r)| (!r.is_empty()).then_some((g, r.len())))
            .collect::<Vec<_>>();
        let dimension = nonempty.iter().map(|(_, d)| d).sum::<usize>();
        let (phases, moves) = self.maps[0].statistics();
        let [h, k] = &self.names;
        let defect = self
            .defect
            .map_or_else(|| "graded".to_owned(), |x| x.to_string());
        format!("rvDoubleCentralData = boundaries=[{h},{k}], bulk={}, defect={defect}, dimension={dimension}, conductor={}, grades={nonempty:?}, nonzero_phases={phases}, moved_columns={moves}",
            self.bulk,self.conductor)
    }

    fn lookup(&self, inverse: bool, g: usize, c: usize) -> Result<Image> {
        self.maps[usize::from(inverse)]
            .0
            .get(g)
            .and_then(|r| r.get(c))
            .copied()
            .ok_or_else(|| Error::plain("central map grade or column is out of range"))
    }

    pub(crate) fn image(&self, [g, i]: [usize; 2]) -> Result<String> {
        let x = self.lookup(false, g, i)?;
        Ok(format!(
            "rvDoubleCentralImage[{g},{i}] = target={}, phase={}, conductor={}",
            x.target, x.phase, self.conductor
        ))
    }

    pub(crate) fn basis(&self, [side, g, i]: [usize; 3]) -> Result<String> {
        let p = self
            .bases
            .get(side)
            .and_then(|v| v.get(g))
            .and_then(|r| r.get(i))
            .ok_or_else(|| {
                Error::plain("central map basis side, grade or column is out of range")
            })?;
        Ok(format!(
            "rvDoubleCentralBasis[{side},{g},{i}] = grades={:?}, columns={:?}",
            p.grades, p.columns
        ))
    }

    pub(crate) fn entry(
        &self,
        model: &mut Model,
        inverse: bool,
        [g, r, c]: [usize; 3],
    ) -> Result<Value> {
        self.lookup(inverse, g, r)?;
        let x = self.lookup(inverse, g, c)?;
        if x.target != r {
            return Ok(rv_lens::zero());
        }
        let mut field = Field::new(model, self.conductor);
        let root = field.phase(x.phase)?;
        rv_lens::rational_value(model, root)
    }
}
