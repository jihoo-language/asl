use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_lens::{self, Rational};
use crate::syntax::{Label, ParsedExpr};

use super::matrix::Field;
use super::Double;

mod catalog;
mod chain;
mod fourier;
pub(super) use catalog::Catalog;
pub(crate) use chain::Move;

const MAX_GENUS: usize = 64;
const MAX_SUPPORT: usize = 65_536;

pub(crate) type Marking = Vec<usize>;

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct State(Vec<(Marking, Rational)>);

pub(crate) enum Command {
    Dimensions(Label, usize),
    Sewing(Label, [usize; 2]),
    Charge(Label, Marking),
    Chain(Label, Move, Marking),
    Modular(Label),
}

pub(crate) enum Expression {
    Basis(Label, usize, Marking),
    Zero(Label, [usize; 2]),
    Coefficient(Marking, Box<ParsedExpr>),
    Move(Move, Box<ParsedExpr>),
    Cut([usize; 2], Box<ParsedExpr>),
    Gauge(Label, Vec<i128>, Box<ParsedExpr>),
    ToTube(Label, Box<ParsedExpr>),
    FromTube(Box<ParsedExpr>),
    Add(Box<[ParsedExpr; 2]>),
    Scale(Box<[ParsedExpr; 2]>),
    Sew(Box<[ParsedExpr; 2]>),
}

pub(crate) fn check_genus(g: usize) -> Result<()> {
    if g > MAX_GENUS {
        return Err(Error::plain(
            "explicit surface states support at most 64 handles",
        ));
    }
    Ok(())
}

impl State {
    pub(crate) fn zero() -> Self {
        Self(Vec::new())
    }

    pub(crate) fn basis(
        model: &mut Model,
        double: &Double,
        q: usize,
        marks: Marking,
    ) -> Result<Self> {
        check_genus(marks.len() / 2)?;
        double.surface_require_marking(model, q, &marks)?;
        Ok(Self(vec![(marks, Rational::integer(1))]))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        if self.0.len() > MAX_SUPPORT {
            return Err(Error::plain(
                "surface sparse support exceeds 65,536 markings",
            ));
        }
        let cells = self
            .0
            .iter()
            .map(|(m, c)| 5 + m.len() + 2 * c.numerator.len())
            .sum::<usize>();
        model.allocate(cells + 4)?;
        Ok(Value::DoubleSurface(Rc::new(self)))
    }

    pub(crate) fn format(&self) -> String {
        self.0
            .iter()
            .map(|(marks, x)| format!("({marks:?},{:?}/{})", x.numerator, x.denominator))
            .collect::<Vec<_>>()
            .join(", ")
    }

    pub(crate) fn coefficient(&self, marks: &[usize]) -> Rational {
        self.0
            .binary_search_by(|(m, _)| m.as_slice().cmp(marks))
            .ok()
            .map(|i| self.0[i].1.clone())
            .unwrap_or_else(Rational::zero)
    }

    pub(crate) fn add(&self, model: &mut Model, right: &Self) -> Result<Self> {
        let mut values = self.0.iter().cloned().collect::<BTreeMap<_, _>>();
        for (marks, value) in &right.0 {
            let x = values.entry(marks.clone()).or_insert_with(Rational::zero);
            *x = rv_lens::sum(model, x, value)?;
        }
        Ok(Self(
            values.into_iter().filter(|(_, x)| !x.is_zero()).collect(),
        ))
    }

    pub(crate) fn scale(
        &self,
        model: &mut Model,
        double: &Double,
        scale: &Rational,
    ) -> Result<Self> {
        if scale.is_zero() {
            return Ok(Self::zero());
        }
        let mut result = Vec::with_capacity(self.0.len());
        for (m, x) in &self.0 {
            let p = rv_lens::multiply(model, double.conductor, x, scale)?;
            if !p.is_zero() {
                result.push((m.clone(), p));
            }
        }
        Ok(Self(result))
    }

    pub(crate) fn moved(
        &self,
        model: &mut Model,
        double: &Double,
        word: Move,
        q: usize,
    ) -> Result<Self> {
        let catalog = double.surface_catalog(model)?;
        let mut values = Vec::new();
        for (marks, x) in &self.0 {
            model.spend(1 + 4 * marks.len())?;
            let (m, _, phase) = word.apply(double, marks)?;
            if catalog.charge(&m)? != q {
                return Err(Error::plain(
                    "certified integral marking move changes the electric sector",
                ));
            }
            values.push((
                m,
                rv_lens::shifted(model, double.conductor, x, phase as i128)?,
            ));
        }
        values.sort_by(|a, b| a.0.cmp(&b.0));
        Ok(Self(values))
    }

    pub(crate) fn sewn(&self, model: &mut Model, double: &Double, other: &Self) -> Result<Self> {
        let size = self
            .0
            .len()
            .checked_mul(other.0.len())
            .filter(|&n| n <= MAX_SUPPORT)
            .ok_or_else(|| Error::plain("sewn surface support exceeds 65,536 markings"))?;
        model.spend(size)?;
        let mut out = Vec::with_capacity(size);
        for (a, x) in &self.0 {
            for (b, y) in &other.0 {
                out.push((
                    [a.as_slice(), b.as_slice()].concat(),
                    rv_lens::multiply(model, double.conductor, x, y)?,
                ));
            }
        }
        // The left genus is fixed, so concatenation is injective and in lexicographic order.
        Ok(Self(out))
    }

    pub(crate) fn cut(
        &self,
        model: &mut Model,
        double: &Double,
        [cut, q]: [usize; 2],
    ) -> Result<Self> {
        let catalog = double.surface_catalog(model)?;
        catalog.check_charge(q)?;
        let mut out = Vec::new();
        for (m, x) in &self.0 {
            model.spend(m.len() + 1)?;
            if catalog.charge(&m[..2 * cut])? == q {
                out.push((m.clone(), x.clone()));
            }
        }
        Ok(Self(out))
    }

    pub(crate) fn gauge(
        &self,
        model: &mut Model,
        from: &Double,
        to: &Double,
        phases: &[i128],
    ) -> Result<Self> {
        let catalog = from.surface_catalog(model)?;
        let target = to.surface_catalog(model)?;
        let n = catalog.order;
        let conductor = from.conductor as i128;
        if from.conductor != to.conductor
            || n != target.order
            || catalog.characters != target.characters
        {
            return Err(Error::plain(
                "surface gauge requires the same marked group, conductor and electric characters",
            ));
        }
        if phases.len() != (n - 1).pow(2) {
            return Err(Error::plain(
                "surface gauge needs the normalized two-cochain entries",
            ));
        }
        model.spend(n.pow(3))?;
        let mu = |a: usize, b: usize| {
            if a == 0 || b == 0 {
                0
            } else {
                phases[(a - 1) * (n - 1) + (b - 1)].rem_euclid(conductor)
            }
        };
        for a in 0..n {
            for b in 0..n {
                let ab = from.group.mul(a, b);
                if ab != to.group.mul(a, b) {
                    return Err(Error::plain(
                        "surface gauge requires identical marked group tables",
                    ));
                }
                for c in 0..n {
                    let bc = from.group.mul(b, c);
                    let boundary = mu(b, c) - mu(ab, c) + mu(a, bc) - mu(a, b);
                    let delta =
                        to.cocycle_phase([a, b, c]) as i128 - from.cocycle_phase([a, b, c]) as i128;
                    if (boundary - delta).rem_euclid(conductor) != 0 {
                        return Err(Error::plain(
                            "surface gauge cochain does not relate the two bulk cocycles",
                        ));
                    }
                }
            }
        }
        let mut result = Vec::with_capacity(self.0.len());
        for (marks, x) in &self.0 {
            model.spend(marks.len() + 1)?;
            if catalog.charge(marks)? != target.charge(marks)? {
                return Err(Error::plain(
                    "surface gauge changed a residual electric character",
                ));
            }
            let phase = chain::cycle(marks)
                .iter()
                .map(|(&[a, b], &s)| s * mu(a, b))
                .sum::<i128>();
            result.push((
                marks.clone(),
                rv_lens::shifted(model, from.conductor, x, phase)?,
            ));
        }
        Ok(Self(result))
    }
}

impl Double {
    pub(crate) fn surface_require_marking(
        &self,
        model: &mut Model,
        q: usize,
        marks: &[usize],
    ) -> Result<()> {
        let catalog = self.surface_catalog(model)?;
        catalog.check_charge(q)?;
        let actual = catalog.charge(marks)?;
        if actual != q {
            return Err(Error::plain(format!(
                "surface marking has gauge charge {actual}, expected {q}"
            )));
        }
        Ok(())
    }

    pub(crate) fn surface_inverse_charge(&self, model: &mut Model, q: usize) -> Result<usize> {
        let catalog = self.surface_catalog(model)?;
        catalog.check_charge(q)?;
        Ok(catalog.inverse(q))
    }

    pub(crate) fn surface_chain(
        &self,
        model: &mut Model,
        op: Move,
        marks: &Marking,
    ) -> Result<String> {
        check_genus(marks.len() / 2)?;
        let catalog = self.surface_catalog(model)?;
        let q = catalog.charge(marks)?;
        let (out, chain, phase) = op.apply(self, marks)?;
        if catalog.charge(&out)? != q {
            return Err(Error::plain("marking chain changes its electric sector"));
        }
        Ok(format!("rvDoubleSurfaceChain = output={out:?}, charge={q}, conductor={}, phase={phase}, chain={chain:?}, boundary_verified=true",self.conductor))
    }
}
