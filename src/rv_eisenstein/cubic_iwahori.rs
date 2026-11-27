use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::cubic_tree::{Tree, Vertex};

mod algebra;
pub(crate) use algebra::Word as GalleryWord;
use algebra::{Basis, Iwahori, Laurent, Spherical, Vector, Word};

const Q: i128 = 27;
const Q_MINUS_ONE: i128 = Q - 1;
const MAX_VECTOR_SUPPORT: usize = 8_192;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) enum Sort {
    Affine,
    Spherical,
    Laurent,
    Vector,
}

impl Sort {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Affine => "RVCubicIw",
            Self::Spherical => "RVCubicSph",
            Self::Laurent => "RVCubicSatake",
            Self::Vector => "RVCubicChamberVec",
        }
    }
}

pub(crate) enum Expression {
    Base(Label),
    Chamber(Box<[ParsedExpr; 2]>),
    Endpoint(usize, Box<ParsedExpr>),
    Panel(usize, usize, Box<ParsedExpr>),
    Omega(i64, Box<ParsedExpr>),
    Constant(Label, Constant),
    Unary(Unary, Box<ParsedExpr>),
    Multiply(Box<[ParsedExpr; 2]>),
    Act(Box<[ParsedExpr; 2]>),
}

pub(crate) enum Constant {
    Word(Word),
    Omega(i64),
    Bernstein,
    Projector,
    Spherical(usize),
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Delta,
    Satake,
    Compress,
    Scale(Rational),
}

pub(crate) enum Report {
    Gallery(Box<[ParsedExpr; 2]>),
    Word(Word),
    Kernel(Box<[ParsedExpr; 3]>),
    AffineCoefficient(i64, Word, Box<ParsedExpr>),
    SphericalCoefficient(usize, Box<ParsedExpr>),
    SatakeCoefficient(i64, Box<ParsedExpr>),
    VectorCoefficient(Box<[ParsedExpr; 2]>),
    Busemann(Box<[ParsedExpr; 2]>),
    Counts(Box<ParsedExpr>),
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Chamber([Vertex; 2]);

impl Chamber {
    pub(crate) fn new(tree: &Tree<'_>, a: Vertex, b: Vertex) -> Result<Self> {
        if tree.distance(&a, &b)? != 1 {
            return Err(Error::plain(
                "cubic Iwahori chamber requires adjacent lattice vertices",
            ));
        }
        Ok(Self::oriented(a, b))
    }

    fn oriented(a: Vertex, b: Vertex) -> Self {
        if a.level().is_multiple_of(2) {
            Self([a, b])
        } else {
            Self([b, a])
        }
    }

    pub(crate) fn base(model: &mut Model, tree: &Tree<'_>) -> Result<Self> {
        Ok(Self([tree.root(), tree.at(model, 1, 0)?]))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(97)?;
        Ok(Value::CubicChamber(Rc::new(self)))
    }

    pub(crate) fn description(&self) -> String {
        format!(
            "even=({}), odd=({})",
            self.0[0].description(),
            self.0[1].description()
        )
    }

    pub(crate) fn endpoint(&self, panel: usize) -> Result<Vertex> {
        self.0.get(panel).cloned().ok_or_else(panel_index)
    }

    fn alternatives(&self, model: &mut Model, tree: &Tree<'_>, panel: usize) -> Result<Vec<Self>> {
        let a = self.0.get(panel).ok_or_else(panel_index)?;
        let b = &self.0[1 - panel];
        let mut output = Vec::with_capacity(Q as usize);
        for next in tree.neighbors(model, a)? {
            if next != *b {
                output.push(Self::oriented(a.clone(), next));
            }
        }
        if output.len() != Q as usize {
            return Err(inconsistent());
        }
        Ok(output)
    }

    pub(crate) fn panel(
        &self,
        model: &mut Model,
        tree: &Tree<'_>,
        panel: usize,
        index: usize,
    ) -> Result<Self> {
        if index >= Q as usize {
            return Err(Error::plain("cubic Iwahori panel index must be in 0..27"));
        }
        self.alternatives(model, tree, panel)?
            .into_iter()
            .nth(index)
            .ok_or_else(inconsistent)
    }

    pub(crate) fn omega(&self, model: &mut Model, tree: &Tree<'_>, power: i64) -> Result<Self> {
        let a = tree.omega(model, power, &self.0[0])?;
        let b = tree.omega(model, power, &self.0[1])?;
        Self::new(tree, a, b)
    }

    fn relative(&self, tree: &Tree<'_>, target: &Self) -> Result<Word> {
        if self == target {
            return Ok(Word::ID);
        }
        let distance = |a: &Vertex| -> Result<usize> {
            Ok(tree
                .distance(a, &target.0[0])?
                .min(tree.distance(a, &target.0[1])?))
        };
        let a = distance(&self.0[0])?;
        let b = distance(&self.0[1])?;
        if a.abs_diff(b) != 1 {
            return Err(inconsistent());
        }
        Ok(Word {
            length: 1 + a.min(b),
            first: u8::from(b < a),
        })
    }

    pub(crate) fn gallery(&self, tree: &Tree<'_>, target: &Self) -> Result<String> {
        Ok(format!(
            "rvCubicGallery = {}",
            word_description(self.relative(tree, target)?)
        ))
    }
}

pub(crate) fn word_description(word: Word) -> String {
    match word.last() {
        Some(last) => format!(
            "length={}, left_descent={}, right_descent={last}",
            word.length, word.first
        ),
        None => "length=0, left_descent=none, right_descent=none".to_owned(),
    }
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) enum Object {
    Affine(Iwahori),
    Spherical(Spherical),
    Laurent(Laurent),
    Vector(Vector),
}

impl Object {
    pub(crate) fn zero(sort: Sort) -> Self {
        match sort {
            Sort::Affine => Self::Affine(Iwahori::zero()),
            Sort::Spherical => Self::Spherical(Spherical::zero()),
            Sort::Laurent => Self::Laurent(Laurent::zero()),
            Sort::Vector => Self::Vector(Vector::zero()),
        }
    }

    pub(crate) fn constant(model: &mut Model, c: Constant) -> Result<Self> {
        let result = match c {
            Constant::Word(word) => Self::Affine(Iwahori::basis(Basis { word, omega: 0 })),
            Constant::Omega(omega) => Self::Affine(Iwahori::basis(Basis {
                word: Word::ID,
                omega,
            })),
            Constant::Bernstein => Self::Affine(algebra::bernstein(model)?),
            Constant::Projector => Self::Affine(algebra::compression(model, &Spherical::basis(0))?),
            Constant::Spherical(n) => Self::Spherical(Spherical::basis(n)),
        };
        Ok(result)
    }

    pub(crate) fn delta(chamber: Chamber) -> Self {
        Self::Vector(Vector::basis(chamber))
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        let cells = match &self {
            Self::Affine(a) => 1 + 9 * a.entries().len(),
            Self::Spherical(a) => 1 + 6 * a.entries().len(),
            Self::Laurent(a) => 1 + 6 * a.entries().len(),
            Self::Vector(a) => 1 + 104 * a.entries().len(),
        };
        model.allocate(cells)?;
        Ok(Value::CubicHecke(Rc::new(self)))
    }

    pub(crate) fn description(&self) -> String {
        match self {
            Self::Affine(a) => {
                let rows: Vec<_> = a
                    .entries()
                    .iter()
                    .map(|(k, v)| (k.omega, k.word.first, k.word.length, fraction(*v)))
                    .collect();
                format!("terms={rows:?}")
            }
            Self::Spherical(a) => {
                let rows: Vec<_> = a
                    .entries()
                    .iter()
                    .map(|(k, v)| (*k, fraction(*v)))
                    .collect();
                format!("radii={rows:?}")
            }
            Self::Laurent(a) => {
                let rows: Vec<_> = a
                    .entries()
                    .iter()
                    .map(|(k, v)| (*k, fraction(*v)))
                    .collect();
                format!("terms={rows:?}")
            }
            Self::Vector(a) => format!("support={}", a.entries().len()),
        }
    }

    pub(crate) fn add(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        Ok(match (self, rhs) {
            (Self::Affine(a), Self::Affine(b)) => Self::Affine(a.add(model, b)?),
            (Self::Spherical(a), Self::Spherical(b)) => Self::Spherical(a.add(model, b)?),
            (Self::Laurent(a), Self::Laurent(b)) => Self::Laurent(a.add(model, b)?),
            (Self::Vector(a), Self::Vector(b)) => Self::Vector(a.add(model, b)?),
            _ => unreachable!(),
        })
    }

    pub(crate) fn scale(&self, model: &mut Model, c: Rational) -> Result<Self> {
        Ok(match self {
            Self::Affine(a) => Self::Affine(a.scale(model, c)?),
            Self::Spherical(a) => Self::Spherical(a.scale(model, c)?),
            Self::Laurent(a) => Self::Laurent(a.scale(model, c)?),
            Self::Vector(a) => Self::Vector(a.scale(model, c)?),
        })
    }

    pub(crate) fn multiply(&self, model: &mut Model, rhs: &Self) -> Result<Self> {
        Ok(match (self, rhs) {
            (Self::Affine(a), Self::Affine(b)) => Self::Affine(algebra::combine(model, a, b)?),
            (Self::Spherical(a), Self::Spherical(b)) => {
                Self::Spherical(algebra::radial_product(model, a, b)?)
            }
            (Self::Laurent(a), Self::Laurent(b)) => {
                Self::Laurent(algebra::laurent_product(model, a, b)?)
            }
            _ => unreachable!(),
        })
    }

    pub(crate) fn satake(&self, model: &mut Model) -> Result<Self> {
        let Self::Spherical(a) = self else {
            unreachable!()
        };
        Ok(Self::Laurent(algebra::satake(model, a)?))
    }

    pub(crate) fn compress(&self, model: &mut Model) -> Result<Self> {
        let Self::Spherical(a) = self else {
            unreachable!()
        };
        Ok(Self::Affine(algebra::compression(model, a)?))
    }

    pub(crate) fn coefficients(&self, query: Coefficient<'_>) -> String {
        let (name, coefficient) = match (self, query) {
            (Self::Affine(a), Coefficient::Affine(omega, word)) => {
                ("rvCubicIwCoeff", a.coefficient(&Basis { omega, word }))
            }
            (Self::Spherical(a), Coefficient::Spherical(n)) => {
                ("rvCubicSphCoeff", a.coefficient(&n))
            }
            (Self::Laurent(a), Coefficient::Laurent(k)) => {
                ("rvCubicSatakeCoeff", a.coefficient(&k))
            }
            (Self::Vector(a), Coefficient::Vector(v)) => ("rvCubicIwVectorCoeff", a.coefficient(v)),
            _ => unreachable!(),
        };
        format!("{name} = {}", fraction(coefficient))
    }

    pub(crate) fn kernel(
        &self,
        model: &mut Model,
        tree: &Tree<'_>,
        source: &Chamber,
        target: &Chamber,
    ) -> Result<String> {
        let Self::Affine(a) = self else {
            unreachable!()
        };
        let mut total = Rational::ZERO;
        let mut targets = BTreeMap::new();
        for (&basis, &weight) in a.entries() {
            model.spend(1)?;
            let power = -basis.omega.rem_euclid(6);
            let word = match targets.entry(power) {
                std::collections::btree_map::Entry::Vacant(slot) => {
                    let preimage = target.omega(model, tree, power)?;
                    *slot.insert(source.relative(tree, &preimage)?)
                }
                std::collections::btree_map::Entry::Occupied(slot) => *slot.get(),
            };
            if word == basis.word {
                total = total.add(weight).map_err(|_| arithmetic())?;
            }
        }
        Ok(format!("rvCubicIwKernel = {}", fraction(total)))
    }

    pub(crate) fn action(&self, model: &mut Model, tree: &Tree<'_>, vector: &Self) -> Result<Self> {
        let (Self::Affine(operator), Self::Vector(source)) = (self, vector) else {
            unreachable!()
        };
        let mut output = Vector::zero();
        for (&basis, &weight) in operator.entries() {
            let mut current = source.clone();
            for i in 0..basis.word.length {
                model.spend(1)?;
                if current.entries().is_empty() {
                    break;
                }
                let panel = usize::from(basis.word.panel(i));
                let mut next = Vector::zero();
                for (c, &w) in current.entries() {
                    for target in c.alternatives(model, tree, panel)? {
                        next.insert(model, target, w)?;
                    }
                    require_vector_capacity(&next)?;
                }
                current = next;
            }
            if basis.omega.rem_euclid(6) == 0 {
                output.add_scaled(model, &current, weight)?;
                require_vector_capacity(&output)?;
            } else {
                for (c, &w) in current.entries() {
                    let target = c.omega(model, tree, basis.omega)?;
                    output.insert(model, target, w.multiply(weight).map_err(|_| arithmetic())?)?;
                    require_vector_capacity(&output)?;
                }
            }
        }
        Ok(Self::Vector(output))
    }
}

pub(crate) enum Coefficient<'a> {
    Affine(i64, Word),
    Spherical(usize),
    Laurent(i64),
    Vector(&'a Chamber),
}

pub(crate) fn busemann(tree: &Tree<'_>, x: &Vertex, ray: &Vertex) -> Result<i64> {
    if x.level() > ray.level() {
        return Err(Error::plain(
            "Busemann ray requires at least the target vertex radius",
        ));
    }
    let meet = i64::try_from(tree.lca(x, ray)?.level()).map_err(|_| arithmetic())?;
    let level = i64::try_from(x.level()).map_err(|_| arithmetic())?;
    Ok(2 * meet - level)
}

pub(crate) fn satake_counts(model: &mut Model, tree: &Tree<'_>, ray: &Vertex) -> Result<String> {
    let radius = ray.level();
    if radius > 3 {
        return Err(Error::plain(
            "explicit Busemann sphere enumeration supports radius at most 3",
        ));
    }
    let count = if radius == 0 {
        1
    } else {
        (Q as usize + 1) * (Q as usize).pow(radius as u32 - 1)
    };
    let mut heights = Laurent::zero();
    for i in 0..count {
        let x = tree.at(model, radius, i)?;
        let k = busemann(tree, &x, ray)?;
        heights.insert(model, k, Rational::integer(1))?;
    }
    if heights != algebra::satake(model, &Spherical::basis(radius))? {
        return Err(inconsistent());
    }
    let rows: Vec<_> = heights
        .entries()
        .iter()
        .map(|(&k, c)| (k, c.numerator))
        .collect();
    Ok(format!(
        "rvCubicSatakeCounts = radius={radius}, heights={rows:?}"
    ))
}

fn require_vector_capacity(vector: &Vector) -> Result<()> {
    if vector.entries().len() > MAX_VECTOR_SUPPORT {
        return Err(Error::plain("explicit cubic Iwahori vector exceeds 8192 chambers; use the endpoint kernel for individual gallery coefficients"));
    }
    Ok(())
}

pub(crate) fn fraction(x: Rational) -> String {
    if x.denominator == 1 {
        x.numerator.to_string()
    } else {
        format!("{}/{}", x.numerator, x.denominator)
    }
}

pub(crate) fn arithmetic() -> Error {
    Error::plain("cubic Hecke coefficients exceed the signed 128-bit exact arithmetic range")
}

fn panel_index() -> Error {
    Error::plain("cubic Iwahori panel must be 0 or 1")
}

fn inconsistent() -> Error {
    Error::plain("inconsistent cubic Iwahori lattice, gallery or Satake identity")
}
