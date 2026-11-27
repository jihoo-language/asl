use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::Model;
use crate::syntax::{Label, ParsedExpr};

use super::cubic_tree::{Operator, Tree, Vertex};

mod boundary;
mod knapp;

pub(crate) use boundary::{mass, Function};
pub(crate) use knapp::Parameter as Knapp;

use boundary::{complementary, Node};

pub(crate) type Rat = Rational;
pub(super) const Q: i128 = 27;

pub(crate) enum Expression {
    Cylinder(Box<ParsedExpr>),
    Scalar(Label, Rat),
    Unary(Unary, Box<ParsedExpr>),
    Compact(Rat, Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Refine(usize),
    Coarsen(usize),
    Average(usize),
    Difference(usize),
    Scale(Rat),
    Knapp(Rat),
    Move(Move, Rat),
}

#[derive(Clone, Copy)]
pub(crate) enum Move {
    Translate,
    TranslateBack,
    Omega(i64),
}

pub(crate) enum Report {
    Mass(Box<ParsedExpr>),
    Boundary(BoundaryQuery, Box<[ParsedExpr; 2]>),
    CompactCovariance(Rat, Box<[ParsedExpr; 3]>),
    KnappSpectrum(Rat, usize),
    KnappKernel(Rat, Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Copy)]
pub(crate) enum BoundaryQuery {
    Value,
    Integral,
    Poisson(Rat),
    Eigen(Rat),
    Covariance(Move, Rat),
}

#[derive(Clone, Copy)]
pub(crate) enum Motion<'a> {
    Compact(&'a Operator),
    Move(Move),
}

impl Motion<'_> {
    fn image(self, model: &mut Model, tree: &Tree<'_>, v: &Vertex) -> Result<Vertex> {
        match self {
            Self::Compact(operator) => tree.apply(model, operator, v),
            Self::Move(Move::Translate) => tree.translate(model, v),
            Self::Move(Move::TranslateBack) => tree.translate_back(model, v),
            Self::Move(Move::Omega(power)) => tree.omega(model, power, v),
        }
    }

    fn preimage(self, model: &mut Model, tree: &Tree<'_>, v: &Vertex) -> Result<Vertex> {
        match self {
            Self::Compact(operator) => {
                let inverse = tree.inverse(model, operator)?;
                tree.apply(model, &inverse, v)
            }
            Self::Move(Move::Translate) => tree.translate_back(model, v),
            Self::Move(Move::TranslateBack) => tree.translate(model, v),
            Self::Move(Move::Omega(power)) => tree.omega(model, -power.rem_euclid(6), v),
        }
    }

    pub(crate) fn act(
        self,
        model: &mut Model,
        tree: &Tree<'_>,
        z: Rat,
        f: &Function,
    ) -> Result<Function> {
        let mut transported = Node::Constant(Rat::ZERO);
        f.visit(&mut |key, c| {
            if key.is_empty() {
                return transported.add_indicator(model, key, 0, c);
            }
            let endpoint = tree.vertex_at_address(model, key)?;
            let parent = tree.ancestor(&endpoint, key.len() - 1)?;
            let target = self.image(model, tree, &endpoint)?;
            let nearer = self.image(model, tree, &parent)?;
            if tree.distance(&target, &nearer)? != 1 {
                return Err(Error::plain(
                    "principal-series motion failed to preserve a cylinder edge",
                ));
            }
            let complement = target.level() < nearer.level();
            let bound = if complement { &nearer } else { &target };
            let key = tree.address(model, bound)?;
            complementary(&mut transported, model, complement, &key, c)
        })?;
        let transported = Function::from_node(transported, f.presentation());
        let root = self.image(model, tree, &tree.root())?;
        let path = tree.address(model, &root)?;
        if path.is_empty() {
            return Ok(transported);
        }
        let depth = i64::try_from(path.len()).map_err(|_| arithmetic())?;
        let mut previous = power(z, -depth)?;
        let mut weight = Node::Constant(previous);
        for k in 1..=path.len() {
            let j = i64::try_from(k).map_err(|_| arithmetic())?;
            let value = power(z, 2 * j - depth)?;
            weight.add_indicator(model, &path[..k], 0, subtract(value, previous)?)?;
            previous = value;
        }
        transported.product(model, &Function::from_node(weight, path.len()))
    }

    pub(crate) fn covariance(
        self,
        model: &mut Model,
        tree: &Tree<'_>,
        t: Rat,
        f: &Function,
        v: &Vertex,
    ) -> Result<String> {
        let z = divide(q(), t)?;
        let transported = self.act(model, tree, z, f)?;
        let lhs = poisson(model, tree, &transported, v, t)?;
        let preimage = self.preimage(model, tree, v)?;
        let rhs = poisson(model, tree, f, &preimage, t)?;
        if lhs != rhs {
            return Err(Error::plain(
                "Poisson transform failed principal-series covariance",
            ));
        }
        Ok(format!(
            "rvCubicPoissonCovariance = z={}, value={}, inverse_radius={}",
            fraction(z),
            fraction(lhs),
            preimage.level()
        ))
    }
}

pub(crate) fn poisson(
    model: &mut Model,
    tree: &Tree<'_>,
    f: &Function,
    v: &Vertex,
    t: Rat,
) -> Result<Rat> {
    let path = tree.address(model, v)?;
    f.poisson(model, &path, t)
}

pub(crate) fn poisson_eigen(
    model: &mut Model,
    tree: &Tree<'_>,
    f: &Function,
    v: &Vertex,
    t: Rat,
) -> Result<String> {
    let eigenvalue = add(t, divide(q(), t)?)?;
    let center = poisson(model, tree, f, v, t)?;
    let mut adjacent = Rat::ZERO;
    for neighbor in tree.neighbors(model, v)? {
        adjacent = add(adjacent, poisson(model, tree, f, &neighbor, t)?)?;
    }
    if adjacent != multiply(eigenvalue, center)? {
        return Err(Error::plain(
            "Poisson transform failed the tree adjacency eigenvalue identity",
        ));
    }
    Ok(format!(
        "rvCubicPoissonEigen = eigenvalue={}, value={}, adjacency={}",
        fraction(eigenvalue),
        fraction(center),
        fraction(adjacent)
    ))
}

pub(crate) fn rational(a: i128, b: i128) -> Result<Rat> {
    if b <= 0 {
        return Err(Error::plain(
            "boundary rational denominator must be positive",
        ));
    }
    Rat::ratio(a, b).map_err(|_| arithmetic())
}

pub(crate) fn spectral(a: i128, b: i128) -> Result<Rat> {
    if a == 0 {
        return Err(Error::plain(
            "principal-series spectral parameter must be nonzero",
        ));
    }
    rational(a, b)
}

pub(crate) fn fraction(x: Rat) -> String {
    super::cubic_iwahori::fraction(x)
}

pub(super) const fn integer(x: i128) -> Rat {
    Rat::integer(x)
}
pub(super) const fn q() -> Rat {
    integer(Q)
}

pub(super) fn add(a: Rat, b: Rat) -> Result<Rat> {
    a.add(b).map_err(|_| arithmetic())
}

pub(super) fn neg(a: Rat) -> Result<Rat> {
    Ok(Rat {
        numerator: a.numerator.checked_neg().ok_or_else(arithmetic)?,
        denominator: a.denominator,
    })
}

pub(super) fn subtract(a: Rat, b: Rat) -> Result<Rat> {
    add(a, neg(b)?)
}

pub(super) fn multiply(a: Rat, b: Rat) -> Result<Rat> {
    a.multiply(b).map_err(|_| arithmetic())
}

pub(super) fn divide(a: Rat, b: Rat) -> Result<Rat> {
    if b.numerator == 0 {
        return Err(Error::plain(
            "principal-series rational denominator vanishes",
        ));
    }
    multiply(
        a,
        Rat::ratio(b.denominator, b.numerator).map_err(|_| arithmetic())?,
    )
}

pub(super) fn power(mut base: Rat, exponent: i64) -> Result<Rat> {
    if exponent < 0 {
        base = divide(integer(1), base)?;
    }
    let mut exponent = exponent.unsigned_abs();
    let mut result = integer(1);
    while exponent != 0 {
        if !exponent.is_multiple_of(2) {
            result = multiply(result, base)?;
        }
        exponent /= 2;
        if exponent != 0 {
            base = multiply(base, base)?;
        }
    }
    Ok(result)
}

fn arithmetic() -> Error {
    Error::plain("principal-series rational arithmetic exceeds the checked signed 128-bit range")
}
