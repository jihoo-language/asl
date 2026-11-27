use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::sheaf::{self, Matrix, Sheaf};
use crate::sheaf_hyper::{Complex, Map as SheafMap, NamedSheaf};
use crate::sheaf_leray;
use crate::syntax::Label;

pub(crate) mod braid;
pub(crate) mod convolution;
mod hom;
mod linear;
pub(crate) mod parabolic;
mod realize;

pub(crate) use hom::{Hom, Pair};
pub(crate) use realize::Realization;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Extension {
    Shriek,
    Star,
    Middle,
}

pub(crate) enum Origin {
    Declared,
    Extended(Rc<str>, Extension),
    Dual(Rc<str>),
}

pub(crate) struct Disk {
    pub(crate) prime: Prime,
    nearby: usize,
    vanishing: usize,
    u: Rc<Matrix>,
    v: Rc<Matrix>,
    t: Rc<Matrix>,
    s: Rc<Matrix>,
    t_inverse: Matrix,
    pub(crate) origin: Origin,
}

impl Disk {
    pub(crate) fn declared(
        model: &mut Model,
        characteristic: usize,
        nearby: usize,
        vanishing: usize,
        u: Vec<i64>,
        v: Vec<i64>,
    ) -> Result<Self> {
        let p = sheaf::field(model, characteristic)?;
        sheaf::check_cells(nearby.saturating_add(vanishing))?;
        let u = linear::literal(model, p, vanishing, nearby, u)?;
        let v = linear::literal(model, p, nearby, vanishing, v)?;
        Self::new(model, p, nearby, vanishing, u, v, Origin::Declared)
    }

    fn new(
        model: &mut Model,
        p: Prime,
        nearby: usize,
        vanishing: usize,
        u: Matrix,
        v: Matrix,
        origin: Origin,
    ) -> Result<Self> {
        let vu = Matrix::compose(model, p, &v, &u)?;
        let uv = Matrix::compose(model, p, &u, &v)?;
        let t = linear::one_minus(model, p, nearby, &vu)?;
        let s = linear::one_minus(model, p, vanishing, &uv)?;
        let t_inverse = linear::inverse(model, p, nearby, &t)?;
        model.allocate(u.cells() + v.cells() + t.cells() + s.cells() + t_inverse.cells() + 6)?;
        Ok(Self {
            prime: p,
            nearby,
            vanishing,
            u: Rc::new(u),
            v: Rc::new(v),
            t: Rc::new(t),
            s: Rc::new(s),
            t_inverse,
            origin,
        })
    }

    fn extension(&self) -> Option<Extension> {
        match self.origin {
            Origin::Extended(_, extension) => Some(extension),
            _ => None,
        }
    }

    pub(crate) fn extend(
        &self,
        model: &mut Model,
        source: Rc<str>,
        extension: Extension,
    ) -> Result<Self> {
        let p = self.prime;
        let n = self.nearby;
        let k = linear::one_minus(model, p, n, &self.t)?;
        let (u, v, m) = match extension {
            Extension::Shriek => (Matrix::identity(model, n)?, k, n),
            Extension::Star => (k, Matrix::identity(model, n)?, n),
            Extension::Middle => linear::image_factorization(model, p, n, &k)?,
        };
        let result = Self::new(model, p, n, m, u, v, Origin::Extended(source, extension))?;
        debug_assert!(*result.t == *self.t);
        Ok(result)
    }

    pub(crate) fn dual(&self, model: &mut Model, source: Rc<str>) -> Result<Self> {
        let p = self.prime;
        let u = linear::transpose(model, &self.v)?;
        let u_transpose = linear::transpose(model, &self.u)?;
        let dual_nearby = linear::transpose(model, &self.t_inverse)?;
        let product = Matrix::compose(model, p, &dual_nearby, &u_transpose)?;
        let v = linear::negate(model, p, &product)?;
        let result = Self::new(
            model,
            p,
            self.nearby,
            self.vanishing,
            u,
            v,
            Origin::Dual(source),
        )?;
        debug_assert!(*result.t == dual_nearby);
        Ok(result)
    }

    pub(crate) fn canonical(&self, model: &mut Model, target: &Self) -> Result<Pair> {
        if self.prime != target.prime || self.nearby != target.nearby || self.t != target.t {
            return Err(Error::plain(
                "canonical disk map requires the same based nearby monodromy",
            ));
        }
        let b = if self.extension() == Some(Extension::Shriek) {
            target.u.clone()
        } else if target.extension() == Some(Extension::Star) {
            self.v.clone()
        } else {
            return Err(Error::plain(
                "canonical disk map requires a shriek source or a star target",
            ));
        };
        Ok(Pair {
            a: Rc::new(Matrix::identity(model, self.nearby)?),
            b,
        })
    }

    pub(crate) fn bidual(&self, model: &mut Model, inverse: bool) -> Result<Pair> {
        let t = if inverse { &self.t_inverse } else { &self.t };
        Ok(Pair {
            a: Rc::new(linear::negate(model, self.prime, t)?),
            b: Rc::new(Matrix::identity(model, self.vanishing)?),
        })
    }

    pub(crate) fn extended_map(
        &self,
        model: &mut Model,
        target: &Self,
        map: Pair,
        extension: Extension,
    ) -> Result<Pair> {
        let b = if extension == Extension::Middle {
            let composition = Matrix::compose(model, self.prime, &map.a, &self.v)?;
            Rc::new(linear::factor_through(
                model,
                self.prime,
                &target.v,
                self.vanishing,
                &composition,
            )?)
        } else {
            map.a.clone()
        };
        Ok(Pair { a: map.a, b })
    }

    pub(crate) fn complex(
        &self,
        model: &mut Model,
        name: &Label,
        costalk: bool,
    ) -> Result<Complex> {
        let (ranks, map_degree, d) = if costalk {
            (vec![0, self.vanishing, self.nearby], 1, self.v.clone())
        } else {
            (vec![self.nearby, self.vanishing], 0, self.u.clone())
        };
        let one = model.named("One").expect("builtin singleton");
        let mut components = Vec::new();
        for (degree, rank) in ranks.into_iter().enumerate() {
            components.push(NamedSheaf {
                name: sheaf_leray::component_name(model, name, degree)?,
                sheaf: Rc::new(constant_on_point(model, one, self.prime, rank)?),
            });
        }
        let mut maps = Vec::new();
        for (degree, pair) in components.windows(2).enumerate() {
            let matrix = if degree == map_degree {
                d.clone()
            } else {
                Rc::new(Matrix::default())
            };
            let edge =
                SheafMap::from_matrices(model, pair[0].clone(), pair[1].clone(), vec![matrix])?;
            let label = Label {
                text: format!("{} differential at {degree}", name.text),
                span: name.span,
            };
            maps.push((label, Rc::new(edge)));
        }
        Complex::build(model, components, maps)
    }

    pub(crate) fn dimensions(&self) -> (usize, usize) {
        (self.nearby, self.vanishing)
    }

    pub(crate) fn cohomology_ranks(&self, model: &mut Model) -> Result<([usize; 2], [usize; 2])> {
        let u = linear::rank(model, self.prime, &self.u)?;
        let v = linear::rank(model, self.prime, &self.v)?;
        Ok((
            [self.nearby - u, self.vanishing - u],
            [self.vanishing - v, self.nearby - v],
        ))
    }

    pub(crate) fn report(&self, model: &mut Model, operator: Operator) -> Result<String> {
        let (matrix, h, w) = self.operator(operator);
        linear::report(model, matrix, h, w)
    }

    pub(crate) fn native_operator(
        &self,
        model: &mut Model,
        source: NamedSheaf,
        target: NamedSheaf,
        operator: Operator,
    ) -> Result<SheafMap> {
        let (matrix, h, w) = self.operator(operator);
        let one = model.named("One").expect("builtin singleton");
        if source.sheaf.base() != one
            || target.sheaf.base() != one
            || source.sheaf.prime != self.prime
            || target.sheaf.prime != self.prime
            || source.sheaf.dims() != [w]
            || target.sheaf.dims() != [h]
        {
            return Err(Error::plain(
                "disk operator needs matching vector sheaves on One over its prime field",
            ));
        }
        SheafMap::from_matrices(model, source, target, vec![matrix.clone()])
    }

    fn operator(&self, op: Operator) -> (&Rc<Matrix>, usize, usize) {
        match op {
            Operator::U => (&self.u, self.vanishing, self.nearby),
            Operator::V => (&self.v, self.nearby, self.vanishing),
            Operator::Nearby => (&self.t, self.nearby, self.nearby),
            Operator::Vanishing => (&self.s, self.vanishing, self.vanishing),
        }
    }
}

fn constant_on_point(
    model: &mut Model,
    base: crate::model::TypeId,
    prime: Prime,
    rank: usize,
) -> Result<Sheaf> {
    Sheaf::from_cover_maps(model, base, prime, vec![rank], |_, _, _| {
        unreachable!("singleton")
    })
}

#[derive(Clone, Copy)]
pub(crate) enum Operator {
    U,
    V,
    Nearby,
    Vanishing,
}

impl Operator {
    pub(crate) fn from_label(label: &Label) -> Result<Self> {
        match label.text.as_str() {
            "u" => Ok(Self::U),
            "v" => Ok(Self::V),
            "T" => Ok(Self::Nearby),
            "S" => Ok(Self::Vanishing),
            _ => Err(Error::new(
                label.span,
                "expected a disk operator u, v, T or S",
            )),
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::U => "u",
            Self::V => "v",
            Self::Nearby => "T",
            Self::Vanishing => "S",
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Bridge {
    Canonical,
    Bidual,
    BidualInverse,
    DualMap,
    ExtensionMap,
}

impl Bridge {
    pub(crate) fn detect(name: &str) -> Option<Self> {
        Some(match name {
            "shDiskCanonical" => Self::Canonical,
            "shDiskBidual" => Self::Bidual,
            "shDiskBidualInv" => Self::BidualInverse,
            "shDiskDualMap" => Self::DualMap,
            "shDiskExtMap" => Self::ExtensionMap,
            _ => return None,
        })
    }

    pub(crate) fn has_input(self) -> bool {
        matches!(self, Self::DualMap | Self::ExtensionMap)
    }

    pub(crate) fn constant(
        self,
        model: &mut Model,
        from: Rc<str>,
        to: Rc<str>,
    ) -> Result<(crate::model::TypeId, crate::model::Value)> {
        let a = model.disk_named(&from)?;
        let b = model.disk_named(&to)?;
        let (from, to, raw) = match self {
            Self::Canonical => {
                let raw = a.canonical(model, &b)?;
                (from, to, raw)
            }
            Self::Bidual | Self::BidualInverse => {
                let Origin::Dual(middle) = &b.origin else {
                    return Err(Error::plain(
                        "disk biduality requires the declared double dual",
                    ));
                };
                let middle = model.disk_named(middle)?;
                if !matches!(&middle.origin,Origin::Dual(name) if name==&from) {
                    return Err(Error::plain(
                        "disk biduality requires the declared double dual",
                    ));
                }
                let inverse = matches!(self, Self::BidualInverse);
                let raw = a.bidual(model, inverse)?;
                if inverse {
                    (to, from, raw)
                } else {
                    (from, to, raw)
                }
            }
            _ => unreachable!("map variable"),
        };
        let (ty, hom) = model.disk_hom_type(from, to)?;
        Ok((ty, hom.encode(model, &raw)?))
    }

    pub(crate) fn map(
        self,
        model: &mut Model,
        from: Rc<str>,
        to: Rc<str>,
        (source, target): (Rc<str>, Rc<str>),
    ) -> Result<(crate::model::TypeId, Unary)> {
        let a = model.disk_named(&from)?;
        let b = model.disk_named(&to)?;
        let action = match self {
            Self::DualMap => {
                if !matches!(&a.origin,Origin::Dual(original) if original==&target)
                    || !matches!(&b.origin,Origin::Dual(original) if original==&source)
                {
                    return Err(Error::plain(
                        "disk dual map requires the target dual followed by the source dual",
                    ));
                }
                UnaryKind::Dual
            }
            Self::ExtensionMap => match (&a.origin, &b.origin) {
                (Origin::Extended(left, f), Origin::Extended(right, g))
                    if left == &source && right == &target && f == g =>
                {
                    UnaryKind::Extend(*f)
                }
                _ => return Err(Error::plain(
                    "disk extension map requires the same kind of extension of each map endpoint",
                )),
            },
            _ => unreachable!("constant"),
        };
        let (_, input) = model.disk_hom_type(source, target)?;
        let (ty, output) = model.disk_hom_type(from, to)?;
        Ok((
            ty,
            Unary {
                input,
                output,
                action,
            },
        ))
    }
}

enum UnaryKind {
    Dual,
    Extend(Extension),
}

pub(crate) struct Unary {
    input: Rc<Hom>,
    output: Rc<Hom>,
    action: UnaryKind,
}

impl Unary {
    pub(crate) fn apply(
        &self,
        model: &mut Model,
        coordinates: &[crate::prime_linear::Coordinate],
    ) -> Result<crate::model::Value> {
        let input = self.input.decode(model, coordinates)?;
        let raw = match self.action {
            UnaryKind::Dual => input.dual(model)?,
            UnaryKind::Extend(extension) => {
                self.output
                    .source
                    .extended_map(model, &self.output.target, input, extension)?
            }
        };
        self.output.encode(model, &raw)
    }
}

pub(crate) struct Composition {
    pub(crate) outer: Rc<Hom>,
    pub(crate) inner: Rc<Hom>,
    pub(crate) output: Rc<Hom>,
}

impl Composition {
    pub(crate) fn apply(
        &self,
        model: &mut Model,
        outer: &[crate::prime_linear::Coordinate],
        inner: &[crate::prime_linear::Coordinate],
    ) -> Result<crate::model::Value> {
        let after = self.outer.decode(model, outer)?;
        let before = self.inner.decode(model, inner)?;
        let raw = after.compose(model, self.output.source.prime, &before)?;
        self.output.encode(model, &raw)
    }
}

pub(crate) enum Command {
    Declare(Label, [usize; 3], [Vec<i64>; 2]),
    Extend(Label, Label, Extension),
    Dual(Label, Label),
    Map(Label, [Label; 2], [Vec<i64>; 2]),
    Matrix(Label, Operator),
    MapMatrix(Box<crate::syntax::ParsedExpr>),
    Ranks(Label),
    HomDimension([Label; 2]),
    Complex(Label, Label, bool),
    Operator(Label, Label, Operator, [Label; 2]),
    Realize(Label, Box<[Label; 7]>),
}

pub(crate) enum Expression {
    Id(Label),
    Generator([Label; 2], usize),
    Compose(Box<[crate::syntax::ParsedExpr; 2]>),
    Bridge(Bridge, [Label; 2], Option<Box<crate::syntax::ParsedExpr>>),
    Lift(Label, Box<crate::syntax::ParsedExpr>),
}
