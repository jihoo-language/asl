use std::rc::Rc;

use super::composition::{self, Coordinates};
use super::Cycle;
use super::Hom;
use crate::ec_morphism::Morphism;
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_qp_quad::{self, Element};
use crate::syntax::{Label, ParsedExpr};

pub(crate) mod period;
pub(crate) mod theta;

pub(crate) enum Statement {
    Period(Box<period::Statement>),
    Theta(Box<theta::Statement>),
    Define(Label, [Label; 4]),
    Report(Report, Box<ParsedExpr>),
    Pair(Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Copy)]
pub(crate) enum Report {
    Data,
    Determinant,
    Map,
}

pub(crate) enum Expression {
    Period(Box<period::Expression>),
    Theta(Box<theta::Expression>),
    Divisor(Label, [i128; 2], Box<ParsedExpr>),
    Mukai([i128; 2], Box<ParsedExpr>),
    Map(Label, Box<[ParsedExpr; 4]>),
    FromGraph(Label, Box<ParsedExpr>),
    Scale([i128; 2], Box<ParsedExpr>),
    Block(usize, Box<ParsedExpr>),
    Unary(Unary, Box<ParsedExpr>),
    Binary(Binary, Box<[ParsedExpr; 2]>),
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Fourier,
    Line,
    FirstChern,
    Graph,
    Dagger,
    Adjugate,
}

#[derive(Clone, Copy)]
pub(crate) enum Binary {
    Tensor,
    Pull,
    Push,
    Compose,
}

type Matrix = [[Element; 2]; 2];

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Divisor {
    a: Element,
    d: Element,
    z: Element,
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Mukai {
    rank: Element,
    divisor: Divisor,
    last: Element,
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Endomorphism {
    matrix: Matrix,
    determinant: Element,
    degree: i128,
}

pub(crate) struct Abelian {
    pub(crate) hom_name: Label,
    frames: [Rc<Hom>; 4],
    bases: [[Element; 2]; 4],
    field: Quadratic,
    weight: i128,
}

#[derive(Clone, Copy)]
struct Quadratic(i128);

impl Abelian {
    pub(crate) fn new(
        model: &mut Model,
        hom_name: Label,
        [forward, reverse, end_e, end_f]: [Rc<Hom>; 4],
    ) -> Result<Self> {
        let [f, g] = &forward.maps;
        let [fd, gd] = &reverse.maps;
        let (e, ff) = (&f.source, &f.target);
        if [&forward, &reverse, &end_e, &end_f]
            .into_iter()
            .any(|frame| {
                !frame.maximal || frame.index != 1 || frame.fundamental != forward.fundamental
            })
        {
            return Err(Error::plain("CM abelian marking requires four full certified maximal-order Hom frames in a common field"));
        }
        if e == ff
            || fd.source != *ff
            || fd.target != *e
            || end_e.maps[0].source != *e
            || end_e.maps[0].target != *e
            || end_f.maps[0].source != *ff
            || end_f.maps[0].target != *ff
        {
            return Err(Error::plain("CM abelian marking needs forward, reverse and two endomorphism frames on the named factors"));
        }
        let [domain, range] = forward.endpoints.each_ref().map(Rc::as_ref);
        composition::verify_dual(model, f, fd, [domain, range])?;
        composition::verify_dual(model, g, gd, [domain, range])?;
        let alpha_e = fd.compose(model, domain, g)?;
        let alpha_f = g.compose(model, range, fd)?;
        let first = f.compose(model, domain, &alpha_e)?;
        let second = alpha_f.compose(model, domain, f)?;
        if !first.compare(model, &second)?.is_eq() {
            return Err(inconsistent());
        }
        let [a, b, _] = forward.form;
        let disc = forward.fundamental;
        let (field, alpha) = if disc % 4 == 0 {
            (Quadratic(disc / 4), Element::new(b, 2, 2)?)
        } else {
            (Quadratic(disc), Element::new(b, 1, 2)?)
        };
        let q = field.scale(model, alpha, 1, a)?;
        let reverse_alpha = field.conjugate(model, alpha)?;
        let e_basis = endomorphism_basis(model, field, &end_e, alpha_e, alpha)?;
        let f_basis = endomorphism_basis(model, field, &end_f, alpha_f, alpha)?;
        let frames = [end_e, reverse, forward, end_f];
        let bases = [
            e_basis,
            [integer(a), reverse_alpha],
            [Element::ONE, q],
            f_basis,
        ];
        let atlas = Self {
            hom_name,
            frames,
            bases,
            field,
            weight: a,
        };
        for (i, basis) in atlas.bases.iter().enumerate() {
            let frame = &atlas.frames[i];
            let [u, v] = *basis;
            let sum = field.add(model, u, v)?;
            if atlas.block_degree(model, i, u)? != integer(frame.form[0])
                || atlas.block_degree(model, i, v)? != integer(frame.form[2])
                || atlas.block_degree(model, i, sum)? != integer(frame.plus as i128)
            {
                return Err(inconsistent());
            }
        }
        model.allocate(44 + atlas.hom_name.text.len())?;
        Ok(atlas)
    }

    pub(crate) fn frame(&self, index: usize) -> Result<&Hom> {
        self.frames
            .get(index)
            .map(Rc::as_ref)
            .ok_or_else(|| Error::plain("CM abelian block index must be between zero and three"))
    }

    pub(crate) fn divisor(
        &self,
        model: &mut Model,
        [a, d]: [i128; 2],
        coordinate: &Coordinates,
    ) -> Result<Divisor> {
        let z = self.embed(model, 2, coordinate)?;
        Ok(Divisor {
            a: integer(a),
            d: integer(d),
            z,
        })
    }

    fn embed(&self, model: &mut Model, index: usize, coordinate: &Coordinates) -> Result<Element> {
        if coordinate.denominator != 1 {
            return Err(Error::plain(
                "CM abelian block must belong to the certified integral Hom frame",
            ));
        }
        let [x, y] = coordinate.numerator;
        let [u, v] = self.bases[index];
        let first = self.field.scale(model, u, x, 1)?;
        let second = self.field.scale(model, v, y, 1)?;
        let value = self.field.add(model, first, second)?;
        if self.block_degree(model, index, value)? != integer(coordinate.map.degree() as i128) {
            return Err(inconsistent());
        }
        Ok(value)
    }

    fn block_degree(&self, model: &mut Model, index: usize, entry: Element) -> Result<Element> {
        let degree = self.field.norm(model, entry)?;
        match index {
            1 => self.field.scale(model, degree, 1, self.weight),
            2 => self.field.scale(model, degree, self.weight, 1),
            _ => Ok(degree),
        }
    }

    pub(crate) fn endomorphism(
        &self,
        model: &mut Model,
        entries: [Rc<Morphism>; 4],
    ) -> Result<Endomorphism> {
        let mut matrix = [[Element::ZERO; 2]; 2];
        for (index, map) in entries.into_iter().enumerate() {
            let value = Coordinates::read(model, &self.frames[index], map)?;
            matrix[index / 2][index % 2] = self.embed(model, index, &value)?;
        }
        self.make_endomorphism(model, matrix)
    }

    fn make_endomorphism(&self, model: &mut Model, matrix: Matrix) -> Result<Endomorphism> {
        for (index, &entry) in matrix.iter().flatten().enumerate() {
            let coordinates = self.field.coordinates(model, self.bases[index], entry)?;
            for x in coordinates {
                require_integer(x)?;
            }
        }
        let determinant = self.field.det(model, matrix)?;
        let degree = require_integer(self.field.norm(model, determinant)?)?;
        if degree < 0 {
            return Err(inconsistent());
        }
        Ok(Endomorphism {
            matrix,
            determinant,
            degree,
        })
    }

    pub(crate) fn block(
        &self,
        model: &mut Model,
        endomorphism: Endomorphism,
        index: usize,
    ) -> Result<Morphism> {
        let frame = self.frame(index)?;
        let v = endomorphism.matrix[index / 2][index % 2];
        let [a, b] = self.field.coordinates(model, self.bases[index], v)?;
        frame.combination(model, [require_integer(a)?, require_integer(b)?])
    }

    pub(crate) fn dagger(
        &self,
        model: &mut Model,
        endomorphism: Endomorphism,
    ) -> Result<Endomorphism> {
        let matrix = self.dagger_matrix(model, endomorphism.matrix)?;
        self.make_endomorphism(model, matrix)
    }

    fn dagger_matrix(&self, model: &mut Model, [[a, b], [c, d]]: Matrix) -> Result<Matrix> {
        let b = self.field.conjugate(model, b)?;
        let c = self.field.conjugate(model, c)?;
        Ok([
            [
                self.field.conjugate(model, a)?,
                self.field.scale(model, c, self.weight, 1)?,
            ],
            [
                self.field.scale(model, b, 1, self.weight)?,
                self.field.conjugate(model, d)?,
            ],
        ])
    }

    pub(crate) fn adjugate(
        &self,
        model: &mut Model,
        endomorphism: Endomorphism,
    ) -> Result<Endomorphism> {
        let matrix = self.field.adj(model, endomorphism.matrix)?;
        self.make_endomorphism(model, matrix)
    }

    pub(crate) fn compose(
        &self,
        model: &mut Model,
        a: Endomorphism,
        b: Endomorphism,
    ) -> Result<Endomorphism> {
        let matrix = self.field.multiply(model, a.matrix, b.matrix)?;
        let result = self.make_endomorphism(model, matrix)?;
        if a.degree.checked_mul(b.degree) != Some(result.degree) {
            return Err(inconsistent());
        }
        Ok(result)
    }

    fn hermitian(&self, model: &mut Model, matrix: Matrix) -> Result<Divisor> {
        if matrix != self.dagger_matrix(model, matrix)?
            || matrix[0][0].b != 0
            || matrix[1][1].b != 0
        {
            return Err(inconsistent());
        }
        Ok(Divisor {
            a: matrix[0][0],
            d: matrix[1][1],
            z: matrix[1][0],
        })
    }

    fn matrix(&self, model: &mut Model, v: Divisor) -> Result<Matrix> {
        let adjoint = self.field.conjugate(model, v.z)?;
        let upper = self.field.scale(model, adjoint, self.weight, 1)?;
        Ok([[v.a, upper], [v.z, v.d]])
    }

    pub(crate) fn determinant(&self, model: &mut Model, v: Divisor) -> Result<Element> {
        let conjugate = self.field.conjugate(model, v.z)?;
        rv_qp_quad::product_sum(
            model,
            self.field.0,
            [(1, v.a, v.d), (-self.weight, v.z, conjugate)],
        )
    }

    pub(crate) fn add_divisor(&self, model: &mut Model, a: Divisor, b: Divisor) -> Result<Divisor> {
        Ok(Divisor {
            a: self.field.add(model, a.a, b.a)?,
            d: self.field.add(model, a.d, b.d)?,
            z: self.field.add(model, a.z, b.z)?,
        })
    }

    pub(crate) fn add_mukai(&self, model: &mut Model, a: Mukai, b: Mukai) -> Result<Mukai> {
        Ok(Mukai {
            rank: self.field.add(model, a.rank, b.rank)?,
            divisor: self.add_divisor(model, a.divisor, b.divisor)?,
            last: self.field.add(model, a.last, b.last)?,
        })
    }

    pub(crate) fn scale_divisor(
        &self,
        model: &mut Model,
        v: Divisor,
        scalar: Element,
    ) -> Result<Divisor> {
        Ok(Divisor {
            a: self.field.mul(model, v.a, scalar)?,
            d: self.field.mul(model, v.d, scalar)?,
            z: self.field.mul(model, v.z, scalar)?,
        })
    }

    pub(crate) fn scale_mukai(
        &self,
        model: &mut Model,
        v: Mukai,
        scalar: Element,
    ) -> Result<Mukai> {
        Ok(Mukai {
            rank: self.field.mul(model, v.rank, scalar)?,
            divisor: self.scale_divisor(model, v.divisor, scalar)?,
            last: self.field.mul(model, v.last, scalar)?,
        })
    }

    fn intersection_terms(
        &self,
        model: &mut Model,
        a: Divisor,
        b: Divisor,
    ) -> Result<[(i128, Element, Element); 4]> {
        Ok([
            (1, a.a, b.d),
            (1, b.a, a.d),
            (-self.weight, a.z, self.field.conjugate(model, b.z)?),
            (-self.weight, b.z, self.field.conjugate(model, a.z)?),
        ])
    }

    pub(crate) fn intersection(
        &self,
        model: &mut Model,
        a: Divisor,
        b: Divisor,
    ) -> Result<Element> {
        let terms = self.intersection_terms(model, a, b)?;
        rv_qp_quad::product_sum(model, self.field.0, terms)
    }

    pub(crate) fn pair(&self, model: &mut Model, a: Mukai, b: Mukai) -> Result<Element> {
        let [aa, dd, zz, ww] = self.intersection_terms(model, a.divisor, b.divisor)?;
        rv_qp_quad::product_sum(
            model,
            self.field.0,
            [aa, dd, zz, ww, (-1, a.rank, b.last), (-1, b.rank, a.last)],
        )
    }

    pub(crate) fn line(&self, model: &mut Model, v: Divisor) -> Result<Mukai> {
        let last = self.determinant(model, v)?;
        Ok(Mukai {
            rank: Element::ONE,
            divisor: v,
            last,
        })
    }

    pub(crate) fn fourier(&self, model: &mut Model, v: Mukai) -> Result<Mukai> {
        let divisor = Divisor {
            a: self.field.negate(model, v.divisor.d)?,
            d: self.field.negate(model, v.divisor.a)?,
            z: v.divisor.z,
        };
        Ok(Mukai {
            rank: v.last,
            divisor,
            last: v.rank,
        })
    }

    pub(crate) fn tensor(&self, model: &mut Model, v: Mukai, d: Divisor) -> Result<Mukai> {
        let product = self.scale_divisor(model, d, v.rank)?;
        let divisor = self.add_divisor(model, v.divisor, product)?;
        let h = self.intersection(model, v.divisor, d)?;
        let square = self.determinant(model, d)?;
        let shift = self.field.mul(model, square, v.rank)?;
        let last = self.field.add(model, v.last, h)?;
        let last = self.field.add(model, last, shift)?;
        Ok(Mukai {
            rank: v.rank,
            divisor,
            last,
        })
    }

    pub(crate) fn pull_push(
        &self,
        model: &mut Model,
        v: Mukai,
        map: Endomorphism,
        push: bool,
    ) -> Result<Mukai> {
        let (rank, actual, last) = if push {
            if map.degree == 0 {
                return Err(Error::plain(
                    "CM abelian push requires a nonsingular product isogeny",
                ));
            }
            (
                self.field.scale(model, v.rank, map.degree, 1)?,
                self.field.adj(model, map.matrix)?,
                v.last,
            )
        } else {
            (
                v.rank,
                map.matrix,
                self.field.scale(model, v.last, map.degree, 1)?,
            )
        };
        let h = self.matrix(model, v.divisor)?;
        let dagger = self.dagger_matrix(model, actual)?;
        let left = self.field.multiply(model, dagger, h)?;
        let matrix = self.field.multiply(model, left, actual)?;
        let divisor = self.hermitian(model, matrix)?;
        Ok(Mukai {
            rank,
            divisor,
            last,
        })
    }

    pub(crate) fn graph(&self, model: &mut Model, v: Divisor) -> Result<Cycle> {
        let [u, w] = self.field.coordinates(model, self.bases[2], v.z)?;
        let first = self.field.negate(model, u)?;
        let second = self.field.negate(model, w)?;
        let [f, _, g] = self.frames[2].form;
        let right_f = self.field.scale(model, u, f, 1)?;
        let right_g = self.field.scale(model, w, g, 1)?;
        let delta = self.field.add(model, right_f, right_g)?;
        let right = self.field.add(model, v.a, delta)?;
        let delta = self.field.add(model, u, w)?;
        let left = self.field.add(model, v.d, delta)?;
        Ok([
            require_integer(left)?,
            require_integer(right)?,
            require_integer(first)?,
            require_integer(second)?,
        ])
    }

    pub(crate) fn divisor_from_graph(
        &self,
        model: &mut Model,
        [left, right, f, g]: Cycle,
    ) -> Result<Divisor> {
        let [u, _, w] = self.frames[2].form;
        let first = self.field.scale(model, self.bases[2][0], f, 1)?;
        let second = self.field.scale(model, self.bases[2][1], g, 1)?;
        let z = self.field.add(model, first, second)?;
        let z = self.field.negate(model, z)?;
        let a = self.field.scale(model, integer(f), u, 1)?;
        let b = self.field.scale(model, integer(g), w, 1)?;
        let a = self.field.add(model, a, b)?;
        let a = self.field.add(model, a, integer(right))?;
        let d = self.field.add(model, integer(left), integer(f))?;
        let d = self.field.add(model, d, integer(g))?;
        Ok(Divisor { a, d, z })
    }

    pub(crate) fn report(&self, model: &mut Model, v: Mukai) -> Result<String> {
        let matrix = entries(self.matrix(model, v.divisor)?);
        let square = scalar_text(self.pair(model, v, v)?);
        Ok(format!(
            "rvCMAbelianData = rank={}, matrix={matrix:?}, last={}, square={square}",
            scalar_text(v.rank),
            scalar_text(v.last)
        ))
    }

    pub(crate) fn report_map(&self, map: Endomorphism) -> String {
        let matrix = entries(map.matrix);
        let determinant = triple(map.determinant);
        format!("rvCMAbelianMapData = radicand={}, adjoint_weight={}, matrix={matrix:?}, determinant={determinant:?}, degree={}",self.field.0,self.weight,map.degree)
    }
}

impl Divisor {
    pub(crate) const ZERO: Self = Self {
        a: Element::ZERO,
        d: Element::ZERO,
        z: Element::ZERO,
    };
    pub(crate) fn coordinates(self) -> [[i128; 3]; 3] {
        [triple(self.a), triple(self.d), triple(self.z)]
    }
}

impl Mukai {
    pub(crate) const ZERO: Self = Self {
        rank: Element::ZERO,
        divisor: Divisor::ZERO,
        last: Element::ZERO,
    };
    pub(crate) fn new([r, s]: [i128; 2], divisor: Divisor) -> Self {
        Self {
            rank: integer(r),
            divisor,
            last: integer(s),
        }
    }
    pub(crate) fn divisor(self) -> Divisor {
        self.divisor
    }
    pub(crate) fn coordinates(self) -> [[i128; 3]; 5] {
        let [a, d, z] = self.divisor.coordinates();
        [triple(self.rank), a, d, z, triple(self.last)]
    }
}

impl Endomorphism {
    pub(crate) fn coordinates(self) -> [[[i128; 3]; 2]; 2] {
        entries(self.matrix)
    }
}

fn endomorphism_basis(
    model: &mut Model,
    field: Quadratic,
    frame: &Hom,
    alpha_map: Morphism,
    alpha: Element,
) -> Result<[Element; 2]> {
    let id = Morphism::identity(model, frame.maps[0].source.clone(), &frame.endpoints[0])?;
    let id = Coordinates::read(model, frame, Rc::new(id))?;
    let alpha_map = Coordinates::read(model, frame, Rc::new(alpha_map))?;
    if id.denominator != 1 || alpha_map.denominator != 1 {
        return Err(inconsistent());
    }
    let [a, b] = id.numerator;
    let [c, d] = alpha_map.numerator;
    let determinant = a
        .checked_mul(d)
        .and_then(|x| b.checked_mul(c).and_then(|y| x.checked_sub(y)))
        .ok_or_else(inconsistent)?;
    if determinant.unsigned_abs() != 1 {
        return Err(Error::plain("CM abelian endomorphism frame does not have the expected integral basis of identity and anchored CM action"));
    }
    let basis = [[d, -b], [-c, a]];
    let [u, v] = basis.map(|[real, imag]| -> Result<Element> {
        let component = field.scale(model, alpha, imag, 1)?;
        let total = field.add(model, integer(real), component)?;
        field.scale(model, total, determinant, 1)
    });
    Ok([u?, v?])
}

impl Quadratic {
    fn add(self, model: &mut Model, a: Element, b: Element) -> Result<Element> {
        rv_qp_quad::affine_combination(model, a, b, false)
    }
    fn sub(self, model: &mut Model, a: Element, b: Element) -> Result<Element> {
        rv_qp_quad::affine_combination(model, a, b, true)
    }
    fn mul(self, model: &mut Model, a: Element, b: Element) -> Result<Element> {
        rv_qp_quad::product(model, self.0, a, b)
    }
    fn negate(self, model: &mut Model, a: Element) -> Result<Element> {
        self.sub(model, Element::ZERO, a)
    }
    fn conjugate(self, model: &mut Model, a: Element) -> Result<Element> {
        model.spend(1)?;
        Element::new(a.a, a.b.checked_neg().ok_or_else(inconsistent)?, a.c)
    }
    fn scale(
        self,
        model: &mut Model,
        a: Element,
        numerator: i128,
        denominator: i128,
    ) -> Result<Element> {
        self.mul(model, a, Element::new(numerator, 0, denominator)?)
    }
    fn norm(self, model: &mut Model, a: Element) -> Result<Element> {
        let conjugate = self.conjugate(model, a)?;
        self.mul(model, a, conjugate)
    }
    fn det(self, model: &mut Model, [[a, b], [c, d]]: Matrix) -> Result<Element> {
        rv_qp_quad::product_sum(model, self.0, [(1, a, d), (-1, b, c)])
    }

    fn adj(self, model: &mut Model, [[a, b], [c, d]]: Matrix) -> Result<Matrix> {
        Ok([[d, self.negate(model, b)?], [self.negate(model, c)?, a]])
    }
    fn multiply(self, model: &mut Model, a: Matrix, b: Matrix) -> Result<Matrix> {
        let mut result = [[Element::ZERO; 2]; 2];
        for (i, row) in result.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value = rv_qp_quad::product_sum(
                    model,
                    self.0,
                    [(1, a[i][0], b[0][j]), (1, a[i][1], b[1][j])],
                )?;
            }
        }
        Ok(result)
    }
    fn coordinates(
        self,
        model: &mut Model,
        [u, v]: [Element; 2],
        w: Element,
    ) -> Result<[Element; 2]> {
        let column =
            |x: Element| Ok::<_, Error>([Element::new(x.a, 0, x.c)?, Element::new(x.b, 0, x.c)?]);
        let [x, y, z] = [column(u)?, column(v)?, column(w)?];
        let det = self.det(model, [x, y])?;
        let inverse = det.rational_inverse(model)?;
        let a = self.det(model, [z, y])?;
        let b = self.det(model, [x, z])?;
        Ok([self.mul(model, a, inverse)?, self.mul(model, b, inverse)?])
    }
}

pub(crate) fn scalar_text(x: Element) -> String {
    debug_assert_eq!(x.b, 0);
    if x.c == 1 {
        x.a.to_string()
    } else {
        format!("{}/{}", x.a, x.c)
    }
}

pub(crate) fn integer(x: i128) -> Element {
    Element { a: x, b: 0, c: 1 }
}
fn triple(x: Element) -> [i128; 3] {
    [x.a, x.b, x.c]
}
fn entries(matrix: Matrix) -> [[[i128; 3]; 2]; 2] {
    matrix.map(|row| row.map(triple))
}

fn require_integer(x: Element) -> Result<i128> {
    if x.b != 0 || x.c != 1 {
        return Err(Error::plain(
            "CM abelian coordinate is outside the certified integral lattice",
        ));
    }
    Ok(x.a)
}

fn inconsistent() -> Error {
    Error::plain("CM abelian certification disagrees with the actual Hom maps or Hermitian form")
}
