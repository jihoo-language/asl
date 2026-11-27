use super::{entries, integer, scalar_text, Abelian, Divisor, Endomorphism, Matrix, Mukai};
use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_qp_quad::{self, Element, Wide};
use crate::syntax::{Label, ParsedExpr};

pub(crate) enum Statement {
    Data(Box<ParsedExpr>),
    Ray(Box<[ParsedExpr; 2]>),
    WallData(Box<ParsedExpr>),
    WallAt([i128; 4], Box<ParsedExpr>),
}

pub(crate) enum Expression {
    New(Box<[ParsedExpr; 2]>),
    ChargeLiteral(Label, [i128; 3]),
    Unary(Unary, Box<ParsedExpr>),
    Binary(Binary, Box<[ParsedExpr; 2]>),
    Wall(Box<[ParsedExpr; 4]>),
}

#[derive(Clone, Copy)]
pub(crate) enum Unary {
    Fourier,
    Real,
    Imaginary,
    Factor,
    Determinant,
}

#[derive(Clone, Copy)]
pub(crate) enum Binary {
    Tensor,
    Pull,
    Pair,
    Multiply,
}

// These are independent of the imaginary quadratic CM field. In particular
// this type still describes a central complex scalar when the CM field is Q(i).
#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Charge(Element);

type PeriodCoordinates = ([[i128; 3]; 3], [[i128; 3]; 3], [i128; 3]);

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Period {
    real: Divisor,
    imaginary: Divisor,
    factor: Charge,
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Wall {
    polynomials: [[Element; 3]; 2],
    coefficients: [i128; 3],
}

impl Period {
    pub(crate) fn new(
        model: &mut Model,
        surface: &Abelian,
        real: Divisor,
        imaginary: Divisor,
    ) -> Result<Self> {
        require_positive(model, surface, imaginary)?;
        Ok(Self {
            real,
            imaginary,
            factor: Charge::ONE,
        })
    }

    pub(crate) fn parts(self) -> [Divisor; 2] {
        [self.real, self.imaginary]
    }

    pub(crate) fn factor(self) -> Charge {
        self.factor
    }

    pub(crate) fn determinant(self, model: &mut Model, surface: &Abelian) -> Result<Charge> {
        let a = surface.determinant(model, self.real)?;
        let b = surface.determinant(model, self.imaginary)?;
        let real = sub(model, a, b)?;
        let imaginary = surface.intersection(model, self.real, self.imaginary)?;
        Charge::parts(model, real, imaginary)
    }

    pub(crate) fn pair(self, model: &mut Model, surface: &Abelian, v: Mukai) -> Result<Charge> {
        let det = self.determinant(model, surface)?;
        let b = surface.intersection(model, self.real, v.divisor)?;
        let w = surface.intersection(model, self.imaginary, v.divisor)?;
        let mixed = Charge::parts(model, b, w)?;
        rv_qp_quad::product_sum(
            model,
            -1,
            [
                (1, mixed.0, Element::ONE),
                (-1, v.last, Element::ONE),
                (-1, v.rank, det.0),
            ],
        )
        .map(Charge)
    }

    pub(crate) fn fourier(self, model: &mut Model, surface: &Abelian) -> Result<Self> {
        let determinant = self.determinant(model, surface)?;
        let inverse = determinant.inverse(model)?;
        let [x, y] = self.parts().map(|d| -> Result<Divisor> {
            Ok(surface.fourier(model, Mukai::new([0, 0], d))?.divisor())
        });
        let [x, y] = [x?, y?];
        let [c, d] = inverse.components()?;
        let negative_d = sub(model, Element::ZERO, d)?;
        let real = combine_divisors(model, surface, [(c, x), (negative_d, y)])?;
        let imaginary = combine_divisors(model, surface, [(d, x), (c, y)])?;
        verify_imaginary(model, surface, self, imaginary)?;
        let mut result = Self::new(model, surface, real, imaginary)?;
        let next_determinant = result.determinant(model, surface)?;
        if determinant.multiply(model, next_determinant)? != Charge::ONE {
            return Err(inconsistent());
        }
        result.factor = self.factor.multiply(model, determinant)?;
        Ok(result)
    }

    pub(crate) fn tensor(
        self,
        model: &mut Model,
        surface: &Abelian,
        divisor: Divisor,
    ) -> Result<Self> {
        let real = surface.add_divisor(model, self.real, divisor)?;
        Ok(Self { real, ..self })
    }

    pub(crate) fn pull(
        self,
        model: &mut Model,
        surface: &Abelian,
        map: Endomorphism,
    ) -> Result<Self> {
        if map.degree == 0 {
            return Err(Error::plain(
                "CM period pull requires a nonsingular product isogeny",
            ));
        }
        let [a, b] = self.parts().map(|v| -> Result<Divisor> {
            Ok(surface
                .pull_push(model, Mukai::new([0, 0], v), map, false)?
                .divisor())
        });
        let mut result = Self::new(model, surface, a?, b?)?;
        let det = self.determinant(model, surface)?;
        if result.determinant(model, surface)? != Charge(mul(model, integer(map.degree), det.0)?) {
            return Err(inconsistent());
        }
        result.factor = self.factor;
        Ok(result)
    }

    pub(crate) fn report(self, model: &mut Model, surface: &Abelian) -> Result<String> {
        let real = entries(surface.matrix(model, self.real)?);
        let imaginary = entries(surface.matrix(model, self.imaginary)?);
        let determinant = self.determinant(model, surface)?.coordinates();
        let positive = surface.determinant(model, self.imaginary)?;
        Ok(format!("rvCMPeriodData = real={real:?}, imaginary={imaginary:?}, determinant={determinant:?}, imaginary_determinant={}, factor={:?}",
            scalar_text(positive), self.factor.coordinates()))
    }

    pub(crate) fn coordinates(self) -> PeriodCoordinates {
        (
            self.real.coordinates(),
            self.imaginary.coordinates(),
            self.factor.coordinates(),
        )
    }
}

fn combine_divisors(
    model: &mut Model,
    surface: &Abelian,
    [(a, x), (b, y)]: [(Element, Divisor); 2],
) -> Result<Divisor> {
    let x = surface.scale_divisor(model, x, a)?;
    let y = surface.scale_divisor(model, y, b)?;
    surface.add_divisor(model, x, y)
}

fn require_positive(model: &mut Model, surface: &Abelian, v: Divisor) -> Result<()> {
    if v.a.a <= 0 || surface.determinant(model, v)?.a <= 0 {
        return Err(Error::plain(
            "CM period imaginary divisor must be positive definite",
        ));
    }
    Ok(())
}

fn verify_imaginary(
    model: &mut Model,
    surface: &Abelian,
    input: Period,
    imaginary: Divisor,
) -> Result<()> {
    let field = surface.field;
    let b = surface.matrix(model, input.real)?;
    let w = surface.matrix(model, input.imaginary)?;
    let v = surface.matrix(model, imaginary)?;
    let left = field.multiply(model, b, v)?;
    let right = field.multiply(model, w, v)?;
    let bb = field.multiply(model, left, b)?;
    let ww = field.multiply(model, right, w)?;
    let bw = field.multiply(model, left, w)?;
    let wb = field.multiply(model, right, b)?;
    if matrix_sum(model, bb, ww, false)? != w
        || matrix_sum(model, bw, wb, true)? != [[Element::ZERO; 2]; 2]
    {
        return Err(inconsistent());
    }
    Ok(())
}

fn matrix_sum(model: &mut Model, a: Matrix, b: Matrix, minus: bool) -> Result<Matrix> {
    let mut result = [[Element::ZERO; 2]; 2];
    for (row, (p, q)) in result.iter_mut().zip(a.into_iter().zip(b)) {
        for (entry, (x, y)) in row.iter_mut().zip(p.into_iter().zip(q)) {
            *entry = rv_qp_quad::affine_combination(model, x, y, minus)?;
        }
    }
    Ok(result)
}

impl Charge {
    pub(crate) const ZERO: Self = Self(Element::ZERO);
    const ONE: Self = Self(Element::ONE);

    pub(crate) fn new([a, b, c]: [i128; 3]) -> Result<Self> {
        Ok(Self(Element::new(a, b, c)?))
    }

    fn parts(model: &mut Model, real: Element, imaginary: Element) -> Result<Self> {
        if real.b != 0 || imaginary.b != 0 {
            return Err(inconsistent());
        }
        Ok(Self(add(
            model,
            real,
            Element::new(0, imaginary.a, imaginary.c)?,
        )?))
    }

    fn components(self) -> Result<[Element; 2]> {
        Ok([
            Element::new(self.0.a, 0, self.0.c)?,
            Element::new(self.0.b, 0, self.0.c)?,
        ])
    }

    pub(crate) fn multiply(self, model: &mut Model, other: Self) -> Result<Self> {
        Ok(Self(mul(model, self.0, other.0)?))
    }

    fn crossed(self, other: Self) -> Result<Wide<16>> {
        Wide::product(self.0.b, other.0.a).add(Wide::product(self.0.a, other.0.b).negate())
    }

    fn inverse(self, model: &mut Model) -> Result<Self> {
        if self == Self::ZERO {
            return Err(Error::plain(
                "positive CM period has a singular complex determinant",
            ));
        }
        let n = Wide::<16>::product(self.0.a, self.0.a).add(Wide::product(self.0.b, self.0.b))?;
        rv_qp_quad::normalize_wide(
            model,
            Wide::product(self.0.a, self.0.c),
            Wide::product(self.0.b, self.0.c).negate(),
            n,
        )
        .map(Self)
    }

    pub(crate) fn coordinates(self) -> [i128; 3] {
        [self.0.a, self.0.b, self.0.c]
    }

    pub(crate) fn ray(self, model: &mut Model, other: Self) -> Result<String> {
        let [zero_self, zero_other] = [self == Self::ZERO, other == Self::ZERO];
        if zero_self || zero_other {
            return Ok(format!("orientation=zero, first_zero={zero_self}, second_zero={zero_other}, ratio=undefined"));
        }
        model.spend(8)?;
        if !self.crossed(other)?.is_zero() {
            return Ok(
                "orientation=different_rays, first_zero=false, second_zero=false, ratio=undefined"
                    .into(),
            );
        }
        let dot =
            Wide::<16>::product(self.0.a, other.0.a).add(Wide::product(self.0.b, other.0.b))?;
        let norm =
            Wide::<16>::product(other.0.a, other.0.a).add(Wide::product(other.0.b, other.0.b))?;
        let label = if dot.negative {
            "opposite_ray"
        } else {
            "same_ray"
        };
        let ratio = rv_qp_quad::normalize_wide(
            model,
            dot.multiply(other.0.c)?,
            Wide::from_i128(0),
            norm.multiply(self.0.c)?,
        )?;
        Ok(format!(
            "orientation={label}, first_zero=false, second_zero=false, ratio={}",
            scalar_text(ratio)
        ))
    }
}

impl Wall {
    pub(crate) fn new(
        model: &mut Model,
        surface: &Abelian,
        l: Divisor,
        beta: Divisor,
        v: [Mukai; 2],
    ) -> Result<Self> {
        require_positive(model, surface, l)?;
        let det_l = surface.determinant(model, l)?;
        let bisector = surface.intersection(model, beta, l)?;
        let det_beta = surface.determinant(model, beta)?;
        let [x, y] = v.map(|v| -> Result<[Element; 3]> {
            let a = mul(model, v.rank, det_l)?;
            let a = sub(model, Element::ZERO, a)?;
            let b = surface.intersection(model, l, v.divisor)?;
            let c = surface.intersection(model, beta, v.divisor)?;
            let b =
                rv_qp_quad::product_sum(model, -1, [(1, b, Element::ONE), (-1, v.rank, bisector)])?;
            let c = rv_qp_quad::product_sum(
                model,
                -1,
                [
                    (1, c, Element::ONE),
                    (-1, v.last, Element::ONE),
                    (-1, v.rank, det_beta),
                ],
            )?;
            Ok([a, b, c])
        });
        let [[a, b, c], [d, e, f]] = [x?, y?];
        let terms = [
            rv_qp_quad::product_sum(model, -1, [(1, a, e), (-1, b, d)])?,
            rv_qp_quad::product_sum(model, -1, [(2, a, f), (-2, c, d)])?,
            rv_qp_quad::product_sum(model, -1, [(1, b, f), (-1, c, e)])?,
        ];
        let coefficients = primitive(model, terms)?;
        Ok(Self {
            polynomials: [[a, b, c], [d, e, f]],
            coefficients,
        })
    }

    pub(crate) fn report(self, model: &mut Model) -> Result<String> {
        let [a, b, c] = self.coefficients;
        let geometry = if a != 0 {
            let q = divide(model, integer(b), integer(a))?;
            let r = divide(model, integer(c), integer(a))?;
            let center = mul(model, Element::new(-1, 0, 2)?, q)?;
            let radius =
                rv_qp_quad::product_sum(model, -1, [(1, center, center), (-1, r, Element::ONE)])?;
            format!(
                "kind=quadratic, center={}, radius_squared={}, positive_height_locus={}",
                scalar_text(center),
                scalar_text(radius),
                radius.a > 0
            )
        } else if b != 0 {
            let location = divide(model, integer(c), integer(b))?;
            let location = sub(model, Element::ZERO, location)?;
            format!("kind=linear, real_coordinate={}", scalar_text(location))
        } else if c == 0 {
            "kind=identically_collinear".into()
        } else {
            "kind=empty".into()
        };
        let polynomials = self.polynomials.map(|p| p.map(scalar_text));
        let quadratics = polynomials.map(|p| format!("[{}, {}, {}]", p[0], p[1], p[2]));
        Ok(format!("rvCMPeriodWallData = {geometry}, coefficients={:?}, charge_quadratics=[{}, {}], orientation_test=sign(Re_v*Re_w+t^2*(2*A_v*b+B_v)*(2*A_w*b+B_w))",
            self.coefficients,quadratics[0],quadratics[1]))
    }

    pub(crate) fn at(self, model: &mut Model, [bn, bd, tn, td]: [i128; 4]) -> Result<String> {
        let b = Element::new(bn, 0, bd)?;
        let t = Element::new(tn, 0, td)?;
        if t.a <= 0 {
            return Err(Error::plain("CM numerical wall height must be positive"));
        }
        let z = Charge::parts(model, b, t)?;
        let [x, y] = self.polynomials.map(|[a, b, c]| -> Result<Charge> {
            let first = mul(model, a, z.0)?;
            let next = add(model, first, b)?;
            let second = mul(model, next, z.0)?;
            Ok(Charge(add(model, second, c)?))
        });
        let [x, y] = [x?, y?];
        let [p, q, r] = self.coefficients;
        let evaluate = rv_qp_quad::product_sum(
            model,
            -1,
            [
                (p, b, b),
                (p, t, t),
                (q, b, Element::ONE),
                (r, Element::ONE, Element::ONE),
            ],
        )?;
        let on_locus = evaluate.is_zero();
        if on_locus != x.crossed(y)?.is_zero() {
            return Err(inconsistent());
        }
        let orientation = x.ray(model, y)?;
        Ok(format!(
            "rvCMPeriodWallAt = on_locus={on_locus}, {orientation}, charges=[{:?}, {:?}]",
            x.coordinates(),
            y.coordinates()
        ))
    }

    pub(crate) fn coefficients(self) -> [i128; 3] {
        self.coefficients
    }
}

fn primitive(model: &mut Model, values: [Element; 3]) -> Result<[i128; 3]> {
    type Coefficient = Wide<16>;
    let [a, b, c] = values;
    let numerators = [
        Coefficient::product(a.a, b.c).multiply(c.c)?,
        Coefficient::product(b.a, a.c).multiply(c.c)?,
        Coefficient::product(c.a, a.c).multiply(b.c)?,
    ];
    let mut common = crate::wide_integer::Natural::ZERO;
    for value in numerators {
        common = common.gcd(model, value.magnitude)?;
    }
    if common.is_zero() {
        return Ok([0; 3]);
    }
    let negative = numerators
        .iter()
        .find(|x| !x.is_zero())
        .is_some_and(|x| x.negative);
    let common = Coefficient {
        magnitude: common,
        negative,
    };
    let mut output = [0; 3];
    for (entry, n) in output.iter_mut().zip(numerators) {
        *entry = n
            .divide_exact(model, common)?
            .and_then(Coefficient::signed)
            .ok_or_else(overflow)?;
    }
    Ok(output)
}

fn add(model: &mut Model, a: Element, b: Element) -> Result<Element> {
    rv_qp_quad::affine_combination(model, a, b, false)
}
fn sub(model: &mut Model, a: Element, b: Element) -> Result<Element> {
    rv_qp_quad::affine_combination(model, a, b, true)
}
fn mul(model: &mut Model, a: Element, b: Element) -> Result<Element> {
    rv_qp_quad::product(model, -1, a, b)
}
fn divide(model: &mut Model, a: Element, b: Element) -> Result<Element> {
    let inverse = b.rational_inverse(model)?;
    mul(model, a, inverse)
}
fn inconsistent() -> Error {
    Error::plain(
        "CM period transform disagrees with its exact Hermitian identity or wall certificate",
    )
}
fn overflow() -> Error {
    Error::plain("CM period coordinates exceed the exact signed integer range")
}
