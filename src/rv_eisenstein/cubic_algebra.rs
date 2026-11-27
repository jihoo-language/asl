use std::rc::Rc;

use crate::error::{Error, Result};
use crate::hw_algebra::{self, CoefficientRing, Matrix};
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::cubic_hilbert::{self as cubic, Local};
use super::kummer::inverse_base_unit;
use super::{Algebra as BaseAlgebra, Element as BaseElement, Extension as Base};

pub(crate) struct Declaration(
    pub(crate) Label,
    pub(crate) cubic::LocalSpec,
    pub(crate) Box<[ParsedExpr; 2]>,
);

pub(crate) struct SplitDeclaration(
    pub(crate) Label,
    pub(crate) Label,
    pub(crate) Vec<Vec<i128>>,
);

pub(crate) enum Expression {
    Literal(Label, Vec<Vec<i128>>),
    Basis(Label, usize, usize),
    Term(Label, usize, usize, Box<ParsedExpr>),
    SplitStatus(Label),
    Multiply(Box<[ParsedExpr; 2]>),
    Unary(Unary, Box<ParsedExpr>),
    SplitImage(Label, Box<ParsedExpr>),
    SplitMultiply(Box<[ParsedExpr; 2]>),
    SplitEntry(usize, usize, Box<ParsedExpr>),
    SplitDeterminant(Box<ParsedExpr>),
}

pub(crate) enum Unary {
    Coefficient(usize, usize),
    Trace,
    Norm,
    Adjugate,
    Inverse,
    Negate,
}

pub(crate) enum Report {
    Data(Label),
    SplitData(Label),
    Characteristic(Box<ParsedExpr>),
}

pub(crate) struct Parameter {
    pub(crate) local: Local,
    pub(crate) digits: usize,
}

impl Parameter {
    pub(crate) fn from_integral(
        model: &mut Model,
        field: &Base,
        source: usize,
        value: BaseElement,
    ) -> Result<Self> {
        let (local, digits) = Local::factor_integral(model, field, source, value)?;
        Ok(Self { local, digits })
    }
}

pub(crate) struct Extension {
    pub(crate) base_name: Box<str>,
    pub(crate) precision: usize,
    pub(crate) base: Rc<Base>,
    parameters: [BaseElement; 2],
    shifts: [i128; 2],
    zeta: [BaseElement; 3],
    exponent: u8,
    representation: Vec<Matrix<CubicScalar>>,
}

impl Extension {
    pub(crate) fn new(
        model: &mut Model,
        grade: cubic::Grade,
        base: Rc<Base>,
        [a, b]: [Parameter; 2],
    ) -> Result<Self> {
        let algebra = base.algebra(grade.precision)?;
        let normalize = |model: &mut Model, parameter: Parameter| {
            let valuation = parameter.local.valuation;
            let (order, shift) = if valuation < 0 {
                (valuation.rem_euclid(3) as u128, -valuation.div_euclid(3))
            } else {
                (valuation as u128, 0)
            };
            let r = order.min(grade.precision as u128) as usize;
            if parameter.digits < cubic::CLASS_PRECISION.max(grade.precision.saturating_sub(r)) {
                return Err(Error::plain(
                    "cubic algebra parameter has insufficient relative precision",
                ));
            }
            let class = parameter.local.class(model, &base)?;
            let power = algebra.power(model, algebra.pi(), r as u64)?;
            let unit = algebra.reduce(parameter.local.unit);
            let reduced = algebra.multiply(model, &power, &unit)?;
            Ok((class, reduced, shift))
        };
        let (a_class, a, us) = normalize(model, a)?;
        let (b_class, b, vs) = normalize(model, b)?;
        let exponent = model.cubic_certificate()?.exponent(a_class, b_class);
        let zeta = algebra.add(model, &algebra.scalar(1), &algebra.pi())?;
        let squared = algebra.multiply(model, &zeta, &zeta)?;
        let zeta = [algebra.scalar(1), zeta, squared];
        let etale = Etale {
            base: Coefficients(&algebra),
            a: &a,
        };
        let root = etale.alpha();
        let u = diagonal(
            &etale,
            &[
                root.clone(),
                etale.multiply(model, &root, &etale.constant(zeta[1].clone()))?,
                etale.multiply(model, &root, &etale.constant(zeta[2].clone()))?,
            ],
        )?;
        let v = Matrix::from_entries(
            3,
            vec![
                etale.zero(),
                etale.one(),
                etale.zero(),
                etale.zero(),
                etale.zero(),
                etale.one(),
                etale.constant(b.clone()),
                etale.zero(),
                etale.zero(),
            ],
        )?;
        let representation = crossed_representation(
            model,
            &etale,
            u,
            v,
            [
                &etale.constant(a.clone()),
                &etale.constant(b.clone()),
                &etale.constant(zeta[1].clone()),
            ],
        )?;
        model.allocate(9 * 9 * 3 * 4 + 32)?;
        Ok(Self {
            base_name: grade.field,
            precision: grade.precision,
            base,
            parameters: [a, b],
            shifts: [us, vs],
            zeta,
            exponent,
            representation,
        })
    }

    pub(crate) fn is_split(&self) -> bool {
        self.exponent == 0
    }

    pub(super) fn exponent(&self) -> u8 {
        self.exponent
    }

    pub(super) fn parameters(&self) -> &[BaseElement; 2] {
        &self.parameters
    }

    pub(crate) fn report(&self, name: &str) -> String {
        let kind = if self.is_split() { "split" } else { "division" };
        let coordinates = self.parameters.each_ref().map(coordinates);
        format!("rvCubicAlgData[{name}] = base={}, precision={}, parameters={coordinates:?}, shifts={:?}, invariant={}, kind={kind}",
            self.base_name, self.precision, self.shifts, invariant(self.exponent))
    }

    pub(crate) fn algebra(&self) -> Result<Algebra<'_>> {
        Ok(Algebra {
            extension: self,
            base: self.base.algebra(self.precision)?,
        })
    }
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Element([BaseElement; 9]);

impl Element {
    pub(crate) fn coefficient(&self, i: usize) -> BaseElement {
        self.0[i].clone()
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(45)?;
        Ok(Value::CubicAlg(Rc::new(self)))
    }

    pub(crate) fn rows(&self) -> [[u64; 2]; 9] {
        self.0.each_ref().map(coordinates)
    }
}

pub(crate) fn index(i: usize, j: usize) -> Result<usize> {
    if i >= 3 || j >= 3 {
        return Err(Error::plain(
            "cubic algebra basis and matrix indices must be in 0..=2",
        ));
    }
    Ok(3 * i + j)
}

pub(crate) struct Algebra<'a> {
    extension: &'a Extension,
    base: BaseAlgebra<'a>,
}

impl Algebra<'_> {
    pub(crate) fn zero(&self) -> Element {
        Element(std::array::from_fn(|_| self.base.zero()))
    }

    pub(crate) fn literal(&self, model: &mut Model, rows: Vec<Vec<i128>>) -> Result<Element> {
        if rows.len() > 9 {
            return Err(Error::plain(
                "cubic algebra literals accept up to nine basis coefficients",
            ));
        }
        let mut output = self.zero();
        for (slot, row) in output.0.iter_mut().zip(rows) {
            *slot = self
                .base
                .literal(model, row.into_iter().map(|x| vec![x]).collect())?;
        }
        Ok(output)
    }

    pub(crate) fn term(&self, index: usize, coefficient: BaseElement) -> Element {
        let mut output = self.zero();
        output.0[index] = coefficient;
        output
    }

    pub(crate) fn basis(&self, index: usize) -> Element {
        self.term(index, self.base.scalar(1))
    }

    pub(crate) fn add(&self, model: &mut Model, x: &Element, y: &Element) -> Result<Element> {
        let mut output = self.zero();
        for (out, (a, b)) in output.0.iter_mut().zip(x.0.iter().zip(&y.0)) {
            *out = self.base.add(model, a, b)?;
        }
        Ok(output)
    }

    pub(crate) fn negate(&self, model: &mut Model, x: &Element) -> Result<Element> {
        let mut output = self.zero();
        for (out, a) in output.0.iter_mut().zip(&x.0) {
            *out = self.base.negate(model, a)?;
        }
        Ok(output)
    }

    fn scale(&self, model: &mut Model, x: &Element, a: &BaseElement) -> Result<Element> {
        let mut output = self.zero();
        for (out, b) in output.0.iter_mut().zip(&x.0) {
            *out = self.base.multiply(model, a, b)?;
        }
        Ok(output)
    }

    pub(crate) fn multiply(&self, model: &mut Model, x: &Element, y: &Element) -> Result<Element> {
        model.spend(81)?;
        let mut output = self.zero();
        for (m, left) in x.0.iter().enumerate().filter(|(_, a)| !zero(a)) {
            let (i, j) = (m / 3, m % 3);
            for (n, right) in y.0.iter().enumerate().filter(|(_, a)| !zero(a)) {
                let (k, l) = (n / 3, n % 3);
                let mut term = self.base.multiply(model, left, right)?;
                let phase = j * k % 3;
                if phase != 0 {
                    term = self
                        .base
                        .multiply(model, &term, &self.extension.zeta[phase])?;
                }
                if i + k >= 3 {
                    term = self
                        .base
                        .multiply(model, &term, &self.extension.parameters[0])?;
                }
                if j + l >= 3 {
                    term = self
                        .base
                        .multiply(model, &term, &self.extension.parameters[1])?;
                }
                let index = 3 * ((i + k) % 3) + (j + l) % 3;
                output.0[index] = self.base.add(model, &output.0[index], &term)?;
            }
        }
        Ok(output)
    }

    fn etale(&self) -> Etale<'_, '_> {
        Etale {
            base: Coefficients(&self.base),
            a: &self.extension.parameters[0],
        }
    }

    fn reduced_data(&self, model: &mut Model, x: &Element) -> Result<Reduced> {
        let ring = self.etale();
        let coeffs = x.0.each_ref().map(|a| ring.constant(a.clone()));
        let matrix = image(model, &ring, &self.extension.representation, &coeffs)?;
        let characteristic = hw_algebra::characteristic(model, &ring, &matrix)?;
        let mut base = Vec::with_capacity(4);
        for coefficient in characteristic {
            if coefficient.0[1..].iter().any(|x| !zero(x)) {
                return Err(inconsistent());
            }
            base.push(coefficient.0[0].clone());
        }

        let square = self.multiply(model, x, x)?;
        let linear = self.scale(model, x, &base[1])?;
        let top = self.add(model, &square, &linear)?;
        let adjugate = self.add(model, &top, &self.term(0, base[2].clone()))?;
        let norm = self.base.negate(model, &base[3])?;
        if self.multiply(model, x, &adjugate)? != self.term(0, norm.clone()) {
            return Err(inconsistent());
        }
        Ok(Reduced {
            characteristic: base,
            norm,
            adjugate,
        })
    }

    pub(crate) fn characteristic(&self, model: &mut Model, x: &Element) -> Result<String> {
        let data = self.reduced_data(model, x)?;
        let coefficients: Vec<_> = data.characteristic.iter().map(coordinates).collect();
        Ok(format!(
            "rvCubicReducedChar = precision={}, coefficients={coefficients:?}",
            self.extension.precision
        ))
    }

    pub(crate) fn trace(&self, model: &mut Model, x: &Element) -> Result<BaseElement> {
        let data = self.reduced_data(model, x)?;
        self.base.negate(model, &data.characteristic[1])
    }

    pub(crate) fn norm(&self, model: &mut Model, x: &Element) -> Result<BaseElement> {
        Ok(self.reduced_data(model, x)?.norm)
    }

    pub(crate) fn adjugate(&self, model: &mut Model, x: &Element) -> Result<Element> {
        Ok(self.reduced_data(model, x)?.adjugate)
    }

    pub(crate) fn inverse(&self, model: &mut Model, x: &Element) -> Result<Element> {
        let data = self.reduced_data(model, x)?;
        let valuation = self.base.valuation(model, &data.norm)?;
        if valuation == self.extension.precision {
            return Err(Error::plain(
                "reduced norm vanishes at the available precision. A nonzero norm is uncertified",
            ));
        }
        if valuation != 0 {
            return Err(Error::plain("reduced norm has positive valuation. The inverse leaves the integral cubic algebra chart"));
        }
        let norm_inverse = inverse_base_unit(model, &self.base, &data.norm)?;
        let inverse = self.scale(model, &data.adjugate, &norm_inverse)?;
        if self.multiply(model, x, &inverse)? != self.basis(0)
            || self.multiply(model, &inverse, x)? != self.basis(0)
        {
            return Err(inconsistent());
        }
        Ok(inverse)
    }
}

struct Reduced {
    characteristic: Vec<BaseElement>,
    norm: BaseElement,
    adjugate: Element,
}

pub(crate) struct Split {
    pub(crate) name_of_algebra: Box<str>,
    pub(crate) extension: Rc<Extension>,
    matrices: Vec<Matrix<BaseElement>>,
    norm: BaseElement,
}

impl Split {
    pub(crate) fn new(
        model: &mut Model,
        name_of_algebra: Box<str>,
        extension: Rc<Extension>,
        rows: Vec<Vec<i128>>,
    ) -> Result<Self> {
        if !extension.is_split() {
            return Err(Error::plain(
                "cubic division algebra cannot admit a split matrix witness",
            ));
        }
        if rows.len() > 3 {
            return Err(Error::plain(
                "a cubic split witness takes up to three base coefficients in alpha",
            ));
        }
        let algebra = extension.algebra()?;
        let etale = algebra.etale();
        let mut c = etale.zero();
        for (slot, row) in c.0.iter_mut().zip(rows) {
            *slot = algebra
                .base
                .literal(model, row.into_iter().map(|x| vec![x]).collect())?;
        }
        let multiplication = etale.multiplication_matrix(model, &c)?;
        let norm = determinant(model, &etale.base, &multiplication)?;
        if norm != extension.parameters[1] {
            return Err(Error::plain(
                "cubic split witness norm differs from b at the chart precision",
            ));
        }
        let u = etale.multiplication_matrix(model, &etale.alpha())?;
        let sigma = diagonal(&etale.base, &extension.zeta)?;
        let v = hw_algebra::multiply(model, &etale.base, &multiplication, &sigma)?;
        let matrices = crossed_representation(
            model,
            &etale.base,
            u,
            v,
            [
                &extension.parameters[0],
                &extension.parameters[1],
                &extension.zeta[1],
            ],
        )?;
        model.allocate(9 * 9 * 5 + 10)?;
        Ok(Self {
            name_of_algebra,
            extension,
            matrices,
            norm,
        })
    }

    pub(crate) fn report(&self, name: &str) -> String {
        format!(
            "rvCubicSplitData[{name}] = algebra={}, precision={}, norm_c={:?}, relations=verified",
            self.name_of_algebra,
            self.extension.precision,
            coordinates(&self.norm)
        )
    }

    pub(crate) fn image(&self, model: &mut Model, x: &Element) -> Result<SplitMatrix> {
        let algebra = self.extension.algebra()?;
        let ring = Coefficients(&algebra.base);
        let matrix = image(model, &ring, &self.matrices, &x.0)?;
        let characteristic = hw_algebra::characteristic(model, &ring, &matrix)?;
        if characteristic != algebra.reduced_data(model, x)?.characteristic {
            return Err(inconsistent());
        }
        Ok(SplitMatrix(matrix))
    }

    pub(crate) fn multiply(
        &self,
        model: &mut Model,
        x: &SplitMatrix,
        y: &SplitMatrix,
    ) -> Result<SplitMatrix> {
        let algebra = self.extension.algebra()?;
        Ok(SplitMatrix(hw_algebra::multiply(
            model,
            &Coefficients(&algebra.base),
            &x.0,
            &y.0,
        )?))
    }

    pub(crate) fn determinant(&self, model: &mut Model, x: &SplitMatrix) -> Result<BaseElement> {
        let algebra = self.extension.algebra()?;
        determinant(model, &Coefficients(&algebra.base), &x.0)
    }
}

pub(crate) struct SplitMatrix(Matrix<BaseElement>);

impl SplitMatrix {
    pub(crate) fn coefficient(&self, index: usize) -> BaseElement {
        self.0.entries()[index].clone()
    }

    pub(crate) fn rows(&self) -> [[u64; 2]; 9] {
        std::array::from_fn(|i| coordinates(&self.0.entries()[i]))
    }

    pub(crate) fn entries(&self) -> &[BaseElement] {
        self.0.entries()
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(45)?;
        Ok(Value::CubicSplitMatrix(Rc::new(self)))
    }
}

#[derive(Clone, Copy)]
struct Coefficients<'a, 'b>(&'a BaseAlgebra<'b>);

impl CoefficientRing for Coefficients<'_, '_> {
    type Element = BaseElement;
    fn cells(&self) -> usize {
        4
    }
    fn zero(&self) -> BaseElement {
        self.0.zero()
    }
    fn one(&self) -> BaseElement {
        self.0.scalar(1)
    }
    fn is_zero(&self, x: &BaseElement) -> bool {
        zero(x)
    }
    fn is_unit(&self, x: &BaseElement) -> bool {
        !x.0[0][0].is_multiple_of(3)
    }
    fn add(&self, model: &mut Model, a: &BaseElement, b: &BaseElement) -> Result<BaseElement> {
        self.0.add(model, a, b)
    }
    fn subtract(&self, model: &mut Model, a: &BaseElement, b: &BaseElement) -> Result<BaseElement> {
        let negative = self.0.negate(model, b)?;
        self.0.add(model, a, &negative)
    }
    fn multiply(&self, model: &mut Model, a: &BaseElement, b: &BaseElement) -> Result<BaseElement> {
        self.0.multiply(model, a, b)
    }
    fn inverse_unit(&self, model: &mut Model, x: &BaseElement) -> Result<BaseElement> {
        inverse_base_unit(model, self.0, x)
    }
}

#[derive(Clone, Eq, PartialEq)]
struct CubicScalar([BaseElement; 3]);

struct Etale<'a, 'b> {
    base: Coefficients<'a, 'b>,
    a: &'a BaseElement,
}

impl Etale<'_, '_> {
    fn constant(&self, a: BaseElement) -> CubicScalar {
        CubicScalar([a, self.base.zero(), self.base.zero()])
    }
    fn alpha(&self) -> CubicScalar {
        CubicScalar([self.base.zero(), self.base.one(), self.base.zero()])
    }
    fn multiplication_matrix(
        &self,
        model: &mut Model,
        a: &CubicScalar,
    ) -> Result<Matrix<BaseElement>> {
        let mut entries = vec![self.base.zero(); 9];
        let mut power = self.one();
        for j in 0..3 {
            let column = self.multiply(model, a, &power)?;
            for (i, x) in column.0.into_iter().enumerate() {
                entries[3 * i + j] = x;
            }
            if j < 2 {
                power = self.multiply(model, &power, &self.alpha())?;
            }
        }
        Matrix::from_entries(3, entries)
    }
}

impl CoefficientRing for Etale<'_, '_> {
    type Element = CubicScalar;
    fn cells(&self) -> usize {
        3 * self.base.cells()
    }
    fn zero(&self) -> CubicScalar {
        self.constant(self.base.zero())
    }
    fn one(&self) -> CubicScalar {
        self.constant(self.base.one())
    }
    fn is_zero(&self, x: &CubicScalar) -> bool {
        x.0.iter().all(zero)
    }
    fn is_unit(&self, x: &CubicScalar) -> bool {
        let a = self.a.0[0][0] % 3;
        let [x0, x1, x2] = x.0.each_ref().map(|x| x.0[0][0] % 3);
        !(x0 + a * x1 + a * a * x2).is_multiple_of(3)
    }
    fn add(&self, model: &mut Model, a: &CubicScalar, b: &CubicScalar) -> Result<CubicScalar> {
        let mut out = self.zero();
        for (x, (a, b)) in out.0.iter_mut().zip(a.0.iter().zip(&b.0)) {
            *x = self.base.add(model, a, b)?;
        }
        Ok(out)
    }
    fn subtract(&self, model: &mut Model, a: &CubicScalar, b: &CubicScalar) -> Result<CubicScalar> {
        let mut out = self.zero();
        for (x, (a, b)) in out.0.iter_mut().zip(a.0.iter().zip(&b.0)) {
            *x = self.base.subtract(model, a, b)?;
        }
        Ok(out)
    }
    fn multiply(&self, model: &mut Model, a: &CubicScalar, b: &CubicScalar) -> Result<CubicScalar> {
        let mut out = self.zero();
        model.spend(9)?;
        for (i, x) in a.0.iter().enumerate().filter(|(_, x)| !zero(x)) {
            for (j, y) in b.0.iter().enumerate().filter(|(_, x)| !zero(x)) {
                let mut prod = self.base.multiply(model, x, y)?;
                if i + j >= 3 {
                    prod = self.base.multiply(model, &prod, self.a)?;
                }
                let slot = &mut out.0[(i + j) % 3];
                *slot = self.base.add(model, slot, &prod)?;
            }
        }
        Ok(out)
    }
    fn inverse_unit(&self, model: &mut Model, x: &CubicScalar) -> Result<CubicScalar> {
        let matrix = self.multiplication_matrix(model, x)?;
        let c = hw_algebra::characteristic(model, &self.base, &matrix)?;
        let negative = self.base.subtract(model, &self.base.zero(), &c[3])?;
        let inverse = self.base.inverse_unit(model, &negative)?;
        let square = self.multiply(model, x, x)?;
        let linear = self.multiply(model, x, &self.constant(c[1].clone()))?;
        let polynomial = self.add(model, &square, &linear)?;
        let adj = self.add(model, &polynomial, &self.constant(c[2].clone()))?;
        self.multiply(model, &adj, &self.constant(inverse))
    }
}

fn crossed_representation<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    u: Matrix<R::Element>,
    v: Matrix<R::Element>,
    [a, b, zeta]: [&R::Element; 3],
) -> Result<Vec<Matrix<R::Element>>> {
    let identity = diagonal(ring, &[ring.one(), ring.one(), ring.one()])?;
    let squared_u = hw_algebra::multiply(model, ring, &u, &u)?;
    let squared_v = hw_algebra::multiply(model, ring, &v, &v)?;
    if hw_algebra::multiply(model, ring, &squared_u, &u)?
        != diagonal(ring, &[a.clone(), a.clone(), a.clone()])?
        || hw_algebra::multiply(model, ring, &squared_v, &v)?
            != diagonal(ring, &[b.clone(), b.clone(), b.clone()])?
    {
        return Err(inconsistent());
    }
    let uv = hw_algebra::multiply(model, ring, &u, &v)?;
    let scaled = uv
        .entries()
        .iter()
        .map(|x| ring.multiply(model, x, zeta))
        .collect::<Result<_>>()?;
    if hw_algebra::multiply(model, ring, &v, &u)? != Matrix::from_entries(3, scaled)? {
        return Err(inconsistent());
    }
    let upowers = [identity.clone(), u, squared_u];
    let vpowers = [identity, v, squared_v];
    let mut answer = Vec::with_capacity(9);
    for ui in &upowers {
        for vj in &vpowers {
            answer.push(hw_algebra::multiply(model, ring, ui, vj)?);
        }
    }
    Ok(answer)
}

fn image<R: CoefficientRing>(
    model: &mut Model,
    ring: &R,
    matrices: &[Matrix<R::Element>],
    coefficients: &[R::Element],
) -> Result<Matrix<R::Element>> {
    let mut output = vec![ring.zero(); 9];
    for (matrix, coefficient) in matrices
        .iter()
        .zip(coefficients)
        .filter(|(_, c)| !ring.is_zero(c))
    {
        for (out, x) in output
            .iter_mut()
            .zip(matrix.entries())
            .filter(|(_, x)| !ring.is_zero(x))
        {
            let term = ring.multiply(model, coefficient, x)?;
            *out = ring.add(model, out, &term)?;
        }
    }
    Matrix::from_entries(3, output)
}

fn diagonal<R: CoefficientRing>(
    ring: &R,
    diagonal: &[R::Element; 3],
) -> Result<Matrix<R::Element>> {
    let mut matrix = vec![ring.zero(); 9];
    for (i, x) in diagonal.iter().enumerate() {
        matrix[4 * i] = x.clone();
    }
    Matrix::from_entries(3, matrix)
}

fn determinant(
    model: &mut Model,
    ring: &Coefficients<'_, '_>,
    matrix: &Matrix<BaseElement>,
) -> Result<BaseElement> {
    let c = hw_algebra::characteristic(model, ring, matrix)?;
    ring.subtract(model, &ring.zero(), &c[3])
}

fn coordinates(x: &BaseElement) -> [u64; 2] {
    [x.0[0][0], x.0[1][0]]
}

fn zero(x: &BaseElement) -> bool {
    coordinates(x) == [0, 0]
}

fn invariant(exponent: u8) -> String {
    if exponent == 0 {
        "0".to_owned()
    } else {
        format!("{exponent}/3")
    }
}

fn inconsistent() -> Error {
    Error::plain("inconsistent cubic algebra crossed-product or reduced characteristic identity")
}
