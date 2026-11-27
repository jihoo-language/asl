use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::rv_lens::Rational;
use crate::syntax::{Label, ParsedExpr};

use super::super::super::matrix::{self, Field, Matrix};
use super::{Category, Double};

mod build;
mod stack;
use build::Builder;

const MAX_TUBE: usize = 128;

pub(crate) enum Command {
    Declare(Label, Label),
    Data(Label),
    Annular(Label, [usize; 2]),
    BasisData(Label, usize),
    BlockData(Label, usize),
    Stack(Label, [usize; 2]),
    Decompose(Box<ParsedExpr>),
}

pub(crate) enum Expression {
    Zero(Label),
    Unit(Label),
    Basis(Label, usize),
    Central(Label, usize),
    MatrixUnit(Label, [usize; 3]),
    Coefficient(usize, Box<ParsedExpr>),
    BlockEntry([usize; 3], Box<ParsedExpr>),
    Add(Box<[ParsedExpr; 2]>),
    Product(Box<[ParsedExpr; 2]>),
    Scale(Box<[ParsedExpr; 2]>),
}

pub(crate) struct Tube {
    pub(crate) category: Label,
    pub(crate) conductor: usize,
    forms: Vec<Form>,
    blocks: Vec<Block>,
    forward: Matrix,
    inverse: Matrix,
    unit: Element,
    boundary_simples: usize,
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Element(pub(crate) Vec<(usize, Rational)>);

#[derive(Clone, Copy, Eq, PartialEq)]
struct Form {
    j: usize,
    k: usize,
    i: usize,
    c: usize,
    u: usize,
    v: usize,
}

struct Block {
    start: usize,
    slots: Vec<[usize; 2]>,
}

impl Block {
    fn size(&self) -> usize {
        self.slots.len()
    }

    fn index(&self, r: usize, c: usize) -> usize {
        self.start + r * self.size() + c
    }
}

pub(super) fn bound() -> Error {
    Error::plain(format!(
        "tube algebra supports at most {MAX_TUBE} annular basis elements"
    ))
}

impl Tube {
    pub(crate) fn build(
        model: &mut Model,
        double: &Double,
        category: &Category,
        label: Label,
    ) -> Result<Self> {
        let mut builder = Builder::new(model, double, category)?;
        let forms = builder.forms()?;
        let (blocks, forward, vacuum) = builder.fourier(&forms)?;
        let inverse = annular_inverse(&mut builder.field, &forms, &blocks, &forward)?;
        let unit_dense = builder.geometric_unit(&forms, vacuum)?;
        let unit = Element::from_dense(unit_dense);
        let result = Self {
            category: label,
            conductor: builder.ambient.conductor,
            forms,
            blocks,
            forward,
            inverse,
            unit,
            boundary_simples: category.count(),
        };
        let actual = result.transform(&mut builder.field, &result.forward, &result.unit)?;
        let mut expected = vec![Rational::zero(); result.forms.len()];
        for b in &result.blocks {
            for i in 0..b.size() {
                expected[b.index(i, i)] = Rational::integer(1);
            }
        }
        if actual != Element::from_dense(expected) {
            return Err(Error::plain(
                "geometric annular unit does not act by identity blocks",
            ));
        }
        builder.field.model.allocate(
            16 + result.forms.len() * 8
                + result
                    .blocks
                    .iter()
                    .map(|b| 2 * b.size() + 4)
                    .sum::<usize>()
                + result.forward.cells()
                + result.inverse.cells()
                + result.unit.cells(),
        )?;
        Ok(result)
    }

    fn transform(
        &self,
        field: &mut Field<'_>,
        matrix: &Matrix,
        input: &Element,
    ) -> Result<Element> {
        let mut values = vec![Rational::zero(); self.forms.len()];
        field.model.spend(values.len() * input.0.len())?;
        for &(c, ref value) in &input.0 {
            for (r, v) in values.iter_mut().enumerate() {
                let x = matrix.get(r, c);
                if !x.is_zero() {
                    let term = field.multiply(x, value)?;
                    *v = field.add(v, &term)?;
                }
            }
        }
        Ok(Element::from_dense(values))
    }

    pub(crate) fn block(&self, a: usize) -> Result<usize> {
        Ok(self
            .blocks
            .get(a)
            .ok_or_else(|| Error::plain("tube bulk block index is out of range"))?
            .size())
    }

    fn form(&self, index: usize) -> Result<Form> {
        self.forms
            .get(index)
            .copied()
            .ok_or_else(|| Error::plain("tube annular basis index is out of range"))
    }

    pub(crate) fn basis(&self, index: usize) -> Result<Element> {
        self.form(index)?;
        Ok(Element(vec![(index, Rational::integer(1))]))
    }

    pub(crate) fn unit(&self) -> Element {
        self.unit.clone()
    }

    pub(crate) fn central(&self, model: &mut Model, a: usize) -> Result<Element> {
        let size = self.block(a)?;
        let b = &self.blocks[a];
        let entries = (0..size)
            .map(|i| (b.index(i, i), Rational::integer(1)))
            .collect();
        self.transform(
            &mut Field::new(model, self.conductor),
            &self.inverse,
            &Element(entries),
        )
    }

    pub(crate) fn matrix_unit(&self, model: &mut Model, [a, r, c]: [usize; 3]) -> Result<Element> {
        let i = self.index(a, r, c)?;
        self.transform(
            &mut Field::new(model, self.conductor),
            &self.inverse,
            &Element(vec![(i, Rational::integer(1))]),
        )
    }

    fn index(&self, a: usize, r: usize, c: usize) -> Result<usize> {
        let size = self.block(a)?;
        if r >= size || c >= size {
            return Err(Error::plain(
                "tube matrix coordinate is outside its bulk block",
            ));
        }
        Ok(self.blocks[a].index(r, c))
    }

    pub(crate) fn coefficient(&self, x: &Element, i: usize) -> Result<Rational> {
        self.form(i)?;
        Ok(x.coefficient(i))
    }

    pub(crate) fn block_entry(
        &self,
        model: &mut Model,
        x: &Element,
        [a, r, c]: [usize; 3],
    ) -> Result<Rational> {
        let i = self.index(a, r, c)?;
        Ok(self
            .transform(&mut Field::new(model, self.conductor), &self.forward, x)?
            .coefficient(i))
    }

    pub(crate) fn add(&self, model: &mut Model, a: &Element, b: &Element) -> Result<Element> {
        let mut field = Field::new(model, self.conductor);
        let mut result = vec![Rational::zero(); self.forms.len()];
        for &(i, ref value) in a.0.iter().chain(&b.0) {
            result[i] = field.add(&result[i], value)?;
        }
        Ok(Element::from_dense(result))
    }

    pub(crate) fn scale(
        &self,
        model: &mut Model,
        scalar: &Rational,
        value: &Element,
    ) -> Result<Element> {
        let mut field = Field::new(model, self.conductor);
        let mut coefficients = Vec::with_capacity(value.0.len());
        for &(i, ref v) in &value.0 {
            let p = field.multiply(scalar, v)?;
            if !p.is_zero() {
                coefficients.push((i, p));
            }
        }
        Ok(Element(coefficients))
    }

    pub(crate) fn product(&self, model: &mut Model, a: &Element, b: &Element) -> Result<Element> {
        let mut field = Field::new(model, self.conductor);
        let a = self.transform(&mut field, &self.forward, a)?;
        let b = self.transform(&mut field, &self.forward, b)?;
        let mut result = vec![Rational::zero(); self.forms.len()];
        for block in &self.blocks {
            let n = block.size();
            let mut a0 = Matrix::zero(n, n)?;
            let mut b0 = Matrix::zero(n, n)?;
            for r in 0..n {
                for c in 0..n {
                    let j = block.index(r, c);
                    a0.set(r, c, a.coefficient(j));
                    b0.set(r, c, b.coefficient(j));
                }
            }
            let product = a0.multiply(&mut field, &b0)?;
            for r in 0..n {
                for c in 0..n {
                    result[block.index(r, c)] = product.get(r, c).clone();
                }
            }
        }
        self.transform(&mut field, &self.inverse, &Element::from_dense(result))
    }

    pub(crate) fn stack(
        &self,
        model: &mut Model,
        double: &Double,
        category: &Category,
        key: [usize; 2],
    ) -> Result<String> {
        let [i, j] = key;
        let f = self.form(i)?;
        let g = self.form(j)?;
        let product = self.product(model, &self.basis(i)?, &self.basis(j)?)?;
        Builder::new(model, double, category)?.check_stack(f, g, &product, &self.forms)?;
        Ok(format!(
            "rvDoubleTubeStack[{i},{j}] = terms={}, verified=true",
            product.0.len()
        ))
    }

    pub(crate) fn data(&self, name: &str) -> String {
        let mut counts = BTreeMap::new();
        for block in &self.blocks {
            *counts.entry(block.size()).or_insert(0) += 1;
        }
        format!("rvDoubleTubeData[{name}] = dimension={}, center={}, boundary_simples={}, conductor={}, blocks={counts:?}, fourier_verified=true",
            self.forms.len(),self.blocks.len(),self.boundary_simples,self.conductor)
    }

    pub(crate) fn annular(&self, name: &str, [j, i]: [usize; 2]) -> Result<String> {
        if j >= self.boundary_simples || i >= self.boundary_simples {
            return Err(Error::plain(
                "tube annular endpoints are outside the boundary category",
            ));
        }
        let dimension = self.forms.iter().filter(|f| f.j == j && f.i == i).count();
        Ok(format!(
            "rvDoubleTubeAnnular[{name},{j},{i}] = dimension={dimension}"
        ))
    }

    pub(crate) fn basis_data(&self, name: &str, b: usize) -> Result<String> {
        let Form { j, k, i, c, u, v } = self.form(b)?;
        Ok(format!("rvDoubleTubeBasisData[{name},{b}] = left={j}, loop={k}, right={i}, channel={c}, input_multiplicity={u}, output_multiplicity={v}"))
    }

    pub(crate) fn block_data(&self, name: &str, a: usize) -> Result<String> {
        let d = self.block(a)?;
        Ok(format!(
            "rvDoubleTubeBlockData[{name},{a}] = dimension={d}, channels={:?}",
            self.blocks[a].slots
        ))
    }

    pub(crate) fn decompose(&self, model: &mut Model, x: &Element) -> Result<String> {
        let mut field = Field::new(model, self.conductor);
        let coordinates = self.transform(&mut field, &self.forward, x)?;
        if self.transform(&mut field, &self.inverse, &coordinates)? != *x {
            return Err(Error::plain(
                "tube Fourier and inverse failed to recover the annular coordinates",
            ));
        }
        let mut out = Vec::new();
        for (a, b) in self.blocks.iter().enumerate() {
            for r in 0..b.size() {
                for c in 0..b.size() {
                    let v = coordinates.coefficient(b.index(r, c));
                    if !v.is_zero() {
                        out.push(format!("({a},{r},{c},{:?}/{})", v.numerator, v.denominator));
                    }
                }
            }
        }
        Ok(format!(
            "rvDoubleTubeDecompose = conductor={}, blocks=[{}]",
            self.conductor,
            out.join(", ")
        ))
    }
}

impl Element {
    fn from_dense(values: Vec<Rational>) -> Self {
        Self(
            values
                .into_iter()
                .enumerate()
                .filter(|(_, v)| !v.is_zero())
                .collect(),
        )
    }

    fn coefficient(&self, i: usize) -> Rational {
        self.0
            .binary_search_by_key(&i, |(j, _)| *j)
            .ok()
            .map(|j| self.0[j].1.clone())
            .unwrap_or_else(Rational::zero)
    }

    fn cells(&self) -> usize {
        4 + self
            .0
            .iter()
            .map(|(_, x)| 4 + 2 * x.numerator.len())
            .sum::<usize>()
    }

    pub(crate) fn persist(self, model: &mut Model) -> Result<Value> {
        model.allocate(self.cells())?;
        Ok(Value::DoubleTube(Rc::new(self)))
    }

    pub(crate) fn format(&self) -> String {
        self.0
            .iter()
            .map(|(i, x)| format!("({i},{:?}/{})", x.numerator, x.denominator))
            .collect::<Vec<_>>()
            .join(", ")
    }
}

fn annular_inverse(
    field: &mut Field<'_>,
    forms: &[Form],
    blocks: &[Block],
    forward: &Matrix,
) -> Result<Matrix> {
    type Locations = BTreeMap<[usize; 2], Vec<usize>>;
    let mut domain = Locations::new();
    for (i, f) in forms.iter().enumerate() {
        domain.entry([f.j, f.i]).or_default().push(i);
    }
    let mut codomain = Locations::new();
    for b in blocks {
        for (i, &[left, _]) in b.slots.iter().enumerate() {
            for (j, &[right, _]) in b.slots.iter().enumerate() {
                codomain
                    .entry([left, right])
                    .or_default()
                    .push(b.index(i, j));
            }
        }
    }
    let n = forms.len();
    let mut inverse = Matrix::zero(n, n)?;
    if !domain.keys().eq(codomain.keys()) {
        return Err(Error::plain(
            "tube Fourier transform disagrees with annular support",
        ));
    }
    for (key, cols) in domain {
        let rows = &codomain[&key];
        if rows.len() != cols.len() {
            return Err(Error::plain(
                "tube Fourier transform disagrees with an annular Hom dimension",
            ));
        }
        let mut block = Matrix::zero(rows.len(), cols.len())?;
        for (i, &r) in rows.iter().enumerate() {
            for (j, &c) in cols.iter().enumerate() {
                block.set(i, j, forward.get(r, c).clone());
            }
        }
        let block = block.inverse(field)?;
        for (j, &c) in cols.iter().enumerate() {
            for (i, &r) in rows.iter().enumerate() {
                inverse.set(c, r, block.get(j, i).clone());
            }
        }
    }
    let identity = matrix::Matrix::identity(n)?;
    if forward.multiply(field, &inverse)? != identity
        || inverse.multiply(field, forward)? != identity
    {
        return Err(Error::plain(
            "tube annular Fourier transform is not invertible",
        ));
    }
    Ok(inverse)
}
