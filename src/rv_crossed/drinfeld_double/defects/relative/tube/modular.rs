use crate::error::{Error, Result};
use crate::model::Model;
use crate::rv_lens::{self, Rational, Step};
use crate::syntax::{Label, ParsedExpr};

use super::super::central::{Bulk, Central};
use super::super::linear::{self, GradedMap};
use super::super::tensor::Layout;
use super::super::{decompose, Ambient, Image};
use super::matrix::{Field, Matrix};
use super::{Category, Double, Element, Tube};

mod rotation;
pub(crate) use rotation::Rotation;

pub(crate) enum Command {
    Data(Label),
    Period(Label, usize),
    RotationData(Box<ParsedExpr>),
}

pub(crate) enum Expression {
    SEntry(Label, [usize; 2]),
    Twist(Label, usize),
    Lens(Label, i128),
    Indicator(Label, [usize; 2]),
    Rotate(Label, [usize; 2]),
    RotationEntry([usize; 2], Box<ParsedExpr>),
    RotationTrace(Box<ParsedExpr>),
    Word(Vec<Step>, Box<ParsedExpr>),
}

pub(crate) struct Torus {
    conductor: usize,
    s: Matrix,
    twists: Vec<usize>,
    dimensions: Vec<usize>,
    order: usize,
    exponent: usize,
}

impl Tube {
    pub(crate) fn torus<'a>(
        &'a self,
        model: &mut Model,
        double: &Double,
        category: &Category,
    ) -> Result<&'a Torus> {
        if self.modular.get().is_none() {
            let t = Torus::build(model, double, category)?;
            if t.conductor != self.conductor || t.s.rows != self.blocks.len() {
                return Err(Error::plain(
                    "boundary modular data disagree with the tube center",
                ));
            }
            model.allocate(t.s.cells() + 2 * t.twists.len() + 2 * t.dimensions.len() + 12)?;
            self.modular.set(t).map_err(|_| {
                Error::plain("could not store the verified boundary modular action")
            })?;
        }
        Ok(self.modular.get().expect("verified boundary modular data"))
    }

    fn center(&self, field: &mut Field<'_>, x: &Element) -> Result<Vec<Rational>> {
        let entries = self.transform(field, &self.forward, x)?;
        let mut values = Vec::with_capacity(self.blocks.len());
        for b in &self.blocks {
            let lead = entries.coefficient(b.start);
            for r in 0..b.size() {
                for c in 0..b.size() {
                    let value = entries.coefficient(b.index(r, c));
                    if if r == c {
                        value != lead
                    } else {
                        !value.is_zero()
                    } {
                        return Err(Error::plain(
                            "modular action requires an element of the tube center",
                        ));
                    }
                }
            }
            values.push(lead);
        }
        Ok(values)
    }

    pub(crate) fn extract_center(&self, model: &mut Model, x: &Element) -> Result<Vec<Rational>> {
        self.center(&mut Field::new(model, self.conductor), x)
    }

    pub(crate) fn lift_center(&self, model: &mut Model, values: &[Rational]) -> Result<Element> {
        self.central_element(&mut Field::new(model, self.conductor), values)
    }

    fn central_element(&self, field: &mut Field<'_>, values: &[Rational]) -> Result<Element> {
        if values.len() != self.blocks.len() {
            return Err(Error::plain(
                "tube center coordinate count disagrees with its bulk colors",
            ));
        }
        let mut result = vec![Rational::zero(); self.forms.len()];
        for (x, b) in values.iter().zip(&self.blocks) {
            for r in 0..b.size() {
                result[b.index(r, r)] = x.clone();
            }
        }
        self.transform(field, &self.inverse, &Element::from_dense(result))
    }

    pub(crate) fn modular_word(
        &self,
        model: &mut Model,
        torus: &Torus,
        x: &Element,
        steps: &[Step],
    ) -> Result<Element> {
        let mut field = Field::new(model, self.conductor);
        let values = self.center(&mut field, x)?;
        let values = torus.act(&mut field, values, steps)?;
        self.central_element(&mut field, &values)
    }

    pub(crate) fn indicator(
        &self,
        model: &mut Model,
        torus: &Torus,
        [i, n]: [usize; 2],
    ) -> Result<Rational> {
        if n == 0 {
            return Err(Error::plain(
                "higher Frobenius–Schur indicator index must be positive",
            ));
        }
        let multiplicities = self.branches(i)?;
        let mut field = Field::new(model, self.conductor);
        let mut numerator = Rational::zero();
        for (a, m) in multiplicities.into_iter().enumerate() {
            if m == 0 {
                continue;
            }
            let p = torus.power(a, n as i128);
            let root = field.phase(p)?;
            let w = Rational::integer((m * torus.dimensions[a]) as i128);
            let contribution = field.multiply(&root, &w)?;
            numerator = field.add(&numerator, &contribution)?;
        }
        let unit = rv_lens::fraction(field.model, &[1], torus.order as u128)?;
        field.multiply(&numerator, &unit)
    }

    pub(crate) fn period(&self, torus: &Torus, i: usize) -> Result<String> {
        let mut period = 1;
        for (a, m) in self.branches(i)?.into_iter().enumerate() {
            if m != 0 {
                period = lcm(period, torus.root_order(torus.twists[a]));
            }
        }
        Ok(format!(
            "rvDoubleTubeIndicatorPeriod[{i}] = period={period}, exponent={}",
            torus.exponent
        ))
    }

    fn branches(&self, i: usize) -> Result<Vec<usize>> {
        if i >= self.boundary_simples {
            return Err(Error::plain(
                "indicator simple is outside the boundary category",
            ));
        }
        Ok(self
            .blocks
            .iter()
            .map(|b| b.slots.iter().filter(|x| x[0] == i).count())
            .collect())
    }
}

impl Torus {
    fn build(model: &mut Model, double: &Double, category: &Category) -> Result<Self> {
        let ambient = Ambient::new(model, double, &[category], [])?;
        let order = ambient.order();
        let h = category.boundaries[0].members()?.len();
        let modules = (0..double.count())
            .map(|i| Bulk::new(model, &ambient, category, i))
            .collect::<Result<Vec<_>>>()?;
        let mut field = Field::new(model, ambient.conductor);
        let mut dimensions = Vec::with_capacity(modules.len());
        let mut twists = Vec::with_capacity(modules.len());
        for (a, bulk) in modules.iter().enumerate() {
            let dim = bulk.module.dimension() / h;
            if dim != double.simples[a].dimension {
                return Err(Error::plain(
                    "boundary spherical dimension differs from bulk",
                ));
            }
            let phase = double.twists[a] * (ambient.conductor / double.conductor);
            let theta = checked_twist(&mut field, &ambient, category, bulk)?;
            if theta != field.phase(phase)? {
                return Err(Error::plain(
                    "boundary spherical twist disagrees with the chosen bulk orientation",
                ));
            }
            dimensions.push(dim);
            twists.push(phase);
        }
        let scale = ambient.conductor / double.conductor;
        let mut s = Matrix::zero(modules.len(), modules.len())?;
        let inv = rv_lens::fraction(field.model, &[1], (h * order) as u128)?;
        for (a, aa) in modules.iter().enumerate() {
            for (b, bb) in modules.iter().enumerate().skip(a) {
                let ab = Central::new(field.model, &ambient, [aa, aa], &bb.module)?;
                let ba = Central::new(field.model, &ambient, [bb, bb], &aa.module)?;
                if ab.products[0].bases != ba.products[1].bases
                    || ab.products[1].bases != ba.products[0].bases
                {
                    return Err(Error::plain(
                        "boundary double crossing does not return to its initial balanced basis",
                    ));
                }
                // The forward central moves use pi^{-1}; this is the
                // conjugated modular orientation used by the bulk table.
                let rows = ab.maps[0]
                    .0
                    .iter()
                    .zip(&ba.maps[0].0)
                    .map(|(l, r)| {
                        l.iter()
                            .map(|x| {
                                let y = r[x.target];
                                Image {
                                    target: y.target,
                                    phase: ambient.sum([x.phase, y.phase]),
                                }
                            })
                            .collect()
                    })
                    .collect::<Vec<Vec<_>>>();
                let raw = trace(&mut field, &rows)?;
                let coefficient = field.multiply(&raw, &inv)?;
                let direct = double.s[a][b];
                let phase = field.phase(direct.phase * scale)?;
                let weight = rv_lens::fraction(field.model, &[direct.weight], order as u128)?;
                if coefficient != field.multiply(&phase, &weight)? {
                    return Err(Error::plain(
                        "boundary spherical S disagrees with the chosen bulk orientation",
                    ));
                }
                if a == 0
                    && coefficient
                        != rv_lens::fraction(field.model, &[dimensions[b] as i128], order as u128)?
                {
                    return Err(Error::plain(
                        "boundary spherical S has the wrong vacuum column",
                    ));
                }
                s.set(a, b, coefficient.clone());
                s.set(b, a, coefficient);
            }
        }
        let mut c = Matrix::zero(s.rows, s.cols)?;
        let mut st = s.clone();
        for (j, &d) in double.duals.iter().enumerate() {
            c.set(d, j, Rational::integer(1));
            for i in 0..s.rows {
                st.set(i, j, field.shift(s.get(i, j), twists[j])?);
            }
        }
        let st2 = st.multiply(&mut field, &st)?;
        if s.multiply(&mut field, &s)? != c || st2.multiply(&mut field, &st)? != c {
            return Err(Error::plain(
                "boundary modular operators fail S² = (ST)³ = C",
            ));
        }
        let mut result = Self {
            conductor: ambient.conductor,
            s,
            twists,
            dimensions,
            order,
            exponent: 1,
        };
        result.exponent = result
            .twists
            .iter()
            .fold(1, |m, &p| lcm(m, result.root_order(p)));
        Ok(result)
    }

    fn root_order(&self, phase: usize) -> usize {
        self.conductor / (rv_lens::gcd(self.conductor as u128, phase as u128) as usize)
    }

    fn power(&self, a: usize, n: i128) -> usize {
        (n.rem_euclid(self.conductor as i128) as usize * self.twists[a]) % self.conductor
    }

    pub(crate) fn data(&self, name: &str) -> String {
        format!("rvDoubleTubeModular[{name}] = center={}, conductor={}, exponent={}, orientation=conjugated_double_crossing, verified=true",self.s.rows,self.conductor,self.exponent)
    }

    pub(crate) fn s_entry(&self, key: [usize; 2]) -> Result<Rational> {
        self.s.entry(key)
    }

    pub(crate) fn twist(&self, model: &mut Model, a: usize) -> Result<Rational> {
        let &p = self
            .twists
            .get(a)
            .ok_or_else(|| Error::plain("modular simple index is out of range"))?;
        Field::new(model, self.conductor).phase(p)
    }

    pub(crate) fn lens(&self, model: &mut Model, power: i128) -> Result<Rational> {
        let mut input = vec![Rational::zero(); self.s.rows];
        input[0] = Rational::integer(1);
        let mut field = Field::new(model, self.conductor);
        Ok(self.act(&mut field, input, &[Step::S, Step::T(power), Step::S])?[0].clone())
    }

    fn act(
        &self,
        field: &mut Field<'_>,
        mut values: Vec<Rational>,
        steps: &[Step],
    ) -> Result<Vec<Rational>> {
        for &step in steps {
            match step {
                Step::S => {
                    let mut input = Matrix::zero(values.len(), 1)?;
                    for (i, x) in values.into_iter().enumerate() {
                        input.set(i, 0, x);
                    }
                    let output = self.s.multiply(field, &input)?;
                    values = (0..output.rows).map(|i| output.get(i, 0).clone()).collect();
                }
                Step::T(power) => {
                    for (a, x) in values.iter_mut().enumerate() {
                        *x = field.shift(x, self.power(a, power))?;
                    }
                }
            }
        }
        Ok(values)
    }
}

fn checked_twist(
    field: &mut Field<'_>,
    ambient: &Ambient<'_>,
    category: &Category,
    bulk: &Bulk,
) -> Result<Rational> {
    let self_cross = Central::new(field.model, ambient, [bulk, bulk], &bulk.module)?;
    let [source, target] = &self_cross.products;
    if source.bases != target.bases {
        return Err(Error::plain(
            "spherical twist requires identical graded source and target",
        ));
    }
    let actual = GradedMap::monomial(field, &self_cross.maps[1].0)?;
    let raw = trace(field, &self_cross.maps[1].0)?;
    let inv = rv_lens::fraction(field.model, &[1], bulk.module.dimension() as u128)?;
    let theta = field.multiply(&raw, &inv)?;
    let identity = GradedMap::identity(&bulk.module)?;
    let layout = Layout {
        product: source,
        factors: [&bulk.module, &bulk.module],
    };
    let homs = decompose(field.model, ambient, category, &bulk.module)?;
    for (i, hom) in homs
        .iter()
        .enumerate()
        .filter(|(_, h)| !h.orbits.is_empty())
    {
        let simple = ambient.induce(field.model, category, i)?;
        let g = category.simple(i)?.simple.flux;
        let maps = linear::hom_basis(field, ambient, &simple, &bulk.module, g, hom)?;
        for p in &maps {
            for q in &maps {
                let probe = p.inclusion.compose(field, &q.projection)?;
                let tensor = linear::tensor(field, ambient, layout, layout, [&probe, &identity])?;
                let observed = tensor.inner(field, &actual)?;
                let reference = probe.inner(field, &identity)?;
                if observed != field.multiply(&theta, &reference)? {
                    return Err(Error::plain(
                        "spherical partial trace of the boundary crossing is not scalar",
                    ));
                }
            }
        }
    }
    Ok(theta)
}

fn trace(field: &mut Field<'_>, rows: &[Vec<Image>]) -> Result<Rational> {
    let mut sum = Rational::zero();
    for row in rows {
        field.model.spend(row.len())?;
        for (i, x) in row.iter().enumerate() {
            if i == x.target {
                let root = field.phase(x.phase)?;
                sum = field.add(&sum, &root)?;
            }
        }
    }
    Ok(sum)
}

fn lcm(a: usize, b: usize) -> usize {
    a / (rv_lens::gcd(a as u128, b as u128) as usize) * b
}
