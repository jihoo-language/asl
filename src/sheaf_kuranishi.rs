use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::{Model, Value};
use crate::odd_cobar::Prime;
use crate::prime_linear::{add_scaled, Basis, Coordinate};
use crate::sheaf::{self, Cochain, Cochains, Support, Term};
use crate::sheaf_deform::{Boundaries, Jet};
use crate::sheaf_ext::{self, Hom};
use crate::syntax::Label;

pub(crate) enum Query {
    Coefficient(Label, Vec<usize>, Coefficient),
    Arc(Label, Vec<i64>, Vec<usize>, usize, bool),
}

#[derive(Clone, Copy)]
pub(crate) enum Coefficient {
    Lift,
    Residual,
    Obstruction,
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct Monomial(Box<[(usize, usize)]>);

impl Monomial {
    fn product(&self, model: &mut Model, other: &Self) -> Result<Self> {
        model.spend(self.0.len() + other.0.len())?;
        let mut indices = Vec::with_capacity(self.0.len() + other.0.len());
        let mut rhs = other.0.iter().copied().peekable();
        for &(i, power) in &self.0 {
            while rhs.peek().is_some_and(|&(j, _)| j < i) {
                indices.push(rhs.next().unwrap());
            }
            let power = if rhs.peek().is_some_and(|&(j, _)| j == i) {
                power + rhs.next().unwrap().1
            } else {
                power
            };
            indices.push((i, power));
        }
        indices.extend(rhs);
        Ok(Self(indices.into()))
    }

    fn cells(&self) -> usize {
        4 + 2 * self.0.len()
    }

    pub(crate) fn factors(&self) -> &[(usize, usize)] {
        &self.0
    }
}

type Layer<T> = BTreeMap<Monomial, T>;
type Polynomial<T> = BTreeMap<usize, Layer<T>>;
type CPol = Polynomial<Support>;
type HPol = Polynomial<Rc<[Coordinate]>>;

#[derive(Default)]
struct CochainLayer {
    rows: Layer<Cochain>,
    cells: usize,
}

impl CochainLayer {
    fn add(
        &mut self,
        model: &mut Model,
        p: Prime,
        key: Monomial,
        raw: Vec<Term>,
        factor: u32,
    ) -> Result<()> {
        use std::collections::btree_map::Entry;
        let key_cells = key.cells();
        model.spend(key_cells * (1 + self.rows.len().max(1).ilog2() as usize))?;
        match self.rows.entry(key) {
            Entry::Vacant(entry) => {
                let mut cochain = Cochain::default();
                for (g, c) in raw {
                    cochain.add(model, p, g, p.mul(c, factor))?;
                }
                if cochain.cells() != 0 {
                    self.cells += key_cells + cochain.cells();
                    sheaf::check_cells(self.cells)?;
                    entry.insert(cochain);
                }
            }
            Entry::Occupied(mut entry) => {
                let cochain = entry.get_mut();
                self.cells -= cochain.cells();
                for (g, c) in raw {
                    cochain.add(model, p, g, p.mul(c, factor))?;
                }
                self.cells += cochain.cells();
                sheaf::check_cells(self.cells)?;
                if cochain.cells() == 0 {
                    self.cells -= key_cells;
                    entry.remove();
                }
            }
        }
        Ok(())
    }

    fn into_terms(self) -> Layer<Vec<Term>> {
        self.rows
            .into_iter()
            .map(|(m, c)| (m, c.into_terms()))
            .collect()
    }
}

type Factor<'a> = (&'a CPol, usize);

fn multiply(
    model: &mut Model,
    hom: &Hom,
    after: Factor<'_>,
    before: Factor<'_>,
    degree: usize,
    factor: u32,
    result: &mut CochainLayer,
) -> Result<()> {
    let (after, m) = after;
    let (before, n) = before;
    for (&e, outer) in after.range(..degree) {
        model.spend(1)?;
        let Some(inner) = before.get(&(degree - e)) else {
            continue;
        };
        for (x, a) in outer {
            model.spend(inner.len())?;
            for (y, b) in inner {
                let raw = sheaf_ext::compose(model, (hom, m, a), (hom, n, b))?;
                if !raw.is_empty() {
                    let key = x.product(model, y)?;
                    result.add(model, hom.prime(), key, raw, factor)?;
                }
            }
        }
    }
    Ok(())
}

fn add_raw(
    model: &mut Model,
    p: Prime,
    sum: &mut Cochain,
    raw: &[Term],
    factor: u32,
) -> Result<()> {
    for (g, c) in raw {
        sum.add(model, p, g.clone(), p.mul(*c, factor))?;
    }
    Ok(())
}

fn combine(
    model: &mut Model,
    p: Prime,
    left: &[Term],
    right: &[Term],
    factor: u32,
) -> Result<Vec<Term>> {
    if right.is_empty() {
        model.spend(left.len())?;
        return Ok(left.to_vec());
    }
    let mut sum = Cochain::default();
    add_raw(model, p, &mut sum, left, 1)?;
    add_raw(model, p, &mut sum, right, factor)?;
    Ok(sum.into_terms())
}

#[derive(Default)]
struct Contraction {
    incoming2: Option<Boundaries>,
    incoming3: Option<Boundaries>,
}

impl Contraction {
    fn boundary_preimage(
        &mut self,
        model: &mut Model,
        hom: &Hom,
        degree: usize,
        raw: &[Term],
        factor: u32,
    ) -> Result<Vec<Term>> {
        if raw.is_empty() {
            return Ok(Vec::new());
        }
        let slot = if degree == 2 {
            &mut self.incoming2
        } else {
            &mut self.incoming3
        };
        if slot.is_none() {
            *slot = Some(Boundaries::build(model, hom, degree)?);
        }
        let cells = self.incoming2.as_ref().map_or(0, Boundaries::cells)
            + self.incoming3.as_ref().map_or(0, Boundaries::cells);
        sheaf::check_cells(cells)?;
        let reducer = if degree == 2 {
            &self.incoming2
        } else {
            &self.incoming3
        };
        reducer
            .as_ref()
            .expect("incoming reducer built")
            .preimage(model, hom, raw, factor)?
            .ok_or_else(|| Error::plain("Kuranishi contraction could not solve an exact cochain"))
    }

    fn step(&mut self, model: &mut Model, hom: &Hom, q: &[Term], d_q: &[Term]) -> Result<Step> {
        let p = hom.prime();
        let neg = p.negate(1);
        // Remove the noncycle part using the fixed linear preimage of dQ.
        // That preimage lies in the chosen complement of the cycles.
        let noncycle = self.boundary_preimage(model, hom, 3, d_q, 1)?;
        let closed = combine(model, p, q, &noncycle, neg)?;
        let class = hom.cohomology.class_coordinates(model, hom, 2, &closed)?;
        let representative = hom
            .cohomology
            .representative_support(model, hom, 2, &class)?;
        let exact = combine(model, p, &closed, &representative, neg)?;
        let lift = self.boundary_preimage(model, hom, 2, &exact, neg)?;
        let derivative = hom.differential_support(model, &lift)?;
        let residual = combine(model, p, q, &derivative, 1)?;
        let expected = combine(model, p, &representative, &noncycle, 1)?;
        if !sheaf::compare(model, &residual, &expected)?.is_eq() {
            return Err(Error::plain(
                "Kuranishi cochain contraction identity failed",
            ));
        }
        Ok(Step {
            lift,
            residual,
            class,
        })
    }
}

struct Step {
    lift: Vec<Term>,
    residual: Vec<Term>,
    class: Vec<Coordinate>,
}

struct DegreeCoefficients {
    lift: Layer<Support>,
    residual: Layer<Support>,
    class: Layer<Rc<[Coordinate]>>,
}

pub(crate) struct Chart {
    pub(crate) sheaf: Rc<str>,
    pub(crate) order: usize,
    pub(crate) directions: usize,
    pub(crate) obstructions: usize,
    a: CPol,
    residual: CPol,
    kappa: HPol,
}

impl Chart {
    pub(crate) fn build(
        model: &mut Model,
        hom: &Hom,
        sheaf: Rc<str>,
        order: usize,
        selected: Option<Vec<Rc<[Coordinate]>>>,
    ) -> Result<Self> {
        let p = hom.prime();
        let tangent_rank = hom.cohomology.rank(model, hom, 1)?;
        let obstructions = hom.cohomology.rank(model, hom, 2)?;
        let mut residues = BTreeMap::new();
        let directions = selected.as_ref().map_or(tangent_rank, Vec::len);
        let mut independent = Basis::default();
        let mut cells = directions + 1;
        sheaf::check_cells(cells)?;
        model.spend(directions)?;
        for i in 0..directions {
            let fallback = [(i, 1)];
            let class = selected.as_ref().map_or(&fallback[..], |values| &values[i]);
            if selected.is_some() {
                let old = independent.rank();
                independent.insert(model, p, class.to_vec())?;
                if independent.rank() == old {
                    return Err(Error::plain(
                        "Kuranishi tangent directions must be linearly independent",
                    ));
                }
            }
            if order != 0 {
                let representative = hom
                    .cohomology
                    .representative_support(model, hom, 1, class)?;
                let monomial = Monomial(vec![(i, 1)].into());
                check_support(&mut cells, &monomial, &representative)?;
                residues.insert(monomial, representative.into());
            }
        }

        let mut chart = Self {
            sheaf,
            order,
            directions,
            obstructions,
            a: BTreeMap::new(),
            residual: BTreeMap::new(),
            kappa: BTreeMap::new(),
        };
        let mut pending = BTreeSet::new();
        let mut contraction = Contraction::default();
        if !residues.is_empty() {
            chart.a.insert(1, residues);
            if order >= 2 {
                pending.insert(2);
            }
        }
        while let Some(degree) = pending.pop_first() {
            let DegreeCoefficients {
                lift,
                residual,
                class,
            } = chart.degree(model, hom, &mut contraction, degree, &mut cells)?;
            let a_changed = !lift.is_empty();
            let r_changed = !residual.is_empty();
            if a_changed {
                chart.a.insert(degree, lift);
            }
            if r_changed {
                chart.residual.insert(degree, residual);
            }
            if !class.is_empty() {
                chart.kappa.insert(degree, class);
            }
            // Later possible terms come from A*A, or A*R in the Bianchi identity.
            if a_changed || r_changed {
                schedule(model, &mut pending, order, degree, chart.a.keys().copied())?;
            }
            if a_changed {
                schedule(
                    model,
                    &mut pending,
                    order,
                    degree,
                    chart.residual.keys().copied(),
                )?;
            }
        }
        model.allocate(cells)?;
        Ok(chart)
    }

    fn degree(
        &self,
        model: &mut Model,
        hom: &Hom,
        contraction: &mut Contraction,
        degree: usize,
        cells: &mut usize,
    ) -> Result<DegreeCoefficients> {
        let neg = hom.prime().negate(1);
        let mut quadratic = CochainLayer::default();
        multiply(
            model,
            hom,
            (&self.a, 1),
            (&self.a, 1),
            degree,
            neg,
            &mut quadratic,
        )?;
        let mut bianchi = CochainLayer::default();
        multiply(
            model,
            hom,
            (&self.a, 1),
            (&self.residual, 2),
            degree,
            1,
            &mut bianchi,
        )?;
        multiply(
            model,
            hom,
            (&self.residual, 2),
            (&self.a, 1),
            degree,
            neg,
            &mut bianchi,
        )?;
        sheaf::check_cells(quadratic.cells + bianchi.cells + *cells)?;
        let mut bianchi = bianchi.into_terms();
        let mut a = BTreeMap::new();
        let mut r = BTreeMap::new();
        let mut k = BTreeMap::new();
        for (mono, q) in quadratic.into_terms() {
            let differential = hom.differential_support(model, &q)?;
            let expected = bianchi.remove(&mono).unwrap_or_default();
            if !sheaf::compare(model, &differential, &expected)?.is_eq() {
                return Err(Error::plain("Kuranishi Bianchi identity failed"));
            }
            let Step {
                lift,
                residual,
                class,
            } = contraction.step(model, hom, &q, &differential)?;
            if !lift.is_empty() {
                check_support(cells, &mono, &lift)?;
                a.insert(mono.clone(), lift.into());
            }
            if !residual.is_empty() {
                check_support(cells, &mono, &residual)?;
                r.insert(mono.clone(), residual.into());
            }
            if !class.is_empty() {
                *cells += mono.cells() + 2 * class.len();
                sheaf::check_cells(*cells)?;
                k.insert(mono, class.into());
            }
        }
        if !bianchi.is_empty() {
            return Err(Error::plain(
                "Kuranishi Bianchi identity failed on an absent quadratic term",
            ));
        }
        Ok(DegreeCoefficients {
            lift: a,
            residual: r,
            class: k,
        })
    }

    pub(crate) fn coefficient(
        &self,
        model: &mut Model,
        powers: &[usize],
        coefficient: Coefficient,
    ) -> Result<Value> {
        if powers.len() != self.directions {
            return Err(Error::plain(
                "Kuranishi exponent list has the wrong number of directions",
            ));
        }
        model.spend(powers.len())?;
        let mut degree = 0usize;
        let mut indices = Vec::new();
        for (i, &exponent) in powers.iter().enumerate() {
            if exponent != 0 {
                degree = degree.checked_add(exponent).ok_or_else(outside)?;
                indices.push((i, exponent));
            }
        }
        self.check_order(degree)?;
        let key = Monomial(indices.into());
        Ok(match coefficient {
            Coefficient::Lift => Value::SheafC(lookup(&self.a, degree, &key)),
            Coefficient::Residual => Value::SheafC(lookup(&self.residual, degree, &key)),
            Coefficient::Obstruction => Value::OddClass(lookup(&self.kappa, degree, &key)),
        })
    }

    pub(crate) fn equation(
        &self,
        model: &mut Model,
        index: usize,
    ) -> Result<impl Iterator<Item = (&Monomial, u32)> + '_> {
        if index >= self.obstructions {
            return Err(Error::plain(
                "Kuranishi equation coordinate is out of range",
            ));
        }
        for row in self.kappa.values().flatten().map(|(_, row)| row) {
            model.spend(1 + row.len().max(1).ilog2() as usize)?;
        }
        Ok(self
            .kappa
            .values()
            .flatten()
            .filter_map(move |(mono, values)| {
                values
                    .binary_search_by_key(&index, |(i, _)| *i)
                    .ok()
                    .map(|i| (mono, values[i].1))
            }))
    }

    fn check_order(&self, order: usize) -> Result<()> {
        if order > self.order {
            Err(outside())
        } else {
            Ok(())
        }
    }

    pub(crate) fn specialized(
        &self,
        model: &mut Model,
        hom: &Hom,
        amplitudes: &[i64],
        weights: &[usize],
        order: usize,
        class: bool,
    ) -> Result<Value> {
        self.check_order(order)?;
        let p = hom.prime();
        if amplitudes.len() != self.directions || weights.len() != self.directions {
            return Err(Error::plain(
                "Kuranishi arc requires one amplitude and weight per direction",
            ));
        }
        model.spend(self.directions)?;
        if weights.contains(&0) {
            return Err(Error::plain("Kuranishi arc weights must be positive"));
        }
        let arc = Arc {
            p,
            weights,
            amplitudes: amplitudes.iter().map(|&x| sheaf::residue(p, x)).collect(),
            order,
        };
        if class {
            let support = arc.obstruction(model, &self.kappa, order)?;
            return crate::odd_ext::persist(model, support);
        }
        if let Some(first) = arc.first_obstruction(model, &self.kappa)? {
            return Err(Error::plain(format!(
                "Kuranishi arc has a nonzero obstruction in order {first}"
            )));
        }
        let mut coefficients: BTreeMap<usize, Cochain> = BTreeMap::new();
        let mut cells = 0;
        for (m, raw) in self.a.range(..=order).flat_map(|(_, layer)| layer) {
            let Some((e, c)) = arc.weighted(model, m)? else {
                continue;
            };
            let coefficient = coefficients.entry(e).or_default();
            cells -= coefficient.cells();
            add_raw(model, p, coefficient, raw, c)?;
            cells += coefficient.cells();
            sheaf::check_cells(cells + coefficients.len())?;
        }
        model.spend(order)?;
        let mut jet = Jet::new(model, Vec::new())?;
        for i in 1..=order {
            let cochain = coefficients.remove(&i).unwrap_or_default().into_terms();
            let Value::SheafC(cochain) = sheaf::persist(model, cochain)? else {
                unreachable!()
            };
            jet = jet.append(model, hom, cochain)?;
        }
        Ok(Value::SheafDeformation(jet))
    }
}

fn lookup<T>(polynomial: &Polynomial<Rc<[T]>>, degree: usize, key: &Monomial) -> Rc<[T]> {
    polynomial
        .get(&degree)
        .and_then(|layer| layer.get(key))
        .cloned()
        .unwrap_or_default()
}

fn outside() -> Error {
    Error::plain("Kuranishi coefficient is beyond the computed truncation")
}

fn check_support(cells: &mut usize, mono: &Monomial, raw: &[Term]) -> Result<()> {
    *cells += mono.cells() + raw.iter().map(|(g, _)| g.cells()).sum::<usize>();
    sheaf::check_cells(*cells)
}

fn schedule(
    model: &mut Model,
    pending: &mut BTreeSet<usize>,
    limit: usize,
    degree: usize,
    others: impl Iterator<Item = usize>,
) -> Result<()> {
    for e in others {
        model.spend(1)?;
        if let Some(next) = degree.checked_add(e).filter(|&next| next <= limit) {
            pending.insert(next);
            sheaf::check_cells(pending.len())?;
        }
    }
    Ok(())
}

struct Arc<'a> {
    p: Prime,
    amplitudes: Vec<u32>,
    weights: &'a [usize],
    order: usize,
}

impl Arc<'_> {
    fn weighted(&self, model: &mut Model, mono: &Monomial) -> Result<Option<(usize, u32)>> {
        let mut degree = 0;
        let mut scalar = 1;
        model.spend(mono.0.len())?;
        for &(i, e) in &mono.0 {
            let a = self.amplitudes[i];
            if a == 0 || e > (self.order - degree) / self.weights[i] {
                return Ok(None);
            }
            degree += e * self.weights[i];
            scalar = self.p.mul(scalar, power(model, self.p, a, e)?);
        }
        Ok(Some((degree, scalar)))
    }

    fn obstruction(
        &self,
        model: &mut Model,
        poly: &HPol,
        degree: usize,
    ) -> Result<Vec<Coordinate>> {
        let mut sum = Vec::new();
        for (m, row) in poly.range(..=degree).flat_map(|(_, layer)| layer) {
            if let Some((e, c)) = self.weighted(model, m)? {
                if e == degree {
                    sum = add_scaled(model, self.p, &sum, row, c)?;
                }
            }
        }
        Ok(sum)
    }

    fn first_obstruction(&self, model: &mut Model, poly: &HPol) -> Result<Option<usize>> {
        let mut sum: BTreeMap<usize, Vec<Coordinate>> = BTreeMap::new();
        let mut cells = 0;
        for (m, row) in poly.range(..=self.order).flat_map(|(_, layer)| layer) {
            if let Some((e, c)) = self.weighted(model, m)? {
                let output = sum.entry(e).or_default();
                cells -= 2 * output.len();
                *output = add_scaled(model, self.p, output, row, c)?;
                cells += 2 * output.len();
                sheaf::check_cells(cells + sum.len())?;
            }
        }
        Ok(sum
            .into_iter()
            .find_map(|(degree, row)| (!row.is_empty()).then_some(degree)))
    }
}

fn power(model: &mut Model, p: Prime, mut base: u32, mut exponent: usize) -> Result<u32> {
    model.spend(2 * (usize::BITS - exponent.leading_zeros()) as usize)?;
    let mut result = 1;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = p.mul(result, base);
        }
        base = p.mul(base, base);
        exponent >>= 1;
    }
    Ok(result)
}
