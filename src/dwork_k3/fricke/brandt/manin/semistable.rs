use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::lee_linear::Rational as Q;
use crate::model::{Model, Value};
use crate::syntax::{Label, ParsedExpr};

use super::super::super::hecke;
use super::super::Basis;
use super::linear::{self, Frame, Matrix, ONE, ZERO};
use super::{Operator, Space};
use geometry::Edges;
use lattice::{self as lat, Lattice};

mod geometry;
mod lattice;

pub(crate) enum Command {
    Define(Label, Label, Label),
    Describe(Label),
    WeilDeligne(Label),
    Hecke(Label, usize),
    Equation(usize, usize, usize),
}

impl Command {
    pub(crate) fn execute(self, model: &mut Model) -> Result<Option<String>> {
        match self {
            Self::Define(name, manin, brandt) => {
                let space = model.manin_space(&manin)?;
                let basis = model.fricke_ss(&brandt)?;
                let data =
                    Semistable::new(model, space, basis).map_err(|error| error.at(name.span))?;
                model.define_fricke_bad(name, data)?;
                Ok(None)
            }
            Self::Describe(name) => Ok(Some(model.fricke_bad(&name)?.describe())),
            Self::WeilDeligne(name) => {
                let data = model.fricke_bad(&name)?;
                data.wd(model).map(Some)
            }
            Self::Hecke(name, ell) => {
                let data = model.fricke_bad(&name)?;
                data.hecke_report(model, ell).map(Some)
            }
            Self::Equation(ell, i, j) => {
                let equation = model.fricke_oriented_equation(ell)?;
                Ok(Some(format!(
                    "ordered Γ0(2) Hecke coefficient [{ell},{i},{j}] = {}",
                    equation.coefficient_text(i, j)
                )))
            }
        }
    }
}

pub(crate) enum Expression {
    Vector(Label, bool, Vec<i128>),
    Component(Label, usize),
    Node(Label, usize),
    Action(Action, Box<ParsedExpr>),
}

#[derive(Clone, Copy)]
pub(crate) enum Action {
    Monodromy,
    Class,
    Lift,
    Frobenius,
    Atkin(usize),
    Hecke(usize),
    Negate,
}

#[derive(Clone, Copy, Eq, PartialEq)]
pub(crate) enum Grade {
    Cycle,
    Dual,
    Component,
}

pub(crate) struct Semistable {
    basis: Rc<Basis>,
    space: Rc<Space>,
    edges: Edges,
    lattice: Lattice,
    frobenius: Matrix,
    w2: Matrix,
    quotient: Matrix,
    operators: RefCell<BTreeMap<usize, Rc<HeckeAction>>>,
}

struct HeckeAction {
    edges: Matrix,
    cycles: Matrix,
    component_bridge: Matrix,
}

impl Semistable {
    fn new(model: &mut Model, space: Rc<Space>, basis: Rc<Basis>) -> Result<Self> {
        if space.geometry.prime != basis.p as usize {
            return Err(Error::plain(
                "semistable Manin level must be twice the Brandt characteristic",
            ));
        }
        let edges = Edges::new(model, &basis)?;
        let lattice = Lattice::new(model, &edges.widths)?;
        let frobenius = lat::restrict(&lat::permutation(&edges.frobenius))?;
        let w2 = lat::restrict(&lat::permutation(&edges.w2))?;
        let g = frobenius.len();
        if space.plus.vectors.len() != g {
            return Err(certificate("cycle rank equals the full cuspidal genus"));
        }
        let quotient = quotient_map(basis.rank(), &edges.quotient);
        let mut equations = vec![vec![ZERO; g]; g];
        for (i, row) in equations.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                *entry = linear::sub(w2[i][j], if i == j { ONE } else { ZERO })?;
            }
        }
        let plus = Frame::nullspace(model, g, equations)?;
        let quotient_on_plus =
            linear::multiply(model, &quotient, &linear::columns(&plus.vectors, g))?;
        if plus.vectors.len() + 1 != basis.rank()
            || linear::inverse(model, &quotient_on_plus)?.is_none()
        {
            return Err(certificate("Fricke quotient identifies W2+ cycles over Q"));
        }
        for involution in [&frobenius, &w2] {
            if linear::multiply(model, involution, involution)? != linear::identity(g) {
                return Err(certificate("geometric cycle permutations are involutions"));
            }
            lat::intertwine(
                model,
                involution,
                &lat::transpose(involution),
                &lattice.gram,
            )?;
        }
        lat::intertwine(model, &frobenius, &frobenius, &w2)?;
        let edge_projection = edge_map(basis.rank(), &edges.quotient);
        lat::intertwine(
            model,
            &lat::permutation(&edges.frobenius),
            &lat::permutation(&basis.frobenius),
            &edge_projection,
        )?;
        lat::intertwine(
            model,
            &lat::permutation(&edges.w2),
            &linear::identity(basis.rank()),
            &edge_projection,
        )?;
        let wp = space.operator(model, Operator::Atkin(basis.p as usize))?;
        let minus_wp = lat::negative(&space.plus.restrict(model, &wp)?)?;
        if linear::characteristic(model, &frobenius)? != linear::characteristic(model, &minus_wp)? {
            return Err(certificate(
                "bad Frobenius agrees with −Wp on full Manin cusp homology",
            ));
        }
        model.allocate(28 + 20 * g * g + 12 * g)?;
        let data = Self {
            basis,
            space,
            edges,
            lattice,
            frobenius,
            w2,
            quotient,
            operators: RefCell::new(BTreeMap::new()),
        };
        data.wd_matrices(model)?;
        Ok(data)
    }

    pub(crate) fn rank(&self) -> usize {
        self.frobenius.len()
    }
    pub(crate) fn order(&self) -> usize {
        self.lattice.order
    }
    pub(crate) fn field_label(&self, span: crate::error::Span) -> Label {
        self.basis.field_label(span)
    }
    pub(crate) fn node(&self, index: usize) -> Result<Vec<u64>> {
        self.edges
            .parameters
            .get(index)
            .cloned()
            .ok_or_else(|| Error::plain("semistable node index is out of range"))
    }

    pub(crate) fn check_action(&self, action: Action, input: Grade) -> Result<Grade> {
        match action {
            Action::Monodromy if input == Grade::Cycle => Ok(Grade::Dual),
            Action::Class if input == Grade::Dual => Ok(Grade::Component),
            Action::Lift if input == Grade::Component => Ok(Grade::Dual),
            Action::Monodromy | Action::Class | Action::Lift => Err(Error::plain(
                "semistable monodromy, component, or lift has the wrong input lattice",
            )),
            Action::Atkin(q) if q != 2 && q != self.basis.p as usize => {
                Err(Error::plain("semistable Atkin–Lehner index must be 2 or p"))
            }
            Action::Hecke(ell) if ell == self.basis.p as usize => Err(Error::plain(
                "oriented Hecke level must be prime to the bad characteristic",
            )),
            Action::Hecke(ell) => {
                hecke::prime_level(ell)?;
                Ok(input)
            }
            _ => Ok(input),
        }
    }

    pub(crate) fn apply(
        &self,
        model: &mut Model,
        action: Action,
        grade: Grade,
        input: Value,
    ) -> Result<Value> {
        match (action, input) {
            (Action::Monodromy, Value::FrickeBadVector(value)) => {
                let result = lat::integral_apply(model, &self.lattice.gram, &value)?;
                persist(model, result)
            }
            (Action::Class, Value::FrickeBadVector(value)) => {
                Ok(Value::FrickeComponent(self.lattice.class(&value)))
            }
            (Action::Lift, Value::FrickeComponent(residue)) => {
                persist(model, self.lattice.lift(residue)?)
            }
            (Action::Negate, Value::FrickeBadVector(value)) => {
                let result = value
                    .iter()
                    .map(|v| v.checked_neg().ok_or_else(arithmetic))
                    .collect::<Result<_>>()?;
                persist(model, result)
            }
            (Action::Negate, Value::FrickeComponent(residue)) => Ok(Value::FrickeComponent(
                (self.order() - residue) % self.order(),
            )),
            (action, Value::FrickeBadVector(value)) => {
                let matrix = self.action_matrix(model, action, grade)?;
                let output = lat::integral_apply(model, &matrix, &value)?;
                persist(model, output)
            }
            (action, Value::FrickeComponent(residue)) => {
                let matrix = self.action_matrix(model, action, Grade::Dual)?;
                let output = lat::integral_apply(model, &matrix, &self.lattice.lift(residue)?)?;
                Ok(Value::FrickeComponent(self.lattice.class(&output)))
            }
            _ => Err(certificate(
                "typed lattice operation uses the certified domain",
            )),
        }
    }

    pub(crate) fn link(&self, x: usize, y: usize) -> Result<String> {
        let pairing = self.lattice.link(x, y)?;
        Ok(format!(
            "component linking (mod Z) = {}",
            linear::describe(&[pairing]).trim_matches(['[', ']'])
        ))
    }

    fn action_matrix(&self, model: &mut Model, action: Action, grade: Grade) -> Result<Matrix> {
        let cycle = match action {
            Action::Frobenius => self.frobenius.clone(),
            Action::Atkin(2) => self.w2.clone(),
            Action::Atkin(_) => lat::negative(&self.frobenius)?,
            Action::Hecke(ell) => self.hecke(model, ell)?.cycles.clone(),
            _ => return Err(certificate("linear graph action")),
        };
        // On cycle isometries the dual is inverse transpose; these are
        // involutions. For Hecke the dual action is the adjoint transpose.
        Ok(if grade == Grade::Cycle {
            cycle
        } else {
            lat::transpose(&cycle)
        })
    }

    fn hecke(&self, model: &mut Model, ell: usize) -> Result<Rc<HeckeAction>> {
        self.check_action(Action::Hecke(ell), Grade::Cycle)?;
        if let Some(matrix) = self.operators.borrow().get(&ell) {
            return Ok(matrix.clone());
        }
        let equation = model.fricke_oriented_equation(ell)?;
        let full = self
            .edges
            .oriented_matrix(model, &self.basis, &equation, ell)?;
        for (i, row) in full.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                if c.numerator * self.edges.widths[i] != full[j][i].numerator * self.edges.widths[j]
                {
                    return Err(certificate(
                        "ordered Hecke operator is self-adjoint for edge widths",
                    ));
                }
            }
        }
        for involution in [&self.edges.frobenius, &self.edges.w2] {
            let permutation = lat::permutation(involution);
            lat::intertwine(model, &full, &full, &permutation)?;
        }
        let graph = lat::restrict(&full)?;
        lat::intertwine(model, &graph, &lat::transpose(&graph), &self.lattice.gram)?;
        for involution in [&self.frobenius, &self.w2] {
            lat::intertwine(model, &graph, &graph, involution)?;
        }
        for prior in self.operators.borrow().values() {
            lat::intertwine(model, &graph, &graph, &prior.cycles)?;
        }
        let fricke = self.basis.operator(model, ell)?;
        let fricke_matrix: Matrix = fricke
            .0
            .iter()
            .map(|row| row.iter().copied().map(Q::integer).collect())
            .collect();
        let fricke_push = lat::transpose(&fricke_matrix);
        lat::intertwine(
            model,
            &full,
            &fricke_push,
            &edge_map(self.basis.rank(), &self.edges.quotient),
        )?;
        let fricke_cycle = lat::restrict(&fricke_push)?;
        lat::intertwine(model, &graph, &fricke_cycle, &self.quotient)?;
        let manin = self.space.operator(model, Operator::Hecke(ell))?;
        let manin = self.space.plus.restrict(model, &manin)?;
        if linear::characteristic(model, &graph)? != linear::characteristic(model, &manin)? {
            return Err(certificate(
                "oriented Hecke spectrum equals the full cuspidal J+ Manin spectrum",
            ));
        }
        let mut difference = lat::transpose(&graph);
        for (i, row) in difference.iter_mut().enumerate() {
            row[i] = linear::sub(row[i], Q::integer((ell + 1) as i128))?;
        }
        let bridge = linear::multiply(model, &self.lattice.inverse, &difference)?;
        if bridge.iter().flatten().any(|c| c.denominator != 1) {
            return Err(certificate(
                "Hecke acts by ell+1 on the integral component quotient",
            ));
        }
        model.allocate(8 + 8 * graph.len().pow(2) + 4 * full.len().pow(2))?;
        let action = Rc::new(HeckeAction {
            edges: full,
            cycles: graph,
            component_bridge: bridge,
        });
        self.operators.borrow_mut().insert(ell, action.clone());
        Ok(action)
    }

    fn describe(&self) -> String {
        let exceptional = self
            .edges
            .widths
            .iter()
            .enumerate()
            .filter_map(|(i, &w)| (w == 2).then_some(i))
            .collect::<Vec<_>>();
        format!(
            "X0({}) at bad p={}; two P1 components; {} geometric nodes; toric rank {}\n\
             oriented t = {:?}\n\
             ordered endpoint j pairs = {:?}\n\
             edge widths = {:?}; exceptional edge indices after resolution = {exceptional:?}\n\
             geometric Frobenius permutation = {:?}; W2 permutation = {:?}; Wp reverses the two components and acts by −F on cycles\n\
             Gram = {}; Smith = {:?}; component order = {}",
            self.space.geometry.level,self.basis.p,self.edges.parameters.len(),self.rank(),
            self.edges.parameters,self.edges.endpoints,self.edges.widths,self.edges.frobenius,
            self.edges.w2,linear::describe_matrix(&self.lattice.gram),self.lattice.smith,self.order()
        )
    }

    fn wd_matrices(&self, model: &mut Model) -> Result<(Matrix, Matrix)> {
        let g = self.rank();
        let mut n = vec![vec![ZERO; 2 * g]; 2 * g];
        let mut phi = n.clone();
        for i in 0..g {
            for j in 0..g {
                n[i][g + j] = self.lattice.gram[i][j];
                phi[i][j] = self.frobenius[j][i];
                phi[i + g][j + g] =
                    self.frobenius[i][j].multiply(Q::integer(i128::from(self.basis.p)))?;
            }
        }
        let n_phi = linear::multiply(model, &n, &phi)?;
        let mut other = linear::multiply(model, &phi, &n)?;
        for c in other.iter_mut().flatten() {
            *c = c.multiply(Q::integer(i128::from(self.basis.p)))?;
        }
        if linear::multiply(model, &n, &n)?
            .iter()
            .flatten()
            .any(|c| *c != ZERO)
            || n_phi != other
        {
            return Err(certificate("N²=0 and NΦ=pΦN"));
        }
        Ok((n, phi))
    }

    fn wd(&self, model: &mut Model) -> Result<String> {
        let (n, phi) = self.wd_matrices(model)?;
        let polynomial = linear::characteristic(model, &self.frobenius)?;
        Ok(format!(
            "Weil–Deligne on X∨ ⊕ X(−1), N = {}, Φ = {}; N²=0; NΦ=pΦN\n\
             bad Euler factor (ascending) = {}; conductor exponent = {}",
            linear::describe_matrix(&n),
            linear::describe_matrix(&phi),
            linear::describe(&polynomial),
            self.rank()
        ))
    }

    fn hecke_report(&self, model: &mut Model, ell: usize) -> Result<String> {
        let t = self.hecke(model, ell)?;
        let characteristic = linear::characteristic(model, &t.cycles)?;
        Ok(format!(
            "oriented Γ0(2) T{ell} on cycles = {}; characteristic (descending) = {}; full Manin J+ and Fricke quotient actions certified; component action = {} mod {}\n\
             ordered edge push (target rows) = {}\n\
             integral component bridge = {}",
            linear::describe_matrix(&t.cycles),linear::describe(&characteristic),ell+1,self.order(),
            linear::describe_matrix(&t.edges),linear::describe_matrix(&t.component_bridge)
        ))
    }
}

fn edge_map(n: usize, indices: &[usize]) -> Matrix {
    (0..n)
        .map(|i| {
            indices
                .iter()
                .map(|&j| Q::integer(i128::from(i == j)))
                .collect()
        })
        .collect()
}

fn quotient_map(n: usize, indices: &[usize]) -> Matrix {
    let last = indices[indices.len() - 1];
    (0..n - 1)
        .map(|i| {
            indices[..indices.len() - 1]
                .iter()
                .map(|&j| Q::integer(i128::from(i == j) - i128::from(i == last)))
                .collect()
        })
        .collect()
}

pub(crate) fn persist(model: &mut Model, value: Vec<i128>) -> Result<Value> {
    model.allocate(3 + 2 * value.len())?;
    Ok(Value::FrickeBadVector(value.into()))
}

pub(crate) fn add(model: &mut Model, left: &[i128], right: &[i128]) -> Result<Value> {
    let values = left
        .iter()
        .zip(right)
        .map(|(a, b)| a.checked_add(*b).ok_or_else(arithmetic))
        .collect::<Result<_>>()?;
    persist(model, values)
}

fn arithmetic() -> Error {
    Error::plain("semistable cycle arithmetic exceeds the exact signed coefficient range")
}

fn certificate(identity: &str) -> Error {
    Error::plain(format!("semistable Fricke certificate failed: {identity}"))
}
