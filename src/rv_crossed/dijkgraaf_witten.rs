use std::rc::Rc;

use crate::complex::{for_each_boundary, Simplex};
use crate::error::{Error, Result};
use crate::model::{Model, Poset, TypeId, MAX_SHOWN_BYTES};
use crate::odd_cobar::Prime;
use crate::rv_cyclotomic;
use crate::rv_lens::{self, Rational};
use crate::syntax::{Label, ParsedExpr};

use super::flat::{pullback_layer, Flat};
use super::path_integral::{
    edge_name, order_layers, spanning_forest, FlatPresentation, Orbit, Quotient,
};
use super::postnikov::{Cochain, Postnikov};
use super::{Crossed, Group};

pub(crate) enum Command {
    Coefficients(Label, Label, usize),
    Declare {
        name: Label,
        base: Label,
        post: Label,
        cocycle: Box<ParsedExpr>,
        cycle: Vec<i128>,
    },
    Component(Label, usize),
    Class(Label, [usize; 2]),
    Phase(Label, Box<ParsedExpr>, Vec<i128>),
}

pub(crate) fn coefficients(
    model: &mut Model,
    group: Rc<Group>,
    group_name: Rc<str>,
    p: usize,
) -> Result<Postnikov> {
    if p != 2 {
        Prime::new(model, p)?;
    }
    let kernel = Rc::new(Group::cyclic(model, p)?);
    let action = (0..group.order()).map(|_| (0..p).collect()).collect();
    let crossed = Crossed::new(
        model,
        [
            (Rc::from("DW coefficient group"), kernel.clone()),
            (group_name, group),
        ],
        vec![0; p],
        action,
    )?;
    model.allocate(kernel.cells() + crossed.cells())?;
    Postnikov::new(model, Rc::new(crossed), 1)
}

pub(crate) struct Integral {
    site: Rc<Poset>,
    pub(crate) prime: usize,
    pub(crate) value: Rational,
    components: Vec<Component>,
}

struct Component {
    points: Vec<usize>,
    tree: Vec<[usize; 2]>,
    generators: usize,
    relators: usize,
    solutions: usize,
    classes: Vec<Class>,
    value: Rational,
}

struct Class {
    orbit: Orbit,
    phase: u32,
}

impl Integral {
    pub(crate) fn new(
        model: &mut Model,
        base: TypeId,
        post: &Postnikov,
        k: &Cochain,
        cycle: Vec<i128>,
    ) -> Result<Self> {
        require_twist(model, post, k)?;
        let site = model.poset(base).expect("named finite space");
        let cycle = Cycle::new(model, &site, cycle)?;
        let group = Quotient::new(model, &post.crossed)?;
        let forest = spanning_forest(model, &site)?;
        let prime = post.prime.get() as usize;
        let mut value = Rational::integer(1);
        let mut components = Vec::with_capacity(forest.len());
        for (points, tree) in forest {
            let geometry = FlatPresentation::new(model, &site, points, tree)?;
            let component = Component::new(model, geometry, &group, post, k, &cycle)?;
            value = rv_lens::multiply(model, prime, &value, &component.value)?;
            components.push(component);
        }
        model.allocate(4 + 2 * value.numerator.len() + 9 * components.len())?;
        Ok(Self {
            site,
            prime,
            value,
            components,
        })
    }

    pub(crate) fn component(&self, name: &str, i: usize) -> Result<String> {
        let data = self.get_component(i)?;
        let root = &self.site.names[data.points[0]];
        let tree = data
            .tree
            .iter()
            .map(|&edge| edge_name(&self.site, &data.points, edge))
            .collect::<Vec<_>>()
            .join(",");
        let value = &data.value;
        Ok(format!(
            "rvDWComponent[{name},{i}] = root={root}, vertices={}, tree=[{tree}], generators={}, relators={}, normalized_solutions={}, conjugacy_classes={}, value={:?}/{}",
            data.points.len(), data.generators, data.relators, data.solutions, data.classes.len(), value.numerator, value.denominator
        ))
    }

    pub(crate) fn class(
        &self,
        model: &mut Model,
        name: &str,
        [i, j]: [usize; 2],
    ) -> Result<String> {
        let component = self.get_component(i)?;
        let class = component.classes.get(j).ok_or_else(|| {
            Error::plain("Dijkgraaf–Witten conjugacy-class index is out of range")
        })?;
        let orbit = &class.orbit;
        let edges = orbit.edge_report(model, &self.site, &component.points)?;
        let result = format!(
            "rvDWClass[{name},{i},{j}] = labels=[{edges}], centralizer={}, orbit_size={}, phase={}",
            orbit.centralizer, orbit.size, class.phase
        );
        if result.len() > MAX_SHOWN_BYTES {
            return Err(Error::plain(
                "Dijkgraaf–Witten class report exceeds the display size limit",
            ));
        }
        Ok(result)
    }

    fn get_component(&self, i: usize) -> Result<&Component> {
        self.components
            .get(i)
            .ok_or_else(|| Error::plain("Dijkgraaf–Witten component index is out of range"))
    }
}

impl Component {
    fn new(
        model: &mut Model,
        geometry: FlatPresentation,
        group: &Quotient,
        post: &Postnikov,
        k: &Cochain,
        cycle: &Cycle,
    ) -> Result<Self> {
        let tetrahedra = &geometry.layers[3];
        let restricted = cycle.restrict(model, &geometry.points, tetrahedra, post.prime)?;
        let mut classes = Vec::new();
        let p = post.prime.get() as usize;
        let mut phase_counts = vec![0i128; p];
        let solutions = geometry.visit(model, group, |model, arrows, orbit| {
            let pulled = pullback_layer(model, post, arrows, tetrahedra, 3, k)?;
            let phase = pair(model, post.prime, &restricted, &pulled)?;
            // All representations in the orbit have the same phase.
            // Dividing the raw count by |G| gives weight 1/centralizer.
            phase_counts[phase as usize] += orbit.size as i128;
            model.allocate(6 + 3 * orbit.edges.len())?;
            classes.push(Class { orbit, phase });
            Ok(())
        })?;
        let coefficients = rv_cyclotomic::reduced(model, p, phase_counts)?;
        let value = rv_lens::fraction(model, &coefficients, group.order() as u128)?;
        let generators = geometry.generators();
        let relators = geometry.relators();
        let FlatPresentation { points, tree, .. } = geometry;
        model.allocate(points.len() + 2 * tree.len() + 2 * value.numerator.len())?;
        Ok(Self {
            points,
            tree,
            generators,
            relators,
            solutions,
            classes,
            value,
        })
    }
}

struct Cycle {
    tetrahedra: Vec<Simplex>,
    coefficients: Vec<i128>,
}

impl Cycle {
    fn new(model: &mut Model, site: &Poset, coefficients: Vec<i128>) -> Result<Self> {
        let points: Vec<_> = (0..site.names.len()).collect();
        let [_, _, triangles, tetrahedra] = order_layers(model, site, &points)?;
        if coefficients.len() != tetrahedra.len() {
            return Err(Error::plain(format!(
                "Dijkgraaf–Witten cycle requires {} integral tetrahedron coefficients",
                tetrahedra.len()
            )));
        }
        let mut boundary = vec![0i128; triangles.len()];
        for_each_boundary(model, &triangles, &tetrahedra, |_, index, faces| {
            let coefficient = coefficients[index];
            for (i, face) in faces.into_iter().enumerate() {
                let total = &mut boundary[face];
                *total = if i.is_multiple_of(2) {
                    total.checked_add(coefficient)
                } else {
                    total.checked_sub(coefficient)
                }
                .ok_or_else(|| {
                    Error::plain("integral 3-cycle boundary exceeds the signed 128-bit range")
                })?;
            }
            Ok(())
        })?;
        if let Some((index, _)) = boundary.iter().enumerate().find(|&(_, &x)| x != 0) {
            let vertices = triangles[index]
                .iter()
                .map(|&i| site.names[i].as_str())
                .collect::<Vec<_>>();
            return Err(Error::plain(format!(
                "Dijkgraaf–Witten integral chain is not a 3-cycle at '{}'",
                vertices.join("<")
            )));
        }
        Ok(Self {
            tetrahedra,
            coefficients,
        })
    }

    fn restrict(
        &self,
        model: &mut Model,
        points: &[usize],
        tetrahedra: &[Simplex],
        p: Prime,
    ) -> Result<Vec<u32>> {
        model.spend(tetrahedra.len() * 4 * (1 + self.tetrahedra.len().max(1).ilog2() as usize))?;
        Ok(tetrahedra
            .iter()
            .map(|chain| {
                let [a, b, c, d] = *chain.as_slice() else {
                    unreachable!()
                };
                let tetrahedron = [points[a], points[b], points[c], points[d]];
                let index = self
                    .tetrahedra
                    .binary_search_by(|x| x.as_slice().cmp(&tetrahedron))
                    .expect("component tetrahedron belongs to the full complex");
                residue(self.coefficients[index], p)
            })
            .collect())
    }
}

fn residue(value: i128, p: Prime) -> u32 {
    value.rem_euclid(i128::from(p.get())) as u32
}

fn pair(model: &mut Model, p: Prime, cycle: &[u32], pullback: &Cochain) -> Result<u32> {
    model.spend(pullback.0.len())?;
    Ok(pullback
        .0
        .iter()
        .fold(0, |sum, &(i, value)| p.add(sum, p.mul(cycle[i], value))))
}

fn require_twist(model: &mut Model, post: &Postnikov, k: &Cochain) -> Result<()> {
    if post.action.iter().any(|&x| x != 1) {
        return Err(Error::plain(
            "Dijkgraaf–Witten cocycles require trivial coefficient action",
        ));
    }
    post.require_cocycle(model, k)
}

pub(crate) fn phase(model: &mut Model, flat: &Flat, k: &Cochain, cycle: Vec<i128>) -> Result<u32> {
    let post = &flat.post;
    require_twist(model, post, k)?;
    let site = model.poset(flat.base).expect("flat base space");
    let cycle = Cycle::new(model, &site, cycle)?;
    let cycle: Vec<_> = cycle
        .coefficients
        .into_iter()
        .map(|x| residue(x, post.prime))
        .collect();
    let pulled = flat.pullback(model, 3, k)?;
    pair(model, post.prime, &cycle, &pulled)
}
