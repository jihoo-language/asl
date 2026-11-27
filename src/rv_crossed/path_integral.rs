use std::collections::VecDeque;
use std::rc::Rc;

use crate::complex::{OrderComplex, Simplex};
use crate::error::{Error, Result};
use crate::lee_linear::Rational;
use crate::model::{Model, Poset, TypeId, MAX_SHOWN_BYTES};
use crate::prime_linear::Basis;
use crate::set_family::Diagram;
use crate::syntax::{Label, ParsedExpr};

use super::flat::{pullback_layer, twisted_differential, Arrows};
use super::postnikov::{Cochain, Postnikov};

mod presentation;
use presentation::Presentation;
pub(super) use presentation::Quotient;

const MAX_ASSIGNMENTS: usize = 250_000;

pub(crate) enum Command {
    Declare(Label, Label, Label, Box<ParsedExpr>),
    Cardinality(Label),
    Component(Label, usize),
    Class(Label, [usize; 2]),
}

pub(crate) struct Integral {
    site: Rc<Poset>,
    components: Vec<Component>,
    total: Rational,
}

struct Component {
    points: Vec<usize>,
    tree: Vec<[usize; 2]>,
    generators: usize,
    relators: usize,
    solutions: usize,
    classes: Vec<Class>,
    total: Rational,
}

struct Class {
    orbit: Orbit,
    betti: [usize; 3],
    obstruction: bool,
    fixed: Rational,
    weight: Rational,
}

pub(super) struct Orbit {
    pub(super) edges: Vec<([usize; 2], usize)>,
    pub(super) centralizer: usize,
    pub(super) size: usize,
}

impl Orbit {
    pub(super) fn edge_report(
        &self,
        model: &mut Model,
        site: &Poset,
        points: &[usize],
    ) -> Result<String> {
        model.spend(self.edges.len())?;
        Ok(self
            .edges
            .iter()
            .map(|&(edge, g)| format!("{}:{g}", edge_name(site, points, edge)))
            .collect::<Vec<_>>()
            .join(","))
    }
}

pub(super) struct FlatPresentation {
    pub(super) points: Vec<usize>,
    pub(super) tree: Vec<[usize; 2]>,
    pub(super) layers: [Vec<Simplex>; 4],
    presentation: Presentation,
}

impl Integral {
    pub(crate) fn new(
        model: &mut Model,
        base: TypeId,
        post: &Postnikov,
        k: &Cochain,
    ) -> Result<Self> {
        post.require_cocycle(model, k)?;
        let site = model.poset(base).expect("named finite space");
        let group = Quotient::new(model, &post.crossed)?;
        let forest = spanning_forest(model, &site)?;
        let mut total = Rational::integer(1);
        let mut components = Vec::with_capacity(forest.len());
        for (points, tree) in forest {
            let component = Component::new(model, &site, points, tree, &group, post, k)?;
            total = exact(total.multiply(component.total))?;
            components.push(component);
        }
        model.allocate(1 + 10 * components.len())?;
        Ok(Self {
            site,
            components,
            total,
        })
    }

    pub(crate) fn cardinality(&self, name: &str) -> String {
        format!(
            "rvHomotopyCardinality[{name}] = {}, components={}",
            show(self.total),
            self.components.len()
        )
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
        Ok(format!(
            "rvHomotopyComponent[{name},{i}] = root={root}, vertices={}, tree=[{tree}], generators={}, relators={}, normalized_solutions={}, conjugacy_classes={}, cardinality={}",
            data.points.len(), data.generators, data.relators, data.solutions, data.classes.len(), show(data.total),
        ))
    }

    pub(crate) fn class(
        &self,
        model: &mut Model,
        name: &str,
        [i, j]: [usize; 2],
    ) -> Result<String> {
        let component = self.get_component(i)?;
        let class = component
            .classes
            .get(j)
            .ok_or_else(|| Error::plain("homotopy conjugacy-class index is out of range"))?;
        let edges = class
            .orbit
            .edge_report(model, &self.site, &component.points)?;
        let result = format!(
            "rvHomotopyClass[{name},{i},{j}] = labels=[{edges}], centralizer={}, orbit_size={}, betti={:?}, obstructed={}, fixed_cardinality={}, weight={}",
            class.orbit.centralizer, class.orbit.size, class.betti, class.obstruction, show(class.fixed), show(class.weight),
        );
        if result.len() > MAX_SHOWN_BYTES {
            return Err(Error::plain(
                "homotopy class report exceeds the display size limit",
            ));
        }
        Ok(result)
    }

    fn get_component(&self, i: usize) -> Result<&Component> {
        self.components
            .get(i)
            .ok_or_else(|| Error::plain("homotopy component index is out of range"))
    }
}

pub(super) fn edge_name(site: &Poset, points: &[usize], [a, b]: [usize; 2]) -> String {
    format!("{}<{}", site.names[points[a]], site.names[points[b]])
}

impl Component {
    fn new(
        model: &mut Model,
        site: &Poset,
        points: Vec<usize>,
        tree: Vec<[usize; 2]>,
        group: &Quotient,
        post: &Postnikov,
        k: &Cochain,
    ) -> Result<Self> {
        let geometry = FlatPresentation::new(model, site, points, tree)?;
        let layers = geometry.layers.each_ref().map(Vec::as_slice);
        let generators = geometry.generators();
        let relators = geometry.relators();
        let mut classes = Vec::new();
        let mut total = Rational::ZERO;
        let solutions = geometry.visit(model, group, |model, arrows, orbit| {
            let (betti, obstruction) = local_weight(model, post, k, arrows, layers)?;
            let fixed = if obstruction {
                Rational::ZERO
            } else {
                let [h0, h1, h2] = betti;
                prime_power(post.prime.get(), h0 as i64 - h1 as i64 + h2 as i64)?
            };
            let weight = exact(fixed.multiply(Rational::new(1, orbit.centralizer as i128)?))?;
            total = exact(total.add(weight))?;
            model.allocate(12 + orbit.edges.len() * 3)?;
            classes.push(Class {
                orbit,
                betti,
                obstruction,
                fixed,
                weight,
            });
            Ok(())
        })?;
        let FlatPresentation { points, tree, .. } = geometry;
        model.allocate(points.len() + 2 * tree.len())?;
        Ok(Self {
            points,
            tree,
            generators,
            relators,
            solutions,
            classes,
            total,
        })
    }
}

impl FlatPresentation {
    pub(super) fn new(
        model: &mut Model,
        site: &Poset,
        points: Vec<usize>,
        tree: Vec<[usize; 2]>,
    ) -> Result<Self> {
        let layers = order_layers(model, site, &points)?;
        let presentation = Presentation::new(model, points.len(), &tree, &layers[1], &layers[2])?;
        Ok(Self {
            points,
            tree,
            layers,
            presentation,
        })
    }

    pub(super) fn generators(&self) -> usize {
        self.presentation.generators.len()
    }

    pub(super) fn relators(&self) -> usize {
        self.presentation.relators.len()
    }

    pub(super) fn visit(
        &self,
        model: &mut Model,
        group: &Quotient,
        mut visitor: impl FnMut(&mut Model, &Arrows, Orbit) -> Result<()>,
    ) -> Result<usize> {
        let presentation = &self.presentation;
        let [_, edges, triangles, _] = &self.layers;
        let generators = self.generators();
        let n = group.order();
        let assignments = n.checked_pow(generators as u32).filter(|&x| x <= MAX_ASSIGNMENTS)
            .ok_or_else(|| Error::plain("homotopy enumeration exceeds the normalized assignment limit after elimination"))?;
        let mut values = vec![0; presentation.original_generators + 1];
        let mut arrows = vec![vec![None; self.points.len()]; self.points.len()];
        let mut solutions = 0;
        let mut orbit_sizes = 0;
        for mut sequence in 0..assignments {
            model.spend(1 + generators)?;
            for &g in presentation.generators.iter().rev() {
                values[g] = sequence % n;
                sequence /= n;
            }
            if !presentation.satisfies(model, group, &values)? {
                continue;
            }
            solutions += 1;
            let Some(centralizer) = presentation.centralizer(model, group, &values)? else {
                continue;
            };
            let size = n / centralizer;
            orbit_sizes += size;
            presentation.recover(model, group, &mut values)?;
            let mut nonidentity = Vec::new();
            model.spend(edges.len() + triangles.len())?;
            for (edge, &variable) in edges.iter().zip(&presentation.edge_generators) {
                let &[a, b] = edge.as_slice() else {
                    unreachable!()
                };
                let g = values[variable];
                arrows[a][b] = Some(g);
                if g != 0 {
                    nonidentity.push(([a, b], g));
                }
            }
            for triangle in triangles {
                let &[a, b, c] = triangle.as_slice() else {
                    unreachable!()
                };
                if group.mul(arrows[a][b].unwrap(), arrows[b][c].unwrap()) != arrows[a][c].unwrap()
                {
                    return Err(Error::plain(
                        "recovered edge labels failed a triangle relation",
                    ));
                }
            }
            visitor(
                model,
                &arrows,
                Orbit {
                    edges: nonidentity,
                    centralizer,
                    size,
                },
            )?;
        }
        if solutions != orbit_sizes {
            return Err(Error::plain(
                "homotopy conjugacy classes failed the orbit-stabilizer check",
            ));
        }
        Ok(solutions)
    }
}

pub(super) fn order_layers(
    model: &mut Model,
    site: &Poset,
    points: &[usize],
) -> Result<[Vec<Simplex>; 4]> {
    let mut order = OrderComplex::induced(model, site, points)?;
    let vertices = order.vertices();
    let edges = order.extend(model, &vertices)?;
    let triangles = order.extend(model, &edges)?;
    let tetrahedra = order.extend(model, &triangles)?;
    Ok([vertices, edges, triangles, tetrahedra])
}

fn local_weight(
    model: &mut Model,
    post: &Postnikov,
    k: &Cochain,
    arrows: &Arrows,
    layers: [&[Simplex]; 4],
) -> Result<([usize; 3], bool)> {
    let mut ranks = [0; 4];
    let mut obstruction = false;
    for i in 0..3 {
        let (_, columns) = twisted_differential(model, post, arrows, layers[i], layers[i + 1])?;
        let mut basis = Basis::default();
        for column in columns {
            basis.insert(model, post.prime, column)?;
        }
        ranks[i + 1] = basis.rank();
        if i == 2 {
            let pulled = pullback_layer(model, post, arrows, layers[3], 3, k)?;
            obstruction = !basis.reduce(model, post.prime, pulled.0)?.is_empty();
        }
    }
    let mut betti = [0; 3];
    for i in 0..3 {
        betti[i] = layers[i]
            .len()
            .checked_sub(ranks[i] + ranks[i + 1])
            .ok_or_else(|| Error::plain("twisted cochain ranks are inconsistent"))?;
    }
    Ok((betti, obstruction))
}

type SpanningTree = (Vec<usize>, Vec<[usize; 2]>);

pub(super) fn spanning_forest(model: &mut Model, site: &Poset) -> Result<Vec<SpanningTree>> {
    let (_, covers) = Diagram::site_orders(model, site)?;
    let mut neighbors = vec![Vec::new(); covers.len()];
    for (u, row) in covers.iter().enumerate() {
        for &v in row {
            neighbors[u].push(v);
            neighbors[v].push(u);
        }
    }
    for row in &mut neighbors {
        row.sort_unstable();
    }
    let mut visited = vec![false; covers.len()];
    let mut forest = Vec::new();
    for root in 0..covers.len() {
        if visited[root] {
            continue;
        }
        let mut points = Vec::new();
        let mut tree = Vec::new();
        let mut queue = VecDeque::from([root]);
        visited[root] = true;
        while let Some(u) = queue.pop_front() {
            points.push(u);
            model.spend(1 + neighbors[u].len())?;
            for &v in &neighbors[u] {
                if !std::mem::replace(&mut visited[v], true) {
                    let edge = if site.strict_less(u, v) {
                        [u, v]
                    } else {
                        [v, u]
                    };
                    tree.push(edge);
                    queue.push_back(v);
                }
            }
        }
        points.sort_unstable();
        for edge in &mut tree {
            *edge = edge.map(|i| points.binary_search(&i).expect("tree endpoint"));
        }
        tree.sort_unstable();
        forest.push((points, tree));
    }
    Ok(forest)
}

fn rational_limit() -> Error {
    Error::plain("homotopy cardinality exceeds the exact signed 128-bit rational range")
}

fn exact<T>(value: Result<T>) -> Result<T> {
    value.map_err(|_| rational_limit())
}

fn prime_power(prime: u32, exponent: i64) -> Result<Rational> {
    let value = i128::from(prime)
        .checked_pow(exponent.unsigned_abs() as u32)
        .ok_or_else(rational_limit)?;
    if exponent < 0 {
        Rational::new(1, value)
    } else {
        Ok(Rational::integer(value))
    }
}

fn show(value: Rational) -> String {
    if value.denominator == 1 {
        value.numerator.to_string()
    } else {
        format!("{}/{}", value.numerator, value.denominator)
    }
}
