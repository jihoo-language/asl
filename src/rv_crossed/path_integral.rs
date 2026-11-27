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
use presentation::{Presentation, Quotient};

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
    edges: Vec<([usize; 2], usize)>,
    centralizer: usize,
    orbit_size: usize,
    betti: [usize; 3],
    obstruction: bool,
    fixed: Rational,
    weight: Rational,
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
            .map(|&edge| self.edge_name(&data.points, edge))
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
        model.spend(class.edges.len())?;
        let edges = class
            .edges
            .iter()
            .map(|&(edge, g)| format!("{}:{g}", self.edge_name(&component.points, edge)))
            .collect::<Vec<_>>()
            .join(",");
        let result = format!(
            "rvHomotopyClass[{name},{i},{j}] = labels=[{edges}], centralizer={}, orbit_size={}, betti={:?}, obstructed={}, fixed_cardinality={}, weight={}",
            class.centralizer, class.orbit_size, class.betti, class.obstruction, show(class.fixed), show(class.weight),
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

    fn edge_name(&self, points: &[usize], [a, b]: [usize; 2]) -> String {
        format!(
            "{}<{}",
            self.site.names[points[a]], self.site.names[points[b]]
        )
    }
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
        let mut order = OrderComplex::induced(model, site, &points)?;
        let v = order.vertices();
        let e = order.extend(model, &v)?;
        let f = order.extend(model, &e)?;
        let t = order.extend(model, &f)?;
        let presentation = Presentation::new(model, points.len(), &tree, &e, &f)?;
        let generators = presentation.generators.len();
        let relators = presentation.relators.len();
        let n = group.order();
        let assignments = n.checked_pow(generators as u32).filter(|&x| x <= MAX_ASSIGNMENTS)
            .ok_or_else(|| Error::plain("homotopy enumeration exceeds the normalized assignment limit after elimination"))?;
        let mut values = vec![0; presentation.original_generators + 1];
        let mut arrows = vec![vec![None; points.len()]; points.len()];
        let mut classes = Vec::new();
        let mut solutions = 0;
        let mut orbit_sizes = 0;
        let mut total = Rational::ZERO;
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
            let orbit_size = n / centralizer;
            orbit_sizes += orbit_size;
            presentation.recover(model, group, &mut values)?;
            let mut edges = Vec::new();
            model.spend(e.len() + f.len())?;
            for (edge, &variable) in e.iter().zip(&presentation.edge_generators) {
                let &[a, b] = edge.as_slice() else {
                    unreachable!()
                };
                let g = values[variable];
                arrows[a][b] = Some(g);
                if g != 0 {
                    edges.push(([a, b], g));
                }
            }
            for triangle in &f {
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
            let (betti, obstruction) = local_weight(model, post, k, &arrows, [&v, &e, &f, &t])?;
            let fixed = if obstruction {
                Rational::ZERO
            } else {
                let [h0, h1, h2] = betti;
                prime_power(post.prime.get(), h0 as i64 - h1 as i64 + h2 as i64)?
            };
            let weight = exact(fixed.multiply(Rational::new(1, centralizer as i128)?))?;
            total = exact(total.add(weight))?;
            model.allocate(12 + edges.len() * 3)?;
            classes.push(Class {
                edges,
                centralizer,
                orbit_size,
                betti,
                obstruction,
                fixed,
                weight,
            });
        }
        if solutions != orbit_sizes {
            return Err(Error::plain(
                "homotopy conjugacy classes failed the orbit-stabilizer check",
            ));
        }
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

fn spanning_forest(model: &mut Model, site: &Poset) -> Result<Vec<SpanningTree>> {
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
