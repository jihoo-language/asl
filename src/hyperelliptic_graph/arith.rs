use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::integer_smith::{self, arithmetic, Matrix, Track, Vector};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;
use crate::rv_qp_quad::{gauge::Coefficients, Element};
use crate::stable_marks::Tree;

use super::{period_inner, Cover, Cycle, Grade};

pub(crate) type Rational = (i128, i128);

#[derive(Eq, Ord, PartialEq, PartialOrd)]
pub(super) struct Chart {
    polynomial: Vec<u64>,
    pub(super) genus: Option<usize>,
    pub(super) nonsquare: bool,
}

impl Chart {
    pub(super) fn normalized(
        model: &mut Model,
        prime: Prime,
        leading: u32,
        roots: impl ExactSizeIterator<Item = u32>,
        odd_at_infinity: bool,
        nonsquare: bool,
    ) -> Result<Self> {
        let ramified = roots.len() + usize::from(odd_at_infinity);
        if !ramified.is_multiple_of(2) {
            return Err(Error::plain(
                "normalized hyperelliptic chart has an odd branch count",
            ));
        }
        let genus = (ramified != 0).then(|| ramified / 2 - 1);
        let mut coefficients = vec![leading];
        for root in roots {
            model.spend(1 + coefficients.len())?;
            let mut next = vec![0; coefficients.len() + 1];
            for (i, c) in coefficients.into_iter().enumerate() {
                next[i] = prime.add(next[i], prime.negate(prime.mul(root, c)));
                next[i + 1] = prime.add(next[i + 1], c);
            }
            coefficients = next;
        }
        Ok(Self {
            polynomial: coefficients.into_iter().map(u64::from).collect(),
            genus,
            nonsquare,
        })
    }

    pub(super) fn cells(&self) -> usize {
        4 + self.polynomial.len()
    }
}

#[derive(Clone, Copy)]
struct Contact {
    depth: i64,
    unit: u32,
}

struct BranchData {
    prime: Prime,
    leading: u32,
    pairs: Vec<Vec<Contact>>,
}

impl BranchData {
    fn new(model: &mut Model, prime: Prime, leading: Rational, roots: &[Rational]) -> Result<Self> {
        let leading = Element::new(leading.0, 0, leading.1)?;
        if leading.rational_valuation(model, prime)? != Some(0) {
            return Err(Error::plain(
                "hyperelliptic leading coefficient must be a p-adic unit",
            ));
        }
        let leading = leading.rational_angle(model, prime)?;
        let mut points = Vec::with_capacity(roots.len());
        let mut inverse = Vec::with_capacity(roots.len());
        let mut residues = BTreeSet::new();
        let modulus = i128::from(prime.get());
        for (i, &(a, c)) in roots.iter().enumerate() {
            let root = Element::new(a, 0, c)?;
            if root.c % modulus == 0 {
                return Err(Error::plain(format!(
                    "hyperelliptic branch root {i} must be p-integral"
                )));
            }
            let reciprocal = invert(model, prime, root.c.rem_euclid(modulus) as u32)?;
            model.spend(1 + residues.len().max(1).ilog2() as usize)?;
            residues.insert(prime.mul(root.a.rem_euclid(modulus) as u32, reciprocal));
            points.push([Element::new(root.a, 0, 1)?, Element::new(root.c, 0, 1)?]);
            inverse.push(reciprocal);
        }
        if residues.len() < 3 {
            return Err(Error::plain(
                "arithmetic hyperelliptic branch configuration requires at least three depth-zero residue directions",
            ));
        }
        let count = roots.len();
        integer_smith::storage(3 * count * (1 + count))?;
        let mut pairs = vec![vec![Contact { depth: 0, unit: 0 }; count]; count];
        let local = Coefficients::Rational(prime);
        for i in 0..count {
            for j in i + 1..count {
                let Some((depth, angular)) = local.bracket(model, points[i], points[j])? else {
                    return Err(Error::plain(format!(
                        "hyperelliptic branch roots {i} and {j} coincide"
                    )));
                };
                let unit = prime.mul(angular, prime.mul(inverse[i], inverse[j]));
                pairs[i][j] = Contact { depth, unit };
                pairs[j][i] = Contact {
                    depth,
                    unit: prime.negate(unit),
                };
            }
        }
        Ok(Self {
            prime,
            leading,
            pairs,
        })
    }

    fn tree(&self, model: &mut Model) -> Result<Rc<Tree>> {
        let n = self.pairs.len();
        model.spend(n * n)?;
        let orders = self
            .pairs
            .iter()
            .flatten()
            .map(|x| i128::from(x.depth))
            .collect::<Vec<_>>();
        Tree::from_pairs(model, n, &orders)
    }

    fn cluster_depth(&self, model: &mut Model, side: &[usize]) -> Result<i64> {
        // In an ultrametric the minimum from any center in the cluster
        // is the same as the minimum over all pairs.
        model.spend(1 + side.len())?;
        Ok(side[1..]
            .iter()
            .map(|&other| self.pairs[side[0]][other].depth)
            .min()
            .expect("a stable split has at least two marks on each side"))
    }

    fn nonsquare(&self, model: &mut Model, mut value: u32) -> Result<bool> {
        let mut power = (self.prime.get() - 1) / 2;
        let mut character = 1;
        while power != 0 {
            model.spend(1)?;
            if power & 1 == 1 {
                character = self.prime.mul(character, value);
            }
            value = self.prime.mul(value, value);
            power >>= 1;
        }
        debug_assert!(character == 1 || character == self.prime.get() - 1);
        Ok(character != 1)
    }

    fn chart(
        &self,
        model: &mut Model,
        inside: &[usize],
        depth: i64,
        leading: u32,
        nonsquare: &mut Option<u32>,
    ) -> Result<Chart> {
        let center = inside[0];
        let mut odd = BTreeSet::new();
        let mut scale = None;
        for &root in inside {
            let contact = self.pairs[root][center];
            if root != center && contact.depth < depth {
                return Err(Error::plain("normalized root is outside its cluster chart"));
            }
            let residue = if contact.depth == depth {
                contact.unit
            } else {
                0
            };
            if residue != 0 && scale.is_none() {
                scale = Some(residue);
            }
            model.spend(1 + odd.len().max(1).ilog2() as usize)?;
            if !odd.remove(&residue) {
                odd.insert(residue);
            }
        }
        let scale = scale.ok_or_else(|| {
            Error::plain("normalized cluster chart has only one finite direction")
        })?;
        let inverse = invert(model, self.prime, scale)?;
        // The first two finite residue directions give coordinates 0 and 1.
        // The discarded even powers are squares, so only odd factors
        // contribute to the leading unit under this coordinate change.
        model.spend(1 + odd.len())?;
        let unit = (0..odd.len()).fold(leading, |unit, _| self.prime.mul(unit, scale));
        let sign = self.nonsquare(model, unit)?;
        let leading = if sign {
            if nonsquare.is_none() {
                for residue in 2..self.prime.get() {
                    if self.nonsquare(model, residue)? {
                        *nonsquare = Some(residue);
                        break;
                    }
                }
            }
            nonsquare.ok_or_else(|| Error::plain("odd prime residue field lacks a nonsquare"))?
        } else {
            1
        };
        Chart::normalized(
            model,
            self.prime,
            leading,
            odd.into_iter().map(|x| self.prime.mul(x, inverse)),
            !inside.len().is_multiple_of(2),
            sign,
        )
    }

    fn charts(&self, model: &mut Model, tree: &Tree) -> Result<(Vec<Chart>, Vec<bool>)> {
        let mut charts: Vec<Option<Chart>> = (0..tree.vertices().len()).map(|_| None).collect();
        let mut nonsquare = None;
        let mut edge_signs = Vec::with_capacity(tree.splits().len());
        let mut incidences = BTreeMap::new();
        for (i, directions) in tree.vertices().iter().enumerate() {
            for direction in directions {
                model
                    .spend(1 + direction.len() * (1 + incidences.len().max(1).ilog2() as usize))?;
                incidences.insert(direction.as_slice(), i);
            }
        }
        for (a, b, _) in tree.splits() {
            let da = self.cluster_depth(model, a)?;
            let db = self.cluster_depth(model, b)?;
            let (inside, outside, depth) = match (da > 0, db > 0) {
                (true, false) => (a, b, da),
                (false, true) => (b, a, db),
                _ => {
                    return Err(Error::plain(
                        "stable branch split lacks a unique positive-depth cluster",
                    ))
                }
            };
            let center = inside[0];
            let mut nu = (inside.len() as i64) * depth;
            let mut unit = self.leading;
            model.spend(1 + outside.len())?;
            for &other in outside {
                let contact = self.pairs[center][other];
                if contact.depth >= depth {
                    return Err(Error::plain(
                        "stable branch cluster has an incorrect outer depth",
                    ));
                }
                nu += contact.depth;
                unit = self.prime.mul(unit, contact.unit);
            }
            if nu & 1 != 0 {
                return Err(Error::plain(format!(
                    "arithmetic hyperelliptic chart for cluster {inside:?} has odd normalization valuation {nu}",
                )));
            }
            let swap = self.nonsquare(model, unit)?;
            let &vertex = incidences
                .get(outside.as_slice())
                .ok_or_else(|| Error::plain("stable branch cluster has no child vertex"))?;
            let chart = self.chart(model, inside, depth, unit, &mut nonsquare)?;
            if charts[vertex].replace(chart).is_some() {
                return Err(Error::plain(
                    "stable branch tree gives a vertex two parent clusters",
                ));
            }
            edge_signs.push(swap);
        }
        if charts.iter().filter(|chart| chart.is_none()).count() != 1 {
            return Err(Error::plain(
                "stable branch tree lacks a unique depth-zero vertex",
            ));
        }
        let all = (0..tree.marks).collect::<Vec<_>>();
        let vertices = charts
            .into_iter()
            .map(|chart| match chart {
                Some(chart) => Ok(chart),
                None => self.chart(model, &all, 0, self.leading, &mut nonsquare),
            })
            .collect::<Result<Vec<_>>>()?;
        Ok((vertices, edge_signs))
    }
}

pub(crate) struct Arithmetic {
    pub(crate) cover: Rc<Cover>,
    vertices: Vec<usize>,
    edges: Vec<usize>,
    frobenius: Vec<Vector>,
    periods: Vec<Vec<i128>>,
    charts: Vec<Chart>,
}

impl Arithmetic {
    pub(crate) fn build(
        model: &mut Model,
        grade: Grade,
        leading: Rational,
        roots: Vec<Rational>,
    ) -> Result<Rc<Self>> {
        let data = BranchData::new(model, grade.prime, leading, &roots)?;
        let tree = data.tree(model)?;
        let (charts, edge_signs) = data.charts(model, &tree)?;
        let cover = Cover::build(model, grade, &tree)?;
        // At multiplier one every semistable edge, including an odd
        // separating edge, must have integral length.
        let lengths = cover.metric(model, 1)?;
        let vertices: Vec<_> = cover
            .vertices
            .iter()
            .enumerate()
            .map(|(i, v)| {
                model.spend(1)?;
                let chart = &charts[v.base];
                let valid = chart.genus.map_or(v.degree == 1 && v.genus == 0, |g| {
                    v.degree == 2 && v.genus == g
                });
                if !valid {
                    return Err(Error::plain(
                        "normalized chart genus disagrees with the hyperelliptic graph",
                    ));
                }
                Ok(if chart.nonsquare { v.deck } else { i })
            })
            .collect::<Result<_>>()?;
        let edges: Vec<_> = cover
            .edges
            .iter()
            .enumerate()
            .map(|(i, e)| if edge_signs[e.base] { e.deck } else { i })
            .collect();
        Self::certify_permutation(model, &cover, &vertices, &edges)?;
        let frobenius = Self::action(model, &cover, &edges, &lengths)?;
        let periods = cover.monodromy(model, 1)?;
        let cells = 7
            + charts.iter().map(Chart::cells).sum::<usize>()
            + vertices.len()
            + edges.len()
            + periods.len() * (3 + periods.len() * 2)
            + frobenius
                .iter()
                .map(|column| 3 + 4 * column.len())
                .sum::<usize>();
        integer_smith::storage(cells)?;
        model.allocate(cells)?;
        Ok(Rc::new(Self {
            cover,
            vertices,
            edges,
            frobenius,
            periods,
            charts,
        }))
    }

    fn certify_permutation(
        model: &mut Model,
        cover: &Cover,
        vertices: &[usize],
        edges: &[usize],
    ) -> Result<()> {
        for (i, &v) in vertices.iter().enumerate() {
            model.spend(1)?;
            if vertices[v] != i || cover.vertices[v].base != cover.vertices[i].base {
                return Err(Error::plain(
                    "invalid arithmetic Frobenius on cover vertices",
                ));
            }
        }
        for (i, &target) in edges.iter().enumerate() {
            model.spend(1)?;
            let before = cover.edges[i];
            let after = cover.edges[target];
            if edges[target] != i
                || after.from != vertices[before.from]
                || after.to != vertices[before.to]
                || after.length != before.length
                || after.dilation != before.dilation
            {
                return Err(Error::plain(
                    "arithmetic Frobenius does not preserve the cover incidences",
                ));
            }
        }
        Ok(())
    }

    fn action(
        model: &mut Model,
        cover: &Cover,
        edges: &[usize],
        lengths: &[i128],
    ) -> Result<Vec<Vector>> {
        let mut images: Vec<Cycle> = Vec::new();
        let mut result = Vec::new();
        for cycle in &cover.cycles {
            model.spend(cycle.len() * (1 + cycle.len().max(1).ilog2() as usize))?;
            let mut image: Cycle = cycle.iter().map(|&(edge, c)| (edges[edge], c)).collect();
            image.sort_unstable();
            let mut column = Vector::new();
            let mut check = vec![0i32; cover.edges.len()];
            for (i, &chord) in cover.chords.iter().enumerate() {
                model.spend(1 + image.len().max(1).ilog2() as usize)?;
                if let Ok(k) = image.binary_search_by_key(&chord, |entry| entry.0) {
                    let c = image[k].1;
                    column.insert(i, i128::from(c));
                    model.spend(1 + cover.cycles[i].len())?;
                    for &(edge, value) in &cover.cycles[i] {
                        check[edge] += c * value;
                    }
                }
            }
            model.spend(check.len())?;
            if !check
                .into_iter()
                .enumerate()
                .filter(|&(_, c)| c != 0)
                .eq(image.iter().copied())
            {
                return Err(Error::plain(
                    "arithmetic Frobenius does not preserve the integral cycle lattice",
                ));
            }
            for (original, changed) in cover.cycles.iter().zip(images.iter().chain([&image])) {
                if period_inner(model, cycle, original, lengths)?
                    != period_inner(model, &image, changed, lengths)?
                {
                    return Err(Error::plain(
                        "arithmetic Frobenius does not preserve monodromy pairing",
                    ));
                }
            }
            images.push(image);
            result.push(column);
        }
        Ok(result)
    }

    pub(crate) fn compare(&self, model: &mut Model, rhs: &Self) -> Result<Ordering> {
        let geometric = self.cover.compare(model, &rhs.cover)?;
        if geometric != Ordering::Equal {
            return Ok(geometric);
        }
        model.spend(
            1 + self.vertices.len()
                + self.edges.len()
                + self.charts.iter().map(Chart::cells).sum::<usize>(),
        )?;
        Ok((&self.vertices, &self.edges, &self.charts).cmp(&(
            &rhs.vertices,
            &rhs.edges,
            &rhs.charts,
        )))
    }

    pub(crate) fn summary(&self, model: &mut Model) -> Result<String> {
        model.spend(1)?;
        Ok(format!(
            "RVHyperArithmetic({}, {}){{genus={}, toric={}}}",
            self.cover.grade.prime.get(),
            self.cover.grade.marks,
            self.cover.grade.genus(),
            self.cover.cycles.len(),
        ))
    }

    fn action_text(&self, model: &mut Model) -> Result<String> {
        let n = self.frobenius.len();
        model.spend(1 + self.vertices.len() + self.edges.len() + n * n)?;
        let mut matrix = vec![vec![0i128; n]; n];
        let mut dual = vec![vec![0i128; n]; n];
        for (j, column) in self.frobenius.iter().enumerate() {
            for (&i, &x) in column {
                matrix[i][j] = x;
                dual[j][i] = x;
            }
        }
        Ok(format!(
            "rvHyperFrobenius = vertices={:?}, edges={:?}, cycles={matrix:?}, dual={dual:?}",
            self.vertices, self.edges,
        ))
    }

    fn fixed_lattice(&self, model: &mut Model) -> Result<Vec<Vector>> {
        let n = self.frobenius.len();
        let mut matrix = Matrix::new(model, n, 2 * n)?;
        for i in 0..n {
            for j in 0..n {
                let x = self.frobenius[i].get(&j).copied().unwrap_or(0) - i128::from(i == j);
                matrix.write(model, i, j, x)?;
                matrix.write(model, i, n + j, -self.periods[i][j])?;
            }
        }
        // Since Q is nonsingular, projecting ker[F^T-I|-Q] to
        // the first n coordinates identifies it with the fixed lattice.
        let smith = integer_smith::reduce(model, matrix, Track::None, Track::Basis)?;
        if smith.factors.len() != n {
            return Err(Error::plain(
                "arithmetic fixed-lattice equation lacks full rank",
            ));
        }
        let mut columns = Vec::with_capacity(n);
        for j in n..2 * n {
            let kernel = smith.columns.basis(j);
            model.spend(1 + kernel.len())?;
            columns.push(kernel.range(..n).map(|(&i, &v)| (i, v)).collect());
        }
        Ok(columns)
    }

    fn fixed_text(&self, model: &mut Model) -> Result<String> {
        let n = self.frobenius.len();
        let columns = self.fixed_lattice(model)?;
        let fixed = reduce_fixed(model, &self.periods, &columns)?;
        let mut factors = Vec::new();
        let mut generators = Vec::new();
        let mut order = 1i128;
        for (i, &d) in fixed.factors.iter().enumerate() {
            if d > 1 {
                order = order.checked_mul(d).ok_or_else(arithmetic)?;
                factors.push(d);
                let mut representative = Vector::new();
                for (&j, &c) in fixed.rows.basis(i) {
                    integer_smith::add_scaled(model, &mut representative, &columns[j], c)?;
                }
                let mut entries = vec![0i128; n];
                for (j, c) in representative {
                    entries[j] = c;
                }
                integer_smith::storage((3 + 2 * n) * (generators.len() + 1))?;
                generators.push(entries);
            }
        }
        Ok(format!(
            "rvHyperFixed = factors={factors:?}, order={order}, representatives={generators:?}"
        ))
    }
}

fn reduce_fixed(
    model: &mut Model,
    periods: &[Vec<i128>],
    columns: &[Vector],
) -> Result<integer_smith::Smith> {
    let n = periods.len();
    let mut matrix = Matrix::new(model, n, n)?;
    for (j, column) in columns.iter().enumerate() {
        for (&i, &x) in column {
            matrix.write(model, i, j, x)?;
        }
    }
    let frame = integer_smith::reduce(model, matrix, Track::Both, Track::Basis)?;
    if frame.factors.len() != n {
        return Err(Error::plain("arithmetic fixed lattice lacks full rank"));
    }
    let mut inclusion = Matrix::new(model, n, n)?;
    // Q is symmetric, so each row also lists a monodromy relation.
    for (j, relation) in periods.iter().enumerate() {
        let mut column = Vector::new();
        for (i, &d) in frame.factors.iter().enumerate() {
            let row = frame.rows.coordinate(i);
            model.spend(1 + row.len())?;
            let mut image = 0i128;
            for (&k, &x) in row {
                let component = x.checked_mul(relation[k]).ok_or_else(arithmetic)?;
                image = image.checked_add(component).ok_or_else(arithmetic)?;
            }
            if image % d != 0 {
                return Err(Error::plain(
                    "monodromy relation is outside the arithmetic fixed lattice",
                ));
            }
            integer_smith::add_scaled(model, &mut column, frame.columns.basis(i), image / d)?;
        }
        for (i, x) in column {
            inclusion.write(model, i, j, x)?;
        }
    }
    let result = integer_smith::reduce(model, inclusion, Track::Basis, Track::None)?;
    if result.factors.len() != n {
        return Err(Error::plain("arithmetic component group is not finite"));
    }
    Ok(result)
}

pub(crate) enum Query {
    Frobenius,
    Fixed,
}

impl Query {
    pub(crate) fn report(self, model: &mut Model, cover: &Arithmetic) -> Result<String> {
        match self {
            Self::Frobenius => cover.action_text(model),
            Self::Fixed => cover.fixed_text(model),
        }
    }
}
