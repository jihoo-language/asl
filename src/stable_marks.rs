use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::fmt;
use std::rc::Rc;

use crate::error::{Error, Result};
use crate::model::Model;
use crate::odd_cobar::Prime;
use crate::prime_linear::invert;
pub(crate) use crate::rv_qp_quad::gauge::Coefficients;
use crate::rv_qp_quad::{Element, Field};
use crate::tropical_matroid::{Grade, Matroid};

pub(crate) fn check_count(count: usize) -> Result<()> {
    if count < 3 {
        return Err(Error::plain(
            "stable marked reduction requires at least three marks",
        ));
    }
    // The same coordinate bound also applies when starting from Plucker data.
    Grade::new(2, count, 2)?;
    Ok(())
}

type Split = (Vec<usize>, Vec<usize>, u128);
type Directions = Vec<Vec<usize>>;

pub(crate) struct Tree {
    pub(crate) marks: usize,
    splits: Vec<Split>,
    vertices: Vec<Directions>,
}

impl Tree {
    pub(crate) fn from_matroid(model: &mut Model, matroid: &Matroid) -> Result<Rc<Self>> {
        let count = matroid.grade.ground;
        check_count(count)?;
        if matroid.grade.rank != 2 || matroid.entries().len() != count * (count - 1) / 2 {
            return Err(Error::plain(
                "stable marked tree requires a full-support rank-two valuated matroid",
            ));
        }
        model.spend(count * count)?;
        let mut values = vec![0i128; count * count];
        for (basis, valuation) in matroid.entries() {
            let i = basis[0];
            let j = basis[1];
            values[i * count + j] = i128::from(*valuation);
            values[j * count + i] = i128::from(*valuation);
        }
        Self::from_pairs(model, count, &values)
    }

    fn from_pairs(model: &mut Model, count: usize, values: &[i128]) -> Result<Rc<Self>> {
        let anchor = count - 1;
        model.spend(count * count)?;
        let mut depths = vec![0i128; count * count];
        for i in 0..anchor {
            for j in (i + 1)..anchor {
                let c =
                    values[i * count + j] - values[i * count + anchor] - values[j * count + anchor];
                depths[i * count + j] = c;
                depths[j * count + i] = c;
            }
        }
        let mut tree = Self {
            marks: count,
            splits: Vec::new(),
            vertices: Vec::new(),
        };
        let root: Vec<_> = (0..anchor).collect();
        tree.add_cluster(model, &depths, &root, None)?;
        model.spend(tree.splits.len() * count * (1 + tree.splits.len().max(1).ilog2() as usize))?;
        tree.splits.sort_unstable();
        model.allocate(tree.cells())?;
        Ok(Rc::new(tree))
    }

    fn add_cluster(
        &mut self,
        model: &mut Model,
        depths: &[i128],
        group: &[usize],
        parent_height: Option<i128>,
    ) -> Result<()> {
        let count = self.marks;
        model.spend(count + 4 * group.len() * group.len())?;
        let depth = |i, j| depths[i * count + j];
        let mut height = i128::MAX;
        for (index, &i) in group.iter().enumerate() {
            for &j in &group[index + 1..] {
                height = height.min(depth(i, j));
            }
        }
        let mut children: Directions = Vec::new();
        let mut ids = vec![0; count];
        for &mark in group {
            let child = children
                .iter()
                .position(|branch| depth(branch[0], mark) > height);
            let id = child.unwrap_or(children.len());
            if let Some(child) = child {
                children[child].push(mark);
            } else {
                children.push(vec![mark]);
            }
            ids[mark] = id;
        }
        for (index, &i) in group.iter().enumerate() {
            for &j in &group[index + 1..] {
                if (ids[i] == ids[j]) != (depth(i, j) > height) {
                    return Err(Error::plain(
                        "rank-two Plücker distances fail the ultrametric cluster condition",
                    ));
                }
            }
        }
        let outside = complement(count, group);
        if let Some(parent) = parent_height {
            let length = height
                .checked_sub(parent)
                .and_then(|x| u128::try_from(x).ok())
                .filter(|&x| x > 0)
                .ok_or_else(|| {
                    Error::plain("stable marked tree edge does not have positive length")
                })?;
            let side = group.to_vec();
            self.splits.push(if side < outside {
                (side, outside.clone(), length)
            } else {
                (outside.clone(), side, length)
            });
        }

        // A cluster has at least two children, so adding the direction
        // towards its parent or anchor always yields a stable component.
        let mut directions = children.clone();
        let insert = directions.partition_point(|d| d[0] < outside[0]);
        directions.insert(insert, outside);
        self.vertices.push(directions);
        for child in children.iter().filter(|child| child.len() > 1) {
            self.add_cluster(model, depths, child, Some(height))?;
        }
        Ok(())
    }

    pub(crate) fn three(model: &mut Model) -> Result<Rc<Self>> {
        Self::from_pairs(model, 3, &[0; 9])
    }

    fn cells(&self) -> usize {
        4 + self.splits.len() * (9 + self.marks)
            + self
                .vertices
                .iter()
                .map(|d| self.marks + 3 * d.len() + 3)
                .sum::<usize>()
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        model.spend(1 + self.cells().min(other.cells()))?;
        Ok((&self.splits, &self.vertices).cmp(&(&other.splits, &other.vertices)))
    }

    pub(crate) fn splits(&self) -> &[Split] {
        &self.splits
    }

    pub(crate) fn vertices(&self) -> &[Directions] {
        &self.vertices
    }

    fn directions(&self, vertex: usize) -> Result<&Directions> {
        self.vertices
            .get(vertex)
            .ok_or_else(|| Error::plain(format!("unknown stable vertex {vertex}")))
    }

    pub(crate) fn summary(&self, model: &mut Model) -> Result<String> {
        model.spend(self.cells())?;
        Ok(format!("RVStableTree({}){{{:?}}}", self.marks, self.splits))
    }
}

fn complement(count: usize, group: &[usize]) -> Vec<usize> {
    let mut marks = group.iter().copied().peekable();
    (0..count)
        .filter(|i| {
            if marks.peek() == Some(i) {
                marks.next();
                false
            } else {
                true
            }
        })
        .collect()
}

#[derive(Clone, Copy)]
struct Local {
    order: i64,
    angular: u32,
}

struct Pairs {
    marks: usize,
    prime: Prime,
    values: Vec<Option<Local>>,
}

impl Pairs {
    fn build(model: &mut Model, field: Coefficients, marks: &[[Element; 2]]) -> Result<Self> {
        let count = marks.len();
        let prime = field.prime();
        model.spend(count * count)?;
        let mut values = vec![None; count * count];
        for (i, &point) in marks.iter().enumerate() {
            for (j, &later) in marks.iter().enumerate().skip(i + 1) {
                let Some((order, angular)) = field.bracket(model, point, later)? else {
                    return Err(Error::plain(format!(
                        "projective marks {i} and {j} collide exactly"
                    )));
                };
                values[i * count + j] = Some(Local { order, angular });
                values[j * count + i] = Some(Local {
                    order,
                    angular: prime.negate(angular),
                });
            }
        }
        Ok(Self {
            marks: count,
            prime,
            values,
        })
    }

    fn orders(&self, model: &mut Model) -> Result<Vec<i128>> {
        model.spend(self.values.len())?;
        Ok(self
            .values
            .iter()
            .map(|x| x.map_or(0, |v| i128::from(v.order)))
            .collect())
    }

    fn get(&self, i: usize, j: usize) -> Option<Local> {
        self.values[i * self.marks + j]
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Residue {
    Infinity,
    Finite(u32),
}

impl fmt::Display for Residue {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Infinity => write!(out, "infinity"),
            Self::Finite(value) => write!(out, "{value}"),
        }
    }
}

struct Normalization<'a> {
    pairs: &'a Pairs,
    infinity: usize,
    zero: usize,
    order: i128,
    angular: u32,
}

impl<'a> Normalization<'a> {
    fn new(model: &mut Model, pairs: &'a Pairs, [infinity, zero, one]: [usize; 3]) -> Result<Self> {
        let numerator = pairs
            .get(one, infinity)
            .expect("distinct normalization marks");
        let denominator = pairs.get(one, zero).expect("distinct normalization marks");
        let unit_inverse = invert(model, pairs.prime, denominator.angular)?;
        Ok(Self {
            pairs,
            infinity,
            zero,
            order: i128::from(numerator.order) - i128::from(denominator.order),
            angular: pairs.prime.mul(numerator.angular, unit_inverse),
        })
    }

    fn reduce(&self, model: &mut Model, mark: usize) -> Result<Residue> {
        model.spend(1)?;
        let Some(numerator) = self.pairs.get(mark, self.zero) else {
            return Ok(Residue::Finite(0));
        };
        let Some(denominator) = self.pairs.get(mark, self.infinity) else {
            return Ok(Residue::Infinity);
        };
        // [mark, zero] [one, infinity] / ([mark, infinity] [one, zero]).
        let total = self.order + i128::from(numerator.order) - i128::from(denominator.order);
        match total.cmp(&0) {
            Ordering::Greater => Ok(Residue::Finite(0)),
            Ordering::Less => Ok(Residue::Infinity),
            Ordering::Equal => {
                let prime = self.pairs.prime;
                let unit_inverse = invert(model, prime, denominator.angular)?;
                Ok(Residue::Finite(prime.mul(
                    self.angular,
                    prime.mul(numerator.angular, unit_inverse),
                )))
            }
        }
    }
}

pub(crate) struct Reduction {
    pub(crate) prime: u32,
    pub(crate) tree: Rc<Tree>,
    residues: Vec<Vec<Residue>>,
}

impl Reduction {
    pub(crate) fn rational(
        model: &mut Model,
        prime: Prime,
        raw: Vec<[i128; 2]>,
    ) -> Result<Rc<Self>> {
        check_count(raw.len())?;
        let mut points = Vec::with_capacity(raw.len());
        for [mut x, mut y] in raw {
            model.spend(1)?;
            let divisor = crate::rv_ore::gcd(x.unsigned_abs(), y.unsigned_abs());
            if divisor > 1 {
                if x != 0 {
                    x = crate::rv_phodge::divide_content(x, divisor);
                }
                if y != 0 {
                    y = crate::rv_phodge::divide_content(y, divisor);
                }
            }
            points.push([Element::new(x, 0, 1)?, Element::new(y, 0, 1)?]);
        }
        Self::build(model, Coefficients::Rational(prime), &points)
    }

    pub(crate) fn quadratic(
        model: &mut Model,
        field: Field,
        points: &[[Element; 2]],
    ) -> Result<Rc<Self>> {
        Self::build(model, Coefficients::Quadratic(field), points)
    }

    fn build(model: &mut Model, field: Coefficients, points: &[[Element; 2]]) -> Result<Rc<Self>> {
        check_count(points.len())?;
        // Detect a zero coordinate pair before testing pairwise determinants.
        model.spend(points.len())?;
        if let Some((mark, _)) = points
            .iter()
            .enumerate()
            .find(|(_, pair)| pair.iter().all(|x| x.is_zero()))
        {
            return Err(Error::plain(format!(
                "mark {mark} has zero homogeneous coordinates"
            )));
        }
        let pairs = Pairs::build(model, field, points)?;
        let orders = pairs.orders(model)?;
        let tree = Tree::from_pairs(model, points.len(), &orders)?;
        let mut residues = Vec::with_capacity(tree.vertices.len());
        let mut cells = 4;
        for (vertex, directions) in tree.vertices.iter().enumerate() {
            let normalization = Normalization::new(
                model,
                &pairs,
                [directions[0][0], directions[1][0], directions[2][0]],
            )?;
            let mut seen = BTreeSet::new();
            let mut row = Vec::with_capacity(directions.len());
            for branch in directions {
                let residue = normalization.reduce(model, branch[0])?;
                for &i in &branch[1..] {
                    if normalization.reduce(model, i)? != residue {
                        return Err(Error::plain(format!("stable residue is not constant on the branch of mark {i} at vertex {vertex}")));
                    }
                }
                model.spend(1 + seen.len().max(1).ilog2() as usize)?;
                if !seen.insert(residue) {
                    return Err(Error::plain(format!(
                        "distinct directions collide on the residue component at vertex {vertex}"
                    )));
                }
                row.push(residue);
            }
            cells += 3 + row.len();
            residues.push(row);
        }
        model.allocate(cells)?;
        Ok(Rc::new(Self {
            prime: pairs.prime.get(),
            tree,
            residues,
        }))
    }

    pub(crate) fn three(model: &mut Model, prime: u32) -> Result<Rc<Self>> {
        let tree = Tree::three(model)?;
        model.allocate(10)?;
        Ok(Rc::new(Self {
            prime,
            tree,
            residues: vec![vec![
                Residue::Infinity,
                Residue::Finite(0),
                Residue::Finite(1),
            ]],
        }))
    }

    pub(crate) fn compare(&self, model: &mut Model, other: &Self) -> Result<Ordering> {
        let tree = self.tree.compare(model, &other.tree)?;
        if tree != Ordering::Equal {
            return Ok(tree);
        }
        model.spend(self.residues.iter().map(|r| 1 + r.len()).sum())?;
        Ok(self.residues.cmp(&other.residues))
    }

    fn components(&self, model: &mut Model) -> Result<String> {
        let mut components = Vec::new();
        for (directions, row) in self.tree.vertices.iter().zip(&self.residues) {
            model.spend(self.tree.marks + 1 + row.len())?;
            let branches: Vec<_> = directions
                .iter()
                .zip(row)
                .map(|(branch, residue)| format!("({branch:?}, {residue})"))
                .collect();
            components.push(format!("[{}]", branches.join(", ")));
        }
        Ok(format!("[{}]", components.join(", ")))
    }

    pub(crate) fn summary(&self, model: &mut Model) -> Result<String> {
        Ok(format!(
            "RVStableMarks({}, {}){{splits={:?}, components={}}}",
            self.prime,
            self.tree.marks,
            self.tree.splits,
            self.components(model)?
        ))
    }

    fn at(&self, model: &mut Model, vertex: usize, mark: usize) -> Result<Residue> {
        if mark >= self.tree.marks {
            return Err(Error::plain("stable mark is outside the ground set"));
        }
        let directions = self.tree.directions(vertex)?;
        model.spend(1 + self.tree.marks)?;
        let i = directions
            .iter()
            .position(|branch| branch.binary_search(&mark).is_ok())
            .expect("the directions partition the marks");
        Ok(self.residues[vertex][i])
    }
}

pub(crate) enum Query {
    Splits,
    Directions(usize),
    Components,
    Residue(usize, usize),
}

impl Query {
    pub(crate) fn requires_tree(&self) -> bool {
        matches!(self, Self::Splits | Self::Directions(_))
    }

    pub(crate) fn tree(self, model: &mut Model, value: &Tree) -> Result<String> {
        model.spend(value.cells())?;
        Ok(match self {
            Self::Splits => format!("rvStableSplits = {:?}", value.splits),
            Self::Directions(vertex) => format!(
                "rvStableDirections[{vertex}] = {:?}",
                value.directions(vertex)?
            ),
            _ => unreachable!("checked stable tree operation"),
        })
    }

    pub(crate) fn curve(self, model: &mut Model, value: &Reduction) -> Result<String> {
        Ok(match self {
            Self::Components => format!("rvStableComponents = {}", value.components(model)?),
            Self::Residue(vertex, mark) => format!(
                "rvStableResidue[{vertex}, {mark}] = {}",
                value.at(model, vertex, mark)?
            ),
            _ => unreachable!("checked stable residue operation"),
        })
    }
}
